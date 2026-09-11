#!/usr/bin/env python3

import argparse
import logging
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

from .version import __version__


PRESETS = ["ont-legacy", "ont-q20", "pb-clr", "pb-hifi", "asm"]


def add_common_arguments(parser):
    parser.add_argument("-r", "--ref", required=True, help="Reference genome (FASTA or GBK)")
    parser.add_argument("-t", "--threads", type=int, default=2, help="Number of threads (default: 2)")
    parser.add_argument("-m", "--model-name", default="r1041_e82_400bps_sup_v430_bacteria_finetuned", help="Clair3 model name")
    parser.add_argument(
        "-p", "--model-path", default=None,
        help="Path to Clair3 model directory (default: found automatically from the model name)",
    )
    parser.add_argument("--preset", choices=PRESETS, default="ont-q20", help="Alignment preset (default: ont-q20)")
    parser.add_argument("--aligner", choices=["minimap2", "winnowmap"], default="minimap2", help="Read aligner (default: minimap2)")


def build_parser():
    parser = argparse.ArgumentParser(description="The DNA Variant Identification using ONT (dviONT) Pipeline.")
    parser.add_argument("-v", "--version", action="version", version=f"dviONT v{__version__}")
    subparsers = parser.add_subparsers(dest="command", required=True)

    call = subparsers.add_parser("call", help="Run the standard single-sample DviONT workflow")
    add_common_arguments(call)
    call.add_argument("-o", "--output-dir", required=True, help="Output directory for results")
    call.add_argument("-i", "--reads", required=True, help="Reads file (FASTQ)")
    call.add_argument("-s", "--sample", default="SAMPLE", help="Sample name (default: SAMPLE)")

    cohort = subparsers.add_parser("cohort", help="Call multiple samples and build a SNP alignment")
    add_common_arguments(cohort)
    cohort.add_argument("--reads-list", required=True, help="Tab-separated sample_id and reads_path file")
    cohort.add_argument("--out", required=True, help="Cohort output directory")
    return parser


def parse_reads_list(path):
    samples = []
    seen = set()
    with open(path) as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 2 or not all(fields):
                raise ValueError(f"{path}:{line_number}: expected sample_id<TAB>reads_path")
            sample, reads = fields
            if sample in seen:
                raise ValueError(f"{path}:{line_number}: duplicate sample ID {sample}")
            if not os.path.isfile(reads):
                raise FileNotFoundError(f"Reads file for {sample} does not exist: {reads}")
            seen.add(sample)
            samples.append((sample, reads))
    if not samples:
        raise ValueError(f"No samples found in {path}")
    return samples


def run_call(args):
    """Run the standard single-sample workflow."""
    from .clair3_module import run_clair3
    from .directory_management import PipelineManager
    from .extract_fasta_and_gbk import extract_fasta_and_gbk
    from .minimap2 import run_minimap2_alignment
    from .ref_format import determine_ref_format
    from .snpEff_module import run_snpEff
    from .vcf_processor import VCFProcessor
    from .winnowmap import run_winnowmap_alignment

    start_time = time.time()
    pipeline_manager = PipelineManager(args.output_dir, args.sample, model_path=args.model_path)
    ref_dir = pipeline_manager.create_output_directory()
    log_file = pipeline_manager.get_log_file()

    ref_fmt = determine_ref_format(args.ref)
    logging.info("Reference format determined: %s", ref_fmt)
    fasta_out = extract_fasta_and_gbk(args.ref, ref_dir, ref_fmt, args.output_dir)

    if args.aligner == "winnowmap":
        bam_output = run_winnowmap_alignment(
            fasta_out, args.reads, args.threads, args.output_dir, args.sample, args.preset
        )
    else:
        bam_output = run_minimap2_alignment(
            fasta_out, args.reads, args.threads, args.output_dir, args.sample, args.preset
        )
    if not bam_output:
        raise RuntimeError(f"{args.aligner} alignment failed")

    result = run_clair3(
        args.output_dir, fasta_out, bam_output, args.threads,
        args.model_name, args.sample, args.model_path,
    )
    if not result:
        raise RuntimeError("Clair3 calling failed")
    final_vcf, _ = result

    if ref_fmt == "genbank":
        vcf_file_to_process = run_snpEff(
            args.output_dir, final_vcf, fasta_out, args.sample
        )
    else:
        logging.warning("Reference is FASTA; skipping SnpEff annotation.")
        vcf_file_to_process = final_vcf

    processor = VCFProcessor(
        vcf_file=vcf_file_to_process,
        ref_fmt=ref_fmt,
        output_dir=args.output_dir,
        sample=args.sample,
        genbank_file=args.ref if ref_fmt == "genbank" else None,
    )
    processor.parse_vcf()
    logging.info("DviONT completed for '%s' in %.2f seconds", args.sample, time.time() - start_time)
    logging.info("Logs available at: %s", log_file)
    return final_vcf, fasta_out


def require_executable(name):
    if shutil.which(name) is None:
        raise RuntimeError(f"Required executable not found on PATH: {name}")


def fasta_sequence_length(path):
    """Return the total number of bases in a FASTA file."""
    length = 0
    records = 0
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                records += 1
            else:
                length += len(line.strip())
    if records == 0:
        raise ValueError(f"No FASTA records found in {path}")
    return length


def write_consensus(vcf, reference, sample, output):
    """Write one full-reference-length consensus record for a cohort sample."""
    command = [
        "bcftools", "consensus", "-f", str(reference), "-s", sample, str(vcf)
    ]
    result = subprocess.run(command, check=True, text=True, capture_output=True)
    sequence = "".join(
        line.strip() for line in result.stdout.splitlines() if not line.startswith(">")
    )
    if not sequence:
        raise ValueError(f"bcftools consensus produced no sequence for {sample}")
    with open(output, "w") as handle:
        handle.write(f">{sample}\n")
        for offset in range(0, len(sequence), 60):
            handle.write(sequence[offset:offset + 60] + "\n")
    return len(sequence)


def run_cohort(args):
    """Run ordinary DviONT calls and combine their final processed VCFs."""
    require_executable("bcftools")
    require_executable("snp-dists")
    samples = parse_reads_list(args.reads_list)
    out = Path(args.out)
    calls_dir = out / "calls"
    alignments_dir = out / "alignments"
    consensus_dir = alignments_dir / "consensus_snps"
    cohort_vcfs_dir = out / "cohort_vcfs"
    distances_dir = out / "distances"
    calls_dir.mkdir(parents=True, exist_ok=True)
    consensus_dir.mkdir(parents=True, exist_ok=True)
    cohort_vcfs_dir.mkdir(parents=True, exist_ok=True)
    distances_dir.mkdir(parents=True, exist_ok=True)

    vcfs = []
    cohort_reference = None
    for sample, reads in samples:
        call_args = argparse.Namespace(**vars(args))
        call_args.output_dir = str(calls_dir / sample)
        call_args.reads = reads
        call_args.sample = sample
        final_vcf, sample_reference = run_call(call_args)
        if not os.path.isfile(final_vcf):
            raise FileNotFoundError(f"Final VCF not found: {final_vcf}")
        subprocess.run(["bcftools", "index", "-f", final_vcf], check=True)
        vcfs.append(final_vcf)
        if cohort_reference is None:
            cohort_reference = sample_reference

    merged = cohort_vcfs_dir / "cohort_merged.vcf.gz"
    normalized = cohort_vcfs_dir / "cohort_merged.norm.vcf.gz"
    subprocess.run(["bcftools", "merge", "-m", "none", "-Oz", "-o", str(merged), *vcfs], check=True)
    subprocess.run(["bcftools", "index", "-f", str(merged)], check=True)
    subprocess.run(["bcftools", "norm", "-f", cohort_reference, "-m", "-any", "-Oz", "-o", str(normalized), str(merged)], check=True)
    subprocess.run(["bcftools", "index", "-f", str(normalized)], check=True)
    filtered = cohort_vcfs_dir / "cohort_merged.snps.vcf.gz"
    subprocess.run(["bcftools", "view", "-m2", "-M2", "-v", "snps", "-Oz", "-o", str(filtered), str(normalized)], check=True)
    subprocess.run(["bcftools", "index", "-f", str(filtered)], check=True)

    sample_names = [sample for sample, _ in samples]
    reference_length = fasta_sequence_length(cohort_reference)
    consensus_fastas = []
    for sample in sample_names:
        consensus = consensus_dir / f"{sample}.fasta"
        consensus_length = write_consensus(filtered, cohort_reference, sample, consensus)
        if consensus_length != reference_length:
            raise ValueError(
                f"Consensus length for {sample} is {consensus_length}; "
                f"expected reference length {reference_length}"
            )
        consensus_fastas.append(consensus)

    alignment = alignments_dir / "cohort.snp_alignment.fasta"
    with open(alignment, "wb") as destination:
        for consensus in consensus_fastas:
            with open(consensus, "rb") as source:
                shutil.copyfileobj(source, destination)

    distances = distances_dir / "cohort.snp_distance_matrix.tsv"
    with open(distances, "w") as handle:
        subprocess.run(["snp-dists", str(alignment)], check=True, stdout=handle)


def main():
    parser = build_parser()
    args = parser.parse_args()
    try:
        from .clair3_models import resolve_model_path
        # Resolve the Clair3 model up front so a missing model fails before alignment starts.
        args.model_path = resolve_model_path(args.model_name, args.model_path)
        if args.command == "call":
            run_call(args)
        else:
            run_cohort(args)
    except Exception as error:
        logging.error("Error occurred: %s", error, exc_info=True)
        print(f"An error occurred: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
