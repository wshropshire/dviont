#!/usr/bin/env python3

import argparse
import os
import shutil
import logging
import time
import sys
from .directory_management import PipelineManager
from .ref_format import determine_ref_format
from .extract_fasta_and_gbk import extract_fasta_and_gbk
from .minimap2 import run_minimap2_alignment
from .clair3_module import run_clair3
from .snpEff_module import run_snpEff
from .vcf_processor import VCFProcessor
from .version import __version__

# Ensure relative imports work in a module-based execution context
if __name__ == "__main__" and not hasattr(sys, 'argv'):
    sys.argv = [__file__]
    __package__ = "dviont"

def get_arguments():
    """Parse assembler arguments"""
    parser = argparse.ArgumentParser(description="The DNA Variant Identification using ONT (dviONT) Pipeline.")
    parser.add_argument("-o", "--output_dir", required=True, help="Output directory for results")
    parser.add_argument("-r", "--ref", required=True, help="Path to reference genome file (FASTA or GBK)")
    parser.add_argument("-i", "--reads", required=True, help="Path to reads file (FASTQ)")
    parser.add_argument("-t", "--threads", type=int, default=2, help="Number of threads to use (default: 2)")
    parser.add_argument("-m", "--model_name", default="r1041_e82_400bps_sup_v430_bacteria_finetuned", help="Model name for Clair3 (default: r1041_e82_400bps_sup_v430_bacteria_finetuned)")
    parser.add_argument("-s", "--sample", default="SAMPLE", help="Sample name")
    parser.add_argument("-p", "--model_path", default=None, help="Path to Clair3 model (optional)")
    parser.add_argument("-v", "--version", action="version", version=f"dviONT v{__version__}", help="Show version and exit")

    args = parser.parse_args()
    return args 

def main():
    """Main workflow: Extract FASTA/GBK first, run minimap2 followed by Clair3 and SnpEff (if GBK provided)."""
    # Record log time
    start_time = time.time()  # Record start time for wall time calculation

    # Get arguments for script
    args = get_arguments()
    output_dir = args.output_dir  
    ref = args.ref
    reads = args.reads
    model_path = args.model_path
    model_name = args.model_name
    threads = args.threads
    sample = args.sample

    try:
        # Step 1: Create output directories and initialize logging using PipelineManager
        pipeline_manager = PipelineManager(output_dir, sample)
        ref_dir = pipeline_manager.create_output_directory()
        log_file = pipeline_manager.get_log_file()

        # Step 2: Determine reference file format (FASTA or GenBank)
        ref_fmt = determine_ref_format(ref)
        logging.info(f"Reference format determined: {ref_fmt}")

        # Step 3: Extract FASTA and GFF from the reference file
        fasta_out = extract_fasta_and_gbk(ref, ref_dir, ref_fmt, output_dir)
        logging.info(f"FASTA file extracted: {fasta_out}")

        # Step 4: Run Minimap2 to align reads
        bam_output = run_minimap2_alignment(fasta_out, reads, threads, output_dir, sample)

        # Step 5: Run Clair3 to generate a VCF file
        vcf_out, consensus_path = run_clair3(output_dir, fasta_out, bam_output, threads, model_name, sample, model_path)

        # Step 6: Optionally run SnpEff for annotation (if reference is GenBank)
        if ref_fmt == "genbank":
            snp_eff_annotated_vcf = run_snpEff(output_dir, vcf_out, fasta_out, sample)
            vcf_file_to_process = snp_eff_annotated_vcf
        else:
            logging.warning("Reference is FASTA; skipping SnpEff annotation.")
            vcf_file_to_process = vcf_out

        # Step 7: Process the VCF file using VCFProcessor
        logging.info(f"Processing VCF file: {vcf_file_to_process}")
        vcf_processor = VCFProcessor(
            vcf_file=vcf_file_to_process,
            ref_fmt=ref_fmt,
            output_dir=output_dir,
            sample=sample,
            genbank_file=ref if ref_fmt == "genbank" else None
        )
        vcf_processor.parse_vcf()
        logging.info(f"✅ VCF processing complete. Results saved in {output_dir}")

        # End of pipeline message
        end_time = time.time()
        wall_time = end_time - start_time
        logging.info(f"🎉 dviONT pipeline completed successfully for sample '{sample}'.")
        logging.info(f"🕒 Total compute (wall time): {wall_time:.2f} seconds.")
        logging.info(f"📜 Logs available at: {log_file}")

    except Exception as e:
        logging.error(f"Error occurred: {e}", exc_info=True)
        print(f"❌ An error occurred. Check the log file for details: {log_file}")
        sys.exit(1)


if __name__ == "__main__":main()
