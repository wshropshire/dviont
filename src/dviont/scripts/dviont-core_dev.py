#!/usr/bin/env python3

import argparse
import os
import subprocess
import logging
import sys
from directory_management import PipelineManager, run_command

def get_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Core genome alignment pipeline.")
    parser.add_argument("-o", "--output_dir", required=True, help="Output directory")
    parser.add_argument("-r", "--ref", required=True, help="Reference genome")
    parser.add_argument("-v", "--vcf", required=True, help="Input VCF file")
    parser.add_argument("-b", "--bam", required=True, help="Input BAM file")
    parser.add_argument("-t", "--threads", type=int, default=2, help="Number of threads")
    parser.add_argument("-s", "--sample", default="SAMPLE", help="Sample name")
    parser.add_argument("--mincovfrac", type=float, default=0.2, help="Minimum coverage fraction threshold (default: 0.2)")
#    parser.add_argument("--snpfrac", type=float, default=0.9, help="Minimum allele frequency for SNP filtering (default: 0.9)")
    parser.add_argument("--auto_mask", action="store_true", help="Automatically identify repeat regions using NUCmer and merge with coverage-based masking")
    return parser.parse_args()

def extract_snps_and_dels(vcf, output_vcf): # remove snpfrac
    """Extract SNPs and deletions from the VCF file."""
 #   command = f"bcftools view -i 'TYPE=\"snp\" || (strlen(REF) > strlen(ALT))' {vcf} | bcftools filter -i 'AF>={snpfrac}' -Oz -o {output_vcf}"
    command = f"bcftools view -i 'TYPE=\"snp\" || (strlen(REF) > strlen(ALT))' {vcf} | bcftools view -e 'GT=\"1/0\"' -Oz -o {output_vcf}"
    run_command(command)
    run_command(f"bcftools index {output_vcf}")

def calculate_average_depth(bam_file):
    """Calculate the average sequencing depth."""
    command = f"samtools depth -a {bam_file}"
    depth_data = run_command(command).splitlines()
    
    total_depth, total_positions = 0, 0
    for line in depth_data:
        _, _, depth = line.split()
        depth = int(depth)
        if depth > 0:
            total_depth += depth
            total_positions += 1

    avg_depth = total_depth / total_positions if total_positions else 0
    logging.info(f"Calculated average depth: {avg_depth:.2f}")
    return avg_depth

def generate_lowcov_masking_bed(bam_file, avg_depth, output_bed, mincovfrac):
    """Generate BED file marking missing and low-coverage regions with N."""
    try:
        threshold = avg_depth * mincovfrac
        logging.info(f"Identifying missing regions and positions with coverage depths less than {mincovfrac} threshold")
        command = f"samtools depth -a {bam_file}"
        depth_data = run_command(command)
        if depth_data is None:
            return

        masked_regions = []
        current_region = None

        for line in depth_data.splitlines():
            try:
                chrom, pos, depth = line.split()
                pos, depth = int(pos) - 1, int(depth)
                if depth == 0 or depth < threshold:
                    if current_region is None:
                        current_region = [chrom, pos, pos]
                    else:
                        current_region[2] = pos
                else:
                    if current_region:
                        masked_regions.append(current_region)
                        current_region = None
            except ValueError:
                logging.warning(f"Skipping malformed depth line: {line}")
                continue

        if current_region:
            masked_regions.append(current_region)

        with open(output_bed, "w") as bed:
            for region in masked_regions:
                bed.write(f"{region[0]}\t{region[1]}\t{region[2]+1}\n")
        
        logging.info(f"Generated correctly formatted BED file: {output_bed}")
    except Exception as e:
        logging.error(f"Error generating BED file: {e}")

def run_nucmer(reference, output_dir, prefix):
    """Run NUCmer for All-vs-All comparison to identify repeat regions."""
    nucmer_cmd = f"nucmer --maxmatch --nosimplify {reference} {reference} -p {output_dir}/{prefix}"
    run_command(nucmer_cmd)
    logging.info("NUCmer All-vs-All comparison completed.")

def generate_auto_masking_bed(delta_file, bed_file):
    """Generate BED file of repeat regions, removing the first two lines and duplicates."""
    temp_bed_file = bed_file + ".tmp"
    show_coords_cmd = f"show-coords -r -T {delta_file} | awk '{{if ($1 != $3 && $2 != $4) print $0}}' | awk '{{print $8\"\t\"$1\"\t\"$2}}' > {temp_bed_file}"
    run_command(show_coords_cmd)
    
    try:
        with open(temp_bed_file, "r") as infile:
            lines = infile.readlines()[2:]
        unique_lines = list(set(lines))
        
        with open(bed_file, "w") as outfile:
            outfile.writelines(unique_lines)
        
        os.remove(temp_bed_file)
        logging.info(f"Filtered and generated auto-mask BED file {bed_file}.")
    except Exception as e:
        logging.error(f"Error processing auto-mask BED file: {e}")

def merge_bed_files(bed1, bed2, output_bed):
    """Merge two BED files using bedtools to take the union of overlapping and unique regions."""
    merge_cmd = f"cat {bed1} {bed2} | sort -k1,1 -k2,2n | bedtools merge > {output_bed}"
    run_command(merge_cmd)
    logging.info(f"Merged BED files into {output_bed}.")

def generate_consensus(sample, ref, output_vcf, output_fasta, output_bed):
    """Generate consensus FASTA from SNP VCF with proper deletion and BED masking."""
    temp_fasta = output_fasta + ".tmp"
    command = f"bcftools consensus -f {ref} -o {temp_fasta} --mark-del N -m {output_bed} {output_vcf}"
    run_command(command)
    
    with open(temp_fasta, "r") as infile, open(output_fasta, "w") as outfile:
        for line in infile:
            if line.startswith(">"):
                outfile.write(f">{sample}\n")
            else:
                outfile.write(line)
    
    os.remove(temp_fasta)
    logging.info(f"Generated consensus FASTA with header replaced to sample name: {output_fasta}")

def main():
    """Main pipeline execution."""
    args = get_arguments()
    
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Initialize logging using PipelineManager
    manager = PipelineManager(args.output_dir, args.sample)
    manager.setup_logging()
    
    snps_and_dels_vcf = os.path.join(args.output_dir, f"{args.sample}_snps_and_dels.vcf.gz")
    extract_snps_and_dels(args.vcf, snps_and_dels_vcf) # args.snpfrac removed
    
    avg_depth = calculate_average_depth(args.bam)
    bed_file = os.path.join(args.output_dir, f"{args.sample}_lowcov_mask.bed")
    generate_lowcov_masking_bed(args.bam, avg_depth, bed_file, args.mincovfrac)
    
    if args.auto_mask:
        nucmer_prefix = f"{args.sample}_nucmer"
        delta_file = f"{args.output_dir}/{nucmer_prefix}.delta"
        repeat_bed_file = os.path.join(args.output_dir, f"{args.sample}_repeat_regions.bed")
        merged_bed_file = os.path.join(args.output_dir, f"{args.sample}_merged_mask.bed")
        
        run_nucmer(args.ref, args.output_dir, nucmer_prefix)
        generate_auto_masking_bed(delta_file, repeat_bed_file)
        merge_bed_files(bed_file, repeat_bed_file, merged_bed_file)
        bed_file = merged_bed_file
    
    consensus_fasta = os.path.join(args.output_dir, f"{args.sample}_consensus.fasta")
    generate_consensus(args.sample, args.ref, snps_and_dels_vcf, consensus_fasta, bed_file)
    
    logging.info("Pipeline completed successfully.")

if __name__ == "__main__":
    main()
