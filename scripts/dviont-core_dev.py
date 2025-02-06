#!/usr/bin/env python3

import argparse
import os
import subprocess
import logging
import time
import re

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', handlers=[
    logging.StreamHandler(),
    logging.FileHandler('pipeline.log')
])

def run_command(command):
    """Runs a shell command and logs output."""
    try:
        result = subprocess.run(command, shell=True, check=True, text=True, capture_output=True)
        logging.info(f"Command succeeded: {command}")
        return result.stdout
    except subprocess.CalledProcessError as e:
        logging.error(f"Error running command: {command}\n{e.stderr}")
        raise

def extract_snps(vcf, output_vcf):
    """Extract SNPs from the VCF file and index it."""
    command = f"bcftools view -v snps -Oz -o {output_vcf} {vcf}"
    run_command(command)
    
    # Index the VCF file
    index_command = f"bcftools index {output_vcf}"
    run_command(index_command)

def calculate_average_depth(bam_file):
    """Calculate the average depth excluding positions with zero coverage."""
    command = f"samtools depth -a {bam_file}"
    depth_data = run_command(command).splitlines()
    
    total_depth, total_positions = 0, 0
    for line in depth_data:
        _, _, depth = line.split()
        depth = int(depth)
        if depth > 0:
            total_depth += depth
            total_positions += 1
    
    return total_depth / total_positions if total_positions else 0

def generate_masking_bed(bam_file, avg_depth, output_bed):
    """Generate BED file with positions of 0 or <10% coverage."""
    command = f"samtools depth -a {bam_file}"
    depth_data = run_command(command).splitlines()
    
    threshold = avg_depth * 0.1
    with open(output_bed, 'w') as bed:
        for line in depth_data:
            chrom, pos, depth = line.split()
            if int(depth) == 0 or int(depth) < threshold:
                bed.write(f"{chrom}\t"-"\t{pos}\n")

def generate_consensus(ref, snp_vcf, output_fasta):
    """Generate consensus FASTA from SNP VCF."""
    command = f"bcftools consensus -f {ref} -o {output_fasta} {snp_vcf}"
    run_command(command)

def remove_insertions(consensus_fasta, cleaned_fasta):
    """Remove lowercase (insertions) from the consensus FASTA."""
    with open(consensus_fasta, 'r') as infile, open(cleaned_fasta, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                outfile.write(line)
            else:
                outfile.write(re.sub(r'[a-z]', '', line))

def mask_consensus(consensus_fasta, bed_file, masked_fasta):
    """Mask the consensus file using a BED file."""
    command = f"bedtools maskfasta -fi {consensus_fasta} -bed {bed_file} -fo {masked_fasta}"
    run_command(command)

def main(output_dir, ref, vcf, bam, threads, sample):
    """Main workflow."""
    os.makedirs(output_dir, exist_ok=True)
    start_time = time.time()
    
    snp_vcf = os.path.join(output_dir, f"{sample}_snps.vcf.gz")
    extract_snps(vcf, snp_vcf)
    
    avg_depth = calculate_average_depth(bam)
    bed_file = os.path.join(output_dir, f"{sample}_mask.bed")
    generate_masking_bed(bam, avg_depth, bed_file)
    
    consensus_fasta = os.path.join(output_dir, f"{sample}_consensus.fasta")
    generate_consensus(ref, snp_vcf, consensus_fasta)
    
    cleaned_fasta = os.path.join(output_dir, f"{sample}_cleaned.fasta")
    remove_insertions(consensus_fasta, cleaned_fasta)
    
    masked_fasta = os.path.join(output_dir, f"{sample}_masked_consensus.fasta")
    mask_consensus(cleaned_fasta, bed_file, masked_fasta)
    
    logging.info(f"Pipeline completed for {sample} in {time.time() - start_time:.2f} seconds.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Core genome alignment pipeline.")
    parser.add_argument("-o", "--output_dir", required=True)
    parser.add_argument("-r", "--ref", required=True)
    parser.add_argument("-v", "--vcf", required=True)
    parser.add_argument("-b", "--bam", required=True)
    parser.add_argument("-t", "--threads", type=int, default=2)
    parser.add_argument("-s", "--sample", default="SAMPLE")
    args = parser.parse_args()
    main(args.output_dir, args.ref, args.vcf, args.bam, args.threads, args.sample)
