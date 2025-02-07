import os
import shutil
import logging
from Bio import SeqIO
import subprocess
import gzip

def extract_fasta_and_gbk(reference, ref_dir, ref_fmt, output_dir):
    """Extract FASTA and GenBank files, ensuring correct formatting."""
    try:
        if not os.path.exists(reference):
            logging.error(f"Reference file not found: {reference}")
        
        os.makedirs(os.path.join(ref_dir, "genomes"), exist_ok=True)
        os.makedirs(os.path.join(ref_dir, "ref"), exist_ok=True)

        if ref_fmt == "genbank":
            # Define paths
            genes_gbk_path = os.path.join(ref_dir, "ref", "genes.gbk")
            genes_gbk_gz_path = genes_gbk_path + ".gz"

            # Copy and gzip the GenBank file
            shutil.copy(reference, genes_gbk_path)
            with open(genes_gbk_path, 'rb') as f_in, gzip.open(genes_gbk_gz_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
            
            # Remove uncompressed file after gzipping
            os.remove(genes_gbk_path)
            logging.info(f"Copied and gzipped GenBank file: {genes_gbk_gz_path}")

            # Convert GenBank to FASTA
            fasta_out = os.path.join(ref_dir, "ref.fa")
            with open(fasta_out, "w") as fasta_file:
                records = SeqIO.parse(reference, "genbank")
                for record in records:
                    SeqIO.write(record, fasta_file, "fasta")
            
            logging.info(f"Converted GenBank to FASTA: {fasta_out}")

        elif ref_fmt == "fasta":
            # Copy FASTA file
            fasta_out = os.path.join(ref_dir, "ref.fa")
            shutil.copy(reference, fasta_out)
            logging.info(f"Copied FASTA file to {fasta_out}")

        # Ensure `fasta_out` is symlinked to `genomes/`
        genomes_fasta_out = os.path.join(ref_dir, "genomes", "ref.fa")
        if not os.path.exists(genomes_fasta_out):
            # Create a relative symlink pointing to the ref.fa in the parent directory
            os.symlink(os.path.relpath(fasta_out, start=os.path.join(ref_dir, "genomes")), genomes_fasta_out)
            logging.info(f"Created relative symlink: {genomes_fasta_out} -> {fasta_out}")
        else:
            logging.info(f"Symlink already exists: {genomes_fasta_out}")

        # Run samtools faidx
        try:
            result = subprocess.run(["samtools", "faidx", fasta_out], capture_output=True, text=True, check=True)
        except subprocess.CalledProcessError as e:
            logging.error(f"❌ Error running samtools faidx: {e.stderr}")
            return None
            
        return fasta_out

    except Exception as e:
        logging.error(f"❌ Error in extract_fasta_and_gbk: {e}")
        return None
