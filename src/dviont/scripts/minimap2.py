import subprocess
import os
import logging

def run_minimap2_alignment(fasta_file, reads, threads, output_dir, sample):
    """
    Run Minimap2 alignment and sort the output BAM file with sample name.

    Args:
        fasta_file (str): Path to the reference FASTA file (symlink in output_dir).
        reads (str): Path to the input reads (FASTQ file).
        threads (int): Number of threads to use.
        output_dir (str): Directory where the output BAM file will be saved.
        sample (str): Sample name to append to the output BAM file name.
    
    Returns:
        str: Path to the sorted BAM file if successful, None if failed.
    """
    bam_output = os.path.join(output_dir, f"{sample}_aln_sort.bam")

    minimap_cmd = ["minimap2", "-t", str(threads), "-ax", "lr:hq", "-O", "20,60", fasta_file, reads]
    samtools_sort_cmd = ["samtools", "sort", "-@", str(threads), "-o", bam_output]

    try:
        # Print the Minimap2 command in a single print statement
        logging.info(f"🔹 Running minimap2 with command:{' '.join(minimap_cmd)}")


        # Run minimap2 and capture both stdout and stderr
        with subprocess.Popen(minimap_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True) as minimap_proc, \
             subprocess.Popen(samtools_sort_cmd, stdin=minimap_proc.stdout, stderr=subprocess.PIPE, text=True) as sort_proc:

            # Print Minimap2's stderr in real-time
            for line in minimap_proc.stderr:
                logging.info(f"Minimap2: {line.strip()}")

            minimap_proc.stdout.close()  # Allow sort_proc to get EOF

            # Print the Samtools sort command
            logging.info(f"🔹 Running samtools sort with command:{' '.join(samtools_sort_cmd)}")

            # Wait for minimap2 and samtools to complete
            minimap_return = minimap_proc.wait()
            sort_stderr = sort_proc.communicate()[1]  # Capture Samtools stderr

            if minimap_return != 0:
                logging.error(f"❌ Minimap2 failed with error code {minimap_return}. Check logs above.")
                return None
            
            if sort_proc.returncode != 0:
                logging.error(f"❌ Samtools sorting failed: {sort_stderr}")
                return None

        # Index the sorted BAM file
        subprocess.run(["samtools", "index", bam_output], check=True)
        logging.info(f"✅ Minimap2 alignment and Samtools sort complete: {bam_output}")
        return bam_output

    except subprocess.CalledProcessError as e:
        logginf.error(f"❌ Error running Minimap2 or Samtools: {e.stderr}")
        return None
    except Exception as e:
        logging.error(f"❌ Unexpected error: {e}")
        return None
