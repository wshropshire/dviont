import os
import logging
import subprocess
import shutil
import pysam
from .merge_vcfs import merge_vcfs  # Import merge function

class Clair3Pipeline:
    def __init__(self, output_dir, ref, bam_output, sample, threads=2, model_name="r1041_e82_400bps_sup_v500", model_path=None):
        """Initializes the Clair3Pipeline class."""
        self.output_dir = os.path.abspath(output_dir)
        self.ref = os.path.abspath(ref)
        self.bam_output = os.path.abspath(bam_output)
        self.sample = sample
        self.threads = threads
        self.model_name = model_name

        # Locate dviont package root directory
        package_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # Moves up one level

        # Correct Model Path
        self.model_path = model_path or os.path.join(package_root, "models", model_name)

        self.clair3_output_dir = os.path.join(self.output_dir, "clair3")

    def run_pipeline(self):
        """Runs Clair3, merges VCFs, normalizes, and applies multiallelic filtering."""
        logging.info(f"🔹 Running Clair3 pipeline for sample: {self.sample}")
        logging.info(f"🔹 Using Clair3 Model Path: {self.model_path}")

        try:
            os.makedirs(self.clair3_output_dir, exist_ok=True)

            # Run Clair3
            clair3_cmd = [
                "run_clair3.sh",
                f"--bam_fn={self.bam_output}",
                f"--ref_fn={self.ref}",
                f"--threads={self.threads}",
                "--platform=ont",
                f"--model_path={self.model_path}",
                f"--output={self.clair3_output_dir}",
                f"--sample_name={self.sample}",
                "--include_all_ctgs",
                "--haploid_precise",
                "--no_phasing_for_fa",
                "--snp_min_af=0.02", # Testing minimum AF at lower proportion (incraese SNP sensitivity)
                "--enable_long_indel"
            ]

            logging.info(f"Executing Clair3 command: {' '.join(clair3_cmd)}")
            subprocess.run(clair3_cmd, check=True)

            # Remove 'tmp' directory from Clair3 output
            self.cleanup_tmp_folder()

            # Define VCF file paths
            pileup_vcf = os.path.join(self.clair3_output_dir, "pileup.vcf.gz")
            full_vcf = os.path.join(self.clair3_output_dir, "full_alignment.vcf.gz")
            merged_vcf = os.path.join(self.output_dir, f"{self.sample}_merged.vcf.gz")

            # Ensure both VCFs exist before merging
            if not os.path.exists(pileup_vcf) or not os.path.exists(full_vcf):
                logging.error("❌ Clair3 did not produce expected VCFs.")
                return None

            # Merge the VCFs
            logging.info("🔹 Merging pileup and full-alignment VCFs...")
            merge_vcfs(pileup_vcf, full_vcf, merged_vcf)

            # Normalize the merged VCF
            norm_vcf = os.path.join(self.output_dir, f"{self.sample}_merged.norm.vcf.gz")
            norm_cmd = [
                "bcftools", "norm", "-a", "-m-both", "-f", self.ref, merged_vcf, "-O", "z", "-o", norm_vcf
            ]
            logging.info(f"🔹 Running bcftools norm: {' '.join(norm_cmd)}")
            subprocess.run(norm_cmd, check=True)

            # Apply filtering for multiallelic sites and sort the VCF
            final_sorted_vcf = os.path.join(self.output_dir, f"{self.sample}_filtered.sorted.vcf.gz")
            self.filter_multiallelic_sites(norm_vcf, final_sorted_vcf)

            logging.info(f"✅ Clair3 pipeline complete. Final sorted VCF: {final_sorted_vcf}")
            return final_sorted_vcf

        except subprocess.CalledProcessError as e:
            logging.error(f"❌ Clair3 execution failed: {e}")
            return None
        except Exception as e:
            logging.error(f"❌ Unexpected error: {e}", exc_info=True)
            return None

    def cleanup_tmp_folder(self):
        """Removes the tmp directory inside Clair3 output if it exists."""
        tmp_dir = os.path.join(self.clair3_output_dir, "tmp")
        if os.path.exists(tmp_dir):
            shutil.rmtree(tmp_dir)

    def filter_multiallelic_sites(self, input_vcf, output_vcf):
        """
        Filters multiallelic sites to **keep only the variant with the highest allele frequency (AF)**
        and ensures the final VCF file is sorted.
        """
        logging.info(f"🔹 Filtering multiallelic sites in: {input_vcf}")

        # Parse VCF and store variants by (CHROM, POS)
        variants = {}
        with pysam.VariantFile(input_vcf) as vcf:
            for record in vcf:
                chrom, pos = record.chrom, record.pos
                sample = list(record.samples.keys())[0]  # Assuming single sample

                # Extract AF (Allele Frequency) safely
                allele_freqs = record.samples[sample].get("AF", [None])
                allele_freqs = [af for af in allele_freqs if af is not None]  # Remove None values

                allele_freq = max(allele_freqs) if allele_freqs else 0.0  # Default to 0.0 if AF missing

                key = (chrom, pos)

                # Store the variant at this position
                if key not in variants:
                    variants[key] = []

                variants[key].append({
                    "record": record,
                    "allele_freq": allele_freq
                })

        # Select only the variant with the highest AF per position
        filtered_variants = []
        for key, records in variants.items():
            best_call = max(records, key=lambda x: x["allele_freq"])
            filtered_variants.append(best_call["record"])

        # Write filtered VCF
        with pysam.VariantFile(input_vcf) as template_vcf, pysam.VariantFile(output_vcf, "w", header=template_vcf.header) as output:
            for variant in filtered_variants:
                output.write(variant)

        # Sort the filtered VCF in-place
        temp_sorted_vcf = output_vcf + ".tmp.gz"
        sort_cmd = [
            "bcftools", "sort", "-O", "z", "-o", temp_sorted_vcf, output_vcf
        ]

        try:
            subprocess.run(sort_cmd, check=True)
            shutil.move(temp_sorted_vcf, output_vcf) 
        except subprocess.CalledProcessError as e:
            logging.error(f"❌ Sorting VCF failed: {e}")

        return output_vcf


def run_clair3(output_dir, ref, bam_output, threads, model_name, sample, model_path=None):
    """
    Wrapper function to initialize and run Clair3Pipeline.

    Returns:
        str: Path to final sorted VCF file.
    """
    pipeline = Clair3Pipeline(output_dir, ref, bam_output, sample, threads, model_name, model_path)
    return pipeline.run_pipeline()
