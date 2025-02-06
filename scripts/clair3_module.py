import subprocess
import os
import shutil
import logging

def run_clair3(output_dir, ref, bam_output, threads=2, model_name="r1041_e82_400bps_sup_v500", sample="", model_path=None):
    """
    Run Clair3 for long-read alignment and variant calling, followed by variant normalization and filtering out LowQual variants.
    """
    try:
        # Determine absolute path to the 'models' directory relative to the script's location
        script_dir = os.path.dirname(os.path.abspath(__file__))
        base_dir = os.path.dirname(script_dir)
        default_model_path = os.path.join(base_dir, "models", model_name)
        parallel_path = os.path.join(base_dir, "binaries", "parallel")

        # Set model path if not provided or explicitly set to None
        model_path = model_path if model_path is not None else default_model_path

        # Use absolute path for output directory
        output_dir = os.path.abspath(output_dir)

        # Create the clair3 subdirectory within the output directory
        clair3_output_dir = os.path.join(output_dir, "clair3")
        os.makedirs(clair3_output_dir, exist_ok=True)

        # Clair3 command
        clair3_cmd = [
            "run_clair3.sh",
            f"--bam_fn={os.path.abspath(bam_output)}",
            f"--ref_fn={os.path.abspath(ref)}",
            f"--threads={threads}",
            "--platform=ont",
            f"--parallel={parallel_path}",
            f"--model_path={os.path.abspath(model_path)}",
            f"--output={clair3_output_dir}",
            f"--sample_name={sample}",
            "--include_all_ctgs",
            "--haploid_precise",
            "--no_phasing_for_fa",
            "--enable_long_indel"
        ]

        # Print the command for debugging
        logging.info(f"🔹 Running Clair3 with command: {' '.join(clair3_cmd)}")

        # Check if model path exists
        if not os.path.exists(model_path):
            logging.error(f"❌ ERROR: Model path does not exist: {model_path}")
            return None
        else:
            logging.info(f"✅ Clair3 model path exists: {model_path}")

        # Run the Clair3 command and capture stdout and stderr
        with subprocess.Popen(clair3_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True) as clair3_proc:
            # Print Clair3's stdout and stderr in real-time
            for line in clair3_proc.stdout:
                logging.info(f"Clair3: {line.strip()}")  # Capture and log stdout in real-time
            for line in clair3_proc.stderr:
                logging.error(f"Clair3: {line.strip()}")  # Capture and log stderr in real-time

            # Wait for Clair3 to complete
            clair3_return = clair3_proc.wait()

            # If the process didn't complete successfully, return None
            if clair3_return != 0:
                logging.error(f"❌ Clair3 failed with error code {clair3_return}. Check logs above.")
                return None

        # Remove 'tmp' directory within clair3_output_dir if it exists
        tmp_dir = os.path.join(clair3_output_dir, "tmp")
        if os.path.exists(tmp_dir):
            shutil.rmtree(tmp_dir)
            logging.info(f"✅ Removed 'tmp' directory from {clair3_output_dir}")

        # Move the result files to the top output directory and rename them
        merge_vcf = os.path.join(clair3_output_dir, "merge_output.vcf.gz")
        merge_vcf_tbi = os.path.join(clair3_output_dir, "merge_output.vcf.gz.tbi")
        
        if os.path.exists(merge_vcf) and os.path.exists(merge_vcf_tbi):
            # Define output VCF path
            vcf_out = os.path.join(output_dir, f"{sample}_output_raw.vcf.gz")
            vcf_tbi_out = os.path.join(output_dir, f"{sample}_output_raw.vcf.gz.tbi")
            
            # Move and rename the files
            shutil.move(merge_vcf, vcf_out)
            shutil.move(merge_vcf_tbi, vcf_tbi_out)
            
            logging.info(f"✅ Moved and renamed files to {output_dir}:")
            logging.info(f"{vcf_out}")
            logging.info(f"{vcf_tbi_out}")

            # Filter out LowQual variants using bcftools view
            vcf_out_filt = os.path.join(output_dir, f"{sample}_output.filt.vcf")
            bcftools_view_cmd = [
                "bcftools", "view", "-f", "PASS", vcf_out, "-O", "v", "-o", vcf_out_filt
            ]
            logging.info(f"🔹 Running bcftools view to filter LowQual variants: {' '.join(bcftools_view_cmd)}")
            subprocess.run(bcftools_view_cmd, check=True)

            # Normalize the VCF using bcftools norm (left-align, right-align, split multi-allelic variants)
            norm_vcf = os.path.join(output_dir, f"{sample}_output.filt.norm.vcf.gz")
            bcftools_norm_cmd = [
                "bcftools", "norm", "-a", "-m-both", "-f", os.path.abspath(ref), vcf_out_filt, "-O", "z", "-o", norm_vcf
            ]
            logging.info(f"🔹 Running bcftools norm to normalize variants: {' '.join(bcftools_norm_cmd)}")
            subprocess.run(bcftools_norm_cmd, check=True)
            logging.info(f"✅ Clair3 variant calling and VCF processing complete: {vcf_out}")

            # Return the filtered and normalized VCF file path
            return vcf_out_filt

        else:
            logging.error(f"❌ ERROR: Expected VCF files not found in {clair3_output_dir}")
            return None

    except subprocess.CalledProcessError as e:
        logging.error(f"❌ Error running Clair3: {e}")
        return None
    except Exception as e:
        logging.error(f"❌ Unexpected error: {e}")
        return None
