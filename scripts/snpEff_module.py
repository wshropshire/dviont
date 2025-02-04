import os
import subprocess
import shutil
import logging
from extract_headers import extract_headers_from_fasta


def run_snpEff(output_dir, vcf_out, ref, sample=""):
    """
    Run snpEff to predict and annotate SNP changes.

    Args:
        output_dir (str): Directory where output will be saved.
        vcf_out (str): Path to the input VCF file (output from Clair3).
        ref (str): Path to the reference genome fasta file.
        sample (str, optional): Sample name used for output naming.

    Returns:
        str: Path to the snpEff annotated VCF file if successful, None otherwise.
    """

    try:
        # Determine paths for configuration and reference files
        script_dir = os.path.dirname(os.path.abspath(__file__))
        base_dir = os.path.dirname(script_dir)
        snpeff_config_src = os.path.join(base_dir, "etc", "snpeff.config")

        # Point to reference directory inside output_dir, ensuring it is absolute
        reference_dir = os.path.abspath(os.path.join(output_dir, "reference"))

        # Copy the snpeff_config file to the reference directory
        snpeff_config_work = os.path.join(reference_dir, "snpeff.config")
        shutil.copy(snpeff_config_src, snpeff_config_work)

        # Extract chromosome headers from the reference genome
        chromosomes = extract_headers_from_fasta(ref)

        # Append genome information to snpEff config (instead of overwriting)
        with open(snpeff_config_work, "a") as config_file:  # Use append mode
            config_file.write(f"\n# Custom genome configuration\n")
            config_file.write(f"{'ref'}.genome : dviONT Reference\n")
            config_file.write(f"{'ref'}.chromosome : {', '.join(chromosomes)}\n")
            for chrom in chromosomes:
                config_file.write(f"{'ref'}.{chrom}.codonTable : Bacterial_and_Plant_Plastid\n")

        # Define snpEff database name
        snpeff_db = "ref"

        # Run snpEff build to create an index file
        snpeff_build_cmd = [
            "snpEff", "build", "-c", snpeff_config_work, "-dataDir", reference_dir, "-genbank", snpeff_db
        ]
        logging.info(f"🔹 Running snpEff build: {' '.join(snpeff_build_cmd)}")
        subprocess.run(snpeff_build_cmd, check=True)

        # Define snpEff annotated VCF output file
        snpeff_annotated_vcf = os.path.join(output_dir, f"{sample}_annotated.vcf")

        # Run snpEff annotation with correct options
        snpeff_cmd = [
            "snpEff", "ann", "-noLog", "-noStats",
            "-no-downstream", "-no-upstream", "-no-utr",
            "-c", snpeff_config_work, "-dataDir", reference_dir,
            snpeff_db, vcf_out
        ]

        logging.info(f"🔹 Running snpEff annotation: {' '.join(snpeff_cmd)}")
        with open(snpeff_annotated_vcf, "w") as output_file:
            subprocess.run(snpeff_cmd, stdout=output_file, check=True)

        logging.info(f"✅ SnpEff annotation completed: {snpeff_annotated_vcf}")
        return snpeff_annotated_vcf

    except subprocess.CalledProcessError as e:
        logging.error(f"❌ Error running snpEff: {e}")
        return None
    except Exception as e:
        logging.error(f"❌ Unexpected error: {e}")
        return None
