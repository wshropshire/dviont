# dviONT Pipeline

## Overview

dviONT (DNA Variant Identification using ONT) is a variant calling pipeline designed specifically for Q20+ Oxford Nanopore Technologies sequencing data. This pipeline was heavily inspired by (1) Torsten Seemann's [short-read variant calling tool Snippy](https://github.com/tseemann/snippy) and (2) The best practices of ONT long-read variant calling as described in [this Michael Hall et al. eLife 2024 paper](https://doi.org/10.7554/eLife.98300). The pipeline facilitates the following:

1. Alignment of ONT sequencing reads that have been basecalled using dorado (tested with v.0.9.1) and the super accurate model (r1041_e82_400bps_sup_v500) to provided reference genome using Minimap2.
2. Variant calling using Clair3 parameters following best practices as described in aforementioned eLife journal article.
3. Optional annotation of variants with SnpEff (if a GenBank reference is provided).
4. Post-processing of variant call files (VCF), which is intended to be used for further downstream analyses.
5. A readable tab separated variant calling report heavily inspired by Torsten Seemann and the output from the short-read variant calling tool, Snippy.

---

## Features

- **Reference Support:** Handles both FASTA and GenBank formats.
- **Variant Calling:** Leverages Clair3 which has been shown to have precision/recall for variant calling using Q20+ ONT reads comparable or even better than some short-read variant calling workflows.
- **Annotation:** SnpEff integration for functional annotation of variants.
- **Core Genome Alignment:** Currently in production; can use the normalized vcf output to potentially create a core genome alignment that can be used to infer a phylogeny.

---

## Installation

**NOTE** I've only tested this with a Linux, RHEL 7.9 operating system. I am not sure how this will function in other OS environments. The easiest way to install is:

(1) Clone the GitHub repository and then (2) create a conda environment. **Make sure to modify the env name and prefix in the yaml file to point to correct conda environment**

```bash
git clone https://github.com/wshropshire/dviont
cd dviont
conda env create -f dviont.yaml
```

Ensure all required executables (e.g., Minimap2, Clair3, SnpEff) are in your PATH or specify their paths accordingly. Specifically, the dviont executable will be in the `bin` directory so either create a symlink to conda environment or modify your bashrc/zshrc script as necessary:

```bash
echo 'export PATH="$PATH:/path/to/dviont/bin"' >> ~/.bashrc
source ~/.bashrc
```

I don't know how to make conda packages, but open to help from others 😀

---

## Usage

The pipeline can be executed using the following command:

```bash
./dviont/bin/dviont \
    -o <output_directory> \
    -r <reference_genome> \
    -i <reads_file> \
    -t <threads> \
    -m <clair3_model_name> \
    -p <clair3_model_path> \
    -s <sample_name>
```

### Required Arguments

- `-o`, `--output_dir`: Path to the directory where results will be saved.
- `-r`, `--ref`: Reference genome file (FASTA or GenBank).
- `-i`, `--reads`: ONT Q20+ reads(FASTQ).

### Optional Arguments

- `-t`, `--threads`: Number of threads to use (default: 2).
- `-m`, `--model_name`: Clair3 model name (default: `r1041_e82_400bps_sup_v500`).
- `-p`, `--model_path`: Path to the Clair3 model (optional).
- `-s`, `--sample`: Prefix for output (default: `SAMPLE`).
- `-v`, `--version`: Display the version of the dviONT pipeline.

---

## Example

Note that I have included example fasta/GenBank reference files as well as ONT Q20+ reads in the `data` directory

```bash
python3 dviONT.py \
    -o ./data/dviont_results \
    -r ./data/test.gb \
    -i ./data/reads.fastq.gz \
    -t 4 \
    -m r1041_e82_400bps_sup_v500 \
    -s SAMPLE1
```

---

## Outputs

- **Aligned Sorted Reads:** `<output_dir>/<sample_name>_aligned_reads.bam`
- **Filtered Variants:** `<output_dir>/<sample_name>_output.filt.vcf`
- **Normalized Variants (For core genome alignment):**`<output_dir>/<sample_name>_output.filt.norm.vcf.gz`
- **Filtered Annotated Variants (if GenBank):** `<output_dir>/<sample_name>_annotated.vcf`
- **dviONT Variant Calling Report:** `<output_dir>/<sample_name>_dviont_report.tsv`

---

## Columns in the dviONT Report

Name | Description
-----|------------
CHROM | The header of the sequence the variant was detected
POS | Position in the sequence (1-based)
TYPE | The variant type: snp mnp indel
REF | The nucleotide(s) in the reference
ALT | The alternate nucleotide(s) supported by the reads
EVIDENCE | Fraction of support for ALT and REF respectively

If you supply a Genbank file as the `--ref` rather than a FASTA
file, dviONT provides further annotation.

Name | Description
-----|------------
ANNOT | Mutation type (i.e., intergenic vs. coding mutation)
IMPACT | snpEff predicted impact of mutation
GENE | The `/gene` tag of the feature (if it existed)
LOCUS_TAG | The `/locus_tag` of the feature (if it existed)
HGVS.c | Nucleotide mutation following HGVS nomenclature
HGVS.p | Protein mutation following HGVS nomenclautre 
PRODUCT | The `/product` tag of the feature (if it existed)
EFFECT | The `snpEff` annotated consequence of this variant (ANN tag in .vcf)

---

## License

[MIT License](LICENSE)

---

## Contributing

Feel free to contribute to the project by submitting issues or pull requests.

---

## Version

dviONT v0.0.1