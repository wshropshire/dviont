<img src="src/dviont/docs/DVIONT.png" width="50%">

# DNA Variant Identification using ONT (dviONT) Pipeline

## Overview

dviONT (DNA Variant Identification using ONT) is a bacteria variant calling pipeline designed specifically for Q20+ Oxford Nanopore Technologies sequencing data. This pipeline was heavily inspired by (1) Torsten Seemann's [short-read variant calling tool Snippy](https://github.com/tseemann/snippy) and (2) The best practices of ONT long-read variant calling as described in [this Michael Hall et al. eLife 2024 paper](https://doi.org/10.7554/eLife.98300). The pipeline facilitates the following:

1. Alignment of ONT sequencing reads to a provided reference genome using Minimap2 (default) or Winnowmap.
2. Variant calling using Clair3 parameters following best practices described in the aforementioned eLife article.
3. Optional annotation of variants with SnpEff when a GenBank reference is provided.
4. Post-processing of variant call files (VCFs) for reporting and downstream analyses.
5. A readable tab-separated variant calling report heavily inspired by Torsten Seemann's short-read variant calling tool Snippy.
6. Cohort analysis of multiple ONT samples against a common reference, including generation of a combined SNP alignment and pairwise SNP distance matrix.

**dviONT has been developed and tested extensively using *Escherichia coli* datasets. Performance on other bacterial species has not been evaluated as extensively and should be validated by the user.**

---

## Features

- **Reference Support:** Handles both FASTA and GenBank formats.
- **Read Alignment:** Supports Minimap2 (default) and Winnowmap.
- **Variant Calling:** Uses Clair3 for ONT variant calling.
- **Final Variant Callset:** By default dviONT merges Clair3's pileup and full-alignment VCFs, normalizes the calls with `bcftools norm`, and resolves multiallelic sites. This retains high-confidence pileup-only calls (typically INDELs) that Clair3's native merge drops when the full-alignment model does not confirm them. Pass `--final-vcf clair3` to use Clair3's native `clair3/merge_output.vcf.gz` instead.
- **Annotation:** SnpEff integration for functional annotation of variants when a GenBank reference is supplied.
- **Consensus FASTA:** Generates a consensus FASTA file based on Clair3 variant calls, which can be used for genome assembly QC purposes.
- **Cohort SNP Alignment:** `dviont cohort` processes multiple samples against a shared reference and generates a combined SNP alignment.
- **SNP Distance Matrix:** Cohort mode uses `snp-dists` to calculate pairwise SNP distances from the cohort SNP alignment.

---

## Installation

> [!WARNING]
> dviONT has been developed and validated primarily on Linux/HPC with Clair3 v1.2.0 (TensorFlow). Clair3 v1.2.0 is not built for macOS, so the macOS environment uses Clair3 v2.0.3 (PyTorch) instead. the MacOS VERSION IS IN DEVELOPMENT. Variant calls may differ slightly between Clair3 versions, so validate macOS results against the included test data before relying on them. Clair3 does not support Intel Macs.

Clone the repository:

```bash
git clone https://github.com/wshropshire/dviont
cd dviont
```

### Linux / HPC

```bash
# Create the dviONT environment (Clair3 v1.x, Python 3.10)
mamba env create -f ./src/dviont/build/dviont_env.yaml   # or: conda env create -f ...
conda activate dviont_env

# Install dviONT
pip install .

# Download the Clair3 models used by dviONT (run once, on a node with internet access)
download_clair3_models
```

### macOS (Apple Silicon)

```bash
# Confirm your conda/mamba install is native Apple Silicon; this must report osx-arm64
conda info | grep platform

# Create the dviONT environment (Clair3 v2.x, Python 3.11)
mamba env create -f ./src/dviont/build/dviont_env_macOS.yaml
conda activate dviont_env

# Install dviONT
pip install .

# Clair3 v2 ships its PyTorch models inside the environment; dviONT finds them automatically
ls "$CONDA_PREFIX/bin/models"
```

> [!NOTE]
> You do not need `download_clair3_models` on macOS. It fetches Clair3 v1 (TensorFlow) models, which Clair3 v2 cannot load, so it exits without downloading when Clair3 v2 is installed.

## Clair3 models

dviONT finds the model named by `-m/--model-name` automatically. It looks in the `models` sub-directory of the installed dviONT package and in the models bundled with Clair3 (`$CONDA_PREFIX/bin/models`), and only uses a model whose format matches the installed Clair3 version (TensorFlow for Clair3 v1, PyTorch for Clair3 v2). Use `-p/--model-path` to point at a specific model directory instead.

- **Linux/HPC:** `download_clair3_models` saves models into the `models` sub-directory of the installed dviONT package. The default `r1041_e82_400bps_sup_v430_bacteria_finetuned` model is not bundled with Clair3 v1, so run it once after installing. Without arguments it downloads `r1041_e82_400bps_sup_v430_bacteria_finetuned`, `r1041_e82_400bps_sup_v500`, `r1041_e82_400bps_sup_v420` and `r941_prom_sup_g5014`. To download specific models or use another location: `download_clair3_models [--output_dir DIR] [model_name ...]` (then pass `-p DIR/<model_name>`).
- **macOS:** all models listed by `ls "$CONDA_PREFIX/bin/models"` work out of the box; nothing to download.

Choose the model that matches the Dorado basecalling model used for your ONT data.

---

## Usage

The `call` command runs the standard dviONT single-isolate workflow:

```bash
dviont call \
    -o <output_directory> \
    -r <reference_genome> \
    -i <reads_file> \
    -t <threads> \
    -m <clair3_model_name> \
    -p <clair3_model_path> \
    -s <sample_name> \
    --preset <alignment_preset> \
    --aligner <minimap2_or_winnowmap>
```

### Required Arguments

- `-o`, `--output-dir`: Path to the directory where results will be saved.
- `-r`, `--ref`: Reference genome file (FASTA or GenBank).
- `-i`, `--reads`: ONT Q20+ reads(FASTQ).

### Optional Arguments

- `-t`, `--threads`: Number of threads to use (default: 2).
- `-m`, `--model-name`: Clair3 model name (default: `r1041_e82_400bps_sup_v430_bacteria_finetuned`).
- `-p`, `--model-path`: Path to the Clair3 model (optional).
- `-s`, `--sample`: Prefix for output (default: `SAMPLE`).
- `--preset`: Minimap2 alignment preset (default: ont-q20).
    - `ont-legacy`: map-ont (ONT R9.x Guppy HAC)
    - `ont-q20`: lr:hq (ONT R10 Q20+ / Dorado SUP or duplex)
    - `pb-clr`: map-pb (PacBio CLR)
    - `pb-hifi`: map-hifi (PacBio HiFi/CCS)
    - `asm`: asm5 (assembly-to-assembly alignment)
- `--aligner`: Read aligner: `minimap2` (default) or `winnowmap`. Winnowmap uses a weighted repeat k-mer list generated by meryl.
- `--final-vcf`: Callset used for the final VCF, report, and consensus FASTA (default: `dviont`).
    - `dviont`: merge Clair3 pileup and full-alignment VCFs (a PASS call in either model is retained; when both PASS, the higher-QUAL record wins), normalize with `bcftools norm`, and keep the highest-AF allele at multiallelic sites.
    - `clair3`: use Clair3's native `clair3/merge_output.vcf.gz` as-is.
- `-v`, `--version`: Display the version of the dviONT pipeline.

---

## Example

Example GenBank/FASTA references and ONT Q20+ reads are included in `src/dviont/data`. From the repository root:

```bash
dviont call \
    -o ./dviont_test_results \
    -r ./src/dviont/data/test.gb \
    -i ./src/dviont/data/test_sup_v500_dorado091.fastq.gz \
    -t 4 \
    -m r1041_e82_400bps_sup_v500 \
    -s SAMPLE1 \
    --preset ont-q20
```

---

## Cohort mode

`dviont cohort` runs the ordinary dviONT workflow for multiple ONT read sets against one reference and preserves each sample's outputs. Cohort aggregation uses the final VCF returned by each completed call. It then uses `bcftools consensus` to build a full reference-length pseudoalignment before calculating the pairwise SNP distance matrix.

Provide a tab-separated samples file with one sample and reads path per line:

```text
SAMPLE1<TAB>/path/to/SAMPLE1.fastq.gz
SAMPLE2<TAB>/path/to/SAMPLE2.fastq.gz
```

Run cohort mode with:

```bash
dviont cohort \
    --ref ref.fasta \
    --reads-list samples.tsv \
    --out cohort_out \
    --threads 16 \
    --model-name r1041_e82_400bps_sup_v430_bacteria_finetuned \
    --model-path /path/to/clair3/models/r1041_e82_400bps_sup_v430_bacteria_finetuned \
    --preset ont-q20 \
    --aligner minimap2
```

The cohort output is organized as follows:

```text
cohort_out/
├── calls/                         # Standard per-sample dviONT output directories
│   ├── SAMPLE1/
│   └── SAMPLE2/
├── alignments/
│   ├── consensus_snps/            # Full reference-length sample consensuses
│   │   ├── SAMPLE1.fasta
│   │   └── SAMPLE2.fasta
│   └── cohort.snp_alignment.fasta
├── cohort_vcfs/
│   ├── cohort_merged.vcf.gz
│   ├── cohort_merged.vcf.gz.csi
│   ├── cohort_merged.norm.vcf.gz
│   ├── cohort_merged.norm.vcf.gz.csi
│   ├── cohort_merged.snps.vcf.gz
│   └── cohort_merged.snps.vcf.gz.csi
└── distances/
    └── cohort.snp_distance_matrix.tsv
```

All consensus sequences are checked against the processed cohort reference length. The final SNP alignment retains reference bases at nonvariant positions and genomic spacing between variants.

---

## Outputs

- **Aligned Sorted Reads:** `<output_dir>/<sample>_aln_sort.bam`
- **Clair3 Raw Calls:** `<output_dir>/clair3/pileup.vcf.gz`, `<output_dir>/clair3/full_alignment.vcf.gz`, `<output_dir>/clair3/merge_output.vcf.gz`
- **Merged Variants (`--final-vcf dviont`, default):** `<output_dir>/<sample>_merged.vcf.gz`
- **Normalized Variants (`--final-vcf dviont`, default):** `<output_dir>/<sample>_merged.norm.vcf.gz`
- **Final Filtered Variants (`--final-vcf dviont`, default):** `<output_dir>/<sample>_filtered.sorted.vcf.gz`
- **Final Variants (`--final-vcf clair3`):** `<output_dir>/clair3/merge_output.vcf.gz`
- **Consensus FASTA:** `<output_dir>/<sample>_consensus.fasta`
- **Annotated Variants (GenBank references only):** `<output_dir>/<sample>_annotated.vcf`
- **dviONT Variant Calling Report:** `<output_dir>/<sample>_dviont_report.tsv`

The final dviONT VCF used for reporting, the consensus FASTA, and downstream analysis is `<sample>_filtered.sorted.vcf.gz` by default, or the native `clair3/merge_output.vcf.gz` when `--final-vcf clair3` is given. Cohort mode uses the same final callset returned by the corresponding single-sample workflow.

---

## Columns in the dviONT Report

| Name | Description |
| --- | --- |
| CHROM | Reference sequence or contig containing the variant |
| POS | Variant position (1-based) |
| TYPE | Variant type: SNP, MNP, or INDEL |
| REF | Reference allele |
| ALT | Alternate allele |
| EVIDENCE | Read-depth evidence for the alternate and reference alleles derived from the VCF DP and AD fields |

If a GenBank file is supplied with `--ref`, dviONT also reports functional annotation from SnpEff:

| Name | Description |
| --- | --- |
| ANNOT | Predicted variant annotation |
| IMPACT | SnpEff predicted impact |
| GENE | Gene name, when available |
| LOCUS_TAG | GenBank locus tag |
| HGVS.c | Coding DNA-level HGVS annotation |
| HGVS.p | Protein-level HGVS annotation |
| PRODUCT_ID | Protein identifier from the GenBank annotation, when available |
| PRODUCT | Product description from the GenBank annotation, when available |
---

## License

[MIT License](LICENSE.txt)

---

## Contributing

Feel free to contribute to the project by submitting issues or pull requests.

---

## Version

dviONT v0.6.0
