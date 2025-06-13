import pysam
import logging

def parse_vcf(vcf_file):
    """ Parses a VCF file and returns a dictionary of variants. """
    variants = {}
    with pysam.VariantFile(vcf_file) as vcf:
        for record in vcf:
            alts = tuple(record.alts) if record.alts is not None else ("<NO_ALT>",)  # Placeholder for no ALT
            key = (record.chrom, record.pos, record.ref, alts)  # Unique key for variant

            sample = list(record.samples.keys())[0]  # Assuming single sample per VCF
            genotype = record.samples[sample]["GT"]

            variants[key] = {
                "record": record,
                "qual": record.qual,
                "filter": record.filter.keys(),
                "genotype": genotype,
            }
    return variants

def merge_vcfs(pileup_vcf, full_vcf, output_vcf):
    """ Merges two Clair3 VCFs (pileup and full-alignment) based on defined rules. """
    logging.info(f"🔹 Parsing pileup VCF: {pileup_vcf}")
    pileup_variants = parse_vcf(pileup_vcf)

    logging.info(f"🔹 Parsing full-alignment VCF: {full_vcf}")
    full_variants = parse_vcf(full_vcf)

    merged_variants = {}
    all_keys = set(pileup_variants.keys()).union(set(full_variants.keys()))

    for key in all_keys:
        pileup = pileup_variants.get(key)
        full = full_variants.get(key)

        if pileup and full:
            if "PASS" in pileup["filter"] and "PASS" in full["filter"]:
                # Rule 1: If one model calls 1/1 or 0/1, and the other calls 1/2 → Take 1/1 or 0/1
                if (pileup["genotype"] == (1, 2) and full["genotype"] in [(1, 1), (0, 1)]) or \
                   (full["genotype"] == (1, 2) and pileup["genotype"] in [(1, 1), (0, 1)]):
                    chosen_variant = pileup if pileup["genotype"] in [(1, 1), (0, 1)] else full
                    merged_variants[key] = chosen_variant["record"]
                else:
                    merged_variants[key] = pileup["record"] if pileup["qual"] > full["qual"] else full["record"]
            # Select pileup records that are exclusive to pileup model
            elif "PASS" in pileup["filter"]:
                merged_variants[key] = pileup["record"]
            # Select full alignment records that are exclusive to full alignment model
            elif "PASS" in full["filter"]:
                merged_variants[key] = full["record"]
        elif pileup and "PASS" in pileup["filter"]:
            merged_variants[key] = pileup["record"]
        elif full and "PASS" in full["filter"]:
            merged_variants[key] = full["record"]

    # Write merged VCF file
    with pysam.VariantFile(pileup_vcf) as template_vcf, pysam.VariantFile(output_vcf, "w", header=template_vcf.header) as output:
        for variant in merged_variants.values():
            output.write(variant)

    logging.info(f"✅ Merged VCF saved to: {output_vcf}")
    return output_vcf
