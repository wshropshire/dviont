import csv
from Bio import SeqIO
import os
import gzip

class VCFProcessor:
    def __init__(self, vcf_file, ref_fmt, output_dir, sample, genbank_file=None):
        self.vcf_file = vcf_file
        self.ref_fmt = ref_fmt.lower()
        self.output_dir = output_dir
        self.sample = sample
        self.genbank_file = genbank_file
        self.genbank_dict = self.load_genbank(self.genbank_file) if self.genbank_file else {}

    def load_genbank(self, genbank_file):
        """Parse the GenBank file to create a dictionary linking LOCUS_TAG to product_id and product."""
        genbank_dict = {}
        for record in SeqIO.parse(genbank_file, "genbank"):
            for feature in record.features:
                if feature.type == "CDS" and 'locus_tag' in feature.qualifiers:
                    locus_tag = feature.qualifiers['locus_tag'][0]
                    product_id = feature.qualifiers.get('protein_id', [''])[0]
                    product = feature.qualifiers.get('product', [''])[0]
                    genbank_dict[locus_tag] = {'product_id': product_id, 'product': product}
        return genbank_dict

    def parse_vcf(self):
        """Parse the VCF file and generate the summary report."""
        header = ["CHROM", "POS", "TYPE", "REF", "ALT", "EVIDENCE"]
        if self.ref_fmt == "genbank":
            header += ["ANNOT", "IMPACT", "GENE", "LOCUS_TAG", "HGVS.c", "HGVS.p", "PRODUCT_ID", "PRODUCT"]

        rows = []

        open_func = gzip.open if self.vcf_file.endswith(".gz") else open
        with open_func(self.vcf_file, 'rt') as vcf:
            for line in vcf:
                if line.startswith("#"):
                    continue

                fields = line.strip().split("\t")
                # Minimum 8 columns for a valid VCF
                if len(fields) < 8:
                    continue

                chrom, pos, _, ref, alt, _, _, info, *rest = fields

                # Determine TYPE
                if len(ref) == 1 and len(alt) == 1:
                    variant_type = "SNP"
                elif len(ref) == len(alt):
                    variant_type = "MNP"
                else:
                    variant_type = "INDEL"

                # Extract EVIDENCE
                evidence = "NA"
                if rest:
                    sample_fields = rest[-1].split(":")  # FORMAT: GT:GQ:DP:AD:AF
                    if len(sample_fields) >= 4:
                        dp = sample_fields[2]
                        ad_values = sample_fields[3].split(",")  # e.g. "15,19"
                        if len(ad_values) >= 2:
                            ref_count, alt_count = ad_values[0], ad_values[1]
                            evidence = f"ALT:{dp}/{alt_count};REF:{dp}/{ref_count}"

                row = [chrom, pos, variant_type, ref, alt, evidence]

                # Only add GenBank annotations if requested
                if self.ref_fmt == "genbank" and "ANN=" in info:
                    info_fields = info.split("|")
                    annot = info_fields[1] if len(info_fields) > 1 else ""
                    impact = info_fields[2] if len(info_fields) > 2 else ""
                    gene = info_fields[3] if len(info_fields) > 3 else ""
                    locus_tag = info_fields[4] if len(info_fields) > 4 else ""

                    if annot == "intergenic_region":
                        gene = ""
                    elif gene == locus_tag:
                        gene = "hyp"

                    hgvs_c = info_fields[9] if len(info_fields) > 9 else ""
                    hgvs_p = info_fields[10] if len(info_fields) > 10 else ""

                    product_id, product = "", ""
                    if locus_tag in self.genbank_dict:
                        product_id = self.genbank_dict[locus_tag]["product_id"]
                        product = self.genbank_dict[locus_tag]["product"]

                    row += [annot, impact, gene, locus_tag, hgvs_c, hgvs_p, product_id, product]

                rows.append(row)

        # Write to TSV file
        self.write_tsv(header, rows)

    def write_tsv(self, header, rows):
        """Write the summary data to a TSV file."""
        os.makedirs(self.output_dir, exist_ok=True)
        output_path = os.path.join(self.output_dir, f"{self.sample}_dviont_report.tsv")

        with open(output_path, "w", newline="") as tsvfile:
            writer = csv.writer(tsvfile, delimiter='\t')
            writer.writerow(header)
            writer.writerows(rows)
