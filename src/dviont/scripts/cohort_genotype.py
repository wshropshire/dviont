#!/usr/bin/env python3
"""Cohort consensus by re-genotyping candidate SNP sites from each isolate's reads.

Discovery is left to the variant caller (Clair3/dviONT): every SNP called in ANY isolate
becomes a candidate site. Each isolate is then genotyped at every candidate site directly
from its own BAM, and every other position is REF only where that isolate's reads cover it:

  callable position      >= --min-depth reads with MAPQ >= --min-mapq (and base quality
                         >= --min-baseq); optionally <= --max-depth-factor x median depth
  candidate site         base = the majority A/C/G/T base if it has >= --min-af of the reads,
                         otherwise N
  other callable sites   REF
  not callable           N

So an absent call is never silently turned into REF: a site is REF only when the reads say so.
N never creates a difference; it can only remove one, and every N is counted in the stats.

Used by the st_snp_pipeline (03_build_aln.py) and as a drop-in for the dviONT cohort step
(`dviont cohort`, which previously used `bcftools consensus` with no mask).

CLI:
  cohort_genotype.py --ref ref.fa --samples samples.tsv --out-prefix out/cohort \
      [--sites-vcf cohort_merged.snps.vcf.gz] [--min-depth 10 --min-af 0.8 --min-mapq 5]
      [--sensitivity-af 0.7,0.9]
samples.tsv: sample<TAB>bam[<TAB>vcf]   (vcf = that isolate's own calls; used for candidate sites
              when --sites-vcf is not given, and for the per-isolate concordance stats)

Outputs:
  <prefix>.aln.fasta              reference-length alignment, one record per isolate
  <prefix>.genotype_stats.tsv     per isolate: callable bp, N counts, Clair3 calls confirmed /
                                  overturned to REF / masked, ALT recovered at sites it did not call
  <prefix>.sites.tsv              per candidate site: how many isolates ALT / REF / N
  <prefix>.sensitivity_af.tsv     pairwise distances at each --sensitivity-af threshold
                                  (unmasked, candidate sites only) -- with --sensitivity-af
"""
import argparse, gzip, itertools, logging, os, sys
from collections import OrderedDict

import numpy as np

BASES = np.frombuffer(b"ACGT", dtype=np.uint8)
N = ord("N")


def read_fasta(path):
    seqs, name, buf = OrderedDict(), None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf).upper()
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf).upper()
    return seqs


def _open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def vcf_snps(path, keep_filters=("PASS", "."), sample_col=0):
    """{(contig, pos): alt} for SNPs whose genotype carries an ALT allele (any sample column if
    sample_col is None; MNPs of equal length are split into single-base SNPs)."""
    out = {}
    with _open(path) as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if keep_filters and f[6] not in keep_filters:
                continue
            ref, alts = f[3], f[4].split(",")
            if len(f) > 9:
                cols = range(9, len(f)) if sample_col is None else [9 + sample_col]
                fmt = f[8].split(":")
                carried = set()
                for c in cols:
                    gt = dict(zip(fmt, f[c].split(":"))).get("GT", ".")
                    carried |= {int(x) for x in gt.replace("|", "/").split("/") if x.isdigit() and x != "0"}
            else:
                carried = set(range(1, len(alts) + 1))
            for k in carried:
                alt = alts[k - 1]
                if alt in ("*", ".") or len(alt) != len(ref):
                    continue
                for off, (r, a) in enumerate(zip(ref, alt)):
                    if r != a and a in "ACGT":
                        out[(f[0], int(f[1]) + off)] = a
    return out


def base_counts(bam_path, contig, length, min_mapq, min_baseq):
    """4 x length array of A/C/G/T counts from reads passing MAPQ/flag filters."""
    import pysam

    def ok(read):
        return (not read.is_unmapped and not read.is_secondary and not read.is_qcfail
                and not read.is_duplicate and read.mapping_quality >= min_mapq)

    with pysam.AlignmentFile(bam_path) as bam:
        if contig not in bam.references:
            return np.zeros((4, length), dtype=np.uint32)
        cov = bam.count_coverage(contig, 0, length, quality_threshold=min_baseq, read_callback=ok)
    return np.vstack([np.frombuffer(c, dtype=np.uint64 if c.itemsize == 8 else np.uint32).astype(np.uint32)
                      for c in cov])


def genotype_sample(bam, ref, sites, a):
    """Return (sequence dict contig -> uint8 array, per-site (majority base idx, fraction, callable),
    stats dict)."""
    seqs, site_info = {}, {}
    st = {"callable_bp": 0, "ref_bp": 0, "n_lowdepth": 0, "n_highdepth": 0,
          "sites_alt": 0, "sites_ref": 0, "sites_N_mixed": 0, "sites_N_lowdepth": 0, "pileup_only_variants": 0}
    for contig, rs in ref.items():
        L = len(rs)
        r = np.frombuffer(rs.encode(), dtype=np.uint8).copy()
        cnt = base_counts(bam, contig, L, a.min_mapq, a.min_baseq)
        depth = cnt.sum(axis=0)
        callable_ = depth >= a.min_depth
        st["n_lowdepth"] += int((~callable_).sum())
        if a.max_depth_factor and a.max_depth_factor > 0:
            med = np.median(depth[depth > 0]) if (depth > 0).any() else 0
            hi = depth > a.max_depth_factor * med
            st["n_highdepth"] += int((hi & callable_).sum())
            callable_ &= ~hi
        maj = cnt.argmax(axis=0)
        frac = np.divide(cnt.max(axis=0), depth, out=np.zeros(L), where=depth > 0)
        seq = np.where(callable_, r, N).astype(np.uint8)
        # positions where the reads clearly disagree with REF but no isolate had a call (QC only)
        clear = callable_ & (frac >= a.min_af) & (BASES[maj] != r) & np.isin(r, BASES)
        pos = np.asarray(sites.get(contig, []), dtype=np.int64) - 1
        pos = pos[(pos >= 0) & (pos < L)]
        if len(pos):
            m_ok = callable_[pos] & (frac[pos] >= a.min_af)
            seq[pos] = np.where(m_ok, BASES[maj[pos]], N)
            clear[pos] = False
            site_info[contig] = (pos, maj[pos].astype(np.uint8), frac[pos].astype(np.float32), callable_[pos])
            st["sites_alt"] += int((m_ok & (BASES[maj[pos]] != r[pos])).sum())
            st["sites_ref"] += int((m_ok & (BASES[maj[pos]] == r[pos])).sum())
            st["sites_N_mixed"] += int((callable_[pos] & ~m_ok).sum())
            st["sites_N_lowdepth"] += int((~callable_[pos]).sum())
        st["pileup_only_variants"] += int(clear.sum())
        st["callable_bp"] += int(callable_.sum())
        seqs[contig] = seq
    return seqs, site_info, st


_W = None


def _work(item):
    s, bam = item
    if not os.path.exists(bam):
        return s, None, None, None
    ref, sites, a = _W
    seqs, info, st = genotype_sample(bam, ref, sites, a)
    return s, seqs, info, st


def pair_distances(names, site_tables, ref, thresholds):
    """Pairwise distances at candidate sites for each AF threshold (site_tables: name -> contig ->
    (pos, maj, frac, callable))."""
    rows = []
    contigs = list(ref)
    for t in thresholds:
        calls = {}
        for n in names:
            parts = []
            for c in contigs:
                if c not in site_tables[n]:
                    continue
                pos, maj, frac, cal = site_tables[n][c]
                parts.append(np.where(cal & (frac >= t), maj, 255).astype(np.uint8))
            calls[n] = np.concatenate(parts) if parts else np.zeros(0, np.uint8)
        for x, y in itertools.combinations(names, 2):
            a_, b_ = calls[x], calls[y]
            both = (a_ != 255) & (b_ != 255)
            rows.append((t, x, y, int((both & (a_ != b_)).sum()), int(both.sum())))
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--samples", required=True, help="sample<TAB>bam[<TAB>vcf]")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--sites-vcf", help="candidate SNP sites (e.g. dviONT cohort_merged.snps.vcf.gz); default = union of sample VCFs")
    ap.add_argument("--keep-filters", default="PASS,.", help="FILTER values accepted when reading candidate sites from VCFs")
    ap.add_argument("--min-depth", type=int, default=10)
    ap.add_argument("--min-af", type=float, default=0.8, help="majority-base fraction required to genotype a site")
    ap.add_argument("--min-mapq", type=int, default=5)
    ap.add_argument("--min-baseq", type=int, default=0)
    ap.add_argument("--max-depth-factor", type=float, default=0, help="N where depth > factor x median (collapsed repeats); 0 = off")
    ap.add_argument("--sensitivity-af", default="", help="comma-separated extra AF thresholds for the sensitivity table, e.g. 0.7,0.9")
    ap.add_argument("--threads", type=int, default=1, help="isolates genotyped in parallel")
    ap.add_argument("--extra-record", default="", help="name for an extra record holding the unmodified reference (e.g. REF_MB10992)")
    a = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    ref = read_fasta(a.ref)
    samples = []
    for line in open(a.samples):
        f = line.rstrip("\n").split("\t")
        if not f[0] or f[0].startswith("#"):
            continue
        samples.append((f[0], f[1], f[2] if len(f) > 2 and f[2] else None))
    keep = tuple(x for x in a.keep_filters.split(",") if x)

    own = {s: vcf_snps(v, keep) if v and os.path.exists(v) else {} for s, _, v in samples}
    if a.sites_vcf:
        cand = vcf_snps(a.sites_vcf, keep, sample_col=None)
    else:
        cand = {}
        for d in own.values():
            cand.update(d)
    sites = {}
    for (c, p) in cand:
        sites.setdefault(c, []).append(p)
    for c in sites:
        sites[c] = sorted(set(sites[c]))
    n_sites = sum(len(v) for v in sites.values())
    logging.info("%d isolates, %d candidate SNP sites", len(samples), n_sites)

    os.makedirs(os.path.dirname(os.path.abspath(a.out_prefix)), exist_ok=True)
    stats, site_tables, site_counts = [], {}, {}
    total = sum(len(s) for s in ref.values())
    with open(f"{a.out_prefix}.aln.fasta", "w") as aln:
        if a.extra_record:
            aln.write(f">{a.extra_record}\n" + "".join(ref.values()) + "\n")
        global _W
        _W = (ref, sites, a)
        todo = [(s, bam) for s, bam, _ in samples]
        if a.threads > 1:
            import multiprocessing as mp
            pool = mp.get_context("fork").Pool(a.threads)
            results = pool.imap(_work, todo)
        else:
            pool, results = None, map(_work, todo)
        for (s, bam, vcf), (s2, seqs, info, st) in zip(samples, results):
            if seqs is None:
                logging.warning("%s: BAM missing (%s) -- written as all N", s, bam)
                aln.write(f">{s}\n" + "N" * total + "\n")
                stats.append({"sample": s, "note": "BAM missing"})
                continue
            site_tables[s] = info
            full = b"".join(seqs[c].tobytes() for c in ref)
            aln.write(f">{s}\n" + full.decode() + "\n")
            # concordance with this isolate's own Clair3 calls
            o = own.get(s, {})
            conf = rej = msk = 0
            for (c, p), alt in o.items():
                if c not in seqs or p > len(seqs[c]):
                    continue
                b = chr(seqs[c][p - 1])
                if b == alt:
                    conf += 1
                elif b == "N":
                    msk += 1
                elif b == ref[c][p - 1]:
                    rej += 1
            recovered = 0
            for c, (pos, maj, frac, cal) in info.items():
                called = {p for (cc, p) in o if cc == c}
                ok = cal & (frac >= a.min_af) & (BASES[maj] != np.frombuffer(ref[c].encode(), np.uint8)[pos])
                recovered += int(sum(1 for p in (pos[ok] + 1) if int(p) not in called))
                for p, alt_ok, m in zip(pos + 1, ok, cal & (frac >= a.min_af)):
                    k = (c, int(p))
                    sc = site_counts.setdefault(k, [0, 0, 0])
                    if alt_ok:
                        sc[0] += 1
                    elif m:
                        sc[1] += 1
                    else:
                        sc[2] += 1
            n_bases = sum(len(x) for x in seqs.values())
            nN = sum(int((x == N).sum()) for x in seqs.values())
            st.update({"sample": s, "reference_bp": n_bases, "frac_N": round(nN / n_bases, 4),
                       "clair3_snps": len(o), "clair3_confirmed": conf, "clair3_overturned_to_ref": rej,
                       "clair3_masked_N": msk, "alt_recovered_not_called": recovered})
            stats.append(st)
            logging.info("%s: callable %.1f%%, ALT %d, Clair3 SNPs %d (confirmed %d, REF %d, N %d), recovered %d",
                         s, 100 * st["callable_bp"] / n_bases, st["sites_alt"], len(o), conf, rej, msk, recovered)

    if pool is not None:
        pool.close(); pool.join()
    cols = ["sample", "reference_bp", "callable_bp", "frac_N", "n_lowdepth", "n_highdepth", "sites_alt", "sites_ref",
            "sites_N_mixed", "sites_N_lowdepth", "clair3_snps", "clair3_confirmed", "clair3_overturned_to_ref",
            "clair3_masked_N", "alt_recovered_not_called", "pileup_only_variants", "note"]
    with open(f"{a.out_prefix}.genotype_stats.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for st in stats:
            fh.write("\t".join(str(st.get(k, "")) for k in cols) + "\n")
    with open(f"{a.out_prefix}.sites.tsv", "w") as fh:
        fh.write("contig\tpos\tref\tn_alt\tn_ref\tn_N\n")
        for c in ref:
            for p in sites.get(c, []):
                v = site_counts.get((c, p), [0, 0, 0])
                fh.write(f"{c}\t{p}\t{ref[c][p - 1]}\t{v[0]}\t{v[1]}\t{v[2]}\n")
    if a.sensitivity_af:
        th = sorted({a.min_af} | {float(x) for x in a.sensitivity_af.split(",") if x})
        names = [s for s, _, _ in samples if s in site_tables]
        with open(f"{a.out_prefix}.sensitivity_af.tsv", "w") as fh:
            fh.write("min_af\tsample1\tsample2\tsnps_at_candidate_sites\tsites_called_in_both\n")
            for r in pair_distances(names, site_tables, ref, th):
                fh.write("\t".join(map(str, r)) + "\n")
    logging.info("wrote %s.aln.fasta (+ genotype_stats, sites%s)", a.out_prefix, ", sensitivity_af" if a.sensitivity_af else "")
    return 0


if __name__ == "__main__":
    sys.exit(main())
