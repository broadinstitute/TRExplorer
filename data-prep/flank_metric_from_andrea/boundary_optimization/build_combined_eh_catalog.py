"""Build one combined ExpansionHunter catalog covering three locus populations:
  1. The 16,264 Andrea-catalog loci (exact TRExplorer v2 matches) that the gap-purity rule extends
     -- both __original and __extended definitions (from data/gap_purity_scan_v2_exact_matches.tsv).
  2. ALL 331,440 Andrea-catalog loci (exact TRExplorer v2 matches) that the rule does NOT extend --
     single __original definition each (same source file, left_ext==0 and right_ext==0).
  3. The known disease-associated loci (variant_catalog_without_offtargets.GRCh38.json, single-region
     only, excluding any IUPAC-degenerate-motif locus e.g. GCN/NGC) -- __original + __extended (deduped
     when the rule doesn't extend that locus; only EP400 among these actually extends).

LocusId convention throughout: "{source locus's own stable LocusId}__{original|extended}" -- keyed to
each locus's own original identity (not the resulting region's own coordinates), so the original and
extended variant of one locus always share a stable base id, and two different loci extending to the
same resulting region can never collide.

Usage:
    python3 build_combined_eh_catalog.py
"""

import json
import sys
from pathlib import Path

import pandas as pd
import pysam

sys.path.insert(0, str(Path(__file__).parent))
from gap_purity_extension import compute_candidate_definitions, MAX_OFFSET

DATA_DIR = Path(__file__).parent / "data"
ANDREA_SCAN_TSV = DATA_DIR / "gap_purity_scan_v2_exact_matches.tsv"
DISEASE_CATALOG_PATH = Path("~/code/str-analysis/str_analysis/variant_catalogs/variant_catalog_without_offtargets.GRCh38.json").expanduser()
FASTA_PATH = "/Users/weisburd/hg38.fa"
OUTPUT_PATH = DATA_DIR / "combined_eh_catalog.json"


def make_entry(base_locus_id, label, chrom, start, end, motif):
    return {
        "LocusId": f"{base_locus_id}__{label}",
        "ReferenceRegion": f"{chrom}:{start}-{end}",
        "LocusStructure": f"({motif})*",
        "VariantType": "Repeat",
    }


def andrea_entries():
    df = pd.read_csv(ANDREA_SCAN_TSV, sep="\t")
    n_extended = (df["left_ext"].gt(0) | df["right_ext"].gt(0)).sum()
    n_not_extended = len(df) - n_extended
    print(f"{len(df):,} total Andrea v2-exact-match loci ({n_extended:,} extended, {n_not_extended:,} not extended)")

    entries = []
    for row in df.itertuples():
        extended = row.left_ext > 0 or row.right_ext > 0
        entries.append(make_entry(row.LocusId, "original", row.chrom, row.start_0based, row.end_1based, row.motif))
        if extended:
            new_start, new_end = row.start_0based - row.left_ext, row.end_1based + row.right_ext
            entries.append(make_entry(row.LocusId, "extended", row.chrom, new_start, new_end, row.motif))
    print(f"  -> {len(entries):,} catalog entries")
    return entries


def disease_locus_entries():
    d = json.load(open(DISEASE_CATALOG_PATH))
    single = [e for e in d if not isinstance(e["ReferenceRegion"], list)]
    # exclude any locus whose repeat unit contains an IUPAC ambiguity code (e.g. GCN, NGC) -- this
    # project's earlier rounds found these need special handling; keep them out of this pipeline
    # entirely rather than risk passing an ambiguity-coded motif through incorrectly. Filter on the
    # LocusStructure-derived motif, not RepeatUnit -- RepeatUnit can be a simplified/literal form
    # (e.g. RFC1: RepeatUnit="AAAAG" but LocusStructure="(AARRG)*", the real degenerate repeat unit)
    # that would silently let a degenerate locus through the filter.
    non_degenerate = [e for e in single if set(e["LocusStructure"].strip("()*").upper()) <= set("ACGT")]
    print(f"{len(single)} single-region disease loci, {len(non_degenerate)} after excluding IUPAC-degenerate-motif loci")

    fasta = pysam.FastaFile(FASTA_PATH)
    entries = []
    for e in non_degenerate:
        chrom, rng = e["ReferenceRegion"].split(":")
        start, end = (int(x) for x in rng.split("-"))
        # use the source's own LocusStructure-derived motif, not RepeatUnit -- RepeatUnit can be a
        # simplified/literal form (e.g. RFC1: RepeatUnit="AAAAG" but LocusStructure="(AARRG)*", the
        # real IUPAC-degenerate repeat unit) and substituting RepeatUnit would silently narrow the motif.
        motif = e["LocusStructure"].strip("()*")

        result = compute_candidate_definitions(chrom, start, end, motif, fasta)
        entries.append(make_entry(e["LocusId"], "original", chrom, start, end, motif))
        if result["left_extension"] > 0 or result["right_extension"] > 0:
            new_start = start - result["left_extension"]
            new_end = end + result["right_extension"]
            entries.append(make_entry(e["LocusId"], "extended", chrom, new_start, new_end, motif))
            print(f"  {e['LocusId']}: extends -> {chrom}:{new_start}-{new_end}")
    print(f"  -> {len(entries):,} catalog entries")
    return entries


def main():
    andrea = andrea_entries()
    disease = disease_locus_entries()

    all_entries = andrea + disease
    seen_ids = set()
    deduped = []
    n_dupes = 0
    for e in all_entries:
        if e["LocusId"] in seen_ids:
            n_dupes += 1
            continue
        seen_ids.add(e["LocusId"])
        deduped.append(e)
    if n_dupes:
        print(f"{n_dupes} duplicate LocusIds dropped (disease locus coincides with an Andrea-catalog locus)")

    with open(OUTPUT_PATH, "w") as f:
        json.dump(deduped, f, indent=2)
    print(f"\nWrote {len(deduped):,} total catalog entries to {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
