"""Generate a single BED file with the original and rule-recommended-extension definition of every
locus that gets genuinely extended by the motif-reanchoring rule, across three locus sets:
  - Round 2: the 12 known disease loci (+EP400) extended at the original threshold=0.333
    (data/reanchor_validation_catalog_2.json; EP400 not in that catalog, added directly from
    report_v2 section 2's verified narrower/wider coordinates).
  - Round 3: the TRExplorer v1 threshold-disagreement loci that genuinely extend at threshold=0.15
    (data/threshold_sweep_trexplorer_v1.tsv).
  - Round 4: the flank-repetitiveness-ranked top-20k loci that genuinely extend at threshold=0.15
    (data/top20k_flank_repetitive_extensions.tsv).

All source coordinates are already 0-based half-open (same convention used throughout this project's
ReferenceRegion/fasta.fetch calls), so no conversion is needed for BED output.

Usage:
    python3 generate_extended_definitions_bed.py
"""

import json
from pathlib import Path

import pandas as pd

DATA_DIR = Path(__file__).parent / "data"
OUTPUT_PATH = DATA_DIR / "extended_definitions.bed"

EP400_ORIGINAL = ("chr12", 132062548, 132062611)
EP400_NEW = ("chr12", 132062524, 132062611)


def round2_rows():
    catalog = json.load(open(DATA_DIR / "reanchor_validation_catalog_2.json"))
    by_locus = {}
    for entry in catalog:
        base_id, variant = entry["LocusId"].rsplit("_", 1)
        chrom, rng = entry["ReferenceRegion"].split(":")
        start, end = (int(x) for x in rng.split("-"))
        motif = entry["LocusStructure"].strip("()*")
        by_locus.setdefault(base_id, {})[variant] = (chrom, start, end, motif)

    rows = []
    for locus_id, variants in by_locus.items():
        chrom, ostart, oend, motif = variants["baseline"]
        _, nstart, nend, _ = variants["reanchored"]
        rows.append(("round2", locus_id, chrom, ostart, oend, nstart, nend, motif))

    chrom, ostart, oend = EP400_ORIGINAL
    _, nstart, nend = EP400_NEW
    rows.append(("round2", "EP400", chrom, ostart, oend, nstart, nend, "CAG"))
    return rows


def round3_rows():
    df = pd.read_csv(DATA_DIR / "threshold_sweep_trexplorer_v1.tsv", sep="\t")
    disagree = df[~df["agree"]]
    extended = disagree[(disagree["left_ext_015"] > 0) | (disagree["right_ext_015"] > 0)]
    rows = []
    for row in extended.itertuples():
        nstart = row.start_0based - row.left_ext_015
        nend = row.end + row.right_ext_015
        rows.append(("round3", row.locus_id, row.chrom, row.start_0based, row.end, nstart, nend, row.motif))
    return rows


def round4_rows():
    df = pd.read_csv(DATA_DIR / "top20k_flank_repetitive_extensions.tsv", sep="\t")
    extended = df[(df["left_ext_015"] > 0) | (df["right_ext_015"] > 0)]
    rows = []
    for row in extended.itertuples():
        nstart = row.start_0based - row.left_ext_015
        nend = row.end_1based + row.right_ext_015
        rows.append(("round4", row.LocusId, row.chrom, row.start_0based, row.end_1based, nstart, nend, row.motif))
    return rows


def main():
    all_rows = round2_rows() + round3_rows() + round4_rows()
    print(f"round2: {sum(1 for r in all_rows if r[0] == 'round2')} loci")
    print(f"round3: {sum(1 for r in all_rows if r[0] == 'round3')} loci")
    print(f"round4: {sum(1 for r in all_rows if r[0] == 'round4')} loci")
    print(f"total: {len(all_rows)} loci, {2 * len(all_rows)} BED entries")

    bed_rows = []
    for round_name, locus_id, chrom, ostart, oend, nstart, nend, motif in all_rows:
        bed_rows.append((chrom, ostart, oend, f"{locus_id}__{round_name}__original", 0, "."))
        bed_rows.append((chrom, nstart, nend, f"{locus_id}__{round_name}__new", 0, "."))

    bed_rows.sort(key=lambda r: (r[0], r[1], r[2]))

    with open(OUTPUT_PATH, "w") as f:
        for chrom, start, end, name, score, strand in bed_rows:
            f.write(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\n")

    print(f"\nWrote {len(bed_rows)} BED entries to {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
