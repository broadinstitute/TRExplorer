"""Sort Andrea's flank-metrics catalog by how repetitive the flanks look (flank_motif_similarity_left/
right -- how closely the flank matches a perfect in-frame repeat of the locus's own annotated motif),
most-repetitive-flank first. These are exactly the loci most likely to be under-extended: a highly
repetitive flank is direct evidence the boundary stopped short of real repeat sequence.

Ranking key: max(flank_motif_similarity_left, flank_motif_similarity_right) per locus -- a locus
qualifies if EITHER side looks repetitive, since the boundary-extension rule evaluates each side
independently anyway.

Usage:
    python3 sort_catalog_by_flank_repetitiveness.py
"""

from pathlib import Path

import pandas as pd

CATALOG_PATH = Path(__file__).parent.parent / "TR_catalog_2_6bp_with_flank_metrics.tsv.gz"
OUTPUT_PATH = Path(__file__).parent / "data" / "TR_catalog_2_6bp_sorted_by_flank_repetitiveness.tsv.gz"
TOP_N_PATH = Path(__file__).parent / "data" / "top20k_by_flank_repetitiveness.tsv.gz"
TOP_N = 20_000


def main():
    df = pd.read_csv(CATALOG_PATH, sep="\t")
    print(f"loaded {len(df):,} loci from {CATALOG_PATH.name}")

    df["max_flank_motif_similarity"] = df[["flank_motif_similarity_left", "flank_motif_similarity_right"]].max(axis=1)
    df = df.sort_values("max_flank_motif_similarity", ascending=False)

    cols = ["max_flank_motif_similarity"] + [c for c in df.columns if c != "max_flank_motif_similarity"]
    df = df[cols]

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUTPUT_PATH, sep="\t", index=False, compression="gzip")
    print(f"wrote {len(df):,} rows to {OUTPUT_PATH}")

    df.head(TOP_N).to_csv(TOP_N_PATH, sep="\t", index=False, compression="gzip")
    print(f"wrote top {TOP_N:,} rows to {TOP_N_PATH} (input for scan_top20k_flank_repetitive_extensions.py)")

    preview_cols = [
        "LocusId", "ReferenceMotif", "MotifSize", "NumRepeatsInReference", "ReferenceRepeatPurity",
        "flank_motif_similarity_left", "flank_motif_similarity_right", "max_flank_motif_similarity",
    ]
    print(f"\n=== Top 30 loci by flank repetitiveness ===")
    print(df[preview_cols].head(30).to_string(index=False))

    print(f"\n=== max_flank_motif_similarity distribution ===")
    print(df["max_flank_motif_similarity"].describe())
    for q in [0.999, 0.99, 0.95, 0.90]:
        n_above = (df["max_flank_motif_similarity"] >= df["max_flank_motif_similarity"].quantile(q)).sum()
        print(f"  top {(1-q)*100:.1f}% ({n_above:,} loci) have max_flank_motif_similarity >= {df['max_flank_motif_similarity'].quantile(q):.4f}")


if __name__ == "__main__":
    main()
