"""How many of the first 50,000 loci in Andrea's flank-metrics catalog would get extended by the
motif-reanchoring rule (motif_reanchor_extension.py) at gap-mismatch thresholds 0.15, 0.20, and 0.33
(the original spec)? Andrea's catalog isn't sorted by any obvious metric (checked by eye -- chrom
order is shuffled), so "top 50k" here means the first 50k rows as given, i.e. a representative slice
of the 415k-locus 2-6bp-motif catalog.

Usage:
    python3 scan_top50k_threshold_rates.py
"""

import sys
from pathlib import Path

import pandas as pd
import pysam

sys.path.insert(0, str(Path(__file__).parent))
from motif_reanchor_extension import compute_candidate_definitions

CATALOG_PATH = Path(__file__).parent.parent / "TR_catalog_2_6bp_with_flank_metrics.tsv.gz"
FASTA_PATH = "/Users/weisburd/hg38.fa"
N_LOCI = 50_000
THRESHOLDS = [0.15, 0.20, 0.33]
LABELS = {0.15: "015", 0.20: "020", 0.33: "033"}


def main():
    df = pd.read_csv(
        CATALOG_PATH,
        sep="\t",
        nrows=N_LOCI,
        usecols=["start_0based", "end_1based", "ReferenceRegion", "LocusId", "ReferenceMotif"],
    )
    print(f"loaded {len(df):,} loci from {CATALOG_PATH.name}")
    df["chrom"] = df["ReferenceRegion"].str.split(":").str[0]

    fasta = pysam.FastaFile(FASTA_PATH)

    counts = {t: 0 for t in THRESHOLDS}
    both_sides_counts = {t: 0 for t in THRESHOLDS}
    total_extension_bp = {t: 0 for t in THRESHOLDS}
    rows = []

    for i, row in enumerate(df.itertuples(index=False)):
        if i % 10_000 == 0:
            print(f"  {i:,} / {len(df):,}")
        per_threshold = {}
        for t in THRESHOLDS:
            result = compute_candidate_definitions(
                row.chrom, row.start_0based, row.end_1based, row.ReferenceMotif, fasta, threshold=t
            )
            left, right = result["left_extension"], result["right_extension"]
            per_threshold[t] = (left, right)
            if left > 0 or right > 0:
                counts[t] += 1
                total_extension_bp[t] += left + right
            if left > 0 and right > 0:
                both_sides_counts[t] += 1
        rows.append(
            {
                "LocusId": row.LocusId,
                "chrom": row.chrom,
                "start_0based": row.start_0based,
                "end_1based": row.end_1based,
                "motif": row.ReferenceMotif,
                **{f"left_ext_{LABELS[t]}": per_threshold[t][0] for t in THRESHOLDS},
                **{f"right_ext_{LABELS[t]}": per_threshold[t][1] for t in THRESHOLDS},
            }
        )

    out_df = pd.DataFrame(rows)
    out_path = Path(__file__).parent / "data" / "top50k_threshold_scan.tsv"
    out_df.to_csv(out_path, sep="\t", index=False)
    print(f"\nWrote {len(out_df):,} rows to {out_path}")

    print(f"\n=== Extension rate across first {N_LOCI:,} catalog loci ===")
    for t in THRESHOLDS:
        pct = 100 * counts[t] / len(df)
        both_pct = 100 * both_sides_counts[t] / len(df)
        mean_bp = total_extension_bp[t] / len(df)
        print(
            f"threshold={t:.2f}: {counts[t]:,} / {len(df):,} loci extended ({pct:.2f}%), "
            f"{both_sides_counts[t]:,} extended on both sides ({both_pct:.2f}%), "
            f"mean extension = {mean_bp:.2f}bp/locus, total = {total_extension_bp[t]:,}bp"
        )

    print("\n=== Overlap: how much of 0.15's extended set is a subset of 0.20's, and 0.20's of 0.33's ===")
    extended_015 = {r["LocusId"] for r in rows if r["left_ext_015"] > 0 or r["right_ext_015"] > 0}
    extended_020 = {r["LocusId"] for r in rows if r["left_ext_020"] > 0 or r["right_ext_020"] > 0}
    extended_033 = {r["LocusId"] for r in rows if r["left_ext_033"] > 0 or r["right_ext_033"] > 0}
    print(f"0.15 subset of 0.20: {extended_015 <= extended_020} ({len(extended_015 - extended_020)} loci in 0.15 but not 0.20)")
    print(f"0.20 subset of 0.33: {extended_020 <= extended_033} ({len(extended_020 - extended_033)} loci in 0.20 but not 0.33)")
    print(f"loci extended at 0.15 AND NOT at 0.33 (should be 0, since 0.15 is stricter): {len(extended_015 - extended_033)}")


if __name__ == "__main__":
    main()
