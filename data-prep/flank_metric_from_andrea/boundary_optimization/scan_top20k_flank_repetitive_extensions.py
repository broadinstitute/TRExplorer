"""For the top 20,000 loci ranked by flank repetitiveness (data/top20k_by_flank_repetitiveness.tsv.gz,
from sort_catalog_by_flank_repetitiveness.py), compute the motif-reanchoring rule's recommended
extension at thresholds 0.15 and 0.20. This is the free/local pre-check before deciding how large a
Hail Batch validation catalog to build -- only loci that actually extend need real-read testing (a
narrow-vs-wide pOk comparison is meaningless for a locus the rule doesn't touch).

Usage:
    python3 scan_top20k_flank_repetitive_extensions.py
"""

from pathlib import Path
import sys

import pandas as pd
import pysam

sys.path.insert(0, str(Path(__file__).parent))
from motif_reanchor_extension import compute_candidate_definitions

CATALOG_PATH = Path(__file__).parent / "data" / "top20k_by_flank_repetitiveness.tsv.gz"
FASTA_PATH = "/Users/weisburd/hg38.fa"
THRESHOLDS = [0.15, 0.20]
LABELS = {0.15: "015", 0.20: "020"}


def main():
    df = pd.read_csv(
        CATALOG_PATH, sep="\t",
        usecols=["start_0based", "end_1based", "ReferenceRegion", "LocusId", "ReferenceMotif", "max_flank_motif_similarity"],
    )
    print(f"loaded {len(df):,} loci")
    df["chrom"] = df["ReferenceRegion"].str.split(":").str[0]

    fasta = pysam.FastaFile(FASTA_PATH)

    rows = []
    for i, row in enumerate(df.itertuples(index=False)):
        if i % 5_000 == 0:
            print(f"  {i:,} / {len(df):,}")
        rec = {
            "LocusId": row.LocusId, "chrom": row.chrom, "start_0based": row.start_0based,
            "end_1based": row.end_1based, "motif": row.ReferenceMotif,
            "max_flank_motif_similarity": row.max_flank_motif_similarity,
        }
        for t in THRESHOLDS:
            result = compute_candidate_definitions(row.chrom, row.start_0based, row.end_1based, row.ReferenceMotif, fasta, threshold=t)
            rec[f"left_ext_{LABELS[t]}"] = result["left_extension"]
            rec[f"right_ext_{LABELS[t]}"] = result["right_extension"]
        rows.append(rec)

    out_df = pd.DataFrame(rows)
    out_path = Path(__file__).parent / "data" / "top20k_flank_repetitive_extensions.tsv"
    out_df.to_csv(out_path, sep="\t", index=False)
    print(f"\nWrote {len(out_df):,} rows to {out_path}")

    for t in THRESHOLDS:
        lab = LABELS[t]
        extended = (out_df[f"left_ext_{lab}"] > 0) | (out_df[f"right_ext_{lab}"] > 0)
        print(f"threshold={t}: {extended.sum():,} / {len(out_df):,} extended ({100*extended.mean():.1f}%)")

    extended_015 = (out_df["left_ext_015"] > 0) | (out_df["right_ext_015"] > 0)
    extended_020 = (out_df["left_ext_020"] > 0) | (out_df["right_ext_020"] > 0)
    print(f"union (0.15 or 0.20): {(extended_015 | extended_020).sum():,}")
    print(f"0.15 subset of 0.20: {(extended_015 & ~extended_020).sum() == 0}")


if __name__ == "__main__":
    main()
