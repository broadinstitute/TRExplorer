"""Apply the gap-purity boundary-extension rule (gap_purity_extension.py) to every locus in Andrea's
catalog whose original definition exactly matches a TRExplorer v2 locus (chrom/start/end/motif) --
data/andrea_v2_exact_matches.sorted.txt, built earlier by intersecting data/andrea_locus_ids.sorted.txt
and data/trexplorer_v2_locus_ids.sorted.txt. Restricting to exact matches avoids re-litigating the
external-source boundary-convention mismatches found earlier (see conversation) -- these are loci
where the starting definition itself is trustworthy.

Usage:
    python3 scan_gap_purity_on_v2_exact_matches.py
"""

from pathlib import Path

import pandas as pd
import pysam

from gap_purity_extension import extend_by_gap_purity, MAX_OFFSET

DATA_DIR = Path(__file__).parent / "data"
CATALOG_PATH = Path(__file__).parent.parent / "TR_catalog_2_6bp_with_flank_metrics.tsv.gz"
MATCH_IDS_PATH = DATA_DIR / "andrea_v2_exact_matches.sorted.txt"
FASTA_PATH = "/Users/weisburd/hg38.fa"
OUTPUT_PATH = DATA_DIR / "gap_purity_scan_v2_exact_matches.tsv"


def main():
    match_ids = set(MATCH_IDS_PATH.read_text().split())
    print(f"{len(match_ids):,} loci to scan (exact TRExplorer v2 matches)")

    df = pd.read_csv(
        CATALOG_PATH, sep="\t",
        usecols=["LocusId", "ReferenceRegion", "start_0based", "end_1based", "ReferenceMotif"],
    )
    df = df[df["LocusId"].isin(match_ids)].reset_index(drop=True)
    print(f"{len(df):,} rows after filtering catalog")
    df["chrom"] = df["ReferenceRegion"].str.split(":").str[0]

    fasta = pysam.FastaFile(FASTA_PATH)

    rows = []
    for i, row in enumerate(df.itertuples(index=False)):
        if i % 50_000 == 0:
            print(f"  {i:,} / {len(df):,}")
        core_seq = fasta.fetch(row.chrom, row.start_0based, row.end_1based)
        core_length = len(core_seq)
        left_flank = fasta.fetch(row.chrom, max(0, row.start_0based - MAX_OFFSET), row.start_0based)[::-1]
        right_flank = fasta.fetch(row.chrom, row.end_1based, row.end_1based + MAX_OFFSET)
        left_ext = extend_by_gap_purity(left_flank, row.ReferenceMotif, core_length, "left")
        right_ext = extend_by_gap_purity(right_flank, row.ReferenceMotif, core_length, "right")
        # extend_by_gap_purity only sees MAX_OFFSET bases of flank, so an extension that lands exactly
        # on MAX_OFFSET can't be distinguished from a genuine rule-determined stop -- it may just mean
        # the flank window ran out before the rule did. Flag it rather than silently treat it as final.
        rows.append((row.LocusId, row.chrom, row.start_0based, row.end_1based, row.ReferenceMotif, left_ext, right_ext,
                     left_ext == MAX_OFFSET, right_ext == MAX_OFFSET))

    out_df = pd.DataFrame(rows, columns=["LocusId", "chrom", "start_0based", "end_1based", "motif", "left_ext", "right_ext",
                                          "left_ext_capped", "right_ext_capped"])
    out_df.to_csv(OUTPUT_PATH, sep="\t", index=False)
    print(f"\nWrote {len(out_df):,} rows to {OUTPUT_PATH}")
    n_capped = (out_df["left_ext_capped"] | out_df["right_ext_capped"]).sum()
    if n_capped:
        print(f"WARNING: {n_capped} locus/loci hit the {MAX_OFFSET}bp flank-fetch cap on at least one side -- "
              f"their extension may be truncated, not a genuine rule-determined stop. See left/right_ext_capped.")

    extended = (out_df["left_ext"] > 0) | (out_df["right_ext"] > 0)
    both_sides = (out_df["left_ext"] > 0) & (out_df["right_ext"] > 0)
    print(f"\nExtended: {extended.sum():,} / {len(out_df):,} ({100*extended.mean():.2f}%)")
    print(f"Extended on both sides: {both_sides.sum():,} ({100*both_sides.mean():.2f}%)")
    print(f"Mean extension (extended loci only): {(out_df.loc[extended,'left_ext']+out_df.loc[extended,'right_ext']).mean():.1f}bp")


if __name__ == "__main__":
    main()
