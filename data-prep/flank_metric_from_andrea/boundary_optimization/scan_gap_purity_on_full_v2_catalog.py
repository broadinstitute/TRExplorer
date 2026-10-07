"""Apply the gap-purity boundary-extension rule (gap_purity_extension.py) to every locus in the full
TRExplorer v2 catalog (not just Andrea's 2-6bp-motif subset scanned by
scan_gap_purity_on_v2_exact_matches.py).

Usage:
    python3 scan_gap_purity_on_full_v2_catalog.py
"""

import csv
import gzip
from pathlib import Path

import pysam

from gap_purity_extension import scan_side, MAX_OFFSET, MAX_RESCAN_ATTEMPTS

CATALOG_BED_PATH = Path(
    "~/code/tandem-repeat-catalogs/results__2026-02-01/release_draft_2026-02-01/"
    "TRExplorer.repeat_catalog_v2.hg38.1_to_1000bp_motifs.bed.gz"
).expanduser()
FASTA_PATH = "/Users/weisburd/hg38.fa"
DATA_DIR = Path(__file__).parent / "data"
OUTPUT_PATH = DATA_DIR / "gap_purity_scan_full_v2_catalog.tsv"


def main():
    fasta = pysam.FastaFile(FASTA_PATH)

    n_loci = 0
    n_extended = 0
    n_extended_both_sides = 0
    n_capped = 0

    with gzip.open(CATALOG_BED_PATH, "rt") as bed_file, open(OUTPUT_PATH, "w", newline="") as out_file:
        writer = csv.writer(out_file, delimiter="\t")
        writer.writerow(["chrom", "start_0based", "end_1based", "motif", "left_ext", "right_ext",
                          "left_ext_capped", "right_ext_capped"])

        for line in bed_file:
            chrom, start_0based, end_1based, motif, _ = line.rstrip("\n").split("\t")
            start_0based, end_1based = int(start_0based), int(end_1based)
            core_length = end_1based - start_0based

            left_ext, left_capped = scan_side(fasta, chrom, "left", start_0based, end_1based, motif, core_length)
            right_ext, right_capped = scan_side(fasta, chrom, "right", start_0based, end_1based, motif, core_length)

            writer.writerow([chrom, start_0based, end_1based, motif, left_ext, right_ext, left_capped, right_capped])

            n_loci += 1
            if left_ext > 0 or right_ext > 0:
                n_extended += 1
            if left_ext > 0 and right_ext > 0:
                n_extended_both_sides += 1
            if left_capped or right_capped:
                n_capped += 1

            if n_loci % 500_000 == 0:
                print(f"  {n_loci:,} loci scanned...")

    print(f"\nWrote {n_loci:,} rows to {OUTPUT_PATH}")
    if n_capped:
        print(f"WARNING: {n_capped:,} locus/loci still hit the fetch-window edge on at least one side "
              f"after {MAX_RESCAN_ATTEMPTS} rescan attempts (the last window is "
              f"{3 ** (MAX_RESCAN_ATTEMPTS - 1)}x the first, which is at least {MAX_OFFSET:,}bp "
              f"and grows with the motif) -- their "
              f"extension may still be truncated, not a genuine rule-determined stop. See left/right_ext_capped.")

    print(f"\nTotal loci in TRExplorer v2: {n_loci:,}")
    print(f"Extended: {n_extended:,} / {n_loci:,} ({100*n_extended/n_loci:.2f}%)")
    print(f"Extended on both sides: {n_extended_both_sides:,} ({100*n_extended_both_sides/n_loci:.2f}%)")


if __name__ == "__main__":
    main()
