"""Subset extended_definitions.bed to only the loci whose ORIGINAL definition (chrom/start/end/motif)
exactly matches an atomic locus in TRExplorer v2's own bed (data/trexplorer_v2_locus_ids.sorted.txt,
built by the earlier ad hoc comparison -- see conversation). Loci from external-source annotations
(e.g. TRExplorerV2:PolymorphicTRsInT2TAssembliesV2) mostly fail this check and are dropped.

Usage:
    python3 subset_bed_to_v2_exact_matches.py
"""

from pathlib import Path

from generate_extended_definitions_bed import round2_rows, round3_rows, round4_rows

DATA_DIR = Path(__file__).parent / "data"
V2_IDS_PATH = DATA_DIR / "trexplorer_v2_locus_ids.sorted.txt"
OUTPUT_PATH = DATA_DIR / "extended_definitions.subset_to_v2_exact_matches.bed"


def main():
    v2_ids = set(V2_IDS_PATH.read_text().split())

    all_rows = round2_rows() + round3_rows() + round4_rows()
    kept, dropped = [], []
    for round_name, locus_id, chrom, ostart, oend, nstart, nend, motif in all_rows:
        original_id = f"{chrom.removeprefix('chr')}-{ostart}-{oend}-{motif}"
        (kept if original_id in v2_ids else dropped).append(
            (round_name, locus_id, chrom, ostart, oend, nstart, nend)
        )

    print(f"{len(kept)} / {len(all_rows)} loci have an original definition that exactly matches TRExplorer v2")
    for round_name in ["round2", "round3", "round4"]:
        n_kept = sum(1 for r in kept if r[0] == round_name)
        n_total = sum(1 for r in all_rows if r[0] == round_name)
        print(f"  {round_name}: {n_kept} / {n_total}")

    bed_rows = []
    for round_name, locus_id, chrom, ostart, oend, nstart, nend in kept:
        bed_rows.append((chrom, ostart, oend, f"{locus_id}__{round_name}__original", 0, "."))
        bed_rows.append((chrom, nstart, nend, f"{locus_id}__{round_name}__new", 0, "."))
    bed_rows.sort(key=lambda r: (r[0], r[1], r[2]))

    with open(OUTPUT_PATH, "w") as f:
        for chrom, start, end, name, score, strand in bed_rows:
            f.write(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\n")

    print(f"\nWrote {len(bed_rows)} BED entries ({len(kept)} loci) to {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
