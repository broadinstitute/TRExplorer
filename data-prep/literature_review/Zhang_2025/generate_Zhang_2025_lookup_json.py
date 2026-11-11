"""Process the STR MPRA results from Zhang et al. 2025 and create a lookup JSON for TRExplorer.

Zhang et al. 2025, "Systematic evaluation of the impact of promoter proximal short tandem repeats on
expression" (bioRxiv, doi:10.1101/2025.09.14.676153), tested 30,516 hg38 STRs within 3.6 kb of a TSS in a
massively parallel reporter assay (MPRA) in HEK293T cells. Each STR was tested at its hg38 length and at
-5, +3 (some loci) and +5 repeat copies. For each STR with enough data, a linear regression of reporter
expression (RNA/DNA ratio) on the change in repeat copy number gives a slope (beta_1), a Pearson r and a
Benjamini-Hochberg adjusted p-value (padj). The paper calls an STR significant at FDR < 0.1.

The same oligo pool was cloned into two libraries: hSTR1 (the paper's primary dataset) and hSTR2 (more
barcodes per variant but fewer reads per barcode, so less power). The LengthEffect fields come from hSTR1
when it had enough data to test the STR, and otherwise from hSTR2. The two libraries' slopes are on
different scales (among STRs significant in both, hSTR2 slopes are about half of hSTR1's), so Details
records which library the LengthEffect fields came from.

Inputs (see run.sh):
  - GSE306816_hSTR1_HEK293_linear_regression.csv.gz: per-STR regression results for the hSTR1 library
    (19,818 STRs with enough data)
  - GSE306816_hSTR2_HEK293_linear_regression.csv.gz: the same for the hSTR2 library (15,632 STRs). Despite
    its .csv name, this file is tab-separated.
  - tss_str_pairs.tsv and array_probes.tsv from https://github.com/gymreklab/str_mpra_design: the hg38
    coordinates (str_pos is 1-based) and motif of every STR in the library (30,516)

Each STR is matched to the TRExplorer catalog the same way as in the Manigbas 2024 script (by exact
locus ID, or else by interval overlap with the same canonical motif and the highest Jaccard similarity, using
../trexplorer_catalog_matching.find_exact_or_best_jaccard_match), but with the Tanudisastro 2024 script's lower
Jaccard cutoff of 0.2 (see MIN_JACCARD_SIMILARITY). A match with the same interval and canonical motif is
"exact" too: Zhang et al. write each motif starting at its alphabetically first base, so the same locus can
have a motif that starts at a different base than in its TRExplorer locus ID.

Output is a JSON lookup table keyed by TRExplorer locus ID. It has an entry for every TRExplorer locus
that matched an STR in the library. A locus without an entry was either not tested by this study, or was
tested by one of the few STRs that could not be matched to the catalog (43 of 30,516 with the v2.1
catalog, for example because the design motif differs from the catalog motif at the same interval). The
unmatched STRs are printed when the script runs.
  - LengthEffectTestResult: "significant" (FDR < 0.1), "not significant", or "not enough data" (in the
    library, but too few barcodes or variants to fit the regression in either hSTR1 or hSTR2)
  - LengthEffectSlope: change in RNA/DNA expression ratio per added repeat copy (null if not enough data)
  - LengthEffectFDR: Benjamini-Hochberg adjusted p-value of the slope (null if not enough data)
  - Details: zhangLocusId, zhangInterval (0-based start), zhangMotif, matchType ("exact" for the same interval and
    canonical motif, or "fuzzy"),
    lengthEffectSourceLibrary ("hSTR1", "hSTR2", or null if not enough data), and an hSTR1 and/or hSTR2
    object (lengthEffectSlope, lengthEffectFDR, pearsonR) for each library that had enough data
"""

import argparse
import collections
import gzip
import json
from pathlib import Path
import sys

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from trexplorer_catalog_matching import build_catalog_interval_trees, find_exact_or_best_jaccard_match  # noqa: E402

SIGNIFICANCE_FDR_THRESHOLD = 0.1
# Same cutoff as the Tanudisastro 2024 script. Zhang's HipSTR intervals are often about twice as long as
# TRExplorer's for the same repeat (they extend into imperfect flanking repeat), so the 0.66 cutoff of the
# Manigbas 2024 script left 16% of the STRs unmatched.
MIN_JACCARD_SIMILARITY = 0.2
TEST_RESULT_SIGNIFICANT = "significant"
TEST_RESULT_NOT_SIGNIFICANT = "not significant"
TEST_RESULT_NOT_ENOUGH_DATA = "not enough data"
# The LengthEffect fields come from the first library in this list that had enough data to test the STR
LIBRARIES_IN_ORDER_OF_PREFERENCE = ["hSTR1", "hSTR2"]
REGRESSION_RESULTS_COLUMNS = ["locus", "repeat_unit", "r", "beta_1", "padj"]


def parse_args():
    parser = argparse.ArgumentParser(description="Generate Zhang 2025 STR MPRA lookup JSON")
    parser.add_argument("--hstr1-regression-results", default="GSE306816_hSTR1_HEK293_linear_regression.csv.gz",
                        help="Per-STR linear regression results for the hSTR1 library from GEO GSE306816 (comma-separated)")
    parser.add_argument("--hstr2-regression-results", default="GSE306816_hSTR2_HEK293_linear_regression.csv.gz",
                        help="Per-STR linear regression results for the hSTR2 library from GEO GSE306816 (tab-separated)")
    parser.add_argument("--tss-str-pairs", default="tss_str_pairs.tsv",
                        help="tss_str_pairs.tsv from the str_mpra_design repo (STR coordinates)")
    parser.add_argument("--array-probes", default="array_probes.tsv",
                        help="array_probes.tsv from the str_mpra_design repo (STR motifs)")
    parser.add_argument("--trexplorer-catalog", required=True, help="Path to TRExplorer v2.1 JSON catalog")
    parser.add_argument("-o", "--output", default="Zhang_2025_lookup.json.gz",
                        help="Output JSON path (default: Zhang_2025_lookup.json.gz)")
    return parser.parse_args()


def read_zhang_strs(tss_str_pairs_path, array_probes_path):
    """Returns one row per human STR in the MPRA library, with its hg38 interval and motif.

    Args:
        tss_str_pairs_path: path of tss_str_pairs.tsv (several rows per STR when it is near several TSSs)
        array_probes_path: path of array_probes.tsv (one row per probe, several probes per STR)

    Returns:
        DataFrame with columns zhang_locus_number, zhang_locus_id, chrom (without "chr"), start_0based,
        end_1based, motif.
    """
    pairs = pd.read_csv(tss_str_pairs_path, sep="\t", usecols=["organism", "str_chr", "str_pos", "str_end", "str_id"])
    strs = pairs[pairs["organism"] == "hg38"].drop_duplicates("str_id")

    probes = pd.read_csv(array_probes_path, sep="\t", usecols=["organism", "id", "motif"])
    motifs = probes[probes["organism"] == "hg38"].drop_duplicates("id").rename(columns={"id": "str_id"})

    strs = strs.merge(motifs[["str_id", "motif"]], on="str_id", how="left", validate="one_to_one")
    if strs["motif"].isna().any():
        raise ValueError(f"{strs['motif'].isna().sum():,d} STRs in {tss_str_pairs_path} have no motif in {array_probes_path}")

    return pd.DataFrame({
        "zhang_locus_number": strs["str_id"].str.rsplit("_", n=1).str[1].astype(int),
        "zhang_locus_id": strs["str_id"],
        "chrom": strs["str_chr"].str.replace("chr", "", regex=False),
        "start_0based": strs["str_pos"].astype(int) - 1,
        "end_1based": strs["str_end"].astype(int),
        "motif": strs["motif"],
    })


def merge_regression_results(strs, results, library):
    """Adds one library's regression results to the table of STRs.

    Args:
        strs: DataFrame returned by read_zhang_strs.
        results: DataFrame with one row of regression results per STR that the library had enough data to test.
        library: library name ("hSTR1" or "hSTR2"), used as the prefix of the added columns.

    Returns:
        A copy of strs with the columns <library>_slope, <library>_fdr and <library>_pearson_r added, which are NaN
        for STRs the library didn't have enough data to test.

    Raises:
        ValueError: if results is missing a column, or has a row that doesn't match an STR in strs.
    """
    missing_columns = set(REGRESSION_RESULTS_COLUMNS) - set(results.columns)
    if missing_columns:
        raise ValueError(f"{library} regression results are missing the columns {sorted(missing_columns)}. "
                         f"Found: {list(results.columns)}")

    results = results[REGRESSION_RESULTS_COLUMNS].rename(columns={
        "locus": "zhang_locus_number",
        "repeat_unit": f"{library}_repeat_unit",
        "beta_1": f"{library}_slope",
        "padj": f"{library}_fdr",
        "r": f"{library}_pearson_r",
    })
    merged = strs.merge(results, on="zhang_locus_number", how="left", validate="one_to_one")
    analyzed = merged[f"{library}_repeat_unit"].notna()
    if analyzed.sum() != len(results):
        raise ValueError(f"Only {analyzed.sum():,d} of {len(results):,d} {library} regression results matched an STR "
                         f"in the library")
    n_motif_mismatches = (merged.loc[analyzed, f"{library}_repeat_unit"] != merged.loc[analyzed, "motif"]).sum()
    print(f"  {n_motif_mismatches:,d} STRs analyzed in {library} have a regression repeat_unit that differs from the "
          f"design motif")

    return merged.drop(columns=f"{library}_repeat_unit")


def classify_test_result(padj):
    """Returns the LengthEffectTestResult for an STR's adjusted p-value (None or NaN if it wasn't analyzed)."""
    if padj is None or pd.isna(padj):
        return TEST_RESULT_NOT_ENOUGH_DATA
    return TEST_RESULT_SIGNIFICANT if padj < SIGNIFICANCE_FDR_THRESHOLD else TEST_RESULT_NOT_SIGNIFICANT


def build_lookup_entry(results_by_library, details):
    """Returns the lookup entry for one Zhang STR.

    Args:
        results_by_library: {library: {"lengthEffectSlope": float, "lengthEffectFDR": float, "pearsonR": float}}
            for each library that had enough data to test the STR.
        details: dict with the STR's zhangLocusId, zhangInterval, zhangMotif and matchType.

    Returns:
        dict whose LengthEffectTestResult, LengthEffectSlope and LengthEffectFDR come from the first library in
        LIBRARIES_IN_ORDER_OF_PREFERENCE that tested the STR, and whose Details add that library's name
        (lengthEffectSourceLibrary) and the results of every library that tested the STR.
    """
    source_library = next(
        (library for library in LIBRARIES_IN_ORDER_OF_PREFERENCE if library in results_by_library), None)
    source_library_results = results_by_library.get(source_library, {})
    return {
        "LengthEffectTestResult": classify_test_result(source_library_results.get("lengthEffectFDR")),
        "LengthEffectSlope": source_library_results.get("lengthEffectSlope"),
        "LengthEffectFDR": source_library_results.get("lengthEffectFDR"),
        "Details": {**details, "lengthEffectSourceLibrary": source_library, **results_by_library},
    }


def pick_best_entry(entries):
    """Returns the entry to keep when several Zhang STRs match the same TRExplorer locus.

    Prefers STRs whose LengthEffect fields come from hSTR1, then from hSTR2, then STRs without enough data, and
    within each of these groups the lowest adjusted p-value.
    """
    def sort_key(entry):
        source_library = entry["Details"]["lengthEffectSourceLibrary"]
        if source_library is None:
            library_rank = len(LIBRARIES_IN_ORDER_OF_PREFERENCE)
        else:
            library_rank = LIBRARIES_IN_ORDER_OF_PREFERENCE.index(source_library)
        return library_rank, entry["LengthEffectFDR"] or 0

    return min(entries, key=sort_key)


def main():
    args = parse_args()

    strs = read_zhang_strs(args.tss_str_pairs, args.array_probes)
    print(f"Read {len(strs):,d} human STRs in the MPRA library")

    # GEO's hSTR1 file is comma-separated, but its hSTR2 file is tab-separated despite its .csv name.
    # merge_regression_results raises an error if a file is read with the wrong separator.
    for library, results_path, separator in [
        ("hSTR1", args.hstr1_regression_results, ","),
        ("hSTR2", args.hstr2_regression_results, "\t"),
    ]:
        results = pd.read_csv(results_path, sep=separator)
        print(f"Read {library} regression results for {len(results):,d} STRs from {results_path}")
        strs = merge_regression_results(strs, results, library)

    print(f"Reading TRExplorer catalog: {args.trexplorer_catalog}")
    # Only the catalog loci that overlap a Zhang STR are kept, so memory stays small
    catalog_interval_trees = build_catalog_interval_trees(
        args.trexplorer_catalog, only_loci_that_overlap=zip(strs["chrom"], strs["start_0based"], strs["end_1based"]),
        show_progress_bar=False)

    entries_by_trexplorer_locus = collections.defaultdict(list)
    match_type_counts = collections.Counter()
    unmatched_test_results = collections.Counter()
    unmatched_strs = []
    for row in strs.itertuples():
        results_by_library = {}
        for library in LIBRARIES_IN_ORDER_OF_PREFERENCE:
            if pd.isna(getattr(row, f"{library}_fdr")):
                continue
            results_by_library[library] = {
                "lengthEffectSlope": float(getattr(row, f"{library}_slope")),
                "lengthEffectFDR": float(getattr(row, f"{library}_fdr")),
                "pearsonR": float(getattr(row, f"{library}_pearson_r")),
            }

        trexplorer_locus_id, match_type = find_exact_or_best_jaccard_match(
            catalog_interval_trees, row.chrom, row.start_0based, row.end_1based, row.motif,
            MIN_JACCARD_SIMILARITY)
        entry = build_lookup_entry(results_by_library, {
            "zhangLocusId": row.zhang_locus_id,
            "zhangInterval": f"chr{row.chrom}:{row.start_0based}-{row.end_1based}",
            "zhangMotif": row.motif,
            "matchType": match_type,
        })
        if trexplorer_locus_id is None:
            match_type_counts["no match"] += 1
            unmatched_test_results[entry["LengthEffectTestResult"]] += 1
            unmatched_strs.append(f"{row.zhang_locus_id} chr{row.chrom}:{row.start_0based}-{row.end_1based} "
                                  f"{row.motif} ({entry['LengthEffectTestResult']}, "
                                  f"from {entry['Details']['lengthEffectSourceLibrary']})")
            continue
        match_type_counts[match_type] += 1
        entries_by_trexplorer_locus[trexplorer_locus_id].append(entry)

    output_lookup = {locus_id: pick_best_entry(entries) for locus_id, entries in entries_by_trexplorer_locus.items()}
    n_loci_with_multiple_strs = sum(len(entries) > 1 for entries in entries_by_trexplorer_locus.values())

    print(f"\nZhang 2025 STR matching statistics:")
    for match_type in ["exact", "fuzzy", "no match"]:
        print(f"  {match_type + ':':10s} {match_type_counts[match_type]:,d}")
    print(f"  Unmatched STRs by test result: {dict(unmatched_test_results)}")
    for unmatched_str in unmatched_strs:
        print(f"    unmatched: {unmatched_str}")
    print(f"  TRExplorer loci matched by more than one STR (kept the best): {n_loci_with_multiple_strs:,d}")

    test_result_counts = collections.Counter(
        (entry["LengthEffectTestResult"], entry["Details"]["lengthEffectSourceLibrary"])
        for entry in output_lookup.values())
    print(f"\nOutput: {len(output_lookup):,d} TRExplorer loci")
    for test_result in [TEST_RESULT_SIGNIFICANT, TEST_RESULT_NOT_SIGNIFICANT, TEST_RESULT_NOT_ENOUGH_DATA]:
        counts_by_library = {
            library: count for (result, library), count in test_result_counts.items() if result == test_result}
        print(f"  {test_result + ':':18s} {sum(counts_by_library.values()):,d}  "
              f"(LengthEffect fields from: {counts_by_library})")
    n_tested_in_both_libraries = sum("hSTR1" in entry["Details"] and "hSTR2" in entry["Details"]
                                     for entry in output_lookup.values())
    print(f"  Loci tested in both hSTR1 and hSTR2: {n_tested_in_both_libraries:,d}")

    fopen = gzip.open if args.output.endswith("gz") else open
    with fopen(args.output, "wt") as f:
        json.dump(output_lookup, f, indent=2)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
