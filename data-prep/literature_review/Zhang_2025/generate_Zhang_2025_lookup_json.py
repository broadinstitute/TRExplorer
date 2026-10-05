"""Process the STR MPRA results from Zhang et al. 2025 and create a lookup JSON for TRExplorer.

Zhang et al. 2025, "Systematic evaluation of the impact of promoter proximal short tandem repeats on
expression" (bioRxiv, doi:10.1101/2025.09.14.676153), tested 30,516 hg38 STRs within 3.6 kb of a TSS in a
massively parallel reporter assay (MPRA) in HEK293T cells. Each STR was tested at its hg38 length and at
-5, +3 (some loci) and +5 repeat copies. For each STR with enough data, a linear regression of reporter
expression (RNA/DNA ratio) on the change in repeat copy number gives a slope (beta_1) and a
Benjamini-Hochberg adjusted p-value (padj). The paper calls an STR significant at FDR < 0.1.

Inputs (see run.sh):
  - GSE306816_hSTR1_HEK293_linear_regression.csv.gz: per-STR regression results for the hSTR1 library,
    the paper's primary dataset (19,818 STRs with enough data)
  - tss_str_pairs.tsv and array_probes.tsv from https://github.com/gymreklab/str_mpra_design: the hg38
    coordinates (str_pos is 1-based) and motif of every STR in the library (30,516)

Each STR is matched to the TRExplorer catalog the same way as in the Manigbas 2024 script (by exact
locus ID, or else by interval overlap with the same canonical motif and the highest Jaccard similarity),
but with the Tanudisastro 2024 script's lower Jaccard cutoff of 0.2 (see MIN_JACCARD_SIMILARITY).

Output is a JSON lookup table keyed by TRExplorer locus ID. It has an entry for every TRExplorer locus
that matched an STR in the library. A locus without an entry was either not tested by this study, or was
tested by one of the few STRs that could not be matched to the catalog (43 of 30,516 with the v2.1
catalog, for example because the design motif differs from the catalog motif at the same interval). The
unmatched STRs are printed when the script runs.
  - LengthEffectTestResult: "significant" (FDR < 0.1), "not significant", or "not enough data" (in the
    library, but too few barcodes or variants to fit the regression)
  - LengthEffectSlope: change in RNA/DNA expression ratio per added repeat copy (null if not enough data)
  - LengthEffectFDR: Benjamini-Hochberg adjusted p-value of the slope (null if not enough data)
  - Details: zhangLocusId, zhangInterval (0-based start), zhangMotif, matchType ("exact" or "fuzzy")
"""

import argparse
import bisect
import collections
import gzip
import json
import re

import pandas as pd

from str_analysis.utils.canonical_repeat_unit import compute_canonical_motif
from str_analysis.utils.eh_catalog_utils import get_variant_catalog_iterator
from str_analysis.utils.misc_utils import parse_interval

SIGNIFICANCE_FDR_THRESHOLD = 0.1
# Same cutoff as the Tanudisastro 2024 script. Zhang's HipSTR intervals are often about twice as long as
# TRExplorer's for the same repeat (they extend into imperfect flanking repeat), so the 0.66 cutoff of the
# Manigbas 2024 script left 16% of the STRs unmatched.
MIN_JACCARD_SIMILARITY = 0.2
TEST_RESULT_SIGNIFICANT = "significant"
TEST_RESULT_NOT_SIGNIFICANT = "not significant"
TEST_RESULT_NOT_ENOUGH_DATA = "not enough data"


def parse_args():
    parser = argparse.ArgumentParser(description="Generate Zhang 2025 STR MPRA lookup JSON")
    parser.add_argument("--regression-results", default="GSE306816_hSTR1_HEK293_linear_regression.csv.gz",
                        help="Per-STR linear regression results from GEO GSE306816")
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


def classify_test_result(padj):
    """Returns the LengthEffectTestResult for an STR's adjusted p-value (None or NaN if it wasn't analyzed)."""
    if padj is None or pd.isna(padj):
        return TEST_RESULT_NOT_ENOUGH_DATA
    return TEST_RESULT_SIGNIFICANT if padj < SIGNIFICANCE_FDR_THRESHOLD else TEST_RESULT_NOT_SIGNIFICANT


def collect_catalog_candidates(catalog_path, str_intervals_by_chrom, max_str_length):
    """Returns the TRExplorer loci that overlap at least one Zhang STR.

    Only these loci are kept, so memory stays small even though the catalog has millions of loci.

    Args:
        catalog_path: path of the TRExplorer JSON catalog.
        str_intervals_by_chrom: {chrom: sorted list of (start_0based, end_1based)} of the Zhang STRs.
        max_str_length: length of the longest Zhang STR, which bounds the search window.

    Returns:
        (locus_ids, candidates_by_chrom) where locus_ids is the set of candidate TRExplorer locus IDs and
        candidates_by_chrom is {chrom: [dict(locus_id, start_0based, end_1based, canonical_motif, motif_length)]}.
    """
    locus_structure_regex = re.compile(r"^[(]([A-Z]+)[)][+*]")
    str_starts_by_chrom = {chrom: [start for start, _ in intervals] for chrom, intervals in str_intervals_by_chrom.items()}

    locus_ids = set()
    candidates_by_chrom = collections.defaultdict(list)
    for record in get_variant_catalog_iterator(catalog_path, show_progress_bar=False):
        reference_region = record["ReferenceRegion"]
        if isinstance(reference_region, list):
            reference_region = reference_region[0]
        chrom, start_0based, end_1based = parse_interval(reference_region)
        chrom = chrom.replace("chr", "")
        if chrom not in str_starts_by_chrom:
            continue

        # STRs that can overlap this locus start in (start_0based - max_str_length, end_1based).
        starts = str_starts_by_chrom[chrom]
        i = bisect.bisect_right(starts, start_0based - max_str_length)
        j = bisect.bisect_left(starts, end_1based)
        if not any(str_end > start_0based for _, str_end in str_intervals_by_chrom[chrom][i:j]):
            continue

        if "ReferenceMotif" in record:
            motif = record["ReferenceMotif"]
        else:
            match = locus_structure_regex.match(record["LocusStructure"])
            if not match:
                continue
            motif = match.group(1)

        locus_ids.add(record["LocusId"])
        candidates_by_chrom[chrom].append({
            "locus_id": record["LocusId"],
            "start_0based": start_0based,
            "end_1based": end_1based,
            "canonical_motif": compute_canonical_motif(motif),
            "motif_length": len(motif),
        })
    return locus_ids, candidates_by_chrom


def find_trexplorer_match(chrom, start_0based, end_1based, motif, catalog_locus_ids, candidates_by_chrom):
    """Returns (trexplorer_locus_id, match_type) for a Zhang STR, or (None, None) if nothing matches.

    Uses the same rules as the Manigbas 2024 script, apart from the Jaccard cutoff: an exact locus ID match, or else the overlapping
    locus with the same canonical motif (or the same motif length for motifs over 6 bp) and the highest
    Jaccard similarity, if that similarity is at least MIN_JACCARD_SIMILARITY.

    Args:
        chrom: chromosome without "chr".
        start_0based: 0-based start.
        end_1based: 1-based end.
        motif: repeat motif.
        catalog_locus_ids: set of TRExplorer locus IDs to check for an exact match.
        candidates_by_chrom: {chrom: [candidate dict]} returned by collect_catalog_candidates.

    Returns:
        (locus ID, "exact" or "fuzzy"), or (None, None).
    """
    zhang_locus_id = f"{chrom}-{start_0based}-{end_1based}-{motif}"
    if zhang_locus_id in catalog_locus_ids:
        return zhang_locus_id, "exact"

    canonical_motif = compute_canonical_motif(motif)
    best_match = None
    best_jaccard = 0
    for candidate in candidates_by_chrom.get(chrom, []):
        intersection = min(end_1based, candidate["end_1based"]) - max(start_0based, candidate["start_0based"])
        if intersection <= 0:
            continue
        if len(motif) <= 6:
            if canonical_motif != candidate["canonical_motif"]:
                continue
        elif len(motif) != candidate["motif_length"]:
            continue

        union = (end_1based - start_0based) + (candidate["end_1based"] - candidate["start_0based"]) - intersection
        jaccard = intersection / union
        if jaccard >= MIN_JACCARD_SIMILARITY and jaccard > best_jaccard:
            best_jaccard = jaccard
            best_match = candidate["locus_id"]

    return (best_match, "fuzzy") if best_match else (None, None)


def pick_best_entry(entries):
    """Returns the entry to keep when several Zhang STRs match the same TRExplorer locus.

    Prefers STRs with regression results over those without, then the lowest adjusted p-value.
    """
    return min(entries, key=lambda entry: (entry["LengthEffectFDR"] is None, entry["LengthEffectFDR"] or 0))


def main():
    args = parse_args()

    strs = read_zhang_strs(args.tss_str_pairs, args.array_probes)
    print(f"Read {len(strs):,d} human STRs in the MPRA library")

    results = pd.read_csv(args.regression_results)
    print(f"Read regression results for {len(results):,d} STRs from {args.regression_results}")
    strs = strs.merge(results[["locus", "repeat_unit", "beta_1", "padj"]],
                      left_on="zhang_locus_number", right_on="locus", how="left", validate="one_to_one")
    if strs["locus"].notna().sum() != len(results):
        raise ValueError(f"Only {strs['locus'].notna().sum():,d} of {len(results):,d} regression results "
                         f"matched an STR in the library")
    analyzed = strs["locus"].notna()
    n_motif_mismatches = (strs.loc[analyzed, "repeat_unit"] != strs.loc[analyzed, "motif"]).sum()
    print(f"  {n_motif_mismatches:,d} analyzed STRs have a regression repeat_unit that differs from the design motif")

    str_intervals_by_chrom = {
        chrom: sorted(zip(group["start_0based"], group["end_1based"])) for chrom, group in strs.groupby("chrom")
    }
    max_str_length = int((strs["end_1based"] - strs["start_0based"]).max())
    print(f"Reading TRExplorer catalog: {args.trexplorer_catalog}")
    catalog_locus_ids, candidates_by_chrom = collect_catalog_candidates(
        args.trexplorer_catalog, str_intervals_by_chrom, max_str_length)
    print(f"  Kept {len(catalog_locus_ids):,d} TRExplorer loci that overlap a Zhang STR")

    entries_by_trexplorer_locus = collections.defaultdict(list)
    match_type_counts = collections.Counter()
    unmatched_test_results = collections.Counter()
    unmatched_strs = []
    for row in strs.itertuples():
        test_result = classify_test_result(row.padj)
        trexplorer_locus_id, match_type = find_trexplorer_match(
            row.chrom, row.start_0based, row.end_1based, row.motif, catalog_locus_ids, candidates_by_chrom)
        if trexplorer_locus_id is None:
            match_type_counts["no match"] += 1
            unmatched_test_results[test_result] += 1
            unmatched_strs.append(f"{row.zhang_locus_id} chr{row.chrom}:{row.start_0based}-{row.end_1based} "
                                  f"{row.motif} ({test_result})")
            continue
        match_type_counts[match_type] += 1
        entries_by_trexplorer_locus[trexplorer_locus_id].append({
            "LengthEffectTestResult": test_result,
            "LengthEffectSlope": None if pd.isna(row.beta_1) else float(row.beta_1),
            "LengthEffectFDR": None if pd.isna(row.padj) else float(row.padj),
            "Details": {
                "zhangLocusId": row.zhang_locus_id,
                "zhangInterval": f"chr{row.chrom}:{row.start_0based}-{row.end_1based}",
                "zhangMotif": row.motif,
                "matchType": match_type,
            },
        })

    output_lookup = {locus_id: pick_best_entry(entries) for locus_id, entries in entries_by_trexplorer_locus.items()}
    n_loci_with_multiple_strs = sum(len(entries) > 1 for entries in entries_by_trexplorer_locus.values())

    print(f"\nZhang 2025 STR matching statistics:")
    for match_type in ["exact", "fuzzy", "no match"]:
        print(f"  {match_type + ':':10s} {match_type_counts[match_type]:,d}")
    print(f"  Unmatched STRs by test result: {dict(unmatched_test_results)}")
    for unmatched_str in unmatched_strs:
        print(f"    unmatched: {unmatched_str}")
    print(f"  TRExplorer loci matched by more than one STR (kept the best): {n_loci_with_multiple_strs:,d}")

    test_result_counts = collections.Counter(entry["LengthEffectTestResult"] for entry in output_lookup.values())
    print(f"\nOutput: {len(output_lookup):,d} TRExplorer loci")
    for test_result in [TEST_RESULT_SIGNIFICANT, TEST_RESULT_NOT_SIGNIFICANT, TEST_RESULT_NOT_ENOUGH_DATA]:
        print(f"  {test_result + ':':18s} {test_result_counts[test_result]:,d}")

    fopen = gzip.open if args.output.endswith("gz") else open
    with fopen(args.output, "wt") as f:
        json.dump(output_lookup, f, indent=2)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
