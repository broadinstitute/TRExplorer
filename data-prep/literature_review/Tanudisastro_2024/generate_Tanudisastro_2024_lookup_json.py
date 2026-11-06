"""Match the sc-eTRs from Tanudisastro et al. 2026 (Nature) to TRExplorer loci and write a lookup JSON.

The paper (doi 10.1038/s41586-026-10438-2, first posted on bioRxiv in 2024 as 10.1101/2024.11.02.621562) tested
1.7 million TRs for association with gene expression in 28 immune cell types of 1,925 TenK10K individuals. This script
reads three of its supplementary tables, from Nature's 41586_2026_10438_MOESM3_ESM.zip:

  - Table S1: the sc-eTRs (single-cell expression TR loci), one row per TR, cell type and eGene, at gene-level FDR < 5%
  - Table S6: the candidate causal sc-eTRs, with their SuSiE fine-mapping posterior inclusion probability (PIP >= 0.7).
    498 of its 2,541 rows have no matching Table S1 row, even though all of them have p <= 5e-8. They are added to
    the Table S1 associations, and 185 of their loci are not in Table S1 at all.
  - Table S4: the sc-eTRs that are also associated with chromatin accessibility at a scATAC-seq peak within 10 kb of
    the TR (caQTLs)

Each Tanudisastro locus is matched to the TRExplorer catalog by locus ID (exact) or by interval overlap and motif
(fuzzy). Table S6 gives only start coordinates, so a Table S6 locus that's not in Table S1 is matched to the TRExplorer
locus with the same motif that contains its start, preferring the closest start. When several Tanudisastro loci match
the same TRExplorer locus, the one with the smallest p-value is kept.

The output JSON is keyed by TRExplorer locus ID, and each entry has:
  - SignificantCellTypes: comma-separated cell types, ordered by their smallest p-value
  - MinPvalue: the smallest p-value across the locus's associations
  - FineMappingMaxPIP: the largest fine-mapping PIP, or None if the locus wasn't fine-mapped
  - ChromatinAccessibilityCellTypes: comma-separated cell types in which the TR is a caQTL, ordered by the smallest
    caQTL p-value, or None
  - Details: {"associations": [...]}, one association per cell type and eGene, most significant first, each with
    cellType, eGene (gene symbol, or Ensembl ID for genes without one), pvalue, effectSize, and fineMappingPIP for the
    fine-mapped ones. Details also has tanudisastroLocusId when the paper's locus has a different interval or
    canonical motif than the TRExplorer locus (matched by overlap). It's chrom-start-end-motif, or chrom-start-motif for
    a locus that is only in Table S6.
"""

import argparse
import collections
import gzip
import json
import math
from pathlib import Path
import re

import intervaltree
import pandas as pd
import tqdm

from str_analysis.utils.canonical_repeat_unit import compute_canonical_motif
from str_analysis.utils.eh_catalog_utils import get_variant_catalog_iterator
from str_analysis.utils.misc_utils import parse_interval

DEFAULT_SUPPLEMENTARY_TABLES_DIR = Path(__file__).parent / "41586_2026_10438_MOESM3_ESM"


def parse_args():
    parser = argparse.ArgumentParser(description="Generate the Tanudisastro et al. 2026 sc-eTR lookup JSON")
    parser.add_argument("--catalog-path", required=True,
                        help="Path to TRExplorer catalog JSON (e.g., TRExplorer.repeat_catalog_v2.1.hg38.1_to_1000bp_motifs.EH.json.gz)")
    parser.add_argument("--supplementary-tables-dir", default=str(DEFAULT_SUPPLEMENTARY_TABLES_DIR),
                        help="Directory with the paper's TableS1.xlsx, TableS4.xlsx and TableS6.xlsx "
                             f"(default: {DEFAULT_SUPPLEMENTARY_TABLES_DIR.name} next to this script)")
    parser.add_argument("--output", default="Tanudisastro_2024_lookup.json.gz",
                        help="Output JSON path (default: Tanudisastro_2024_lookup.json.gz)")
    return parser.parse_args()


def read_supplementary_table(supplementary_tables_dir, table_name):
    """Reads one of the paper's supplementary tables (e.g. TableS1) from its xlsx file, as strings."""
    return pd.read_excel(Path(supplementary_tables_dir) / f"{table_name}.xlsx", sheet_name=table_name, dtype=str)


def parse_table_s1(df):
    """Returns Table S1 (one row per sc-eTR, cell type and eGene) with chrom (no 'chr' prefix), start_0based,
    end_1based, motif, cell_type, egene_id (Ensembl ID), egene (gene symbol, or the Ensembl ID for genes without one),
    pvalue and effect_size columns."""
    return pd.DataFrame({
        "chrom": df["chromosome"].str.replace("chr", "", regex=False),
        "start_0based": df["start coordinate (hg38)"].astype(int),
        "end_1based": df["end coordinate (hg38)"].astype(int),
        "motif": df["motif"],
        "cell_type": df["cell_type"],
        "egene_id": df["eGene"],
        "egene": df["eGene_symbol"],
        "pvalue": df["nominal p-value"].astype(float),
        "effect_size": df["effect size"].astype(float),
    })


def parse_table_s6(df):
    """Returns Table S6 (one row per fine-mapped sc-eTR, cell type and eGene; start coordinates only) with chrom,
    start_0based, motif, cell_type, egene_id, egene, pvalue, effect_size and pip columns."""
    return pd.DataFrame({
        "chrom": df["chromosome"].str.replace("chr", "", regex=False),
        "start_0based": df["start coordinate (hg38)"].astype(int),
        "motif": df["motif"],
        "cell_type": df["cell_type"],
        "egene_id": df["eGene"],
        "egene": df["eGene_symbol"],
        "pvalue": df["nominal p-value"].astype(float),
        "effect_size": df["effect size"].astype(float),
        "pip": df["pip"].astype(float),
    })


def parse_table_s4(df):
    """Returns Table S4 (one row per caQTL: sc-eTR, cell type, eGene and scATAC-seq peak) with chrom, start_0based,
    end_1based, motif, cell_type and atac_pvalue columns."""
    return pd.DataFrame({
        "chrom": df["chromosome"].str.replace("chr", "", regex=False),
        "start_0based": df["start coordinate (hg38)"].astype(int),
        "end_1based": df["end coordinate (hg38)"].astype(int),
        "motif": df["motif"],
        "cell_type": df["cell_type"],
        "atac_pvalue": df["atac_pval"].astype(float),
    })


def build_tanudisastro_loci(s1, s6, s4):
    """Combines the parsed Tables S1, S6 and S4 into one record per Tanudisastro locus.

    A Table S6 row is added to the Table S1 loci with the same chrom, start and motif. Table S6 has no end coordinate to
    tell apart several such loci (with different ends), so its PIP goes to every one of them that has the row's cell
    type and eGene, or if none does, the row is added to the first one. Table S6 rows of loci that are not in Table S1
    form loci of their own, with end_1based None.

    Args:
        s1 (pd.DataFrame): parse_table_s1 output
        s6 (pd.DataFrame): parse_table_s6 output
        s4 (pd.DataFrame): parse_table_s4 output

    Returns:
        dict: (chrom, start_0based, end_1based, motif) -> {"associations": {(cell_type, egene_id): association dict
            with cell_type, egene, pvalue, effect_size and pip (None unless fine-mapped)},
            "caqtl_min_pvalue_by_cell_type": {cell_type: smallest caQTL p-value}}
    """
    loci = collections.defaultdict(lambda: {"associations": {}, "caqtl_min_pvalue_by_cell_type": {}})
    for row in s1.itertuples(index=False):
        locus = loci[(row.chrom, row.start_0based, row.end_1based, row.motif)]
        key = (row.cell_type, row.egene_id)
        if key in locus["associations"]:
            raise ValueError(f"Table S1 has more than one row for {row.chrom}:{row.start_0based} {row.motif} {key}")
        locus["associations"][key] = {"cell_type": row.cell_type, "egene": row.egene, "pvalue": row.pvalue,
                                      "effect_size": row.effect_size, "pip": None}

    s1_locus_keys_by_start_and_motif = collections.defaultdict(list)
    for locus_key in sorted(loci):
        chrom, start_0based, _, motif = locus_key
        s1_locus_keys_by_start_and_motif[(chrom, start_0based, motif)].append(locus_key)

    for row in s6.itertuples(index=False):
        key = (row.cell_type, row.egene_id)
        candidate_locus_keys = s1_locus_keys_by_start_and_motif.get((row.chrom, row.start_0based, row.motif), [])
        locus_keys = ([k for k in candidate_locus_keys if key in loci[k]["associations"]]
                      or candidate_locus_keys[:1] or [(row.chrom, row.start_0based, None, row.motif)])
        for locus_key in locus_keys:
            associations = loci[locus_key]["associations"]
            if key not in associations:
                associations[key] = {"cell_type": row.cell_type, "egene": row.egene, "pvalue": row.pvalue,
                                     "effect_size": row.effect_size, "pip": None}
            associations[key]["pip"] = max(row.pip, associations[key]["pip"] or 0)

    for row in s4.itertuples(index=False):
        locus_key = (row.chrom, row.start_0based, row.end_1based, row.motif)
        if locus_key not in loci:
            raise ValueError(f"Table S4 locus {locus_key} is not in Table S1")
        min_pvalues = loci[locus_key]["caqtl_min_pvalue_by_cell_type"]
        min_pvalues[row.cell_type] = min(row.atac_pvalue, min_pvalues.get(row.cell_type, math.inf))

    return dict(loci)


def get_min_pvalue(tanudisastro_locus):
    """Returns the smallest p-value among a Tanudisastro locus's associations."""
    return min(a["pvalue"] for a in tanudisastro_locus["associations"].values())


def is_same_locus(tanudisastro_locus_key, trexplorer_locus_id):
    """Returns whether a Tanudisastro locus is the TRExplorer locus: the same start, the same end (unless the
    Tanudisastro locus has none, as for a locus that is only in Table S6) and the same canonical motif, so that a motif
    written starting at a different base or on the other strand still counts (e.g. the paper's
    9-27573528-27573546-CCCCGG is TRExplorer's 9-27573528-27573546-GGCCCC)."""
    chrom, start_0based, end_1based, motif = tanudisastro_locus_key
    trexplorer_chrom, trexplorer_start_0based, trexplorer_end_1based, trexplorer_motif = trexplorer_locus_id.split("-")
    return (chrom == trexplorer_chrom and start_0based == int(trexplorer_start_0based)
            and end_1based in (None, int(trexplorer_end_1based))
            and compute_canonical_motif(motif) == compute_canonical_motif(trexplorer_motif))


def pick_tanudisastro_locus_for_trexplorer_locus(candidate_keys, tanudisastro_loci, trexplorer_locus_id):
    """Picks which of the Tanudisastro loci that matched a TRExplorer locus to show for it: the one with the smallest
    p-value, and among those with the same smallest p-value (often 0, for p-values too small for a double-precision
    number), the one that is the TRExplorer locus itself (is_same_locus), then the one that sorts first.

    Args:
        candidate_keys (list): the (chrom, start_0based, end_1based, motif) keys of the matching Tanudisastro loci
        tanudisastro_loci (dict): build_tanudisastro_loci output
        trexplorer_locus_id (str): the TRExplorer locus

    Returns:
        tuple: the key of the picked Tanudisastro locus
    """
    return min(candidate_keys, key=lambda k: (
        get_min_pvalue(tanudisastro_loci[k]), not is_same_locus(k, trexplorer_locus_id), str(k)))


def build_lookup_entry(tanudisastro_locus_key, tanudisastro_locus, trexplorer_locus_id):
    """Returns the lookup JSON entry of a TRExplorer locus from the Tanudisastro locus it was matched to.

    Args:
        tanudisastro_locus_key (tuple): (chrom, start_0based, end_1based, motif), with end_1based None for a locus that
            is only in Table S6
        tanudisastro_locus (dict): the build_tanudisastro_loci record of the locus
        trexplorer_locus_id (str): the TRExplorer locus it was matched to

    Returns:
        dict: SignificantCellTypes, MinPvalue, FineMappingMaxPIP, ChromatinAccessibilityCellTypes and Details
    """
    associations = sorted(tanudisastro_locus["associations"].values(),
                          key=lambda a: (a["pvalue"], a["cell_type"], a["egene"]))
    significant_cell_types = list(dict.fromkeys(a["cell_type"] for a in associations))
    pips = [a["pip"] for a in associations if a["pip"] is not None]
    caqtl_min_pvalues = tanudisastro_locus["caqtl_min_pvalue_by_cell_type"]
    caqtl_cell_types = sorted(caqtl_min_pvalues, key=lambda cell_type: (caqtl_min_pvalues[cell_type], cell_type))

    details = {}
    if not is_same_locus(tanudisastro_locus_key, trexplorer_locus_id):
        details["tanudisastroLocusId"] = "-".join(str(x) for x in tanudisastro_locus_key if x is not None)
    details["associations"] = []
    for a in associations:
        association = {
            "cellType": a["cell_type"],
            "eGene": a["egene"],
            "pvalue": float(f"{a['pvalue']:.2e}"),  # 3 significant digits
            "effectSize": round(a["effect_size"], 4),
        }
        if a["pip"] is not None:
            association["fineMappingPIP"] = round(a["pip"], 3)
        details["associations"].append(association)

    return {
        "SignificantCellTypes": ",".join(significant_cell_types),
        "MinPvalue": associations[0]["pvalue"],
        "FineMappingMaxPIP": round(max(pips), 3) if pips else None,
        "ChromatinAccessibilityCellTypes": ",".join(caqtl_cell_types) if caqtl_cell_types else None,
        "Details": details,
    }


def build_trexplorer_lookups(catalog_path):
    """Builds the set of locus IDs and per-chromosome interval trees of the TRExplorer catalog.

    Args:
        catalog_path: Path to TRExplorer catalog JSON

    Returns:
        tuple: (locus_ids set, interval_trees dict)
            - locus_ids: the catalog's locus IDs
            - interval_trees: {chrom: IntervalTree}
    """
    locus_ids = set()
    interval_trees = collections.defaultdict(intervaltree.IntervalTree)

    locus_structure_regex = re.compile(r"^[(]([A-Z]+)[)][+*]")

    print(f"Building lookup structures from TRExplorer catalog: {catalog_path}")
    for record in tqdm.tqdm(get_variant_catalog_iterator(catalog_path, show_progress_bar=False),
                            unit=" records", unit_scale=True):
        locus_id = record["LocusId"]
        reference_region = record["ReferenceRegion"]
        if isinstance(reference_region, list):
            reference_region = reference_region[0]

        chrom, start_0based, end_1based = parse_interval(reference_region)
        chrom = chrom.replace("chr", "")

        # Extract motif from ReferenceMotif field, or parse from LocusStructure
        if "ReferenceMotif" in record:
            motif = record["ReferenceMotif"]
        else:
            match = locus_structure_regex.match(record["LocusStructure"])
            if match:
                motif = match.group(1)
            else:
                # Skip complex locus structures
                continue

        locus_ids.add(locus_id)

        interval_trees[chrom].add(intervaltree.Interval(
            start_0based, end_1based,
            data={
                "locus_id": locus_id,
                "start_0based": start_0based,
                "end_1based": end_1based,
                "canonical_motif": compute_canonical_motif(motif),
                "motif_length": len(motif),
            }
        ))

    print(f"  Found {len(locus_ids):,d} loci")
    print(f"  Built interval trees for {len(interval_trees):,d} chromosomes")

    return locus_ids, interval_trees


def find_trexplorer_match(chrom, start_0based, end_1based, motif, exact_lookup, interval_trees):
    """Find matching TRExplorer locus for a Tanudisastro locus.

    Args:
        chrom: Chromosome (without 'chr' prefix)
        start_0based: 0-based start coordinate
        end_1based: 1-based end coordinate
        motif: Repeat motif
        exact_lookup: The TRExplorer locus IDs (a set, or a dict keyed by locus ID)
        interval_trees: Dict of {chrom: IntervalTree}

    Returns:
        tuple: (trexplorer_locus_id, match_type) or (None, None)
    """
    tanudisastro_locus_id = f"{chrom}-{start_0based}-{end_1based}-{motif}"

    # Exact match first (O(1))
    if tanudisastro_locus_id in exact_lookup:
        return tanudisastro_locus_id, "exact"

    # Fuzzy match via interval tree
    if chrom not in interval_trees:
        return None, None

    candidates = interval_trees[chrom].overlap(start_0based, end_1based)
    tanudisastro_canonical = compute_canonical_motif(motif)

    best_match = None
    best_jaccard = 0

    for interval in candidates:
        data = interval.data
        if not is_matching_motif(motif, tanudisastro_canonical, data):
            continue

        # Jaccard similarity
        intersection = max(0, min(end_1based, data['end_1based']) - max(start_0based, data['start_0based']))
        union = (end_1based - start_0based) + (data['end_1based'] - data['start_0based']) - intersection
        jaccard = intersection / union if union > 0 else 0

        if jaccard >= 0.2 and jaccard > best_jaccard:
            best_jaccard = jaccard
            best_match = data['locus_id']

    if best_match:
        return best_match, "fuzzy"
    return None, None


def is_matching_motif(motif, canonical_motif, trexplorer_interval_data):
    """Returns whether a TRExplorer locus has the same motif as a Tanudisastro locus: the same canonical motif for
    motifs up to 6 bp, and the same motif length for longer ones.

    Args:
        motif: The Tanudisastro locus's motif
        canonical_motif: compute_canonical_motif(motif)
        trexplorer_interval_data: the data of the TRExplorer locus's interval in the interval trees
    """
    if len(motif) <= 6:
        return canonical_motif == trexplorer_interval_data["canonical_motif"]
    return len(motif) == trexplorer_interval_data["motif_length"]


def find_trexplorer_match_by_start_and_motif(chrom, start_0based, motif, interval_trees):
    """Finds the TRExplorer locus for a Table S6 locus, which has no end coordinate: the locus that contains its start
    coordinate, has the same motif (is_matching_motif), and has the closest start.

    Args:
        chrom: Chromosome (without 'chr' prefix)
        start_0based: 0-based start coordinate
        motif: Repeat motif
        interval_trees: Dict of {chrom: IntervalTree}

    Returns:
        tuple: (trexplorer_locus_id, match_type), where match_type is "same start" or "overlapping start", or
            (None, None)
    """
    if chrom not in interval_trees:
        return None, None
    canonical_motif = compute_canonical_motif(motif)
    candidates = sorted(
        (abs(interval.data["start_0based"] - start_0based), interval.data["locus_id"])
        for interval in interval_trees[chrom].at(start_0based)
        if is_matching_motif(motif, canonical_motif, interval.data))
    if not candidates:
        return None, None
    distance, trexplorer_locus_id = candidates[0]
    return trexplorer_locus_id, "same start" if distance == 0 else "overlapping start"


def main():
    args = parse_args()

    locus_ids, interval_trees = build_trexplorer_lookups(args.catalog_path)

    print(f"\nReading supplementary tables from {args.supplementary_tables_dir}")
    s1 = parse_table_s1(read_supplementary_table(args.supplementary_tables_dir, "TableS1"))
    s6 = parse_table_s6(read_supplementary_table(args.supplementary_tables_dir, "TableS6"))
    s4 = parse_table_s4(read_supplementary_table(args.supplementary_tables_dir, "TableS4"))
    print(f"  Read {len(s1):,d} Table S1 rows, {len(s6):,d} Table S6 rows and {len(s4):,d} Table S4 rows")

    tanudisastro_loci = build_tanudisastro_loci(s1, s6, s4)
    print(f"  Found {len(tanudisastro_loci):,d} unique Tanudisastro loci "
          f"({sum(k[2] is None for k in tanudisastro_loci):,d} of them only in Table S6)")

    print("\nMatching Tanudisastro loci to TRExplorer catalog...")
    # Multiple Tanudisastro loci may match the same TRExplorer locus: {trexplorer_locus_id: [tanudisastro_locus_key, ...]}
    trexplorer_candidates = collections.defaultdict(list)
    stats = collections.Counter()
    unmatched_loci = []
    for tanudisastro_locus_key in tqdm.tqdm(sorted(tanudisastro_loci, key=str), unit=" loci"):
        chrom, start_0based, end_1based, motif = tanudisastro_locus_key
        if end_1based is None:
            trexplorer_locus_id, match_type = find_trexplorer_match_by_start_and_motif(
                chrom, start_0based, motif, interval_trees)
        else:
            trexplorer_locus_id, match_type = find_trexplorer_match(
                chrom, start_0based, end_1based, motif, locus_ids, interval_trees)

        if trexplorer_locus_id is None:
            stats["no match"] += 1
            unmatched_loci.append(tanudisastro_locus_key)
            continue

        stats[match_type] += 1
        trexplorer_candidates[trexplorer_locus_id].append(tanudisastro_locus_key)

    output_lookup = {}
    multi_match_count = 0
    for trexplorer_locus_id, candidate_keys in sorted(trexplorer_candidates.items()):
        if len(candidate_keys) > 1:
            multi_match_count += 1
        best_key = pick_tanudisastro_locus_for_trexplorer_locus(candidate_keys, tanudisastro_loci, trexplorer_locus_id)
        output_lookup[trexplorer_locus_id] = build_lookup_entry(best_key, tanudisastro_loci[best_key], trexplorer_locus_id)

    print(f"\nTanudisastro locus matching statistics:")
    print(f"  Total Tanudisastro loci:  {len(tanudisastro_loci):,d}")
    for match_type in ["exact", "fuzzy", "same start", "overlapping start", "no match"]:
        print(f"  {match_type + ':':25s} {stats[match_type]:,d}")
    print(f"  TRExplorer loci with multiple Tanudisastro matches: {multi_match_count:,d} (kept best by min p-value)")
    print(f"  TRExplorer loci with a fine-mapping PIP: "
          f"{sum(e['FineMappingMaxPIP'] is not None for e in output_lookup.values()):,d}")
    print(f"  TRExplorer loci that are caQTLs: "
          f"{sum(e['ChromatinAccessibilityCellTypes'] is not None for e in output_lookup.values()):,d}")

    if unmatched_loci:
        unmatched_loci.sort(key=lambda k: get_min_pvalue(tanudisastro_loci[k]))
        print(f"\nWARNING: {len(unmatched_loci):,d} Tanudisastro loci have no TRExplorer match. The most significant:")
        for key in unmatched_loci[:20]:
            print(f"  {key}: min p-value = {get_min_pvalue(tanudisastro_loci[key]):.2e}")

    print(f"\nWriting output to: {args.output}")
    fopen = gzip.open if args.output.endswith("gz") else open
    with fopen(args.output, "wt") as f:
        json.dump(output_lookup, f, indent=1)

    print(f"Done. Wrote {len(output_lookup):,d} matched loci to {args.output}")


if __name__ == "__main__":
    main()
