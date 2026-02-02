"""Process sc-eTR data from Tanudisastro et al. 2024 and match to TRExplorer loci.

This script reads the source data from TableS1v0.1.tsv and matches each locus to the
TRExplorer catalog using exact matching (by locus ID) or fuzzy matching (by interval
overlap and motif similarity).

Output is a JSON lookup table keyed by TRExplorer locus ID containing:
- Tanudisastro2024_LocusId: the original locus ID from the source data
- MatchType: "exact" or "fuzzy"
- SignificantCellTypes: comma-separated list of cell types, ordered by p-value
- MinPvalue: minimum p-value across all cell types
- Details: per-cell-type information (pvalue, effectSize, eGene, rank, totalInCellType)
"""

import argparse
import collections
import gzip
import json
from pathlib import Path
import re
import sys

import intervaltree
import pandas as pd
import tqdm

from str_analysis.utils.canonical_repeat_unit import compute_canonical_motif
from str_analysis.utils.eh_catalog_utils import get_variant_catalog_iterator
from str_analysis.utils.misc_utils import parse_interval


def parse_args():
    parser = argparse.ArgumentParser(description="Generate Tanudisastro 2024 lookup JSON")
    parser.add_argument("--catalog-path", required=True,
                        help="Path to TRExplorer catalog JSON (e.g., TRExplorer.repeat_catalog_v2.hg38.1_to_1000bp_motifs.EH.json.gz)")
    parser.add_argument("--output", default="Tanudisastro_2024_lookup.json.gz",
                        help="Output JSON path (default: Tanudisastro_2024_lookup.json.gz)")
    return parser.parse_args()


def build_trexplorer_lookups(catalog_path):
    """Build exact lookup dict and interval trees from TRExplorer catalog.

    Args:
        catalog_path: Path to TRExplorer catalog JSON

    Returns:
        tuple: (exact_lookup dict, interval_trees dict)
            - exact_lookup: {locus_id: record}
            - interval_trees: {chrom: IntervalTree}
    """
    exact_lookup = {}
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

        canonical_motif = compute_canonical_motif(motif)

        exact_lookup[locus_id] = record

        interval_trees[chrom].add(intervaltree.Interval(
            start_0based, end_1based,
            data={
                "locus_id": locus_id,
                "start_0based": start_0based,
                "end_1based": end_1based,
                "canonical_motif": canonical_motif,
                "motif_length": len(motif),
            }
        ))

    print(f"  Built exact lookup with {len(exact_lookup):,d} loci")
    print(f"  Built interval trees for {len(interval_trees):,d} chromosomes")

    return exact_lookup, interval_trees


def find_trexplorer_match(chrom, start_0based, end_1based, motif, exact_lookup, interval_trees):
    """Find matching TRExplorer locus for a Tanudisastro locus.

    Args:
        chrom: Chromosome (without 'chr' prefix)
        start_0based: 0-based start coordinate
        end_1based: 1-based end coordinate
        motif: Repeat motif
        exact_lookup: Dict of {locus_id: record}
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
    tanudisastro_motif_len = len(motif)

    best_match = None
    best_jaccard = 0

    for interval in candidates:
        data = interval.data
        # Motif check: canonical match for 1-6bp, length match for >6bp
        if tanudisastro_motif_len <= 6:
            if tanudisastro_canonical != data['canonical_motif']:
                continue
        else:
            if tanudisastro_motif_len != data['motif_length']:
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


def main():
    args = parse_args()

    # Build lookup structures from TRExplorer catalog
    exact_lookup, interval_trees = build_trexplorer_lookups(args.catalog_path)

    # Read and process TableS1v0.1.tsv
    script_dir = Path(__file__).parent
    source_tsv = script_dir / "TableS1v0.1.tsv"
    print(f"\nReading source data: {source_tsv}")
    df = pd.read_table(source_tsv)
    print(f"  Read {len(df):,d} rows")

    # Deduplicate identical rows
    before_dedup = len(df)
    df = df.drop_duplicates()
    print(f"  After deduplication: {len(df):,d} rows ({before_dedup - len(df):,d} duplicates removed)")

    # Strip 'chr' prefix
    df['chromosome'] = df['chromosome'].str.replace('chr', '', regex=False)

    # Compute p-value rank within each cell type
    df['rank'] = df.groupby('cell_type')['nominal p-value'].rank(method='min', ascending=True).astype(int)
    df['total_in_cell_type'] = df.groupby('cell_type')['nominal p-value'].transform('count')

    # Build output lookup
    print("\nMatching Tanudisastro loci to TRExplorer catalog...")

    # Group source data by Tanudisastro locus ID
    tanudisastro_loci = collections.defaultdict(list)
    for _, row in df.iterrows():
        locus_id = f"{row['chromosome']}-{row['start coordinate (hg38)']}-{row['end coordinate (hg38)']}-{row['motif']}"
        tanudisastro_loci[locus_id].append(row)

    print(f"  Found {len(tanudisastro_loci):,d} unique Tanudisastro loci")

    # Match each locus to TRExplorer and collect candidates
    # (multiple Tanudisastro loci may map to the same TRExplorer locus)
    trexplorer_candidates = collections.defaultdict(list)  # {trexplorer_id: [(tanudisastro_id, match_type, rows, min_pvalue), ...]}
    stats = {"exact": 0, "fuzzy": 0, "no_match": 0}
    unmatched_significant = []  # Track unmatched loci with significant p-values

    for tanudisastro_locus_id, rows in tqdm.tqdm(tanudisastro_loci.items(), unit=" loci"):
        # Parse locus ID
        parts = tanudisastro_locus_id.split("-")
        chrom = parts[0]
        start_0based = int(parts[1])
        end_1based = int(parts[2])
        motif = parts[3]

        # Find TRExplorer match
        trexplorer_locus_id, match_type = find_trexplorer_match(
            chrom, start_0based, end_1based, motif, exact_lookup, interval_trees
        )

        min_pvalue = min(row['nominal p-value'] for row in rows)

        if trexplorer_locus_id is None:
            stats["no_match"] += 1
            if min_pvalue < 0.05:
                unmatched_significant.append((tanudisastro_locus_id, min_pvalue))
            continue

        stats[match_type] += 1
        trexplorer_candidates[trexplorer_locus_id].append(
            (tanudisastro_locus_id, match_type, rows, min_pvalue)
        )

    # Build output lookup, keeping only the Tanudisastro locus with smallest min p-value
    # for each TRExplorer locus
    output_lookup = {}
    multi_match_count = 0

    for trexplorer_locus_id, candidates in trexplorer_candidates.items():
        # Sort by min p-value and keep the best one
        candidates.sort(key=lambda x: x[3])
        if len(candidates) > 1:
            multi_match_count += 1

        tanudisastro_locus_id, match_type, rows, min_pvalue = candidates[0]

        # Build details dict, sorted by p-value
        sorted_rows = sorted(rows, key=lambda r: r['nominal p-value'])
        details = {}
        for row in sorted_rows:
            cell_type = row['cell_type']
            details[cell_type] = {
                "pvalue": row['nominal p-value'],
                "effectSize": row['effect size'],
                "eGene": row['eGene'],
                "rank": int(row['rank']),
                "totalInCellType": int(row['total_in_cell_type']),
            }

        output_lookup[trexplorer_locus_id] = {
            "Tanudisastro2024_LocusId": tanudisastro_locus_id,
            "MatchType": match_type,
            "SignificantCellTypes": ",".join(details.keys()),
            "MinPvalue": min_pvalue,
            "Details": details,
        }

    # Print statistics
    total = sum(stats.values())
    print(f"\nTanudisastro 2024 locus matching statistics:")
    print(f"  Total Tanudisastro loci:  {total:,d}")
    print(f"  Exact matches:            {stats['exact']:,d}")
    print(f"  Fuzzy matches (J>=0.2):   {stats['fuzzy']:,d}")
    print(f"  No match found:           {stats['no_match']:,d}")
    print(f"  TRExplorer loci with multiple Tanudisastro matches: {multi_match_count:,d} (kept best by min p-value)")

    # Print errors for unmatched significant loci
    if unmatched_significant:
        unmatched_significant.sort(key=lambda x: x[1])  # Sort by p-value
        print(f"\nWARNING: {len(unmatched_significant):,d} unmatched loci with p < 0.05:")
        for locus_id, min_pvalue in unmatched_significant[:50]:
            print(f"  {locus_id}: min p-value = {min_pvalue:.2e}")
        if len(unmatched_significant) > 50:
            print(f"  ... and {len(unmatched_significant) - 50:,d} more")

    # Write output JSON
    print(f"\nWriting output to: {args.output}")
    fopen = gzip.open if args.output.endswith("gz") else open
    with fopen(args.output, "wt") as f:
        json.dump(output_lookup, f, indent=2)

    print(f"Done. Wrote {len(output_lookup):,d} matched loci to {args.output}")


if __name__ == "__main__":
    main()
