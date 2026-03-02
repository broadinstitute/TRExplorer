"""Process PheWAS data from Manigbas et al. 2024 and create a lookup JSON for TRExplorer.

This script reads the supplementary data from Manigbas et al. 2024 and matches each locus
to the TRExplorer catalog using exact matching (by locus ID) or fuzzy matching (by interval
overlap and motif similarity).

Output is a JSON lookup table keyed by TRExplorer locus ID containing:
- AssociatedTraits: comma-separated list of traits, ordered by p-value (most significant first)
- MinPvalue: minimum p-value across all traits (null if no associations)
- Details: per-trait information (traitType, pvalue, log10Pvalue, sampleSize, caviarRank,
           caviarProbability, replicatedInAoU, manigbasLocusId)

For loci that were genotyped but have no associations, output includes:
- AssociatedTraits: "" (empty string)
- MinPvalue: null
- Details: {} (empty dict)
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
    parser = argparse.ArgumentParser(description="Generate Manigbas 2024 lookup JSON")
    parser.add_argument("--excel-file", default="41467_2024_54678_MOESM3_ESM.xlsx",
                        help="Path to Manigbas 2024 supplementary Excel file (default: 41467_2024_54678_MOESM3_ESM.xlsx)")
    parser.add_argument("--trexplorer-catalog", required=True,
                        help="Path to TRExplorer v2 JSON catalog")
    parser.add_argument("-o", "--output", default="Manigbas_2024_lookup.json.gz",
                        help="Output JSON path (default: Manigbas_2024_lookup.json.gz)")
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
    """Find matching TRExplorer locus for a Manigbas locus.

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
    manigbas_locus_id = f"{chrom}-{start_0based}-{end_1based}-{motif}"

    # Exact match first (O(1))
    if manigbas_locus_id in exact_lookup:
        return manigbas_locus_id, "exact"

    # Fuzzy match via interval tree
    if chrom not in interval_trees:
        return None, None

    candidates = interval_trees[chrom].overlap(start_0based, end_1based)
    manigbas_canonical = compute_canonical_motif(motif)
    manigbas_motif_len = len(motif)

    best_match = None
    best_jaccard = 0

    for interval in candidates:
        data = interval.data
        # Motif check: canonical match for 1-6bp, length match for >6bp
        if manigbas_motif_len <= 6:
            if manigbas_canonical != data['canonical_motif']:
                continue
        else:
            if manigbas_motif_len != data['motif_length']:
                continue

        # Jaccard similarity
        intersection = max(0, min(end_1based, data['end_1based']) - max(start_0based, data['start_0based']))
        union = (end_1based - start_0based) + (data['end_1based'] - data['start_0based']) - intersection
        jaccard = intersection / union if union > 0 else 0

        if jaccard >= 0.66 and jaccard > best_jaccard:
            best_jaccard = jaccard
            best_match = data['locus_id']

    if best_match:
        return best_match, "fuzzy"
    return None, None


def main():
    args = parse_args()

    # Resolve excel file path
    script_dir = Path(__file__).parent
    if not Path(args.excel_file).is_absolute():
        excel_path = script_dir / args.excel_file
    else:
        excel_path = Path(args.excel_file)

    if not excel_path.exists():
        sys.exit(f"ERROR: Excel file not found: {excel_path}")

    # Build lookup structures from TRExplorer catalog
    exact_lookup, interval_trees = build_trexplorer_lookups(args.trexplorer_catalog)

    # Read Data S1: All genotyped TRs
    print(f"\nReading Data S1 (genotyped TRs) from: {excel_path}")
    df_s1 = pd.read_excel(excel_path, sheet_name="Data S1_TRs genotyped", header=2)
    print(f"  Read {len(df_s1):,d} rows")

    # Filter to loci used in PheWAS
    df_s1_phewas = df_s1[df_s1['TR used in PheWAS?'] == 1].copy()
    print(f"  Filtered to {len(df_s1_phewas):,d} loci used in PheWAS")

    # Read Data S5: Causal analysis (high-confidence associations)
    print(f"\nReading Data S5 (causal analysis) from: {excel_path}")
    df_s5 = pd.read_excel(excel_path, sheet_name="Data S5_Causal analysis", header=2)
    print(f"  Read {len(df_s5):,d} rows")

    # Filter to high-confidence associations (top ranked by CAVIAR)
    df_s5_high = df_s5[df_s5['TR is top ranked by CAVIAR and significant after conditioning on lead SNV?'] == "Yes"].copy()
    print(f"  Filtered to {len(df_s5_high):,d} high-confidence associations")

    # Read Data S7: All of Us replication
    print(f"\nReading Data S7 (All of Us replication) from: {excel_path}")
    df_s7 = pd.read_excel(excel_path, sheet_name="Data S7_All of Us replication", header=2)
    print(f"  Read {len(df_s7):,d} rows")

    # Build replication lookup from Data S7
    # Key: (chrom, start, end, normalized_trait_name)
    replication_lookup = {}
    nan_pvalue_count = 0
    for _, row in df_s7.iterrows():
        chrom = str(row['Chr']).replace('chr', '')
        start = int(row['TR start, hg38'])
        end = int(row['TR end, hg38'])
        trait = row['UK Biobank, trait description'].strip().lower()
        pvalue = row['All of Us, meta-analysis p-value']

        if pd.isna(pvalue):
            nan_pvalue_count += 1
            replicated = False
        else:
            replicated = pvalue <= 0.05

        key = (chrom, start, end, trait)
        replication_lookup[key] = replicated

    if nan_pvalue_count > 0:
        print(f"  WARNING: {nan_pvalue_count:,d} rows in Data S7 have NaN p-values")

    print(f"  Built replication lookup with {len(replication_lookup):,d} entries")

    # Strip 'chr' prefix from chromosomes in Data S1 and S5
    df_s1_phewas['Chr'] = df_s1_phewas['Chr'].astype(str).str.replace('chr', '', regex=False)
    df_s5_high['Chr'] = df_s5_high['Chr'].astype(str).str.replace('chr', '', regex=False)

    # Build Manigbas locus IDs for all PheWAS loci
    print("\nMatching Manigbas loci to TRExplorer catalog...")

    # First, create locus ID for all PheWAS loci from Data S1
    phewas_loci = set()
    for _, row in df_s1_phewas.iterrows():
        chrom = str(row['Chr'])
        start = int(row['TR start, hg38'])
        end = int(row['TR end, hg38'])
        motif = row['TR motif']
        phewas_loci.add((chrom, start, end, motif))

    print(f"  Found {len(phewas_loci):,d} unique PheWAS loci from Data S1")

    # Group associations by Manigbas locus from Data S5
    manigbas_associations = collections.defaultdict(list)
    for _, row in df_s5_high.iterrows():
        chrom = str(row['Chr'])
        start = int(row['TR start, hg38'])
        end = int(row['TR end, hg38'])
        motif = row['TR motif']
        locus_key = (chrom, start, end, motif)

        # Build association record
        trait = row['Associated trait']
        trait_normalized = trait.strip().lower()
        repl_key = (chrom, start, end, trait_normalized)
        replicated = replication_lookup.get(repl_key, False)

        association = {
            "trait": trait,
            "traitType": row['Associated trait type'],
            "pvalue": row['METAL nominal p-value'],
            "log10Pvalue": row['METAL -log10 p-value'],
            "sampleSize": int(row['Total sample size']),
            "caviarRank": int(row['CAVIAR Rank of TR']),
            "caviarProbability": row['Probability of TR as being causal variant from CAVIAR'],
            "replicatedInAoU": replicated,
        }
        manigbas_associations[locus_key].append(association)

    print(f"  Found {len(manigbas_associations):,d} Manigbas loci with associations in Data S5")

    # Match all PheWAS loci to TRExplorer
    # Track which TRExplorer loci have associations vs negative results
    trexplorer_with_associations = collections.defaultdict(list)  # {trexplorer_id: [(manigbas_key, match_type, associations), ...]}
    trexplorer_negative = collections.defaultdict(list)  # {trexplorer_id: [(manigbas_key, match_type), ...]}
    stats = {"exact": 0, "fuzzy": 0, "no_match": 0}
    unmatched_with_associations = []

    for locus_key in tqdm.tqdm(phewas_loci, unit=" loci"):
        chrom, start, end, motif = locus_key
        manigbas_locus_id = f"{chrom}-{start}-{end}-{motif}"

        # Find TRExplorer match
        trexplorer_locus_id, match_type = find_trexplorer_match(
            chrom, start, end, motif, exact_lookup, interval_trees
        )

        if trexplorer_locus_id is None:
            stats["no_match"] += 1
            if locus_key in manigbas_associations:
                unmatched_with_associations.append((manigbas_locus_id, manigbas_associations[locus_key]))
            continue

        stats[match_type] += 1

        if locus_key in manigbas_associations:
            trexplorer_with_associations[trexplorer_locus_id].append(
                (manigbas_locus_id, match_type, manigbas_associations[locus_key])
            )
        else:
            trexplorer_negative[trexplorer_locus_id].append(
                (manigbas_locus_id, match_type)
            )

    # Build output lookup
    output_lookup = {}
    multi_match_association_count = 0
    multi_match_negative_count = 0

    # Process loci WITH associations
    for trexplorer_locus_id, candidates in trexplorer_with_associations.items():
        if len(candidates) > 1:
            multi_match_association_count += 1
            print(f"  WARNING: {len(candidates)} Manigbas loci fuzzy-match to TRExplorer locus {trexplorer_locus_id} - merging associations")

        # Merge all associations from all matching Manigbas loci
        all_associations = []
        manigbas_ids = []
        for manigbas_locus_id, match_type, associations in candidates:
            manigbas_ids.append(manigbas_locus_id)
            for assoc in associations:
                assoc_copy = assoc.copy()
                assoc_copy["manigbasLocusId"] = manigbas_locus_id
                all_associations.append(assoc_copy)

        # Sort by p-value (most significant first)
        all_associations.sort(key=lambda x: x["pvalue"])

        # Build Details dict and ordered traits list
        details = {}
        ordered_traits = []
        for assoc in all_associations:
            trait = assoc["trait"]
            if trait not in details:
                ordered_traits.append(trait)
            details[trait] = {
                "traitType": assoc["traitType"],
                "pvalue": assoc["pvalue"],
                "log10Pvalue": assoc["log10Pvalue"],
                "sampleSize": assoc["sampleSize"],
                "caviarRank": assoc["caviarRank"],
                "caviarProbability": assoc["caviarProbability"],
                "replicatedInAoU": assoc["replicatedInAoU"],
                "manigbasLocusId": assoc["manigbasLocusId"],
            }

        output_lookup[trexplorer_locus_id] = {
            "AssociatedTraits": ",".join(ordered_traits),
            "MinPvalue": min(a["pvalue"] for a in details.values()),
            "Details": details,
        }

    # Process loci WITHOUT associations (negative results)
    for trexplorer_locus_id, candidates in trexplorer_negative.items():
        # Skip if this locus already has associations
        if trexplorer_locus_id in output_lookup:
            continue

        if len(candidates) > 1:
            multi_match_negative_count += 1

        output_lookup[trexplorer_locus_id] = {
            "AssociatedTraits": "",
            "MinPvalue": None,
            "Details": {},
        }

    # Print statistics
    total = sum(stats.values())
    print(f"\nManigbas 2024 locus matching statistics:")
    print(f"  Total PheWAS loci:        {total:,d}")
    print(f"  Exact matches:            {stats['exact']:,d}")
    print(f"  Fuzzy matches (J>=0.66):  {stats['fuzzy']:,d}")
    print(f"  No match found:           {stats['no_match']:,d}")
    print(f"  TRExplorer loci with multiple Manigbas matches (with associations): {multi_match_association_count:,d}")
    print(f"  TRExplorer loci with multiple Manigbas matches (negative): {multi_match_negative_count:,d}")

    # Count output categories
    loci_with_associations = sum(1 for v in output_lookup.values() if v["AssociatedTraits"])
    loci_negative = sum(1 for v in output_lookup.values() if not v["AssociatedTraits"])
    total_traits = sum(len(v["Details"]) for v in output_lookup.values())
    replicated_count = sum(
        1 for v in output_lookup.values()
        for d in v["Details"].values()
        if d.get("replicatedInAoU")
    )

    print(f"\nOutput statistics:")
    print(f"  TRExplorer loci with associations:    {loci_with_associations:,d}")
    print(f"  TRExplorer loci without associations: {loci_negative:,d}")
    print(f"  Total trait associations:             {total_traits:,d}")
    print(f"  Associations replicated in AoU:       {replicated_count:,d}")

    # Print warnings for unmatched loci with associations
    if unmatched_with_associations:
        print(f"\nWARNING: {len(unmatched_with_associations):,d} Manigbas loci with associations could not be matched to TRExplorer:")
        for manigbas_id, associations in unmatched_with_associations[:20]:
            traits = ", ".join(a["trait"] for a in associations)
            min_p = min(a["pvalue"] for a in associations)
            print(f"  {manigbas_id}: {len(associations)} trait(s) [{traits}], min p={min_p:.2e}")
        if len(unmatched_with_associations) > 20:
            print(f"  ... and {len(unmatched_with_associations) - 20:,d} more")

    # Write output JSON
    print(f"\nWriting output to: {args.output}")
    fopen = gzip.open if args.output.endswith("gz") else open
    with fopen(args.output, "wt") as f:
        json.dump(output_lookup, f, indent=2)

    print(f"Done. Wrote {len(output_lookup):,d} loci to {args.output}")


if __name__ == "__main__":
    main()
