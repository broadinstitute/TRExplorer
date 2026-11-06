"""Add or update the HPRC256_ReferenceLPSRepeatCount and HPRC256_MethylationInterval columns of an existing catalog
table.

  HPRC256_ReferenceLPSRepeatCount: the longest pure segment (LPS) of each locus's motif in the hg38 reference sequence
      of its HPRC256_NarrowestInterval, which is where the HPRC256:LPS allele histograms and genotype heatmaps mark the
      reference allele
  HPRC256_MethylationInterval: the TRGT interval that each locus's HPRC256_AlleleSizeAndMethylationDistribution comes
      from, which the methylation heatmap's Population and Sex queries are pinned to

The regular way to get these columns is to rebuild the catalog table with `load_bigquery_main_table.py`, which computes
them the same way. This script fills them in the live table instead, without a rebuild:

  1. reads every locus's HPRC256_NarrowestInterval and HPRC256_ModeAllele from the table
  2. computes HPRC256_ReferenceLPSRepeatCount from the reference genome
     (bigquery_loader_utils.compute_reference_lps_repeat_count_in_interval)
  3. picks each locus's methylation interval from the HPRC256 methylation TSV, the way the loader picks the interval of
     HPRC256_AlleleSizeAndMethylationDistribution (bigquery_loader_utils.parse_all_samples_distribution_column)
  4. fills both columns with bigquery_loader_utils.fill_catalog_table_columns: one MERGE statement that either applies
     fully or leaves the table unchanged

The website's charts read these columns, so run this script on the table the website queries before publishing a
regenerated index.html and locus.html.

With --check-only, nothing is changed: the script reports the counts it would write, and how often the reference LPS
equals the HPRC256 mode allele, as a check that it is measured the same way as the alleles.

Usage:
    python3 add_or_update_hprc256_reference_lps_and_methylation_interval_columns_in_catalog_table.py --check-only
    python3 add_or_update_hprc256_reference_lps_and_methylation_interval_columns_in_catalog_table.py
"""

import argparse
import collections
import os

import pyfaidx
import tqdm
from google.cloud import bigquery

from bigquery_loader_utils import DATASET_ID, PROJECT_ID, compute_reference_lps_repeat_count_in_interval, \
    count_catalog_table_rows_with_a_value_in_each_column, count_staging_rows_with_a_value_in_each_column, \
    fill_catalog_table_columns, parse_all_samples_distribution_column, read_website_table_id
from global_constants import HPRC256_REFERENCE_LPS_AND_METHYLATION_INTERVAL_COLUMN_NAMES

DEFAULT_HPRC256_METHYLATION_STRATIFIED_TSV = \
    "../data-prep/hprc-lps/trgt-hprc.methylation.stratified.by_population.by_sex.256_samples.tsv.gz"


def read_hprc256_loci_from_catalog_table(client, catalog_table):
    """Returns the LocusId, HPRC256_NarrowestInterval and HPRC256_ModeAllele of every catalog table row that has an
    HPRC256_NarrowestInterval, sorted by chromosome."""
    query = (f"SELECT LocusId, chrom, HPRC256_NarrowestInterval, HPRC256_ModeAllele FROM `{catalog_table}` "
             f"WHERE HPRC256_NarrowestInterval IS NOT NULL ORDER BY chrom_index")
    return [dict(row) for row in client.query(query).result()]


def compute_reference_lps_by_locus(hprc256_loci, reference_fasta):
    """Returns {LocusId: HPRC256_ReferenceLPSRepeatCount} for the loci from read_hprc256_loci_from_catalog_table, which
    must be sorted by chromosome."""
    fasta = pyfaidx.Fasta(os.path.expanduser(reference_fasta), one_based_attributes=False, as_raw=True)
    reference_lps_by_locus = {}
    chrom = chrom_sequence = None
    for locus in tqdm.tqdm(hprc256_loci, unit=" loci", unit_scale=True):
        if locus["chrom"] != chrom:
            chrom = locus["chrom"]
            # The table's chrom has no "chr" prefix, unlike the reference genome's chromosome names
            chrom_sequence = str(fasta[f"chr{chrom}"])
        reference_lps_by_locus[locus["LocusId"]] = compute_reference_lps_repeat_count_in_interval(
            chrom_sequence, locus["HPRC256_NarrowestInterval"], locus["LocusId"].split("-")[-1])
    return reference_lps_by_locus


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--reference-fasta", default="~/hg38.fa", help="hg38 reference genome FASTA, indexed")
    parser.add_argument("--hprc256-methylation-stratified-tsv", default=DEFAULT_HPRC256_METHYLATION_STRATIFIED_TSV,
                        help="the HPRC256 methylation TSV that the catalog table was loaded from")
    parser.add_argument("--table-id", help="catalog table to fill (default: the one the website queries, i.e. the "
                                           "TABLE_ID in ../website/header_template.html)")
    parser.add_argument("--check-only", action="store_true",
                        help="Only report the counts that would be written and the reference LPS check")
    args = parser.parse_args()

    client = bigquery.Client(project=PROJECT_ID)
    table_id = args.table_id or read_website_table_id()
    catalog_table = f"{PROJECT_ID}.{DATASET_ID}.{table_id}"
    print(f"Catalog table: {catalog_table}")

    hprc256_loci = read_hprc256_loci_from_catalog_table(client, catalog_table)
    print(f"Read {len(hprc256_loci):,} loci with an HPRC256_NarrowestInterval from {table_id}")
    reference_lps_by_locus = compute_reference_lps_by_locus(hprc256_loci, args.reference_fasta)

    methylation_lookup = parse_all_samples_distribution_column(
        args.hprc256_methylation_stratified_tsv, "AlleleSizeAndMethylationDistribution",
        narrowest_interval_by_locus={locus["LocusId"]: locus["HPRC256_NarrowestInterval"] for locus in hprc256_loci})

    rows = [{
        "LocusId": locus_id,
        "HPRC256_ReferenceLPSRepeatCount": reference_lps_by_locus.get(locus_id),
        "HPRC256_MethylationInterval": methylation_lookup[locus_id]["interval"] if locus_id in methylation_lookup else None,
    } for locus_id in sorted(set(reference_lps_by_locus) | set(methylation_lookup))]
    staging_counts = count_staging_rows_with_a_value_in_each_column(
        rows, HPRC256_REFERENCE_LPS_AND_METHYLATION_INTERVAL_COLUMN_NAMES)
    print(f"Counts to write: {staging_counts}")

    # Where most samples carry the reference allele, the mode allele should equal the reference LPS
    num_loci_where_reference_lps_differs_from_mode = collections.Counter(
        locus["HPRC256_ModeAllele"] - reference_lps_by_locus[locus["LocusId"]] for locus in hprc256_loci
        if locus["HPRC256_ModeAllele"] is not None)
    num_loci_with_a_mode = sum(num_loci_where_reference_lps_differs_from_mode.values())
    print(f"Mode allele - reference LPS, for {num_loci_with_a_mode:,} loci:")
    for difference, n in sorted(num_loci_where_reference_lps_differs_from_mode.items(), key=lambda x: -x[1])[:10]:
        print(f"  {difference:+d}: {n:10,d} ({n / num_loci_with_a_mode:.1%})")

    if args.check_only:
        print("Nothing was changed")
        return

    fill_catalog_table_columns(client, catalog_table, rows, HPRC256_REFERENCE_LPS_AND_METHYLATION_INTERVAL_COLUMN_NAMES,
                               "hprc256_reference_lps_and_methylation_interval_staging")
    table_counts = count_catalog_table_rows_with_a_value_in_each_column(
        client, catalog_table, HPRC256_REFERENCE_LPS_AND_METHYLATION_INTERVAL_COLUMN_NAMES)
    print(f"Counts in {table_id} after the update: {table_counts}")
    if table_counts != staging_counts:
        print("WARNING: these differ from the counts that were written")


if __name__ == "__main__":
    main()
