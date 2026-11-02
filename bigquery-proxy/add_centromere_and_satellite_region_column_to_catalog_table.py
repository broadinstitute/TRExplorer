"""Add the CentromereAndSatelliteRegion column to an existing catalog table and fill it from an annotated catalog.

The regular way to get this column is to rebuild the catalog table with `load_bigquery_main_table.py` from an
annotated catalog that has the field. This script fills it in the live table instead, without a rebuild
(bigquery_loader_utils.fill_catalog_table_columns):

  1. adds the CentromereAndSatelliteRegion column if the table does not have it yet
  2. loads the LocusId and CentromereAndSatelliteRegion of every annotated catalog locus that has a value into a
     temporary staging table
  3. in one MERGE statement, copies the staging table's values into the rows with the same LocusId and sets the
     column to NULL in every other row that has a value (so loci without a value end up NULL, even after an earlier
     fill); being one statement, it either applies fully or leaves the table unchanged
  4. deletes the staging table

The annotated catalog is the TRExplorer.repeat_catalog_v2.1.hg38.1_to_1000bp_motifs.EH.with_annotations.json.gz
written by tandem-repeat-catalogs/run_all_steps_to_generate_the_hg38_catalog.py, whose step 16 adds the field. After
filling the column, add allowCustomFilter, allowExport and a group to its entry in global_constants.py and regenerate
the website (`cd ../website && python3 generate_website.py`) so its filters and export dialog list it.

With --check-only, nothing is changed: the script only reports how many of the first CHECK_ONLY_MAX_LOCI LocusIds
that have a value are in the table.

Usage:
    python3 add_centromere_and_satellite_region_column_to_catalog_table.py <annotated catalog .json.gz> --check-only
    python3 add_centromere_and_satellite_region_column_to_catalog_table.py <annotated catalog .json.gz>   # the table the website queries
    python3 add_centromere_and_satellite_region_column_to_catalog_table.py <annotated catalog .json.gz> --table-id catalog_20261018_153653
"""

import argparse
import collections
import gzip
import os

import ijson
from google.cloud import bigquery

from bigquery_loader_utils import count_locus_ids_in_catalog_table, fill_catalog_table_columns, read_website_table_id

PROJECT_ID = "cmg-analysis"
DATASET_ID = "tandem_repeat_explorer"
CHECK_ONLY_MAX_LOCI = 10000
COLUMN_NAME = "CentromereAndSatelliteRegion"


def read_locus_ids_and_values(catalog_json_path):
    """Returns the LocusId and CentromereAndSatelliteRegion of each annotated catalog locus that has a value.

    Args:
        catalog_json_path (str): path of the annotated catalog JSON file, which can be gzipped

    Returns:
        list: (LocusId, CentromereAndSatelliteRegion) tuples, in catalog order
    """
    rows = []
    num_records = 0
    num_records_with_field = 0
    with (gzip.open if catalog_json_path.endswith("gz") else open)(catalog_json_path, "rb") as f:
        for record in ijson.items(f, "item"):
            num_records += 1
            if COLUMN_NAME not in record:
                continue
            num_records_with_field += 1
            if record[COLUMN_NAME]:
                rows.append((record["LocusId"], record[COLUMN_NAME]))

    if num_records_with_field < num_records:
        raise SystemExit(f"{num_records - num_records_with_field:,} of the {num_records:,} records in "
                         f"{catalog_json_path} have no {COLUMN_NAME} field. Pass an annotated catalog built after "
                         f"the field was added.")

    duplicate_locus_ids = [locus_id for locus_id, count in collections.Counter(
        locus_id for locus_id, _ in rows).items() if count > 1]
    if duplicate_locus_ids:
        raise SystemExit(f"{len(duplicate_locus_ids):,} LocusIds appear more than once in {catalog_json_path}, "
                         f"eg. {duplicate_locus_ids[:3]}")

    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("catalog_json", help="Annotated catalog JSON file that has the CentromereAndSatelliteRegion field")
    parser.add_argument("--table-id", help="catalog table to fill (default: the one the website queries, i.e. the "
                                           "TABLE_ID in ../website/header_template.html)")
    parser.add_argument("--check-only", action="store_true",
                        help="Only report how many of the catalog's LocusIds that have a value are in the table")
    args = parser.parse_args()

    rows = read_locus_ids_and_values(os.path.expanduser(args.catalog_json))
    print(f"Read {len(rows):,} loci that have a {COLUMN_NAME} value from {args.catalog_json}")

    client = bigquery.Client(project=PROJECT_ID)
    table_id = args.table_id or read_website_table_id()
    catalog_table = f"{PROJECT_ID}.{DATASET_ID}.{table_id}"
    print(f"Catalog table: {catalog_table}")

    if args.check_only:
        locus_ids = [locus_id for locus_id, _ in rows[:CHECK_ONLY_MAX_LOCI]]
        n_found, gb_scanned = count_locus_ids_in_catalog_table(client, catalog_table, locus_ids)
        print(f"{n_found:,} of the first {len(locus_ids):,} LocusIds that have a value are in {table_id} "
              f"({gb_scanned:.2f} GB scanned); nothing was changed")
        return

    fill_catalog_table_columns(client, catalog_table, [{"LocusId": locus_id, COLUMN_NAME: value} for locus_id, value in rows],
                               [COLUMN_NAME], "centromere_and_satellite_region_staging")
    counts = client.query(f"SELECT {COLUMN_NAME} AS value, COUNT(*) AS n FROM `{catalog_table}` "
                          f"WHERE {COLUMN_NAME} IS NOT NULL GROUP BY 1 ORDER BY 2 DESC").result()
    counts = [(row["value"], row["n"]) for row in counts]
    print(f"Filled {COLUMN_NAME} in {sum(n for _, n in counts):,} rows of {table_id} and set it to NULL in the "
          f"rest (the catalog has {len(rows):,} loci with a value; a difference means some LocusIds are not in "
          f"the table):")
    for value, n in counts:
        print(f"  {n:10,d}  {value}")
    print("Next: add allowCustomFilter, allowExport and a group to the column's entry in global_constants.py, then "
          "regenerate the website so its filters and export dialog list it: cd ../website && python3 generate_website.py")


if __name__ == "__main__":
    main()
