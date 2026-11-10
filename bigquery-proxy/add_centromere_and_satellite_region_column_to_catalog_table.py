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
written by tandem-repeat-catalogs/run_all_steps_to_generate_the_hg38_catalog.py, whose step 16 adds the field.

With --check-only, nothing is changed: the script only reports how many of the first CHECK_ONLY_MAX_LOCI LocusIds
that have a value are in the table, and how many rows of the table have each class now.

Usage:
    python3 add_centromere_and_satellite_region_column_to_catalog_table.py <annotated catalog .json.gz> --check-only
    python3 add_centromere_and_satellite_region_column_to_catalog_table.py <annotated catalog .json.gz>   # the table the website queries
    python3 add_centromere_and_satellite_region_column_to_catalog_table.py <annotated catalog .json.gz> --table-id catalog_20261018_153653
"""

import collections

from bigquery_loader_utils import run_script_that_fills_catalog_table_column_from_annotated_catalog_field

CHECK_ONLY_MAX_LOCI = 10000
COLUMN_NAME = "CentromereAndSatelliteRegion"


def count_staging_rows_by_region_class(rows):
    """Returns {CentromereAndSatelliteRegion class: number of staging rows with that class}."""
    return dict(collections.Counter(row[COLUMN_NAME] for row in rows))


def count_catalog_table_rows_by_region_class(client, catalog_table):
    """Returns {CentromereAndSatelliteRegion class: number of rows of the catalog table with that class}."""
    return {row["value"]: row["n"] for row in client.query(
        f"SELECT {COLUMN_NAME} AS value, COUNT(*) AS n FROM `{catalog_table}` "
        f"WHERE {COLUMN_NAME} IS NOT NULL GROUP BY 1").result()}


def print_counts_by_region_class(title, counts):
    """Prints {CentromereAndSatelliteRegion class: number of rows}, the most common class first."""
    print(f"{title}:")
    for value, n in sorted(counts.items(), key=lambda item: -item[1]):
        print(f"  {n:10,d}  {value}")


def main():
    run_script_that_fills_catalog_table_column_from_annotated_catalog_field(
        __doc__, COLUMN_NAME, "centromere_and_satellite_region_staging",
        count_values_in_staging_rows=count_staging_rows_by_region_class,
        count_values_in_catalog_table=count_catalog_table_rows_by_region_class,
        print_counts=print_counts_by_region_class, max_locus_ids_to_check=CHECK_ONLY_MAX_LOCI)


if __name__ == "__main__":
    main()
