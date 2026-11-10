"""Add the PseudoautosomalRegion column to an existing catalog table and fill it from an annotated catalog.

The regular way to get this column is to rebuild the catalog table with `load_bigquery_main_table.py` from an
annotated catalog that has the field. This script fills it in the live table instead, without a rebuild
(bigquery_loader_utils.fill_catalog_table_columns):

  1. adds the PseudoautosomalRegion column if the table does not have it yet
  2. loads the LocusId and PseudoautosomalRegion of every annotated catalog locus that has a value into a temporary
     staging table
  3. in one MERGE statement, copies the staging table's values into the rows with the same LocusId and sets the
     column to NULL in every other row that has a value (so loci without a value end up NULL, even after an earlier
     fill); being one statement, it either applies fully or leaves the table unchanged
  4. deletes the staging table

Run it on the table the website queries before publishing website pages that use this column: the locus page selects
it for every locus, and the ExpansionHunter export dialog counts the loci that have a value, so their queries fail on
a table that doesn't have the column.

The annotated catalog is the TRExplorer.repeat_catalog_v2.1.hg38.1_to_1000bp_motifs.EH.with_annotations.json.gz
written by tandem-repeat-catalogs/run_all_steps_to_generate_the_hg38_catalog.py, whose step 16 adds the field (with
tandem-repeat-catalogs/scripts/add_pseudoautosomal_region_annotations.py, which can also add it to an existing
annotated catalog).

With --check-only, nothing is changed: the script only reports how many of the LocusIds that have a value are in the
table, and how many rows of the table have a value now.

Usage:
    python3 add_pseudoautosomal_region_column_to_catalog_table.py <annotated catalog .json.gz> --check-only
    python3 add_pseudoautosomal_region_column_to_catalog_table.py <annotated catalog .json.gz>   # the table the website queries
    python3 add_pseudoautosomal_region_column_to_catalog_table.py <annotated catalog .json.gz> --table-id catalog_20261030_170523
"""

from bigquery_loader_utils import run_script_that_fills_catalog_table_column_from_annotated_catalog_field

COLUMN_NAME = "PseudoautosomalRegion"


def main():
    run_script_that_fills_catalog_table_column_from_annotated_catalog_field(
        __doc__, COLUMN_NAME, "pseudoautosomal_region_staging")


if __name__ == "__main__":
    main()
