"""Add or update the Tanudisastro2024_* columns of an existing catalog table, filling them from the Tanudisastro lookup
JSON (the sc-eTRs of Tanudisastro et al. 2026, Nature).

The regular way to get these columns is to rebuild the catalog table with `load_bigquery_main_table.py`, which reads
the same lookup (--tanudisastro-2024-lookup-json). This script fills them in the live table instead, without a rebuild
(bigquery_loader_utils.add_or_update_catalog_table_columns_from_lookup_json):

  1. adds any Tanudisastro2024_* columns (global_constants.TANUDISASTRO_2024_COLUMN_NAMES) the table does not have yet
  2. loads the lookup into a temporary staging table, one row per LocusId, with the same values
     load_bigquery_main_table.py writes (bigquery_loader_utils.get_tanudisastro_2024_column_values)
  3. in one MERGE statement, copies the staging table's values into the rows with the same LocusId and sets the
     Tanudisastro2024_* columns to NULL in every other row that has a value (so loci missing from the lookup end up
     NULL); being one statement, it either applies fully or leaves the table unchanged
  4. deletes the staging table

The lookup is written by ../data-prep/literature_review/Tanudisastro_2024/generate_Tanudisastro_2024_lookup_json.py.
The website's locus page and results table read the Tanudisastro2024_Details format that this lookup writes, so run
this script on the table the website queries right before publishing a regenerated index.html and locus.html. After filling the columns, regenerate the
website (`cd ../website && python3 generate_website.py`) if their descriptions in global_constants.py changed.

With --check-only, nothing is changed: the script reports how many of the lookup's first CHECK_ONLY_MAX_LOCI LocusIds
are in the table, and the table's current number of rows with a value in each Tanudisastro2024_* column.

Usage:
    python3 add_or_update_tanudisastro_2024_columns_in_catalog_table.py --check-only
    python3 add_or_update_tanudisastro_2024_columns_in_catalog_table.py   # the table the website queries
    python3 add_or_update_tanudisastro_2024_columns_in_catalog_table.py --table-id catalog_20261030_170523
"""

import argparse

from bigquery_loader_utils import add_or_update_catalog_table_columns_from_lookup_json, \
    get_tanudisastro_2024_column_values
from global_constants import TANUDISASTRO_2024_COLUMN_NAMES

CHECK_ONLY_MAX_LOCI = 10000
DEFAULT_TANUDISASTRO_2024_LOOKUP_JSON = "../data-prep/literature_review/Tanudisastro_2024/Tanudisastro_2024_lookup.json.gz"


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--tanudisastro-2024-lookup-json", default=DEFAULT_TANUDISASTRO_2024_LOOKUP_JSON,
                        help="Tanudisastro lookup JSON written by ../data-prep/literature_review/Tanudisastro_2024/"
                             "generate_Tanudisastro_2024_lookup_json.py")
    parser.add_argument("--table-id", help="catalog table to fill (default: the one the website queries, i.e. the "
                                           "TABLE_ID in ../website/header_template.html)")
    parser.add_argument("--check-only", action="store_true",
                        help="Only report how many of the lookup's LocusIds are in the table, and the table's current "
                             "number of rows with a value in each Tanudisastro2024_* column")
    args = parser.parse_args()

    add_or_update_catalog_table_columns_from_lookup_json(
        args.tanudisastro_2024_lookup_json, args.table_id, args.check_only, TANUDISASTRO_2024_COLUMN_NAMES,
        get_tanudisastro_2024_column_values, "tanudisastro_2024_staging", max_locus_ids_to_check=CHECK_ONLY_MAX_LOCI)


if __name__ == "__main__":
    main()
