"""Add or update the Garg2019_* and Garg2020_* columns of an existing catalog table, filling them from the Garg 2019 and
Garg 2020 DMR lookup JSON.

The regular way to get these columns is to rebuild the catalog table with `load_bigquery_main_table.py`, which reads
the same lookup (--garg-2019-and-2020-dmr-lookup-json). This script fills them in the live table instead, without a
rebuild (bigquery_loader_utils.add_or_update_catalog_table_columns_from_lookup_json):

  1. adds any Garg2019_* and Garg2020_* columns (global_constants.GARG_2019_AND_2020_COLUMN_NAMES) the table does not
     have yet
  2. loads the lookup into a temporary staging table, one row per LocusId, with the same values
     load_bigquery_main_table.py writes (bigquery_loader_utils.get_garg_2019_and_2020_column_values)
  3. in one MERGE statement, copies the staging table's values into the rows with the same LocusId and sets the
     Garg2019_* and Garg2020_* columns to NULL in every other row that has a value (so loci missing from the lookup end
     up NULL, even after an earlier fill); being one statement, it either applies fully or leaves the table unchanged
  4. deletes the staging table

The lookup is written by ../data-prep/literature_review/Garg_2019_and_2020/run.sh. The website's locus page and the
results table columns for these data query the Garg2019_* and Garg2020_* columns, so run this script on the table the
website queries before publishing a regenerated index.html and locus.html. After filling the columns, regenerate the
website (`cd ../website && python3 generate_website.py`) if their descriptions in global_constants.py changed.

With --check-only, nothing is changed: the script reports how many of the lookup's LocusIds are in the table, and the
table's current number of rows with a value in each Garg2019_* and Garg2020_* column.

Usage:
    python3 add_or_update_garg_2019_and_2020_columns_in_catalog_table.py --check-only
    python3 add_or_update_garg_2019_and_2020_columns_in_catalog_table.py   # the table the website queries
    python3 add_or_update_garg_2019_and_2020_columns_in_catalog_table.py --table-id catalog_20261030_170523
"""

import argparse

from bigquery_loader_utils import add_or_update_catalog_table_columns_from_lookup_json, \
    get_garg_2019_and_2020_column_values
from global_constants import GARG_2019_AND_2020_COLUMN_NAMES

DEFAULT_GARG_2019_AND_2020_DMR_LOOKUP_JSON = \
    "../data-prep/literature_review/Garg_2019_and_2020/Garg_2019_and_2020_DMR_lookup.json.gz"


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--garg-2019-and-2020-dmr-lookup-json", default=DEFAULT_GARG_2019_AND_2020_DMR_LOOKUP_JSON,
                        help="Garg 2019 and Garg 2020 DMR lookup JSON written by "
                             "../data-prep/literature_review/Garg_2019_and_2020/run.sh")
    parser.add_argument("--table-id", help="catalog table to fill (default: the one the website queries, i.e. the "
                                           "TABLE_ID in ../website/header_template.html)")
    parser.add_argument("--check-only", action="store_true",
                        help="Only report how many of the lookup's LocusIds are in the table, and the table's current "
                             "counts of rows with a value in each Garg2019_* and Garg2020_* column")
    args = parser.parse_args()

    add_or_update_catalog_table_columns_from_lookup_json(
        args.garg_2019_and_2020_dmr_lookup_json, args.table_id, args.check_only, GARG_2019_AND_2020_COLUMN_NAMES,
        get_garg_2019_and_2020_column_values, "garg_2019_and_2020_staging")


if __name__ == "__main__":
    main()
