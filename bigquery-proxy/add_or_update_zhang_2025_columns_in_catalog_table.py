"""Add or update the Zhang2025_* columns of an existing catalog table, filling them from the Zhang 2025 lookup JSON.

The regular way to get these columns is to rebuild the catalog table with `load_bigquery_main_table.py`, which reads
the same lookup (--zhang-2025-lookup-json). This script fills them in the live table instead, without a rebuild
(bigquery_loader_utils.add_or_update_catalog_table_columns_from_lookup_json):

  1. adds any Zhang2025_* columns (global_constants.ZHANG_2025_COLUMN_NAMES) the table does not have yet
  2. loads the lookup into a temporary staging table, one row per LocusId, with the same values
     load_bigquery_main_table.py writes (bigquery_loader_utils.get_zhang_2025_column_values)
  3. in one MERGE statement, copies the staging table's values into the rows with the same LocusId and sets the
     Zhang2025_* columns to NULL in every other row that has a value (so loci missing from the lookup end up NULL,
     even after an earlier fill); being one statement, it either applies fully or leaves the table unchanged
  4. deletes the staging table

The lookup is written by ../data-prep/literature_review/Zhang_2025/run.sh. After filling the columns, regenerate the
website (`cd ../website && python3 generate_website.py`) if their descriptions in global_constants.py changed.

With --check-only, nothing is changed: the script reports how many of the lookup's first CHECK_ONLY_MAX_LOCI LocusIds
are in the table, and the table's current counts of each Zhang2025_LengthEffectTestResult.

Usage:
    python3 add_or_update_zhang_2025_columns_in_catalog_table.py --check-only
    python3 add_or_update_zhang_2025_columns_in_catalog_table.py   # the table the website queries
    python3 add_or_update_zhang_2025_columns_in_catalog_table.py --table-id catalog_20261018_153653
"""

import argparse
import collections
import json

from bigquery_loader_utils import add_or_update_catalog_table_columns_from_lookup_json, get_zhang_2025_column_values
from global_constants import ZHANG_2025_COLUMN_NAMES

CHECK_ONLY_MAX_LOCI = 10000
DEFAULT_ZHANG_2025_LOOKUP_JSON = "../data-prep/literature_review/Zhang_2025/Zhang_2025_lookup.json.gz"


def count_staging_rows_by_test_result_and_source_library(rows):
    """Returns {(Zhang2025_LengthEffectTestResult, lengthEffectSourceLibrary in Zhang2025_Details): row count} for the
    staging rows."""
    return collections.Counter(
        (row["Zhang2025_LengthEffectTestResult"], json.loads(row["Zhang2025_Details"])["lengthEffectSourceLibrary"])
        for row in rows)


def count_rows_by_test_result_and_source_library(client, catalog_table):
    """Returns {(Zhang2025_LengthEffectTestResult, lengthEffectSourceLibrary in Zhang2025_Details): row count} for the
    rows of the table that have a test result. The library is None for 'not enough data' rows, and for every row
    filled before the lookup recorded it."""
    result = client.query(
        f"SELECT Zhang2025_LengthEffectTestResult AS test_result, "
        f"JSON_VALUE(Zhang2025_Details, '$.lengthEffectSourceLibrary') AS source_library, COUNT(*) AS n "
        f"FROM `{catalog_table}` WHERE Zhang2025_LengthEffectTestResult IS NOT NULL GROUP BY 1, 2").result()
    return {(row["test_result"], row["source_library"]): row["n"] for row in result}


def print_counts(title, counts):
    """Prints {(test result, source library): count} sorted by test result, then source library."""
    print(f"{title} ({sum(counts.values()):,} loci in total):")
    for (test_result, source_library), n in sorted(counts.items(), key=lambda item: (item[0][0], str(item[0][1]))):
        print(f"  {n:10,d}  {test_result} (from {source_library})")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--zhang-2025-lookup-json", default=DEFAULT_ZHANG_2025_LOOKUP_JSON,
                        help="Zhang 2025 lookup JSON written by ../data-prep/literature_review/Zhang_2025/run.sh")
    parser.add_argument("--table-id", help="catalog table to fill (default: the one the website queries, i.e. the "
                                           "TABLE_ID in ../website/header_template.html)")
    parser.add_argument("--check-only", action="store_true",
                        help="Only report how many of the lookup's LocusIds are in the table, and the table's current "
                             "Zhang2025_* counts")
    args = parser.parse_args()

    add_or_update_catalog_table_columns_from_lookup_json(
        args.zhang_2025_lookup_json, args.table_id, args.check_only, ZHANG_2025_COLUMN_NAMES,
        get_zhang_2025_column_values, "zhang_2025_staging", count_staging_rows_by_test_result_and_source_library,
        count_rows_by_test_result_and_source_library, print_counts, max_locus_ids_to_check=CHECK_ONLY_MAX_LOCI)


if __name__ == "__main__":
    main()
