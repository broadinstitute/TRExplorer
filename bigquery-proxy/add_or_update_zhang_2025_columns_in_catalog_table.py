"""Add or update the Zhang2025_* columns of an existing catalog table, filling them from the Zhang 2025 lookup JSON.

The regular way to get these columns is to rebuild the catalog table with `load_bigquery_main_table.py`, which reads
the same lookup (--zhang-2025-lookup-json). This script fills them in the live table instead, without a rebuild:

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
import datetime
import gzip
import json
import os

from google.cloud import bigquery

from bigquery_loader_utils import get_zhang_2025_column_values, read_website_table_id
from global_constants import MAIN_BIGQUERY_TABLE_COLUMNS, ZHANG_2025_COLUMN_NAMES

PROJECT_ID = "cmg-analysis"
DATASET_ID = "tandem_repeat_explorer"
CHECK_ONLY_MAX_LOCI = 10000
DEFAULT_ZHANG_2025_LOOKUP_JSON = "../data-prep/literature_review/Zhang_2025/Zhang_2025_lookup.json.gz"


def zhang_2025_schema_fields():
    """The SchemaFields of the Zhang2025_* columns, as declared in MAIN_BIGQUERY_TABLE_COLUMNS."""
    types = {c["name"]: c["type"] for c in MAIN_BIGQUERY_TABLE_COLUMNS}
    return [bigquery.SchemaField(name, types[name]) for name in ZHANG_2025_COLUMN_NAMES]


def build_staging_rows(zhang_2025_lookup):
    """Returns one staging table row per locus in the Zhang 2025 lookup: its LocusId and Zhang2025_* column values.

    Args:
        zhang_2025_lookup (dict): LocusId -> lookup entry, as read from the lookup JSON

    Returns:
        list: dicts with LocusId and every column in ZHANG_2025_COLUMN_NAMES

    Raises:
        SystemExit: if the values of an entry don't cover exactly the columns in ZHANG_2025_COLUMN_NAMES
    """
    rows = []
    for locus_id, entry in zhang_2025_lookup.items():
        column_values = get_zhang_2025_column_values(entry)
        if set(column_values) != set(ZHANG_2025_COLUMN_NAMES):
            raise SystemExit(f"The Zhang2025_* values {sorted(column_values)} don't match the columns declared in "
                             f"global_constants.py: {ZHANG_2025_COLUMN_NAMES}")
        rows.append({"LocusId": locus_id, **column_values})
    return rows


def build_fill_sql(catalog_table, staging_table):
    """Returns one MERGE statement that copies the staging table's Zhang2025_* values into the rows with the same
    LocusId and sets them to NULL in every other row that has a value. A single statement either applies fully or not
    at all, so a failure never leaves the live table with its Zhang2025_* columns cleared but not refilled."""
    copy = ", ".join(f"T.{c} = S.{c}" for c in ZHANG_2025_COLUMN_NAMES)
    clear = ", ".join(f"T.{c} = NULL" for c in ZHANG_2025_COLUMN_NAMES)
    has_value = " OR ".join(f"T.{c} IS NOT NULL" for c in ZHANG_2025_COLUMN_NAMES)
    return (f"MERGE `{catalog_table}` T USING `{staging_table}` S ON T.LocusId = S.LocusId "
            f"WHEN MATCHED THEN UPDATE SET {copy} "
            f"WHEN NOT MATCHED BY SOURCE AND ({has_value}) THEN UPDATE SET {clear}")


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

    lookup_path = os.path.expanduser(args.zhang_2025_lookup_json)
    with (gzip.open if lookup_path.endswith("gz") else open)(lookup_path, "rt") as f:
        zhang_2025_lookup = json.load(f)
    rows = build_staging_rows(zhang_2025_lookup)
    lookup_counts = collections.Counter(
        (entry["LengthEffectTestResult"], entry["Details"]["lengthEffectSourceLibrary"])
        for entry in zhang_2025_lookup.values())
    print(f"Read {len(rows):,} loci from {lookup_path}")
    print_counts("Counts in the lookup", lookup_counts)

    client = bigquery.Client(project=PROJECT_ID)
    dataset_ref = f"{PROJECT_ID}.{DATASET_ID}"
    table_id = args.table_id or read_website_table_id()
    catalog_table = f"{dataset_ref}.{table_id}"
    print(f"Catalog table: {catalog_table}")

    table = client.get_table(catalog_table)
    existing_columns = {field.name for field in table.schema}

    if args.check_only:
        locus_ids = [row["LocusId"] for row in rows[:CHECK_ONLY_MAX_LOCI]]
        job = client.query(f"SELECT COUNT(*) AS n FROM `{catalog_table}` WHERE LocusId IN UNNEST(@locus_ids)",
                           job_config=bigquery.QueryJobConfig(query_parameters=[
                               bigquery.ArrayQueryParameter("locus_ids", "STRING", locus_ids)]))
        n_found = list(job.result())[0]["n"]
        print(f"{n_found:,} of the lookup's first {len(locus_ids):,} LocusIds are in {table_id} "
              f"({job.total_bytes_processed / 1e9:.2f} GB scanned)")
        if set(ZHANG_2025_COLUMN_NAMES) <= existing_columns:
            print_counts(f"Current counts in {table_id}",
                         count_rows_by_test_result_and_source_library(client, catalog_table))
        else:
            print(f"{table_id} doesn't have all the Zhang2025_* columns yet")
        print("Nothing was changed")
        return

    new_fields = [field for field in zhang_2025_schema_fields() if field.name not in existing_columns]
    if new_fields:
        table.schema = list(table.schema) + new_fields
        client.update_table(table, ["schema"])
        print(f"Added columns: {', '.join(field.name for field in new_fields)}")

    staging_table = f"{dataset_ref}.zhang_2025_staging_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}"
    try:
        client.load_table_from_json(
            rows, staging_table,
            job_config=bigquery.LoadJobConfig(
                schema=[bigquery.SchemaField("LocusId", "STRING", mode="REQUIRED")] + zhang_2025_schema_fields()),
        ).result()
        print(f"Loaded {client.get_table(staging_table).num_rows:,} of {len(rows):,} rows into {staging_table}")

        client.query(build_fill_sql(catalog_table, staging_table)).result()
        table_counts = count_rows_by_test_result_and_source_library(client, catalog_table)
        print_counts(f"Counts in {table_id} after the update", table_counts)
        if table_counts != lookup_counts:
            print(f"WARNING: these differ from the counts in the lookup, which means some of its LocusIds are not "
                  f"in {table_id}")
    finally:
        client.delete_table(staging_table, not_found_ok=True)


if __name__ == "__main__":
    main()
