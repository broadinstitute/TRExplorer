"""Update an existing BigQuery table schema to match the current BIGQUERY_COLUMNS in global_constants.py.

Adds any columns present in global_constants.py but missing from the table.
New columns will have NULL values for all existing rows until data is reloaded.

Usage:
    python update_bigquery_schema.py                         # updates the most recent catalog_* table
    python update_bigquery_schema.py catalog_20260306_055435 # updates a specific table
"""

import argparse
from google.cloud import bigquery

from global_constants import BIGQUERY_COLUMNS

PROJECT_ID = "cmg-analysis"
DATASET_ID = "tandem_repeat_explorer"

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("table_id", nargs="?", help="BigQuery table ID (default: most recent catalog_* table)")
args = parser.parse_args()

client = bigquery.Client()
dataset_ref = f"{PROJECT_ID}.{DATASET_ID}"

if args.table_id:
    table_id = args.table_id
else:
    tables = [t.table_id for t in client.list_tables(dataset_ref) if t.table_id.startswith("catalog_")]
    if not tables:
        parser.error(f"No catalog_* tables found in {dataset_ref}")
    table_id = sorted(tables)[-1]
    print(f"Using most recent table: {table_id}")

table = client.get_table(f"{dataset_ref}.{table_id}")

existing_columns = {field.name for field in table.schema}
desired_schema = [
    bigquery.SchemaField(col["name"], col["type"], mode=col.get("mode", "NULLABLE"))
    for col in BIGQUERY_COLUMNS
]

new_columns = [field for field in desired_schema if field.name not in existing_columns]
if not new_columns:
    print("No new columns to add. Schema is already up to date.")
else:
    print(f"Adding {len(new_columns)} new column(s): {', '.join(f.name for f in new_columns)}")
    table.schema = list(table.schema) + new_columns
    client.update_table(table, ["schema"])
    print("Done.")
