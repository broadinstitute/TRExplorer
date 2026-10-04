"""Add the SpliceAI_* columns to an existing catalog table and fill them from the per-locus SpliceAI summary TSV.

The regular way to get these columns is to rebuild the catalog table with
`load_bigquery_main_table.py --spliceai-summary-tsv <tsv>`. This script fills them in the live table
instead, without a rebuild:

  1. adds any SpliceAI_* columns (global_constants.SPLICEAI_COLUMN_NAMES) the table does not have yet
  2. loads the TSV into a temporary staging table
  3. in one MERGE statement, copies the staging table's values into the rows with the same LocusId and
     sets the SpliceAI_* columns to NULL in every other row (so loci missing from the TSV end up NULL, even
     after an earlier fill); being one statement, it either applies fully or leaves the table unchanged
  4. deletes the staging table

The TSV comes from summarize_full_run_spliceai_scores_per_locus.py in ../data-prep/splicing_prediction/
(see the README.md there). After
filling the columns, regenerate the website (`cd ../website && python3 generate_website.py`) so its
filters and export dialog list them.

With --check-only, nothing is changed: the script only reports how many of the TSV's first
CHECK_ONLY_MAX_LOCI LocusIds are in the table.

Usage:
    python3 add_spliceai_columns_to_catalog_table.py TRExplorer_SpliceAI_summary.basic.tsv.gz --check-only
    python3 add_spliceai_columns_to_catalog_table.py TRExplorer_SpliceAI_summary.basic.tsv.gz             # the table the website queries
    python3 add_spliceai_columns_to_catalog_table.py TRExplorer_SpliceAI_summary.basic.tsv.gz --table-id catalog_20260918_132438
"""

import argparse
import datetime
import gzip
import os
import re
import shutil
import tempfile

from google.cloud import bigquery

from global_constants import MAIN_BIGQUERY_TABLE_COLUMNS, SPLICEAI_COLUMN_NAMES

PROJECT_ID = "cmg-analysis"
DATASET_ID = "tandem_repeat_explorer"
CHECK_ONLY_MAX_LOCI = 10000


def spliceai_schema_fields():
    """The SchemaFields of the SpliceAI_* columns, as declared in MAIN_BIGQUERY_TABLE_COLUMNS."""
    types = {c["name"]: c["type"] for c in MAIN_BIGQUERY_TABLE_COLUMNS}
    return [bigquery.SchemaField(name, types[name]) for name in SPLICEAI_COLUMN_NAMES]


def build_fill_sql(catalog_table, staging_table):
    """Returns one MERGE statement that copies the staging table's SpliceAI_* values into the rows with the same
    LocusId and sets them to NULL in every other row. A single statement either applies fully or not at all, so a
    failure never leaves the live table with its SpliceAI columns cleared but not refilled."""
    copy = ", ".join(f"T.{c} = S.{c}" for c in SPLICEAI_COLUMN_NAMES)
    clear = ", ".join(f"T.{c} = NULL" for c in SPLICEAI_COLUMN_NAMES)
    return (f"MERGE `{catalog_table}` T USING `{staging_table}` S ON T.LocusId = S.LocusId "
            f"WHEN MATCHED THEN UPDATE SET {copy} "
            f"WHEN NOT MATCHED BY SOURCE THEN UPDATE SET {clear}")


def read_website_table_id():
    """Returns the catalog table the website queries: the TABLE_ID constant in ../website/header_template.html."""
    header_template = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "website", "header_template.html")
    with open(header_template) as f:
        match = re.search(r"const TABLE_ID\s*=\s*'(catalog[^']*)'", f.read())
    if not match:
        raise SystemExit(f"No TABLE_ID constant found in {header_template}; pass --table-id")
    return match.group(1)


def read_tsv_locus_ids(tsv_path, max_loci=None):
    """Returns the TSV's LocusIds (the first max_loci of them, if given), checking its header names every SpliceAI column."""
    fopen = gzip.open if tsv_path.endswith("gz") else open
    with fopen(tsv_path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        if header != ["LocusId"] + SPLICEAI_COLUMN_NAMES:
            raise SystemExit(f"{tsv_path} header {header} does not match LocusId + {SPLICEAI_COLUMN_NAMES}")
        locus_ids = []
        for line in f:
            locus_ids.append(line.split("\t", 1)[0])
            if max_loci and len(locus_ids) >= max_loci:
                break
    return locus_ids


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("tsv", help="Per-locus SpliceAI summary TSV (LocusId plus the SpliceAI_* columns)")
    parser.add_argument("--table-id", help="catalog table to fill (default: the one the website queries, i.e. the "
                                           "TABLE_ID in ../website/header_template.html)")
    parser.add_argument("--check-only", action="store_true", help="Only report how many of the TSV's LocusIds are in the table")
    args = parser.parse_args()

    tsv_path = os.path.expanduser(args.tsv)
    client = bigquery.Client(project=PROJECT_ID)
    dataset_ref = f"{PROJECT_ID}.{DATASET_ID}"
    table_id = args.table_id or read_website_table_id()
    catalog_table = f"{dataset_ref}.{table_id}"
    print(f"Catalog table: {catalog_table}")

    if args.check_only:
        locus_ids = read_tsv_locus_ids(tsv_path, CHECK_ONLY_MAX_LOCI)
        job = client.query(f"SELECT COUNT(*) AS n FROM `{catalog_table}` WHERE LocusId IN UNNEST(@locus_ids)",
                           job_config=bigquery.QueryJobConfig(query_parameters=[
                               bigquery.ArrayQueryParameter("locus_ids", "STRING", locus_ids)]))
        n_found = list(job.result())[0]["n"]
        print(f"{n_found:,} of the TSV's first {len(locus_ids):,} LocusIds are in {table_id} "
              f"({job.total_bytes_processed / 1e9:.2f} GB scanned); nothing was changed")
        return

    n_tsv_loci = len(read_tsv_locus_ids(tsv_path))
    table = client.get_table(catalog_table)
    existing = {field.name for field in table.schema}
    new_fields = [field for field in spliceai_schema_fields() if field.name not in existing]
    if new_fields:
        table.schema = list(table.schema) + new_fields
        client.update_table(table, ["schema"])
        print(f"Added columns: {', '.join(field.name for field in new_fields)}")

    staging_table = f"{dataset_ref}.spliceai_summary_staging_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}"
    try:
        with tempfile.NamedTemporaryFile(suffix=".tsv") as uncompressed:
            fopen = gzip.open if tsv_path.endswith("gz") else open
            with fopen(tsv_path, "rb") as f:
                shutil.copyfileobj(f, uncompressed)
            uncompressed.flush()
            with open(uncompressed.name, "rb") as f:
                client.load_table_from_file(f, staging_table, job_config=bigquery.LoadJobConfig(
                    source_format=bigquery.SourceFormat.CSV, field_delimiter="\t", skip_leading_rows=1,
                    schema=[bigquery.SchemaField("LocusId", "STRING", mode="REQUIRED")] + spliceai_schema_fields(),
                )).result()
        print(f"Loaded {client.get_table(staging_table).num_rows:,} of {n_tsv_loci:,} TSV rows into {staging_table}")

        client.query(build_fill_sql(catalog_table, staging_table)).result()
        # SpliceAI_MaxDeltaScore is never empty in the TSV, so it is non-NULL exactly in the rows that were filled.
        n_filled = list(client.query(
            f"SELECT COUNTIF(SpliceAI_MaxDeltaScore IS NOT NULL) AS n FROM `{catalog_table}`").result())[0]["n"]
        print(f"Filled the SpliceAI columns of {n_filled:,} rows of {table_id} and set them to NULL in the rest "
              f"(the TSV has {n_tsv_loci:,} loci; a difference means some LocusIds are not in the table)")
    finally:
        client.delete_table(staging_table, not_found_ok=True)
    print("Next: regenerate the website so its filters and export dialog list the new columns: "
          "cd ../website && python3 generate_website.py")


if __name__ == "__main__":
    main()
