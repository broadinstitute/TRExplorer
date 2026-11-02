"""Helpers shared by the scripts that load or fill tables in BigQuery."""

import datetime
import gzip
import json
import os
import re

import tqdm
from google.cloud import bigquery

from global_constants import MAIN_BIGQUERY_TABLE_COLUMNS

PROJECT_ID = "cmg-analysis"
DATASET_ID = "tandem_repeat_explorer"


def read_website_table_id():
    """Returns the catalog table the website queries: the TABLE_ID constant in ../website/header_template.html."""
    header_template = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "website", "header_template.html")
    with open(header_template) as f:
        match = re.search(r"const TABLE_ID\s*=\s*'(catalog[^']*)'", f.read())
    if not match:
        raise SystemExit(f"No TABLE_ID constant found in {header_template}; pass --table-id")
    return match.group(1)


def get_catalog_column_schema_fields(column_names):
    """Returns the SchemaFields of these catalog table columns, with the types declared in MAIN_BIGQUERY_TABLE_COLUMNS."""
    types = {c["name"]: c["type"] for c in MAIN_BIGQUERY_TABLE_COLUMNS}
    return [bigquery.SchemaField(name, types[name]) for name in column_names]


def build_staging_rows(lookup, column_names, get_column_values):
    """Returns one staging table row per locus in a lookup JSON: its LocusId and its values for these columns.

    Args:
        lookup (dict): LocusId -> lookup entry, as read from the lookup JSON
        column_names (list): the catalog table columns to fill
        get_column_values (function): returns {column name: value} for a lookup entry, the same values
            load_bigquery_main_table.py writes

    Returns:
        list: dicts with LocusId and every column in column_names

    Raises:
        SystemExit: if the values of an entry don't cover exactly the columns in column_names
    """
    rows = []
    for locus_id, entry in lookup.items():
        column_values = get_column_values(entry)
        if set(column_values) != set(column_names):
            raise SystemExit(f"The values {sorted(column_values)} don't match the columns declared in "
                             f"global_constants.py: {column_names}")
        rows.append({"LocusId": locus_id, **column_values})
    return rows


def build_fill_sql(catalog_table, staging_table, column_names):
    """Returns one MERGE statement that copies the staging table's values of these columns into the rows with the same
    LocusId and sets the columns to NULL in every other row that has a value. A single statement either applies fully
    or not at all, so a failure never leaves the live table with the columns cleared but not refilled."""
    copy = ", ".join(f"T.{c} = S.{c}" for c in column_names)
    clear = ", ".join(f"T.{c} = NULL" for c in column_names)
    has_value = " OR ".join(f"T.{c} IS NOT NULL" for c in column_names)
    return (f"MERGE `{catalog_table}` T USING `{staging_table}` S ON T.LocusId = S.LocusId "
            f"WHEN MATCHED THEN UPDATE SET {copy} "
            f"WHEN NOT MATCHED BY SOURCE AND ({has_value}) THEN UPDATE SET {clear}")


def count_locus_ids_in_catalog_table(client, catalog_table, locus_ids):
    """Returns how many of these LocusIds are in the catalog table, and the GB that the query scanned."""
    job = client.query(f"SELECT COUNT(*) AS n FROM `{catalog_table}` WHERE LocusId IN UNNEST(@locus_ids)",
                       job_config=bigquery.QueryJobConfig(query_parameters=[
                           bigquery.ArrayQueryParameter("locus_ids", "STRING", locus_ids)]))
    return list(job.result())[0]["n"], job.total_bytes_processed / 1e9


def add_missing_catalog_table_columns(client, catalog_table, column_names):
    """Adds the columns that the catalog table does not have yet, with the types declared in MAIN_BIGQUERY_TABLE_COLUMNS.

    Args:
        client (google.cloud.bigquery.Client): the BigQuery client
        catalog_table (str): the catalog table, as project.dataset.table
        column_names (list): the columns the table should have
    """
    table = client.get_table(catalog_table)
    existing_columns = {field.name for field in table.schema}
    new_fields = [field for field in get_catalog_column_schema_fields(column_names) if field.name not in existing_columns]
    if new_fields:
        table.schema = list(table.schema) + new_fields
        client.update_table(table, ["schema"])
        print(f"Added columns: {', '.join(field.name for field in new_fields)}")


def fill_catalog_table_columns(client, catalog_table, rows, column_names, staging_table_prefix):
    """Fills these columns of an existing catalog table from staging rows, without rebuilding the table:

      1. adds the columns the table does not have yet
      2. loads the rows into a temporary staging table
      3. runs the build_fill_sql MERGE, which copies the staging table's values into the rows with the same LocusId and
         sets the columns to NULL in every other row that has a value
      4. deletes the staging table

    Args:
        client (google.cloud.bigquery.Client): the BigQuery client
        catalog_table (str): the catalog table to fill, as project.dataset.table
        rows (list): staging rows returned by build_staging_rows
        column_names (list): the columns to fill
        staging_table_prefix (str): name of the staging table, before the timestamp that makes it unique
    """
    add_missing_catalog_table_columns(client, catalog_table, column_names)

    dataset_ref = catalog_table.rsplit(".", 1)[0]
    staging_table = f"{dataset_ref}.{staging_table_prefix}_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}"
    try:
        client.load_table_from_json(
            rows, staging_table,
            job_config=bigquery.LoadJobConfig(
                schema=[bigquery.SchemaField("LocusId", "STRING", mode="REQUIRED")]
                       + get_catalog_column_schema_fields(column_names)),
        ).result()
        print(f"Loaded {client.get_table(staging_table).num_rows:,} of {len(rows):,} rows into {staging_table}")
        client.query(build_fill_sql(catalog_table, staging_table, column_names)).result()
    finally:
        client.delete_table(staging_table, not_found_ok=True)


def add_or_update_catalog_table_columns_from_lookup_json(
        lookup_path, table_id, check_only, column_names, get_column_values, staging_table_prefix,
        count_values_in_staging_rows, count_values_in_catalog_table, print_counts, max_locus_ids_to_check=None):
    """The steps of the add_or_update_*_columns_in_catalog_table.py scripts: fills these columns of a catalog table from
    a lookup JSON with fill_catalog_table_columns, then checks that the table's counts match the lookup's. With
    check_only, nothing is changed: it only reports how many of the lookup's LocusIds are in the table, and the table's
    current counts.

    Args:
        lookup_path (str): the lookup JSON (LocusId -> lookup entry), which can be gzipped
        table_id (str): the catalog table to fill, or None for the one the website queries (read_website_table_id)
        check_only (bool): only report on the table, without changing it
        column_names (list): the columns to fill
        get_column_values (function): returns {column name: value} for a lookup entry
        staging_table_prefix (str): name of the staging table, before the timestamp that makes it unique
        count_values_in_staging_rows (function): returns a dict of counts that summarizes the staging rows
        count_values_in_catalog_table (function): (client, catalog_table) -> the same counts for the table's values
        print_counts (function): (title, counts) -> prints the counts
        max_locus_ids_to_check (int): with check_only, how many of the lookup's LocusIds to look for in the table
            (default: all of them)
    """
    lookup_path = os.path.expanduser(lookup_path)
    with (gzip.open if lookup_path.endswith("gz") else open)(lookup_path, "rt") as f:
        lookup = json.load(f)
    rows = build_staging_rows(lookup, column_names, get_column_values)
    lookup_counts = count_values_in_staging_rows(rows)
    print(f"Read {len(rows):,} loci from {lookup_path}")
    print_counts("Counts in the lookup", lookup_counts)

    client = bigquery.Client(project=PROJECT_ID)
    table_id = table_id or read_website_table_id()
    catalog_table = f"{PROJECT_ID}.{DATASET_ID}.{table_id}"
    print(f"Catalog table: {catalog_table}")

    if check_only:
        locus_ids = [row["LocusId"] for row in rows[:max_locus_ids_to_check]]
        n_found, gb_scanned = count_locus_ids_in_catalog_table(client, catalog_table, locus_ids)
        print(f"{n_found:,} of the lookup's first {len(locus_ids):,} LocusIds are in {table_id} "
              f"({gb_scanned:.2f} GB scanned)")
        if set(column_names) <= {field.name for field in client.get_table(catalog_table).schema}:
            print_counts(f"Current counts in {table_id}", count_values_in_catalog_table(client, catalog_table))
        else:
            print(f"{table_id} doesn't have all of these columns yet: {', '.join(column_names)}")
        print("Nothing was changed")
        return

    fill_catalog_table_columns(client, catalog_table, rows, column_names, staging_table_prefix)
    table_counts = count_values_in_catalog_table(client, catalog_table)
    print_counts(f"Counts in {table_id} after the update", table_counts)
    if table_counts != lookup_counts:
        print(f"WARNING: these differ from the counts in the lookup, which means some of its LocusIds are not in "
              f"{table_id}")


def create_or_replace_keep_best_non_overlapping_loci_function(client, project_id, dataset_id):
    """Creates the stored BigQuery function that the website's "longest" and "purest" export options call.

    The function gets the loci of one overlap group (one OverlapGroupId), ordered best first, and returns the ids of
    the loci to export: each locus that doesn't overlap a locus it already kept. An overlap group is a chain of loci
    (A overlaps B and B overlaps C, even when A and C don't touch), so keeping just the best locus of each group would
    also leave out loci that overlap nothing that was kept. The export query can't define this function itself,
    because the proxy (main.py) wraps that query in EXPORT DATA ... AS <query>, which accepts only a single query.

    Args:
        client (google.cloud.bigquery.Client): the BigQuery client
        project_id (str): the BigQuery project
        dataset_id (str): the dataset that holds the catalog tables
    """
    # BigQuery passes INT64 values to JavaScript as strings, hence Number(), which is exact for integers up to 2^53
    client.query(f'''
        CREATE OR REPLACE FUNCTION `{project_id}.{dataset_id}.keep_best_non_overlapping_loci`(
            loci ARRAY<STRUCT<id INT64, start_0based INT64, end_1based INT64>>)
        RETURNS ARRAY<INT64>
        LANGUAGE js AS r"""
            const kept = [];
            for (const locus of loci) {{
                const start = Number(locus.start_0based), end = Number(locus.end_1based);
                if (kept.every(k => start >= k.end || end <= k.start)) {{
                    kept.push({{id: Number(locus.id), start: start, end: end}});
                }}
            }}
            return kept.map(k => k.id);
        """
    ''').result()


def get_zhang_2025_column_values(zhang_2025_lookup_entry):
    """Returns the Zhang2025_* column values for one locus of the Zhang 2025 lookup JSON.

    Both load_bigquery_main_table.py and add_or_update_zhang_2025_columns_in_catalog_table.py use this, so a table
    filled either way holds the same values.

    Args:
        zhang_2025_lookup_entry (dict): the lookup's entry for the locus, written by
            ../data-prep/literature_review/Zhang_2025/generate_Zhang_2025_lookup_json.py

    Returns:
        dict: column name -> value, with Zhang2025_Details as a JSON string
    """
    return {
        "Zhang2025_LengthEffectTestResult": zhang_2025_lookup_entry["LengthEffectTestResult"],
        "Zhang2025_LengthEffectSlope": zhang_2025_lookup_entry["LengthEffectSlope"],
        "Zhang2025_LengthEffectFDR": zhang_2025_lookup_entry["LengthEffectFDR"],
        "Zhang2025_Details": json.dumps(zhang_2025_lookup_entry["Details"]),
    }


def get_garg_2019_and_2020_column_values(garg_lookup_entry):
    """Returns the Garg2019_* and Garg2020_* column values for one locus of the Garg 2019 and 2020 DMR lookup JSON.

    Both load_bigquery_main_table.py and add_or_update_garg_2019_and_2020_columns_in_catalog_table.py use this, so a
    table filled either way holds the same values.

    Args:
        garg_lookup_entry (dict): the lookup's entry for the locus, written by
            ../data-prep/literature_review/Garg_2019_and_2020/generate_Garg_2019_and_2020_DMR_lookup_json.py, with a
            list of DMRs for each of Garg2019_AutismDMRs, Garg2019_SchizophreniaDMRs and Garg2020_DMRs that the locus
            matched

    Returns:
        dict: column name -> value, with the *Details columns as JSON strings, and None for the columns of a dataset
            the locus didn't match
    """
    def get_methylation_change(dmrs):
        return ",".join(sorted({dmr["methylationChange"] for dmr in dmrs})) if dmrs else None

    autism_dmrs = garg_lookup_entry.get("Garg2019_AutismDMRs")
    schizophrenia_dmrs = garg_lookup_entry.get("Garg2019_SchizophreniaDMRs")
    garg_2020_dmrs = garg_lookup_entry.get("Garg2020_DMRs")
    return {
        "Garg2019_AutismDMRMethylationChange": get_methylation_change(autism_dmrs),
        "Garg2019_AutismDMRDetails": json.dumps(autism_dmrs) if autism_dmrs else None,
        "Garg2019_SchizophreniaDMRMethylationChange": get_methylation_change(schizophrenia_dmrs),
        "Garg2019_SchizophreniaDMRDetails": json.dumps(schizophrenia_dmrs) if schizophrenia_dmrs else None,
        "Garg2020_NumIndividualsWithDMR": max(
            dmr["numIndividualsWithDMR"] for dmr in garg_2020_dmrs) if garg_2020_dmrs else None,
        "Garg2020_DMRDetails": json.dumps(garg_2020_dmrs) if garg_2020_dmrs else None,
    }


def scan_for_smallest_trid(input_tsv, col_idx, locus_id_column="locus_id", interval_column="interval",
                           trid_column="trid"):
    """Find the row to keep for each (locus id, interval) in a per-locus TSV.

    Two VCF records can describe the same repeat over the same interval under different TRIDs, so
    the upstream TSVs now emit a row for each. These BigQuery schemas have no TRID column to tell
    those rows apart, and the website fetches them with an unordered LIMIT 1, so a loader keeps only
    one row per (locus id, interval): the one with the smallest TRID, which is the same row
    load_bigquery_main_table.py picks when it builds the main catalog table.

    Only the key and the TRID are held in memory, since the histogram and distribution columns these
    TSVs carry are large.

    Args:
        input_tsv: path of the gzipped TSV.
        col_idx: column name -> index, from the TSV header.
        locus_id_column: name of the locus id column ("locus_id" or "LocusId").
        interval_column: name of the interval column ("interval" or "Interval").
        trid_column: name of the TRID column ("trid" or "TRID").

    Returns:
        A dictionary of "<locus id>\\t<interval>" -> the smallest TRID seen for that key.
    """
    locus_i, interval_i = col_idx[locus_id_column], col_idx[interval_column]
    trid_i = col_idx[trid_column]
    smallest_trid_by_key = {}
    with gzip.open(input_tsv, "rt") as f:
        f.readline()
        for line in tqdm.tqdm(f, unit=" rows", unit_scale=True, desc="Finding the smallest TRID per locus"):
            fields = line.rstrip("\n").split("\t")
            key = f"{fields[locus_i]}\t{fields[interval_i]}"
            trid = fields[trid_i]
            if key not in smallest_trid_by_key or trid < smallest_trid_by_key[key]:
                smallest_trid_by_key[key] = trid
    return smallest_trid_by_key


def is_smallest_trid_row(fields, col_idx, smallest_trid_by_key, locus_id_column="locus_id",
                         interval_column="interval", trid_column="trid"):
    """Return whether this row is the one scan_for_smallest_trid chose for its (locus id, interval).

    Always True when `smallest_trid_by_key` is empty, which is how a TSV written before the TRID
    tie-break was added (no TRID column, already one row per key) loads unchanged.
    """
    if not smallest_trid_by_key:
        return True
    key = f"{fields[col_idx[locus_id_column]]}\t{fields[col_idx[interval_column]]}"
    return fields[col_idx[trid_column]] == smallest_trid_by_key[key]
