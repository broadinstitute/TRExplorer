"""Create and populate a BigQuery table containing HPRC256 stratified LPS allele-frequency
histograms (per-population, per-sex, and Pop_Sex strata).

Input: a wide-format TSV produced by
    str_analysis/convert_multisample_LPS_table_to_allele_frequency_histograms.py
when run with `--stratify-by-population --stratify-by-sex`. Default path:
    ../data-prep/hprc-lps/hprc-lps.per_locus_and_motif.by_population.by_sex.256_samples.tsv.gz

The new table mirrors `HPRC256_LPS_STRATIFIED_BIGQUERY_COLUMNS` from global_constants.py:
LocusId + 17 strata x 2 histogram types = 35 columns. The unstratified base histograms
and per-locus statistics already live on the main `catalog_*` table (HPRC256_AlleleHistogram
etc.), so they are not duplicated here. The table is clustered on LocusId for fast
single-locus lookups; no partitioning.

On completion, the new table id is written to the HTML files in `../website/header_template.html`,
`../index.html`, and `../locus.html` so the frontend can target it from the dropdown UI.
"""

import argparse
import datetime
import gzip
import os
import re
import sys
import time

import tqdm
from google.cloud import bigquery

from global_constants import HPRC256_LPS_STRATIFIED_BIGQUERY_COLUMNS, HPRC256_STRATA_LABELS

PROJECT_ID = "cmg-analysis"
DATASET_ID = "tandem_repeat_explorer"
TABLE_ID = "hprc256_stratified_lps"
HTML_TABLE_ID_CONST_NAME = "HPRC256_LPS_STRATIFIED_TABLE_ID"

DEFAULT_INPUT = "../data-prep/hprc-lps/hprc-lps.per_locus_and_motif.by_population.by_sex.256_samples.tsv.gz"


def insert_with_retries(client, table_ref, rows_to_insert, batch_size=1000, max_retries=5):
    for i in range(0, len(rows_to_insert), batch_size):
        batch = rows_to_insert[i:i + batch_size]
        for retries in range(max_retries):
            try:
                errors = client.insert_rows_json(table_ref, batch)
                if errors:
                    raise RuntimeError(f"BigQuery insert errors: {errors}")
                break
            except Exception as e:
                print(f"Error inserting batch (attempt {retries + 1}/{max_retries}): {e}")
                if retries + 1 == max_retries:
                    raise
                time.sleep(5)


def update_html_table_id(html_paths, const_name, new_table_id):
    pattern = re.compile(rf"const {re.escape(const_name)}[\s]*=[\s]*'[^']*'")
    for path in html_paths:
        if not os.path.isfile(path):
            print(f"Skipping HTML update; file not found: {path}")
            continue
        with open(path, "r") as f:
            content = f.read()
        replacement = f"const {const_name} = '{new_table_id}'"
        if pattern.search(content):
            content = pattern.sub(replacement, content)
        else:
            # Insert the new constant immediately after the existing TABLE_ID line.
            content = re.sub(
                r"(const TABLE_ID[\s]*=[\s]*'catalog[^']*')",
                rf"\1\n    {replacement}",
                content,
                count=1,
            )
        with open(path, "wt") as f:
            f.write(content)
        print(f"Updated {const_name} in {path}")


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--input-tsv", default=DEFAULT_INPUT,
                        help="Stratified LPS wide-format TSV (gzipped).")
    parser.add_argument("--skip-html-update", action="store_true",
                        help="Don't rewrite TABLE_ID constants in the website HTML files.")
    args = parser.parse_args()

    input_tsv = os.path.expanduser(args.input_tsv)
    if not os.path.isfile(input_tsv):
        parser.error(f"Input file not found: {input_tsv}")

    schema = [
        bigquery.SchemaField(c["name"], c["type"], mode=c.get("mode", "NULLABLE"),
                             description=c.get("description", ""))
        for c in HPRC256_LPS_STRATIFIED_BIGQUERY_COLUMNS
    ]
    schema_field_names = [c["name"] for c in HPRC256_LPS_STRATIFIED_BIGQUERY_COLUMNS]
    # Per-column caster: TSV values arrive as strings, but BigQuery's
    # insert_rows_json rejects empty strings for FLOAT64/INT64 columns. Cast
    # explicitly here; empty strings become None (NULL in BQ).
    _NUMERIC_CASTERS = {"FLOAT64": float, "FLOAT": float, "INT64": int, "INTEGER": int}
    column_casters = {
        c["name"]: _NUMERIC_CASTERS.get(c["type"]) for c in HPRC256_LPS_STRATIFIED_BIGQUERY_COLUMNS
    }

    def coerce_cell(name, raw):
        caster = column_casters.get(name)
        if caster is None:
            return raw
        if raw is None or raw == "":
            return None
        return caster(raw)

    client = bigquery.Client(project=PROJECT_ID)
    dataset_ref = client.dataset(DATASET_ID)

    new_table_id = f"{TABLE_ID}_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}"
    new_table_ref = dataset_ref.table(new_table_id)
    new_table = bigquery.Table(new_table_ref, schema=schema)
    # Cluster on (LocusId, Interval) so per-locus + per-interval lookups still
    # hit a small cluster after the table grew to multiple rows per LocusId.
    new_table.clustering_fields = ["LocusId", "Interval"]
    client.create_table(new_table)
    print(f"Created BigQuery table {DATASET_ID}.{new_table_id}")

    # Read the wide TSV and pull out only the schema columns (LocusId + 34 stratified histograms).
    # Required input columns: LocusId plus an AlleleSizeHistogram__<stratum> and
    # BiallelicHistogram__<stratum> for every stratum in HPRC256_STRATA_LABELS.
    print(f"Reading {input_tsv}")
    with gzip.open(input_tsv, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        col_idx = {name: i for i, name in enumerate(header)}

        if "LocusId" not in col_idx:
            parser.error(f"Input is missing required 'LocusId' column. Header: {header[:5]}...")

        missing = [name for name in schema_field_names if name not in col_idx]
        if missing:
            parser.error(
                f"Input is missing {len(missing)} required column(s) (this script expects the "
                f"output of convert_multisample_LPS_table_to_allele_frequency_histograms.py "
                f"run with --stratify-by-population --stratify-by-sex). First missing: {missing[:5]}"
            )

        rows_to_insert = []
        rows_loaded = 0
        for line_num, line in enumerate(tqdm.tqdm(f, unit=" rows", unit_scale=True), start=2):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != len(header):
                raise ValueError(
                    f"Field count mismatch on line {line_num}: got {len(fields)}, expected {len(header)}"
                )
            row = {name: coerce_cell(name, fields[col_idx[name]]) for name in schema_field_names}
            rows_to_insert.append(row)
            if len(rows_to_insert) >= 1000:
                insert_with_retries(client, new_table_ref, rows_to_insert)
                rows_loaded += len(rows_to_insert)
                rows_to_insert = []
        if rows_to_insert:
            insert_with_retries(client, new_table_ref, rows_to_insert)
            rows_loaded += len(rows_to_insert)

    print(f"Loaded {rows_loaded:,d} rows into {DATASET_ID}.{new_table_id}")

    if not args.skip_html_update:
        update_html_table_id(
            ["../website/header_template.html", "../index.html", "../locus.html"],
            HTML_TABLE_ID_CONST_NAME,
            new_table_id,
        )

    print("Done!")
    return new_table_id


if __name__ == "__main__":
    main()
