"""Create and populate two BigQuery tables containing HPRC256 stratified joint distributions:

  1. Allele size x purity:
     hprc256_allele_purity_<timestamp>
     Source: trgt-hprc.allele_size_purity.stratified.by_population.by_sex.<N>_samples.tsv.gz

  2. Allele size x methylation:
     hprc256_methylation_<timestamp>
     Source: trgt-hprc.methylation.stratified.by_population.by_sex.<N>_samples.tsv.gz

Both inputs are produced by
    data-prep/hprc-lps/compute_allele_size_purity_and_methylation_distributions_from_vcf.py
when run with `--stratify-by-population --stratify-by-sex`. Each table is keyed on
(locus_id, interval) and clustered on (locus_id, interval) (no partitioning); the schemas live in
global_constants.py as HPRC256_ALLELE_PURITY_BIGQUERY_TABLE_COLUMNS and
HPRC256_METHYLATION_BIGQUERY_TABLE_COLUMNS.

On completion, the new table ids are written to the HTML files in
`../website/header_template.html`, `../index.html`, and `../locus.html` so the frontend
can target them from the dropdown UI.
"""

import argparse
import atexit
import datetime
import gzip
import json
import os
import re
import tempfile

import tqdm
from google.cloud import bigquery

from global_constants import (
    HPRC256_ALLELE_PURITY_BIGQUERY_TABLE_COLUMNS,
    HPRC256_METHYLATION_BIGQUERY_TABLE_COLUMNS,
)

PROJECT_ID = "cmg-analysis"
DATASET_ID = "tandem_repeat_explorer"

DEFAULT_PURITY_INPUT = "../data-prep/hprc-lps/trgt-hprc.allele_size_purity.stratified.by_population.by_sex.256_samples.tsv.gz"
DEFAULT_METHYLATION_INPUT = "../data-prep/hprc-lps/trgt-hprc.methylation.stratified.by_population.by_sex.256_samples.tsv.gz"

PURITY_TABLE_ID = "hprc256_allele_purity"
METHYLATION_TABLE_ID = "hprc256_methylation"
PURITY_HTML_CONST = "HPRC256_ALLELE_PURITY_TABLE_ID"
METHYLATION_HTML_CONST = "HPRC256_METHYLATION_TABLE_ID"


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
            content, n = re.subn(
                r"(const TABLE_ID[\s]*=[\s]*'catalog[^']*')",
                rf"\1\n    {replacement}",
                content,
                count=1,
            )
            if n == 0:
                print(f"Warning: could not find anchor 'const TABLE_ID = ...' in {path}; skipping.")
                continue
        with open(path, "wt") as f:
            f.write(content)
        print(f"Updated {const_name} in {path}")


def load_table(client, dataset_ref, input_tsv, table_prefix, schema_columns, parser):
    """Create a new table from `schema_columns` and bulk-load rows from `input_tsv` into it.
    Returns the new table id."""
    if not os.path.isfile(input_tsv):
        parser.error(f"Input file not found: {input_tsv}")

    schema_field_names = [c["name"] for c in schema_columns]
    schema = [
        bigquery.SchemaField(c["name"], c["type"], mode=c.get("mode", "NULLABLE"),
                             description=c.get("description", ""))
        for c in schema_columns
    ]

    new_table_id = f"{table_prefix}_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}"
    new_table_ref = dataset_ref.table(new_table_id)
    new_table = bigquery.Table(new_table_ref, schema=schema)
    # Cluster on (locus_id, interval) so per-locus + per-interval lookups still
    # hit a small cluster after the table grew to multiple rows per locus_id.
    new_table.clustering_fields = ["locus_id", "interval"]
    client.create_table(new_table)
    print(f"Created BigQuery table {DATASET_ID}.{new_table_id}")

    print(f"Reading {input_tsv}")
    ndjson_path = os.path.join(tempfile.gettempdir(), f"{new_table_id}.jsonl")
    atexit.register(lambda: os.path.exists(ndjson_path) and os.remove(ndjson_path))
    rows_loaded = 0
    with gzip.open(input_tsv, "rt") as f, open(ndjson_path, "wt") as ndjson_file:
        header = f.readline().rstrip("\n").split("\t")
        col_idx = {name: i for i, name in enumerate(header)}

        if "locus_id" not in col_idx:
            parser.error(f"Input {input_tsv} is missing required 'locus_id' column. Header: {header[:5]}...")

        missing = [name for name in schema_field_names if name not in col_idx]
        if missing:
            parser.error(
                f"Input {input_tsv} is missing {len(missing)} required column(s) (run the upstream "
                f"compute_allele_size_purity_and_methylation_distributions_from_vcf.py with "
                f"--stratify-by-population --stratify-by-sex). First missing: {missing[:5]}"
            )

        for line_num, line in enumerate(tqdm.tqdm(f, unit=" rows", unit_scale=True), start=2):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != len(header):
                raise ValueError(
                    f"Field count mismatch on line {line_num}: got {len(fields)}, expected {len(header)}"
                )
            row = {name: fields[col_idx[name]] for name in schema_field_names}
            ndjson_file.write(json.dumps(row) + "\n")
            rows_loaded += 1

    print(f"Loading {DATASET_ID}.{new_table_id} from {ndjson_path}")
    load_job_config = bigquery.LoadJobConfig(
        schema=schema,
        source_format=bigquery.SourceFormat.NEWLINE_DELIMITED_JSON,
        write_disposition=bigquery.WriteDisposition.WRITE_APPEND,
    )
    with open(ndjson_path, "rb") as ndjson_file_obj:
        load_job = client.load_table_from_file(ndjson_file_obj, new_table_ref, job_config=load_job_config)
    load_job.result()
    os.remove(ndjson_path)

    print(f"Loaded {rows_loaded:,d} rows into {DATASET_ID}.{new_table_id}")
    return new_table_id


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--purity-tsv", default=DEFAULT_PURITY_INPUT,
                        help="Stratified allele-size-and-purity TSV (gzipped).")
    parser.add_argument("--methylation-tsv", default=DEFAULT_METHYLATION_INPUT,
                        help="Stratified allele-size-and-methylation TSV (gzipped).")
    parser.add_argument("--skip-html-update", action="store_true",
                        help="Don't rewrite TABLE_ID constants in the website HTML files.")
    args = parser.parse_args()

    client = bigquery.Client(project=PROJECT_ID)
    dataset_ref = client.dataset(DATASET_ID)

    purity_table_id = load_table(
        client, dataset_ref, os.path.expanduser(args.purity_tsv),
        PURITY_TABLE_ID, HPRC256_ALLELE_PURITY_BIGQUERY_TABLE_COLUMNS, parser,
    )
    methylation_table_id = load_table(
        client, dataset_ref, os.path.expanduser(args.methylation_tsv),
        METHYLATION_TABLE_ID, HPRC256_METHYLATION_BIGQUERY_TABLE_COLUMNS, parser,
    )

    if not args.skip_html_update:
        html_paths = ["../website/header_template.html", "../index.html", "../locus.html"]
        update_html_table_id(html_paths, PURITY_HTML_CONST, purity_table_id)
        update_html_table_id(html_paths, METHYLATION_HTML_CONST, methylation_table_id)

    print("Done!")


if __name__ == "__main__":
    main()
