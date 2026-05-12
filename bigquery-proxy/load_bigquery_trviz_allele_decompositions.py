"""Load an HPRC256 trviz allele-decomposition parquet (produced by
``data-prep/run_decompose_hprc_alleles.py``) into a new dated sidecar
table at ``cmg-analysis:tandem_repeat_explorer.hprc256_decomposed_alleles_<ts>``.

The table is range-partitioned on ``chrom_index`` (1..26) and clustered on
``start_0based`` so the per-locus website query can prune to a single partition
and scan a small cluster. Schema is inferred from the parquet metadata.

On completion, the new fully-qualified backticked table id is written to the
``HPRC256_DECOMPOSED_ALLELES_TABLE_ID`` constant in ``website/header_template.html``,
``locus.html`` and ``index.html`` so the frontend picks it up on next reload.

Replaces the previous ``data-prep/upload_decomposed_alleles.sh`` script.
"""

import argparse
import datetime
import os
import re
import sys

from google.cloud import bigquery

PROJECT_ID = "cmg-analysis"
DATASET_ID = "tandem_repeat_explorer"
TABLE_PREFIX = "hprc256_decomposed_alleles"
HTML_TABLE_ID_CONST_NAME = "HPRC256_DECOMPOSED_ALLELES_TABLE_ID"

# Pattern matches `decomposed_alleles.HPRC256.YYYYMMDD_HHMMSS.parquet` (anywhere in the basename).
_TIMESTAMP_FROM_FILENAME_RE = re.compile(r"\.(\d{8}_\d{6})\.parquet$")


def derive_timestamp(gcs_path):
    """Returns the timestamp embedded in the parquet filename, or current local time."""
    basename = gcs_path.rsplit("/", 1)[-1]
    match = _TIMESTAMP_FROM_FILENAME_RE.search(basename)
    if match:
        return match.group(1)
    fallback = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    print(f"Note: could not parse timestamp from filename '{basename}', using current time: {fallback}")
    return fallback


def update_html_table_id(html_paths, new_table_fqn_backticked):
    """Rewrites the ``HPRC256_DECOMPOSED_ALLELES_TABLE_ID`` const in each HTML file.

    Unlike other table-id constants in the same files (which store the bare table id),
    this constant stores the fully-qualified backticked name because the SQL uses it
    directly in ``FROM ${...}`` without re-assembling the project/dataset prefix.

    Args:
        html_paths: Iterable of HTML file paths to update.
        new_table_fqn_backticked: The single-quoted string contents to substitute,
            e.g. ``"'`cmg-analysis.tandem_repeat_explorer.hprc256_decomposed_alleles_<ts>`'"``.
    """
    # Match the assignment regardless of whether the value sits on the same line or wraps.
    value_pattern = re.compile(
        rf"const {re.escape(HTML_TABLE_ID_CONST_NAME)}\s*=\s*'`[^`]*`'",
        re.DOTALL,
    )

    for path in html_paths:
        if not os.path.isfile(path):
            print(f"Skipping HTML update; file not found: {path}")
            continue
        with open(path, "r") as f:
            content = f.read()
        if not value_pattern.search(content):
            print(f"Warning: did not find {HTML_TABLE_ID_CONST_NAME} in {path}; skipping.")
            continue
        replacement = f"const {HTML_TABLE_ID_CONST_NAME} = {new_table_fqn_backticked}"
        content = value_pattern.sub(replacement, content)
        with open(path, "wt") as f:
            f.write(content)
        print(f"Updated {HTML_TABLE_ID_CONST_NAME} in {path}")


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--input-parquet", required=True,
        help="GCS URI of the merged decomposed-alleles parquet (gs://...).",
    )
    parser.add_argument(
        "--force", action="store_true",
        help="Replace the target table if it already exists (write_disposition=WRITE_TRUNCATE).",
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Print the load configuration without submitting the job.",
    )
    parser.add_argument(
        "--skip-html-update", action="store_true",
        help="Don't rewrite the HPRC256_DECOMPOSED_ALLELES_TABLE_ID constant in the website HTML files.",
    )
    args = parser.parse_args()

    if not args.input_parquet.startswith("gs://"):
        parser.error(f"--input-parquet must be a gs:// URI (got: '{args.input_parquet}')")

    timestamp = derive_timestamp(args.input_parquet)
    new_table_id = f"{TABLE_PREFIX}_{timestamp}"
    table_dotted = f"{PROJECT_ID}.{DATASET_ID}.{new_table_id}"
    new_table_fqn_backticked = f"'`{table_dotted}`'"

    print(f"Source parquet: {args.input_parquet}")
    print(f"Target table:   {table_dotted}")
    print(f"Force overwrite: {args.force}")
    print(f"Dry run:        {args.dry_run}")
    print()

    client = bigquery.Client(project=PROJECT_ID)

    # Idempotency guard: refuse to clobber an existing table unless --force was passed.
    table_exists = False
    try:
        client.get_table(table_dotted)
        table_exists = True
    except Exception:
        table_exists = False

    if table_exists and not args.force:
        print(f"Error: target table '{table_dotted}' already exists.", file=sys.stderr)
        print("       Re-run with --force to replace it.", file=sys.stderr)
        sys.exit(1)
    if table_exists:
        print("Target table already exists; --force was provided so it will be replaced.")

    job_config = bigquery.LoadJobConfig(
        source_format=bigquery.SourceFormat.PARQUET,
        range_partitioning=bigquery.RangePartitioning(
            field="chrom_index",
            range_=bigquery.PartitionRange(start=1, end=26, interval=1),
        ),
        clustering_fields=["start_0based"],
        write_disposition=(
            bigquery.WriteDisposition.WRITE_TRUNCATE if args.force
            else bigquery.WriteDisposition.WRITE_EMPTY
        ),
    )

    if args.dry_run:
        print("Dry run: would submit load job with:")
        print(f"  source_format       = PARQUET")
        print(f"  range_partitioning  = chrom_index, 1..26, interval 1")
        print(f"  clustering_fields   = ['start_0based']")
        print(f"  write_disposition   = {job_config.write_disposition}")
        print(f"  source_uri          = {args.input_parquet}")
        print(f"  destination         = {table_dotted}")
        return

    job = client.load_table_from_uri(
        args.input_parquet,
        table_dotted,
        job_config=job_config,
    )
    print(f"Submitted load job {job.job_id}; waiting for completion...")
    job.result()
    destination = client.get_table(table_dotted)
    print(f"Loaded {destination.num_rows:,d} rows into {table_dotted}")

    if not args.skip_html_update:
        update_html_table_id(
            ["../website/header_template.html", "../locus.html", "../index.html"],
            new_table_fqn_backticked,
        )

    print("Done!")


if __name__ == "__main__":
    main()
