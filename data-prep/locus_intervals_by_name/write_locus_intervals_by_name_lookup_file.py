"""Writes the lookup file that the website's gene and transcript name searches need, and uploads it to Google Cloud
Storage, where the lookup_locus_intervals_by_name Cloud Function (locus-intervals-by-name-function/main.py) reads it.

The search box matches a gene or transcript name or id against the gene and transcript columns of the catalog table,
ignoring case (see GENE_AND_TRANSCRIPT_SEARCH_COLUMNS in website/index_page_template.html). For each name, the lookup
file has one row per chromosome, with the smallest start and the largest end of the loci that have the name in any of
those columns. The website searches for a name only within these intervals. That finds the same loci as a search of
the whole table, at about 1% of the cost, since the table is partitioned by chrom_index and clustered by start_0based.
The website has no other way to search for a name, so every catalog table that it points to needs a lookup file.
load_bigquery_main_table.py runs this script for each catalog table it creates.

Example:
    python3 write_locus_intervals_by_name_lookup_file.py --catalog-table-id catalog_20261030_170523
"""

import argparse
import gzip
import re
from pathlib import Path

from google.cloud import bigquery, storage

BIGQUERY_PROJECT = "cmg-analysis"
BIGQUERY_DATASET = "tandem_repeat_explorer"

# The Cloud Function reads the lookup files over HTTPS from this public bucket, at the path given by
# LOOKUP_FILE_PATH_TEMPLATE. It has its own copies of these constants, of LOOKUP_FILE_COLUMNS, SUPPORTED_GENOME_VERSIONS
# and of CATALOG_TABLE_ID_REGEXP (locus-intervals-by-name-function/main.py), and the tests check that they match.
LOOKUP_FILE_BUCKET = "tandem-repeat-explorer"
LOOKUP_FILE_PATH_TEMPLATE = "locus_intervals_by_name/{genome_version}/{catalog_table_id}.locus_intervals_by_name.tsv.gz"

SUPPORTED_GENOME_VERSIONS = ("hg38",)

# Catalog table ids as load_bigquery_main_table.py names them. The Cloud Function accepts only these.
CATALOG_TABLE_ID_REGEXP = re.compile(r"^catalog_\d{8}_\d{6}$")

# The columns that the search box matches a name against, besides LocusId (the website's
# GENE_AND_TRANSCRIPT_SEARCH_COLUMNS)
GENE_AND_TRANSCRIPT_NAME_COLUMNS = [
    "GencodeGeneName", "GencodeGeneId", "GencodeTranscriptId",
    "RefseqGeneName", "RefseqGeneId", "RefseqTranscriptId",
    "ManeGeneName", "ManeGeneId", "ManeTranscriptId",
]

# The search box also matches a name against LocusId. A search term that LOCUS_ID_REGEXP in website/header_template.html
# matches is searched for as a LocusId instead, so only the LocusIds that it doesn't match can be found by a name search,
# and only those go into the lookup file (none in catalog_20261030_170523). This pattern matches a subset of the strings
# that LOCUS_ID_REGEXP matches, so that no LocusId that a name search can find is left out: it spells out both cases
# rather than ignoring case, since BigQuery's case-insensitive matching treats a few non-ASCII letters (such as the
# Kelvin sign) as ASCII ones and JavaScript's doesn't, and BigQuery's \s matches fewer characters than JavaScript's.
LOCUS_IDS_FOUND_BY_LOCUS_ID_SEARCH_SQL_REGEXP = r"^\s*(?:[cC][hH][rR])?\w+\s*-[\d,]+-[\d,]+-[ACGTNRYMKSWHBVDacgtnrymkswhbvd]+\s*$"

LOOKUP_FILE_COLUMNS = ["name_upper", "chrom", "start_0based", "end_1based"]


def build_sql_for_locus_intervals_by_name(complete_table_id):
    """Returns the query that computes the rows of the lookup file.

    Args:
        complete_table_id (str): The catalog table, as project.dataset.table.

    Returns:
        str: A query that returns, for each uppercased name and chromosome, the smallest start_0based and largest
            end_1based of the loci with that name, sorted by name and chrom_index. The names are uppercased by
            BigQuery's UPPER, which is what the website's search uses to compare them.
    """
    name_expressions = GENE_AND_TRANSCRIPT_NAME_COLUMNS + [
        f"IF(REGEXP_CONTAINS(LocusId, r'{LOCUS_IDS_FOUND_BY_LOCUS_ID_SEARCH_SQL_REGEXP}'), NULL, LocusId)"]
    return f"""
        SELECT
            UPPER(name) AS name_upper,
            CONCAT('chr', chrom) AS chrom,
            MIN(start_0based) AS start_0based,
            MAX(end_1based) AS end_1based
        FROM `{complete_table_id}`, UNNEST([{", ".join(name_expressions)}]) AS name
        WHERE name IS NOT NULL AND name != ''
        GROUP BY name_upper, chrom_index, chrom
        ORDER BY name_upper, chrom_index
    """


def check_locus_intervals_by_name_rows(rows):
    """Raises a ValueError if the rows can't be written to the lookup file as they are.

    Args:
        rows (list): (name_upper, chrom, start_0based, end_1based) tuples.
    """
    if not rows:
        raise ValueError("The catalog table has no gene or transcript names")
    seen_name_and_chrom_pairs = set()
    for name_upper, chrom, start_0based, end_1based in rows:
        if not name_upper or re.search(r"[\t\r\n]", name_upper):
            raise ValueError(f"Name {name_upper!r} is empty or contains a tab or line break, which the TSV can't hold")
        if not re.fullmatch(r"chr\w+", chrom):
            raise ValueError(f"Unexpected chromosome {chrom!r} for name {name_upper!r}")
        if not 0 <= start_0based <= end_1based:
            raise ValueError(f"Invalid interval {chrom}:{start_0based}-{end_1based} for name {name_upper!r}")
        if (name_upper, chrom) in seen_name_and_chrom_pairs:
            raise ValueError(f"Name {name_upper!r} has more than one interval on {chrom}")
        seen_name_and_chrom_pairs.add((name_upper, chrom))


def format_locus_intervals_by_name_tsv(rows):
    """Returns the text of the lookup file.

    Args:
        rows (list): (name_upper, chrom, start_0based, end_1based) tuples, in the order to write them.

    Returns:
        str: Tab-separated lines with a header line.
    """
    lines = ["\t".join(LOOKUP_FILE_COLUMNS)]
    lines += [f"{name_upper}\t{chrom}\t{start_0based}\t{end_1based}" for name_upper, chrom, start_0based, end_1based in rows]
    return "\n".join(lines) + "\n"


def write_gzipped_text(text, output_path):
    """Writes the text gzipped, with no timestamp or file name in the gzip header, so the same text gives the same file.

    Args:
        text (str): The text to write.
        output_path (Path): The file to write.
    """
    with open(output_path, "wb") as output_file, \
            gzip.GzipFile(filename="", fileobj=output_file, mode="wb", mtime=0) as gzip_file:
        gzip_file.write(text.encode("utf-8"))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--catalog-table-id", required=True,
                        help="The catalog table in the cmg-analysis.tandem_repeat_explorer dataset, eg. "
                             "catalog_20261030_170523")
    parser.add_argument("--genome-version", default="hg38", choices=SUPPORTED_GENOME_VERSIONS,
                        help="The reference genome of the catalog table")
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent,
                        help="The directory to write the lookup file to before uploading it")
    parser.add_argument("--skip-upload", action="store_true",
                        help="Write the lookup file without uploading it to Google Cloud Storage")
    args = parser.parse_args()

    if not CATALOG_TABLE_ID_REGEXP.fullmatch(args.catalog_table_id):
        parser.error(f"--catalog-table-id {args.catalog_table_id} doesn't look like catalog_YYYYMMDD_HHMMSS, "
                     f"which is what the Cloud Function accepts")

    complete_table_id = f"{BIGQUERY_PROJECT}.{BIGQUERY_DATASET}.{args.catalog_table_id}"
    print(f"Computing the intervals of the loci with each gene and transcript name in {complete_table_id}")
    query_job = bigquery.Client(project=BIGQUERY_PROJECT).query(build_sql_for_locus_intervals_by_name(complete_table_id))
    rows = [(row.name_upper, row.chrom, row.start_0based, row.end_1based) for row in query_job.result()]
    print(f"Got {len(rows):,d} rows for {len({row[0] for row in rows}):,d} names "
          f"({query_job.total_bytes_billed / 1e9:.2f} GB billed)")
    check_locus_intervals_by_name_rows(rows)

    lookup_file_path = LOOKUP_FILE_PATH_TEMPLATE.format(
        genome_version=args.genome_version, catalog_table_id=args.catalog_table_id)
    local_path = args.output_dir / Path(lookup_file_path).name
    write_gzipped_text(format_locus_intervals_by_name_tsv(rows), local_path)
    print(f"Wrote {local_path}")

    if args.skip_upload:
        print("Skipping the upload since --skip-upload was passed")
        return

    blob = storage.Client(project=BIGQUERY_PROJECT).bucket(LOOKUP_FILE_BUCKET).blob(lookup_file_path)
    # Not cached, so that an instance of the Cloud Function that starts after a rewrite reads the new file. Instances
    # that are already running keep the copy they read, so after rewriting the file of a catalog table that the website
    # uses, redeploy the function (make deploy in locus-intervals-by-name-function) to restart them.
    blob.cache_control = "no-cache"
    blob.upload_from_filename(str(local_path), content_type="application/gzip")
    print(f"Uploaded it to gs://{LOOKUP_FILE_BUCKET}/{lookup_file_path}")


if __name__ == "__main__":
    main()
