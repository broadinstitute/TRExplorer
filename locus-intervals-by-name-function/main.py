"""Cloud Function that returns where the loci with a given gene or transcript name or id are in a TRExplorer catalog
table, so that the website can search for the name within those intervals instead of across the whole table.

Request (POST, JSON):
    {"names": ["ATXN1", "ENSG00000124788"], "genome_version": "hg38", "catalog_table_id": "catalog_20261030_170523"}

Response (JSON):
    {"intervals_by_name": {"ATXN1": [{"chrom": "chr6", "start_0based": 16300565, "end_1based": 16761730}],
                           "ENSG00000124788": [...]}}

Names are compared ignoring case. A name's intervals, one per chromosome, span the catalog loci that have the name in
any of the gene and transcript name and id columns (or as an unusual LocusId; see
data-prep/locus_intervals_by_name/write_locus_intervals_by_name_lookup_file.py, which writes the lookup files). An empty
list means that no locus in the catalog table has the name.
"""

import bisect
import functools
import gzip
import re
import urllib.error
import urllib.request

from flask import jsonify

# The lookup files are public, so they're downloaded over HTTPS with no credentials. This, LOOKUP_FILE_HEADER,
# SUPPORTED_GENOME_VERSIONS and CATALOG_TABLE_ID_REGEXP must match the constants of
# data-prep/locus_intervals_by_name/write_locus_intervals_by_name_lookup_file.py, which its tests check.
LOOKUP_FILE_URL_TEMPLATE = (
    "https://storage.googleapis.com/tandem-repeat-explorer/locus_intervals_by_name/{genome_version}/"
    "{catalog_table_id}.locus_intervals_by_name.tsv.gz")
LOOKUP_FILE_HEADER = "name_upper\tchrom\tstart_0based\tend_1based"
LOOKUP_FILE_DOWNLOAD_TIMEOUT_SECONDS = 30

SUPPORTED_GENOME_VERSIONS = ("hg38",)
CATALOG_TABLE_ID_REGEXP = re.compile(r"^catalog_\d{8}_\d{6}$")

# Well above any search that can succeed: the BigQuery query that searches for the names can hold only about 2,000 of
# them, and the website refuses an uploaded file with more than 2,500 (MAX_NAMES_PER_UPLOADED_FILE in
# website/index_page_template.html) before looking them up
MAX_NAMES_PER_REQUEST = 10_000
MAX_NAME_LENGTH = 1_000

# A lookup file takes about 35 MB of memory once read. Two of them cover the live catalog table plus the one before it,
# which is kept as a fallback.
MAX_LOOKUP_FILES_IN_MEMORY = 2


class LookupFileNotFoundError(Exception):
    """Raised when there's no lookup file for the requested genome version and catalog table."""


def parse_lookup_file_text(text):
    """Reads the lines of a lookup file into a sorted list, so that a name's lines can be found by binary search.

    This takes about a tenth of the time and a third of the memory of parsing every line into a dict, which matters
    because it happens when an instance of the function starts, which a search then waits for.

    Args:
        text (str): The decompressed text of a lookup file.

    Returns:
        list: The lines after the header, without empty lines, sorted.
    """
    # Split on line breaks only, since splitlines() would also split a name at characters like a form feed
    lines = text.split("\n")
    if lines[0] != LOOKUP_FILE_HEADER:
        raise ValueError(f"Unexpected lookup file header: {lines[0][:200]!r}")
    return sorted(line for line in lines[1:] if line)


def find_intervals_of_name(name_upper, sorted_lookup_file_lines):
    """Finds the intervals of an uppercased name.

    Args:
        name_upper (str): The name, uppercased.
        sorted_lookup_file_lines (list): See parse_lookup_file_text.

    Returns:
        list: [{"chrom": ..., "start_0based": ..., "end_1based": ...}, ...], or an empty list if no locus has the name.
    """
    # The lookup file's names have no tabs or line breaks, and one with a tab would match the start of other lines
    if "\t" in name_upper or "\n" in name_upper:
        return []
    # The lines that start with the name followed by a tab come one after another in the sorted list, starting at the
    # first line that sorts at or after that prefix
    line_prefix = f"{name_upper}\t"
    intervals = []
    line_index = bisect.bisect_left(sorted_lookup_file_lines, line_prefix)
    while line_index < len(sorted_lookup_file_lines) and sorted_lookup_file_lines[line_index].startswith(line_prefix):
        _, chrom, start_0based, end_1based = sorted_lookup_file_lines[line_index].split("\t")
        intervals.append({"chrom": chrom, "start_0based": int(start_0based), "end_1based": int(end_1based)})
        line_index += 1
    return intervals


def download_and_parse_lookup_file(lookup_file_url):
    """Downloads and parses a lookup file.

    Args:
        lookup_file_url (str): The file's URL.

    Returns:
        list: See parse_lookup_file_text.

    Raises:
        LookupFileNotFoundError: If there's no file at the URL.
    """
    try:
        with urllib.request.urlopen(lookup_file_url, timeout=LOOKUP_FILE_DOWNLOAD_TIMEOUT_SECONDS) as response:
            compressed_bytes = response.read()
    except urllib.error.HTTPError as e:
        if e.code == 404:
            raise LookupFileNotFoundError(f"There's no lookup file at {lookup_file_url}") from e
        raise
    return parse_lookup_file_text(gzip.decompress(compressed_bytes).decode("utf-8"))


# Kept between requests, so that a running instance of the function reads each lookup file only once. A failed read
# isn't kept, so it's retried on the next request.
@functools.lru_cache(maxsize=MAX_LOOKUP_FILES_IN_MEMORY)
def get_sorted_lookup_file_lines(genome_version, catalog_table_id):
    """Returns the parsed lookup file of a catalog table.

    Args:
        genome_version (str): The catalog table's reference genome.
        catalog_table_id (str): The catalog table.

    Returns:
        list: See parse_lookup_file_text.
    """
    return download_and_parse_lookup_file(
        LOOKUP_FILE_URL_TEMPLATE.format(genome_version=genome_version, catalog_table_id=catalog_table_id))


def parse_request_json(request_json):
    """Reads and checks the parameters of a request.

    Args:
        request_json: The request body, parsed from JSON.

    Returns:
        tuple: (names, genome_version, catalog_table_id).

    Raises:
        ValueError: If a parameter is missing or invalid. The message says which one, and why.
    """
    if not isinstance(request_json, dict):
        raise ValueError("The request body must be a JSON object")
    names = request_json.get("names")
    if not isinstance(names, list) or not names or not all(isinstance(name, str) for name in names):
        raise ValueError("names must be a non-empty list of strings")
    if len(names) > MAX_NAMES_PER_REQUEST:
        raise ValueError(f"names has {len(names):,d} names, more than the {MAX_NAMES_PER_REQUEST:,d} allowed per request")
    if any(len(name) > MAX_NAME_LENGTH for name in names):
        raise ValueError(f"names has a name longer than {MAX_NAME_LENGTH:,d} characters")
    genome_version = request_json.get("genome_version")
    if genome_version not in SUPPORTED_GENOME_VERSIONS:
        raise ValueError(f"genome_version must be one of: {', '.join(SUPPORTED_GENOME_VERSIONS)}")
    catalog_table_id = request_json.get("catalog_table_id")
    if not isinstance(catalog_table_id, str) or not CATALOG_TABLE_ID_REGEXP.fullmatch(catalog_table_id):
        raise ValueError("catalog_table_id must be a catalog table id like catalog_20261030_170523")
    return names, genome_version, catalog_table_id


def look_up_intervals_of_names(names, sorted_lookup_file_lines):
    """Looks up each name's intervals.

    Args:
        names (list): The names to look up, as searched for.
        sorted_lookup_file_lines (list): See parse_lookup_file_text.

    Returns:
        dict: {name: [{"chrom": ..., "start_0based": ..., "end_1based": ...}, ...]} for each of the names, with an empty
            list for a name that no locus has. The lookup file's names were uppercased by BigQuery's UPPER, and Python's
            upper() gives the same result for them, since they're all ASCII.
    """
    return {name: find_intervals_of_name(name.upper(), sorted_lookup_file_lines) for name in names}


def lookup_locus_intervals_by_name(request):
    """Cloud Function entry point: returns the intervals of the loci that have each of the requested names.

    Args:
        request (flask.Request): A POST request whose JSON body is described in this module's docstring.

    Returns:
        tuple: (response body, status code, response headers).
    """
    response_headers = {
        "Access-Control-Allow-Origin": "*",
        "Access-Control-Allow-Methods": "POST",
        "Access-Control-Allow-Headers": "Content-Type",
    }
    if request.method == "OPTIONS":
        return "", 204, response_headers
    if request.method != "POST":
        return jsonify({"error": "Only POST requests are supported"}), 405, response_headers

    try:
        names, genome_version, catalog_table_id = parse_request_json(request.get_json(silent=True))
    except ValueError as e:
        return jsonify({"error": str(e)}), 400, response_headers

    try:
        sorted_lookup_file_lines = get_sorted_lookup_file_lines(genome_version, catalog_table_id)
    except LookupFileNotFoundError as e:
        print(f"ERROR: {e}")
        return jsonify({"error": f"There's no name lookup file for {genome_version} catalog table {catalog_table_id}"}), \
            404, response_headers
    except Exception as e:
        print(f"ERROR: Reading the lookup file of {genome_version} catalog table {catalog_table_id} failed: {e}")
        return jsonify({"error": f"Reading the name lookup file failed: {e}"}), 500, response_headers

    return jsonify({"intervals_by_name": look_up_intervals_of_names(names, sorted_lookup_file_lines)}), 200, \
        response_headers
