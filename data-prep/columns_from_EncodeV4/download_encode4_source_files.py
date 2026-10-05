"""Downloads the ENCODE4 and IGVF files used to compute the validated-elements column, and writes a manifest of them.

The files are chosen by querying the ENCODE portal (https://www.encodeproject.org/) and IGVF portal
(https://data.igvf.org/) REST APIs:

  validated:  files with element-level results of functional assays:
              - ENCODE GRCh38 CRISPR screen "element quantifications" TSVs (these have a Significant
                column), plus the harmonized K562 CRISPR dataset (ENCSR998YDI) used to train ENCODE-rE2G
              - ENCODE GRCh38 MPRA "element quantifications" BEDs in the "element enrichments" format
              - IGVF CRISPR screen "differential element quantifications" files. These come in many
                layouts; compute_encode4_columns_for_catalog_loci.py decides which can be used.

The manifest (encode4_source_files_manifest.tsv) lists one row per file, and
compute_encode4_columns_for_catalog_loci.py reads it to know which file is which.

Files already present with the right md5 are not downloaded again.

Usage:
    python3 download_encode4_source_files.py --dry-run         # print file counts and sizes only
    python3 download_encode4_source_files.py                   # download everything
    python3 download_encode4_source_files.py --source validated
"""

import argparse
import concurrent.futures
import hashlib
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

ENCODE_PORTAL_URL = "https://www.encodeproject.org"
IGVF_PORTAL_API_URL = "https://api.data.igvf.org"
DOWNLOAD_ATTEMPTS = 3
HARMONIZED_K562_CRISPR_ANNOTATION_ACCESSION = "ENCSR998YDI"
SOURCES = ["validated"]
MANIFEST_COLUMNS = [
    "source", "data_source", "file_accession", "dataset_accession", "assay_label", "biosample_term_name",
    "output_type", "file_type", "file_size", "md5sum", "download_url", "local_path",
]
FILE_FIELDS = [
    "accession", "dataset", "assay_title", "biosample_ontology.term_name", "output_type", "file_type",
    "file_size", "md5sum", "href",
]


def query_encode_portal(path_and_query):
    """Returns the parsed JSON response of an ENCODE portal API request.

    The portal answers a search with no results with HTTP 404, so that case returns {"@graph": []}.

    Args:
        path_and_query: URL path plus query string, e.g. "/search/?type=File&...".

    Returns:
        The parsed JSON response as a dict.
    """
    request = urllib.request.Request(ENCODE_PORTAL_URL + path_and_query, headers={"Accept": "application/json"})
    try:
        with urllib.request.urlopen(request, timeout=300) as response:
            return json.load(response)
    except urllib.error.HTTPError as e:
        if e.code == 404 and path_and_query.startswith("/search/"):
            return {"@graph": []}
        raise


def search_files(**filters):
    """Returns the released GRCh38 files matching the given ENCODE search filters.

    Args:
        **filters: search filters, e.g. output_type="peaks". A list value adds the filter once per item.

    Returns:
        List of file dicts with the FILE_FIELDS fields.
    """
    query = [("type", "File"), ("status", "released"), ("assembly", "GRCh38"), ("limit", "all")]
    for key, value in filters.items():
        for v in value if isinstance(value, list) else [value]:
            query.append((key, v))
    query += [("field", field) for field in FILE_FIELDS]
    return query_encode_portal("/search/?" + urllib.parse.urlencode(query))["@graph"]


def accession_from_path(encode_path):
    """Returns the accession from an ENCODE object path like /annotations/ENCSR800VNX/."""
    return encode_path.rstrip("/").rsplit("/", 1)[-1]


def make_manifest_row(source, file_record, assay_label):
    """Returns a manifest row dict (without local_path) for an ENCODE file record.

    Args:
        source: one of SOURCES.
        file_record: file dict returned by search_files.
        assay_label: short assay name stored in the manifest, e.g. "MPRA".

    Returns:
        Dict with the MANIFEST_COLUMNS keys except local_path.
    """
    biosample_term_name = (file_record.get("biosample_ontology") or {}).get("term_name", "")
    return {
        "source": source,
        "data_source": "ENCODE",
        "file_accession": file_record["accession"],
        "dataset_accession": accession_from_path(file_record["dataset"]),
        "assay_label": assay_label,
        "biosample_term_name": biosample_term_name,
        "output_type": file_record["output_type"],
        "file_type": file_record["file_type"],
        "file_size": file_record.get("file_size", 0),
        "md5sum": file_record["md5sum"],
        "download_url": ENCODE_PORTAL_URL + file_record["href"],
    }


def select_igvf_crispr_element_files():
    """Returns manifest rows for the released IGVF CRISPR screen "differential element quantifications" files."""
    query = [
        ("type", "TabularFile"), ("status", "released"), ("content_type", "differential element quantifications"),
        ("limit", "all"),
    ] + [("field", field) for field in [
        "accession", "file_set.accession", "file_set.samples.sample_terms.term_name", "preferred_assay_titles",
        "content_type", "file_format", "file_size", "md5sum", "href", "controlled_access",
    ]]
    request = urllib.request.Request(IGVF_PORTAL_API_URL + "/search/?" + urllib.parse.urlencode(query),
                                     headers={"Accept": "application/json"})
    with urllib.request.urlopen(request, timeout=300) as response:
        file_records = json.load(response)["@graph"]

    rows = []
    for f in file_records:
        if f.get("controlled_access"):
            continue
        assay_titles = f.get("preferred_assay_titles") or []
        if not any("CRISPR" in title or title in ("Perturb-seq", "TAP-seq") for title in assay_titles):
            continue
        sample_term_names = sorted({
            term["term_name"] for sample in f["file_set"].get("samples", []) for term in sample.get("sample_terms", [])
        })
        rows.append({
            "source": "validated",
            "data_source": "IGVF",
            "file_accession": f["accession"],
            "dataset_accession": f["file_set"]["accession"],
            "assay_label": "CRISPR",
            "biosample_term_name": ";".join(sample_term_names),
            "output_type": f["content_type"],
            "file_type": f["file_format"],
            "file_size": f.get("file_size", 0),
            "md5sum": f["md5sum"],
            "download_url": IGVF_PORTAL_API_URL + f["href"],
        })
    return rows


def select_validated_element_files():
    """Returns manifest rows for files with element-level CRISPR and MPRA results from ENCODE and IGVF."""
    rows = select_igvf_crispr_element_files()

    crispr_files = search_files(output_type="element quantifications", file_format="tsv")
    crispr_files += search_files(dataset=f"/annotations/{HARMONIZED_K562_CRISPR_ANNOTATION_ACCESSION}/")
    for f in {f["accession"]: f for f in crispr_files}.values():
        if "CRISPR" in (f.get("assay_title") or "") or f["dataset"].endswith(f"/{HARMONIZED_K562_CRISPR_ANNOTATION_ACCESSION}/"):
            rows.append(make_manifest_row("validated", f, "CRISPR"))

    for f in search_files(assay_title="MPRA", output_type="element quantifications",
                          file_type="bed element enrichments"):
        rows.append(make_manifest_row("validated", f, "MPRA"))

    return sorted(rows, key=lambda row: (row["data_source"], row["assay_label"], row["file_accession"]))


def compute_md5(path):
    """Returns the hex md5 of a file, read in 1 MB chunks."""
    md5 = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            md5.update(chunk)
    return md5.hexdigest()


def download_file(download_url, local_path, expected_md5):
    """Downloads a file unless it already exists with the expected md5, then checks the md5.

    The file is streamed to <local_path>.partial and renamed only after its md5 matches.

    Args:
        download_url: full download URL.
        local_path: Path to write.
        expected_md5: md5 reported by the ENCODE portal.

    Returns:
        True if the file was downloaded, False if it was already present.
    """
    if local_path.exists() and compute_md5(local_path) == expected_md5:
        return False

    partial_path = local_path.with_name(local_path.name + ".partial")
    for attempt in range(1, DOWNLOAD_ATTEMPTS + 1):
        try:
            with urllib.request.urlopen(download_url, timeout=600) as response, open(partial_path, "wb") as f:
                for chunk in iter(lambda: response.read(1 << 20), b""):
                    f.write(chunk)
            break
        except (urllib.error.URLError, ConnectionError, TimeoutError) as e:
            if attempt == DOWNLOAD_ATTEMPTS:
                raise
            print(f"Retrying {download_url} after attempt {attempt} failed: {e}")
            time.sleep(10 * attempt)

    actual_md5 = compute_md5(partial_path)
    if actual_md5 != expected_md5:
        partial_path.unlink()
        raise ValueError(f"md5 mismatch for {download_url}: expected {expected_md5}, got {actual_md5}")
    partial_path.rename(local_path)
    return True


def write_manifest(rows, manifest_path):
    """Writes the manifest rows as a TSV with MANIFEST_COLUMNS."""
    with open(manifest_path, "wt") as f:
        f.write("\t".join(MANIFEST_COLUMNS) + "\n")
        for row in rows:
            f.write("\t".join(str(row[column]) for column in MANIFEST_COLUMNS) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--source", choices=SOURCES, action="append",
                        help="Only download this source. Can be repeated. Default: all sources.")
    parser.add_argument("--output-dir", type=Path, default=Path("encode4_source_files"),
                        help="Directory to download into. Each source gets its own subdirectory.")
    parser.add_argument("--manifest", type=Path, default=Path("encode4_source_files_manifest.tsv"),
                        help="Manifest TSV to write.")
    parser.add_argument("--n-parallel-downloads", type=int, default=4)
    parser.add_argument("--dry-run", action="store_true", help="Print file counts and sizes without downloading.")
    args = parser.parse_args()

    selectors = {
        "validated": select_validated_element_files,
    }
    rows = []
    for source in args.source or SOURCES:
        source_rows = selectors[source]()
        total_gb = sum(int(row["file_size"]) for row in source_rows) / 1e9
        print(f"{source}: {len(source_rows):,d} files, {total_gb:.2f} GB")
        for row in source_rows:
            file_name = row["download_url"].rsplit("/", 1)[-1]
            row["local_path"] = str(args.output_dir / source / file_name)
        rows += source_rows

    total_gb = sum(int(row["file_size"]) for row in rows) / 1e9
    print(f"Total: {len(rows):,d} files, {total_gb:.2f} GB")
    if args.dry_run:
        return

    for source in {row["source"] for row in rows}:
        (args.output_dir / source).mkdir(parents=True, exist_ok=True)

    n_downloaded = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.n_parallel_downloads) as executor:
        futures = {
            executor.submit(download_file, row["download_url"], Path(row["local_path"]), row["md5sum"]): row
            for row in rows
        }
        for i, future in enumerate(concurrent.futures.as_completed(futures), start=1):
            n_downloaded += future.result()
            if i % 100 == 0 or i == len(futures):
                print(f"  {i:,d} of {len(futures):,d} files done ({n_downloaded:,d} downloaded)")

    # Keep the existing manifest's rows for the other SOURCES that were not downloaded in this run.
    selected_sources = {row["source"] for row in rows}
    if args.manifest.exists():
        with open(args.manifest, "rt") as f:
            header = f.readline().rstrip("\n").split("\t")
            rows += [row for row in (dict(zip(header, line.rstrip("\n").split("\t"))) for line in f)
                     if row["source"] in SOURCES and row["source"] not in selected_sources]

    write_manifest(rows, args.manifest)
    print(f"Wrote {len(rows):,d} rows to {args.manifest}")


if __name__ == "__main__":
    main()
