"""Tests for export/download functionality.

Searches for a query, shows all columns for the formats that export the shown columns (TSV
and JSON), then exports to each available format and verifies the file downloads
successfully and is non-empty. For BED, TSV, and
JSON formats, also verifies that the number of rows matches the total result count.
"""

import gzip
import json
import re
import urllib.request

import pytest

from test_basic_search import REGION_QUERY, assert_no_errors, do_search, enable_all_columns

# Formats whose exports have the columns shown in the results table
FORMATS_WITH_SHOWN_COLUMNS = {"tsv", "json"}

# All export format values from the dropdown
EXPORT_FORMATS = [
    "bed",
    "tsv",
    "json",
    "expansion_hunter",
    "gangstr",
    "hipstr",
    "trgt",
    "vamos",
    "longtr",
    "medaka",
    "strdust",
    "inquistr",
    "atarva",
]

# Formats that go directly to export without showing the options dialog
FORMATS_WITHOUT_DIALOG = {"vamos"}

# Formats where the exported row count should match the total search results
FORMATS_WITH_EXACT_ROW_COUNT = {"bed", "tsv", "json"}

# A region whose TR loci (OverlapGroupId 383854 in the 10/30/2026 catalog) form a chain: each overlaps the next, but
# the first and last don't touch, so the 'longest' and 'purest' export options have to keep more than one of them
OVERLAP_HANDLING_EXPORT_QUERY = "chr5:180560940-180560980"



def get_total_results(page):
    """Parse the total result count from 'Showing X-Y out of Z'."""
    text = page.locator(".bottom-pagination .total-results").inner_text()
    match = re.search(r"out of\s+([\d,]+)", text)
    return int(match.group(1).replace(",", "")) if match else 0


def select_export_format(page, fmt):
    """Select an export format using Semantic UI's dropdown API."""
    page.evaluate(
        "fmt => $('.bottom-pagination .export-dropdown').first().dropdown('set selected', fmt)",
        fmt,
    )


def do_export(page, fmt, overlap_handling=None):
    """Trigger an export and return the list of download URLs.

    Monkey-patches anchor click to capture URLs instead of navigating away.
    Handles alert/confirm dialogs and the export options dialog as needed.

    Args:
        overlap_handling: which "overlapHandling" option of the export options dialog to select ('all', 'longest' or
            'purest'), or None to leave the dialog's default
    """
    page.evaluate("""() => {
        window.__capturedDownloadUrls = [];
        const origClick = HTMLAnchorElement.prototype.click;
        HTMLAnchorElement.prototype.click = function() {
            if (this.href && this.href.includes('storage')) {
                window.__capturedDownloadUrls.push(this.href);
            } else {
                origClick.call(this);
            }
        };
    }""")

    # Removed once this export is done, so exporting twice from one page doesn't accept each dialog twice
    accept_dialog = lambda dialog: dialog.accept()
    page.on("dialog", accept_dialog)

    select_export_format(page, fmt)

    if fmt not in FORMATS_WITHOUT_DIALOG:
        page.wait_for_selector(
            "#export-column-dialog", state="visible", timeout=10000
        )
        if overlap_handling is not None:
            page.check(f"input[name='overlapHandling'][value='{overlap_handling}']", force=True)
        page.click("#export-dialog-download-button")

    page.wait_for_function(
        "() => window.__capturedDownloadUrls && window.__capturedDownloadUrls.length > 0",
        timeout=60000,
    )

    # performExport waits 2 s after each download link and then resets the cursor. Wait for that, so the next export
    # from this page starts after this one has finished.
    page.wait_for_function("() => document.body.style.cursor !== 'wait'", timeout=60000)

    page.remove_listener("dialog", accept_dialog)
    return page.evaluate("() => window.__capturedDownloadUrls")


def fetch_export_content(url):
    """Fetch exported file content from a GCS URL, handling gzip decompression."""
    with urllib.request.urlopen(url) as resp:
        data = resp.read()
    # Detect gzip by its magic bytes (0x1f 0x8b) rather than the URL suffix: signed GCS URLs
    # append query params (e.g. ...bed.gz?X-Goog-Algorithm=...), so a `.gz` suffix check misses.
    if data[:2] == b"\x1f\x8b":
        data = gzip.decompress(data)
    return data.decode("utf-8")


def count_data_rows(content, fmt):
    """Count the number of data entries in exported content."""
    if fmt == "json":
        return len(json.loads(content))
    lines = [line for line in content.strip().split("\n") if line.strip()]
    if fmt == "tsv" and lines:
        return len(lines) - 1  # TSV has a header row
    return len(lines)


@pytest.mark.parametrize("fmt", EXPORT_FORMATS)
def test_export_format(index_page, console_errors, fmt):
    """Export to a given format and verify the download is non-empty and error-free."""
    page = index_page
    # A search with a few results, to keep exports fast. Only the formats with the shown columns need all of them shown,
    # which reruns the search.
    do_search(page, REGION_QUERY)
    if fmt in FORMATS_WITH_SHOWN_COLUMNS:
        enable_all_columns(page)

    total_results = get_total_results(page)
    download_urls = do_export(page, fmt)

    assert len(download_urls) > 0, f"No download URL captured for format '{fmt}'"

    for url in download_urls:
        content = fetch_export_content(url)
        assert len(content) > 0, f"Downloaded file is empty for format '{fmt}'"

        if fmt in FORMATS_WITH_EXACT_ROW_COUNT:
            row_count = count_data_rows(content, fmt)
            assert row_count == total_results, (
                f"Expected {total_results} rows for '{fmt}' export, got {row_count}"
            )

    assert_no_errors(page, console_errors, f"export '{fmt}'")


def export_bed_intervals(page, overlap_handling):
    """Exports the current search results as BED and returns each row's (chrom, start, end, motif).

    The motif tells apart loci that have the same coordinates.
    """
    intervals = []
    for url in do_export(page, "bed", overlap_handling=overlap_handling):
        for line in fetch_export_content(url).strip().split("\n"):
            if line.strip():
                fields = line.split("\t")
                intervals.append((fields[0], int(fields[1]), int(fields[2]), fields[3]))
    return intervals


def intervals_overlap(a, b):
    return a[0] == b[0] and a[1] < b[2] and b[1] < a[2]


@pytest.mark.parametrize("overlap_handling", ["longest", "purest"])
def test_export_overlap_handling(index_page, console_errors, overlap_handling):
    """Export with the 'longest' or 'purest' option and check that no two exported loci overlap, and that each
    locus left out overlaps an exported one (for 'longest', one at least as long)."""
    page = index_page
    do_search(page, OVERLAP_HANDLING_EXPORT_QUERY)

    all_loci = export_bed_intervals(page, "all")
    kept_loci = export_bed_intervals(page, overlap_handling)

    assert len(kept_loci) > 1, f"Expected the '{overlap_handling}' export of {OVERLAP_HANDLING_EXPORT_QUERY} to keep several loci"
    assert set(kept_loci) <= set(all_loci)
    for i, a in enumerate(kept_loci):
        for b in kept_loci[i + 1:]:
            assert not intervals_overlap(a, b), f"The '{overlap_handling}' export kept overlapping loci {a} and {b}"

    left_out_loci = [locus for locus in all_loci if locus not in kept_loci]
    assert len(left_out_loci) > 0, f"Expected {OVERLAP_HANDLING_EXPORT_QUERY} to have overlapping loci to leave out"
    for locus in left_out_loci:
        overlapping_kept_loci = [k for k in kept_loci if intervals_overlap(locus, k)]
        assert overlapping_kept_loci, f"The '{overlap_handling}' export left out {locus}, which overlaps no exported locus"
        if overlap_handling == "longest":
            assert max(k[2] - k[1] for k in overlapping_kept_loci) >= locus[2] - locus[1], (
                f"The 'longest' export left out {locus}, which is longer than every exported locus it overlaps")

    assert_no_errors(page, console_errors, f"'{overlap_handling}' export")
