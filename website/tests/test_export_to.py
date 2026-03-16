"""Tests for export/download functionality.

Searches for a query, enables all columns, then exports to each available format
and verifies the file downloads successfully and is non-empty. For BED, TSV, and
JSON formats, also verifies that the number of rows matches the total result count.
"""

import gzip
import json
import re
import urllib.request

import pytest

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
]

# Formats that go directly to export without showing the options dialog
FORMATS_WITHOUT_DIALOG = {"vamos"}

# Formats where the exported row count should match the total search results
FORMATS_WITH_EXACT_ROW_COUNT = {"bed", "tsv", "json"}

# A query with a small number of results to keep exports fast
EXPORT_QUERY = "chr4:3074000-3075000"


def do_search(page, query):
    """Fill the search box, click Search, and wait for results."""
    page.fill("#search-query", query)
    page.click("#search-button")
    page.wait_for_selector("#results-table tbody tr", timeout=30000)


def enable_all_columns(page):
    """Open column settings dialog, check all boxes, and click Apply."""
    page.click("#results-table-settings-button")
    page.wait_for_selector(
        "#results-table-column-settings-dialog", state="visible", timeout=5000
    )

    for cb in page.query_selector_all(
        ".results-table-column-settings-checkboxes input[type='checkbox']"
    ):
        if not cb.is_checked():
            cb.check(force=True)

    page.click("#results-table-column-settings-apply-button")
    page.wait_for_selector("#results-table tbody tr", timeout=30000)


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


def do_export(page, fmt):
    """Trigger an export and return the list of download URLs.

    Monkey-patches anchor click to capture URLs instead of navigating away.
    Handles alert/confirm dialogs and the export options dialog as needed.
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

    page.on("dialog", lambda dialog: dialog.accept())

    select_export_format(page, fmt)

    if fmt not in FORMATS_WITHOUT_DIALOG:
        page.wait_for_selector(
            "#export-column-dialog", state="visible", timeout=10000
        )
        page.click("#export-dialog-download-button")

    page.wait_for_function(
        "() => window.__capturedDownloadUrls && window.__capturedDownloadUrls.length > 0",
        timeout=60000,
    )

    return page.evaluate("() => window.__capturedDownloadUrls")


def fetch_export_content(url):
    """Fetch exported file content from a GCS URL, handling gzip decompression.

    Converts storage.cloud.google.com URLs to storage.googleapis.com for
    direct public access.
    """
    url = url.replace("storage.cloud.google.com/", "storage.googleapis.com/")
    with urllib.request.urlopen(url) as resp:
        data = resp.read()
    if url.endswith(".gz"):
        data = gzip.decompress(data)
    return data.decode("utf-8")


def count_data_rows(content, fmt):
    """Count the number of data entries in exported content."""
    if fmt in ("json", "expansion_hunter"):
        return len(json.loads(content))
    lines = [line for line in content.strip().split("\n") if line.strip()]
    if fmt == "tsv" and lines:
        return len(lines) - 1  # TSV has a header row
    return len(lines)


@pytest.mark.parametrize("fmt", EXPORT_FORMATS)
def test_export_format(index_page, console_errors, fmt):
    """Export to a given format and verify the download is non-empty and error-free."""
    page = index_page
    do_search(page, EXPORT_QUERY)
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

    assert console_errors == [], (
        f"JavaScript errors during export '{fmt}': {console_errors}"
    )
