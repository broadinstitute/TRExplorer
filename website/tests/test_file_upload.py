"""Tests for the panel-file upload feature on the index page.

A (+) button in the search box opens an upload modal (drag-drop / browse) that accepts
TXT/TSV files (first column = TR LocusIds, hg38 regions, or gene/transcript identifiers)
and BED files (chrom, start, end), optionally gzip-compressed. Uploading restricts the
results table to the union of the file's loci/regions.
"""

import gzip

import pytest

from test_basic_search import real_errors, do_search, click_next_page, get_total_pages


def wait_for_results(page):
    """Wait for a search to populate the results table."""
    page.wait_for_selector("#results-table tbody tr", timeout=30000)


def upload_files(page, paths):
    """Set the hidden file input, which triggers parsing + search."""
    page.set_input_files("#tr-file-input", [str(p) for p in paths])
    page.wait_for_selector("#tr-file-chip-container .ui.label", timeout=10000)


def test_upload_button_opens_modal(index_page):
    page = index_page
    page.click("#tr-file-upload-button")
    page.wait_for_selector("#tr-file-upload-modal", state="visible", timeout=5000)
    assert page.is_visible("#tr-file-dropzone")


def test_upload_txt_regions_and_genes(index_page, console_errors, tmp_path):
    """A TXT with a region and a gene symbol filters the table and shows a chip."""
    page = index_page
    f = tmp_path / "panel.txt"
    f.write_text("chr4:3074000-3075000\nATXN1\n")

    upload_files(page, [f])
    wait_for_results(page)

    chip = page.inner_text("#tr-file-chip-container")
    assert "panel.txt" in chip
    assert "2 loci" in chip
    rows = page.query_selector_all("#results-table tbody tr")
    assert len(rows) > 0
    assert real_errors(console_errors) == [], f"JS errors: {real_errors(console_errors)}"


def test_upload_txt_with_header_row(index_page, tmp_path):
    """A leading header row (e.g. 'LocusId') is skipped, not searched."""
    page = index_page
    f = tmp_path / "with_header.txt"
    f.write_text("region\nchr4:3074000-3075000\n")

    upload_files(page, [f])
    wait_for_results(page)
    # header 'region' is skipped, leaving exactly one item
    assert "1 locus" in page.inner_text("#tr-file-chip-container")


def test_upload_locus_id(index_page, console_errors, tmp_path):
    """A TXT containing a real LocusId (chr-start-end-motif) matches that locus."""
    page = index_page
    # First find a real LocusId by running a region search and reading a result row.
    page.fill("#search-query", "chr4:3074000-3075000")
    page.click("#search-button")
    wait_for_results(page)
    locus_id = page.get_attribute(".row-locus-page-link", "data-locus-id")
    assert locus_id and "-" in locus_id

    # Clear the search box so only the uploaded file drives the query.
    page.fill("#search-query", "")
    f = tmp_path / "one_locus.txt"
    f.write_text(locus_id + "\n")

    upload_files(page, [f])
    wait_for_results(page)
    assert "1 locus" in page.inner_text("#tr-file-chip-container")
    locus_ids = [
        el.get_attribute("data-locus-id")
        for el in page.query_selector_all(".row-locus-page-link")
    ]
    assert locus_id in locus_ids
    assert real_errors(console_errors) == [], f"JS errors: {real_errors(console_errors)}"


def test_upload_locus_id_with_chr_prefix(index_page, console_errors, tmp_path):
    """A "chr"-prefixed LocusId (as commonly pasted/exported, e.g. chr4-...-CAG) must still match,
    even though the LocusId is stored in the catalog without the "chr" prefix."""
    page = index_page
    page.fill("#search-query", "chr4:3074000-3075000")
    page.click("#search-button")
    wait_for_results(page)
    locus_id = page.get_attribute(".row-locus-page-link", "data-locus-id")
    assert locus_id and "-" in locus_id
    assert not locus_id.lower().startswith("chr")  # catalog LocusIds have no "chr" prefix

    page.fill("#search-query", "")
    f = tmp_path / "one_locus_chr_prefixed.txt"
    f.write_text(f"chr{locus_id}\n")

    upload_files(page, [f])
    wait_for_results(page)
    assert "1 locus" in page.inner_text("#tr-file-chip-container")
    locus_ids = [
        el.get_attribute("data-locus-id")
        for el in page.query_selector_all(".row-locus-page-link")
    ]
    assert locus_id in locus_ids
    assert real_errors(console_errors) == [], f"JS errors: {real_errors(console_errors)}"


def test_drop_file_on_main_search_box(index_page, console_errors, tmp_path):
    """Dropping a file directly on the main search box (not just in the upload modal) behaves the
    same as dropping it into the modal's dropzone."""
    page = index_page
    data_transfer = page.evaluate_handle(
        """() => {
            const dt = new DataTransfer()
            const file = new File(['chr4:3074000-3075000\\n'], 'dropped_on_search_box.txt', { type: 'text/plain' })
            dt.items.add(file)
            return dt
        }"""
    )
    page.dispatch_event("#tr-main-search-dropzone", "drop", {"dataTransfer": data_transfer})
    wait_for_results(page)
    chip = page.inner_text("#tr-file-chip-container")
    assert "dropped_on_search_box.txt" in chip
    assert "1 locus" in chip
    assert len(page.query_selector_all("#results-table tbody tr")) > 0
    assert real_errors(console_errors) == [], f"JS errors: {real_errors(console_errors)}"


def test_upload_unknown_contig_does_not_break_panel(index_page, console_errors, tmp_path):
    """A token on an out-of-catalog contig (e.g. chr23) must match nothing rather than emitting
    invalid SQL that errors the whole panel; valid loci in the same file still return results."""
    page = index_page
    f = tmp_path / "panel.txt"
    f.write_text("chr23:100-200\nchr4:3074000-3075000\n")  # chr23 unknown; HTT region valid
    upload_files(page, [f])
    wait_for_results(page)
    assert not page.is_visible("#error-message")
    assert len(page.query_selector_all("#results-table tbody tr")) > 0
    assert real_errors(console_errors) == [], f"JS errors: {real_errors(console_errors)}"


def test_upload_bed(index_page, console_errors, tmp_path):
    """A BED file's regions are matched as 0-based half-open overlaps."""
    page = index_page
    f = tmp_path / "panel.bed"
    f.write_text("chr4\t3074000\t3075000\n")

    upload_files(page, [f])
    wait_for_results(page)
    assert "1 region" in page.inner_text("#tr-file-chip-container")
    assert len(page.query_selector_all("#results-table tbody tr")) > 0
    assert real_errors(console_errors) == [], f"JS errors: {real_errors(console_errors)}"


def test_upload_gzip(index_page, tmp_path):
    """A gzip-compressed .txt.gz is transparently decompressed."""
    page = index_page
    f = tmp_path / "panel.txt.gz"
    with gzip.open(f, "wt") as fh:
        fh.write("chr4:3074000-3075000\n")

    upload_files(page, [f])
    wait_for_results(page)
    assert len(page.query_selector_all("#results-table tbody tr")) > 0


def test_upload_after_paging_shows_results(index_page, tmp_path):
    """Regression: uploading a panel after navigating to page >=2 must reset to page 1 and still
    show the panel's matches, not an empty page (currentPage was left stale before the fix)."""
    page = index_page
    do_search(page, "ATXN1")
    if get_total_pages(page) > 1:
        click_next_page(page)

    f = tmp_path / "panel.txt"
    f.write_text("chr4:3074000-3075000\n")  # HTT region: far fewer rows than a page-2 offset
    upload_files(page, [f])
    wait_for_results(page)
    # Uploading clears the leftover "ATXN1" from the search box (panel is the query) and resets to
    # page 1, so the panel's matches are shown instead of an empty page-2 offset.
    assert page.input_value("#search-query") == ""
    assert len(page.query_selector_all("#results-table tbody tr")) > 0


def test_upload_token_is_sql_injection_safe(index_page, console_errors, tmp_path):
    """A backslash-quote payload must stay contained in the SQL literal, not inject `OR TRUE`.

    If the escaping were unsafe, the payload would break out of the string literal and match the
    whole catalog (a full page of rows); safely escaped it matches nothing.
    """
    page = index_page
    f = tmp_path / "evil.txt"
    f.write_text("x\\') OR TRUE--\n")
    page.set_input_files("#tr-file-input", str(f))
    page.wait_for_selector("#tr-file-chip-container .ui.label", timeout=10000)
    page.wait_for_load_state("networkidle")

    assert not page.is_visible("#error-message")
    # Injection would return a full page of unfiltered rows; a safely-escaped token matches nothing.
    assert len(page.query_selector_all("#results-table tbody tr")) < 50
    assert real_errors(console_errors) == [], f"JS errors: {real_errors(console_errors)}"


def test_upload_bed_empty_interval_ignored(index_page, tmp_path):
    """A zero-length BED interval (start == end) is skipped rather than matching spanning loci."""
    page = index_page
    f = tmp_path / "empty.bed"
    f.write_text("chr4\t3074000\t3074000\n")  # empty half-open interval -> no regions parsed
    # Open the modal first so its inline error message is visible when the parse yields nothing.
    page.click("#tr-file-upload-button")
    page.wait_for_selector("#tr-file-upload-modal", state="visible", timeout=5000)
    page.set_input_files("#tr-file-input", str(f))
    # No valid regions -> the "Nothing to search" error is shown and no chip appears.
    page.wait_for_selector("#tr-file-upload-error", state="visible", timeout=10000)
    assert not page.is_visible("#tr-file-chip-container .ui.label")


def test_remove_chip_clears_filter(index_page, tmp_path):
    """Removing the chip clears the uploaded filter."""
    page = index_page
    f = tmp_path / "panel.txt"
    f.write_text("chr4:3074000-3075000\n")

    upload_files(page, [f])
    wait_for_results(page)
    page.click("#tr-file-chip-remove")
    page.wait_for_selector("#tr-file-chip-container .ui.label", state="hidden", timeout=5000)
    assert not page.is_visible("#tr-file-chip-container .ui.label")
