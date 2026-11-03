"""Tests for basic search functionality on the index page.

Each example query type is tested: genomic region, gene name, Ensembl gene ID,
and Ensembl transcript ID. Tests confirm that results load without JS errors.
After basic search validation, all columns are enabled and each is sorted
ascending and descending. Multi-page queries also test pagination.
"""

import re

import pytest


# The tests query the live site's BigQuery table through its public proxy, so every query counts toward the site's
# daily query quota. A region search scans only the table partition and blocks its region falls in (tens of MB), while a
# gene or transcript name search scans whole columns of the table (about 0.5 GiB). So the tests search names only where
# name search itself is tested, and search a region everywhere else.

# Queries from the example links on the page, with expected minimum result counts. Each one is searched once.
SEARCH_QUERIES = [
    ("chr4:3074000-3075000", 1),
    ("ATXN1", 1),
    ("ENSG00000165060", 1),
    ("ENST00000370475", 1),
]

# A region search with several pages of results: the span of the ATXN1 gene
MULTI_PAGE_REGION_QUERY = "chr6:16299112-16761490"

# Columns whose sorting the tests check, one of each kind of sort: by position, by a text column, by a column with
# its own sort columns (the allele histogram sorts by its stdev), and by a number
COLUMNS_TO_TEST_SORTING = ["ReferenceRegion", "GencodeGeneName", "HPRC256_AlleleHistogram", "AoU1027_OE_LengthPercentile"]


def do_search(page, query):
    """Fill the search box, click Search, and wait for results."""
    page.fill("#search-query", query)
    page.click("#search-button")
    page.wait_for_selector("#results-table tbody tr", timeout=30000)


def wait_for_results_refresh(page):
    """Wait for a sort, page change or column change, which re-runs the search on screen, to finish.

    The old rows stay visible while the new query runs, so this waits for the page to set the cursor back to default,
    which it does when a search finishes. The search sets it to 'wait' before the action's click returns.
    """
    page.wait_for_function("() => document.body.style.cursor === 'default'", timeout=30000)
    page.wait_for_selector("#results-table tbody tr", timeout=30000)


def real_errors(console_errors):
    """Return console errors excluding blocked-CDN / analytics network failures.

    Google Analytics / googletagmanager and other third-party CDNs are often unreachable in
    CI/sandbox environments, surfacing as "Failed to load resource: net::ERR_*" console errors
    that are unrelated to the page's own logic. Real JS exceptions do not match these patterns
    and are still caught.
    """
    return [
        e for e in console_errors
        if "failed to load resource" not in e.lower()
        and "google-analytics" not in e.lower()
        and "googletagmanager" not in e.lower()
        and "net::err" not in e.lower()
    ]


def assert_no_errors(page, console_errors, context):
    """Assert no error message is shown and no (non-environmental) JS errors occurred."""
    assert not page.is_visible("#error-message"), (
        f"Error message displayed ({context}): {page.inner_text('#error-text')}"
    )
    assert real_errors(console_errors) == [], (
        f"JavaScript errors ({context}): {real_errors(console_errors)}"
    )


def enable_all_columns(page):
    """Open column settings dialog, check all boxes, and click Apply."""
    page.click("#results-table-settings-button")
    page.wait_for_selector("#results-table-column-settings-dialog", state="visible", timeout=5000)

    for cb in page.query_selector_all(
        ".results-table-column-settings-checkboxes input[type='checkbox']"
    ):
        if not cb.is_checked():
            cb.check(force=True)

    page.click("#results-table-column-settings-apply-button")
    wait_for_results_refresh(page)


def get_sortable_column_names(page):
    """Return list of column data-column attribute values for all sortable headers."""
    return [
        h.get_attribute("data-column")
        for h in page.query_selector_all("#results-table th.results-table-sortable-header")
    ]


def get_total_pages(page):
    """Parse total page count from the pagination 'of N' text."""
    text = page.locator(".bottom-pagination .total-pages").inner_text()
    match = re.search(r"of\s+([\d,]+)", text)
    return int(match.group(1).replace(",", "")) if match else 1


def click_next_page(page):
    """Click the next-page button using bottom pagination and wait for refresh."""
    page.locator(".bottom-pagination .next-page").click(force=True)
    wait_for_results_refresh(page)


def click_last_page(page):
    """Click the last-page button using bottom pagination and wait for refresh."""
    page.locator(".bottom-pagination .last-page").click(force=True)
    wait_for_results_refresh(page)


@pytest.mark.parametrize("query,min_results", SEARCH_QUERIES)
def test_search_returns_results(index_page, console_errors, query, min_results):
    """Search for a query and verify results appear with no JS errors."""
    do_search(index_page, query)

    rows = index_page.query_selector_all("#results-table tbody tr")
    assert len(rows) >= min_results, (
        f"Expected at least {min_results} result(s) for '{query}', got {len(rows)}"
    )
    assert_no_errors(index_page, console_errors, f"search '{query}'")


def test_search_chr4_region_has_htt_disease_locus(index_page, console_errors):
    """Search chr4:3074000-3075000 and verify it contains the HTT disease-associated locus."""
    query = "chr4:3074000-3075000"
    do_search(index_page, query)

    disease_rows = index_page.query_selector_all("#results-table tbody tr.disease-associated-row")
    assert len(disease_rows) == 1, (
        f"Expected exactly 1 disease-associated row, got {len(disease_rows)}"
    )

    row_text = disease_rows[0].inner_text()
    assert "HTT" in row_text, f"Expected gene 'HTT' in disease row, got: {row_text}"
    assert "4:3074876-3074933" in row_text, (
        f"Expected region '4:3074876-3074933' in disease row, got: {row_text}"
    )
    assert_no_errors(index_page, console_errors, f"search '{query}'")


def test_all_columns_and_sorting(index_page, console_errors):
    """Enable all columns, verify they load, then sort a few columns of different kinds in both directions."""
    page = index_page
    query = SEARCH_QUERIES[0][0]
    do_search(page, query)
    enable_all_columns(page)
    assert_no_errors(page, console_errors, f"all columns for '{query}'")
    assert set(COLUMNS_TO_TEST_SORTING) <= set(get_sortable_column_names(page))

    for column_name in COLUMNS_TO_TEST_SORTING:
        selector = f"th.results-table-sortable-header[data-column='{column_name}']"

        # The first click sorts in the column's default direction, and the second in the other direction
        page.click(selector)
        wait_for_results_refresh(page)
        assert_no_errors(page, console_errors, f"sort '{column_name}' click 1 for '{query}'")
        first_direction = page.evaluate("GLOBAL_PAGE_STATE['sd']")

        page.click(selector)
        wait_for_results_refresh(page)
        assert_no_errors(page, console_errors, f"sort '{column_name}' click 2 for '{query}'")
        assert page.evaluate("[GLOBAL_PAGE_STATE['sc'], GLOBAL_PAGE_STATE['sd']]") == [
            column_name, "DESC" if first_direction == "ASC" else "ASC"]


def test_pagination_all_columns(index_page, console_errors):
    """Enable all columns, navigate to pages 2, 3, and the last page without errors."""
    page = index_page
    query = MULTI_PAGE_REGION_QUERY
    do_search(page, query)
    enable_all_columns(page)
    assert_no_errors(page, console_errors, f"all columns page 1 for '{query}'")

    # Verify there are multiple pages
    total_pages = get_total_pages(page)
    assert total_pages > 1, (
        f"Expected more than 1 page of results for '{query}', got {total_pages}"
    )

    # Navigate to page 2
    click_next_page(page)
    assert_no_errors(page, console_errors, f"all columns page 2 for '{query}'")

    # Navigate to page 3 (if it exists)
    if total_pages >= 3:
        click_next_page(page)
        assert_no_errors(page, console_errors, f"all columns page 3 for '{query}'")

    # Jump to last page (if not already there)
    if total_pages > 3:
        click_last_page(page)
        assert_no_errors(page, console_errors, f"all columns last page ({total_pages}) for '{query}'")
