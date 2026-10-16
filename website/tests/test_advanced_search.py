"""Tests for advanced search and filter functionality.

The custom filter tests load filters through the URL's "cf" parameter, where each filter is
"[!]column|operator|URL-encoded value|AND or OR" ('!' marks a disabled filter), filters are joined by ';', and an 'AND'
entry starts the next group of filters.
"""

import urllib.parse

import pytest


@pytest.fixture(autouse=True)
def fail_on_uncaught_javascript_errors(page):
    """Fails the test if the page throws an uncaught JavaScript error while it runs."""
    errors = []
    page.on("pageerror", lambda error: errors.append(str(error)))
    yield
    assert errors == []


GET_CUSTOM_FILTER_ROWS_SHOWN_PER_GROUP_JS = """() =>
    [...document.querySelectorAll('#custom-filter-rows .custom-filter-group')].map(group =>
        [...group.querySelectorAll('.custom-filter-row')].map(row => ({
            column: $(row).find('.custom-filter-field').dropdown('get value'),
            operator: $(row).find('.custom-filter-op').dropdown('get value'),
            value: $(row).find('.custom-filter-value').val(),
            enabled: $(row).find('.custom-filter-enable-checkbox input').is(':checked'),
            conjunction: $(row).find('.custom-filter-conjunction').dropdown('get value'),
        })))
"""

GET_ADD_CUSTOM_FILTER_LINKS_STATE_JS = """() => {
    const container = document.getElementById('add-custom-filter-link-container')
    let position = 'elsewhere'
    if (container.parentElement.id === 'source-catalog-and-add-custom-filter-link-row') {
        position = 'right of catalog selector'
    } else if (container.parentElement === [...document.querySelectorAll('#custom-filter-rows .custom-filter-group-links')].pop()) {
        position = "after the last group's link"
    }
    return {
        position,
        addCustomFilterLinkVisible: $('#add-custom-filter-link').is(':visible'),
        groupLinkVisible: $('#add-custom-filter-group-link').is(':visible'),
        helpPopupIsSetUp: $(container).find('.question').data('module-popup') !== undefined,
        linkTextPerGroup: [...document.querySelectorAll('#custom-filter-rows .custom-filter-group-links > a')].map(a => a.textContent),
    }
}"""


def open_index_page_with_custom_filters(page, base_url, custom_filters=""):
    """Opens the index page with the advanced search form shown and the given custom filters in the URL.

    Args:
        page: Playwright page.
        base_url: URL of the local server that serves the website.
        custom_filters: value for the URL's "cf" parameter.
    """
    url_params = {"showAdvFilters": 1}
    if custom_filters:
        url_params["cf"] = custom_filters
    page.goto(f"{base_url}/index.html#{urllib.parse.urlencode(url_params)}")
    # The page sets the cursor back to default once it has applied the URL's state to the form
    page.wait_for_function("() => document.body.style.cursor === 'default'")


def get_custom_filter_groups(page):
    """Returns the page's custom filter model (the customFilterGroups list)."""
    return page.evaluate("() => customFilterGroups")


def get_custom_filters_in_url(page):
    """Returns the value of the "cf" parameter in the page's URL."""
    return page.evaluate("() => new URLSearchParams(location.hash.slice(1)).get('cf')")


def get_custom_filter_container_children(page):
    """Returns the class of each element in #custom-filter-rows, which shows where the "AND" separators are."""
    return page.evaluate("() => [...document.getElementById('custom-filter-rows').children].map(el => el.className)")


def assert_custom_filter_rows_match_model(page):
    """Checks that the custom filter rows on the page show exactly what the model holds."""
    expected_rows_per_group = [
        [{**custom_filter, "conjunction": group["conjunction"]} for custom_filter in group["filters"]]
        for group in get_custom_filter_groups(page)
    ]
    assert page.evaluate(GET_CUSTOM_FILTER_ROWS_SHOWN_PER_GROUP_JS) == expected_rows_per_group


def add_custom_filters_with_dialog(page, link_selector, column_names):
    """Clicks a link that opens the custom filter column dialog, checks the given columns, and clicks Add Filters."""
    page.click(link_selector)
    page.wait_for_selector("#custom-filter-column-dialog", state="visible")
    assert page.locator("#custom-filter-column-checkboxes input:checked").count() == 0
    for column_name in column_names:
        page.click(f'#custom-filter-column-checkboxes input[data-column-name="{column_name}"] + label')
    page.click("#custom-filter-dialog-add-button")
    page.wait_for_selector("#custom-filter-column-dialog", state="hidden")


def get_add_another_custom_filter_link_selector(group_index):
    """Returns a selector for the "Add another custom filter..." link at the bottom of one custom filter group."""
    return f"#custom-filter-rows .custom-filter-group-links > a >> nth={group_index}"


def get_custom_filter_row(page, group_index, row_index):
    """Returns a locator for one custom filter row."""
    return page.locator("#custom-filter-rows .custom-filter-group").nth(group_index).locator(".custom-filter-row").nth(row_index)


def select_dropdown_item(row, dropdown_class, value):
    """Opens one of a row's dropdowns and clicks the item with the given value."""
    row.locator(f".{dropdown_class}").click()
    row.locator(f'.{dropdown_class} .item[data-value="{value}"]').click()


def capture_sql_of_searches(page):
    """Intercepts the page's queries to the server, so no search reaches BigQuery, and returns the list their SQL is
    appended to."""
    sql_sent = []

    def capture_sql_and_abort(route):
        sql_sent.append(route.request.post_data_json["sql"])
        route.abort()

    page.route("**/query_db", capture_sql_and_abort)
    return sql_sent


def test_custom_filter_links_move_below_the_filters_and_back(page, base_url):
    links_state_with_no_filters = {
        "position": "right of catalog selector", "addCustomFilterLinkVisible": True, "groupLinkVisible": False,
        "helpPopupIsSetUp": True, "linkTextPerGroup": [],
    }
    open_index_page_with_custom_filters(page, base_url)
    assert page.evaluate(GET_ADD_CUSTOM_FILTER_LINKS_STATE_JS) == links_state_with_no_filters

    add_custom_filters_with_dialog(page, "#add-custom-filter-link", ["HPRC256_Stdev"])
    assert page.evaluate(GET_ADD_CUSTOM_FILTER_LINKS_STATE_JS) == {
        "position": "after the last group's link", "addCustomFilterLinkVisible": False, "groupLinkVisible": True,
        "helpPopupIsSetUp": True, "linkTextPerGroup": ["Add another custom filter..."],
    }

    add_custom_filters_with_dialog(page, "#add-custom-filter-group-link", ["TenK10K_Stdev"])
    assert page.evaluate(GET_ADD_CUSTOM_FILTER_LINKS_STATE_JS) == {
        "position": "after the last group's link", "addCustomFilterLinkVisible": False, "groupLinkVisible": True,
        "helpPopupIsSetUp": True, "linkTextPerGroup": ["Add another custom filter...", "Add another custom filter..."],
    }

    page.click("#custom-filter-rows .custom-filter-delete >> nth=1")
    page.click("#custom-filter-rows .custom-filter-delete")
    assert get_custom_filter_groups(page) == []
    assert page.evaluate(GET_ADD_CUSTOM_FILTER_LINKS_STATE_JS) == links_state_with_no_filters

    # The links still work after being moved around
    add_custom_filters_with_dialog(page, "#add-custom-filter-link", ["HPRC256_Stdev"])
    assert [group["filters"][0]["column"] for group in get_custom_filter_groups(page)] == ["HPRC256_Stdev"]


def test_link_from_before_filter_groups_loads_as_one_group(page, base_url):
    open_index_page_with_custom_filters(
        page, base_url, "AoU1027_Stdev|>=|0.2|OR;!HPRC256_Stdev|>=|0.2|OR;ReferenceRepeatPurity|||OR")
    assert get_custom_filter_groups(page) == [{"conjunction": "OR", "filters": [
        {"column": "AoU1027_Stdev", "operator": ">=", "value": "0.2", "enabled": True},
        {"column": "HPRC256_Stdev", "operator": ">=", "value": "0.2", "enabled": False},
        {"column": "ReferenceRepeatPurity", "operator": "", "value": "", "enabled": True},
    ]}]
    assert_custom_filter_rows_match_model(page)
    assert get_custom_filter_container_children(page) == ["custom-filter-group"]


def test_add_group_of_filters_with_and(page, base_url):
    open_index_page_with_custom_filters(page, base_url, "AoU1027_Stdev|>=|0.2|OR;HPRC256_Stdev|>=|0.2|OR")
    first_group = get_custom_filter_groups(page)[0]

    # The dialog adds the checked columns in the order it lists them
    add_custom_filters_with_dialog(page, "#add-custom-filter-group-link", ["TenK10K_Stdev", "ReferenceRepeatPurity"])
    assert get_custom_filter_groups(page) == [first_group, {"conjunction": "AND", "filters": [
        {"column": "ReferenceRepeatPurity", "operator": "", "value": "", "enabled": True},
        {"column": "TenK10K_Stdev", "operator": "", "value": "", "enabled": True},
    ]}]
    assert_custom_filter_rows_match_model(page)
    assert get_custom_filter_container_children(page) == [
        "custom-filter-group", "custom-filter-group-separator", "custom-filter-group"]
    assert page.locator(".custom-filter-group-separator").inner_text() == "AND"
    assert get_custom_filters_in_url(page) == (
        "AoU1027_Stdev|>=|0.2|OR;HPRC256_Stdev|>=|0.2|OR;AND;ReferenceRepeatPurity|||AND;TenK10K_Stdev|||AND")



def test_each_groups_add_another_custom_filter_link_adds_to_that_group(page, base_url):
    open_index_page_with_custom_filters(page, base_url, "AoU1027_Stdev|>=|0.2|OR;AND;TenK10K_Stdev|<=|5|AND")
    add_custom_filters_with_dialog(page, get_add_another_custom_filter_link_selector(0), ["GencodeGeneName"])
    add_custom_filters_with_dialog(page, get_add_another_custom_filter_link_selector(1), ["ReferenceMotif"])
    assert [[f["column"] for f in group["filters"]] for group in get_custom_filter_groups(page)] == [
        ["AoU1027_Stdev", "GencodeGeneName"], ["TenK10K_Stdev", "ReferenceMotif"]]
    assert_custom_filter_rows_match_model(page)
    # A filter added to an OR group shows the group's OR
    assert get_custom_filter_row(page, 0, 1).locator(".custom-filter-conjunction .text").inner_text() == "OR"


def test_new_group_dialog_with_nothing_checked_adds_no_group(page, base_url):
    open_index_page_with_custom_filters(page, base_url, "AoU1027_Stdev|>=|0.2|AND")
    groups_before = get_custom_filter_groups(page)
    add_custom_filters_with_dialog(page, "#add-custom-filter-group-link", [])
    assert get_custom_filter_groups(page) == groups_before
    assert get_custom_filter_container_children(page) == ["custom-filter-group"]


def test_and_or_selectors_are_linked_within_a_group_only(page, base_url):
    open_index_page_with_custom_filters(
        page, base_url,
        "AoU1027_Stdev|>=|0.2|AND;HPRC256_Stdev|>=|0.2|AND;AND;TenK10K_Stdev|<=|5|AND;ReferenceRepeatPurity|>=|0.5|AND")
    select_dropdown_item(get_custom_filter_row(page, 1, 1), "custom-filter-conjunction", "OR")
    assert [group["conjunction"] for group in get_custom_filter_groups(page)] == ["AND", "OR"]
    assert_custom_filter_rows_match_model(page)
    assert get_custom_filters_in_url(page) == (
        "AoU1027_Stdev|>=|0.2|AND;HPRC256_Stdev|>=|0.2|AND;AND;TenK10K_Stdev|<=|5|OR;ReferenceRepeatPurity|>=|0.5|OR")


def test_deleting_the_last_filter_of_a_group_deletes_the_group(page, base_url):
    open_index_page_with_custom_filters(
        page, base_url, "AoU1027_Stdev|>=|0.2|AND;AND;HPRC256_Stdev|>=|0.2|OR;AND;TenK10K_Stdev|<=|5|AND")
    get_custom_filter_row(page, 1, 0).locator(".custom-filter-delete").click()
    assert [group["filters"][0]["column"] for group in get_custom_filter_groups(page)] == [
        "AoU1027_Stdev", "TenK10K_Stdev"]
    assert get_custom_filter_container_children(page) == [
        "custom-filter-group", "custom-filter-group-separator", "custom-filter-group"]
    assert get_custom_filters_in_url(page) == "AoU1027_Stdev|>=|0.2|AND;AND;TenK10K_Stdev|<=|5|AND"

    get_custom_filter_row(page, 0, 0).locator(".custom-filter-delete").click()
    assert get_custom_filter_container_children(page) == ["custom-filter-group"]
    assert get_custom_filters_in_url(page) == "TenK10K_Stdev|<=|5|AND"
    # The AND/OR selector is hidden on the first row of a group
    assert not get_custom_filter_row(page, 0, 0).locator(".custom-filter-conjunction").is_visible()
    assert_custom_filter_rows_match_model(page)


def test_editing_a_filter_row_updates_the_model_and_url(page, base_url):
    open_index_page_with_custom_filters(page, base_url, "AoU1027_Stdev|>=|0.2|AND;AND;TenK10K_Stdev|||AND")
    row = get_custom_filter_row(page, 1, 0)
    select_dropdown_item(row, "custom-filter-op", "<=")
    row.locator(".custom-filter-value").fill("5")
    row.locator(".custom-filter-value").press("Tab")
    get_custom_filter_row(page, 0, 0).locator(".custom-filter-enable-checkbox label").click()
    assert get_custom_filters_in_url(page) == "!AoU1027_Stdev|>=|0.2|AND;AND;TenK10K_Stdev|<=|5|AND"
    assert_custom_filter_rows_match_model(page)

    # A null check hides the value input, and changing the column clears the operator and shows the input again
    select_dropdown_item(row, "custom-filter-op", "IS NULL")
    assert not row.locator(".custom-filter-value").is_visible()
    select_dropdown_item(row, "custom-filter-field", "GencodeGeneName")
    assert get_custom_filter_groups(page)[1]["filters"] == [
        {"column": "GencodeGeneName", "operator": "", "value": "5", "enabled": True}]
    assert row.locator(".custom-filter-value").is_visible()
    assert row.locator('.custom-filter-op .item[data-value="LIKE"]').count() == 1
    assert_custom_filter_rows_match_model(page)


def test_custom_filter_url_format_round_trip(page, base_url):
    open_index_page_with_custom_filters(page, base_url)
    serialized = "AoU1027_Stdev|>=|0.2|OR;!GencodeGeneName|=|a%7Cb%3BAND|OR;AND;ReferenceMotif|IS NULL||AND"
    assert page.evaluate("s => serializeCustomFilterGroups(parseCustomFilterGroups(s))", serialized) == serialized
    assert page.evaluate("s => parseCustomFilterGroups(s)", serialized)[0]["filters"][1]["value"] == "a|b;AND"

    # Stray 'AND' entries don't create empty groups, and an operator the column's type doesn't offer is dropped
    assert page.evaluate("s => parseCustomFilterGroups(s)", "AND;HPRC256_Stdev|LIKE|5|AND;AND;AND") == [
        {"conjunction": "AND", "filters": [{"column": "HPRC256_Stdev", "operator": "", "value": "5", "enabled": True}]}]


def test_search_joins_custom_filter_groups_with_and(page, base_url):
    sql_sent = capture_sql_of_searches(page)
    open_index_page_with_custom_filters(
        page, base_url,
        "AoU1027_Stdev|>=|0.2|OR;HPRC256_Stdev|>=|0.2|OR;AND;TenK10K_Stdev|<=|5|AND;"
        "AND;!ReferenceRepeatPurity|>=|0.5|AND;AND;GencodeGeneName|LIKE|SCN%25|AND")
    with page.expect_request("**/query_db"):
        page.click("#search-button")
    assert len(sql_sent) == 1
    assert ("(AoU1027_Stdev >= 0.2 OR HPRC256_Stdev >= 0.2) AND TenK10K_Stdev <= 5 AND GencodeGeneName LIKE 'SCN%'"
            in sql_sent[0])
    # A group whose filters are all disabled adds no condition
    assert "ReferenceRepeatPurity >=" not in sql_sent[0]


def test_search_reports_an_incomplete_filter_in_a_later_group(page, base_url):
    sql_sent = capture_sql_of_searches(page)
    open_index_page_with_custom_filters(
        page, base_url, "AoU1027_Stdev|>=|0.2|AND;AND;TenK10K_Stdev|<=|5|AND;HPRC256_Stdev|||AND")
    page.click("#search-button")
    assert "Please select an operator for the HPRC256:LPS Stdev custom filter" in page.inner_text("#form-error-message")
    rows_marked_as_errors = page.evaluate("""() => [...document.querySelectorAll('#custom-filter-rows .custom-filter-group')]
        .flatMap((group, groupIndex) => [...group.querySelectorAll('.custom-filter-row')]
            .map((row, rowIndex) => row.classList.contains('error') ? [groupIndex, rowIndex] : null)
            .filter(Boolean))""")
    assert rows_marked_as_errors == [[1, 1]]
    assert sql_sent == []
