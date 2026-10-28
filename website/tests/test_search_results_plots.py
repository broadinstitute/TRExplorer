"""Tests for the search results plots ("Add plot..." below the results table).

The first tests call the plot code's functions in the page. The last ones add, reload, edit and remove plots
through the page, so they query BigQuery through the live query proxy like the other tests here.
"""

import re

from test_basic_search import do_search, real_errors


def evaluate_with_settings(page, expression, settings):
    """Evaluates expression with `s` set to DEFAULT_SEARCH_RESULTS_PLOT_SETTINGS overridden by settings."""
    return page.evaluate(f"(c) => {{ const s = {{...DEFAULT_SEARCH_RESULTS_PLOT_SETTINGS, ...c}}; return {expression} }}", settings)


def wait_for_plots_to_finish_loading(page, num_plots):
    page.wait_for_function(
        f"document.querySelectorAll('.search-results-plot').length === {num_plots} && "
        "[...document.querySelectorAll('.search-results-plot-loading')].every(e => e.style.display === 'none')",
        timeout=90000)


def select_dropdown_item(page, dropdown_id, value):
    page.click(f"#{dropdown_id}")
    page.click(f"#{dropdown_id} .menu .item[data-value='{value}']")


def test_plot_settings_survive_the_url_round_trip(index_page):
    settings = {"type": "box", "x": "GencodeGeneRegion", "y": "HPRC256_Stdev", "yScale": "symlog", "color": "Source",
                "topN": 5, "sort": "MotifSize", "dir": "asc", "missing": 1}
    serialized = evaluate_with_settings(index_page, "serializeSearchResultsPlotSettingsList([s])", settings)
    assert serialized == "t:box,x:GencodeGeneRegion,y:HPRC256_Stdev,c:Source,ys:symlog,n:5,m:1,s:MotifSize,d:asc"
    assert evaluate_with_settings(
        index_page, "JSON.stringify(parseSearchResultsPlotSettingsList(serializeSearchResultsPlotSettingsList([s]))[0]) === JSON.stringify(s)",
        settings)


def test_plots_from_the_url_that_cant_be_plotted_are_dropped(index_page):
    parsed = index_page.evaluate("parseSearchResultsPlotSettingsList('t:bar,x:NoSuchColumn;t:scatter,x:GencodeGeneRegion,y:HPRC256_Stdev;t:bar,x:Source')")
    assert [p["x"] for p in parsed] == ["Source"]


def test_problems_with_plot_settings(index_page):
    assert evaluate_with_settings(index_page, "findProblemWithSearchResultsPlotSettings(s)", {}) == "Select an X-axis column."
    assert evaluate_with_settings(index_page, "findProblemWithSearchResultsPlotSettings(s)", {"type": "box", "x": "Source"}) == "Select a Y-axis column."
    assert "isn't numeric" in evaluate_with_settings(index_page, "findProblemWithSearchResultsPlotSettings(s)", {"type": "bar", "x": "Source", "y": "GencodeGeneRegion"})
    assert evaluate_with_settings(index_page, "findProblemWithSearchResultsPlotSettings(s)", {"type": "bar", "x": "Source"}) == ""


def test_coordinates_can_only_be_scatter_plot_y_axes(index_page):
    for plot_type in ["bar", "box"]:
        assert evaluate_with_settings(index_page, "findProblemWithSearchResultsPlotSettings(s)", {"type": plot_type, "x": "Source", "y": "start_0based"}) == (
            "Start (0-based) can only be a scatter plot Y-axis.")
    assert evaluate_with_settings(index_page, "findProblemWithSearchResultsPlotSettings(s)", {"type": "scatter", "x": "MotifSize", "y": "start_0based"}) == ""

    index_page.evaluate("showSearchResultsPlotDialog(null)")
    index_page.wait_for_selector("#search-results-plot-dialog", state="visible")
    y_menu_offers_start_and_end = """() => ['start_0based', 'end_1based'].map(
        name => !!document.querySelector(`#plot-y-column-dropdown .menu .item[data-value="${name}"]`))"""
    for plot_type, offered in [("bar", False), ("box", False), ("scatter", True)]:
        index_page.click(f"#search-results-plot-dialog .ui.radio.checkbox:has(input[value='{plot_type}'])")
        assert index_page.evaluate(y_menu_offers_start_and_end) == [offered, offered], plot_type


def test_plot_queries_use_the_search_where_clause(index_page):
    where_clause = "WHERE (UPPER(GencodeGeneName) = 'ATXN1')"
    for settings in [
        {"type": "bar", "x": "GencodeGeneRegion", "color": "Source"},
        {"type": "bar", "x": "MotifSize", "xMode": "bins", "xScale": "log"},
        {"type": "scatter", "x": "MotifSize", "y": "HPRC256_Stdev"},
        {"type": "box", "x": "GencodeGeneRegion", "y": "HPRC256_Stdev", "sort": "MotifSize"},
    ]:
        sql = index_page.evaluate(
            "([c, w]) => buildSearchResultsPlotQuery({...DEFAULT_SEARCH_RESULTS_PLOT_SETTINGS, ...c}, COMPLETE_TABLE_ID, w)",
            [settings, where_clause])
        assert where_clause in sql
        assert sql.startswith("WITH base AS (SELECT ")


def test_symlog_axis_labels_prefer_powers_of_ten_and_dont_crowd(index_page):
    assert index_page.evaluate("computeAxisTickValuesForScale(0, 100, 'symlog')") == [0, 1, 2, 5, 10, 20, 50, 100]
    assert index_page.evaluate("computeAxisTickValuesForScale(-10, 10, 'symlog')") == [-10, -5, -2, -1, 0, 1, 2, 5, 10]
    ticks = index_page.evaluate("computeAxisTickValuesForScale(0, 80, 'symlog')")
    assert 1 in ticks and 10 in ticks


def test_log_axis_labels_cover_small_values(index_page):
    assert index_page.evaluate("computeAxisTickValuesForScale(1e-10, 1, 'log')") == [1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1]


def test_plot_numbers_keep_integers_exact_and_small_values_nonzero(index_page):
    assert index_page.evaluate("[formatPlotNumber(94418421), formatPlotNumber(1e-180), formatPlotNumber(0.000123), formatPlotNumber(1234.567)]") == [
        "94,418,421", "1.00e-180", "0.000123", "1,230"]


def test_sorting_text_categories_by_value_puts_9_before_10(index_page):
    sorted_labels = evaluate_with_settings(
        index_page, "sortPlotCategories(s, ['10', '9', 'X', '2'], {}, {})", {"type": "bar", "x": "chrom", "sort": "value", "dir": "asc"})
    assert sorted_labels == ["2", "9", "10", "X"]


def test_add_plot_link_appears_for_any_number_of_results(index_page):
    do_search(index_page, "chr4:3074876-3074933")
    assert 0 < index_page.evaluate("window.searchResults.total_results") < 5
    index_page.wait_for_selector("#add-search-results-plot-link", state="visible", timeout=10000)
    # The first 10 kb of chr1 are all N's, so no TR loci are there
    index_page.fill("#search-query", "chr1:1-5000")
    index_page.click("#search-button")
    index_page.wait_for_function("window.searchResults.total_results === 0", timeout=30000)
    assert not index_page.is_visible("#add-search-results-plot-link")


def test_a_bad_value_in_a_hidden_field_does_not_block_ok(index_page):
    do_search(index_page, "ATXN1")
    index_page.click("#add-search-results-plot-link")
    index_page.wait_for_selector("#search-results-plot-dialog", state="visible")
    select_dropdown_item(index_page, "plot-x-column-dropdown", "Source")
    index_page.fill("#plot-top-n-input", "200")
    # A numerical-only x column makes the bar plot a histogram, which hides the number of categories
    select_dropdown_item(index_page, "plot-x-column-dropdown", "start_0based")
    assert not index_page.is_visible("#plot-top-n-field")
    index_page.click("#search-results-plot-dialog-ok-button")
    wait_for_plots_to_finish_loading(index_page, 1)
    assert index_page.inner_text(".search-results-plot-title") == "Distribution of Start (0-based)"


def click_plot_bar(page, dataset_index, bar_index):
    """Clicks the center of a bar (or box) in the first plot's chart and returns the value it shows."""
    bar = page.evaluate(f"""() => {{
        const chart = Chart.getChart(document.querySelector('.search-results-plot canvas'))
        const element = chart.getDatasetMeta({dataset_index}).data[{bar_index}]
        const {{x, y, base}} = element.getProps(['x', 'y', 'base'], true)
        return {{x, y: (y + base) / 2, value: chart.data.datasets[{dataset_index}].rawValues?.[{bar_index}]}}
    }}""")
    page.locator(".search-results-plot canvas").first.click(position={"x": bar["x"], "y": bar["y"]})
    return bar["value"]


def get_custom_filters(page):
    return page.evaluate("customFilterGroups.flatMap(g => g.filters.map(f => [f.column, f.operator, f.value]))")


def open_search_with_plots(page, query, serialized_plots):
    """Loads the page for a search with the given plots in its URL, and waits for the plots to be drawn."""
    url = f"{page.url.split('#')[0]}#showRs=1&q={query}&plots={serialized_plots.replace(':', '%3A').replace(',', '%2C')}"
    # Going through a blank page makes this a fresh page load, rather than a hash change on the current page
    page.goto("about:blank")
    page.goto(url)
    wait_for_plots_to_finish_loading(page, len(serialized_plots.split(';')))


def test_clicking_a_stacked_bar_segment_filters_to_its_loci(index_page, console_errors):
    open_search_with_plots(index_page, "ATXN1", "t:bar,x:GencodeGeneRegion,c:MotifSize,n:5")
    # The tallest bar is intron, and the first color is motif size 1
    num_loci_in_segment = click_plot_bar(index_page, 0, 0)
    index_page.wait_for_function(f"window.searchResults && window.searchResults.total_results === {num_loci_in_segment}", timeout=60000)
    assert get_custom_filters(index_page) == [["GencodeGeneRegion", "=", "intron"], ["MotifSize", "=", "1"]]
    assert index_page.evaluate("GLOBAL_PAGE_STATE.showAdvFilters") == 1
    assert "cf=" in index_page.url
    assert real_errors(console_errors) == []


def test_clicking_a_histogram_bin_filters_to_its_range(index_page):
    open_search_with_plots(index_page, "ATXN1", "t:bar,x:MotifSize,xm:bins,b:10")
    num_loci_in_second_bin = click_plot_bar(index_page, 0, 1)
    index_page.wait_for_function(f"window.searchResults && window.searchResults.total_results === {num_loci_in_second_bin}", timeout=60000)
    filters = get_custom_filters(index_page)
    assert [f[:2] for f in filters] == [["MotifSize", ">="], ["MotifSize", "<="]]


def get_plot_bar_index(page, display_label):
    return page.evaluate(f"Chart.getChart(document.querySelector('.search-results-plot canvas')).data.labels.indexOf({display_label!r})")


def test_clicking_a_bar_whose_label_has_a_quote_filters_to_its_loci(index_page, console_errors):
    open_search_with_plots(index_page, "ATXN1", "t:bar,x:GencodeGeneRegion")
    num_loci_in_bar = click_plot_bar(index_page, 0, get_plot_bar_index(index_page, "3' UTR"))
    index_page.wait_for_function(f"window.searchResults && window.searchResults.total_results === {num_loci_in_bar}", timeout=60000)
    assert get_custom_filters(index_page) == [["GencodeGeneRegion", "=", "3' UTR"]]
    assert real_errors(console_errors) == []


def test_clicking_a_bar_resets_settings_left_in_the_hidden_advanced_search_panel(index_page):
    # keepCodingLoci is set but doesn't apply, since the advanced search panel is hidden
    url = f"{index_page.url.split('#')[0]}#showRs=1&q=ATXN1&keepCodingLoci=1&plots=t%3Abar%2Cx%3AGencodeGeneRegion"
    index_page.goto("about:blank")
    index_page.goto(url)
    wait_for_plots_to_finish_loading(index_page, 1)
    num_loci_in_bar = click_plot_bar(index_page, 0, get_plot_bar_index(index_page, "intron"))
    index_page.wait_for_function(f"window.searchResults && window.searchResults.total_results === {num_loci_in_bar}", timeout=60000)
    assert index_page.evaluate("[GLOBAL_PAGE_STATE.showAdvFilters, GLOBAL_PAGE_STATE.keepCodingLoci]") == [1, 0]
    assert not index_page.is_checked("input[name='keepCodingLoci']")


def test_clicking_a_bar_ignores_unsubmitted_edits_and_back_hides_the_panel_again(index_page):
    open_search_with_plots(index_page, "ATXN1", "t:bar,x:GencodeGeneRegion")
    # An edit to the search box that wasn't searched for doesn't change what the click selects
    index_page.fill("#search-query", "HTT")
    index_page.locator("#search-query").blur()
    num_loci_in_bar = click_plot_bar(index_page, 0, get_plot_bar_index(index_page, "intron"))
    index_page.wait_for_function(f"window.searchResults && window.searchResults.total_results === {num_loci_in_bar}", timeout=60000)
    assert index_page.input_value("#search-query") == "ATXN1"
    assert index_page.is_visible("#advanced-search-form")
    # Back returns to the search before the click, whose advanced search panel was hidden
    index_page.go_back()
    index_page.wait_for_function("GLOBAL_PAGE_STATE.showAdvFilters === 0")
    index_page.wait_for_selector("#advanced-search-form", state="hidden")
    assert index_page.input_value("#search-query") == "ATXN1"


def test_clicking_a_fractional_histogram_bin_filters_to_exactly_its_loci(index_page):
    open_search_with_plots(index_page, "ATXN1", "t:bar,x:ReferenceRepeatPurity")
    num_loci_in_bin = click_plot_bar(index_page, 0, 9)
    index_page.wait_for_function(f"window.searchResults && window.searchResults.total_results === {num_loci_in_bin}", timeout=60000)


def test_clicking_a_box_keeps_only_loci_with_a_y_value(index_page):
    open_search_with_plots(index_page, "ATXN1", "t:box,x:GencodeGeneRegion,y:HPRC256_Stdev")
    index_page.wait_for_function("Chart.getChart(document.querySelector('.search-results-plot canvas'))")
    num_loci_in_box = index_page.evaluate(
        "Chart.getChart(document.querySelector('.search-results-plot canvas')).data.datasets[0].boxRows[0].n_y")
    click_plot_bar(index_page, 0, 0)
    index_page.wait_for_function(f"window.searchResults && window.searchResults.total_results === {num_loci_in_box}", timeout=60000)
    assert get_custom_filters(index_page)[-1] == ["HPRC256_Stdev", "IS NOT NULL", ""]


def test_a_box_with_no_height_can_be_clicked(index_page):
    # Most intron loci have a purity of 1, so the intron box's 25th and 75th percentiles are both 1
    open_search_with_plots(index_page, "ATXN1", "t:box,x:GencodeGeneRegion,y:ReferenceRepeatPurity")
    intron_index = get_plot_bar_index(index_page, "intron")
    assert index_page.evaluate(
        f"(() => {{ const row = Chart.getChart(document.querySelector('.search-results-plot canvas')).data.datasets[0].boxRows[{intron_index}]; return row.q1 === row.q3 }})()")
    click_plot_bar(index_page, 0, intron_index)
    index_page.wait_for_function("customFilterGroups.length === 1", timeout=60000)


def test_log_histogram_click_filters_of_an_integer_column_use_whole_numbers(index_page):
    # Bin 15 of 30 log bins over 1..25 starts at log10(5), whose inverse rounds to 5.000000000000001
    bins = evaluate_with_settings(
        index_page, "getHistogramBins(s, [{bin_range_lo: 0, bin_range_hi: Math.log10(25), bin_width: Math.log10(25) / 30}])",
        {"type": "bar", "x": "MotifSize", "xMode": "bins", "xScale": "log"})
    assert bins[15]["rawStart"] == 5
    assert all(isinstance(b["rawStart"], int) for b in bins)


def test_symlog_histogram_click_filter_edges_match_the_query(index_page):
    # Bin 15 of 30 symlog bins over 0..63 starts at log10(8) = symlog(7), so 7 belongs to bin 15
    bins = evaluate_with_settings(
        index_page, "getHistogramBins(s, [{bin_range_lo: 0, bin_range_hi: symlogTransform(63), bin_width: symlogTransform(63) / 30}])",
        {"type": "bar", "x": "NsInFlanks", "xMode": "bins", "xScale": "symlog"})
    assert bins[15]["rawStart"] == 7
    # The JavaScript transform matches the query's LOG10(1 + ABS(x)) exactly, and tiny values aren't rounded to 0
    assert index_page.evaluate("symlogTransform(7) === Math.log10(8) && symlogTransform(1e-180) > 0")
    # A non-integer column's edge at symlog(9) is 9, not the 9.000000000000002 that the inverse gives
    assert index_page.evaluate("roundedValueAtOrAfterPlottedEdge(1, 'symlog')") == 9


def test_percent_stacks_end_at_100(index_page):
    open_search_with_plots(index_page, "ATXN1", "t:bar,x:GencodeGeneRegion,c:MotifSize,p:1,n:4,o:1")
    assert index_page.evaluate("Chart.getChart(document.querySelector('.search-results-plot canvas')).scales.y.max") == 100


def test_last_integer_histogram_bin_label_stops_at_the_largest_value(index_page):
    bins = evaluate_with_settings(
        index_page, "getHistogramBins(s, [{bin_range_lo: 1, bin_range_hi: 200, bin_width: 7}])", {"type": "bar", "x": "MotifSize", "xMode": "bins"})
    assert bins[-1]["displayLabel"] == "197-200"


def test_other_and_missing_categories_have_no_click_filter(index_page):
    settings = {"type": "bar", "x": "GencodeGeneRegion", "color": "Source"}
    for x_label, color_label in [("__other__", "TRExplorerV2:HipSTRCatalog"), ("intron", "__missing__")]:
        assert evaluate_with_settings(
            index_page, f"getCustomFiltersForClickedPlotElement(s, {{xLabel: '{x_label}'}}, '{color_label}')", settings) is None
    # The Missing category of a numeric column is exactly the loci with no value
    numeric_missing = evaluate_with_settings(
        index_page, "getCustomFiltersForClickedPlotElement(s, {xLabel: 'intron'}, '__missing__')", {"type": "bar", "x": "GencodeGeneRegion", "color": "MotifSize"})
    assert [[f["column"], f["operator"]] for f in numeric_missing] == [["GencodeGeneRegion", "="], ["MotifSize", "IS NULL"]]
    # The first bin of a log histogram keeps only positive values, and has no lower bound otherwise
    first_log_bin = evaluate_with_settings(
        index_page,
        "getCustomFiltersForClickedPlotElement(s, {xLabel: '0', rawStart: 0.04, rawEnd: 0.050200000000000001, isFirstBin: true, isLastBin: false, isIntegerBins: false}, '')",
        {"type": "bar", "x": "HPRC256_Stdev", "xScale": "log"})
    assert [[f["operator"], f["value"]] for f in first_log_bin] == [[">", "0"], ["<", "0.0502"]]


def test_small_scatter_plot_is_drawn_as_points(index_page):
    # The heatmap query runs first, and finds few enough loci for the plot to show each one
    open_search_with_plots(index_page, "ATXN1", "t:scatter,x:MotifSize,y:HPRC256_Stdev")
    assert index_page.evaluate("Chart.getChart(document.querySelector('.search-results-plot canvas')).config.type") == "scatter"
    assert "Based on" in index_page.inner_text(".search-results-plot-caption")


def test_large_scatter_plot_is_drawn_as_a_heatmap(index_page):
    open_search_with_plots(index_page, "chr21:1-46700000", "t:scatter,x:HPRC256_Stdev,y:AoU1027_Stdev,xs:symlog,ys:symlog")
    assert index_page.evaluate("Chart.getChart(document.querySelector('.search-results-plot canvas')).config.type") == "matrix"
    assert "Each cell's color shows how many" in index_page.inner_text(".search-results-plot-caption")


def test_add_reload_edit_and_remove_plots(index_page, console_errors):
    do_search(index_page, "ATXN1")
    index_page.click("#add-search-results-plot-link")
    index_page.wait_for_selector("#search-results-plot-dialog", state="visible")

    # Clicking OK without an x-axis column keeps the dialog open and says why
    index_page.click("#search-results-plot-dialog-ok-button")
    index_page.wait_for_selector("#search-results-plot-dialog-error", state="visible")
    assert index_page.is_visible("#search-results-plot-dialog")

    select_dropdown_item(index_page, "plot-x-column-dropdown", "GencodeGeneRegion")
    select_dropdown_item(index_page, "plot-color-column-dropdown", "MotifSize")
    index_page.fill("#plot-top-n-input", "5")
    index_page.click("#search-results-plot-dialog-ok-button")
    wait_for_plots_to_finish_loading(index_page, 1)
    assert index_page.is_visible(".search-results-plot canvas")
    assert "plots=t%3Abar%2Cx%3AGencodeGeneRegion%2Cc%3AMotifSize%2Cn%3A5" in index_page.url
    assert re.match(r"Based on [\d,]+ loci in the search results", index_page.inner_text(".search-results-plot-caption"))

    # Reloading the page shows the same plot
    url = index_page.url
    index_page.goto("about:blank")
    index_page.goto(url)
    wait_for_plots_to_finish_loading(index_page, 1)
    assert index_page.inner_text(".search-results-plot-title") == "Number of loci by Gencode Gene Region, colored by Motif Size"

    # The edit dialog starts out with the plot's settings
    index_page.click(".edit-search-results-plot-link")
    index_page.wait_for_selector("#search-results-plot-dialog", state="visible")
    assert index_page.eval_on_selector("#plot-x-column-dropdown", "e => $(e).dropdown('get value')") == "GencodeGeneRegion"
    assert index_page.eval_on_selector("#plot-color-column-dropdown", "e => $(e).dropdown('get value')") == "MotifSize"
    assert index_page.input_value("#plot-top-n-input") == "5"
    index_page.check("#plot-percent-checkbox", force=True)
    index_page.click("#search-results-plot-dialog-ok-button")
    wait_for_plots_to_finish_loading(index_page, 1)
    assert "p%3A1" in index_page.url

    # Going to the next results page doesn't rerun the plot query
    num_cached_queries = index_page.evaluate("searchResultsPlotQueryResultsCache.size")
    index_page.locator(".bottom-pagination .next-page").first.click()
    index_page.wait_for_function("GLOBAL_PAGE_STATE.currentPage === 2")
    index_page.wait_for_selector("#results-table tbody tr", timeout=30000)
    assert index_page.evaluate("searchResultsPlotQueryResultsCache.size") == num_cached_queries

    index_page.click(".remove-search-results-plot-link")
    assert index_page.locator(".search-results-plot").count() == 0
    assert "plots=" not in index_page.url
    assert index_page.is_visible("#add-search-results-plot-link")

    assert real_errors(console_errors) == []
