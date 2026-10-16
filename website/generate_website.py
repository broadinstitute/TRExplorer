"""Parses JINJA templates and converts them to static HTML pages.

This script processes all *_page_template.html files in the current directory,
renders them using Jinja2, and outputs the final HTML files to the parent directory.
"""

import glob
import json
import jinja2
import os
import re
import sys

# Add bigquery-proxy directory to path to import global_constants
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'bigquery-proxy'))
from global_constants import MAIN_BIGQUERY_TABLE_COLUMNS, GROUP_ROWS_IN_COLUMN_DIALOGS, GENE_REGION_PRIORITY, get_column_descriptions, get_custom_filter_columns, get_exportable_columns, find_problems_in_column_descriptions_shown_on_website

# Each What's New entry starts with its date in a bold table cell, e.g. <td style="..."><b>10/15/2026</b></td>
WHATS_NEW_ENTRY_DATE_CELL_REGEX = re.compile(r"<td[^>]*><b>((\d{1,2})/(\d{1,2})/(\d{4}))</b></td>")


def find_latest_entry_date_in_whats_new(whats_new_html):
    """Returns the newest entry date in the What's New page, written the way the page writes it.

    Args:
        whats_new_html (str): contents of whats_new.html.

    Returns:
        str: the newest entry date, e.g. "10/15/2026".

    Raises:
        ValueError: if the page has no entry dates.
    """
    date_cells = list(WHATS_NEW_ENTRY_DATE_CELL_REGEX.finditer(whats_new_html))
    if not date_cells:
        raise ValueError("Found no entry dates like <td><b>10/15/2026</b></td> in whats_new.html")

    latest_date_cell = max(date_cells, key=lambda m: (int(m.group(4)), int(m.group(2)), int(m.group(3))))
    return latest_date_cell.group(1)


def main():
    # The website shows the main table's column descriptions in help text popups, so stop before writing any page if
    # one of them would render wrong there
    column_description_problems = find_problems_in_column_descriptions_shown_on_website(MAIN_BIGQUERY_TABLE_COLUMNS)
    if column_description_problems:
        sys.exit("Fix these column descriptions in bigquery-proxy/global_constants.py:\n  " + "\n  ".join(column_description_problems))

    jinja2_env = jinja2.Environment(loader=jinja2.FileSystemLoader('.'))

    # Get column descriptions to pass to templates
    column_descriptions = get_column_descriptions()

    # Get custom filter columns to pass to templates
    custom_filter_columns = get_custom_filter_columns()
    custom_filter_columns_json = json.dumps(custom_filter_columns)

    # Get exportable columns for the export dialog
    exportable_columns = get_exportable_columns()
    exportable_columns_json = json.dumps(exportable_columns)

    # The header's "What's New (date)" link shows the newest What's New entry date
    with open("../whats_new.html", "r") as f:
        latest_whats_new_date = find_latest_entry_date_in_whats_new(f.read())

    for template_file in glob.glob("*_page_template.html"):
        print(f"Processing {template_file}")
        template_name = os.path.basename(template_file)

        template = jinja2_env.get_template(template_name)
        html_content = template.render(
            column_descriptions=column_descriptions,
            column_group_rows=GROUP_ROWS_IN_COLUMN_DIALOGS,
            gene_region_priority=GENE_REGION_PRIORITY,
            custom_filter_columns_json=custom_filter_columns_json,
            exportable_columns_json=exportable_columns_json,
            latest_whats_new_date=latest_whats_new_date,
        )

        output_html_path = os.path.join("../", template_file.replace("_page_template", ""))
        with open(output_html_path, "wt") as f:
            f.write(html_content)

        print(f"Wrote {output_html_path}")


if __name__ == "__main__":
    main()
