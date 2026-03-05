"""Parses JINJA templates and converts them to static HTML pages.

This script processes all *_page_template.html files in the current directory,
renders them using Jinja2, and outputs the final HTML files to the parent directory.
"""

import glob
import json
import jinja2
import os
import sys

# Add bigquery-proxy directory to path to import global_constants
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'bigquery-proxy'))
from global_constants import BIGQUERY_COLUMNS, GROUP_ORDER, get_column_descriptions, get_custom_filter_columns, get_exportable_columns

jinja2_env = jinja2.Environment(loader=jinja2.FileSystemLoader('.'))

# Get column descriptions to pass to templates
column_descriptions = get_column_descriptions()

# Get custom filter columns to pass to templates
custom_filter_columns = get_custom_filter_columns()
custom_filter_columns_json = json.dumps(custom_filter_columns)

# Get exportable columns for the export dialog
exportable_columns = get_exportable_columns()
exportable_columns_json = json.dumps(exportable_columns)

# Read the data last updated date for use in templates
data_last_updated_path = "data_last_updated_date.json"
with open(data_last_updated_path, "r") as f:
    data_last_updated = json.load(f)
data_last_updated_date = data_last_updated.get("data_last_updated_date", "")

for template_file in glob.glob("*_page_template.html"):
    print(f"Processing {template_file}")
    template_name = os.path.basename(template_file)

    template = jinja2_env.get_template(template_name)
    html_content = template.render(
        column_descriptions=column_descriptions,
        column_groups=GROUP_ORDER,
        custom_filter_columns_json=custom_filter_columns_json,
        exportable_columns_json=exportable_columns_json,
        data_last_updated_date=data_last_updated_date,
    )

    output_html_path = os.path.join("../", template_file.replace("_page_template", ""))
    with open(output_html_path, "wt") as f:
        f.write(html_content)

    print(f"Wrote {output_html_path}")
