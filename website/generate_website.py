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
from global_constants import BIGQUERY_COLUMNS, get_column_descriptions, get_exportable_columns

jinja2_env = jinja2.Environment(loader=jinja2.FileSystemLoader('.'))

# Get column descriptions to pass to templates
column_descriptions = get_column_descriptions()

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
        exportable_columns_json=exportable_columns_json,
        data_last_updated_date=data_last_updated_date,
    )

    output_html_path = os.path.join("../", template_file.replace("_page_template", ""))
    with open(output_html_path, "wt") as f:
        f.write(html_content)

    print(f"Wrote {output_html_path}")

# Update whats_new.html with the latest database update info
whats_new_path = "../whats_new.html"

update_date = data_last_updated_date
update_details = data_last_updated.get("data_last_updated_details", "")

if update_date:
    # Build the description text
    description = "Database updated."
    if update_details:
        description += f" {update_details}"

    # Create the new table row
    new_row = f'''        <tr class="ui table row" style="border: 0px !important;">
          <td style="vertical-align: top; font-size: 1.1em;">{update_date}</td>
          <td style="vertical-align: top; font-size: 1.1em;">
              <ul>
                  <li>{description}</li>
              </ul>
          </td>
        </tr>
'''

    with open(whats_new_path, "r") as f:
        whats_new_content = f.read()

    # Remove any existing "Database updated." entry with the same date (to avoid duplicates on re-runs)
    escaped_date = re.escape(update_date)
    # Match new format (with ul/li)
    whats_new_content = re.sub(
        rf'\s*<tr class="ui table row"[^>]*>\s*<td[^>]*>{escaped_date}</td>\s*<td[^>]*>\s*<ul>\s*<li>Database updated\.[^<]*</li>\s*</ul>\s*</td>\s*</tr>',
        '',
        whats_new_content
    )
    # Match old format (plain text)
    whats_new_content = re.sub(
        rf'\s*<tr class="ui table row"[^>]*>\s*<td[^>]*>{escaped_date}</td>\s*<td[^>]*>Database updated\.[^<]*</td>\s*</tr>',
        '',
        whats_new_content
    )

    # Insert the new row after <tbody class="ui table body">
    whats_new_content = re.sub(
        r'(<tbody class="ui table body">)',
        r'\1\n' + new_row,
        whats_new_content
    )

    with open(whats_new_path, "w") as f:
        f.write(whats_new_content)

    print(f"Updated {whats_new_path} with database update entry: {update_date}")
