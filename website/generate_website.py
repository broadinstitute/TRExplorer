"""Parses the JINJA templates and converts them to static HTML pages."""

import glob
import jinja2
import os
import re

jinja2_env = jinja2.Environment(loader=jinja2.FileSystemLoader('.'))

for template_file in glob.glob("*_page_template.html"):
    print(f"Processing {template_file}")
    template_name = os.path.basename(template_file)

    template = jinja2_env.get_template(template_name)
    html_content = template.render()

    output_html_path = os.path.join("../", template_file.replace("_page_template", ""))
    with open(output_html_path, "wt") as f:
        f.write(html_content)

    print(f"Wrote {output_html_path}")
