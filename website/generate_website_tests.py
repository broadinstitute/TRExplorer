"""Tests for generate_website.py. Run from the website/ directory with: python3 -m unittest generate_website_tests"""

import pathlib
import re
import unittest

from generate_website import find_latest_entry_date_in_whats_new

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent

# The date the generated pages show in the header's "What's New (date)" link
HEADER_WHATS_NEW_DATE_REGEX = re.compile(r'id="whats-new-link"[^>]*>What\'s New<span[^>]*>\(([^)]*)\)</span>')


class FindLatestEntryDateInWhatsNewTests(unittest.TestCase):

    def test_returns_the_newest_date_whatever_the_entry_order(self):
        whats_new_html = """
            <tr><td style="vertical-align: top;"><b>9/30/2026</b></td><td><ul><li>A</li></ul></td></tr>
            <tr><td style="vertical-align: top;"><b>10/5/2026</b></td><td><ul><li>B</li></ul></td></tr>
            <tr><td style="vertical-align: top;"><b>7/10/2025</b></td><td><ul><li>C</li></ul></td></tr>
        """
        self.assertEqual(find_latest_entry_date_in_whats_new(whats_new_html), "10/5/2026")

    def test_ignores_dates_in_the_entry_text(self):
        whats_new_html = '<tr><td><b>10/5/2026</b></td><td><ul><li>Data loaded on <b>12/1/2026</b></li></ul></td></tr>'
        self.assertEqual(find_latest_entry_date_in_whats_new(whats_new_html), "10/5/2026")

    def test_raises_when_there_are_no_entry_dates(self):
        with self.assertRaises(ValueError):
            find_latest_entry_date_in_whats_new("<tr><td>No dates here</td></tr>")


class GeneratedPagesTests(unittest.TestCase):

    def test_header_date_matches_the_newest_whats_new_entry(self):
        """Fails when whats_new.html was edited without rerunning generate_website.py."""
        latest_whats_new_date = find_latest_entry_date_in_whats_new((REPO_ROOT / "whats_new.html").read_text())
        for page_name in ["index.html", "locus.html"]:
            with self.subTest(page_name=page_name):
                self.assertEqual(
                    HEADER_WHATS_NEW_DATE_REGEX.findall((REPO_ROOT / page_name).read_text()),
                    [latest_whats_new_date],
                    f"The header date in {page_name} doesn't match whats_new.html. Rerun generate_website.py.")


if __name__ == "__main__":
    unittest.main()
