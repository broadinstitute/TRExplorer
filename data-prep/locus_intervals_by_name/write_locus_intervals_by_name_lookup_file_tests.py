import gzip
import importlib.util
import re
import tempfile
import unittest
from pathlib import Path

from write_locus_intervals_by_name_lookup_file import (
    CATALOG_TABLE_ID_REGEXP,
    GENE_AND_TRANSCRIPT_NAME_COLUMNS,
    LOCUS_IDS_FOUND_BY_LOCUS_ID_SEARCH_SQL_REGEXP,
    LOOKUP_FILE_BUCKET,
    LOOKUP_FILE_COLUMNS,
    LOOKUP_FILE_PATH_TEMPLATE,
    SUPPORTED_GENOME_VERSIONS,
    build_sql_for_locus_intervals_by_name,
    check_locus_intervals_by_name_rows,
    format_locus_intervals_by_name_tsv,
    write_gzipped_text,
)

REPO_ROOT = Path(__file__).resolve().parents[2]
HEADER_TEMPLATE_PATH = REPO_ROOT / "website" / "header_template.html"
INDEX_PAGE_TEMPLATE_PATH = REPO_ROOT / "website" / "index_page_template.html"
LOOKUP_FUNCTION_MAIN_PATH = REPO_ROOT / "locus-intervals-by-name-function" / "main.py"


def load_lookup_function_module():
    """Returns the Cloud Function's main.py as a module. It's deployed from its own directory, so it can't import this
    script's constants and has its own copies of them."""
    spec = importlib.util.spec_from_file_location("locus_intervals_by_name_function_main", LOOKUP_FUNCTION_MAIN_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def read_website_locus_id_regexp():
    """Returns the website's LOCUS_ID_REGEXP as a Python regexp."""
    pattern, flags = re.search(r"const LOCUS_ID_REGEXP = /(.+)/(\w*)\n", HEADER_TEMPLATE_PATH.read_text()).groups()
    return re.compile(pattern, re.IGNORECASE if "i" in flags else 0)


class WriteLocusIntervalsByNameLookupFileTests(unittest.TestCase):

    def test_sql_reads_every_column_that_the_website_matches_names_against(self):
        website_columns = re.findall(r"'(\w+)'", re.search(
            r"const GENE_AND_TRANSCRIPT_SEARCH_COLUMNS = \[(.*?)\]", INDEX_PAGE_TEMPLATE_PATH.read_text(), re.DOTALL).group(1))
        self.assertEqual(website_columns, GENE_AND_TRANSCRIPT_NAME_COLUMNS + ["LocusId"])
        sql = build_sql_for_locus_intervals_by_name("project.dataset.catalog_20261030_170523")
        for column in GENE_AND_TRANSCRIPT_NAME_COLUMNS:
            self.assertIn(column, sql)
        self.assertIn("FROM `project.dataset.catalog_20261030_170523`", sql)
        self.assertIn(f"IF(REGEXP_CONTAINS(LocusId, r'{LOCUS_IDS_FOUND_BY_LOCUS_ID_SEARCH_SQL_REGEXP}'), NULL, LocusId)", sql)

    def test_sql_locus_id_regexp_matches_only_what_the_website_searches_for_as_a_locus_id(self):
        website_locus_id_regexp = read_website_locus_id_regexp()
        sql_regexp = re.compile(LOCUS_IDS_FOUND_BY_LOCUS_ID_SEARCH_SQL_REGEXP)
        for text in ["1-3074876-3074933-CAG", "chrX-100-200-cag", "CHR4-1,000-2,000-CAG", "MT-5-10-ACGTN",
                     " 1-100-200-RYMKSWHBVD ", "1-100-200-CAGX", "HTT", "1-100-CAG", "chr1:100-200", "1_100_200_CAG",
                     "1-100-200-"]:
            if sql_regexp.search(text):
                self.assertTrue(website_locus_id_regexp.search(text), text)
        self.assertTrue(sql_regexp.search("1-3074876-3074933-CAG"))
        self.assertTrue(sql_regexp.search("chrX-100-200-cag"))
        self.assertFalse(sql_regexp.search("HTT"))
        self.assertFalse(sql_regexp.search("1-100-200-CAGX"))

    def test_cloud_function_reads_the_files_that_this_script_writes(self):
        lookup_function = load_lookup_function_module()
        self.assertEqual(lookup_function.LOOKUP_FILE_URL_TEMPLATE,
                         f"https://storage.googleapis.com/{LOOKUP_FILE_BUCKET}/{LOOKUP_FILE_PATH_TEMPLATE}")
        self.assertEqual(lookup_function.LOOKUP_FILE_HEADER, "\t".join(LOOKUP_FILE_COLUMNS))
        self.assertEqual(lookup_function.CATALOG_TABLE_ID_REGEXP.pattern, CATALOG_TABLE_ID_REGEXP.pattern)
        self.assertEqual(lookup_function.SUPPORTED_GENOME_VERSIONS, SUPPORTED_GENOME_VERSIONS)

    def test_check_rows_accepts_valid_rows(self):
        check_locus_intervals_by_name_rows([("ATXN1", "chr6", 16300565, 16761730), ("5S_RRNA", "chr2", 10, 20),
                                            ("5S_RRNA", "chr7", 30, 30)])

    def test_check_rows_rejects_rows_that_the_lookup_file_cant_hold(self):
        for rows in [
            [],
            [("", "chr6", 1, 2)],
            [("A\tB", "chr6", 1, 2)],
            [("A\nB", "chr6", 1, 2)],
            [("ATXN1", "6", 1, 2)],
            [("ATXN1", "chr6", 3, 2)],
            [("ATXN1", "chr6", -1, 2)],
            [("ATXN1", "chr6", 1, 2), ("ATXN1", "chr6", 5, 9)],
        ]:
            with self.assertRaises(ValueError, msg=str(rows)):
                check_locus_intervals_by_name_rows(rows)

    def test_format_tsv(self):
        self.assertEqual(
            format_locus_intervals_by_name_tsv([("ATXN1", "chr6", 16300565, 16761730), ("HTT", "chr4", 3074876, 3074933)]),
            "name_upper\tchrom\tstart_0based\tend_1based\n"
            "ATXN1\tchr6\t16300565\t16761730\n"
            "HTT\tchr4\t3074876\t3074933\n")

    def test_write_gzipped_text_gives_the_same_bytes_for_the_same_text(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            first_path, second_path = Path(temp_dir) / "first.tsv.gz", Path(temp_dir) / "second.tsv.gz"
            write_gzipped_text("a\tb\n", first_path)
            write_gzipped_text("a\tb\n", second_path)
            self.assertEqual(first_path.read_bytes(), second_path.read_bytes())
            self.assertEqual(gzip.decompress(first_path.read_bytes()).decode(), "a\tb\n")


if __name__ == "__main__":
    unittest.main()
