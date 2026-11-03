import unittest
from unittest import mock

import flask

import main

LOOKUP_FILE_TEXT = (
    "name_upper\tchrom\tstart_0based\tend_1based\n"
    "5S_RRNA\tchr2\t89599771\t89600561\n"
    "5S_RRNA\tchr7\t152593090\t152593109\n"
    "ATXN1\tchr6\t16300565\t16761662\n"
    "ATXN10\tchr22\t45671000\t45845000\n"
    "ATXN1-AS1\tchr6\t16700000\t16720000\n"
    "\n"
)

VALID_REQUEST_JSON = {"names": ["ATXN1"], "genome_version": "hg38", "catalog_table_id": "catalog_20261030_170523"}


class LookupFileTests(unittest.TestCase):

    def test_parse_lookup_file_text_sorts_the_lines_and_drops_the_header_and_empty_lines(self):
        lines = main.parse_lookup_file_text(LOOKUP_FILE_TEXT)
        self.assertEqual(len(lines), 5)
        self.assertEqual(lines, sorted(lines))

    def test_parse_lookup_file_text_rejects_an_unexpected_header(self):
        with self.assertRaises(ValueError):
            main.parse_lookup_file_text("name\tchrom\tstart\tend\nATXN1\tchr6\t1\t2\n")

    def test_find_intervals_of_name(self):
        lines = main.parse_lookup_file_text(LOOKUP_FILE_TEXT)
        self.assertEqual(main.find_intervals_of_name("ATXN1", lines),
                         [{"chrom": "chr6", "start_0based": 16300565, "end_1based": 16761662}])
        self.assertEqual(main.find_intervals_of_name("5S_RRNA", lines), [
            {"chrom": "chr2", "start_0based": 89599771, "end_1based": 89600561},
            {"chrom": "chr7", "start_0based": 152593090, "end_1based": 152593109},
        ])
        self.assertEqual(main.find_intervals_of_name("ATXN1-AS1", lines),
                         [{"chrom": "chr6", "start_0based": 16700000, "end_1based": 16720000}])
        for name_without_loci in ["ATXN", "ATXN2", "ATXN1-AS", "ZZZ", "", "ATXN1\tchr6"]:
            self.assertEqual(main.find_intervals_of_name(name_without_loci, lines), [], name_without_loci)

    def test_look_up_intervals_of_names_ignores_case(self):
        lines = main.parse_lookup_file_text(LOOKUP_FILE_TEXT)
        self.assertEqual(main.look_up_intervals_of_names(["atxn1", "FOOBAR"], lines), {
            "atxn1": [{"chrom": "chr6", "start_0based": 16300565, "end_1based": 16761662}],
            "FOOBAR": [],
        })

    def test_get_sorted_lookup_file_lines_keeps_the_two_most_recently_used_files(self):
        main.get_sorted_lookup_file_lines.cache_clear()
        with mock.patch.object(main, "download_and_parse_lookup_file", side_effect=lambda url: [url]) as download:
            for catalog_table_id in ["catalog_20260101_000000", "catalog_20260101_000000", "catalog_20260202_000000",
                                     "catalog_20260303_000000", "catalog_20260202_000000"]:
                lines = main.get_sorted_lookup_file_lines("hg38", catalog_table_id)
                self.assertTrue(lines[0].endswith(f"/hg38/{catalog_table_id}.locus_intervals_by_name.tsv.gz"))
            self.assertEqual(download.call_count, 3)
            self.assertEqual(main.get_sorted_lookup_file_lines.cache_info().currsize, 2)
        main.get_sorted_lookup_file_lines.cache_clear()

    def test_get_sorted_lookup_file_lines_reads_a_missing_file_again_on_the_next_request(self):
        main.get_sorted_lookup_file_lines.cache_clear()
        with mock.patch.object(main, "download_and_parse_lookup_file",
                               side_effect=main.LookupFileNotFoundError("There's no lookup file at ...")) as download:
            for _ in range(2):
                with self.assertRaises(main.LookupFileNotFoundError):
                    main.get_sorted_lookup_file_lines("hg38", "catalog_20260101_000000")
            self.assertEqual(download.call_count, 2)
        main.get_sorted_lookup_file_lines.cache_clear()


class RequestTests(unittest.TestCase):

    def test_parse_request_json_accepts_a_valid_request(self):
        self.assertEqual(main.parse_request_json(VALID_REQUEST_JSON), (["ATXN1"], "hg38", "catalog_20261030_170523"))

    def test_parse_request_json_rejects_invalid_requests(self):
        for request_json in [
            None,
            ["ATXN1"],
            {**VALID_REQUEST_JSON, "names": []},
            {**VALID_REQUEST_JSON, "names": "ATXN1"},
            {**VALID_REQUEST_JSON, "names": ["ATXN1", 5]},
            {**VALID_REQUEST_JSON, "names": ["A"] * (main.MAX_NAMES_PER_REQUEST + 1)},
            {**VALID_REQUEST_JSON, "names": ["A" * (main.MAX_NAME_LENGTH + 1)]},
            {**VALID_REQUEST_JSON, "genome_version": "chm13"},
            {**VALID_REQUEST_JSON, "genome_version": ["hg38"]},
            {key: value for key, value in VALID_REQUEST_JSON.items() if key != "genome_version"},
            {**VALID_REQUEST_JSON, "catalog_table_id": "catalog"},
            {**VALID_REQUEST_JSON, "catalog_table_id": "../catalog_20261030_170523"},
            {key: value for key, value in VALID_REQUEST_JSON.items() if key != "catalog_table_id"},
        ]:
            with self.assertRaises(ValueError, msg=str(request_json)[:200]):
                main.parse_request_json(request_json)


class EntryPointTests(unittest.TestCase):

    def setUp(self):
        self.app = flask.Flask(__name__)
        main.get_sorted_lookup_file_lines.cache_clear()

    def tearDown(self):
        main.get_sorted_lookup_file_lines.cache_clear()

    def call_entry_point(self, method, request_json=None, data=None):
        with self.app.test_request_context(method=method, json=request_json, data=data):
            body, status_code, headers = main.lookup_locus_intervals_by_name(flask.request)
            return (body.get_json() if hasattr(body, "get_json") else body), status_code, headers

    def test_preflight_and_other_methods(self):
        body, status_code, headers = self.call_entry_point("OPTIONS")
        self.assertEqual((body, status_code), ("", 204))
        self.assertEqual(headers["Access-Control-Allow-Origin"], "*")
        self.assertEqual(self.call_entry_point("GET")[1], 405)

    def test_invalid_request_gives_400(self):
        body, status_code, _ = self.call_entry_point("POST", data="not json")
        self.assertEqual(status_code, 400)
        self.assertIn("JSON object", body["error"])
        self.assertEqual(self.call_entry_point("POST", {**VALID_REQUEST_JSON, "genome_version": "chm13"})[1], 400)

    def test_missing_lookup_file_gives_404(self):
        with mock.patch.object(main, "download_and_parse_lookup_file",
                               side_effect=main.LookupFileNotFoundError("There's no lookup file at ...")):
            body, status_code, _ = self.call_entry_point("POST", VALID_REQUEST_JSON)
        self.assertEqual(status_code, 404)
        self.assertIn("catalog_20261030_170523", body["error"])

    def test_unreadable_lookup_file_gives_500(self):
        with mock.patch.object(main, "download_and_parse_lookup_file", side_effect=OSError("connection reset")):
            body, status_code, _ = self.call_entry_point("POST", VALID_REQUEST_JSON)
        self.assertEqual(status_code, 500)
        self.assertIn("connection reset", body["error"])

    def test_lookup(self):
        with mock.patch.object(main, "download_and_parse_lookup_file",
                               return_value=main.parse_lookup_file_text(LOOKUP_FILE_TEXT)):
            body, status_code, _ = self.call_entry_point("POST", {**VALID_REQUEST_JSON, "names": ["atxn1", "FOOBAR"]})
        self.assertEqual(status_code, 200)
        self.assertEqual(body, {"intervals_by_name": {
            "atxn1": [{"chrom": "chr6", "start_0based": 16300565, "end_1based": 16761662}],
            "FOOBAR": [],
        }})


if __name__ == "__main__":
    unittest.main()
