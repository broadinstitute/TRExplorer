import datetime
import unittest

from delete_previous_bigquery_tables import decide_which_table_versions_to_keep, find_table_ids_in_text, \
    get_table_id_prefix

NOW = datetime.datetime(2026, 11, 10, 12, 0, tzinfo=datetime.timezone.utc)
TWO_DAYS_AGO = NOW - datetime.timedelta(days=2)
ONE_HOUR_AGO = NOW - datetime.timedelta(hours=1)


def decide(table_rows_and_creation_times, website_table_ids, table_ids_the_website_has_used=None):
    """Returns {table ID: True to keep, False to delete} for tables given as {table ID: (rows, created)}."""
    decisions = decide_which_table_versions_to_keep(
        list(table_rows_and_creation_times), website_table_ids,
        {table_id: rows for table_id, (rows, _) in table_rows_and_creation_times.items()},
        {table_id: created for table_id, (_, created) in table_rows_and_creation_times.items()},
        NOW, table_ids_the_website_has_used)
    return {table_id: keep for table_id, keep, _ in decisions}


class DecideWhichTableVersionsToKeepTests(unittest.TestCase):

    def test_an_empty_newest_table_left_by_a_failed_load_is_deleted_and_the_live_one_kept(self):
        self.assertEqual(decide({
            "catalog_20260918_141947": (5_000_000, TWO_DAYS_AGO),
            "catalog_20261005_174918": (5_000_000, TWO_DAYS_AGO),
            "catalog_20261030_170523": (5_660_978, TWO_DAYS_AGO),
            "catalog_20261109_090000": (0, TWO_DAYS_AGO),
        }, {"catalog_20261030_170523"},
            {"catalog_20260918_141947", "catalog_20261005_174918", "catalog_20261030_170523"}), {
            "catalog_20260918_141947": False,
            "catalog_20261005_174918": True,   # fallback
            "catalog_20261030_170523": True,   # the website uses it
            "catalog_20261109_090000": False,  # empty
        })

    def test_a_partial_test_load_never_becomes_the_fallback(self):
        # The 2026-09-18 sequence: a full load, then a 1-row load_bigquery_main_table.py -n 1 test load, which the
        # loader doesn't point the website at, then the next full load
        self.assertEqual(decide({
            "catalog_20260902_151726": (5_599_658, TWO_DAYS_AGO),
            "catalog_20260918_133634": (1, TWO_DAYS_AGO),
            "catalog_20260918_141947": (5_657_793, TWO_DAYS_AGO),
        }, {"catalog_20260918_141947"}, {"catalog_20260902_151726", "catalog_20260918_141947"}), {
            "catalog_20260902_151726": True,   # fallback
            "catalog_20260918_133634": False,
            "catalog_20260918_141947": True,
        })

    def test_an_empty_table_created_less_than_a_day_ago_is_kept_since_its_load_may_still_be_running(self):
        self.assertEqual(decide({
            "catalog_20261030_170523": (5_660_978, TWO_DAYS_AGO),
            "catalog_20261110_110000": (0, ONE_HOUR_AGO),
        }, {"catalog_20261030_170523"}), {"catalog_20261030_170523": True, "catalog_20261110_110000": True})

    def test_a_finished_load_the_website_has_not_switched_to_is_kept(self):
        self.assertEqual(decide({
            "catalog_20261030_170523": (5_660_978, TWO_DAYS_AGO),
            "catalog_20261109_090000": (5_660_978, TWO_DAYS_AGO),
        }, {"catalog_20261030_170523"}), {"catalog_20261030_170523": True, "catalog_20261109_090000": True})

    def test_the_versions_named_by_the_working_copy_and_the_published_site_are_both_kept(self):
        # The working copy already points at the newest load, while the published site still uses the previous one
        self.assertEqual(decide({
            "catalog_20261005_174918": (5_000_000, TWO_DAYS_AGO),
            "catalog_20261016_002129": (5_000_000, TWO_DAYS_AGO),
            "catalog_20261030_170523": (5_660_978, TWO_DAYS_AGO),
        }, {"catalog_20261016_002129", "catalog_20261030_170523"},
            {"catalog_20261005_174918", "catalog_20261016_002129", "catalog_20261030_170523"}), {
            "catalog_20261005_174918": True,   # the fallback is older than both
            "catalog_20261016_002129": True,
            "catalog_20261030_170523": True,
        })

    def test_without_a_fallback_every_older_version_is_deleted(self):
        self.assertEqual(decide({
            "hprc256_methylation_20260101_000000": (100, TWO_DAYS_AGO),
            "hprc256_methylation_20260521_033520": (100, TWO_DAYS_AGO),
        }, {"hprc256_methylation_20260521_033520"}), {
            "hprc256_methylation_20260101_000000": False,
            "hprc256_methylation_20260521_033520": True,
        })

    def test_nothing_is_deleted_from_a_family_the_website_names_no_version_of(self):
        self.assertEqual(decide({
            "catalog_20261005_174918": (5_000_000, TWO_DAYS_AGO),
            "catalog_20261030_170523": (0, TWO_DAYS_AGO),
        }, set()), {"catalog_20261005_174918": True, "catalog_20261030_170523": True})


class TableIdTests(unittest.TestCase):

    def test_find_table_ids_in_text(self):
        header_template_lines = (
            "const TABLE_ID = 'catalog_20261030_170523'\n"
            "const HPRC256_TOTAL_LENGTH_STRATIFIED_TABLE_ID = 'hprc256_stratified_total_allele_length_20261016_004909'\n"
            "const HPRC256_DECOMPOSED_ALLELES_TABLE_ID = "
            "'`cmg-analysis.tandem_repeat_explorer.hprc256_decomposed_alleles_20260520_210656`'\n"
            "const SOMETHING_ELSE = 'catalog_lookup'\n")
        self.assertEqual(find_table_ids_in_text(header_template_lines), {
            "catalog_20261030_170523",
            "hprc256_stratified_total_allele_length_20261016_004909",
            "hprc256_decomposed_alleles_20260520_210656",
        })

    def test_get_table_id_prefix_matches_only_timestamped_versions(self):
        self.assertEqual(get_table_id_prefix("hprc256_stratified_lps_20260902_143605"), "hprc256_stratified_lps")
        self.assertEqual(get_table_id_prefix("hprc256_stratified_total_allele_length_20261016_004909"),
                         "hprc256_stratified_total_allele_length")
        self.assertIsNone(get_table_id_prefix("catalog_20261030_170523_staging"))
        self.assertIsNone(get_table_id_prefix("centromere_and_satellite_region_staging_20261030_170523"))


if __name__ == "__main__":
    unittest.main()
