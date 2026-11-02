import gzip
import json
import os
import tempfile
import unittest

from add_centromere_and_satellite_region_column_to_catalog_table import read_locus_ids_and_values


class ReadLocusIdsAndValuesTests(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.catalog_path = os.path.join(self.temp_dir.name, "catalog.json.gz")

    def tearDown(self):
        self.temp_dir.cleanup()

    def write_catalog(self, records):
        with gzip.open(self.catalog_path, "wt") as f:
            json.dump(records, f)

    def test_returns_the_loci_that_have_a_value(self):
        self.write_catalog([
            {"LocusId": "1-100-130-CAG", "CentromereAndSatelliteRegion": "HSat3", "Source": "Some source"},
            {"LocusId": "1-200-230-CAG", "CentromereAndSatelliteRegion": None},
            {"LocusId": "1-300-330-CAG", "CentromereAndSatelliteRegion": "GAP,HSat3"},
        ])
        self.assertEqual(read_locus_ids_and_values(self.catalog_path),
                         [("1-100-130-CAG", "HSat3"), ("1-300-330-CAG", "GAP,HSat3")])

    def test_exits_when_records_lack_the_field(self):
        self.write_catalog([
            {"LocusId": "1-100-130-CAG", "CentromereAndSatelliteRegion": "HSat3"},
            {"LocusId": "1-200-230-CAG"},
        ])
        with self.assertRaises(SystemExit):
            read_locus_ids_and_values(self.catalog_path)

    def test_exits_when_a_locus_id_appears_twice(self):
        self.write_catalog([
            {"LocusId": "1-100-130-CAG", "CentromereAndSatelliteRegion": "HSat3"},
            {"LocusId": "1-100-130-CAG", "CentromereAndSatelliteRegion": "ct"},
        ])
        with self.assertRaises(SystemExit):
            read_locus_ids_and_values(self.catalog_path)


if __name__ == "__main__":
    unittest.main()
