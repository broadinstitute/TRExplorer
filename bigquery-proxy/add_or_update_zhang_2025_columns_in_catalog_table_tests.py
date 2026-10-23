import json
import unittest

from add_or_update_zhang_2025_columns_in_catalog_table import (
    build_fill_sql, build_staging_rows, zhang_2025_schema_fields)
from global_constants import ZHANG_2025_COLUMN_NAMES


class BuildStagingRowsTests(unittest.TestCase):
    def test_returns_one_row_per_locus_with_details_as_a_json_string(self):
        details = {"zhangLocusId": "Human_STR_7", "lengthEffectSourceLibrary": "hSTR2",
                   "hSTR2": {"lengthEffectSlope": 0.1, "lengthEffectFDR": 0.01, "pearsonR": 0.2}}
        lookup = {
            "1-100-130-CAG": {"LengthEffectTestResult": "significant", "LengthEffectSlope": 0.1,
                              "LengthEffectFDR": 0.01, "Details": details},
            "1-200-230-A": {"LengthEffectTestResult": "not enough data", "LengthEffectSlope": None,
                            "LengthEffectFDR": None, "Details": {"zhangLocusId": "Human_STR_8"}},
        }
        rows = build_staging_rows(lookup)
        self.assertEqual(rows[0], {
            "LocusId": "1-100-130-CAG",
            "Zhang2025_LengthEffectTestResult": "significant",
            "Zhang2025_LengthEffectSlope": 0.1,
            "Zhang2025_LengthEffectFDR": 0.01,
            "Zhang2025_Details": json.dumps(details),
        })
        self.assertEqual(rows[1]["LocusId"], "1-200-230-A")
        self.assertEqual((rows[1]["Zhang2025_LengthEffectSlope"], rows[1]["Zhang2025_LengthEffectFDR"]), (None, None))
        self.assertEqual(set(rows[1]), {"LocusId"} | set(ZHANG_2025_COLUMN_NAMES))


class BuildFillSqlTests(unittest.TestCase):
    def test_copies_matching_rows_and_clears_the_other_rows_that_have_a_value(self):
        self.assertEqual(
            build_fill_sql("project.dataset.catalog_1", "project.dataset.staging_1"),
            "MERGE `project.dataset.catalog_1` T USING `project.dataset.staging_1` S ON T.LocusId = S.LocusId "
            "WHEN MATCHED THEN UPDATE SET "
            "T.Zhang2025_LengthEffectTestResult = S.Zhang2025_LengthEffectTestResult, "
            "T.Zhang2025_LengthEffectSlope = S.Zhang2025_LengthEffectSlope, "
            "T.Zhang2025_LengthEffectFDR = S.Zhang2025_LengthEffectFDR, "
            "T.Zhang2025_Details = S.Zhang2025_Details "
            "WHEN NOT MATCHED BY SOURCE AND (T.Zhang2025_LengthEffectTestResult IS NOT NULL OR "
            "T.Zhang2025_LengthEffectSlope IS NOT NULL OR T.Zhang2025_LengthEffectFDR IS NOT NULL OR "
            "T.Zhang2025_Details IS NOT NULL) "
            "THEN UPDATE SET T.Zhang2025_LengthEffectTestResult = NULL, T.Zhang2025_LengthEffectSlope = NULL, "
            "T.Zhang2025_LengthEffectFDR = NULL, T.Zhang2025_Details = NULL")


class Zhang2025SchemaFieldsTests(unittest.TestCase):
    def test_uses_the_types_declared_in_global_constants(self):
        self.assertEqual([(field.name, field.field_type) for field in zhang_2025_schema_fields()], [
            ("Zhang2025_LengthEffectTestResult", "STRING"),
            ("Zhang2025_LengthEffectSlope", "FLOAT"),
            ("Zhang2025_LengthEffectFDR", "FLOAT"),
            ("Zhang2025_Details", "STRING"),
        ])


if __name__ == "__main__":
    unittest.main()
