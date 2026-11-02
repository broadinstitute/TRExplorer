import json
import unittest

from bigquery_loader_utils import build_fill_sql, build_staging_rows, get_catalog_column_schema_fields, \
    get_garg_2019_and_2020_column_values, get_zhang_2025_column_values
from global_constants import GARG_2019_AND_2020_COLUMN_NAMES, ZHANG_2025_COLUMN_NAMES


class BuildStagingRowsTests(unittest.TestCase):
    def test_zhang_2025_rows_have_details_as_a_json_string(self):
        details = {"zhangLocusId": "Human_STR_7", "lengthEffectSourceLibrary": "hSTR2",
                   "hSTR2": {"lengthEffectSlope": 0.1, "lengthEffectFDR": 0.01, "pearsonR": 0.2}}
        lookup = {
            "1-100-130-CAG": {"LengthEffectTestResult": "significant", "LengthEffectSlope": 0.1,
                              "LengthEffectFDR": 0.01, "Details": details},
            "1-200-230-A": {"LengthEffectTestResult": "not enough data", "LengthEffectSlope": None,
                            "LengthEffectFDR": None, "Details": {"zhangLocusId": "Human_STR_8"}},
        }
        rows = build_staging_rows(lookup, ZHANG_2025_COLUMN_NAMES, get_zhang_2025_column_values)
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

    def test_garg_2019_and_2020_rows_have_details_as_json_strings_and_nulls_for_unmatched_datasets(self):
        autism_dmrs = [{"dmrHg38": "chr1:100-200", "methylationChange": "hypomethylated", "numProbandsWithDMR": 1}]
        schizophrenia_dmrs = [
            {"dmrHg38": "chr1:150-250", "methylationChange": "hypermethylated", "numIndividualsWithDMR": 2},
            {"dmrHg38": "chr1:240-300", "methylationChange": "hypomethylated", "numIndividualsWithDMR": 1},
        ]
        garg_2020_dmrs = [{"trHg38": "chr2:100-130", "methylationChange": "hypermethylated", "numIndividualsWithDMR": 71}]
        lookup = {
            "1-120-260-CGG": {"Garg2019_AutismDMRs": autism_dmrs, "Garg2019_SchizophreniaDMRs": schizophrenia_dmrs},
            "2-100-130-GCC": {"Garg2020_DMRs": garg_2020_dmrs},
        }
        rows = build_staging_rows(lookup, GARG_2019_AND_2020_COLUMN_NAMES, get_garg_2019_and_2020_column_values)
        self.assertEqual(rows, [
            {
                "LocusId": "1-120-260-CGG",
                "Garg2019_AutismDMRMethylationChange": "hypomethylated",
                "Garg2019_AutismDMRDetails": json.dumps(autism_dmrs),
                "Garg2019_SchizophreniaDMRMethylationChange": "hypermethylated,hypomethylated",
                "Garg2019_SchizophreniaDMRDetails": json.dumps(schizophrenia_dmrs),
                "Garg2020_NumIndividualsWithDMR": None,
                "Garg2020_DMRDetails": None,
            },
            {
                "LocusId": "2-100-130-GCC",
                "Garg2019_AutismDMRMethylationChange": None,
                "Garg2019_AutismDMRDetails": None,
                "Garg2019_SchizophreniaDMRMethylationChange": None,
                "Garg2019_SchizophreniaDMRDetails": None,
                "Garg2020_NumIndividualsWithDMR": 71,
                "Garg2020_DMRDetails": json.dumps(garg_2020_dmrs),
            },
        ])

    def test_values_that_dont_match_the_columns_exit(self):
        with self.assertRaises(SystemExit):
            build_staging_rows({"1-100-130-CAG": {}}, ["ColumnA", "ColumnB"], lambda entry: {"ColumnA": 1})


class BuildFillSqlTests(unittest.TestCase):
    def test_copies_matching_rows_and_clears_the_other_rows_that_have_a_value(self):
        self.assertEqual(
            build_fill_sql("project.dataset.catalog_1", "project.dataset.staging_1", ZHANG_2025_COLUMN_NAMES),
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

    def test_one_column(self):
        self.assertEqual(
            build_fill_sql("project.dataset.catalog_1", "project.dataset.staging_1", ["CentromereAndSatelliteRegion"]),
            "MERGE `project.dataset.catalog_1` T USING `project.dataset.staging_1` S ON T.LocusId = S.LocusId "
            "WHEN MATCHED THEN UPDATE SET T.CentromereAndSatelliteRegion = S.CentromereAndSatelliteRegion "
            "WHEN NOT MATCHED BY SOURCE AND (T.CentromereAndSatelliteRegion IS NOT NULL) "
            "THEN UPDATE SET T.CentromereAndSatelliteRegion = NULL")


class GetCatalogColumnSchemaFieldsTests(unittest.TestCase):
    def test_uses_the_types_declared_in_global_constants(self):
        self.assertEqual([(field.name, field.field_type) for field in get_catalog_column_schema_fields(
            ZHANG_2025_COLUMN_NAMES + GARG_2019_AND_2020_COLUMN_NAMES)], [
            ("Zhang2025_LengthEffectTestResult", "STRING"),
            ("Zhang2025_LengthEffectSlope", "FLOAT"),
            ("Zhang2025_LengthEffectFDR", "FLOAT"),
            ("Zhang2025_Details", "STRING"),
            ("Garg2019_AutismDMRMethylationChange", "STRING"),
            ("Garg2019_AutismDMRDetails", "STRING"),
            ("Garg2019_SchizophreniaDMRMethylationChange", "STRING"),
            ("Garg2019_SchizophreniaDMRDetails", "STRING"),
            ("Garg2020_NumIndividualsWithDMR", "INTEGER"),
            ("Garg2020_DMRDetails", "STRING"),
        ])


if __name__ == "__main__":
    unittest.main()
