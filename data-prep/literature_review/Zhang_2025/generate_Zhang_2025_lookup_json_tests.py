import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from generate_Zhang_2025_lookup_json import (
    TEST_RESULT_NOT_ENOUGH_DATA,
    TEST_RESULT_NOT_SIGNIFICANT,
    TEST_RESULT_SIGNIFICANT,
    build_lookup_entry,
    classify_test_result,
    collect_catalog_candidates,
    find_trexplorer_match,
    merge_regression_results,
    pick_best_entry,
    read_zhang_strs,
)


class Zhang2025LookupTests(unittest.TestCase):

    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.dir = Path(self.temp_dir.name)

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_classify_test_result(self):
        self.assertEqual(classify_test_result(0.05), TEST_RESULT_SIGNIFICANT)
        self.assertEqual(classify_test_result(0.1), TEST_RESULT_NOT_SIGNIFICANT)
        self.assertEqual(classify_test_result(float("nan")), TEST_RESULT_NOT_ENOUGH_DATA)
        self.assertEqual(classify_test_result(None), TEST_RESULT_NOT_ENOUGH_DATA)

    def test_read_zhang_strs_uses_one_row_per_human_str_and_converts_to_0based(self):
        pairs = self.dir / "pairs.tsv"
        pd.DataFrame([
            ["hg38", "chr1", 101, 110, "Human_STR_7"],
            ["hg38", "chr1", 101, 110, "Human_STR_7"],   # same STR near a second TSS
            ["mm10", "chr1", 5, 9, "MOUSE_STR_7"],
        ], columns=["organism", "str_chr", "str_pos", "str_end", "str_id"]).to_csv(pairs, sep="\t", index=False)
        probes = self.dir / "probes.tsv"
        pd.DataFrame([
            ["hg38", "Human_STR_7", "CAG", "ref"],
            ["hg38", "Human_STR_7", "CAG", "p5"],
            ["mm10", "MOUSE_STR_7", "A", "ref"],
        ], columns=["organism", "id", "motif", "allele"]).to_csv(probes, sep="\t", index=False)

        strs = read_zhang_strs(pairs, probes)
        self.assertEqual(strs.values.tolist(), [[7, "Human_STR_7", "1", 100, 110, "CAG"]])

    def test_find_trexplorer_match(self):
        candidates = {"1": [
            {"locus_id": "1-100-130-CAG", "start_0based": 100, "end_1based": 130, "canonical_motif": "AGC", "motif_length": 3},
            {"locus_id": "1-100-130-A", "start_0based": 100, "end_1based": 130, "canonical_motif": "A", "motif_length": 1},
        ]}
        locus_ids = {c["locus_id"] for c in candidates["1"]}
        self.assertEqual(find_trexplorer_match("1", 100, 130, "CAG", locus_ids, candidates), ("1-100-130-CAG", "exact"))
        # Rotated motif and slightly different interval: fuzzy match (Jaccard 27/30 = 0.9).
        self.assertEqual(find_trexplorer_match("1", 103, 130, "GCA", locus_ids, candidates), ("1-100-130-CAG", "fuzzy"))
        # Jaccard 10/30 = 0.33 is above the 0.2 cutoff; 5/30 = 0.17 is below it.
        self.assertEqual(find_trexplorer_match("1", 120, 130, "CAG", locus_ids, candidates), ("1-100-130-CAG", "fuzzy"))
        self.assertEqual(find_trexplorer_match("1", 125, 130, "CAG", locus_ids, candidates), (None, None))
        # Different motif.
        self.assertEqual(find_trexplorer_match("1", 100, 130, "GT", locus_ids, candidates), (None, None))

    def test_merge_regression_results(self):
        strs = pd.DataFrame({"zhang_locus_number": [7, 8], "motif": ["CAG", "A"]})
        results = pd.DataFrame({"locus": [8], "repeat_unit": ["A"], "r": [0.4], "beta_0": [1.0], "beta_1": [0.2],
                                "padj": [0.05]})
        merged = merge_regression_results(strs, results, "hSTR2")
        self.assertEqual(list(merged.columns), ["zhang_locus_number", "motif", "hSTR2_pearson_r", "hSTR2_slope", "hSTR2_fdr"])
        self.assertTrue(merged.iloc[0][["hSTR2_pearson_r", "hSTR2_slope", "hSTR2_fdr"]].isna().all())
        self.assertEqual(merged.iloc[1][["hSTR2_pearson_r", "hSTR2_slope", "hSTR2_fdr"]].tolist(), [0.4, 0.2, 0.05])

        # A result for an STR that's not in the library.
        with self.assertRaises(ValueError):
            merge_regression_results(strs, results.assign(locus=[9]), "hSTR2")
        # A tab-separated file read as comma-separated has a single column.
        with self.assertRaises(ValueError):
            merge_regression_results(strs, pd.DataFrame({"locus\trepeat_unit\tr\tbeta_1\tpadj": ["8\tA\t0.4\t0.2\t0.05"]}), "hSTR2")

    def test_build_lookup_entry_uses_hstr1_then_hstr2(self):
        hstr1 = {"lengthEffectSlope": 0.3, "lengthEffectFDR": 0.5, "pearsonR": 0.1}
        hstr2 = {"lengthEffectSlope": 0.1, "lengthEffectFDR": 0.01, "pearsonR": 0.2}
        details = {"zhangLocusId": "Human_STR_7"}

        entry = build_lookup_entry({"hSTR1": hstr1, "hSTR2": hstr2}, details)
        self.assertEqual(entry["LengthEffectTestResult"], TEST_RESULT_NOT_SIGNIFICANT)
        self.assertEqual((entry["LengthEffectSlope"], entry["LengthEffectFDR"]), (0.3, 0.5))
        self.assertEqual(entry["Details"], {"zhangLocusId": "Human_STR_7", "lengthEffectSourceLibrary": "hSTR1",
                                            "hSTR1": hstr1, "hSTR2": hstr2})

        entry = build_lookup_entry({"hSTR2": hstr2}, details)
        self.assertEqual(entry["LengthEffectTestResult"], TEST_RESULT_SIGNIFICANT)
        self.assertEqual((entry["LengthEffectSlope"], entry["LengthEffectFDR"]), (0.1, 0.01))
        self.assertEqual(entry["Details"]["lengthEffectSourceLibrary"], "hSTR2")
        self.assertNotIn("hSTR1", entry["Details"])

        entry = build_lookup_entry({}, details)
        self.assertEqual(entry["LengthEffectTestResult"], TEST_RESULT_NOT_ENOUGH_DATA)
        self.assertEqual((entry["LengthEffectSlope"], entry["LengthEffectFDR"]), (None, None))
        self.assertEqual(entry["Details"], {"zhangLocusId": "Human_STR_7", "lengthEffectSourceLibrary": None})

    def test_pick_best_entry_prefers_hstr1_then_hstr2_then_lowest_fdr(self):
        def entry(source_library, fdr, entry_id):
            return {"LengthEffectFDR": fdr, "Details": {"lengthEffectSourceLibrary": source_library}, "id": entry_id}

        not_analyzed = entry(None, None, "not analyzed")
        weak_hstr1 = entry("hSTR1", 0.5, "weak hSTR1")
        strong_hstr1 = entry("hSTR1", 0.01, "strong hSTR1")
        strong_hstr2 = entry("hSTR2", 0.001, "strong hSTR2")
        self.assertEqual(pick_best_entry([not_analyzed, weak_hstr1, strong_hstr1, strong_hstr2])["id"], "strong hSTR1")
        self.assertEqual(pick_best_entry([not_analyzed, weak_hstr1, strong_hstr2])["id"], "weak hSTR1")
        self.assertEqual(pick_best_entry([not_analyzed, strong_hstr2])["id"], "strong hSTR2")
        self.assertEqual(pick_best_entry([not_analyzed])["id"], "not analyzed")

    def test_collect_catalog_candidates_keeps_only_overlapping_loci(self):
        catalog = self.dir / "catalog.json"
        catalog.write_text(json.dumps([
            {"LocusId": "1-100-130-CAG", "ReferenceRegion": "chr1:100-130", "LocusStructure": "(CAG)*"},
            {"LocusId": "1-5000-5010-A", "ReferenceRegion": "chr1:5000-5010", "LocusStructure": "(A)*"},
            {"LocusId": "2-100-130-A", "ReferenceRegion": "chr2:100-130", "LocusStructure": "(A)*"},
        ]))
        locus_ids, candidates_by_chrom = collect_catalog_candidates(str(catalog), {"1": [(120, 140)]}, max_str_length=20)
        self.assertEqual(locus_ids, {"1-100-130-CAG"})
        self.assertEqual([c["canonical_motif"] for c in candidates_by_chrom["1"]], ["AGC"])


if __name__ == "__main__":
    unittest.main()
