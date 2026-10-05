import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from generate_Zhang_2025_lookup_json import (
    TEST_RESULT_NOT_ENOUGH_DATA,
    TEST_RESULT_NOT_SIGNIFICANT,
    TEST_RESULT_SIGNIFICANT,
    classify_test_result,
    collect_catalog_candidates,
    find_trexplorer_match,
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

    def test_pick_best_entry_prefers_analyzed_then_lowest_fdr(self):
        entries = [
            {"LengthEffectFDR": None, "id": "not analyzed"},
            {"LengthEffectFDR": 0.5, "id": "weak"},
            {"LengthEffectFDR": 0.01, "id": "strong"},
        ]
        self.assertEqual(pick_best_entry(entries)["id"], "strong")
        self.assertEqual(pick_best_entry(entries[:1])["id"], "not analyzed")

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
