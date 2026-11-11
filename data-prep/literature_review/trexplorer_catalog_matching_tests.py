import json
import tempfile
import unittest
from pathlib import Path

from trexplorer_catalog_matching import build_catalog_interval_trees, find_exact_or_best_jaccard_match, \
    is_matching_motif, parse_catalog_record


def get_locus_ids(interval_trees):
    """Returns the set of locus IDs in interval trees returned by build_catalog_interval_trees."""
    return {interval.data["locus_id"] for interval_tree in interval_trees.values() for interval in interval_tree}


class TRExplorerCatalogMatchingTests(unittest.TestCase):

    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.catalog_path = Path(self.temp_dir.name) / "catalog.json"
        self.catalog_path.write_text(json.dumps([
            {"LocusId": "1-100-130-CAG", "ReferenceRegion": "chr1:100-130", "LocusStructure": "(CAG)*"},
            # The same interval and motif as the locus above, so the two tie on Jaccard similarity
            {"LocusId": "1-100-130-AGC", "ReferenceRegion": "chr1:100-130", "LocusStructure": "(AGC)*"},
            {"LocusId": "1-200-260-AAGGAT", "ReferenceRegion": ["chr1:200-260", "chr1:260-270"],
             "LocusStructure": "(AAGGAT)*(A)*", "ReferenceMotif": "AAGGAT"},
            {"LocusId": "1-300-340-complex", "ReferenceRegion": "chr1:300-340", "LocusStructure": "GC(CAG)*"},
            {"LocusId": "2-100-150-ACGTACGT", "ReferenceRegion": "chr2:100-150", "LocusStructure": "(ACGTACGT)*"},
        ]))

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_parse_catalog_record_uses_the_first_region_and_reference_motif(self):
        self.assertEqual(parse_catalog_record({
            "LocusId": "1-200-260-AAGGAT", "ReferenceRegion": ["chr1:200-260", "chr1:260-270"],
            "LocusStructure": "(AAGGAT)*(A)*", "ReferenceMotif": "AAGGAT"}), ("1-200-260-AAGGAT", "1", 200, 260, "AAGGAT"))
        self.assertEqual(parse_catalog_record({
            "LocusId": "1-100-130-CAG", "ReferenceRegion": "chr1:100-130", "LocusStructure": "(CAG)*"}),
            ("1-100-130-CAG", "1", 100, 130, "CAG"))
        self.assertIsNone(parse_catalog_record({
            "LocusId": "1-300-340-complex", "ReferenceRegion": "chr1:300-340", "LocusStructure": "GC(CAG)*"}))

    def test_build_catalog_interval_trees_skips_complex_loci_and_keeps_only_the_loci_asked_for(self):
        interval_trees = build_catalog_interval_trees(str(self.catalog_path), show_progress_bar=False)
        self.assertEqual(get_locus_ids(interval_trees),
                         {"1-100-130-CAG", "1-100-130-AGC", "1-200-260-AAGGAT", "2-100-150-ACGTACGT"})
        # Half-open intervals: chr1:130-200 touches the CAG loci and the AAGGAT locus without overlapping them
        interval_trees = build_catalog_interval_trees(
            str(self.catalog_path), show_progress_bar=False, only_loci_that_overlap=[("1", 130, 200), ("2", 149, 160)])
        self.assertEqual(get_locus_ids(interval_trees), {"2-100-150-ACGTACGT"})
        (interval,) = interval_trees["2"]
        self.assertEqual(interval.data, {"locus_id": "2-100-150-ACGTACGT", "start_0based": 100, "end_1based": 150,
                                         "canonical_motif": "ACGTACGT", "motif_length": 8, "catalog_index": 4})

    def test_is_matching_motif_compares_canonical_motifs_up_to_6_bp_and_lengths_above(self):
        self.assertTrue(is_matching_motif("GCT", "AGC", {"canonical_motif": "AGC", "motif_length": 3}))
        self.assertFalse(is_matching_motif("GT", "AC", {"canonical_motif": "AGC", "motif_length": 3}))
        self.assertTrue(is_matching_motif("AAAAAAAT", "AAAAAAAT", {"canonical_motif": "ACGTACGT", "motif_length": 8}))

    def test_find_exact_or_best_jaccard_match(self):
        interval_trees = build_catalog_interval_trees(str(self.catalog_path), show_progress_bar=False)
        # The same ID wins over a tie that would go to the first locus in the catalog
        self.assertEqual(find_exact_or_best_jaccard_match(interval_trees, "1", 100, 130, "AGC", 0.2),
                         ("1-100-130-AGC", "exact"))
        # Without an ID match, the tie between the two loci goes to the first one in the catalog. Having the same
        # interval and canonical motif, it is the same locus, so an exact match even without the same ID.
        self.assertEqual(find_exact_or_best_jaccard_match(interval_trees, "1", 100, 130, "CTG", 0.2),
                         ("1-100-130-CAG", "exact"))
        self.assertEqual(find_exact_or_best_jaccard_match(interval_trees, "1", 103, 130, "GCA", 0.2),
                         ("1-100-130-CAG", "fuzzy"))
        # Jaccard 10/30 is above a 0.2 cutoff but below 0.66
        self.assertEqual(find_exact_or_best_jaccard_match(interval_trees, "1", 120, 130, "CAG", 0.2),
                         ("1-100-130-CAG", "fuzzy"))
        self.assertEqual(find_exact_or_best_jaccard_match(interval_trees, "1", 120, 130, "CAG", 0.66), (None, None))
        # A different motif, and a chromosome without loci
        self.assertEqual(find_exact_or_best_jaccard_match(interval_trees, "1", 100, 130, "AT", 0.2), (None, None))
        self.assertEqual(find_exact_or_best_jaccard_match(interval_trees, "3", 100, 130, "CAG", 0.2), (None, None))


if __name__ == "__main__":
    unittest.main()
