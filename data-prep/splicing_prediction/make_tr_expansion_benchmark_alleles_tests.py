import unittest

import gzip
import os
import tempfile

from make_tr_expansion_benchmark_alleles import (
    allele_to_variant, canonical_motif, is_inside_a_transcript, load_merged_transcript_spans, percentile_of_histogram,
    TARGET_LABELS, labeled_size_changes_in_repeat_units, short_target_label, target_lengths, target_sort_key)


class TranscriptOverlapTests(unittest.TestCase):

    def test_matches_the_spliceai_annotator_rule(self):
        # TX_START is 0-based: a transcript with TX_START 99, TX_END 200 covers positions 100-200
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "annotation.txt.gz")
            with gzip.open(path, "wt") as f:
                f.write("#NAME\tCHROM\tSTRAND\tTX_START\tTX_END\tEXON_START\tEXON_END\n")
                f.write("T1\tchr1\t+\t99\t200\t99,\t200,\n")
                f.write("T2\tchr1\t+\t150\t300\t150,\t300,\n")
                f.write("T3\tchr1\t+\t500\t600\t500,\t600,\n")
            spans = load_merged_transcript_spans(path)
        self.assertEqual(spans["chr1"], ([100, 501], [300, 600]))
        self.assertFalse(is_inside_a_transcript(spans, "chr1", 99))
        self.assertTrue(is_inside_a_transcript(spans, "chr1", 100))
        self.assertTrue(is_inside_a_transcript(spans, "chr1", 300))
        self.assertFalse(is_inside_a_transcript(spans, "chr1", 301))
        self.assertFalse(is_inside_a_transcript(spans, "chr2", 150))


class PercentileOfHistogramTests(unittest.TestCase):

    def test_nearest_rank(self):
        # 512 haplotypes: 25 at 11 bp, 433 at 12, 46 at 13, 7 at 14 and one 8085 bp outlier
        counts = {11: 25, 12: 433, 13: 46, 14: 7, 8085: 1}
        # rank ceil(0.025 * 512) = 13 falls among the 25 at 11 bp; rank 500 among the 46 at 13 bp; rank 510 at 14 bp
        self.assertEqual([percentile_of_histogram(counts, p) for p in (2.5, 97.5, 99.5)], [11, 13, 14])
        self.assertEqual(percentile_of_histogram(counts, 100), 8085)

    def test_monomorphic(self):
        self.assertEqual(percentile_of_histogram({9: 512}, 2.5), 9)


class TargetLengthsTests(unittest.TestCase):

    def test_observed_percentiles_then_the_99_5th_plus_multiples_of_the_motifs_range(self):
        # a motif range of 11 bp: 14 + 1, 2 and 3 x 11
        self.assertEqual(target_lengths(11, 13, 14, 11), [11, 13, 14, 25, 36, 47])

    def test_short_labels_and_order(self):
        self.assertEqual([short_target_label(t) for t in TARGET_LABELS], ["2.5th", "97.5th", "99.5th", "+1x", "+2x", "+3x"])
        self.assertEqual(sorted(reversed(TARGET_LABELS + ("99.5th + 0.5x motif range",)), key=target_sort_key)[3],
                         "99.5th + 0.5x motif range")


class CanonicalMotifTests(unittest.TestCase):

    def test_rotations_and_reverse_complements_share_one_motif(self):
        self.assertEqual({canonical_motif(m) for m in ("CAG", "AGC", "GCA", "CTG", "TGC", "GCT")}, {"AGC"})
        self.assertEqual(canonical_motif("TTTCA"), canonical_motif("GAAAT"))


class SizeChangesInRepeatUnitsTests(unittest.TestCase):

    def test_rounds_to_whole_units_drops_the_reference_and_merges_labels(self):
        # 12 bp of a 3 bp motif; targets 11, 12, 13, 14 and 20 bp are -1/3, 0, +1/3, +2/3 and +8/3 units
        self.assertEqual(labeled_size_changes_in_repeat_units(12, 3, [11, 12, 13, 14, 20], labels="abcde"),
                         {1: ["d"], 3: ["e"]})
        self.assertEqual(labeled_size_changes_in_repeat_units(12, 3, [15, 16], labels="ab"), {1: ["a", "b"]})

    def test_contractions_keep_at_least_one_unit(self):
        self.assertEqual(list(labeled_size_changes_in_repeat_units(6, 2, [0, 2, 4])), [-2, -1])

    def test_expansions_insert_at_most_5000_bp(self):
        self.assertEqual(list(labeled_size_changes_in_repeat_units(12, 3, [12, 9000])), [1666])

    def test_contractions_delete_at_most_5000_bp(self):
        # 12,375 bp of a 5 bp motif contracted to 7 bp: capped at 5000 // 5 = 1000 units deleted
        self.assertEqual(list(labeled_size_changes_in_repeat_units(12375, 5, [7])), [-1000])

    def test_exact_half_unit_with_a_long_motif_rounds_up(self):
        # chr1-90838466-90838559, a 38 bp motif: reference 93 bp, 99.5th percentile 138 bp, motif range 424 bp;
        # +2x = 138 + 848 = 986 bp is +893 bp = exactly +23.5 units -> +24 (float math gave 23.4999... -> +23)
        self.assertEqual(list(labeled_size_changes_in_repeat_units(93, 38, target_lengths(118, 130, 138, 424)[4:5])), [24])

    def test_exact_halves_round_away_from_zero(self):
        # 20 bp of a 2 bp motif: 23 and 25 bp are +1.5 and +2.5 units -> +2 and +3, not both +2
        self.assertEqual(labeled_size_changes_in_repeat_units(20, 2, [23, 25], labels="ab"), {2: ["a"], 3: ["b"]})
        self.assertEqual(labeled_size_changes_in_repeat_units(20, 2, [19, 17], labels="ab"), {-2: ["b"], -1: ["a"]})


class AlleleToVariantTests(unittest.TestCase):

    def test_expansion_inserts_units_after_the_anchor(self):
        self.assertEqual(allele_to_variant("chr1", 100, "T", "CAGCAGCAG", "CAG", 3, 5), "chr1-100-T-TCAGCAG")

    def test_contraction_deletes_the_first_units_of_the_tract(self):
        self.assertEqual(allele_to_variant("chr1", 100, "T", "CAGCAGCAG", "CAG", 3, 1), "chr1-100-TCAGCAG-T")


if __name__ == "__main__":
    unittest.main()
