import unittest

from summarize_full_run_spliceai_scores_per_locus import allele_size_name, summarize_locus

LABELS = {"2.5": "2.5th percentile", "97.5": "97.5th percentile", "99.5": "99.5th percentile",
          "+1x": "99.5th + 1x motif range", "+2x": "99.5th + 2x motif range", "+3x": "99.5th + 3x motif range"}


def record(ag, al, dg, dl):
    return {"DS_AG": ag, "DS_AL": al, "DS_DG": dg, "DS_DL": dl}


class AlleleSizeNameTests(unittest.TestCase):

    def test_names(self):
        self.assertEqual([allele_size_name(label) for label in LABELS.values()],
                         ["2.5pct", "97.5pct", "99.5pct", "+1x", "+2x", "+3x"])


class SummarizeLocusTests(unittest.TestCase):

    def test_max_score_min_affecting_size_and_repeat_count_string(self):
        # 20 bp of CA = 10 units; alleles: -1 unit (2.5pct), +2 units (97.5pct and 99.5pct merged), +5 (+1x), +15 (+3x)
        alleles = ["chr1-100-ACA-A", "chr1-100-A-ACACA", "chr1-100-A-A" + "CA" * 5, "chr1-100-A-A" + "CA" * 15]
        labels = [[LABELS["2.5"]], [LABELS["97.5"], LABELS["99.5"]], [LABELS["+1x"]], [LABELS["+3x"]]]
        records = {alleles[1]: record("0.005", "0.000", "0.000", "0.000"),
                   alleles[2]: record("0.310", "0.020", "0.000", "0.000"),
                   alleles[3]: record("0.000", "0.880", "0.000", "0.000")}
        row = summarize_locus("chr1-100-120-CA", alleles, labels, records)
        self.assertEqual(row["LocusId"], "1-100-120-CA")
        self.assertEqual(row["SpliceAI_MaxDeltaScore"], "0.880")
        self.assertEqual(row["SpliceAI_MaxDeltaScoreAlleleSize"], "+3x")
        self.assertEqual(row["SpliceAI_MinAlleleSizeThatAffectsSplicing"], "+1x")
        self.assertEqual(row["SpliceAI_DeltaScoreByRepeatCount"], "9:0.000,12:0.005,15:0.310AG,25:0.880AL")

    def test_locus_with_no_stored_alleles(self):
        row = summarize_locus("chr2-5-9-A", ["chr2-5-A-AA"], [[LABELS["+1x"]]], {})
        self.assertEqual((row["SpliceAI_MaxDeltaScore"], row["SpliceAI_MaxDeltaScoreAlleleSize"],
                          row["SpliceAI_MinAlleleSizeThatAffectsSplicing"], row["SpliceAI_DeltaScoreByRepeatCount"]),
                         ("0.000", "", "", "5:0.000"))

    def test_alleles_whose_scoring_returned_an_error_are_left_out(self):
        alleles = ["chr1-100-A-ACA", "chr1-100-A-ACACA"]
        records = {alleles[1]: record("0.300", "0.000", "0.000", "0.000")}
        labels = [[LABELS["+1x"]], [LABELS["+2x"]]]
        row = summarize_locus("chr1-100-120-CA", alleles, labels, records, unscored_variants={alleles[0]})
        self.assertEqual(row["SpliceAI_DeltaScoreByRepeatCount"], "12:0.300AG")
        self.assertIsNone(summarize_locus("chr1-100-120-CA", alleles, labels, records, unscored_variants=set(alleles)))

    def test_ties_go_to_the_smallest_size(self):
        alleles = ["chr1-100-A-ACA", "chr1-100-A-ACACA"]
        records = {v: record("0.500", "0.000", "0.000", "0.000") for v in alleles}
        row = summarize_locus("chr1-100-120-CA", alleles, [[LABELS["+1x"]], [LABELS["+2x"]]], records)
        self.assertEqual(row["SpliceAI_MaxDeltaScoreAlleleSize"], "+1x")


if __name__ == "__main__":
    unittest.main()
