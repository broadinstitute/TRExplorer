"""Tests for compare_eh_to_truth.py."""

import gzip
import json
import shutil
import tempfile
import unittest
from pathlib import Path

import numpy as np

import compare_eh_to_truth


CATALOG_ROWS = [
    # chrom, start_0based, end_1based, motif -- LocusId is "{chrom}-{start}-{end}-{motif}"
    ("chr1", 1000, 1030, "CAG"),   # inside the high-confidence region
    ("chr1", 2000, 2030, "AT"),    # inside the high-confidence region
    ("chr1", 9000, 9030, "CAG"),   # outside every high-confidence region
    ("chrX", 3000, 3030, "GCC"),   # inside, used for the haploid cases
    ("chrX", 4000, 4030, "TTA"),   # inside, used for the haploid cases
]

HIGH_CONFIDENCE_REGIONS = [("chr1", 500, 5000), ("chrX", 2500, 5000)]


def locus_id(row):
    chrom, start, end, motif = row
    return f"{chrom}-{start}-{end}-{motif}"


def write_catalog_bed(path):
    with gzip.open(path, "wt") as f:
        for chrom, start, end, motif in CATALOG_ROWS:
            f.write(f"{chrom}\t{start}\t{end}\t{motif}\t{(end - start) / len(motif):.1f}\n")


def write_high_confidence_bed(path):
    with gzip.open(path, "wt") as f:
        for chrom, start, end in HIGH_CONFIDENCE_REGIONS:
            f.write(f"{chrom}\t{start}\t{end}\n")


def direction_probabilities(pok):
    """Return (pTooShort, pTooLong) splitting the non-OK mass 2:1, so the three sum to 1 as EH's do.

    The two differ from each other and from pOk at every pOk used in these tests, so a test that
    reads the wrong array fails rather than passing on a coincidence.
    """
    too_short = round((1.0 - pok) * 2.0 / 3.0, 3)
    return too_short, round(1.0 - pok - too_short, 3)


def write_eh_json(path, genotypes_and_poks, quick_loci=()):
    """genotypes_and_poks: {LocusId: (genotype string or None, [pOk, ...])}

    quick_loci: the LocusIds ExpansionHunter genotyped on its fast path.
    """
    locus_results = {}
    for lid, (genotype, poks) in genotypes_and_poks.items():
        # Shaped like real ExpansionHunter output: LocusId appears once, at the locus level, and
        # the variant record carries a VariantId instead.
        locus_results[lid] = {
            "LocusId": lid,
            "Coverage": 30.5,
            "Variants": {
                lid: {
                    "VariantId": lid,
                    "Genotype": genotype,
                    "QuickGenotype": lid in quick_loci,
                    "CountsOfSpanningReads": "(10, 12), (11, 3)",
                    "CountsOfFlankingReads": "(9, 4)",
                    "CountsOfInrepeatReads": "()",
                    "AlleleQualityMetrics": {
                        "Alleles": [allele_record(i, pok) for i, pok in enumerate(poks)],
                        "VariantId": lid,
                    },
                }
            },
        }
    with gzip.open(path, "wt") as f:
        json.dump({"LocusResults": locus_results}, f, indent=4)


def allele_record(allele_index, pok):
    """One AlleleQualityMetrics entry carrying all three genotype-quality probabilities."""
    too_short, too_long = direction_probabilities(pok)
    return {"AlleleNumber": allele_index + 1, "pOk": pok,
            "pTooShort": too_short, "pTooLong": too_long}


def write_truth_tsv(path, truth_by_locus_id, reference_copies=10):
    """truth_by_locus_id: {LocusId: (short allele, long allele)}; "" for a missing genotype."""
    columns = ["Chrom", "Start0Based", "End", "LocusId", "Motif", "NumRepeatsInReference",
               "NumRepeatsShortAllele", "NumRepeatsLongAllele"]
    with gzip.open(path, "wt") as f:
        f.write("\t".join(columns) + "\n")
        for row in CATALOG_ROWS:
            lid = locus_id(row)
            if lid not in truth_by_locus_id:
                continue
            short_allele, long_allele = truth_by_locus_id[lid]
            chrom, start, end, motif = row
            f.write("\t".join(str(v) for v in
                              [chrom, start, end, lid, motif, reference_copies,
                               short_allele, long_allele]) + "\n")


class CompareEhToTruthTests(unittest.TestCase):

    def setUp(self):
        self.dir = Path(tempfile.mkdtemp())
        self.catalog_bed = self.dir / "catalog.bed.gz"
        self.high_confidence_bed = self.dir / "sample.dip.bed.gz"
        self.eh_json = self.dir / "sample.json.gz"
        self.truth_tsv = self.dir / "sample.tsv.gz"
        write_catalog_bed(self.catalog_bed)
        write_high_confidence_bed(self.high_confidence_bed)

    def tearDown(self):
        shutil.rmtree(self.dir)

    def compare(self, genotypes_and_poks, truth_by_locus_id, quick_loci=()):
        write_eh_json(self.eh_json, genotypes_and_poks, quick_loci)
        write_truth_tsv(self.truth_tsv, truth_by_locus_id)
        return compare_eh_to_truth.compare_sample(
            self.catalog_bed, self.eh_json, self.truth_tsv, self.high_confidence_bed)

    def test_catalog_row_numbers_follow_file_order(self):
        row_numbers = compare_eh_to_truth.catalog_row_numbers(self.catalog_bed)
        self.assertEqual(row_numbers[locus_id(CATALOG_ROWS[0])], 0)
        self.assertEqual(row_numbers[locus_id(CATALOG_ROWS[4])], 4)
        self.assertEqual(len(row_numbers), len(CATALOG_ROWS))

    def test_exact_match(self):
        arrays = self.compare({locus_id(CATALOG_ROWS[0]): ("10/10", [0.9, 0.9])},
                              {locus_id(CATALOG_ROWS[0]): (10, 10)})
        self.assertEqual(arrays["has_call"][0], 1)
        self.assertEqual(arrays["has_truth"][0], 1)
        self.assertEqual(arrays["is_exact"][0], 1)
        self.assertEqual(arrays["is_within1"][0], 1)
        self.assertEqual(arrays["abs_error"][0], 0)
        self.assertAlmostEqual(float(arrays["pok_sum"][0]), 0.9, places=5)
        self.assertAlmostEqual(float(arrays["pok_sum_compared"][0]), 0.9, places=5)

    def test_off_by_one_is_near_but_not_exact(self):
        arrays = self.compare({locus_id(CATALOG_ROWS[0]): ("10/11", [0.8, 0.7])},
                              {locus_id(CATALOG_ROWS[0]): (10, 10)})
        self.assertEqual(arrays["is_exact"][0], 0)
        self.assertEqual(arrays["is_within1"][0], 1)
        self.assertAlmostEqual(float(arrays["abs_error"][0]), 0.5, places=5)
        self.assertAlmostEqual(float(arrays["pok_sum"][0]), 0.75, places=5)

    def test_off_by_two_is_neither(self):
        arrays = self.compare({locus_id(CATALOG_ROWS[1]): ("8/12", [0.5])},
                              {locus_id(CATALOG_ROWS[1]): (10, 10)})
        self.assertEqual(arrays["is_exact"][1], 0)
        self.assertEqual(arrays["is_within1"][1], 0)
        self.assertAlmostEqual(float(arrays["abs_error"][1]), 2.0, places=5)

    def test_locus_outside_high_confidence_regions_has_no_truth(self):
        outside = locus_id(CATALOG_ROWS[2])
        arrays = self.compare({outside: ("10/10", [0.9])}, {outside: (10, 10)})
        self.assertEqual(arrays["has_call"][2], 1)
        self.assertEqual(arrays["has_truth"][2], 0)
        self.assertEqual(arrays["is_exact"][2], 0)
        # pOk still counts toward the per-call average, but not toward the compared-cohort average.
        self.assertAlmostEqual(float(arrays["pok_sum"][2]), 0.9, places=5)
        self.assertEqual(float(arrays["pok_sum_compared"][2]), 0)

    def test_empty_truth_genotype_is_missing_not_zero(self):
        lid = locus_id(CATALOG_ROWS[0])
        arrays = self.compare({lid: ("10/10", [0.9])}, {lid: ("", "")})
        self.assertEqual(arrays["has_call"][0], 1)
        self.assertEqual(arrays["has_truth"][0], 0)

    def test_haploid_call_scored_against_homozygous_truth(self):
        lid = locus_id(CATALOG_ROWS[3])
        arrays = self.compare({lid: ("7", [0.6])}, {lid: (7, 7)})
        self.assertEqual(arrays["has_truth"][3], 1)
        self.assertEqual(arrays["is_exact"][3], 1)

    def test_haploid_call_skipped_when_truth_is_heterozygous(self):
        lid = locus_id(CATALOG_ROWS[4])
        arrays = self.compare({lid: ("7", [0.6])}, {lid: (7, 9)})
        self.assertEqual(arrays["has_call"][4], 1)
        self.assertEqual(arrays["has_truth"][4], 0)

    def test_fast_path_genotypes_are_flagged(self):
        quick, full = locus_id(CATALOG_ROWS[0]), locus_id(CATALOG_ROWS[1])
        arrays = self.compare({quick: ("10/10", [0.9]), full: ("5/5", [0.7])},
                              {quick: (10, 10), full: (5, 5)}, quick_loci=(quick,))
        self.assertEqual(arrays["is_quick"][0], 1)
        self.assertEqual(arrays["is_quick"][1], 0)

    def test_read_counts_and_coverage_are_captured(self):
        lid = locus_id(CATALOG_ROWS[0])
        arrays = self.compare({lid: ("10/10", [0.9])}, {lid: (10, 10)})
        self.assertAlmostEqual(float(arrays["coverage"][0]), 30.5, places=3)
        self.assertEqual(int(arrays["n_spanning_reads"][0]), 15)   # (10, 12) + (11, 3)
        self.assertEqual(int(arrays["n_flanking_reads"][0]), 4)
        self.assertEqual(int(arrays["n_irr_reads"][0]), 0)

    def test_distance_from_the_reference_is_recorded_for_call_and_truth(self):
        lid = locus_id(CATALOG_ROWS[0])
        # reference is 10 copies; ExpansionHunter calls 10/12, truth is 10/25
        arrays = self.compare({lid: ("10/12", [0.9])}, {lid: (10, 25)})
        self.assertAlmostEqual(float(arrays["eh_diff_from_ref"][0]), 2.0, places=3)
        self.assertAlmostEqual(float(arrays["truth_short_diff_from_ref"][0]), 0.0, places=3)
        self.assertAlmostEqual(float(arrays["truth_long_diff_from_ref"][0]), 15.0, places=3)

    def test_truth_distance_is_negative_for_a_contraction(self):
        lid = locus_id(CATALOG_ROWS[0])
        # reference is 10 copies; both truth alleles sit below it
        arrays = self.compare({lid: ("10/10", [0.9])}, {lid: (3, 4)})
        self.assertAlmostEqual(float(arrays["truth_short_diff_from_ref"][0]), -7.0, places=3)
        self.assertAlmostEqual(float(arrays["truth_long_diff_from_ref"][0]), -6.0, places=3)

    def test_truth_straddling_the_reference_keeps_both_directions(self):
        lid = locus_id(CATALOG_ROWS[0])
        # reference is 10 copies, the short allele contracted and the long one expanded: the case a
        # single signed distance could not represent, since it would have to hide one of the two
        arrays = self.compare({lid: ("10/10", [0.9])}, {lid: (2, 21)})
        self.assertAlmostEqual(float(arrays["truth_short_diff_from_ref"][0]), -8.0, places=3)
        self.assertAlmostEqual(float(arrays["truth_long_diff_from_ref"][0]), 11.0, places=3)

    def test_distance_from_the_reference_is_nan_without_truth(self):
        lid = locus_id(CATALOG_ROWS[0])
        arrays = self.compare({lid: ("10/12", [0.9])}, {})
        self.assertTrue(np.isnan(arrays["eh_diff_from_ref"][0]))
        self.assertTrue(np.isnan(arrays["truth_short_diff_from_ref"][0]))
        self.assertTrue(np.isnan(arrays["truth_long_diff_from_ref"][0]))

    def test_null_genotype_is_not_a_call(self):
        lid = locus_id(CATALOG_ROWS[0])
        arrays = self.compare({lid: (None, [])}, {lid: (10, 10)})
        self.assertEqual(arrays["has_call"][0], 0)
        self.assertEqual(arrays["has_truth"][0], 0)

    def test_locus_id_absent_from_the_catalog_is_an_error(self):
        with self.assertRaises(ValueError):
            self.compare({"chr9-1-2-CAG": ("3/3", [0.9])}, {})

    def test_direction_probabilities_are_captured_separately(self):
        lid = locus_id(CATALOG_ROWS[0])
        arrays = self.compare({lid: ("10/10", [0.9, 0.6])}, {lid: (10, 10)})
        expected_short = np.mean([direction_probabilities(0.9)[0], direction_probabilities(0.6)[0]])
        expected_long = np.mean([direction_probabilities(0.9)[1], direction_probabilities(0.6)[1]])
        self.assertNotAlmostEqual(expected_short, expected_long, places=3)
        for suffix in ("_sum", "_sum_compared"):
            self.assertAlmostEqual(float(arrays["pok" + suffix][0]), 0.75, places=5)
            self.assertAlmostEqual(float(arrays["ptooshort" + suffix][0]), expected_short, places=5)
            self.assertAlmostEqual(float(arrays["ptoolong" + suffix][0]), expected_long, places=5)

    def test_all_three_probabilities_parsed_from_one_line(self):
        """An allele record written on a single line yields all three, not just the first."""
        path = self.dir / "one_line_alleles.json"
        with open(path, "w") as f:
            f.write('{"LocusResults": {"chr1-100-200-CAG": {\n'
                    '"LocusId": "chr1-100-200-CAG",\n'
                    '"Variants": {"v": {\n'
                    '"Genotype": "10/10",\n'
                    '"AlleleQualityMetrics": {"Alleles": [\n'
                    '{"AlleleNumber": 1, "pOk": 0.9, "pTooShort": 0.067, "pTooLong": 0.033}\n'
                    ']}}}}}}\n')
        records = dict(compare_eh_to_truth.parse_eh_json(path))
        record = records["chr1-100-200-CAG"]
        self.assertEqual(record["pok_values"], [0.9])
        self.assertEqual(record["ptooshort_values"], [0.067])
        self.assertEqual(record["ptoolong_values"], [0.033])

    def test_direction_probabilities_are_not_counted_without_truth(self):
        lid = locus_id(CATALOG_ROWS[2])
        arrays = self.compare({lid: ("10/10", [0.9])}, {})
        expected_short, expected_long = direction_probabilities(0.9)
        self.assertAlmostEqual(float(arrays["ptooshort_sum"][2]), expected_short, places=5)
        self.assertEqual(float(arrays["ptooshort_sum_compared"][2]), 0)
        self.assertAlmostEqual(float(arrays["ptoolong_sum"][2]), expected_long, places=5)
        self.assertEqual(float(arrays["ptoolong_sum_compared"][2]), 0)

    def test_loci_with_no_expansion_hunter_record_stay_zero(self):
        arrays = self.compare({locus_id(CATALOG_ROWS[0]): ("10/10", [0.9])},
                              {locus_id(CATALOG_ROWS[0]): (10, 10)})
        self.assertEqual(arrays["has_call"].sum(), 1)
        self.assertEqual(arrays["has_call"][1], 0)
        self.assertEqual(arrays["pok_sum"][1], 0)


if __name__ == "__main__":
    unittest.main()
