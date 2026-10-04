"""Tests for build_accuracy_table.py: how per-sample arrays add up into per-locus rates, and the
cohort checks that stop a partial or double-counted set of samples from passing as a complete one."""

import gzip
import shutil
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

import build_accuracy_table


CATALOG_ROWS = [
    ("chr1", 1000, 1030, "CAG"),
    ("chr1", 2000, 2030, "AT"),
    ("chr1", 994, 1036, "CAG"),
]
ARRAY_NAMES = build_accuracy_table.ARRAY_NAMES


def write_catalog_bed(path):
    with gzip.open(path, "wt") as f:
        for chrom, start, end, motif in CATALOG_ROWS:
            f.write(f"{chrom}\t{start}\t{end}\t{motif}\t{(end - start) / len(motif):.1f}\n")


def write_comparison_npz(path, **per_locus):
    """per_locus: {array name: list of one value per catalog row}."""
    arrays = {name: np.array(per_locus.get(name, [0] * len(CATALOG_ROWS)),
                             dtype=np.float32 if name in build_accuracy_table.FLOAT_ARRAYS
                             else np.uint8)
              for name in ARRAY_NAMES}
    arrays.update({name: np.array(per_locus.get(name, [np.nan] * len(CATALOG_ROWS)), dtype=np.float32)
                   for name in build_accuracy_table.PER_ALLELE_ARRAYS})
    np.savez_compressed(path, **arrays)


def scored_alleles(truth_short, truth_long, eh_short, eh_long):
    """Arrays for one sample scored at the first catalog row only, with the given allele sizes."""
    return {"has_call": [1, 0, 0], "has_truth": [1, 0, 0], "n_alleles_compared": [2, 0, 0],
            "truth_short_allele": [truth_short, np.nan, np.nan], "truth_long_allele": [truth_long, np.nan, np.nan],
            "eh_short_allele": [eh_short, np.nan, np.nan], "eh_long_allele": [eh_long, np.nan, np.nan]}


def read_table(path):
    with gzip.open(path, "rt") as f:
        header = next(f).rstrip("\n").split("\t")
        return [dict(zip(header, line.rstrip("\n").split("\t"))) for line in f]


class BuildAccuracyTableTests(unittest.TestCase):

    def setUp(self):
        self.dir = Path(tempfile.mkdtemp())
        self.catalog_bed = self.dir / "catalog.bed.gz"
        self.comparison_dir = self.dir / "comparisons"
        self.comparison_dir.mkdir()
        write_catalog_bed(self.catalog_bed)

    def tearDown(self):
        shutil.rmtree(self.dir)

    def write_sample_table(self, sample_ids, sample_id_by_label=None):
        path = self.dir / "samples.tsv"
        with open(path, "w") as f:
            if sample_id_by_label:
                f.write("sample_label\tsample_id\tsex\n")
                for label, sample_id in sample_id_by_label.items():
                    f.write(f"{label}\t{sample_id}\tfemale\n")
                return path
            f.write("sample_id\tGender\n")
            for sample_id in sample_ids:
                f.write(f"{sample_id}\tfemale\n")
        return path

    def run_build(self, samples, expected_sample_ids=None, extra_args=(), sample_id_by_label=None):
        """samples: {sample label: {array name: [values per catalog row]}}"""
        for sample_id, arrays in samples.items():
            write_comparison_npz(self.comparison_dir / f"{sample_id}.comparison.npz", **arrays)
        sample_table = self.write_sample_table(
            expected_sample_ids if expected_sample_ids is not None else sorted(samples), sample_id_by_label)
        by_locus = self.dir / "by_locus.tsv.gz"
        argv = sys.argv
        sys.argv = ["build_accuracy_table.py",
                    "--catalog-bed", str(self.catalog_bed),
                    "--comparison-dir", str(self.comparison_dir),
                    "--by-locus-path", str(by_locus),
                    "--sample-table-path", str(sample_table),
                    *extra_args]
        try:
            build_accuracy_table.main()
        finally:
            sys.argv = argv
        return read_table(by_locus)

    def test_rates_are_over_the_samples_that_were_scored(self):
        samples = {
            "S1": {"has_call": [1, 0, 0], "has_truth": [1, 0, 0], "is_exact": [1, 0, 0],
                   "is_within1": [1, 0, 0]},
            "S2": {"has_call": [1, 0, 0], "has_truth": [1, 0, 0], "is_exact": [0, 0, 0],
                   "is_within1": [1, 0, 0], "abs_error": [1.0, 0, 0]},
        }
        row = {r["locus_id"]: r for r in self.run_build(samples)}["chr1-1000-1030-CAG"]
        self.assertEqual(row["n_samples_with_call"], "2")
        self.assertEqual(row["n_compared"], "2")
        self.assertEqual(row["fraction_exact"], "0.5000")
        self.assertEqual(row["fraction_within1"], "1.0000")
        self.assertEqual(row["mean_abs_error"], "0.5000")

    def test_a_call_without_truth_counts_toward_pok_but_not_toward_a_rate(self):
        samples = {
            "S1": {"has_call": [1, 0, 0], "pok_sum": [0.4, 0, 0]},
            "S2": {"has_call": [1, 0, 0], "has_truth": [1, 0, 0], "is_exact": [1, 0, 0],
                   "pok_sum": [0.8, 0, 0], "pok_sum_compared": [0.8, 0, 0]},
        }
        row = {r["locus_id"]: r for r in self.run_build(samples)}["chr1-1000-1030-CAG"]
        self.assertEqual(row["n_samples_with_call"], "2")
        self.assertEqual(row["mean_pOk"], "0.6000")
        self.assertEqual(row["n_compared"], "1")
        self.assertEqual(row["fraction_exact"], "1.0000")
        # The confidence reported next to that rate comes from the scored sample alone.
        self.assertEqual(row["mean_pOk_on_compared"], "0.8000")

    def test_loci_never_genotyped_are_left_out(self):
        rows = self.run_build({"S1": {"has_call": [1, 0, 1], "has_truth": [1, 0, 1],
                                      "is_exact": [1, 0, 0]}})
        self.assertEqual({r["locus_id"] for r in rows},
                         {"chr1-1000-1030-CAG", "chr1-994-1036-CAG"})

    def test_motif_size_comes_from_the_catalog(self):
        rows = {r["locus_id"]: r for r in self.run_build({"S1": {"has_call": [1, 1, 0]}})}
        self.assertEqual(rows["chr1-1000-1030-CAG"]["motif_size"], "3")
        self.assertEqual(rows["chr1-2000-2030-AT"]["motif_size"], "2")

    def test_a_locus_with_no_scored_sample_reports_no_rate(self):
        rows = {r["locus_id"]: r for r in self.run_build({"S1": {"has_call": [1, 0, 0],
                                                                 "pok_sum": [0.3, 0, 0]}})}
        self.assertEqual(rows["chr1-1000-1030-CAG"]["n_compared"], "0")
        self.assertEqual(rows["chr1-1000-1030-CAG"]["fraction_exact"], "")
        self.assertEqual(rows["chr1-1000-1030-CAG"]["mean_pOk"], "0.3000")

    def test_a_missing_sample_stops_the_build_unless_allowed(self):
        samples = {"S1": {"has_call": [1, 0, 1], "has_truth": [1, 0, 1], "is_exact": [1, 0, 1]}}
        with self.assertRaises(SystemExit):
            self.run_build(samples, expected_sample_ids=["S1", "S2"])
        self.assertTrue(self.run_build(samples, expected_sample_ids=["S1", "S2"],
                                       extra_args=("--allow-partial-cohort",)))

    def test_a_sample_with_two_comparison_files_stops_the_build(self):
        (self.comparison_dir / "nested").mkdir()
        write_comparison_npz(self.comparison_dir / "nested" / "S1.comparison.npz", has_call=[1, 0, 1])
        with self.assertRaises(SystemExit):
            self.run_build({"S1": {"has_call": [1, 0, 1]}})

    def test_a_sample_outside_the_cohort_stops_the_build(self):
        write_comparison_npz(self.comparison_dir / "S9.comparison.npz", has_call=[1, 0, 1])
        with self.assertRaises(SystemExit):
            self.run_build({"S1": {"has_call": [1, 0, 1]}}, expected_sample_ids=["S1"])

    def test_fraction_of_alleles_within_10_percent(self):
        samples = {"S1": {"has_call": [1, 0, 0], "has_truth": [1, 0, 0], "n_alleles_compared": [2, 0, 0],
                          "n_alleles_within_10_percent": [1, 0, 0]},
                   "S2": {"has_call": [1, 0, 0], "has_truth": [1, 0, 0], "n_alleles_compared": [1, 0, 0],
                          "n_alleles_within_10_percent": [1, 0, 0]}}
        row = {r["locus_id"]: r for r in self.run_build(samples)}["chr1-1000-1030-CAG"]
        self.assertEqual(row["n_alleles_compared"], "3")
        self.assertEqual(row["fraction_alleles_within_10_percent"], "0.6667")

    def test_score_uses_the_10_longest_truth_alleles_across_samples(self):
        # 12 samples with truth alleles 10/20, 11/21, ... 21/31: the 10 longest are 22..31. EH calls
        # every long allele exactly except the longest (31), which it calls 20.
        samples = {f"S{i:02d}": scored_alleles(10 + i, 20 + i, 10 + i, 20 + i if i < 11 else 20)
                   for i in range(12)}
        row = {r["locus_id"]: r for r in self.run_build(samples)}["chr1-1000-1030-CAG"]
        self.assertEqual(row["n_longest_truth_alleles_scored"], "10")
        self.assertEqual(row["shortest_of_the_longest_truth_alleles"], "22")
        self.assertEqual(row["fraction_of_longest_truth_alleles_within_10_percent"], "0.9000")
        self.assertEqual(row["mean_abs_error_of_longest_truth_alleles"], "1.1000")

    def test_score_uses_all_alleles_when_fewer_than_10_were_scored(self):
        row = {r["locus_id"]: r for r in self.run_build({"S1": scored_alleles(20, 30, 20, 33)})}[
            "chr1-1000-1030-CAG"]
        self.assertEqual(row["n_longest_truth_alleles_scored"], "2")
        self.assertEqual(row["fraction_of_longest_truth_alleles_within_10_percent"], "1.0000")

    def test_haploid_call_counts_once_toward_the_score(self):
        arrays = scored_alleles(7, 7, 9, 9)
        arrays["n_alleles_compared"] = [1, 0, 0]
        row = {r["locus_id"]: r for r in self.run_build({"S1": arrays})}["chr1-1000-1030-CAG"]
        self.assertEqual(row["n_longest_truth_alleles_scored"], "1")
        self.assertEqual(row["fraction_of_longest_truth_alleles_within_10_percent"], "0.0000")

    def test_one_genome_sequenced_at_several_depths_counts_once(self):
        samples = {"HG002": scored_alleles(20, 30, 20, 30), "HG002_10x": scored_alleles(20, 30, 5, 5),
                   "S1": scored_alleles(20, 30, 20, 30)}
        rows = self.run_build(samples, sample_id_by_label={"HG002": "HG002", "HG002_10x": "HG002", "S1": "S1"})
        row = {r["locus_id"]: r for r in rows}["chr1-1000-1030-CAG"]
        self.assertEqual(row["n_compared"], "2")
        self.assertEqual(row["fraction_of_longest_truth_alleles_within_10_percent"], "1.0000")

    def test_comparison_file_from_an_older_script_is_an_error(self):
        path = self.comparison_dir / "S1.comparison.npz"
        np.savez_compressed(path, **{name: np.zeros(len(CATALOG_ROWS), dtype=np.uint8)
                                     for name in ("has_call", "has_truth", "is_exact")})
        with self.assertRaisesRegex(ValueError, "force-comparison-step"):
            self.run_build({}, expected_sample_ids=["S1"])

    def test_mismatched_array_length_is_an_error(self):
        path = self.comparison_dir / "S1.comparison.npz"
        np.savez_compressed(path, **{name: np.zeros(2, dtype=np.uint8)
                                     for name in ARRAY_NAMES + build_accuracy_table.PER_ALLELE_ARRAYS})
        with self.assertRaisesRegex(ValueError, "not built from the same catalog"):
            self.run_build({}, expected_sample_ids=["S1"])


if __name__ == "__main__":
    unittest.main()
