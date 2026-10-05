import gzip
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from compute_encode4_columns_for_catalog_loci import (
    build_catalog_index,
    compute_validated_element_labels,
    find_loci_overlapping_elements,
    find_overlapping_pairs,
    parse_crispr_significant_elements_tsv,
    parse_igvf_crispr_significant_elements,
    parse_mpra_active_elements_bed,
    write_output_tsv,
)

# Locus 0 and 1 are short and adjacent, locus 2 is long (2,000 bp), locus 3 is on another chromosome.
CATALOG_BED_LINES = [
    "chr1\t100\t130\tCAG",
    "chr1\t130\t160\tGT",
    "chr1\t1000\t3000\tAAAAAC",
    "chr2\t500\t520\tA",
]


class Encode4ColumnsTests(unittest.TestCase):

    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.dir = Path(self.temp_dir.name)
        self.catalog_bed = self.write_file("catalog.bed.gz", CATALOG_BED_LINES)
        self.n_loci, self.catalog_index = build_catalog_index(self.catalog_bed)

    def tearDown(self):
        self.temp_dir.cleanup()

    def write_file(self, name, lines):
        path = self.dir / name
        fopen = gzip.open if name.endswith(".gz") else open
        with fopen(path, "wt") as f:
            f.write("\n".join(lines) + "\n")
        return path

    def overlapping_loci(self, intervals):
        elements = pd.DataFrame(intervals, columns=["chrom", "start", "end"])
        element_positions, locus_indices = find_loci_overlapping_elements(self.catalog_index, elements)
        return sorted(zip(element_positions.tolist(), locus_indices.tolist()))

    def test_find_overlapping_pairs_uses_half_open_intervals(self):
        starts, ends = np.array([10, 20]), np.array([20, 30])
        query_positions, interval_positions = find_overlapping_pairs(starts, ends, 10, np.array([20]), np.array([21]))
        self.assertEqual(interval_positions.tolist(), [1])
        query_positions, interval_positions = find_overlapping_pairs(starts, ends, 10, np.array([5]), np.array([10]))
        self.assertEqual(interval_positions.tolist(), [])

    def test_find_loci_overlapping_elements(self):
        self.assertEqual(self.overlapping_loci([("chr1", 129, 131)]), [(0, 0), (0, 1)])
        self.assertEqual(self.overlapping_loci([("chr1", 160, 170)]), [])
        # Inside the long locus, far from its start.
        self.assertEqual(self.overlapping_loci([("chr1", 2900, 2950)]), [(0, 2)])
        # An element that spans a whole short locus.
        self.assertEqual(self.overlapping_loci([("chr2", 0, 10000)]), [(0, 3)])
        # Chromosomes missing from the catalog are ignored.
        self.assertEqual(self.overlapping_loci([("chrUn_x", 0, 10000), ("chr1", 150, 151)]), [(1, 1)])

    def test_crispr_parser_keeps_significant_rows(self):
        path = self.write_file("crispr.tsv", [
            "chrom\tchromStart\tchromEnd\tname\tSignificant",
            "chr1\t100\t110\ta\tTRUE",
            "chr1\t200\t210\tb\tFALSE",
        ])
        self.assertEqual(parse_crispr_significant_elements_tsv(path).values.tolist(), [["chr1", 100, 110]])
        no_significant_column = self.write_file("other.tsv", ["chr\tstart\tend\tlog2_fold_change", "chr1\t1\t2\t3.0"])
        self.assertIsNone(parse_crispr_significant_elements_tsv(no_significant_column))

    def test_mpra_parser_applies_fold_change_and_fdr_thresholds(self):
        path = self.write_file("mpra.bed.gz", [
            "chr1\t100\t110\ta\t0\t+\t1.0\t10\t20\t3.0\t2.0",   # active
            "chr1\t120\t130\tb\t0\t+\t-1.0\t10\t5\t3.0\t2.0",   # negative fold change
            "chr1\t140\t150\tc\t0\t+\t1.0\t10\t20\t1.0\t1.0",   # FDR 0.1
        ])
        self.assertEqual(parse_mpra_active_elements_bed(path).values.tolist(), [["chr1", 100, 110]])
        wrong_format = self.write_file("wrong.bed.gz", ["chr1\t100\t110\ta\t0\t+"])
        self.assertIsNone(parse_mpra_active_elements_bed(wrong_format))
        # Files with the placeholder -1 in every q-value have no FDR, so they are skipped rather than
        # silently yielding no elements.
        no_fdr = self.write_file("no_fdr.bed.gz", [
            "chr1\t100\t110\ta\t0\t+\t1.65\t0.3\t0.8\t-1\t-1",
            "chr1\t120\t130\tb\t0\t+\t2.0\t0.3\t0.8\t-1\t-1",
        ])
        self.assertIsNone(parse_mpra_active_elements_bed(no_fdr))

    def igvf_elements(self, name, lines, file_type="tsv"):
        elements = parse_igvf_crispr_significant_elements(self.write_file(name, lines), file_type)
        return None if elements is None else elements.values.tolist()

    def test_igvf_parser_uses_significant_flag_columns(self):
        self.assertEqual(self.igvf_elements("flowfish.tsv", [
            "name\tname_hg38\tadj.pval.EnhancerEffect.noAux\tSignificant",
            "chr8:1-2\tchr8:100-200\t0.9\tFALSE",
            "chr8:3-4\tchr8:300-400\t0.9\tTRUE",   # the flag wins over the adjusted p-value
        ]), [["chr8", 300, 400]])
        # CSV with a byte order mark, as in the IGVF TeloHAEC files.
        self.assertEqual(self.igvf_elements("telohaec.csv", [
            "﻿intended_target_name,targeting_chr,targeting_start,targeting_end,adj_p_value,significant",
            "x,chr15,100,200,0.01,TRUE",
            "y,chr15,300,400,0.6,FALSE",
        ], file_type="csv"), [["chr15", 100, 200]])
        self.assertEqual(self.igvf_elements("growth.tsv", [
            "dhs\tdhs_coords\tgrowth_significant",
            "93\tchr1:101174581-101175330_93\tTRUE",
            "94\tchr1:5-10_94\tFALSE",
            "95\tNA:NA-NA_ABL1\tTRUE",   # gene-level target without coordinates
        ]), [["chr1", 101174581, 101175330]])

    def test_igvf_parser_applies_published_fdr_cutoffs(self):
        self.assertEqual(self.igvf_elements("fractel.tsv", [
            "FRACTEL_pval_fdr_corr\tintended_target_chr\tintended_target_start\tintended_target_end",
            "0.01\tchr1\t100\t200",
            "0.05\tchr1\t300\t400",   # not below the 0.05 cutoff
        ]), [["chr1", 100, 200]])
        self.assertEqual(self.igvf_elements("cohesin.tsv", [
            "name\tname_hg38\tadj.pval.EnhancerEffect.noAux",
            "chr8:1-2\tchr8:100-200\t0.001",
            "chr8:3-4\tchr8:300-400\t0.2",
        ]), [["chr8", 100, 200]])
        self.assertEqual(self.igvf_elements("tapseq.tsv", [
            "element_hg19\telement_hg38\tminus_auxin_padj\tplus_auxin_padj",
            "chr11:1-2\tchr11:100-200\t0.01\t0.9",
        ]), [["chr11", 100, 200]])

    def test_igvf_parser_skips_layouts_without_a_published_rule(self):
        # IDR-based files, files with an unexplained adjusted p-value, and hg19-only files.
        self.assertIsNone(self.igvf_elements("idr.tsv", [
            "intended_target_chr\tintended_target_start\tintended_target_end\tLFC\tZ-score\tIDR", "chr13\t1\t2\t-2\t-11\t0.00001"]))
        self.assertIsNone(self.igvf_elements("perturbseq.tsv", [
            "p_val_adj\tintended_target_name\tintended_target_chr\tintended_target_start\tintended_target_end", "0.01\tx\tchr1\t1\t2"]))
        self.assertIsNone(self.igvf_elements("hg19.tsv", [
            "minus_auxin_padj\tplus_auxin_padj\tperturbation_name", "0.01\t0.5\tchr11:1-2_ANO1"]))

    def test_validated_element_labels_and_output_tsv(self):
        crispr = self.write_file("crispr.tsv", [
            "chrom\tchromStart\tchromEnd\tname\tSignificant", "chr1\t150\t155\ta\tTRUE"])
        mpra = self.write_file("mpra.bed.gz", ["chr1\t100\t110\ta\t0\t+\t1.0\t10\t20\t3.0\t2.0"])
        igvf = self.write_file("igvf.tsv", [
            "FRACTEL_pval_fdr_corr\tintended_target_chr\tintended_target_start\tintended_target_end", "0.01\tchr1\t100\t101"])
        rows = [
            {"data_source": "ENCODE", "assay_label": "CRISPR", "biosample_term_name": "HepG2", "local_path": crispr},
            {"data_source": "ENCODE", "assay_label": "MPRA", "biosample_term_name": "K562", "local_path": mpra},
            {"data_source": "IGVF", "assay_label": "CRISPR", "biosample_term_name": "HCT116", "local_path": igvf,
             "file_type": "tsv"},
        ]
        validated_labels = compute_validated_element_labels(self.catalog_index, rows)
        self.assertEqual(validated_labels, {0: "ENCODE:MPRA:K562,IGVF:CRISPR:HCT116", 1: "ENCODE:CRISPR:HepG2"})

        output_tsv = self.dir / "out.tsv.gz"
        n_loci_with_value = write_output_tsv(self.catalog_bed, output_tsv, validated_labels)
        output = pd.read_csv(output_tsv, sep="\t", dtype=str, keep_default_na=False)
        # Loci 2 and 3 have no values, so they are left out.
        self.assertEqual(output.values.tolist(), [
            ["1-100-130-CAG", "ENCODE:MPRA:K562,IGVF:CRISPR:HCT116"], ["1-130-160-GT", "ENCODE:CRISPR:HepG2"]])
        self.assertEqual(n_loci_with_value["ENCODE4_and_IGVF_CRISPR_or_MPRA_validated_elements"], 2)


if __name__ == "__main__":
    unittest.main()
