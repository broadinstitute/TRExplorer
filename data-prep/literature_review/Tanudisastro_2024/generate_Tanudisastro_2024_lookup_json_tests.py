import unittest

import intervaltree
import pandas as pd

from str_analysis.utils.canonical_repeat_unit import compute_canonical_motif

from generate_Tanudisastro_2024_lookup_json import (
    build_lookup_entry,
    build_tanudisastro_loci,
    find_trexplorer_match_by_start_and_motif,
    parse_table_s1,
    parse_table_s4,
    parse_table_s6,
    pick_tanudisastro_locus_for_trexplorer_locus,
)


def make_table_s1(rows):
    """Returns a raw Table S1 (as read from the xlsx, all text) from (chrom, start, end, motif, cell type, eGene ID,
    eGene symbol, p-value, effect size) tuples."""
    return pd.DataFrame([
        [chrom, str(start), str(end), motif, cell_type, egene_id, egene_symbol, str(pvalue), str(effect_size)]
        for chrom, start, end, motif, cell_type, egene_id, egene_symbol, pvalue, effect_size in rows
    ], columns=["chromosome", "start coordinate (hg38)", "end coordinate (hg38)", "motif", "cell_type", "eGene",
                "eGene_symbol", "nominal p-value", "effect size"])


def make_table_s6(rows):
    """Returns a raw Table S6 from (chrom, start, motif, cell type, eGene ID, eGene symbol, p-value, effect size, PIP)
    tuples."""
    return pd.DataFrame([
        [chrom, str(start), motif, str(effect_size), str(pvalue), egene_id, egene_symbol, cell_type, str(pip)]
        for chrom, start, motif, cell_type, egene_id, egene_symbol, pvalue, effect_size, pip in rows
    ], columns=["chromosome", "start coordinate (hg38)", "motif", "effect size", "nominal p-value", "eGene",
                "eGene_symbol", "cell_type", "pip"])


def make_table_s4(rows):
    """Returns a raw Table S4 from (chrom, start, end, motif, cell type, ATAC p-value) tuples."""
    return pd.DataFrame([
        [chrom, str(start), str(end), motif, cell_type, str(atac_pvalue)]
        for chrom, start, end, motif, cell_type, atac_pvalue in rows
    ], columns=["chromosome", "start coordinate (hg38)", "end coordinate (hg38)", "motif", "cell_type", "atac_pval"])


class BuildTanudisastroLociTests(unittest.TestCase):

    def setUp(self):
        self.s1 = parse_table_s1(make_table_s1([
            ("chr5", 100, 113, "A", "B_naive", "ENSG1", "CAMK4", 2.8e-28, -0.19),
            ("chr5", 100, 113, "A", "CD4_Naive", "ENSG2", "ENSG2", 2.2e-20, 0.14),
        ]))

    def test_table_s1_rows_become_associations(self):
        loci = build_tanudisastro_loci(self.s1, parse_table_s6(make_table_s6([])), parse_table_s4(make_table_s4([])))
        self.assertEqual(list(loci), [("5", 100, 113, "A")])
        self.assertEqual(loci[("5", 100, 113, "A")]["associations"][("B_naive", "ENSG1")],
                         {"cell_type": "B_naive", "egene": "CAMK4", "pvalue": 2.8e-28, "effect_size": -0.19, "pip": None})

    def test_table_s6_adds_pips_new_associations_and_new_loci(self):
        s6 = parse_table_s6(make_table_s6([
            # the same association as a Table S1 row: only its PIP is added
            ("chr5", 100, "A", "CD4_Naive", "ENSG2", "ENSG2", 2.2e-20, 0.14, 0.99),
            # a Table S1 locus, but a cell type that Table S1 doesn't list for it
            ("chr5", 100, "A", "NK", "ENSG2", "ENSG2", 1e-9, 0.1, 0.8),
            # a locus that's not in Table S1
            ("chr19", 500, "TGC", "NK", "ENSG3", "CACNA1A", 1e-8, 0.2, 0.95),
        ]))
        loci = build_tanudisastro_loci(self.s1, s6, parse_table_s4(make_table_s4([])))
        associations = loci[("5", 100, 113, "A")]["associations"]
        self.assertEqual(associations[("CD4_Naive", "ENSG2")]["pip"], 0.99)
        self.assertEqual(associations[("CD4_Naive", "ENSG2")]["pvalue"], 2.2e-20)
        self.assertEqual(associations[("NK", "ENSG2")],
                         {"cell_type": "NK", "egene": "ENSG2", "pvalue": 1e-9, "effect_size": 0.1, "pip": 0.8})
        self.assertIsNone(associations[("B_naive", "ENSG1")]["pip"])
        self.assertEqual(loci[("19", 500, None, "TGC")]["associations"][("NK", "ENSG3")]["egene"], "CACNA1A")

    def test_table_s6_row_goes_to_the_table_s1_loci_that_have_its_association(self):
        s1 = parse_table_s1(make_table_s1([
            ("chr1", 10, 20, "T", "NK", "ENSG1", "G1", 1e-5, 0.1),
            ("chr1", 10, 30, "T", "B_naive", "ENSG2", "G2", 1e-6, 0.1),
            ("chr1", 10, 31, "T", "B_naive", "ENSG2", "G2", 1e-6, 0.1),
        ]))
        s6 = parse_table_s6(make_table_s6([
            ("chr1", 10, "T", "B_naive", "ENSG2", "G2", 1e-6, 0.1, 0.9),
            # an association that none of the three Table S1 loci with this start and motif has
            ("chr1", 10, "T", "CD4_TCM", "ENSG3", "G3", 1e-9, 0.2, 0.8),
        ]))
        loci = build_tanudisastro_loci(s1, s6, parse_table_s4(make_table_s4([])))
        self.assertEqual(loci[("1", 10, 30, "T")]["associations"][("B_naive", "ENSG2")]["pip"], 0.9)
        self.assertEqual(loci[("1", 10, 31, "T")]["associations"][("B_naive", "ENSG2")]["pip"], 0.9)
        self.assertIsNone(loci[("1", 10, 20, "T")]["associations"][("NK", "ENSG1")]["pip"])
        self.assertEqual(loci[("1", 10, 20, "T")]["associations"][("CD4_TCM", "ENSG3")]["pip"], 0.8)
        self.assertNotIn(("CD4_TCM", "ENSG3"), loci[("1", 10, 30, "T")]["associations"])

    def test_table_s4_gives_the_smallest_caqtl_pvalue_per_cell_type(self):
        s4 = parse_table_s4(make_table_s4([
            ("chr5", 100, 113, "A", "B_naive", 4e-4),
            ("chr5", 100, 113, "A", "B_naive", 2.7e-5),
            ("chr5", 100, 113, "A", "CD4_Naive", 8e-4),
        ]))
        loci = build_tanudisastro_loci(self.s1, parse_table_s6(make_table_s6([])), s4)
        self.assertEqual(loci[("5", 100, 113, "A")]["caqtl_min_pvalue_by_cell_type"], {"B_naive": 2.7e-5, "CD4_Naive": 8e-4})

    def test_table_s4_locus_missing_from_table_s1_raises(self):
        s4 = parse_table_s4(make_table_s4([("chr5", 999, 1013, "A", "B_naive", 4e-4)]))
        with self.assertRaises(ValueError):
            build_tanudisastro_loci(self.s1, parse_table_s6(make_table_s6([])), s4)

    def test_duplicate_table_s1_association_raises(self):
        s1 = parse_table_s1(make_table_s1([
            ("chr5", 100, 113, "A", "B_naive", "ENSG1", "CAMK4", 2.8e-28, -0.19),
            ("chr5", 100, 113, "A", "B_naive", "ENSG1", "CAMK4", 3e-28, -0.19),
        ]))
        with self.assertRaises(ValueError):
            build_tanudisastro_loci(s1, parse_table_s6(make_table_s6([])), parse_table_s4(make_table_s4([])))


class BuildLookupEntryTests(unittest.TestCase):

    def setUp(self):
        s1 = parse_table_s1(make_table_s1([
            ("chr5", 100, 113, "A", "B_naive", "ENSG1", "CAMK4", 2.8121469e-28, -0.188319171),
            ("chr5", 100, 113, "A", "CD4_Naive", "ENSG2", "ENSG2", 2.19735e-20, 0.144942036),
            ("chr5", 100, 113, "A", "CD8_Naive", "ENSG2", "ENSG2", 2.218197e-07, 0.0899370158),
            ("chr5", 100, 113, "A", "CD4_Naive", "ENSG1", "CAMK4", 5e-3, 0.01),
        ]))
        s6 = parse_table_s6(make_table_s6([("chr5", 100, "A", "CD4_Naive", "ENSG2", "ENSG2", 2.19735e-20, 0.144942036, 0.991330078)]))
        s4 = parse_table_s4(make_table_s4([
            ("chr5", 100, 113, "A", "CD4_Naive", 8.42e-4),
            ("chr5", 100, 113, "A", "B_naive", 2.69e-5),
        ]))
        self.loci = build_tanudisastro_loci(s1, s6, s4)

    def test_entry(self):
        entry = build_lookup_entry(("5", 100, 113, "A"), self.loci[("5", 100, 113, "A")], "5-100-113-A")
        self.assertEqual(entry, {
            "SignificantCellTypes": "B_naive,CD4_Naive,CD8_Naive",
            "MinPvalue": 2.8121469e-28,
            "FineMappingMaxPIP": 0.991,
            "ChromatinAccessibilityCellTypes": "B_naive,CD4_Naive",
            "Details": {"associations": [
                {"cellType": "B_naive", "eGene": "CAMK4", "pvalue": 2.81e-28, "effectSize": -0.1883},
                {"cellType": "CD4_Naive", "eGene": "ENSG2", "pvalue": 2.2e-20, "effectSize": 0.1449, "fineMappingPIP": 0.991},
                {"cellType": "CD8_Naive", "eGene": "ENSG2", "pvalue": 2.22e-07, "effectSize": 0.0899},
                {"cellType": "CD4_Naive", "eGene": "CAMK4", "pvalue": 0.005, "effectSize": 0.01},
            ]},
        })

    def test_fuzzy_match_keeps_the_tanudisastro_locus_id(self):
        entry = build_lookup_entry(("5", 100, 113, "A"), self.loci[("5", 100, 113, "A")], "5-99-113-A")
        self.assertEqual(entry["Details"]["tanudisastroLocusId"], "5-100-113-A")

    def test_same_interval_with_a_rotated_or_reverse_complement_motif_is_the_same_locus(self):
        locus = {"associations": {("NK", "ENSG1"): {"cell_type": "NK", "egene": "C9orf72", "pvalue": 1e-60,
                                                    "effect_size": 0.06, "pip": None}},
                 "caqtl_min_pvalue_by_cell_type": {}}
        for trexplorer_locus_id in ["9-27573528-27573546-GGCCCC", "9-27573528-27573546-CCGGGG"]:
            entry = build_lookup_entry(("9", 27573528, 27573546, "CCCCGG"), locus, trexplorer_locus_id)
            self.assertNotIn("tanudisastroLocusId", entry["Details"])
        # A different motif length in the same interval is a different locus
        entry = build_lookup_entry(("9", 27573528, 27573546, "CCCCGG"), locus, "9-27573528-27573546-GGC")
        self.assertEqual(entry["Details"]["tanudisastroLocusId"], "9-27573528-27573546-CCCCGG")

    def test_locus_only_in_table_s6_has_no_caqtl_and_a_tanudisastro_locus_id_only_if_its_start_differs(self):
        s6 = parse_table_s6(make_table_s6([("chr19", 500, "TGC", "NK", "ENSG3", "CACNA1A", 1.04e-8, 0.2, 0.987)]))
        loci = build_tanudisastro_loci(parse_table_s1(make_table_s1([])), s6, parse_table_s4(make_table_s4([])))
        entry = build_lookup_entry(("19", 500, None, "TGC"), loci[("19", 500, None, "TGC")], "19-498-530-GCT")
        self.assertEqual(entry["Details"]["tanudisastroLocusId"], "19-500-TGC")
        entry = build_lookup_entry(("19", 500, None, "TGC"), loci[("19", 500, None, "TGC")], "19-500-530-GCT")
        self.assertEqual(entry, {
            "SignificantCellTypes": "NK",
            "MinPvalue": 1.04e-8,
            "FineMappingMaxPIP": 0.987,
            "ChromatinAccessibilityCellTypes": None,
            "Details": {"associations": [
                {"cellType": "NK", "eGene": "CACNA1A", "pvalue": 1.04e-8, "effectSize": 0.2, "fineMappingPIP": 0.987}]},
        })


class PickTanudisastroLocusForTRExplorerLocusTests(unittest.TestCase):

    def test_smallest_pvalue_then_the_same_locus(self):
        s1 = parse_table_s1(make_table_s1([
            ("chr1", 100, 108, "T", "NK", "ENSG1", "G1", 0, 0.1),
            ("chr1", 100, 110, "T", "NK", "ENSG1", "G1", 0, 0.1),
            ("chr1", 100, 110, "T", "B_naive", "ENSG1", "G1", 1e-5, 0.1),
            ("chr1", 99, 110, "T", "NK", "ENSG1", "G1", 1e-3, 0.1),
        ]))
        loci = build_tanudisastro_loci(s1, parse_table_s6(make_table_s6([])), parse_table_s4(make_table_s4([])))
        candidate_keys = [("1", 99, 110, "T"), ("1", 100, 108, "T"), ("1", 100, 110, "T")]
        # Both loci with p = 0 beat the one with p = 1e-3, and the one that is the TRExplorer locus wins the tie,
        # though ('1', 100, 108, 'T') sorts first
        self.assertEqual(pick_tanudisastro_locus_for_trexplorer_locus(candidate_keys, loci, "1-100-110-T"),
                         ("1", 100, 110, "T"))
        self.assertEqual(pick_tanudisastro_locus_for_trexplorer_locus(candidate_keys, loci, "1-100-109-T"),
                         ("1", 100, 108, "T"))


class FindTRExplorerMatchByStartAndMotifTests(unittest.TestCase):

    def test_match(self):
        tree = intervaltree.IntervalTree()
        for locus_id, start, end, motif in [("19-500-530-GCT", 500, 530, "GCT"), ("19-470-560-GCT", 470, 560, "GCT"),
                                            ("19-490-530-A", 490, 530, "A"), ("19-600-660-AGCGGCCGC", 600, 660, "AGCGGCCGC")]:
            tree.add(intervaltree.Interval(start, end, data={
                "locus_id": locus_id, "start_0based": start, "canonical_motif": compute_canonical_motif(motif),
                "motif_length": len(motif)}))
        interval_trees = {"19": tree}
        # TGC and GCT are rotations of the same motif
        self.assertEqual(find_trexplorer_match_by_start_and_motif("19", 500, "TGC", interval_trees),
                         ("19-500-530-GCT", "same start"))
        # A start inside two loci with the same motif: the one with the closer start
        self.assertEqual(find_trexplorer_match_by_start_and_motif("19", 503, "TGC", interval_trees),
                         ("19-500-530-GCT", "overlapping start"))
        # Motifs over 6 bp only need the same length
        self.assertEqual(find_trexplorer_match_by_start_and_motif("19", 603, "GCTGCAGCG", interval_trees),
                         ("19-600-660-AGCGGCCGC", "overlapping start"))
        # A different motif, a start outside every locus with the motif, and a chromosome with no loci
        self.assertEqual(find_trexplorer_match_by_start_and_motif("19", 495, "AT", interval_trees), (None, None))
        self.assertEqual(find_trexplorer_match_by_start_and_motif("19", 580, "TGC", interval_trees), (None, None))
        self.assertEqual(find_trexplorer_match_by_start_and_motif("2", 500, "TGC", interval_trees), (None, None))


if __name__ == "__main__":
    unittest.main()
