import gzip
import json
import tempfile
import unittest
from pathlib import Path

import openpyxl
import pandas as pd

from generate_Garg_2019_and_2020_DMR_lookup_json import (
    AUTISM,
    GARG_2020,
    GARG_2020_COLUMNS,
    SCHIZOPHRENIA,
    build_garg_2019_dmrs,
    build_garg_2020_dmrs,
    build_interval_trees_by_chrom,
    build_lookup,
    does_repeat_contain_cpg_sites,
    get_methylation_change,
    lift_interval_to_hg38,
    read_fragile_site_statuses,
    split_list,
)


def make_fake_lifter(offset=1000, unmapped=(), reversed_strand=()):
    """Returns a stand-in for liftover.get_lifter(..., one_based=False) whose lifter[chrom][pos] adds offset to pos,
    maps the positions in unmapped to nothing, and maps the positions in reversed_strand to the - strand."""
    class ChromLifter:
        def __init__(self, chrom):
            self.chrom = chrom if chrom.startswith("chr") else f"chr{chrom}"

        def __getitem__(self, pos):
            if pos in unmapped:
                return []
            return [(self.chrom, pos + offset, "-" if pos in reversed_strand else "+")]

    class Lifter:
        def __getitem__(self, chrom):
            return ChromLifter(chrom)

    return Lifter()


def write_catalog(path, records):
    with gzip.open(path, "wt") as f:
        json.dump(records, f)


def catalog_record(chrom, start_0based, end_1based, motif):
    return {
        "LocusId": f"{chrom}-{start_0based}-{end_1based}-{motif}",
        "ReferenceRegion": f"chr{chrom}:{start_0based}-{end_1based}",
        "LocusStructure": f"({motif})*",
        "VariantType": "Repeat",
    }


class DoesRepeatContainCpgSitesTests(unittest.TestCase):

    def test_motifs_with_cg_inside_one_copy_or_across_two_copies(self):
        for motif in ["CG", "GC", "CGG", "GCC", "GAC", "GGGGCC", "ccg"]:
            with self.subTest(motif=motif):
                self.assertTrue(does_repeat_contain_cpg_sites(motif))

    def test_motifs_without_cg_even_if_they_contain_gc(self):
        for motif in ["CAG", "CTG", "AGC", "GCA", "C", "G", "GGCA"]:
            with self.subTest(motif=motif):
                self.assertFalse(does_repeat_contain_cpg_sites(motif))


class ParsingTests(unittest.TestCase):

    def test_split_list(self):
        self.assertEqual(split_list("SAMD11,NOC2L"), ["SAMD11", "NOC2L"])
        self.assertEqual(split_list("."), [])
        self.assertEqual(split_list(float("nan")), [])
        self.assertEqual(split_list("8,10"), ["8", "10"])

    def test_get_methylation_change(self):
        self.assertEqual(get_methylation_change("MAX,MAX"), "hypermethylated")
        self.assertEqual(get_methylation_change("MIN"), "hypomethylated")
        with self.assertRaises(ValueError):
            get_methylation_change("MAX,MIN")
        with self.assertRaises(ValueError):
            get_methylation_change("UP")


class LiftIntervalToHg38Tests(unittest.TestCase):

    def test_lifts_the_first_and_last_base(self):
        self.assertEqual(lift_interval_to_hg38(make_fake_lifter(offset=1000), "chr1", 100, 200), ("1", 1100, 1200))

    def test_reversed_interval_is_sorted(self):
        class ReversingLifter:
            def __getitem__(self, chrom):
                return {100: [("chr2", 5000, "-")], 199: [("chr2", 4901, "-")]}
        self.assertEqual(lift_interval_to_hg38(ReversingLifter(), "chr1", 100, 200), ("2", 4901, 5001))

    def test_unmapped_end_or_strand_change_raises(self):
        with self.assertRaises(ValueError):
            lift_interval_to_hg38(make_fake_lifter(unmapped={199}), "chr1", 100, 200)
        with self.assertRaises(ValueError):
            lift_interval_to_hg38(make_fake_lifter(reversed_strand={100}), "chr1", 100, 200)


class BuildDmrsTests(unittest.TestCase):

    def garg_2019_table(self, **extra_columns):
        return pd.DataFrame({
            "chrom": ["chr1", "chr1"],
            "start": [100, 500],
            "end": [200, 600],
            "individuals": ["11393.s1,11393.p1", "12361.s1"],
            "promoter_genes": ["ZMYM6", "."],
            "gene_body_genes": ["ZMYM6", "PRDM16"],
            "num_familial_controls": [0, 1],
            "num_population_controls": [2, 0],
            "num_intellectual_disability_or_congenital_anomaly_patients": [0, 0],
            "direction": ["MAX,MAX", "MIN"],
            "num_probes": [8, 4],
            "individuals_with_failed_probe_detection": ["11393.p1", "."],
            **extra_columns,
        })

    def test_autism_dmrs_only_found_in_siblings_are_left_out(self):
        table = self.garg_2019_table(statuses=["Sibling,Proband", "Sibling"])
        dmrs = build_garg_2019_dmrs(table, AUTISM, make_fake_lifter(offset=1000))
        self.assertEqual(dmrs, [("1", 1100, 1200, {
            "dmrHg38": "chr1:1100-1200",
            "dmrHg19": "chr1:100-200",
            "methylationChange": "hypermethylated",
            "numProbandsWithDMR": 1,
            "numSiblingsWithDMR": 1,
            "numProbes": 8,
            "numCarriersWithFailedProbeDetection": 1,
            "promoterGenes": ["ZMYM6"],
            "geneBodyGenes": ["ZMYM6"],
            "numPopulationControlsWithDMR": 2,
            "numFamilialControlsWithDMR": 0,
            "numIntellectualDisabilityOrCongenitalAnomalyPatientsWithDMR": 0,
        })])

    def test_autism_statuses_that_dont_pair_up_with_individuals_raise(self):
        table = self.garg_2019_table(statuses=["Proband", "Sibling"])
        with self.assertRaises(ValueError):
            build_garg_2019_dmrs(table, AUTISM, make_fake_lifter())

    def test_schizophrenia_dmrs_count_individuals(self):
        dmrs = build_garg_2019_dmrs(self.garg_2019_table(), SCHIZOPHRENIA, make_fake_lifter(offset=0))
        self.assertEqual([(dmr["dmrHg38"], dmr["numIndividualsWithDMR"], dmr["methylationChange"], dmr["promoterGenes"],
                           dmr["numCarriersWithFailedProbeDetection"]) for _, _, _, dmr in dmrs],
                         [("chr1:100-200", 2, "hypermethylated", ["ZMYM6"], 1),
                          ("chr1:500-600", 1, "hypomethylated", [], 0)])
        self.assertNotIn("numProbandsWithDMR", dmrs[0][3])

    def test_garg_2020_has_one_entry_per_tr_with_a_0based_start(self):
        table = pd.DataFrame({
            "chrom": ["chr10"],
            "start": [1000],
            "end": [3000],
            "promoter_genes": ["DPYSL4"],
            "gene_body_genes": ["DPYSL4"],
            "gene_region": ["5' UTR"],
            "cohorts": ["Infant_GSE82273,MESA_GSE56046"],
            "individuals": ["s409,s809,101103"],
            "direction": ["MAX,MAX,MAX"],
            "trs": ["10:1201-1227,10:2524-2597"],
            "num_tr_alleles": ["8,10"],
            "frequency": ["0.03 - 3.43"],
            "fragile_site": ["FRA10X"],
            "validation": ["Validation with PacBio"],
        })
        dmrs = build_garg_2020_dmrs(table, {"FRA10X": "putative novel"}, make_fake_lifter(offset=0))
        self.assertEqual([interval[:3] for interval in dmrs], [("10", 1200, 1227), ("10", 2523, 2597)])
        self.assertEqual(dmrs[1][3], {
            "trHg38": "chr10:2523-2597",
            "trHg19": "10:2524-2597",
            "dmrHg38": "chr10:1000-3000",
            "dmrHg19": "chr10:1000-3000",
            "methylationChange": "hypermethylated",
            "numIndividualsWithDMR": 3,
            "frequencyPer10000": "0.03 - 3.43",
            "numTRAllelesIn600Genomes": 10,
            "promoterGenes": ["DPYSL4"],
            "geneBodyGenes": ["DPYSL4"],
            "geneRegion": "5' UTR",
            "cohorts": ["Infant_GSE82273", "MESA_GSE56046"],
            "fragileSite": "FRA10X",
            "fragileSiteStatus": "putative novel",
            "validation": "Validation with PacBio",
        })

    def test_garg_2020_trs_that_dont_pair_up_with_allele_counts_raise(self):
        table = pd.DataFrame({
            "chrom": ["chr10"], "start": [1000], "end": [3000], "promoter_genes": ["."], "gene_body_genes": ["."],
            "gene_region": ["exon"], "cohorts": ["X"], "individuals": ["a"], "direction": ["MAX"],
            "trs": ["10:1201-1227,10:2524-2597"], "num_tr_alleles": ["8"], "frequency": ["0 - 1.28"],
            "fragile_site": [float("nan")], "validation": [float("nan")],
        })
        with self.assertRaises(ValueError):
            build_garg_2020_dmrs(table, {}, make_fake_lifter())


class ReadFragileSiteStatusesTests(unittest.TestCase):

    def write_table_s15(self, path, fragile_site_cells):
        """Writes a Table S15 with 2 title rows, a header row, and one row per (fragile site, bold, italic)."""
        workbook = openpyxl.Workbook()
        worksheet = workbook.active
        worksheet.append(["Table S15. 25 epivariations overlapping putative unstable CG-rich repeats"])
        worksheet.append([])
        worksheet.append(["Chr, DMR", GARG_2020_COLUMNS["fragile_site"]])
        for row_index, (fragile_site, bold, italic) in enumerate(fragile_site_cells, start=4):
            worksheet.cell(row_index, 1, "chr1")
            cell = worksheet.cell(row_index, 2, fragile_site)
            cell.font = openpyxl.styles.Font(bold=bold, italic=italic)
        workbook.save(path)

    def test_bold_is_previously_defined_and_italic_is_putative_novel(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            path = str(Path(temp_dir) / "mmc16.xlsx")
            self.write_table_s15(path, [("FRA11A", True, False), (None, False, False), ("FRA1M", False, True)])
            self.assertEqual(read_fragile_site_statuses(path),
                             {"FRA11A": "previously defined", "FRA1M": "putative novel"})

    def test_a_fragile_site_that_is_neither_bold_nor_italic_raises(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            path = str(Path(temp_dir) / "mmc16.xlsx")
            self.write_table_s15(path, [("FRA11A", False, False)])
            with self.assertRaises(ValueError):
                read_fragile_site_statuses(path)


class BuildLookupTests(unittest.TestCase):

    def test_keeps_loci_that_overlap_by_at_least_1bp_and_whose_repeat_contains_cpg_sites(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            catalog_path = str(Path(temp_dir) / "catalog.json.gz")
            write_catalog(catalog_path, [
                catalog_record("1", 120, 150, "CGG"),   # overlaps the autism and both schizophrenia DMRs
                catalog_record("1", 120, 150, "CAG"),   # overlaps them too, but has no CpG sites
                catalog_record("1", 90, 100, "GCC"),    # ends where the autism DMR starts
                catalog_record("1", 300, 330, "GCC"),   # overlaps nothing
                catalog_record("2", 95, 105, "GC"),     # overlaps the Garg 2020 TR by 1 bp
                catalog_record("3", 95, 105, "GC"),     # on a chromosome with no DMRs
            ])
            autism_dmr = {"dmrHg38": "chr1:100-200"}
            schizophrenia_dmrs = [{"dmrHg38": "chr1:140-260"}, {"dmrHg38": "chr1:110-125"}]
            garg_2020_dmr = {"trHg38": "chr2:104-130"}
            lookup = build_lookup(catalog_path, {
                AUTISM: build_interval_trees_by_chrom([("1", 100, 200, autism_dmr)]),
                SCHIZOPHRENIA: build_interval_trees_by_chrom([
                    ("1", 140, 260, schizophrenia_dmrs[0]), ("1", 110, 125, schizophrenia_dmrs[1])]),
                GARG_2020: build_interval_trees_by_chrom([("2", 104, 130, garg_2020_dmr)]),
            })
        self.assertEqual(lookup, {
            # the schizophrenia DMRs are listed by start, whatever their order in the input
            "1-120-150-CGG": {AUTISM: [autism_dmr], SCHIZOPHRENIA: [schizophrenia_dmrs[1], schizophrenia_dmrs[0]]},
            "2-95-105-GC": {GARG_2020: [garg_2020_dmr]},
        })


if __name__ == "__main__":
    unittest.main()
