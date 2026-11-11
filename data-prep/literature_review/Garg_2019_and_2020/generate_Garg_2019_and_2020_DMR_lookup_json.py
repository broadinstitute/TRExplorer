"""Match the DMRs from Garg et al. 2019 and Garg et al. 2020 to TRExplorer loci and create a lookup JSON for TRExplorer.

A DMR (differentially methylated region) here is a rare epivariation: a region where one individual's DNA methylation
in blood is far outside the range seen in controls. Both papers called DMRs from Illumina 450k methylation array data
(hg19), requiring at least 3 probes within 1 kb that are each above the 99.9th (hypermethylated) or below the 0.1th
(hypomethylated) percentile of the controls.

Garg et al. 2019, "Screening for rare epigenetic variations in autism and schizophrenia" (Human Mutation,
doi:10.1002/humu.23740, PMC6801017), studied 297 autism probands and 286 of their unaffected siblings, and 767
individuals with schizophrenia:
  - NIHMS1054049-supplement-Supp_TableS1.xlsx: 150 DMRs found in the autism probands and their siblings
  - NIHMS1054049-supplement-Supp_TableS2.xlsx: 309 DMRs found in the schizophrenia samples

Garg et al. 2020, "A survey of rare epigenetic variation in 23,116 human genomes identifies disease-relevant
epivariations and novel CGG expansions" (American Journal of Human Genetics, doi:10.1016/j.ajhg.2020.08.019,
PMC7536611), found hypermethylated DMRs in 23,116 individuals and intersected them with 180 TRs that they predicted to be
unstable (motif of 3-6 bp, 100% GC, at least 10 different alleles in 600 genomes genotyped by HipSTR):
  - mmc16.xlsx (Table S15): the 25 DMRs within 500 bp of one of these TRs, and the TR (two TRs for one DMR)

Matching to TRExplorer loci (each locus has to overlap the hg38 interval by at least 1 bp):
  - Garg 2019: a locus is matched to a DMR if it overlaps the DMR and its repeat contains CpG sites. Autism DMRs are
    only used if at least one autism proband had them (105 of the 150), since the others were only found in
    unaffected siblings.
  - Garg 2020: a locus is matched to a DMR if it overlaps the DMR's TR named in Table S15 (not the DMR itself, since 6
    of the 25 TRs lie outside their DMR) and its repeat contains CpG sites.
A repeat contains CpG sites if its motif written twice contains CG: GCC (GCCGCC) and CG do, but CAG and CTG (CAGCAG,
CTGCTG) don't, even though their motifs contain GC. Methylation arrays measure CpG sites, so only these TRs can be
methylated themselves.

The hg19 intervals are lifted over to hg38 with the UCSC hg19ToHg38 chain file (the liftover package downloads it). The
DMR coordinates in the three tables are BED-style (0-based start at the first probe's C, end after the last probe's G),
while the Table S15 TR coordinates are HipSTR's (1-based start, inclusive end). Both were checked against the hg19
reference sequence.

Output is a JSON lookup table keyed by TRExplorer locus ID, with an entry for every locus that matched at least one DMR.
Each entry has one or more of these keys, each a list of the matched DMRs (usually one):
  - Garg2019_AutismDMRs: dmrHg38 (0-based start), dmrHg19 (as listed in the paper), methylationChange
    ("hypermethylated" or "hypomethylated"), numProbandsWithDMR, numSiblingsWithDMR, numProbes,
    numCarriersWithFailedProbeDetection (the number of carriers in whom a cluster of the DMR's probes failed detection,
    which can mean a technical problem or a deletion rather than a methylation change), promoterGenes and
    geneBodyGenes (genes whose promoter, +/-2 kb of the TSS, or gene body overlaps the DMR),
    numPopulationControlsWithDMR (of 2,711), numFamilialControlsWithDMR (McRae et al. 2014 families), and
    numIntellectualDisabilityOrCongenitalAnomalyPatientsWithDMR (of 489)
  - Garg2019_SchizophreniaDMRs: the same as Garg2019_AutismDMRs, but with numIndividualsWithDMR instead of
    numProbandsWithDMR and numSiblingsWithDMR
  - Garg2020_DMRs: trHg38 (0-based start), trHg19 (as listed in the paper), dmrHg38, dmrHg19, methylationChange,
    numIndividualsWithDMR (of 23,116), frequencyPer10000 (95% confidence interval, as listed in the paper),
    numTRAllelesIn600Genomes, promoterGenes, geneBodyGenes, geneRegion (e.g. "5' UTR"), cohorts, fragileSite (null if
    none), fragileSiteStatus ("previously defined" for fragile sites already known to be caused by an expansion of the
    TR, or "putative novel" for those the authors propose it causes; the table shows this only through the cell
    formatting), and validation (the paper's notes on long-read validation or the reference that first described the
    expansion; null if none)
"""

import argparse
import collections
import gzip
import json
from pathlib import Path
import re
import sys

import intervaltree
import openpyxl
import pandas as pd
from liftover import get_lifter

from str_analysis.utils.eh_catalog_utils import get_variant_catalog_iterator

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from trexplorer_catalog_matching import parse_catalog_record  # noqa: E402

AUTISM = "Garg2019_AutismDMRs"
SCHIZOPHRENIA = "Garg2019_SchizophreniaDMRs"
GARG_2020 = "Garg2020_DMRs"

METHYLATION_CHANGE_BY_DIRECTION = {"MAX": "hypermethylated", "MIN": "hypomethylated"}

# The columns read from each table. The tables have 2 title rows above the header row.
GARG_2019_COLUMNS = {
    "chrom": "Chr_hg19",
    "start": "Start_hg19",
    "end": "End_hg19",
    "individuals": "Individual with DMR",
    "promoter_genes": "Gene promoter (+/-2kb of TSS) overlap",
    "gene_body_genes": "Gene body overlap",
    "num_familial_controls": "Observed in McRae familial controls?",
    "num_population_controls": "Observed in 2711 population controls?",
    "num_intellectual_disability_or_congenital_anomaly_patients": "Observed in 489 patients with ID/CA?",
    "direction": "Direction of methylation change",
    "num_probes": "Number of probes in DMR",
    "individuals_with_failed_probe_detection": "Cluster of probes with failed detection p-values?",
}
GARG_2019_AUTISM_STATUS_COLUMN = "Individual status"
GARG_2020_COLUMNS = {
    "chrom": "Chr, DMR",
    "start": "Start, DMR",
    "end": "End, DMR",
    "promoter_genes": "Gene Promoter (+/-2kb of TSS)",
    "gene_body_genes": "Gene Body overlap",
    "gene_region": "Location",
    "cohorts": "Cohort",
    "individuals": "Individual IDs with DMR",
    "direction": "Direction",
    "trs": "Putative unstable CG-rich tandem repeat (+/-500bp of DMR)",
    "num_tr_alleles": "Number of TR alleles observed in 600 WGS samples genotyped by hipSTR",
    "frequency": "Frequency of DMR per 10,000 (95% CI)",
    "fragile_site": "Fragile site (Previously defined in bold, putative novel in red italic)",
    "validation": "Reference or validation data with long-read sequencing",
}


def parse_args():
    parser = argparse.ArgumentParser(description="Generate the Garg 2019 and Garg 2020 DMR lookup JSON")
    parser.add_argument("--garg-2019-autism-table", default="NIHMS1054049-supplement-Supp_TableS1.xlsx",
                        help="Garg et al. 2019 Supplementary Table 1 (autism DMRs)")
    parser.add_argument("--garg-2019-schizophrenia-table", default="NIHMS1054049-supplement-Supp_TableS2.xlsx",
                        help="Garg et al. 2019 Supplementary Table 2 (schizophrenia DMRs)")
    parser.add_argument("--garg-2020-table", default="mmc16.xlsx",
                        help="Garg et al. 2020 Table S15 (DMRs near putatively unstable CGG repeats)")
    parser.add_argument("--trexplorer-catalog", required=True, help="Path to TRExplorer v2.1 JSON catalog")
    parser.add_argument("-o", "--output", default="Garg_2019_and_2020_DMR_lookup.json.gz",
                        help="Output JSON path (default: Garg_2019_and_2020_DMR_lookup.json.gz)")
    return parser.parse_args()


def does_repeat_contain_cpg_sites(motif):
    """Returns True if a repeat of this motif contains CpG sites, including ones that span two copies of the motif.

    Args:
        motif (str): repeat motif, eg. "GCC"

    Returns:
        bool: True for motifs like GCC (GCCGCC) or CG, False for motifs like CAG (CAGCAG) that contain GC but no CG
    """
    return "CG" in (motif + motif).upper()


def read_table(path, columns):
    """Reads one of the papers' supplementary tables, which have 2 title rows above the header row.

    Args:
        path (str): path of the .xlsx table
        columns (dict): {short name: column name in the table} of the columns to read

    Returns:
        DataFrame: the table, with its columns renamed to the short names

    Raises:
        ValueError: if the table is missing one of the columns
    """
    table = pd.read_excel(path, header=2)
    missing_columns = set(columns.values()) - set(table.columns)
    if missing_columns:
        raise ValueError(f"{path} is missing the columns {sorted(missing_columns)}. Found: {list(table.columns)}")
    return table[list(columns.values())].rename(columns={name: short for short, name in columns.items()})


def read_fragile_site_statuses(path):
    """Returns the status of each fragile site in Table S15, which the table shows only through the cell formatting.

    Args:
        path (str): path of Table S15 (mmc16.xlsx)

    Returns:
        dict: {fragile site: "previously defined" (in bold in the table) or "putative novel" (in red italic)}

    Raises:
        ValueError: if a fragile site is neither bold nor italic, or both
    """
    worksheet = openpyxl.load_workbook(path).worksheets[0]
    # The table has 2 title rows above the header row
    header = [cell.value for cell in worksheet[3]]
    fragile_site_column_index = header.index(GARG_2020_COLUMNS["fragile_site"])
    statuses = {}
    for row in worksheet.iter_rows(min_row=4):
        cell = row[fragile_site_column_index]
        if cell.value is None:
            continue
        if bool(cell.font.b) == bool(cell.font.i):
            raise ValueError(f"Fragile site {cell.value} in {path} is {'both' if cell.font.b else 'neither'} bold "
                             f"and italic, so its status is unclear")
        statuses[str(cell.value).strip()] = "previously defined" if cell.font.b else "putative novel"
    return statuses


def split_list(value):
    """Splits a comma-separated list from the tables, where "." or an empty cell means an empty list."""
    if pd.isna(value) or str(value).strip() in ("", "."):
        return []
    return [item.strip() for item in str(value).split(",")]


def get_methylation_change(directions):
    """Returns "hypermethylated" or "hypomethylated" for a DMR's comma-separated MAX/MIN directions (one per individual).

    Raises:
        ValueError: if the directions are not all the same, or one is neither MAX nor MIN
    """
    directions = set(split_list(directions))
    if len(directions) != 1 or not directions <= set(METHYLATION_CHANGE_BY_DIRECTION):
        raise ValueError(f"Expected all DMR directions to be MAX, or all to be MIN. Found: {sorted(directions)}")
    return METHYLATION_CHANGE_BY_DIRECTION[directions.pop()]


def lift_interval_to_hg38(lifter, chrom, start_0based, end_1based):
    """Lifts an hg19 interval over to hg38 by lifting its first and last base.

    Prints the interval when its length changes, since that means the two builds differ inside it.

    Args:
        lifter: hg19 to hg38 converter returned by liftover.get_lifter(..., one_based=False)
        chrom (str): hg19 chromosome, with or without "chr"
        start_0based (int): 0-based start
        end_1based (int): 1-based end

    Returns:
        tuple: (chrom without "chr", start_0based, end_1based) in hg38

    Raises:
        ValueError: if either end doesn't lift over to a single position, or the two ends land on different
            chromosomes or strands
    """
    first_base = lifter[chrom][start_0based]
    last_base = lifter[chrom][end_1based - 1]
    hg19_interval = format_interval(chrom.replace("chr", ""), start_0based, end_1based)
    if len(first_base) != 1 or len(last_base) != 1:
        raise ValueError(f"{hg19_interval} (hg19) didn't lift over to hg38")
    (first_base_chrom, first_base_pos, first_base_strand), = first_base
    (last_base_chrom, last_base_pos, last_base_strand), = last_base
    if first_base_chrom != last_base_chrom or first_base_strand != last_base_strand:
        raise ValueError(f"The ends of {hg19_interval} (hg19) lifted over to different chromosomes or strands in hg38")
    lifted_start_0based, lifted_last_base_0based = sorted([first_base_pos, last_base_pos])
    lifted_interval = first_base_chrom.replace("chr", ""), lifted_start_0based, lifted_last_base_0based + 1
    if lifted_interval[2] - lifted_interval[1] != end_1based - start_0based:
        print(f"  {hg19_interval} is {end_1based - start_0based:,d} bp in hg19 and "
              f"{lifted_interval[2] - lifted_interval[1]:,d} bp in hg38 ({format_interval(*lifted_interval)})")
    return lifted_interval


def format_interval(chrom, start_0based, end_1based):
    """Returns the interval as chr<chrom>:<start_0based>-<end_1based>, the way TRExplorer writes intervals."""
    return f"chr{chrom}:{start_0based}-{end_1based}"


def build_garg_2019_dmrs(table, cohort, lifter):
    """Returns the Garg 2019 DMRs to match to TRExplorer loci, lifted over to hg38.

    Args:
        table (DataFrame): Supplementary Table 1 or 2, as returned by read_table
        cohort (str): AUTISM or SCHIZOPHRENIA. For AUTISM, the table also needs a "statuses" column (Proband or Sibling
            for each individual), and DMRs that no proband had are left out.
        lifter: hg19 to hg38 converter returned by liftover.get_lifter(..., one_based=False)

    Returns:
        list: (hg38 chrom, start_0based, end_1based, DMR dict) tuples, with the DMR dict as described in the module
            docstring

    Raises:
        ValueError: if a DMR doesn't lift over to hg38, or an autism DMR's individuals and statuses don't pair up
    """
    dmrs = []
    for row in table.itertuples(index=False):
        dmr_hg19 = f"{row.chrom}:{row.start}-{row.end}"
        individuals = split_list(row.individuals)
        if cohort == AUTISM:
            statuses = split_list(row.statuses)
            if len(statuses) != len(individuals) or not set(statuses) <= {"Proband", "Sibling"}:
                raise ValueError(f"Unexpected individual statuses {statuses} for individuals {individuals} with the "
                                 f"{dmr_hg19} DMR")
            if "Proband" not in statuses:
                continue
            carriers = {"numProbandsWithDMR": statuses.count("Proband"), "numSiblingsWithDMR": statuses.count("Sibling")}
        else:
            carriers = {"numIndividualsWithDMR": len(individuals)}

        lifted_interval = lift_interval_to_hg38(lifter, row.chrom, int(row.start), int(row.end))
        dmr = {
            "dmrHg38": format_interval(*lifted_interval),
            "dmrHg19": dmr_hg19,
            "methylationChange": get_methylation_change(row.direction),
            **carriers,
            "numProbes": int(row.num_probes),
            "numCarriersWithFailedProbeDetection": len(split_list(row.individuals_with_failed_probe_detection)),
            "promoterGenes": split_list(row.promoter_genes),
            "geneBodyGenes": split_list(row.gene_body_genes),
            "numPopulationControlsWithDMR": int(row.num_population_controls),
            "numFamilialControlsWithDMR": int(row.num_familial_controls),
            "numIntellectualDisabilityOrCongenitalAnomalyPatientsWithDMR": int(
                row.num_intellectual_disability_or_congenital_anomaly_patients),
        }
        dmrs.append((*lifted_interval, dmr))
    return dmrs


def build_garg_2020_dmrs(table, fragile_site_statuses, lifter):
    """Returns the Garg 2020 TRs to match to TRExplorer loci, lifted over to hg38, with their DMRs.

    Args:
        table (DataFrame): Table S15, as returned by read_table
        fragile_site_statuses (dict): {fragile site: status}, as returned by read_fragile_site_statuses
        lifter: hg19 to hg38 converter returned by liftover.get_lifter(..., one_based=False)

    Returns:
        list: (hg38 chrom, start_0based, end_1based, DMR dict) tuples, one per TR, with the TR's interval and the DMR
            dict as described in the module docstring

    Raises:
        ValueError: if a DMR or TR doesn't lift over to hg38, or a DMR's TRs and TR allele counts don't pair up
    """
    dmrs = []
    for row in table.itertuples(index=False):
        dmr_hg19 = f"{row.chrom}:{row.start}-{row.end}"
        lifted_dmr_interval = lift_interval_to_hg38(lifter, row.chrom, int(row.start), int(row.end))
        fragile_site = None if pd.isna(row.fragile_site) else str(row.fragile_site).strip()
        trs = split_list(row.trs)
        num_tr_alleles = split_list(row.num_tr_alleles)
        if not trs or len(trs) != len(num_tr_alleles):
            raise ValueError(f"The {dmr_hg19} DMR's TRs {trs} don't pair up with its TR allele counts {num_tr_alleles}")
        for tr_hg19, tr_num_alleles in zip(trs, num_tr_alleles):
            # HipSTR intervals have a 1-based start, eg. 8:104311002-104311038
            tr_chrom, tr_start_1based, tr_end_1based = re.fullmatch(r"(\w+):(\d+)-(\d+)", tr_hg19).groups()
            lifted_tr_interval = lift_interval_to_hg38(
                lifter, f"chr{tr_chrom}", int(tr_start_1based) - 1, int(tr_end_1based))
            dmrs.append((*lifted_tr_interval, {
                "trHg38": format_interval(*lifted_tr_interval),
                "trHg19": tr_hg19,
                "dmrHg38": format_interval(*lifted_dmr_interval),
                "dmrHg19": dmr_hg19,
                "methylationChange": get_methylation_change(row.direction),
                "numIndividualsWithDMR": len(split_list(row.individuals)),
                "frequencyPer10000": str(row.frequency).strip(),
                "numTRAllelesIn600Genomes": int(tr_num_alleles),
                "promoterGenes": split_list(row.promoter_genes),
                "geneBodyGenes": split_list(row.gene_body_genes),
                "geneRegion": str(row.gene_region).strip(),
                "cohorts": split_list(row.cohorts),
                "fragileSite": fragile_site,
                "fragileSiteStatus": fragile_site_statuses[fragile_site] if fragile_site else None,
                "validation": None if pd.isna(row.validation) else str(row.validation).strip(),
            }))
    return dmrs


def build_interval_trees_by_chrom(intervals):
    """Returns {chrom: IntervalTree} of (chrom, start_0based, end_1based, DMR dict) tuples, with the DMR dicts as data."""
    interval_trees = collections.defaultdict(intervaltree.IntervalTree)
    for chrom, start_0based, end_1based, dmr in intervals:
        interval_trees[chrom].addi(start_0based, end_1based, dmr)
    return interval_trees


def build_lookup(catalog_path, interval_trees_by_dataset):
    """Matches each TRExplorer locus whose repeat contains CpG sites to the DMRs (or Garg 2020 TRs) it overlaps.

    Args:
        catalog_path (str): path of the TRExplorer JSON catalog
        interval_trees_by_dataset (dict): {AUTISM, SCHIZOPHRENIA or GARG_2020: interval trees returned by
            build_interval_trees_by_chrom}

    Returns:
        dict: {LocusId: {dataset: list of DMR dicts}}, with a dataset only if the locus matched at least one of its DMRs
    """
    lookup = {}
    for record in get_variant_catalog_iterator(catalog_path, show_progress_bar=True):
        parsed_record = parse_catalog_record(record)
        if parsed_record is None:
            continue
        locus_id, chrom, start_0based, end_1based, motif = parsed_record

        matched_dmrs_by_dataset = {}
        for dataset, interval_trees in interval_trees_by_dataset.items():
            # Sorted so that a locus's DMRs are listed in the same order on every run
            overlapping_intervals = sorted(interval_trees[chrom].overlap(start_0based, end_1based),
                                           key=lambda interval: (interval.begin, interval.end))
            if overlapping_intervals:
                matched_dmrs_by_dataset[dataset] = [interval.data for interval in overlapping_intervals]
        if not matched_dmrs_by_dataset:
            continue

        if not does_repeat_contain_cpg_sites(motif):
            continue
        if locus_id in lookup:
            raise ValueError(f"LocusId {locus_id} appears more than once in {catalog_path}")
        lookup[locus_id] = matched_dmrs_by_dataset
    return lookup


def main():
    args = parse_args()
    lifter = get_lifter("hg19", "hg38", one_based=False)

    autism_table = read_table(args.garg_2019_autism_table,
                              {**GARG_2019_COLUMNS, "statuses": GARG_2019_AUTISM_STATUS_COLUMN})
    schizophrenia_table = read_table(args.garg_2019_schizophrenia_table, GARG_2019_COLUMNS)
    garg_2020_table = read_table(args.garg_2020_table, GARG_2020_COLUMNS)

    print("Lifting the DMRs and TRs over to hg38")
    dmrs_by_dataset = {
        AUTISM: build_garg_2019_dmrs(autism_table, AUTISM, lifter),
        SCHIZOPHRENIA: build_garg_2019_dmrs(schizophrenia_table, SCHIZOPHRENIA, lifter),
        GARG_2020: build_garg_2020_dmrs(garg_2020_table, read_fragile_site_statuses(args.garg_2020_table), lifter),
    }
    print(f"Read {len(autism_table):,d} autism DMRs from {args.garg_2019_autism_table}, and kept the "
          f"{len(dmrs_by_dataset[AUTISM]):,d} that at least one proband had")
    print(f"Read {len(dmrs_by_dataset[SCHIZOPHRENIA]):,d} schizophrenia DMRs from {args.garg_2019_schizophrenia_table}")
    print(f"Read {len(dmrs_by_dataset[GARG_2020]):,d} TRs near {len(garg_2020_table):,d} DMRs from {args.garg_2020_table}")

    print(f"Reading TRExplorer catalog: {args.trexplorer_catalog}")
    lookup = build_lookup(args.trexplorer_catalog, {
        dataset: build_interval_trees_by_chrom(dmrs) for dataset, dmrs in dmrs_by_dataset.items()})

    print(f"\nOutput: {len(lookup):,d} TRExplorer loci")
    for dataset, dmrs in dmrs_by_dataset.items():
        # A Garg 2020 TR is identified by its hg19 interval, and a Garg 2019 DMR by its own
        hg19_interval_key = "trHg19" if dataset == GARG_2020 else "dmrHg19"
        matched_dmrs = [dmr for entry in lookup.values() for dmr in entry.get(dataset, [])]
        unmatched_hg19_intervals = sorted(
            {dmr[hg19_interval_key] for _, _, _, dmr in dmrs} - {dmr[hg19_interval_key] for dmr in matched_dmrs})
        print(f"  {dataset}: {sum(dataset in entry for entry in lookup.values()):,d} loci matched "
              f"{len(dmrs) - len(unmatched_hg19_intervals):,d} of {len(dmrs):,d} "
              f"{'TRs' if dataset == GARG_2020 else 'DMRs'}")
        print(f"    loci by methylation change: "
              f"{dict(collections.Counter(dmr['methylationChange'] for dmr in matched_dmrs))}")
        if dataset == GARG_2020:
            for tr_hg19 in unmatched_hg19_intervals:
                print(f"    WARNING: no TRExplorer locus with CpG sites overlaps the Garg 2020 TR {tr_hg19}")

    fopen = gzip.open if args.output.endswith("gz") else open
    with fopen(args.output, "wt") as f:
        json.dump(lookup, f, indent=2)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
