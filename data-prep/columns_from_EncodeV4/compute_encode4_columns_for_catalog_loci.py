"""Computes the ENCODE4_* columns for every TRExplorer catalog locus from the downloaded ENCODE4 files.

Reads the manifest written by download_encode4_source_files.py and the catalog BED, and writes one row
per locus that has a value in at least one column (loci without any value are left out, the same way
the SpliceAI summary TSV leaves out unscored loci):

  LocusId: chrom (without "chr") - 0-based start - end - motif, e.g. 1-10000-10108-TAACCC
  ENCODE4_and_IGVF_CRISPR_or_MPRA_validated_elements: data source:assay:biosample type triples in which
      an element called functional overlaps the locus, e.g. "ENCODE:MPRA:HepG2,IGVF:CRISPR:HCT116". The
      calls are: ENCODE CRISPR rows with Significant=TRUE, ENCODE MPRA elements with log2 fold change > 0
      and FDR < 0.05, and IGVF CRISPR elements called by the rule published for each file layout (see
      parse_igvf_crispr_significant_elements). An empty value does not mean the locus was tested.

"Overlap" means at least 1 bp in common, with all intervals treated as 0-based half-open.

To keep memory low, the catalog is held as numpy arrays and the ENCODE files are read one at a time.
The peak memory use is printed after each step.

Usage:
    python3 compute_encode4_columns_for_catalog_loci.py \\
        --catalog-bed ../splicing_prediction/TRExplorer.repeat_catalog_v2.1.hg38.1_to_1000bp_motifs.bed.gz
"""

import argparse
import collections
import gzip
import math
import resource
import sys
from pathlib import Path

import numpy as np
import pandas as pd

OUTPUT_COLUMNS = [
    "LocusId",
    "ENCODE4_and_IGVF_CRISPR_or_MPRA_validated_elements",
]

# Loci up to this length are looked up with a fixed search window; the few longer ones are handled
# separately so they don't widen the window for everything else.
SHORT_LOCUS_MAX_LENGTH = 1000

MPRA_MIN_LOG2_FOLD_CHANGE = 0
MPRA_MAX_FDR = 0.05
MPRA_ELEMENT_ENRICHMENTS_COLUMNS = [
    "chrom", "start", "end", "name", "score", "strand", "log2FoldChange", "inputCount", "outputCount",
    "minusLog10PValue", "minusLog10QValue",
]

# FDR cutoff for IGVF CRISPR files without a significant/not column. Both the Engreitz lab Flow-FISH and
# TAP-seq analyses (Guckelberger et al., bioRxiv 2024) and FRACTEL (Gersbach lab) call elements at 0.05.
IGVF_CRISPR_MAX_FDR = 0.05


def print_peak_memory(step_description):
    """Prints the process's peak resident memory so far."""
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    peak_bytes = peak if sys.platform == "darwin" else peak * 1024
    print(f"  Peak memory after {step_description}: {peak_bytes / 1e9:.2f} GB")


def build_catalog_index(catalog_bed):
    """Reads the catalog BED into per-chromosome numpy arrays for overlap lookups.

    Args:
        catalog_bed: path of the catalog BED (chrom, 0-based start, end, motif, ...), optionally gzipped.

    Returns:
        (n_loci, index) where index maps chrom to {"short": (starts, ends, locus_indices),
        "long": (starts, ends, locus_indices)}. Each triple is sorted by start. Locus indices are the
        loci's 0-based row numbers in the catalog BED. "short" loci are at most SHORT_LOCUS_MAX_LENGTH bp.
    """
    catalog = pd.read_csv(catalog_bed, sep="\t", header=None, usecols=[0, 1, 2], names=["chrom", "start", "end"],
                          dtype={"chrom": "category", "start": np.int64, "end": np.int64})
    catalog["locus_index"] = np.arange(len(catalog), dtype=np.int64)

    index = {}
    for chrom, loci in catalog.groupby("chrom", observed=True, sort=False):
        loci = loci.sort_values("start", kind="stable")
        is_short = (loci["end"] - loci["start"]).to_numpy() <= SHORT_LOCUS_MAX_LENGTH
        index[chrom] = {
            subset_name: (loci["start"].to_numpy()[mask], loci["end"].to_numpy()[mask],
                          loci["locus_index"].to_numpy()[mask])
            for subset_name, mask in [("short", is_short), ("long", ~is_short)]
        }
    return len(catalog), index


def find_overlapping_pairs(sorted_starts, ends, max_length, query_starts, query_ends):
    """Returns all (query, interval) pairs that share at least 1 bp.

    Intervals must be sorted by start and be at most max_length bp long. All intervals are 0-based
    half-open. An interval overlaps a query when interval_start < query_end and interval_end > query_start.
    Because interval_end <= interval_start + max_length, only intervals with
    interval_start > query_start - max_length can overlap, which bounds the search window.

    Args:
        sorted_starts: numpy array of interval starts, sorted.
        ends: numpy array of interval ends, in the same order as sorted_starts.
        max_length: the longest interval's length.
        query_starts: numpy array of query starts.
        query_ends: numpy array of query ends.

    Returns:
        (query_positions, interval_positions): numpy arrays of matching positions in the input arrays.
    """
    window_starts = np.searchsorted(sorted_starts, query_starts - max_length, side="right")
    window_ends = np.searchsorted(sorted_starts, query_ends, side="left")
    window_sizes = np.maximum(window_ends - window_starts, 0)

    query_positions = np.repeat(np.arange(len(query_starts)), window_sizes)
    offsets_within_window = np.arange(window_sizes.sum()) - np.repeat(np.cumsum(window_sizes) - window_sizes, window_sizes)
    interval_positions = np.repeat(window_starts, window_sizes) + offsets_within_window

    overlaps = ends[interval_positions] > query_starts[query_positions]
    return query_positions[overlaps], interval_positions[overlaps]


def find_loci_overlapping_elements(catalog_index, elements):
    """Returns the catalog loci that overlap each element.

    Args:
        catalog_index: index returned by build_catalog_index.
        elements: DataFrame with chrom, start, end columns (0-based half-open), any row order.

    Returns:
        (element_row_positions, locus_indices): numpy arrays of matching pairs, where
        element_row_positions are 0-based positions in the elements DataFrame.
    """
    all_element_row_positions = []
    all_locus_indices = []
    elements = elements.reset_index(drop=True)
    for chrom, chrom_elements in elements.groupby("chrom", sort=False):
        if chrom not in catalog_index:
            continue
        element_row_positions = chrom_elements.index.to_numpy()
        element_starts = chrom_elements["start"].to_numpy(dtype=np.int64)
        element_ends = chrom_elements["end"].to_numpy(dtype=np.int64)

        # Short loci: search the loci around each element.
        locus_starts, locus_ends, locus_indices = catalog_index[chrom]["short"]
        if len(locus_starts):
            element_positions, locus_positions = find_overlapping_pairs(
                locus_starts, locus_ends, SHORT_LOCUS_MAX_LENGTH, element_starts, element_ends)
            all_element_row_positions.append(element_row_positions[element_positions])
            all_locus_indices.append(locus_indices[locus_positions])

        # Long loci: search the elements around each locus.
        locus_starts, locus_ends, locus_indices = catalog_index[chrom]["long"]
        if len(locus_starts) and len(element_starts):
            order = np.argsort(element_starts, kind="stable")
            max_element_length = int((element_ends - element_starts).max())
            locus_positions, sorted_element_positions = find_overlapping_pairs(
                element_starts[order], element_ends[order], max_element_length, locus_starts, locus_ends)
            all_element_row_positions.append(element_row_positions[order[sorted_element_positions]])
            all_locus_indices.append(locus_indices[locus_positions])

    if not all_locus_indices:
        return np.array([], dtype=np.int64), np.array([], dtype=np.int64)
    return np.concatenate(all_element_row_positions), np.concatenate(all_locus_indices)


def parse_crispr_significant_elements_tsv(path):
    """Returns the elements with Significant=TRUE (chrom, start, end) from a CRISPR element quantifications TSV.

    Returns None if the file has no Significant column.
    """
    elements = pd.read_csv(path, sep="\t", dtype={"Significant": str})
    if "Significant" not in elements.columns:
        return None
    elements = elements[elements["Significant"].str.upper() == "TRUE"]
    return elements.rename(columns={"chromStart": "start", "chromEnd": "end"})[["chrom", "start", "end"]]


def parse_mpra_active_elements_bed(path):
    """Returns the active elements (chrom, start, end) from an MPRA "element enrichments" BED (BED6+5).

    Active means log2 fold change > MPRA_MIN_LOG2_FOLD_CHANGE and FDR < MPRA_MAX_FDR. Returns None if the
    file doesn't have the 11 columns of the format, or if it has no FDR values: some ENCODE files put the
    placeholder -1 (or 0) in minusLog10QValue on every row, which would otherwise silently yield no elements.
    """
    elements = pd.read_csv(path, sep="\t", header=None, comment="#")
    if len(elements.columns) != len(MPRA_ELEMENT_ENRICHMENTS_COLUMNS):
        return None
    elements.columns = MPRA_ELEMENT_ENRICHMENTS_COLUMNS
    if (elements["minusLog10QValue"] <= 0).all():
        return None
    is_active = ((elements["log2FoldChange"] > MPRA_MIN_LOG2_FOLD_CHANGE)
                 & (elements["minusLog10QValue"] > -math.log10(MPRA_MAX_FDR)))
    return elements.loc[is_active, ["chrom", "start", "end"]]


def parse_interval_strings(interval_strings):
    """Returns (chrom, start, end) columns parsed from strings like "chr1:100-200" or "chr1:100-200_93".

    Strings without coordinates, such as the negative controls ("nc_MYC", "negative_control") and gene-level
    targets ("NA:NA-NA_ABL1") in some IGVF files, are dropped.
    """
    parts = interval_strings.str.extract(r"^(chr[^:]+):(\d+)-(\d+)").dropna()
    return pd.DataFrame({"chrom": parts[0], "start": parts[1].astype(int), "end": parts[2].astype(int)})


def parse_igvf_crispr_significant_elements(path, file_type):
    """Returns the significant elements (chrom, start, end) from an IGVF CRISPR differential element file.

    IGVF CRISPR files come in many layouts. Each layout is recognized by its columns, and its elements
    are called with the rule that the producing lab published:
      - a Significant or significant column (Engreitz lab Flow-FISH and TAP-seq): rows marked TRUE
      - growth_significant or mig_significant (Gersbach lab proliferation and migration screens): TRUE
      - FRACTEL_pval_fdr_corr (Gersbach lab Perturb-seq and FACS screens): FDR < IGVF_CRISPR_MAX_FDR
      - adj.pval.EnhancerEffect.noAux (Engreitz lab cohesin Flow-FISH screens, Guckelberger et al.):
        FDR < IGVF_CRISPR_MAX_FDR in the condition with cohesin present
      - minus_auxin_padj with element_hg38 (Engreitz lab cohesin TAP-seq): FDR < IGVF_CRISPR_MAX_FDR
    Other layouts have no published rule, or no hg38 coordinates, and return None.

    Args:
        path: path of the file.
        file_type: "csv" or "tsv".

    Returns:
        DataFrame with chrom, start, end columns, or None if the file's layout isn't one of the above.
    """
    table = pd.read_csv(path, sep="," if file_type == "csv" else "\t", encoding="utf-8-sig", low_memory=False)
    columns = set(table.columns)

    def is_true(column):
        return table[column].astype(str).str.upper() == "TRUE"

    def fdr_below_cutoff(column):
        return pd.to_numeric(table[column], errors="coerce") < IGVF_CRISPR_MAX_FDR

    flag_column = next((c for c in ["Significant", "significant", "growth_significant", "mig_significant"] if c in columns), None)
    if flag_column:
        is_significant = is_true(flag_column)
    elif "FRACTEL_pval_fdr_corr" in columns:
        is_significant = fdr_below_cutoff("FRACTEL_pval_fdr_corr")
    elif "adj.pval.EnhancerEffect.noAux" in columns:
        is_significant = fdr_below_cutoff("adj.pval.EnhancerEffect.noAux")
    elif {"minus_auxin_padj", "element_hg38"} <= columns:
        is_significant = fdr_below_cutoff("minus_auxin_padj")
    else:
        return None

    significant = table[is_significant]
    for interval_column in ["name_hg38", "element_hg38", "dhs_coords"]:
        if interval_column in columns:
            return parse_interval_strings(significant[interval_column].astype(str))
    for prefix in ["targeting", "intended_target"]:
        if f"{prefix}_chr" in columns:
            return pd.DataFrame({
                "chrom": significant[f"{prefix}_chr"],
                "start": significant[f"{prefix}_start"].astype(int),
                "end": significant[f"{prefix}_end"].astype(int),
            })
    return None


def compute_validated_element_labels(catalog_index, validated_rows):
    """Returns {locus_index: "data source:assay:biosample type,..."} for loci that overlap elements called functional."""
    encode_parsers = {
        "CRISPR": parse_crispr_significant_elements_tsv,
        "MPRA": parse_mpra_active_elements_bed,
    }
    labels_by_locus = collections.defaultdict(set)
    for row in validated_rows:
        if row["data_source"] == "IGVF":
            elements = parse_igvf_crispr_significant_elements(row["local_path"], row["file_type"])
        else:
            elements = encode_parsers[row["assay_label"]](row["local_path"])
        if elements is None:
            print(f"WARNING: skipping {row['local_path']}: its columns don't match a {row['data_source']} "
                  f"{row['assay_label']} format with a known rule for calling elements")
            continue
        label = f"{row['data_source']}:{row['assay_label']}:{row['biosample_term_name'].replace(',', ';')}"
        for locus_index in np.unique(find_loci_overlapping_elements(catalog_index, elements)[1]):
            labels_by_locus[int(locus_index)].add(label)

    return {locus_index: ",".join(sorted(labels)) for locus_index, labels in labels_by_locus.items()}


def read_manifest(manifest_path):
    """Returns the manifest rows grouped by source: {source: [row dict, ...]}."""
    manifest = pd.read_csv(manifest_path, sep="\t", dtype=str, keep_default_na=False)
    return {source: rows.to_dict("records") for source, rows in manifest.groupby("source")}


def write_output_tsv(catalog_bed, output_tsv, validated_labels_by_locus):
    """Streams the catalog BED and writes one row per locus that has a value in any ENCODE4 column.

    Returns:
        Dict mapping each ENCODE4 column name to the number of loci with a value in it.
    """
    n_loci_with_value = collections.Counter()
    fopen = gzip.open if str(catalog_bed).endswith("gz") else open
    with fopen(catalog_bed, "rt") as catalog, gzip.open(output_tsv, "wt") as out:
        out.write("\t".join(OUTPUT_COLUMNS) + "\n")
        for locus_index, line in enumerate(catalog):
            values = [
                validated_labels_by_locus.get(locus_index, ""),
            ]
            if not any(values):
                continue
            chrom, start, end, motif = line.rstrip("\n").split("\t")[:4]
            locus_id = f"{chrom.removeprefix('chr')}-{start}-{end}-{motif}"
            out.write("\t".join([locus_id] + values) + "\n")
            for column, value in zip(OUTPUT_COLUMNS[1:], values):
                n_loci_with_value[column] += bool(value)
    return n_loci_with_value


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--catalog-bed", type=Path, required=True,
                        help="TRExplorer catalog BED (chrom, 0-based start, end, motif), e.g. "
                             "TRExplorer.repeat_catalog_v2.1.hg38.1_to_1000bp_motifs.bed.gz")
    parser.add_argument("--manifest", type=Path, default=Path("encode4_source_files_manifest.tsv"),
                        help="Manifest written by download_encode4_source_files.py")
    parser.add_argument("--output-tsv", type=Path, default=Path("TRExplorer_v2.1_ENCODE4_columns.tsv.gz"))
    args = parser.parse_args()

    rows_by_source = read_manifest(args.manifest)
    if "validated" not in rows_by_source:
        raise ValueError(f"{args.manifest} has no rows for: validated")

    print(f"Reading {args.catalog_bed}")
    n_loci, catalog_index = build_catalog_index(args.catalog_bed)
    print(f"  {n_loci:,d} loci")
    print_peak_memory("reading the catalog")

    print(f"Computing validated elements from {len(rows_by_source['validated']):,d} files")
    validated_labels_by_locus = compute_validated_element_labels(catalog_index, rows_by_source["validated"])
    print_peak_memory("validated elements")

    n_loci_with_value = write_output_tsv(args.catalog_bed, args.output_tsv, validated_labels_by_locus)
    print(f"Wrote {args.output_tsv}")
    for column in OUTPUT_COLUMNS[1:]:
        print(f"  {n_loci_with_value[column]:>10,d} of {n_loci:,d} loci "
              f"({n_loci_with_value[column] / n_loci:.1%}) have {column}")


if __name__ == "__main__":
    main()
