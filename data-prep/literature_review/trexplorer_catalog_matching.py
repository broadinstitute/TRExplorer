"""Helpers shared by the scripts in this directory that match loci from published studies to TRExplorer catalog loci.

The scripts are in subdirectories, so they add this directory to sys.path before importing this module.
"""

import collections
import re

import intervaltree

from str_analysis.utils.canonical_repeat_unit import compute_canonical_motif
from str_analysis.utils.eh_catalog_utils import get_variant_catalog_iterator
from str_analysis.utils.misc_utils import parse_interval

# The motif of the first repeat in a LocusStructure, e.g. CAG in (CAG)*
LOCUS_STRUCTURE_MOTIF_REGEX = re.compile(r"^[(]([A-Z]+)[)][+*]")


def parse_catalog_record(record):
    """Returns the coordinates and motif of a TRExplorer catalog record.

    Args:
        record (dict): a record of the TRExplorer JSON catalog

    Returns:
        tuple: (locus ID, chrom without "chr", start_0based, end_1based, motif) from the record's first ReferenceRegion
            and its ReferenceMotif, or the motif of the first repeat in its LocusStructure. None if it has neither (a
            complex locus structure).
    """
    reference_region = record["ReferenceRegion"]
    if isinstance(reference_region, list):
        reference_region = reference_region[0]
    chrom, start_0based, end_1based = parse_interval(reference_region)
    if "ReferenceMotif" in record:
        motif = record["ReferenceMotif"]
    else:
        match = LOCUS_STRUCTURE_MOTIF_REGEX.match(record["LocusStructure"])
        if not match:
            return None
        motif = match.group(1)
    return record["LocusId"], chrom.replace("chr", ""), start_0based, end_1based, motif


def build_catalog_interval_trees(catalog_path, only_loci_that_overlap=None, show_progress_bar=True):
    """Reads the TRExplorer catalog loci into interval trees.

    Args:
        catalog_path (str): path of the TRExplorer JSON catalog
        only_loci_that_overlap (iterable): optional (chrom without "chr", start_0based, end_1based) intervals of a
            study's loci. If given, only the catalog loci that overlap one of them are kept, so that memory stays small.
        show_progress_bar (bool): whether to show a progress bar while reading the catalog

    Returns:
        dict: chrom without "chr" -> IntervalTree of the kept loci that parse_catalog_record can parse, each with a dict
            of locus_id, start_0based, end_1based, canonical_motif, motif_length and catalog_index (the locus's position
            in the catalog, which find_exact_or_best_jaccard_match uses to break ties) as its data
    """
    study_interval_trees = None
    if only_loci_that_overlap is not None:
        study_interval_trees = collections.defaultdict(intervaltree.IntervalTree)
        for chrom, start_0based, end_1based in only_loci_that_overlap:
            study_interval_trees[chrom].addi(start_0based, end_1based)

    interval_trees = collections.defaultdict(intervaltree.IntervalTree)
    for catalog_index, record in enumerate(
            get_variant_catalog_iterator(catalog_path, show_progress_bar=show_progress_bar)):
        parsed_record = parse_catalog_record(record)
        if parsed_record is None:
            continue
        locus_id, chrom, start_0based, end_1based, motif = parsed_record
        # IntervalTree.overlap bisects, while IntervalTree.overlaps(begin, end) can scan the chromosome's whole tree
        if study_interval_trees is not None and (
                chrom not in study_interval_trees or not study_interval_trees[chrom].overlap(start_0based, end_1based)):
            continue
        interval_trees[chrom].addi(start_0based, end_1based, {
            "locus_id": locus_id,
            "start_0based": start_0based,
            "end_1based": end_1based,
            "canonical_motif": compute_canonical_motif(motif),
            "motif_length": len(motif),
            "catalog_index": catalog_index,
        })
    print(f"Kept {sum(len(interval_tree) for interval_tree in interval_trees.values()):,d} TRExplorer loci from "
          f"{catalog_path}")
    return interval_trees


def is_matching_motif(motif, canonical_motif, catalog_locus_data):
    """Returns whether a TRExplorer locus has the same motif as a study's locus: the same canonical motif for motifs up
    to 6 bp, and the same motif length for longer ones.

    Args:
        motif (str): the study locus's motif
        canonical_motif (str): compute_canonical_motif(motif)
        catalog_locus_data (dict): the data of the TRExplorer locus's interval in build_catalog_interval_trees
    """
    if len(motif) <= 6:
        return canonical_motif == catalog_locus_data["canonical_motif"]
    return len(motif) == catalog_locus_data["motif_length"]


def find_exact_or_best_jaccard_match(interval_trees, chrom, start_0based, end_1based, motif, min_jaccard_similarity):
    """Finds the TRExplorer locus with the study locus's ID (chrom-start_0based-end_1based-motif), or else the one that
    overlaps the study locus, has the same motif (is_matching_motif) and has the highest Jaccard similarity
    (intersection over union of the two intervals), if that is at least min_jaccard_similarity. Ties go to the locus
    that comes first in the catalog.

    Args:
        interval_trees (dict): returned by build_catalog_interval_trees
        chrom (str): chromosome without "chr"
        start_0based (int): the study locus's start
        end_1based (int): the study locus's end
        motif (str): the study locus's motif
        min_jaccard_similarity (float): the lowest Jaccard similarity that counts as a match

    Returns:
        tuple: (locus ID, match type), or (None, None) if no locus matches. The match type is "exact" for the locus with
            the same ID, or for a Jaccard match with the same interval and canonical motif (the same locus with its
            motif written starting at a different base or on the other strand), and "fuzzy" otherwise.
    """
    if chrom not in interval_trees:
        return None, None
    overlapping_intervals = sorted(interval_trees[chrom].overlap(start_0based, end_1based),
                                   key=lambda interval: interval.data["catalog_index"])
    # A TRExplorer locus ID is chrom-start_0based-end_1based-motif, so a locus with the study locus's ID has its interval
    study_locus_id = f"{chrom}-{start_0based}-{end_1based}-{motif}"
    if any(interval.data["locus_id"] == study_locus_id for interval in overlapping_intervals):
        return study_locus_id, "exact"

    canonical_motif = compute_canonical_motif(motif)
    best_match = None
    best_jaccard = 0
    for interval in overlapping_intervals:
        if not is_matching_motif(motif, canonical_motif, interval.data):
            continue
        intersection = min(end_1based, interval.end) - max(start_0based, interval.begin)
        union = (end_1based - start_0based) + (interval.end - interval.begin) - intersection
        jaccard = intersection / union
        if jaccard >= min_jaccard_similarity and jaccard > best_jaccard:
            best_jaccard = jaccard
            best_match = interval.data
    if best_match is None:
        return None, None
    is_same_locus = (best_match["start_0based"], best_match["end_1based"], best_match["canonical_motif"]) == (
        start_0based, end_1based, canonical_motif)
    return best_match["locus_id"], "exact" if is_same_locus else "fuzzy"
