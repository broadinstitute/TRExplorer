"""Compare the simulated allele sizes with the pathogenic thresholds of the STRchive disease loci.

For each STRchive locus with a pathogenic threshold, finds the TRExplorer locus it overlaps (same motif
length, most overlapping bp) in HPRC256_total_allele_length_stats.tsv.gz, computes the target lengths
of make_tr_expansion_benchmark_alleles.py (the 2.5th, 97.5th and 99.5th percentile HPRC256 lengths, and
the 99.5th plus multiples of the motif's range), and prints them as a markdown table in
total repeat units on STRchive's scale (STRchive's reference count plus the change in units from the
hg38 tract, so that differences between STRchive and TRExplorer locus boundaries cancel out), next to
the threshold and the smallest shown size that reaches it.

Loci marked * are not comparable by length: the disease is caused by a different motif inserted into
the repeat, or STRchive's threshold is at or below the reference count.

Usage:
    python3 compare_simulated_sizes_to_pathogenic_thresholds.py
"""
import collections
import csv
import gzip
import json
import os
import sys

# The pipeline modules are in the parent folder
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from make_tr_expansion_benchmark_alleles import (  # noqa: E402
    HIGH_PERCENTILE, HIGHEST_PERCENTILE, HPRC256_ALLELE_LENGTH_STATS_TSV, LOW_PERCENTILE,
    MIN_LOCI_FOR_A_MOTIF_SPECIFIC_RANGE, MOTIF_RANGE_PERCENTILES_TSV, SPLICEAI_ANNOTATION, canonical_motif,
    is_inside_a_transcript, load_merged_transcript_spans, load_motif_ranges, parse_allele_length_histogram,
    TARGET_LABELS, is_polymorphic, labeled_size_changes_in_repeat_units, percentile_of_histogram, short_target_label,
    target_lengths)

SHOWN_TARGETS = {label: short_target_label(label) for label in TARGET_LABELS}

STRCHIVE_LOCI_JSON = "/Users/weisburd/code/STRchive/data/STRchive-loci.json"
NOT_COMPARABLE_BY_LENGTH = {"RAI1", "MIR7-2", "VWA1", "DAB1", "SAMD12", "BEAN1", "RFC1", "STARD7", "YEATS2",
                            "TNRC6A", "MARCHF6", "RAPGEF2", "XYLT1"}


def match_strchive_loci_to_hprc256(strchive_loci):
    """Returns {STRchive id: (locus_id, start, end, motif, {length: count})} for the best-overlapping HPRC256 locus."""
    by_chrom = collections.defaultdict(list)
    for locus in strchive_loci:
        by_chrom[locus["chrom"].replace("chr", "")].append(locus)
    best = {}
    with gzip.open(HPRC256_ALLELE_LENGTH_STATS_TSV, "rt") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            chrom, start, end, motif = row["locus_id"].rsplit("-", 3)
            start, end = int(start), int(end)
            for locus in by_chrom.get(chrom, []):
                overlap = min(end, int(locus["stop_hg38"])) - max(start, int(locus["start_hg38"]))
                if len(motif) == int(locus["motif_len"]) and overlap > best.get(locus["id"], (0,))[0]:
                    best[locus["id"]] = (overlap, row["locus_id"], start, end, motif,
                                         parse_allele_length_histogram(row["allele_length_histogram"]))
    return {strchive_id: match[1:] for strchive_id, match in best.items()}


def main():
    strchive_loci = sorted((x for x in json.load(open(STRCHIVE_LOCI_JSON)) if x.get("pathogenic_min")), key=lambda x: x["gene"])
    matches = match_strchive_loci_to_hprc256(strchive_loci)
    motif_range_bp = load_motif_ranges()
    n_loci_by_motif = {row["group"]: int(row["n_loci"]) for row in csv.DictReader(open(MOTIF_RANGE_PERCENTILES_TSV), delimiter="\t")
                       if row["group_type"] == "canonical_motif"}
    merged_spans = load_merged_transcript_spans(SPLICEAI_ANNOTATION["basic"])
    genes = collections.Counter(x["gene"] for x in strchive_loci)

    header = ["Gene", "Motif", "Motif range (units)", "Ref"] + list(SHOWN_TARGETS.values()) + [
        "Pathogenic", "Smallest size reaching it", "In run"]
    print("| " + " | ".join(header) + " |\n|" + "---|" * len(header))
    seen, tally = collections.Counter(), collections.Counter()
    for x in strchive_loci:
        seen[x["gene"]] += 1
        label = (f"{x['gene']} ({seen[x['gene']]})" if genes[x["gene"]] > 1 else x["gene"]) + (" *" if x["gene"] in NOT_COMPARABLE_BY_LENGTH else "")
        pathogenic = float(x["pathogenic_min"])
        if x["id"] not in matches:
            print(f"| {label} | {x['motif_len']} bp | | no HPRC256 match |" + " |" * len(SHOWN_TARGETS) + f" {pathogenic:.0f} | | |")
            continue
        locus_id, start, end, motif, counts_by_length = matches[x["id"]]
        low_bp, high_bp, highest_bp = (percentile_of_histogram(counts_by_length, p) for p in (LOW_PERCENTILE, HIGH_PERCENTILE, HIGHEST_PERCENTILE))
        reference_bp, motif_size = end - start, len(motif)
        motif_range = motif_range_bp(motif) / motif_size  # in repeat units, for the table
        pooled = n_loci_by_motif.get(canonical_motif(motif), 0) < MIN_LOCI_FOR_A_MOTIF_SPECIFIC_RANGE
        lengths_bp = target_lengths(low_bp, high_bp, highest_bp, motif_range_bp(motif))
        reference_units = float(x["ref_copies"]) if x.get("ref_copies") else reference_bp / motif_size
        # The simulated alleles' changes in whole units, as the allele design rounds and limits them
        changes = labeled_size_changes_in_repeat_units(reference_bp, motif_size, lengths_bp)
        change_by_label = {label: change for change, labels in changes.items() for label in labels}

        polymorphic = is_polymorphic(counts_by_length) and bool(changes)
        if not polymorphic:
            in_run = "not polymorphic"
        elif is_inside_a_transcript(merged_spans, "chr" + locus_id.split("-")[0], start):
            in_run = "yes"
        else:
            in_run = "outside basic transcripts"
        units_by_target = {short: reference_units + change_by_label.get(label, 0) for label, short in SHOWN_TARGETS.items()}
        smallest = next((short for short, units in units_by_target.items() if units >= pathogenic), "none")
        motif_label = canonical_motif(motif) if motif_size <= 12 else f"{motif_size} bp"
        print(f"| {label} | {motif_label} | {motif_range:.1f}{' (pooled)' if pooled else ''} | {reference_units:.0f} | "
              + " | ".join(f"{units:.0f}" for units in units_by_target.values())
              + f" | {pathogenic:.0f} | {smallest} | {in_run} |")
        if x["gene"] not in NOT_COMPARABLE_BY_LENGTH:
            tally["comparable"] += 1
            tally[("smallest", smallest)] += 1
            tally["in run"] += in_run == "yes"
            tally[("smallest in run", smallest)] += in_run == "yes"
    print(f"\n{tally['comparable']} loci comparable by length ({tally['in run']} of them in the run); smallest size that "
          f"reaches the pathogenic threshold:")
    for short in list(SHOWN_TARGETS.values()) + ["none"]:
        print(f"  {short:<7} {tally[('smallest', short)]:>3} loci ({tally[('smallest in run', short)]} in the run)")


if __name__ == "__main__":
    main()
