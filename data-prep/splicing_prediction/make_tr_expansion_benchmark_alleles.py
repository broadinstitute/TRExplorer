"""Make a benchmark set of population-based contraction and expansion alleles at random polymorphic TRExplorer loci.

Samples loci from TRExplorer v2.1 that lie inside a GENCODE v50 transcript of the chosen gene set
(basic or comprehensive), the annotation SpliceAI-lookup scores with (loci outside every transcript
get no scores), optionally limited to a distance range from the nearest GENCODE v50 splice site, and
that are polymorphic in HPRC256: at least MIN_DISTINCT_ALLELE_SIZES_TO_BE_POLYMORPHIC (2) distinct total
allele sizes among the HPRC256 haplotypes (from the TRGT VCF that genotyped all repeats and variation
clusters, using each locus's own single-locus record; the variation-cluster records, whose lengths span
the whole cluster, are skipped, see compute_hprc256_total_allele_length_stats.py), with at least one of the lengths below a whole repeat
unit or more from the reference. The target lengths use percentiles rather than the standard deviation,
which a single extreme allele, such as a mobile element insertion inside the repeat, can inflate a
hundredfold.

Each locus gets up to 6 alleles (TARGET_LABELS), with total lengths:
- the 2.5th, 97.5th and 99.5th percentiles of the HPRC256 lengths (variation seen in people), and
- the 99.5th percentile plus 1x, 2x and 3x its motif's range: the 75th percentile, across every
  HPRC256 locus with the same motif and at least 5 distinct allele lengths, of the 2.5th-97.5th
  percentile range (expansions as large as the motif reaches somewhere in the genome; see
  compute_hprc256_range_percentiles_by_motif.py). A motif with fewer than
  MIN_LOCI_FOR_A_MOTIF_SPECIFIC_RANGE loci uses the value for all motifs of its length. A locus's own
  range says little about how large its rare expansions get (it is unrelated to the distance from
  the 99.5th percentile to the pathogenic threshold at the STRchive disease loci), while +3x the
  motif's range reaches the pathogenic threshold at 53 of 63 of them
  (compare_simulated_sizes_to_pathogenic_thresholds.py).
each written as a whole number of repeat units added to or removed from the hg38 tract:
(target length - reference length) / motif size units, rounded with exact halves away from zero. Lengths
that round to the same change are scored once, and a change of 0 (the reference itself) is dropped.
Contractions keep at least one whole unit and delete at most 5000 bp, and expansions insert at most
5000 bp. Loci whose anchor base or tract contains a base other than A, C, G or T get no alleles, since
the production scorer only accepts ACGT alleles.

Expansions insert whole copies of the tract's most common in-frame unit right before the tract, and
contractions delete the first whole units of the tract; both are anchored on the base before the
tract, as in the earlier SpliceAI analysis of simulated expansions at known disease loci.

Writes benchmark/tr_benchmark_alleles.<gene_set>.json: the alleles grouped by locus, plus the number of loci
the full run would score, which spliceai_l4_cost_and_disk_benchmark.py scales its estimates to.
"""
import argparse
import bisect
import collections
import csv
import gzip
import json
import math
import os
import random

HERE = os.path.dirname(os.path.abspath(__file__))
# Distance from every TRExplorer v2.1 locus to the nearest GENCODE v50 splice site, from
# count_catalog_loci_near_gencode_splice_sites.py (see README.md)
DISTANCES_TSV = os.path.join(HERE, "TRExplorer_v2.1_loci_vs_GENCODE_v50_splice_sites.distances.tsv.gz")
DISTANCE_COLUMN = {
    "basic": "distance_to_nearest_splice_site_bp__gencode.v50.basic.annotation",
    "comprehensive": "distance_to_nearest_splice_site_bp__gencode.v50.annotation",
}
# The GENCODE v50 SpliceAI annotation files the benchmark and the full run score with
SPLICEAI_ANNOTATION = {
    "basic": os.path.join(HERE, "spliceai_lookup_files_gencode_v50", "gencode.v50.basic.annotation.txt.gz"),
    "comprehensive": os.path.join(HERE, "spliceai_lookup_files_gencode_v50", "gencode.v50.annotation.txt.gz"),
}
HPRC256_ALLELE_LENGTH_STATS_TSV = os.path.join(HERE, "HPRC256_total_allele_length_stats.tsv.gz")
MOTIF_RANGE_PERCENTILES_TSV = os.path.join(HERE, "HPRC256_range_percentile_by_motif.tsv")
# Every locus: the farthest from a GENCODE v50 basic splice site is about 2 Mb away.
NO_DISTANCE_LIMIT_BP = 10 ** 9
LOW_PERCENTILE, HIGH_PERCENTILE, HIGHEST_PERCENTILE = 2.5, 97.5, 99.5
MIN_DISTINCT_ALLELE_SIZES_TO_BE_POLYMORPHIC = 2
# A motif's range: this percentile, across its loci, of the LOW-HIGH percentile range
# (75th of loci with 5+ distinct lengths: +3x then reaches the pathogenic threshold at 53 of the 63 STRchive
# disease loci comparable by length, per compare_simulated_sizes_to_pathogenic_thresholds.py; the 99th
# made most +3x alleles 5 kb insertions)
MOTIF_RANGE_PERCENTILE = 75
MIN_LOCI_FOR_A_MOTIF_SPECIFIC_RANGE = 100
# Unseen expansions: the HIGHEST_PERCENTILE length plus these multiples of the motif's range
# (+0.5x was benchmarked too, and found no splice-altering locus that the others missed)
MOTIF_RANGE_MULTIPLES = (1, 2, 3)
# Names of the target lengths, in the order target_lengths returns them
TARGET_LABELS = (f"{LOW_PERCENTILE}th percentile", f"{HIGH_PERCENTILE}th percentile", f"{HIGHEST_PERCENTILE}th percentile") + tuple(
    f"{HIGHEST_PERCENTILE}th + {k}x motif range" for k in MOTIF_RANGE_MULTIPLES)
MAX_INSERTED_BP = 5000
# The production scorer rejects a deletion that reaches past the 10,000 bases scored on either side of
# the variant, so contractions are capped well within that, the same as expansions.
MAX_DELETED_BP = 5000
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def most_common_in_frame_unit(tract, unit_length):
    """Returns the most common unit_length-mer read in frame from the start of the tract.

    Ties go to the unit seen first, so the motif is in the tract's own reading frame.
    """
    units = [tract[i:i + unit_length] for i in range(0, len(tract) - unit_length + 1, unit_length)]
    return collections.Counter(units).most_common(1)[0][0]


def percentile_of_histogram(counts_by_length, percentile):
    """Returns the nearest-rank percentile of the lengths in {length: number of haplotypes}.

    That is the smallest length such that at least percentile% of haplotypes are at most that long.
    """
    rank = max(1, math.ceil(percentile / 100 * sum(counts_by_length.values())))
    n_at_most = 0
    for length in sorted(counts_by_length):
        n_at_most += counts_by_length[length]
        if n_at_most >= rank:
            return length


def parse_allele_length_histogram(histogram):
    """Parses "length:count,..." (compute_hprc256_total_allele_length_stats.py) into {length: count}."""
    return {int(length): int(count) for length, count in (pair.split(":") for pair in histogram.split(","))}


def canonical_motif(motif):
    """Returns the alphabetically first rotation of the motif or its reverse complement, e.g. CAG -> AGC."""
    reverse_complement = motif.translate(COMPLEMENT)[::-1]
    return min(sequence[i:] + sequence[:i] for sequence in (motif, reverse_complement) for i in range(len(sequence)))


def load_motif_ranges():
    """Returns a function that gives a motif's range in whole bp (see MOTIF_RANGE_PERCENTILES_TSV).

    A motif seen at fewer than MIN_LOCI_FOR_A_MOTIF_SPECIFIC_RANGE loci gets the value for all motifs of its
    length, or, for a long motif length with no such loci, the value of the nearest length that has one.
    """
    by_motif, by_length = {}, {}
    with open(MOTIF_RANGE_PERCENTILES_TSV) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            value = int(row["range_percentile_bp"])
            if row["group_type"] == "motif_length":
                by_length[int(row["group"])] = value
            elif int(row["n_loci"]) >= MIN_LOCI_FOR_A_MOTIF_SPECIFIC_RANGE:
                by_motif[row["group"]] = value

    def motif_range_bp(motif):
        if canonical_motif(motif) in by_motif:
            return by_motif[canonical_motif(motif)]
        return by_length[min(by_length, key=lambda length: (abs(length - len(motif)), length))]
    return motif_range_bp


def short_target_label(label):
    """Returns a target label's short form for tables, e.g. "2.5th percentile" -> "2.5th", "99.5th + 1x motif range" -> "+1x"."""
    return label.split(" + ")[1].split()[0].join(["+", ""]) if " + " in label else label.split()[0]


def target_sort_key(label):
    """Sorts target labels as TARGET_LABELS does: percentiles, then increasing multiples of the motif range."""
    return (1, float(short_target_label(label).strip("+x"))) if " + " in label else (0, float(label.split("th")[0]))


def target_lengths(low_bp, high_bp, highest_bp, motif_range_bp):
    """Returns the total allele lengths to simulate at a locus, in the order of TARGET_LABELS.

    Args:
        low_bp, high_bp, highest_bp (int): the locus's LOW, HIGH and HIGHEST percentile HPRC256 lengths
        motif_range_bp (int): its motif's range (see load_motif_ranges), in whole bp, so the lengths are
            whole bp too and labeled_size_changes_in_repeat_units rounds them exactly
    """
    return [low_bp, high_bp, highest_bp] + [highest_bp + k * motif_range_bp for k in MOTIF_RANGE_MULTIPLES]


def labeled_size_changes_in_repeat_units(reference_bp, motif_size, lengths_bp, labels=TARGET_LABELS):
    """Returns {change in whole repeat units: [labels of the lengths it stands for]} for the lengths, sorted by change.

    Each change is (length - reference_bp) / motif_size rounded to a whole number, with exact halves
    rounded away from zero (Python's round() would send them to the even neighbour, merging lengths
    one unit apart, e.g. +1.5 and +2.5 units). The lengths must be whole bp: the rounding is done in
    integer arithmetic, so no floating-point error can move a half-unit change. A contraction that would leave less than one whole unit,
    or delete more than MAX_DELETED_BP, is limited to the largest that does not, and an expansion that
    would insert more than MAX_INSERTED_BP likewise. Lengths that give the same change share one entry,
    and a change of 0 (the reference itself) is left out.
    """
    most_units_removable = min(reference_bp // motif_size - 1, MAX_DELETED_BP // motif_size)
    most_units_insertable = MAX_INSERTED_BP // motif_size
    labels_by_change = collections.defaultdict(list)
    for length, label in zip(lengths_bp, labels):
        difference_bp = length - reference_bp
        # floor(|difference| / motif_size + 1/2), in integers, with the difference's sign
        rounded = (2 * abs(difference_bp) + motif_size) // (2 * motif_size) * (1 if difference_bp >= 0 else -1)
        change = max(-most_units_removable, min(rounded, most_units_insertable))
        if change != 0:
            labels_by_change[change].append(label)
    return dict(sorted(labels_by_change.items()))


def allele_to_variant(chrom, start_0based, anchor, tract, unit, reference_count, repeat_count):
    """Returns the "chrom-pos-ref-alt" string for one simulated allele, anchored on the base before the tract.

    An expansion inserts (repeat_count - reference_count) copies of unit after the anchor; a
    contraction deletes the tract's first (reference_count - repeat_count) whole units.
    """
    if repeat_count > reference_count:
        return f"{chrom}-{start_0based}-{anchor}-{anchor}{unit * (repeat_count - reference_count)}"
    deleted = tract[:len(unit) * (reference_count - repeat_count)]
    return f"{chrom}-{start_0based}-{anchor}{deleted}-{anchor}"


def load_merged_transcript_spans(spliceai_annotation_path):
    """Returns {chrom: (starts, ends)}: the union of transcript spans, 1-based closed, sorted and non-overlapping.

    Uses the SpliceAI annotation's TX_START (0-based) and TX_END columns the way the SpliceAI
    Annotator does: a position is scored when TX_START + 1 <= pos <= TX_END for some transcript.
    """
    spans_by_chrom = collections.defaultdict(list)
    with gzip.open(spliceai_annotation_path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        chrom_i, start_i, end_i = header.index("CHROM"), header.index("TX_START"), header.index("TX_END")
        for line in f:
            fields = line.split("\t")
            spans_by_chrom[fields[chrom_i]].append((int(fields[start_i]) + 1, int(fields[end_i])))
    merged = {}
    for chrom, spans in spans_by_chrom.items():
        starts, ends = [], []
        for start, end in sorted(spans):
            if starts and start <= ends[-1] + 1:
                ends[-1] = max(ends[-1], end)
            else:
                starts.append(start)
                ends.append(end)
        merged[chrom] = (starts, ends)
    return merged


def is_inside_a_transcript(merged_spans, chrom, pos):
    """Returns True if the 1-based pos lies inside some transcript span from load_merged_transcript_spans."""
    starts, ends = merged_spans.get(chrom, ([], []))
    i = bisect.bisect_right(starts, pos) - 1
    return i >= 0 and pos <= ends[i]


def load_loci_near_splice_sites(gene_set, max_distance_bp=NO_DISTANCE_LIMIT_BP, min_distance_bp=0):
    """Returns the TRExplorer loci in a distance range from a splice site of gene_set, and those SpliceAI would score.

    Args:
        gene_set (str): "basic" or "comprehensive"
        max_distance_bp (int): keep loci at most this far from the nearest GENCODE v50 splice site
        min_distance_bp (int): and at least this far (e.g. 501 with max 600 selects the 501-600 bp band)

    Returns:
        tuple: (set of (chrom, start_0based, end, motif) loci in the distance range,
            sorted list of those inside a GENCODE v50 transcript of the same gene set)
    """
    merged_spans = load_merged_transcript_spans(SPLICEAI_ANNOTATION[gene_set])
    # The distances table repeats a few loci; a set counts each locus once.
    loci_near_splice_sites = set()
    with gzip.open(DISTANCES_TSV, "rt") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            distance = row[DISTANCE_COLUMN[gene_set]]
            if distance != "." and min_distance_bp <= int(distance) <= max_distance_bp:
                loci_near_splice_sites.add((row["chrom"], int(row["start_0based"]), int(row["end"]), row["motif"]))
    # The variant is anchored on the base before the tract, whose 1-based position is start_0based.
    scored_loci = sorted(locus for locus in loci_near_splice_sites if is_inside_a_transcript(merged_spans, locus[0], locus[1]))
    print(f"{len(loci_near_splice_sites):,} loci {min_distance_bp}-{max_distance_bp} bp from a {gene_set} splice site, "
          f"{len(scored_loci):,} inside a GENCODE v50 {gene_set} transcript")
    return loci_near_splice_sites, scored_loci


def load_polymorphic_scored_loci(gene_set, max_distance_bp=NO_DISTANCE_LIMIT_BP, min_distance_bp=0):
    """Returns the loci SpliceAI would score that are polymorphic in HPRC256, with their target lengths.

    Streams HPRC256_ALLELE_LENGTH_STATS_TSV (5.6M rows, from compute_hprc256_total_allele_length_stats.py)
    instead of loading it, keeping only the rows of loci SpliceAI would score.

    Returns:
        tuple: (dict of counts at each filtering step, sorted list of
            (chrom, start_0based, end, motif, tuple of target lengths in bp; see target_lengths))
    """
    loci_in_distance_range, scored_loci = load_loci_near_splice_sites(gene_set, max_distance_bp, min_distance_bp)
    scored_loci = set(scored_loci)
    motif_range_bp = load_motif_ranges()
    n_with_stats = 0
    polymorphic = []
    with gzip.open(HPRC256_ALLELE_LENGTH_STATS_TSV, "rt") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            chrom, start, end, motif = row["locus_id"].rsplit("-", 3)
            locus = (f"chr{chrom}", int(start), int(end), motif)
            if locus not in scored_loci:
                continue
            n_with_stats += 1
            counts_by_length = parse_allele_length_histogram(row["allele_length_histogram"])
            if not is_polymorphic(counts_by_length):
                continue
            low_bp, high_bp, highest_bp = (percentile_of_histogram(counts_by_length, p)
                                           for p in (LOW_PERCENTILE, HIGH_PERCENTILE, HIGHEST_PERCENTILE))
            lengths_bp = tuple(target_lengths(low_bp, high_bp, highest_bp, motif_range_bp(motif)))
            # A locus whose every target length rounds to the reference (e.g. a long motif whose
            # percentiles differ by less than half a unit) would get no alleles.
            if labeled_size_changes_in_repeat_units(locus[2] - locus[1], len(motif), lengths_bp):
                polymorphic.append(locus + (lengths_bp,))
    counts = {
        "n_loci_in_distance_range": len(loci_in_distance_range),
        "n_loci_inside_a_transcript": len(scored_loci),
        "n_loci_with_hprc256_allele_lengths": n_with_stats,
        "n_polymorphic_loci": len(polymorphic),
    }
    print(f"{n_with_stats:,} of those have HPRC256 allele lengths, {len(polymorphic):,} with at least "
          f"{MIN_DISTINCT_ALLELE_SIZES_TO_BE_POLYMORPHIC} distinct allele sizes and at least one whole-unit change to simulate")
    return counts, sorted(polymorphic)


def is_polymorphic(counts_by_length):
    """True if the HPRC256 haplotypes have at least MIN_DISTINCT_ALLELE_SIZES_TO_BE_POLYMORPHIC distinct total allele sizes."""
    return len(counts_by_length) >= MIN_DISTINCT_ALLELE_SIZES_TO_BE_POLYMORPHIC


def simulate_alleles_for_locus(fasta, chrom, start, end, motif, lengths_bp):
    """Returns [("chrom-pos-ref-alt", [TARGET_LABELS it stands for])] for the alleles that bring one locus's tract to lengths_bp.

    Returns [] when the anchor base or the tract contains a base other than A, C, G or T (e.g. N), since
    the production scorer rejects such alleles.
    """
    tract = fasta[chrom][start:end].seq.upper()
    anchor = fasta[chrom][start - 1:start].seq.upper()
    if set(anchor + tract) - set("ACGT"):
        return []
    unit = most_common_in_frame_unit(tract, len(motif))
    reference_count = len(tract) // len(motif)
    return [(allele_to_variant(chrom, start, anchor, tract, unit, reference_count, reference_count + change), labels)
            for change, labels in labeled_size_changes_in_repeat_units(len(tract), len(motif), lengths_bp).items()]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--gene-set", choices=["basic", "comprehensive"], required=True)
    parser.add_argument("--n-loci", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--fasta", default=os.path.expanduser("~/hg38.fa"))
    parser.add_argument("--min-distance-bp", type=int, default=0,
                        help="Sample loci at least this far from the nearest splice site (default 0)")
    parser.add_argument("--max-distance-bp", type=int, default=NO_DISTANCE_LIMIT_BP,
                        help="Sample loci at most this far from the nearest splice site (default: no limit)")
    args = parser.parse_args()

    import pyfaidx
    fasta = pyfaidx.Fasta(args.fasta)
    counts, polymorphic_loci = load_polymorphic_scored_loci(args.gene_set, args.max_distance_bp, args.min_distance_bp)
    sampled_loci = random.Random(args.seed).sample(polymorphic_loci, args.n_loci)
    # A sampled locus with no allele (an N in its anchor or tract) is left out, as in the full run.
    sampled_loci_with_alleles = [(locus, simulate_alleles_for_locus(fasta, *locus)) for locus in sampled_loci]
    sampled_loci = [locus for locus, alleles in sampled_loci_with_alleles if alleles]
    labeled_alleles_by_locus = [alleles for _, alleles in sampled_loci_with_alleles if alleles]
    alleles_by_locus = [[variant for variant, _ in alleles] for alleles in labeled_alleles_by_locus]

    # The default range keeps the plain file name, which the cost and disk benchmark reads by default.
    is_default_range = (args.min_distance_bp, args.max_distance_bp) == (0, NO_DISTANCE_LIMIT_BP)
    distance_suffix = "" if is_default_range else f".{args.min_distance_bp}-{args.max_distance_bp}bp"
    output_path = os.path.join(HERE, "benchmark", f"tr_benchmark_alleles.{args.gene_set}{distance_suffix}.json")
    # Imported here: make_full_run_input_chunks imports this module.
    from make_full_run_input_chunks import compute_allele_design_fingerprint
    with open(output_path, "w") as f:
        json.dump({
            "gene_set": args.gene_set,
            # The benchmark refuses a sample whose fingerprint no longer matches the current allele design
            "allele_design_fingerprint": compute_allele_design_fingerprint(args.gene_set),
            "seed": args.seed,
            "distance_range_bp": [args.min_distance_bp, args.max_distance_bp],
            **counts,
            "n_loci_in_full_run": counts["n_polymorphic_loci"],
            "locus_ids": [f"{chrom}-{start}-{end}-{motif}" for chrom, start, end, motif, _ in sampled_loci],
            "alleles_by_locus": alleles_by_locus,
            # For each allele, the target lengths (TARGET_LABELS) it stands for
            "target_labels_by_locus": [[labels for _, labels in alleles] for alleles in labeled_alleles_by_locus],
        }, f, indent=1)
    size_changes = sorted(len(v.split("-")[3]) - len(v.split("-")[2]) for alleles in alleles_by_locus for v in alleles)
    print(f"Wrote {len(size_changes)} alleles at {len(alleles_by_locus)} loci to {output_path} "
          f"({sum(1 for s in size_changes if s < 0)} contractions); size change median "
          f"{size_changes[len(size_changes) // 2]} bp, range {size_changes[0]} to {size_changes[-1]} bp")


if __name__ == "__main__":
    main()
