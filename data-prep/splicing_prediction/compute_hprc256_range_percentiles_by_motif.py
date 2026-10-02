"""For each repeat motif, find how much loci with that motif vary in HPRC256: a high percentile of their ranges.

For every locus in HPRC256_total_allele_length_stats.tsv.gz (compute_hprc256_total_allele_length_stats.py)
with at least MIN_DISTINCT_ALLELE_LENGTHS different total allele lengths among its haplotypes, the range
is the 97.5th minus the 2.5th percentile total allele length, in bp. Less variable loci are
left out: most loci of common motifs barely vary (e.g. short coding CAG repeats), which over all loci
would put the 99th percentile of CAG/CTG at only 3.7 units. Loci are grouped
by canonical motif (the same motif in any rotation or on either strand, e.g. CAG = AGC = GCA = CTG = TGC
= GCT), and separately by motif length. For each group the script writes the MOTIF_RANGE_PERCENTILE of its
loci's ranges, which make_tr_expansion_benchmark_alleles.py uses to size the larger simulated expansions:
how far a motif can vary somewhere in the genome, rather than how far it varies at one locus. Every
locus in a group has the same motif length, so ranking the ranges in bp picks the same locus as ranking
them in repeat units, and the result is a whole number of bp: the design's target lengths and size
changes then stay in exact integer arithmetic (a range stored in repeat units, e.g. 424 bp / 38 =
11.157894736842104, comes back as 423.99999999999994 bp and can round a half-unit target the wrong way).

Writes HPRC256_range_percentile_by_motif.tsv with columns:
    group_type ("canonical_motif" or "motif_length")  group  n_loci  range_percentile_bp
(MOTIF_RANGE_PERCENTILE and the other settings are in make_tr_expansion_benchmark_alleles.py.)

Usage:
    python3 compute_hprc256_range_percentiles_by_motif.py
"""
import collections
import csv
import gzip
import os

from make_tr_expansion_benchmark_alleles import (
    HIGH_PERCENTILE, HPRC256_ALLELE_LENGTH_STATS_TSV, LOW_PERCENTILE, MOTIF_RANGE_PERCENTILE,
    MOTIF_RANGE_PERCENTILES_TSV, canonical_motif, parse_allele_length_histogram, percentile_of_histogram)

MIN_DISTINCT_ALLELE_LENGTHS = 5


def percentile_of_values(values, percentile):
    """Nearest-rank percentile of a list of numbers."""
    values = sorted(values)
    return values[max(0, -(-len(values) * percentile // 100) - 1)]


def main():
    ranges_by_group = collections.defaultdict(list)
    with gzip.open(HPRC256_ALLELE_LENGTH_STATS_TSV, "rt") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            motif = row["locus_id"].rsplit("-", 1)[1]
            counts_by_length = parse_allele_length_histogram(row["allele_length_histogram"])
            if len(counts_by_length) < MIN_DISTINCT_ALLELE_LENGTHS:
                continue
            range_bp = (percentile_of_histogram(counts_by_length, HIGH_PERCENTILE)
                        - percentile_of_histogram(counts_by_length, LOW_PERCENTILE))
            ranges_by_group[("canonical_motif", canonical_motif(motif))].append(range_bp)
            ranges_by_group[("motif_length", str(len(motif)))].append(range_bp)

    with open(MOTIF_RANGE_PERCENTILES_TSV, "w") as f:
        f.write("group_type\tgroup\tn_loci\trange_percentile_bp\n")
        for (group_type, group), ranges in sorted(ranges_by_group.items()):
            f.write(f"{group_type}\t{group}\t{len(ranges)}\t{percentile_of_values(ranges, MOTIF_RANGE_PERCENTILE)}\n")
    print(f"Wrote {len(ranges_by_group):,} groups to {os.path.basename(MOTIF_RANGE_PERCENTILES_TSV)}")


if __name__ == "__main__":
    main()
