"""Same as compare_simulated_sizes_to_pathogenic_thresholds.py, but with HPRC256 longest pure segment (LPS) lengths.

The simulation design uses total allele lengths. This version applies the same rule to the length of
the longest uninterrupted run of the motif (LPS, in repeat units), the measure Danzi et al. found most
variable at disease loci: the 2.5th, 97.5th and 99.5th percentile LPS of each STRchive disease locus,
and the 99.5th plus 1x, 2x and 3x its motif's LPS range (the MOTIF_RANGE_PERCENTILE, across loci with
that canonical motif and at least MIN_DISTINCT_LENGTHS distinct LPS lengths, of the 97.5th-2.5th
percentile LPS range; motifs with fewer than MIN_LOCI_FOR_A_MOTIF_SPECIFIC_RANGE such loci use the
value for all motifs of their length). LPS is already in repeat units, so no conversion to STRchive's
scale is needed; insertions are capped at MAX_INSERTED_BP beyond STRchive's reference count.

Usage:
    python3 compare_simulated_sizes_to_pathogenic_thresholds_using_lps.py
"""
import collections
import gzip
import json

# Also puts the parent folder, which has the pipeline modules, on the import path
from compare_simulated_sizes_to_pathogenic_thresholds import NOT_COMPARABLE_BY_LENGTH, STRCHIVE_LOCI_JSON
from make_tr_expansion_benchmark_alleles import (  # noqa: E402
    HIGH_PERCENTILE, HIGHEST_PERCENTILE, LOW_PERCENTILE, MAX_INSERTED_BP, MIN_LOCI_FOR_A_MOTIF_SPECIFIC_RANGE,
    MOTIF_RANGE_PERCENTILE, canonical_motif, percentile_of_histogram)

LPS_TSV = "/Users/weisburd/code/tandem-repeat-catalogs/hprc_lps.2025_12.per_locus_and_motif.256_samples.tsv.gz"
MIN_DISTINCT_LENGTHS = 5
MOTIF_RANGE_MULTIPLES = (1, 2, 3)


def parse_lps_histogram(histogram):
    """Parses "9x:130,10x:153,..." into {LPS in repeat units: number of haplotypes}."""
    return {int(size.rstrip("x")): int(count) for size, count in (pair.split(":") for pair in histogram.split(","))}


def nearest_rank_percentile(values, percentile):
    values = sorted(values)
    return values[max(0, -(-len(values) * percentile // 100) - 1)]


def read_lps_table(strchive_loci):
    """One pass over the LPS table.

    Returns:
        tuple: ({STRchive id: (motif, {LPS: count})} for the best-matching row, a function giving a
            motif's LPS range in repeat units)
    """
    by_chrom = collections.defaultdict(list)
    for locus in strchive_loci:
        by_chrom[locus["chrom"].replace("chr", "")].append(locus)
    ranges_by_group = collections.defaultdict(list)
    best = {}
    with gzip.open(LPS_TSV, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        locus_i, motif_i, histogram_i = header.index("locus_id"), header.index("motif"), header.index("allele_size_histogram")
        for line in f:
            fields = line.rstrip("\n").split("\t")
            chrom, start, end, _ = fields[locus_i].rsplit("-", 3)
            motif = fields[motif_i]
            counts = parse_lps_histogram(fields[histogram_i])
            if len(counts) >= MIN_DISTINCT_LENGTHS:
                lps_range = percentile_of_histogram(counts, HIGH_PERCENTILE) - percentile_of_histogram(counts, LOW_PERCENTILE)
                ranges_by_group[("motif", canonical_motif(motif))].append(lps_range)
                ranges_by_group[("length", len(motif))].append(lps_range)
            for locus in by_chrom.get(chrom, []):
                overlap = min(int(end), int(locus["stop_hg38"])) - max(int(start), int(locus["start_hg38"]))
                if len(motif) != int(locus["motif_len"]) or overlap <= 0:
                    continue
                # Prefer the row whose motif is the pathogenic one, then the most overlap
                key = (motif in (locus.get("pathogenic_motif_reference_orientation") or []), overlap)
                if key > best.get(locus["id"], ((False, 0),))[0]:
                    best[locus["id"]] = (key, motif, counts)

    by_motif = {group: nearest_rank_percentile(v, MOTIF_RANGE_PERCENTILE) for (kind, group), v in ranges_by_group.items()
                if kind == "motif" and len(v) >= MIN_LOCI_FOR_A_MOTIF_SPECIFIC_RANGE}
    by_length = {group: nearest_rank_percentile(v, MOTIF_RANGE_PERCENTILE) for (kind, group), v in ranges_by_group.items()
                 if kind == "length"}

    def motif_lps_range(motif):
        if canonical_motif(motif) in by_motif:
            return by_motif[canonical_motif(motif)], False
        return by_length[min(by_length, key=lambda length: (abs(length - len(motif)), length))], True
    return {strchive_id: match[1:] for strchive_id, match in best.items()}, motif_lps_range


def main():
    strchive_loci = sorted((x for x in json.load(open(STRCHIVE_LOCI_JSON)) if x.get("pathogenic_min")), key=lambda x: x["gene"])
    matches, motif_lps_range = read_lps_table(strchive_loci)
    genes = collections.Counter(x["gene"] for x in strchive_loci)
    shown = ["2.5th", "97.5th", "99.5th"] + [f"+{k}x" for k in MOTIF_RANGE_MULTIPLES]
    header = ["Gene", "Motif", "Motif LPS range (units)", "Ref"] + shown + ["Pathogenic", "Smallest size reaching it", "LPS polymorphic"]
    print("| " + " | ".join(header) + " |\n|" + "---|" * len(header))
    seen, tally = collections.Counter(), collections.Counter()
    for x in strchive_loci:
        seen[x["gene"]] += 1
        label = (f"{x['gene']} ({seen[x['gene']]})" if genes[x["gene"]] > 1 else x["gene"]) + (" *" if x["gene"] in NOT_COMPARABLE_BY_LENGTH else "")
        pathogenic = float(x["pathogenic_min"])
        reference = f"{float(x['ref_copies']):.0f}" if x.get("ref_copies") else ""
        if x["id"] not in matches:
            print(f"| {label} | {x['motif_len']} bp | | {reference} |" + " |" * len(shown) + f" {pathogenic:.0f} | no LPS match | |")
            continue
        motif, counts = matches[x["id"]]
        low, high, highest = (percentile_of_histogram(counts, p) for p in (LOW_PERCENTILE, HIGH_PERCENTILE, HIGHEST_PERCENTILE))
        motif_range, pooled = motif_lps_range(motif)
        cap = (float(x["ref_copies"]) if x.get("ref_copies") else highest) + MAX_INSERTED_BP // len(motif)
        sizes = [low, high, highest] + [min(highest + k * motif_range, cap) for k in MOTIF_RANGE_MULTIPLES]
        smallest = next((name for name, size in zip(shown, sizes) if size >= pathogenic), "none")
        motif_label = canonical_motif(motif) if len(motif) <= 12 else f"{len(motif)} bp"
        print(f"| {label} | {motif_label} | {motif_range:.0f}{' (pooled)' if pooled else ''} | {reference} | "
              + " | ".join(f"{size:.0f}" for size in sizes)
              + f" | {pathogenic:.0f} | {smallest} | {'yes' if low != high else 'no'} |")
        if x["gene"] not in NOT_COMPARABLE_BY_LENGTH:
            tally["comparable"] += 1
            tally[smallest] += 1
    print(f"\n{tally['comparable']} loci comparable by length; smallest size that reaches the pathogenic threshold:")
    for name in shown + ["none"]:
        print(f"  {name:<7} {tally[name]:>3} loci")


if __name__ == "__main__":
    main()
