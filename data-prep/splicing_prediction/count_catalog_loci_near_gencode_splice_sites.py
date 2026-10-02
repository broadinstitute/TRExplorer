"""Counts tandem repeat catalog loci within various distances of GENCODE splice sites.

A splice site here is the first or last base of an annotated intron (the base next to the exon
boundary) in any transcript with 2 or more exons. The distance from a locus to a splice site is 0
when the site lies inside the repeat tract, and otherwise the number of bp between the site and the
nearest tract end.

The distance cutoffs come from earlier SpliceAI runs on simulated expansions at known disease loci:
changes >= 0.5 all fell within 188 bp of the repeat,
95% of changes >= 0.2 fell within ~500 bp, all changes >= 0.1 fell within ~3.3 kb, and SpliceAI
cannot see anything farther than 5 kb away.

Outputs:
  <output_prefix>.distances.tsv.gz: one row per catalog locus with its distance to the nearest
      splice site in each GTF.
  <output_prefix>.summary.tsv: number and fraction of loci within each distance cutoff.
"""

import argparse
import collections
import gzip

import numpy as np

DISTANCE_CUTOFFS = [0, 10, 50, 188, 500, 1000, 3300, 5000]


def parse_splice_site_positions_from_gtf(gtf_path):
    """Returns {chrom: sorted numpy array of 1-based splice site positions} from a GENCODE GTF.

    Args:
        gtf_path: path of a gzipped GTF file.

    Returns:
        Dict mapping chromosome name to a sorted, deduplicated array of the 1-based positions of the
        first and last base of every intron.
    """
    exons_by_transcript = collections.defaultdict(list)
    with gzip.open(gtf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.split("\t", 9)
            if fields[2] != "exon":
                continue
            transcript_id = fields[8].split('transcript_id "', 1)[1].split('"', 1)[0]
            exons_by_transcript[(fields[0], transcript_id)].append((int(fields[3]), int(fields[4])))

    positions_by_chrom = collections.defaultdict(set)
    for (chrom, _), exons in exons_by_transcript.items():
        for intron_start, intron_end in compute_intron_boundaries(exons):
            positions_by_chrom[chrom].update((intron_start, intron_end))

    return {chrom: np.array(sorted(positions)) for chrom, positions in positions_by_chrom.items()}


def compute_intron_boundaries(exons):
    """Returns (first intron base, last intron base) 1-based pairs for a transcript's exons.

    Args:
        exons: list of (start, end) 1-based inclusive exon coordinates, in any order.

    Returns:
        List of (intron_start, intron_end) tuples.
    """
    exons = sorted(exons)
    return [(previous_end + 1, next_start - 1)
            for (_, previous_end), (next_start, _) in zip(exons, exons[1:])
            if next_start - 1 >= previous_end + 1]


def compute_distances_to_nearest_splice_site(start_0based, end, sorted_positions):
    """Returns the distance in bp from each tract to the nearest splice site position.

    Args:
        start_0based: numpy array of tract start coordinates (0-based).
        end: numpy array of tract end coordinates (1-based inclusive).
        sorted_positions: sorted numpy array of 1-based splice site positions on this chromosome.

    Returns:
        numpy array of distances (0 if a splice site lies inside the tract), or -1 if the
        chromosome has no splice sites.
    """
    if len(sorted_positions) == 0:
        return np.full(len(start_0based), -1)

    # The first splice site at or after the tract's first base, and the last one before it.
    index = np.searchsorted(sorted_positions, start_0based + 1)
    next_site = sorted_positions[np.minimum(index, len(sorted_positions) - 1)]
    previous_site = sorted_positions[np.maximum(index - 1, 0)]

    distance_to_next = np.where(index < len(sorted_positions), np.maximum(next_site - end, 0), np.iinfo(np.int64).max)
    distance_to_previous = np.where(index > 0, start_0based + 1 - previous_site, np.iinfo(np.int64).max)
    return np.minimum(distance_to_next, distance_to_previous)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--catalog-bed", required=True, help="Catalog BED (chrom, 0-based start, end, motif, ...)")
    parser.add_argument("--gtf", action="append", required=True, help="GENCODE GTF; can be given more than once")
    parser.add_argument("--output-prefix", required=True)
    args = parser.parse_args()

    loci_by_chrom = collections.defaultdict(list)
    with gzip.open(args.catalog_bed, "rt") as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            loci_by_chrom[fields[0]].append((int(fields[1]), int(fields[2]), fields[3]))
    total_loci = sum(len(loci) for loci in loci_by_chrom.values())
    print(f"Read {total_loci:,d} loci from {args.catalog_bed}")

    distances_by_gtf = {}
    for gtf_path in args.gtf:
        positions_by_chrom = parse_splice_site_positions_from_gtf(gtf_path)
        print(f"Parsed {sum(len(p) for p in positions_by_chrom.values()):,d} unique splice site positions from {gtf_path}")
        distances_by_gtf[gtf_path] = {
            chrom: compute_distances_to_nearest_splice_site(
                np.array([locus[0] for locus in loci]), np.array([locus[1] for locus in loci]),
                positions_by_chrom.get(chrom, np.array([], dtype=int)))
            for chrom, loci in loci_by_chrom.items()
        }

    gtf_labels = [gtf_path.split("/")[-1].replace(".gtf.gz", "") for gtf_path in args.gtf]
    with gzip.open(f"{args.output_prefix}.distances.tsv.gz", "wt") as f:
        f.write("\t".join(["chrom", "start_0based", "end", "motif"] + [
            f"distance_to_nearest_splice_site_bp__{label}" for label in gtf_labels]) + "\n")
        for chrom, loci in loci_by_chrom.items():
            for i, (start_0based, end, motif) in enumerate(loci):
                f.write("\t".join(map(str, [chrom, start_0based, end, motif] + [
                    distances_by_gtf[gtf_path][chrom][i] for gtf_path in args.gtf])) + "\n")

    with open(f"{args.output_prefix}.summary.tsv", "w") as f:
        f.write("gtf\tmax_distance_bp\tloci_within_distance\tfraction_of_catalog\n")
        for gtf_path, label in zip(args.gtf, gtf_labels):
            all_distances = np.concatenate(list(distances_by_gtf[gtf_path].values()))
            all_distances = all_distances[all_distances >= 0]
            for cutoff in DISTANCE_CUTOFFS:
                count = int((all_distances <= cutoff).sum())
                row = f"{label}\t{cutoff}\t{count}\t{count / total_loci:.4f}"
                f.write(row + "\n")
                print(row)

    print(f"Wrote {args.output_prefix}.distances.tsv.gz")
    print(f"Wrote {args.output_prefix}.summary.tsv")


if __name__ == "__main__":
    main()
