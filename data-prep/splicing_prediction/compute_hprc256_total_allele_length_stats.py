"""Summarize the HPRC256 total allele length distribution of every tandem repeat locus.

Reads the TRGT multisample VCF of the 256 HPRC samples (trgt-hprc.unique_trids.vcf.gz, genotyped
with the TRExplorer catalog) and, for every single-locus record (INFO/STRUC "<TR:...>"), tallies the
allele lengths (FORMAT/AL: the repeat's length in bp, without the padding base) of every called
haplotype. Records for variation clusters ("<VC:...>") are skipped: their AL spans the whole cluster,
and every locus in a cluster also has its own single-locus record in this VCF.

Writes a TSV with one row per locus:
    locus_id  motif_size  n_called_haplotypes  mode_bp  mean_bp  stdev_bp  allele_length_histogram
where stdev_bp is the population standard deviation over called haplotypes, mode_bp is the most
common length (ties go to the length closest to the reference, then the shorter one), and the
histogram is "length:count" pairs sorted by length.

Usage:
    python3 compute_hprc256_total_allele_length_stats.py [--n-processes 8]
"""
import argparse
import collections
import gzip
import math
import multiprocessing
import os
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
VCF = "/Users/weisburd/code/tandem-repeat-explorer/data-prep/hprc-lps_2026-05-19/trgt-hprc.unique_trids.vcf.gz"
OUTPUT = os.path.join(HERE, "HPRC256_total_allele_length_stats.tsv.gz")
CHROMS = [f"chr{c}" for c in list(range(1, 23)) + ["X", "Y", "M"]]


def summarize_allele_lengths(locus_id, allele_lengths):
    """Returns the output row (as a list of strings) for one locus's called allele lengths."""
    _, start, end, motif = locus_id.rsplit("-", 3)
    reference_length = int(end) - int(start)
    counts = collections.Counter(allele_lengths)
    n = len(allele_lengths)
    mean = sum(allele_lengths) / n
    stdev = math.sqrt(sum((x - mean) ** 2 for x in allele_lengths) / n)
    mode = min(counts, key=lambda length: (-counts[length], abs(length - reference_length), length))
    histogram = ",".join(f"{length}:{counts[length]}" for length in sorted(counts))
    return [locus_id, str(len(motif)), str(n), str(mode), f"{mean:.3f}", f"{stdev:.4f}", histogram]


def parse_genotype_fields(sample_fields):
    """Returns the allele lengths of the called haplotypes in "GT:AL" sample fields, e.g. ["0/1:9,12", "./.:."]."""
    allele_lengths = []
    for field in sample_fields:
        genotype, lengths = field.split(":")
        for allele, length in zip(genotype.replace("|", "/").split("/"), lengths.split(",")):
            if allele != "." and length != ".":
                allele_lengths.append(int(length))
    return allele_lengths


def summarize_chromosome(chrom):
    """Returns the output rows for one chromosome's single-locus records."""
    command = ["bcftools", "query", "-r", chrom, "-f", "%INFO/TRID\t%INFO/STRUC[\t%GT:%AL]\n", VCF]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    rows = []
    for line in process.stdout:
        locus_id, structure, *sample_fields = line.rstrip("\n").split("\t")
        if not structure.startswith("<TR:"):
            continue
        allele_lengths = parse_genotype_fields(sample_fields)
        if allele_lengths:
            rows.append(summarize_allele_lengths(locus_id, allele_lengths))
    if process.wait() != 0:
        raise RuntimeError(f"bcftools query failed for {chrom}")
    print(f"{chrom}: {len(rows):,} loci", flush=True)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--n-processes", type=int, default=8)
    args = parser.parse_args()
    with multiprocessing.Pool(args.n_processes) as pool:
        rows_by_chrom = pool.map(summarize_chromosome, CHROMS, chunksize=1)
    with gzip.open(OUTPUT, "wt") as f:
        f.write("\t".join(["locus_id", "motif_size", "n_called_haplotypes", "mode_bp", "mean_bp", "stdev_bp",
                           "allele_length_histogram"]) + "\n")
        for rows in rows_by_chrom:
            for row in rows:
                f.write("\t".join(row) + "\n")
    print(f"Wrote {sum(len(rows) for rows in rows_by_chrom):,} loci to {OUTPUT}")


if __name__ == "__main__":
    main()
