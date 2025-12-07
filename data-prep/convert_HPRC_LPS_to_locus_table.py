"""
This script takes a TSV file with columns:

locus_ids   (example: "10-100000859-100000887-A" or "10-100000859-100000887-A,10-100001413-100001429-T"
motif       (example: "A")
sample_id:  (example: "HG00609")
allele_size (example: 25)

groups the rows by locus_id and motif columns, and then outputs a table with columns:

locus_id    (example: "10-100001413-100001429-T")
motif       (example: "T")
allele_sizes (example: 25,23,24,24,10,15)
"""

import argparse
import collections
import gzip
import os
import numpy as np
import tqdm


def write_line(locus_id, motif, allele_sizes, alleles_for_current_key_by_sample_id, output_file):
    """Process a group of allele sizes and write statistics to output file.
    
    Args:
        locus_id (str): the locus id (can be a comma-separated list of locus ids when it's a variation cluster)
        motif (str): the motif (in variation clusters, this will correspond to the motif at the end of one of the locus ids)
        allele_sizes (list): the allele sizes for all samples
        alleles_for_current_key_by_sample_id (dict): the alleles for the current key by sample id
        output_file (file): the output file
    """
    
    if not allele_sizes:
        return

    if "," in locus_id:
        found_locus_id = None
        for specific_locus_id in locus_id.split(","):
            if specific_locus_id.endswith(f"-{motif}"):
                found_locus_id = specific_locus_id
                break
        else:
            raise ValueError(f"Couldn't resolve locus id for motif {motif} in {locus_id}")

        locus_id = found_locus_id
    
    # Calculate statistics
    allele_sizes = list(sorted(allele_sizes))
    allele_counts = collections.Counter(allele_sizes)
    mode_allele, _ = allele_counts.most_common(1)[0]

    genotype_counts = collections.defaultdict(int)
    for sample_id, allele_list in alleles_for_current_key_by_sample_id.items():
        if len(allele_list) == 1:
            allele_list = allele_list * 2
        elif len(allele_list) != 2:
            raise ValueError(f"Found {len(allele_list)} alleles for {sample_id} in {locus_id} {motif}")

        genotype_counts[tuple(sorted(allele_list))] += 1

    allele_size_histogram = ",".join(f"{allele_size}x:{count}" for allele_size, count in sorted(allele_counts.items()))
    biallelic_histogram = ",".join(f"{genotype[0]}/{genotype[1]}:{count}" for genotype, count in sorted(genotype_counts.items(), key=lambda x: (x[0][1], x[0][0])))

    # Write to output file
    output_file.write("\t".join(map(str, [
        locus_id, 
        motif, 
        allele_size_histogram, 
        biallelic_histogram,
        int(min(allele_sizes)),
        mode_allele, 
        f"{np.mean(allele_sizes):.3f}",
        f"{np.std(allele_sizes):.3f}", 
        int(np.median(allele_sizes)),
        int(np.percentile(allele_sizes, 99)),
        int(max(allele_sizes)),
        len(set(allele_sizes)),
        len(allele_sizes),
    ])) + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-table", default="hprc-lps_2025-12-06/hprc-lps.txt.gz", help="Combine HPRC LPS dataset")
    parser.add_argument("-n", type=int, default=None, help="Number of samples to process")
    args = parser.parse_args()

    if not os.path.isfile(args.input_table):
        parser.error(f"Input file {args.input_table} does not exist")

    if args.n is not None:
        output_path = f"hprc_lps.{args.n}_samples.grouped_by_locus_and_motif.with_biallelic_histogram.tsv.gz"
    else:
        output_path = "hprc_lps.2025_12.grouped_by_locus_and_motif.with_biallelic_histogram.tsv.gz"

    """
    Input table example row:
    $1      trid : 1-92675-92688-AAG
    $2    allele : 1
    $3     motif : AAG
    $4   lps_len : 4
    """

    fopen = gzip.open if args.input_table.endswith("gz") else open
    with fopen(args.input_table, "rt") as infile, gzip.open(output_path, "wt") as outfile:
        header_fields = next(infile).strip().split("\t")
        if header_fields[0:2] != ["trid", "motif"]:
            parser.error(f"Expected header to start with 'trid', 'motif' columns in {args.input_table}. Found: {header_fields}")

        # Write header
        outfile.write("\t".join([
            "locus_id",
            "motif",
            "allele_size_histogram",
            "biallelic_histogram",
            "min_allele",
            "mode_allele",
            "mean",
            "stdev",
            "median",
            "99th_percentile",
            "max_allele",
            "unique_alleles",
            "num_called_alleles",
        ]) + "\n")


        for line_number, line in tqdm.tqdm(enumerate(infile), unit=" lines", unit_scale=True):
            fields = line.strip().split("\t")
            locus_id = fields[0]
            motif = fields[1]

            alleles = []
            alleles_by_sample_id = collections.defaultdict(list)
            for sample_id, allele_sizes in zip(header_fields[2:], fields[2:]):
                for allele_size in allele_sizes.split(","):
                    try:
                        allele_size = int(allele_size)
                    except ValueError:
                        raise ValueError(f"Expected integer allele size, got {allele_size} in line #{line_number + 1}: {line}")
                    alleles.append(allele_size)
                    alleles_by_sample_id[sample_id].append(allele_size)

            write_line(locus_id, motif, alleles, alleles_by_sample_id, outfile)



    print(f"Wrote {line_number + 1:9,d} lines to {output_path}")

if __name__ == "__main__":
    main()
