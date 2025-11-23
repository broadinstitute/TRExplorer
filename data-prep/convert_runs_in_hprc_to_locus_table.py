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
import glob
import gzip
import os
import numpy as np
import tqdm

INPUT_PATHS = "hprc-lps/*.txt.gz"



def process_group(locus_id, motif, allele_sizes, alleles_for_current_key_by_sample_id, output_file):
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
    parser.add_argument("-n", type=int, default=None, help="Number of samples to process")
    args = parser.parse_args()

    input_paths = list(sorted(glob.glob(INPUT_PATHS)))
    
    assert len(input_paths) == 100, f"Expected 100 input paths, got {len(input_paths)}"

    if args.n is not None:
        if args.n > len(input_paths):
            parser.error(f"Requested {args.n} samples, but only {len(input_paths)} input paths available")
        input_paths = input_paths[:args.n]

    
    """
    Input table example row:
    $1      trid : 1-92675-92688-AAG
    $2    allele : 1
    $3     motif : AAG
    $4   lps_len : 4
    """

    alleles_by_locus = collections.defaultdict(list)
    alleles_by_locus_and_sample_id = collections.defaultdict(lambda: collections.defaultdict(list))

    for sample_number, input_path in enumerate(input_paths):
        sample_id = os.path.basename(input_path).split(".")[0]

        print(f"Processing sample #{sample_number+1:3d}: {sample_id}: {os.path.basename(input_path)}")
        infile = gzip.open(input_path, "rt")
        header_line = next(infile) # skip header
        assert header_line.strip() == "\t".join(["trid", "allele", "motif", "lps_len"]), f"Unexpected header line: {header_line.strip()}"

        # iterate through the file, reading two lines at a time
        for line_number, line in enumerate(infile):
            line_number += 1

            # process line #1
            fields = line.strip().split("\t")

            if len(fields) != 4:
                raise ValueError(f"Expected 4 fields, got {len(fields)} in line #{line_number}: {line}")

            locus_id, first_or_second_allele, motif, allele_size = fields

            try:
                allele_size = int(allele_size)
            except ValueError:
                raise ValueError(f"Expected integer allele size, got {allele_size} in line #{line_number}: {line}")

            current_key = (locus_id, motif)

            alleles_by_locus[current_key].append(allele_size)
            alleles_by_locus_and_sample_id[current_key][sample_id].append(allele_size)


        infile.close()


    if args.n is not None:
        output_path = f"hprc_lps.{args.n}_samples.grouped_by_locus_and_motif.with_biallelic_histogram.tsv.gz"
    else:   
        output_path = "hprc_lps.2025_05.grouped_by_locus_and_motif.with_biallelic_histogram.tsv.gz"


    outfile = gzip.open(output_path, "wt")

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


    output_lines_counter = 0
    for (locus_id, motif) in tqdm.tqdm(alleles_by_locus, unit=" loci", unit_scale=True):
        key = (locus_id, motif)
        process_group(locus_id, motif, alleles_by_locus[key], alleles_by_locus_and_sample_id[key], outfile)
        output_lines_counter += 1

    outfile.close()

    print(f"Wrote {output_lines_counter:9,d} lines to {output_path}")

if __name__ == "__main__":
    main()
