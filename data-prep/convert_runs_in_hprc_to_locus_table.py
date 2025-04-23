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
parser = argparse.ArgumentParser()
parser.add_argument("--input-tsv", type=str, default="runs_in_hprc.2025_04.txt.gz")
parser.add_argument("--output-path", type=str, default="runs_in_hprc.2025_04.grouped_by_locus_and_motif.tsv.gz")
args = parser.parse_args()

# check that file exists
if not os.path.exists(args.input_tsv):
    parser.error(f"File {args.input_tsv} does not exist")

def process_group(locus_id, motif, allele_sizes, output_file):
    """Process a group of allele sizes and write statistics to output file."""
    
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

    stdev = np.std(allele_sizes)
    median = np.median(allele_sizes)
    percentile_99 = np.percentile(allele_sizes, 99)
    allele_size_histogram = ",".join(f"{allele_size}x:{count}" for allele_size, count in allele_counts.items())
    
    # Write to output file
    output_file.write("\t".join(map(str, [
        locus_id, 
        motif, 
        allele_size_histogram, 
        mode_allele, 
        f"{stdev:.3f}", 
        int(median), 
        int(percentile_99), 
    ])) + "\n")

def main():

    print(f"Processing {os.path.basename(args.input_tsv)}")
    infile = gzip.open(args.input_tsv, "rt")
    outfile = gzip.open(args.output_path, "wt")

    # Write header
    outfile.write("\t".join([
        "locus_id", 
        "motif", 
        "allele_size_histogram", 
        "mode_allele", 
        "stdev", 
        "median", 
        "99th_percentile",
    ]) + "\n")

    previous_key = None
    previously_seen_locus_ids = set()
    alleles_for_current_key = []
    sample_ids = set()
    counters = collections.Counter()
    for line_number, line in tqdm.tqdm(enumerate(infile, 1), unit=" rows", unit_scale=True, total=896_015_950):
        fields = line.strip().split("\t")
        if len(fields) != 4:
            print(f"WARNING: Skipping malformed line: {line}")
            continue

        current_locus_id, motif, sample_id, allele_size = fields
        
        sample_ids.add(sample_id)

        current_key = (current_locus_id, motif)
        if current_key != previous_key or previous_key is None:
            if current_key in previously_seen_locus_ids:
                parser.error(f"{args.input_tsv} is not sorted by locus id on line #{line_number}: {current_locus_id}") 
            
            previously_seen_locus_ids.add(current_key)
            previous_key = current_key

            process_group(current_locus_id, motif, alleles_for_current_key, outfile)
            counters['output_lines'] += 1

            alleles_for_current_key = []
            
        try:
            allele_size = int(allele_size)
            alleles_for_current_key.append(allele_size)
        except ValueError:
            print(f"Warning: Skipping invalid allele at line #{line_number}: {current_locus_id} {motif} {sample_id} {allele_size}")
            continue
                
    # Process the final group
    if alleles_for_current_key:
        process_group(current_locus_id, motif, alleles_for_current_key, outfile)
        counters['output_lines'] += 1

    infile.close()
    outfile.close()

    print(f"Wrote {counters['output_lines']:9,d} lines to {args.output_path} for {len(sample_ids):9,d} unique sample ids")

if __name__ == "__main__":
    main()