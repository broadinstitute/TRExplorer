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
import pandas as pd
import numpy as np
import tqdm

    
HEADER_FIELDS = [
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
]

def write_line(locus_id, motif, allele_sizes, alleles_by_sample_id, output_file):
    """Process a group of allele sizes and write statistics to output file.
    
    Args:
        locus_id (str): the locus id (can be a comma-separated list of locus ids when it's a variation cluster)
        motif (str): the motif (in variation clusters, this will correspond to the motif at the end of one of the locus ids)
        allele_sizes (list): the allele sizes for all samples
        alleles_by_sample_id (dict): the alleles for the current key by sample id
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
    for sample_id, allele_list in alleles_by_sample_id.items():
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


"""
Example row in 1kgp metadata TSV file:

SampleId                                                            HG00459
Sex                                                                    male
BiosampleId                                                      SAME125269
Subpopulation                                                           CHS
SubpopulationName                                      Southern Han Chinese
Population                                                              EAS
PopulationName                                          East Asian Ancestry
SubpopulationElasticId                                                  CHS
DataSource                1000 Genomes 30x on GRCh38,1000 Genomes phase ...

Population counts:
   1427 AFR
    535 AMR
    623 EAS
    672 EUR
      1 EUR,AFR
    661 SAS


 2343 female
   2653 male

Input table example row:

trid       1-19175-19184-TCC
motif                    TCC
HG00096                  3,3
HG00097                  3,3
HG00099                  3,3
...


"""


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-metadata-tsv", default="hprc-lps_2025-12-06/1kGP_metadata.tsv", help="Sample ancestry metadata TSV file")
    parser.add_argument("--input-table", default="hprc-lps_2025-12-06/hprc-lps.txt.gz", help="Combine HPRC LPS dataset")
    parser.add_argument("--population", choices=["AFR", "AMR", "EAS", "EUR", "SAS"], help="If specified, only process samples from this population")
    parser.add_argument("--sex", choices=["male", "female"], help="If specified, only process samples from this sex")
    parser.add_argument("-n", "--num-samples", type=int, default=None, help="Number of samples to process")
    parser.add_argument("-l", "--num-loci", type=int, default=None, help="Number of loci to process")
    args = parser.parse_args()

    if not os.path.isfile(args.input_table):
        parser.error(f"Input file {args.input_table} does not exist")


    fopen = gzip.open if args.input_table.endswith("gz") else open
    with fopen(args.input_table, "rt") as infile:
        header_fields = next(infile).strip().split("\t")
        if header_fields[0:2] != ["trid", "motif"]:
            parser.error(f"Expected header to start with 'trid', 'motif' columns in {args.input_table}. Found: {header_fields}")
    sample_ids_in_input_table = header_fields[2:]
    print(f"Parsed {len(sample_ids_in_input_table):,d} sample ids from {args.input_table}")

    df_metadata = pd.read_table(args.sample_metadata_tsv)
    all_sample_ids = set(df_metadata.SampleId)
    unexpected_sample_ids = set(sample_ids_in_input_table) - all_sample_ids
    if unexpected_sample_ids:
        parser.error(f"{len(unexpected_sample_ids):,d} sample id(s) not found in {args.sample_metadata_tsv}: {', '.join(sorted(unexpected_sample_ids))}")

    df_metadata = df_metadata[df_metadata.SampleId.isin(set(sample_ids_in_input_table))]
    if args.population:
        df_metadata = df_metadata[df_metadata.Population == args.population]
        print(f"Kept {len(df_metadata):,d} samples from population {args.population}")
    if args.sex:
        df_metadata = df_metadata[df_metadata.Sex == args.sex]
        print(f"Kept {len(df_metadata):,d} samples from sex {args.sex}")

    sample_ids_to_include = [s for s in sample_ids_in_input_table if s in set(df_metadata.SampleId)]
    if args.num_samples is not None and len(sample_ids_to_include) > args.num_samples:
        sample_ids_to_include = sample_ids_to_include[:args.num_samples]

    print(f"Included sample ids: {', '.join(sample_ids_to_include)}")
    output_dir = os.path.dirname(args.sample_metadata_tsv)    
    output_path = os.path.join(output_dir, "hprc_lps.2025_12.per_locus_and_motif")
    if args.population:
        output_path += f".only_{args.population}"
    if args.sex:
        output_path += f".only_{args.sex}"

    output_path += f".{len(sample_ids_to_include)}_samples"
    output_path += ".tsv.gz"
    print(f"Writing data from {len(sample_ids_to_include):,d} samples to {output_path}")
    with fopen(args.input_table, "rt") as infile, gzip.open(output_path, "wt") as outfile:
        next(infile)  # skip header

        # Write header
        outfile.write("\t".join(HEADER_FIELDS) + "\n")


        for line_number, line in tqdm.tqdm(enumerate(infile), unit=" lines", unit_scale=True):
            fields = line.strip().split("\t")
            locus_id = fields[0]
            motif = fields[1]
            if args.num_loci is not None and line_number >= args.num_loci:
                break

            alleles = []
            alleles_by_sample_id = collections.defaultdict(list)
            for sample_id, allele_sizes in zip(header_fields[2:], fields[2:]):
                if sample_id not in sample_ids_to_include:
                    continue
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
