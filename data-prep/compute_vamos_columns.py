"""This script takes the catalog and several files shared by Vamos authors, and outputs a TSV
with Vamos-related columns to include in the BigQuery table.
"""
import argparse
import collections
import gzip
import os
import requests
import tqdm

from str_analysis.utils.canonical_repeat_unit import compute_canonical_motif
from str_analysis.utils.misc_utils import get_json_iterator, parse_interval
from str_analysis.utils.eh_catalog_utils import group_overlapping_loci


def get_reference_region_size(record):
    chrom, start_0based, end_1based = parse_interval(record["ReferenceRegion"])
    return end_1based - start_0based

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--vamos-ori-motifs-tsv", default="../data-prep/oriMotifFinal.all.tsv.gz")
    parser.add_argument("--vamos-eff-motifs-tsv", default="../data-prep/effMotifFinal.all.tsv.gz")
    parser.add_argument("-c", "--catalog-path", help="Path to the annotated catalog JSON file",
                        #default="https://github.com/broadinstitute/tandem-repeat-catalog/releases/download/v1.0/repeat_catalog_v1.hg38.1_to_1000bp_motifs.EH.with_annotations.json.gz")
                        default="~/code/tandem-repeat-catalogs/results__2025-12-07/release_draft_2025-12-07/repeat_catalog_v2.hg38.1_to_1000bp_motifs.EH.with_annotations.json.gz")
    parser.add_argument("-n", type=int, help="Number of records to process from the catalog")
    parser.add_argument("-o", "--output-tsv", default="../data-prep/vamos_ori_and_eff_motif_columns.tsv.gz")
    args = parser.parse_args()

    vamos_data_lookup = parse_vamos_data(args.vamos_ori_motifs_tsv, args.vamos_eff_motifs_tsv)

    # Read catalog data
    print(f"Reading catalog from {args.catalog_path}")
    is_gzipped = args.catalog_path.endswith('gz')
    if args.catalog_path.startswith("http"):
        response = requests.get(args.catalog_path)
        catalog = get_json_iterator(response.content, is_gzipped)
    elif os.path.isfile(os.path.expanduser(args.catalog_path)):
        catalog = get_json_iterator(os.path.expanduser(args.catalog_path), is_gzipped)
    else:
        parser.error(f"Invalid catalog path: {args.catalog_path}")

    counters = collections.Counter()

    output_table_rows = []
    for record_group in group_overlapping_loci(
        catalog,
        only_group_loci_with_similar_motifs=False, # group loci regardless of whether they have the same motif since Vamos doesn't support overlapping locus definitions.
        min_overlap_size=1, # group loci even if they overlap by at least 1bp
        verbose=True,
    ):

        record_with_largest_reference_region = max(record_group, key=lambda record: get_reference_region_size(record))
        for record in record_group:
            counters["total"] += 1
            reference_region = record["ReferenceRegion"]
            output_record = {
                "ReferenceRegion": reference_region,
                "IncludeInVamosCatalog": 1 if record == record_with_largest_reference_region else 0,
            }
            if record == record_with_largest_reference_region:
                counters["records_included_in_vamos_catalog"] += 1

            if reference_region in vamos_data_lookup:
                counters["has_vamos_data"] += 1
                vamos_data = vamos_data_lookup[reference_region]
                output_record.update({
                    "VamosUniqueMotifs": vamos_data["unique_motifs"],
                    "VamosEfficientMotifs": vamos_data["efficient_motifs"],
                    "VamosMotifFrequencies": vamos_data["motif_frequencies"],
                    "VamosNumUniqueMotifs": vamos_data["num_unique_motifs"],
                })
            output_table_rows.append(output_record)

        if args.n and counters["total"] >= args.n:
            break


    print(f"Kept {len(output_table_rows):,d} out of {counters['total']:,d} ({100 * len(output_table_rows) / counters['total']:.2f}%) rows")

    print(f"Writing {len(output_table_rows):,d} rows to {args.output_tsv}")
    header = [
        "ReferenceRegion",
        "IncludeInVamosCatalog",
        "VamosUniqueMotifs",
        "VamosEfficientMotifs",
        "VamosMotifFrequencies",
        "VamosNumUniqueMotifs",
    ]

    fopen = gzip.open if args.output_tsv.endswith(".gz") else open
    with fopen(args.output_tsv, "wt") as f:
        f.write("\t".join(header) + "\n")
        for row in tqdm.tqdm(output_table_rows, unit=" rows", unit_scale=True, total=len(output_table_rows)):
            f.write("\t".join([str(row.get(col, "")) for col in header]) + "\n")

    print(f"Wrote {len(output_table_rows):,d} rows to {args.output_tsv}")

    print(f"{counters['records_included_in_vamos_catalog']:,d} out of {counters['total']:,d} ({counters['records_included_in_vamos_catalog']/counters['total']:0.1%}) were included in the vamos catalog")
    print(f"{counters['has_vamos_data']:,d} out of {counters['total']:,d} ({counters['has_vamos_data']/counters['total']:0.1%}) had vamos data")

def parse_vamos_data(unique_motif_tsv_path, efficient_motif_tsv_path=None):
    """Parse the Vamos motif frequencies TSV file line by line and create a lookup dictionary."""

    efficient_motif_data_available = efficient_motif_tsv_path is not None
    if efficient_motif_data_available:
        vamos_ori_to_efficient_motif_map = parse_vamos_efficient_motif_lookup(efficient_motif_tsv_path)
    else:
        vamos_ori_to_efficient_motif_map = None

    print(f"Parsing unique motif frequencies from {unique_motif_tsv_path}")
    result = {}
    with gzip.open(unique_motif_tsv_path, 'rt') as f:
        header = f.readline().strip().split('\t')

        col_indices = {col: i for i, col in enumerate(header)}
        for line in tqdm.tqdm(f, unit=" lines", unit_scale=True):
            fields = line.strip().split('\t')
            chrom = fields[col_indices['chr']]
            start_0based = int(fields[col_indices['start']])
            end_1based = int(fields[col_indices['end']])
            ori_motifs = fields[col_indices['motifs']]
            ori_motif_counts = fields[col_indices['motifCounts']]
            reference_region = f"{chrom}:{start_0based}-{end_1based}"

            if ori_motifs.count(",") != ori_motif_counts.count(","):
                print(f"WARNING: For {reference_region} number of motif values ({ori_motifs}) != number of motif counts ({ori_motif_counts})")
                continue

            total_counts = 0
            motif_and_count_list = []

            canonical_motif_to_motif = {}
            canonical_motif_to_count = {}
            unique_motif_list = []
            for motif, count in zip(ori_motifs.split(","), ori_motif_counts.split(",")):
                if motif in unique_motif_list:
                    raise ValueError(f"ERROR: Motif {motif} specified more than once in line: {line}")

                count = int(count)
                total_counts += count
                canonical_motif = compute_canonical_motif(motif)
                if canonical_motif not in canonical_motif_to_motif:
                    canonical_motif_to_motif[canonical_motif] = motif
                    canonical_motif_to_count[canonical_motif] = count
                    unique_motif_list.append(motif)
                else:
                    canonical_motif_to_count[canonical_motif] += count

            ori_to_efficient_motif_map = None
            if efficient_motif_data_available:
                if reference_region not in vamos_ori_to_efficient_motif_map:
                    raise ValueError(f"{reference_region} missing from Vamos motif efficient motif table")

                ori_to_efficient_motif_map = vamos_ori_to_efficient_motif_map[reference_region]

            efficient_motif_list = []
            for motif in unique_motif_list:
                canonical_motif = compute_canonical_motif(motif)
                count = canonical_motif_to_count[canonical_motif]
                if efficient_motif_data_available:
                    efficient_motif = ori_to_efficient_motif_map.get(motif)

                    if not efficient_motif:
                        raise ValueError(f"ERROR: No efficient motif found for {motif} in line: {line}")

                    value = f"{motif}:{efficient_motif}:{count}"

                    if efficient_motif not in efficient_motif_list:
                        efficient_motif_list.append(efficient_motif)
                else:
                    value = f"{motif}:{count}"

                motif_and_count_list.append(value)

            result[reference_region] = {
                "unique_motifs": ",".join(unique_motif_list),
                "efficient_motifs": ",".join(efficient_motif_list),
                "motif_frequencies": ",".join(motif_and_count_list),
                "num_unique_motifs": len(motif_and_count_list),
            }

    return result


def parse_vamos_efficient_motif_lookup(efficient_motif_tsv_path):
    """Parse table and return a lookup dictionary mapping the reference region interval to a dictionary
    that maps the original unique motif to its corresponding efficient motif.

    Example row:
        #1      chr     chr7
        #2      start   52646623
        #3      end     52646636
        #4      effMotifs       TA,TAA
        #5      effMotifs_counts        3020,317
        #6      status  OPTIMAL
        #7      effSet/oriSet   2/5
        #8      oriMotifs       TA,TAA,TG,AA,TAG
        #9      effMotif_counterparts   TA,TAA,TA,TA,TA
        #10     mapping_cost    0,0,1,1,1
        #11     effMotif_indicator      1,1,0,0,0
        #12     delta   333
    """
    print(f"Parsing efficient motifs from {efficient_motif_tsv_path}")
    efficient_motif_lookup = {}
    with gzip.open(efficient_motif_tsv_path, 'rt') as f:
        header = f.readline().strip().split('\t')
        col_indices = {col: i for i, col in enumerate(header)}
        for line in tqdm.tqdm(f, unit=" lines", unit_scale=True):
            fields = line.strip().split('\t')
            chrom = fields[col_indices['chr']]
            start_0based = int(fields[col_indices['start']])
            end_1based = int(fields[col_indices['end']])
            ori_motifs = fields[col_indices['oriMotifs']].split(",")
            efficient_motifs_q01 = fields[col_indices['effMotif_counterparts']].split(",")

            assert len(ori_motifs) == len(efficient_motifs_q01), line

            reference_region = f"{chrom}:{start_0based}-{end_1based}"
            efficient_motif_lookup[reference_region] = dict(zip(ori_motifs, efficient_motifs_q01))

    return efficient_motif_lookup


if __name__ == "__main__":
    main()