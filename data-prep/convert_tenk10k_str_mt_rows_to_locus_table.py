"""
This script takes a TSV file. This is an example row:

$1                      locus : chr1:100004721
$2                    alleles : ["A","<STR19>","<STR78>","<STR123>","<STR136>","<STR150>","<STR151>","<STR170>"]
$3                       rsid : NA
$4                       qual : -1.0000e+01
$5                    filters : NA
$6                       info : {"END":100004741,"REF":20,"REPID":"1-100004722-100004741-T","RL":20,"RU":"T","SVTYPE":null,"VARID":"1-100004722-100004741-T"}
$7                 variant_qc : {"AC":[3801,43,1,1,1,2,0,1],"AF":[0.9872727272727273,0.011168831168831168,2.5974025974025974E-4,2.5974025974025974E-4,2.5974025974025974E-4,5.194805194805195E-4,0.0,2.5974025974025974E-4],"AN":3850,"homozygote_count":[1876,0,0,0,0,0,0,0],"call_rate":1.0,"n_called":1925,"n_not_called":0,"n_filtered":0,"n_het":49,"n_non_ref":49,"het_freq_hwe":null,"p_value_hwe":null,"p_value_excess_het":null}
$8                      REPID : 1-100004722-100004741-T
$9         rep_length_alleles : [20,19,78,123,136,150,151,170]
$10              motif_length : 1
$11         bp_length_alleles : [20,19,78,123,136,150,151,170]
$12   aggregated_info_diploid : {"diploid_allele_array_counts":[{"key":"20/150","value":2},{"key":"20/78","value":1},{"key":"20/20","value":1876},{"key":"20/136","value":1},{"key":"20/170","value":1},{"key":"20/19","value":43},{"key":"20/123","value":1}]}
$13           aggregated_info : {"allele_array_counts":[{"key":170,"value":1},{"key":20,"value":3801},{"key":78,"value":1},{"key":123,"value":1},{"key":150,"value":2},{"key":19,"value":43},{"key":136,"value":1}],"mode_allele":20}
$14               num_alleles : 7
$15   sum_alleles_is_not_mode : 49
$16  prop_alleles_is_not_mode : 1.2727e-02
$17                binom_hwep : 8.8440e-01
$18                   obs_het : 2.5455e-02
$19                variant_lc : 4.4834e+01

and outputs a table with columns:

chrom                 (example: "1")
start_0based          (example: 100004721)
end_1based            (example: 100004741)
locus_id              (example: "1-100004721-100004741-T")
motif                 (example: "T")
stdev                 (example: 1.2023e-02)
mean                  (example: 100)
median                (example: 20)
99th_percentile       (example: 170)
allele_size_histogram (example: "19x:51,20x:4766,78x:1,123x:1,136x:1,150x:2,151x:1,170x:1")
"""

import argparse
import collections
import gzip
import ijson
import intervaltree
import json
import numpy as np
import os
import tqdm

from str_analysis.utils.canonical_repeat_unit import compute_canonical_motif
from str_analysis.utils.misc_utils import parse_interval


SAME_AS_IN_CATALOG_LABEL = "the same as"
ALMOST_SAME_AS_IN_CATALOG_LABEL = "almost the same as"
OVERLAP_CATALOG_LABEL = "overlaps"

OUTPUT_HEADER_FIELDS = [
    "chrom", 
    "start_0based", 
    "end_1based", 
    "locus_id", 
    "tenk10k_locus_id", 
    "tenk10k_interval",
    "motif", 
    "allele_size_histogram", 
    "biallelic_histogram",
    "mode_allele",
    "mean",
    "stdev",
    "median",
    "99th_percentile",
    "tenk_10k_vs_catalog_overlap_size",
    "tenk_10k_vs_catalog_size_diff",
    "binom_hwep",
    "variant_lc",
]

def write_to_output(output_row_data, output_tsv, counters):
    fout = gzip.open(output_tsv, "wt")
    fout.write("\t".join(OUTPUT_HEADER_FIELDS) + "\n")

    total = 0

    sorted_output_rows = sorted(output_row_data.values(), key=lambda x: (x['chrom'], x['start_0based'], x['end_1based']))
    for output_row in tqdm.tqdm(sorted_output_rows, unit=" output rows", unit_scale=True):
        found_in_catalog = output_row["found_in_catalog"]
        tenk10k_locus_id = output_row["tenk10k_locus_id"]

        allele_array_counts = [x for x in output_row["allele_array_counts"] if x["key"] is not None]
        allele_array_counts = list(sorted(allele_array_counts, key=lambda x: x["key"]))

        # Extract values and weights for numpy operations
        allele_sizes = np.array([d["key"] for d in allele_array_counts])
        allele_counts = np.array([d["value"] for d in allele_array_counts])

        # Compute mode (most common value)
        mode_idx = np.argmax(allele_counts)
        recomputed_mode_allele = allele_sizes[mode_idx]

        if output_row["mode_allele"] is not None and recomputed_mode_allele != output_row["mode_allele"]:
            print(f"WARNING: Recomputed mode allele = {recomputed_mode_allele} does not match the original mode allele {output_row['mode_allele']} for {tenk10k_locus_id}")
        
        # compute different stats
        output_row["mode_allele"] = int(recomputed_mode_allele)
        np_repeat_result = np.repeat(allele_sizes, allele_counts)
        output_row["median"] = int(np.median(np_repeat_result))
        output_row["99th_percentile"] = int(np.percentile(np_repeat_result, 99))
        mean = np.average(np_repeat_result)
        variance = np.var(np_repeat_result)
        output_row["stdev"] = f"{np.sqrt(variance):.3f}"
        output_row["mean"] = f"{mean:.3f}"
        output_row["min_allele"] = int(np.min(allele_sizes))
        output_row["max_allele"] = int(np.max(allele_sizes))
        output_row["unique_alleles"] = len(set(allele_sizes))
        output_row["num_called_alleles"] = sum(allele_counts)

        output_row["tenk_10k_vs_catalog_overlap_size"] = output_row["found_interval_overlap_size"]
        output_row["tenk_10k_vs_catalog_size_diff"] = output_row["found_interval_size_diff"]

        output_row["allele_size_histogram"] = ""
        output_row["biallelic_histogram"] = ""
        if found_in_catalog is None:
            counters["tenk10k rows not found in catalog"] += 1
        else:
            counters[f"tenk10k rows were {found_in_catalog} catalog entry"] += 1        
            if found_in_catalog in (SAME_AS_IN_CATALOG_LABEL, ALMOST_SAME_AS_IN_CATALOG_LABEL):
                output_row["allele_size_histogram"] = ",".join([f"{key}x:{value}" for key, value in zip(allele_sizes, allele_counts)])
                try:
                    # convert {"diploid_allele_array_counts":[{"key":"20/150","value":2}, ...]} to  biallelic_histogram : 6/6:4,6/12:1,12/15:1,15/15:2,
                    biallelic_histogram = []
                    for diploid_allele_array_count in output_row["diploid_allele_array_counts"]:
                        key = diploid_allele_array_count["key"]
                        if key is None:
                            continue

                        allele_1_allele_2 = key.split("/")
                        allele_1 = int(allele_1_allele_2[0])
                        allele_2 = int(allele_1_allele_2[-1])
                        value = int(diploid_allele_array_count["value"])

                        biallelic_histogram.append(f"{allele_1}/{allele_2}:{value}")
                    output_row["biallelic_histogram"] = ",".join(biallelic_histogram)
                except Exception as e:
                    print(f"WARNING: Unable to parse diploid_allele_array_counts: {output_row['diploid_allele_array_counts']}  {e}")

        output_fields = [output_row[field] for field in OUTPUT_HEADER_FIELDS]
        fout.write("\t".join(map(str, output_fields)) + "\n")
        total += 1

    fout.close()

    print(f"Wrote {total:9,d} lines to {output_tsv}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", type=int, default=None, help="Number of lines to process")
    parser.add_argument("--tenk10k-tsv", default="str_polymorphic_run_mt_bioheart_tob_v1_n1925_rows.tsv.bgz")
    parser.add_argument("--trexplorer-catalog", default="~/code/tandem-repeat-catalogs/results__2024-10-01/release_draft_2024-10-01/repeat_catalog_v1.hg38.1_to_1000bp_motifs.EH.json.gz")
    parser.add_argument("--output-tsv", default="tenk10k_str_mt_rows.reformatted.tsv.gz")
    args = parser.parse_args()



    N = args.n
    input_tsv = args.tenk10k_tsv
    output_tsv = args.output_tsv
    catalog_path = os.path.expanduser(args.trexplorer_catalog)


    # parse catalog
    """
    {
        "LocusId": "1-14069-14081-CCTC",
        "VariantType": "Repeat",
        "ReferenceRegion": "chr1:14069-14081",
        "LocusStructure": "(CCTC)*"
    }
    """

    counters = collections.defaultdict(int)
    catalog_locus_ids = set()
    catalog_intervals = collections.defaultdict(intervaltree.IntervalTree)

    print(f"Parsing catalog: {os.path.basename(catalog_path)}")
    with gzip.open(catalog_path, "rt") as f_catalog:
        catalog_records = json.load(f_catalog)
        for item in tqdm.tqdm(catalog_records, unit=" items", unit_scale=True):
            catalog_locus_id = item["LocusId"]
            catalog_locus_ids.add(catalog_locus_id)

            counters["catalog entries"] += 1
            reference_region = item["ReferenceRegion"]
            if not isinstance(reference_region, str):
                print(f"WARNING: Skipping {catalog_locus_id} in {os.path.basename(catalog_path)} because reference_region is not a string: {reference_region}")
                continue

            catalog_chrom, catalog_start_0based, catalog_end_1based = parse_interval(reference_region)
            catalog_chrom = catalog_chrom.replace("chr", "")
            catalog_motif = item["LocusStructure"].strip("()*+")

            catalog_intervals[catalog_chrom].addi(catalog_start_0based, catalog_end_1based, data=(catalog_locus_id, catalog_motif))

    print(f"Parsed {counters['catalog entries']:9,d} entries from {os.path.basename(catalog_path)}")

    f = gzip.open(input_tsv, "rt")
    header_line = f.readline().strip().split("\t")

    #fout_bed = open(output_tsv.replace(".tsv.gz", ".bed"), "wt")

    print(f"Parsing {os.path.basename(input_tsv)}")
    processed_locus_ids = set()
    output_row_data = {} 
    for i, line in tqdm.tqdm(enumerate(f), unit=" lines", unit_scale=True, total=3_669_907):
        if N is not None and i > N:
            break

        counters["total tenk10k rows"] += 1
        fields = line.strip().split("\t")
        if len(fields) != len(header_line):
            print(f"WARNING: Skipping line #{i+1} because it has {len(fields)} fields")
            continue

        row_data = dict(zip(header_line, fields))
        
        chrom, start_0based = row_data["locus"].split(":")
        chrom = chrom.replace("chr", "")
        start_0based = int(start_0based)

        info_field = row_data["info"]
        info_json = json.loads(info_field)
        end_1based = int(info_json["END"])
        motif = info_json["RU"]
        
        locus_id = f"{chrom}-{start_0based}-{end_1based}-{motif}"
        if locus_id in processed_locus_ids:
            counters["tenk10k rows skipped because the same locus id was already processed"] += 1
            continue

        processed_locus_ids.add(locus_id)

        #fout_bed.write(f"{chrom}\t{start_0based}\t{end_1based}\t{motif}\t.\t.\n")

        aggregated_info = row_data["aggregated_info"]
        aggregated_info_json = json.loads(aggregated_info)

        aggregated_info_diploid = row_data["aggregated_info_diploid"]
        aggregated_info_diploid_json = json.loads(aggregated_info_diploid)

        found_in_catalog = None
        catalog_locus_id = None
        found_interval_size_diff = None
        found_interval_overlap_size = None
        if locus_id in catalog_locus_ids:
            found_in_catalog = SAME_AS_IN_CATALOG_LABEL
            catalog_locus_id = locus_id
            found_interval_size_diff = 0
            found_interval_overlap_size = end_1based - start_0based
        else:
            # if it exists, retrieve the largest overlapping interval in the catalog with the same canonical motif
            found_interval = None
            canonical_motif = None
            overlapping_intervals = catalog_intervals[chrom].overlap(start_0based, end_1based)
            if len(overlapping_intervals) > 0:
                canonical_motif = compute_canonical_motif(motif)
            for interval in overlapping_intervals:
                # make sure the overlap size is at least one repeat unit
                overlap_length = min(end_1based, interval.end) - max(start_0based, interval.begin)
                if overlap_length <= 0:
                    raise ValueError(f"overlap_length is not positive: {overlap_length}")
                    
                if overlap_length < len(canonical_motif):
                    continue

                size_diff = abs((end_1based - start_0based) - interval.length())
                if (found_interval is None 
                    or size_diff < found_interval_size_diff
                    or (
                        size_diff == found_interval_size_diff 
                        and overlap_length > found_interval_overlap_size
                    )
                ):
                    _, catalog_motif = interval.data
                    catalog_canonical_motif = compute_canonical_motif(catalog_motif)
                    if catalog_canonical_motif == canonical_motif:
                        found_interval = interval
                        found_interval_overlap_size = overlap_length
                        found_interval_size_diff = size_diff

            if found_interval is not None:
                catalog_locus_id, _ = found_interval.data
                if found_interval_size_diff == 0:
                    found_in_catalog = SAME_AS_IN_CATALOG_LABEL
                elif found_interval_size_diff < len(canonical_motif):
                    found_in_catalog = ALMOST_SAME_AS_IN_CATALOG_LABEL
                else:
                    found_in_catalog = OVERLAP_CATALOG_LABEL

        if not found_in_catalog:
            counters["tenk10k rows were not found in catalog"] += 1
            continue

        output_row = {
            "chrom": chrom,
            "start_0based": start_0based,
            "end_1based": end_1based,
            "locus_id": locus_id,
            "motif": motif,
            "locus_id": catalog_locus_id,
            "tenk10k_locus_id": info_json["REPID"],
            "tenk10k_interval": f"{chrom}:{start_0based}-{end_1based}",
            "allele_array_counts": aggregated_info_json["allele_array_counts"],
            "diploid_allele_array_counts": aggregated_info_diploid_json["diploid_allele_array_counts"],
            "mode_allele": aggregated_info_json["mode_allele"],
            "binom_hwep": row_data["binom_hwep"],
            "variant_lc": row_data["variant_lc"],

            "found_in_catalog": found_in_catalog,
            "found_interval_overlap_size": found_interval_overlap_size,
            "found_interval_size_diff": found_interval_size_diff,
        }

        if catalog_locus_id not in output_row_data:
            output_row_data[catalog_locus_id] = output_row
        else:
            prev_found_interval_overlap_size = output_row_data[catalog_locus_id]["found_interval_overlap_size"]
            prev_found_interval_size_diff = output_row_data[catalog_locus_id]["found_interval_size_diff"]
            if found_interval_size_diff < prev_found_interval_size_diff or (
                found_interval_size_diff == prev_found_interval_size_diff
                and found_interval_overlap_size > prev_found_interval_overlap_size
            ):
                output_row_data[catalog_locus_id] = output_row

    # fout_bed.close()

    filtered_output_row_data = {}
    for catalog_locus_id, output_row in output_row_data.items():
        if output_row["found_interval_size_diff"] < 3 * len(output_row["motif"]):
            filtered_output_row_data[catalog_locus_id] = output_row
    print(f"Filtered out {len(output_row_data) - len(filtered_output_row_data):,d} out of {len(output_row_data):,d} ({100 * (len(output_row_data) - len(filtered_output_row_data)) / len(output_row_data) if len(output_row_data) > 0 else 0:.1f}%) catalog entries because the interval size difference 3 repeat units or more")

    print(f"Writing tenk10k data for {len(filtered_output_row_data):,d} catalog entries to {output_tsv}")
    write_to_output(filtered_output_row_data, output_tsv, counters)

    # print counters
    print("Counters:")
    for key, count in sorted(counters.items(), key=lambda x: -x[1]):
        print(f"{count:10,d} {key}")

if __name__ == "__main__":
    main()
