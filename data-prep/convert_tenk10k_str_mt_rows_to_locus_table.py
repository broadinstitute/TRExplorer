"""
This script takes a TSV file. This is an example row:

$1                      locus : chr1:100004721
$2                    alleles : ["A","<STR19>","<STR78>","<STR123>","<STR136>","<STR150>","<STR151>","<STR170>"]
$3                       rsid : NA
$4                       qual : -1.0000e+01
$5                    filters : NA
$6                       info : {"END":100004741,"REF":20,"REPID":"1-100004722-100004741-T","RL":20,"RU":"T","SVTYPE":null,"VARID":"1-100004722-100004741-T"}
$7                 variant_qc : {"AC":[4766,51,1,1,1,2,1,1],"AF":[0.9879767827529021,0.010572139303482588,2.0729684908789387E-4,2.0729684908789387E-4,2.0729684908789387E-4,4.1459369817578774E-4,2.0729684908789387E-4,2.0729684908789387E-4],"AN":4824,"homozygote_count":[2354,0,0,0,0,0,0,0],"call_rate":1.0,"n_called":2412,"n_not_called":0,"n_filtered":0,"n_het":58,"n_non_ref":58,"het_freq_hwe":null,"p_value_hwe":null,"p_value_excess_het":null}
$8                      REPID : 1-100004722-100004741-T
$9         rep_length_alleles : [20,19,78,123,136,150,151,170]
$10              motif_length : 1
$11         bp_length_alleles : [20,19,78,123,136,150,151,170]
$12           aggregated_info : {"allele_array_counts":[{"key":170,"value":1},{"key":20,"value":4766},{"key":78,"value":1},{"key":123,"value":1},{"key":150,"value":2},{"key":151,"value":1},{"key":19,"value":51},{"key":136,"value":1}],"mode_allele":20}
$13               num_alleles : 8
$14   sum_alleles_is_not_mode : 58
$15  prop_alleles_is_not_mode : 1.2023e-02
$16                binom_hwep : 8.9373e-01
$17                   obs_het : 2.4046e-02
$18                variant_lc : 4.4575e+01

and outputs a table with columns:

chrom                 (example: "1")
start_0based          (example: 100004721)
end_0based            (example: 100004741)
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

parser = argparse.ArgumentParser()
parser.add_argument("-n", type=int, default=None, help="Number of lines to process")
parser.add_argument("--tenk10k-tsv", default="tenk10k_str_mt_rows.tsv.bgz")
parser.add_argument("--trexplorer-catalog", default="~/code/tandem-repeat-catalogs/results__2024-10-01/release_draft_2024-10-01/repeat_catalog_v1.hg38.1_to_1000bp_motifs.EH.json.gz")
parser.add_argument("--output-tsv", default="tenk10k_str_mt_rows.reformatted.tsv.gz")
args = parser.parse_args()



N = args.n
input_tsv = args.tenk10k_tsv
output_tsv = args.output_tsv
catalog_path = os.path.expanduser(args.trexplorer_catalog)

f = gzip.open(input_tsv, "rt")
header_line = f.readline().strip().split("\t")

fout = gzip.open(output_tsv, "wt")
fout.write("\t".join([
    "chrom", 
    "start_0based", 
    "end_1based", 
    "locus_id", 
    "motif", 
    "stdev",
    "allele_size_histogram", 
    "mode_allele",
]) + "\n")

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
        catalog_motif = item["LocusStructure"].strip("()*+")

        catalog_intervals[catalog_chrom].addi(catalog_start_0based, catalog_end_1based, data=catalog_motif)

print(f"Parsed {counters['catalog entries']:9,d} entries from {os.path.basename(catalog_path)}")

for i, line in tqdm.tqdm(enumerate(f), unit=" lines", unit_scale=True):
    if N is not None and i > N:
        break
    fields = line.strip().split("\t")
    # zip(header_fields, fields)
    if len(fields) != 18:
        print(f"WARNING: Skipping line #{i+1} because it has {len(fields)} fields")
        continue
    
    chrom, start_0based = fields[0].split(":")
    start_0based = int(start_0based)

    info_field = fields[5]
    info_json = json.loads(info_field)
    end_1based = int(info_json["END"])
    motif = info_json["RU"]

    locus_id = f"{chrom}-{start_0based}-{end_1based}-{motif}"
    aggregated_info = fields[11]
    aggregated_info_json = json.loads(aggregated_info)
    allele_array_counts = aggregated_info_json["allele_array_counts"]
    mode_allele = aggregated_info_json["mode_allele"]

    found_in_catalog = None
    if locus_id in catalog_locus_ids:
        found_in_catalog = "same"
    else:
        # if it exists, retrieve the largest overlapping interval in the catalog with the same canonical motif
        canonical_motif = compute_canonical_motif(motif)
        found_interval = None
        for interval in catalog_intervals[chrom].overlap(start_0based, end_1based):
            catalog_canonical_motif = compute_canonical_motif(interval.data)
            if catalog_canonical_motif == canonical_motif and (
                found_interval is None or found_interval.length() < interval.length()):
                found_interval = interval

        if found_interval is not None:
            size_diff = abs((end_1based - start_0based) - found_interval.length())
            if size_diff < len(canonical_motif):
                found_in_catalog = "same"
            else:
                found_in_catalog = "overlaps"

    stdev = ""
    allele_array_counts = ""
    if found_in_catalog is None:
        counters["not_found_in_catalog"] += 1
    else:
        stdev = np.std([d["key"] for d in allele_array_counts for _ in range(d["value"])])
        stdev = f"{stdev:.3f}"
        if found_in_catalog == "same":
            counters["same"] += 1
            allele_array_counts = [f"{d['key']}x:{d['value']}" for d in allele_array_counts]
            allele_array_counts = ",".join(allele_array_counts)
        elif found_in_catalog == "overlaps":
            counters["overlaps"] += 1
        else:
            raise ValueError(f"Unexpected value for found_in_catalog: {found_in_catalog}")

    counters["total"] += 1

    output_row = [chrom, start_0based, end_1based, locus_id, motif, stdev, allele_array_counts, mode_allele]
    fout.write("\t".join(map(str, output_row)) + "\n")

print(f"Wrote {counters['total']:9,d} lines to {output_tsv}")

# print counters
for key, count in sorted(counters.items(), key=lambda x: x[0]):
    print(f"{count:10,d} {key}")
