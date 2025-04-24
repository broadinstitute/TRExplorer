"""This script takes the annotated catalog JSON and loads it into a BigQuery table."""

import argparse
import collections
import datetime
import gzip
import json
import ijson
import os
from pprint import pformat
import re
import requests
from google.cloud import bigquery
import time
import tqdm
import sys

from str_analysis.utils.misc_utils import parse_interval
from str_analysis.utils.canonical_repeat_unit import compute_canonical_motif
from intervaltree import IntervalTree, Interval

PROJECT_ID = "cmg-analysis"
DATASET_ID = "tandem_repeat_explorer"
TABLE_ID = "catalog"

parser = argparse.ArgumentParser(description="Load data into BigQuery from the annotated catalog JSON file.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-c", "--catalog-path", help="Path to the annotated catalog JSON file",
                    default="https://github.com/broadinstitute/tandem-repeat-catalog/releases/download/v1.0/repeat_catalog_v1.hg38.1_to_1000bp_motifs.EH.with_annotations.json.gz")
parser.add_argument("-n", type=int, help="Number of records to read from the catalog")
parser.add_argument("-d", "--known-disease-associated-loci",
                    help="ExpansionHunter catalog with the latest known disease-associated loci",
                    default="https://raw.githubusercontent.com/broadinstitute/str-analysis/refs/heads/main/str_analysis/variant_catalogs/variant_catalog_without_offtargets.GRCh38.json")
parser.add_argument("--tenk10k-tsv", default="../data-prep/tenk10k_str_mt_rows.reformatted.tsv.gz")
parser.add_argument("--hprc100-tsv", default="../data-prep/runs_in_hprc.2025_04.grouped_by_locus_and_motif.tsv.gz")

args = parser.parse_args()

def get_json_iterator(content, is_gzipped=False):
    """Helper function to get an ijson iterator from content."""
    if is_gzipped:
        if isinstance(content, bytes):
            content = gzip.decompress(content)
        else:
            content = gzip.open(content, 'rb')

    return ijson.items(content, "item", use_float=True)

def parse_allele_histograms_from_tsv(tsv_path):
    """Parse the tenk10k TSV file line by line and create a lookup dictionary."""
    
    lookup = {}
    with gzip.open(tsv_path, 'rt') as f:
        header = f.readline().strip().split('\t')
        col_indices = {col: i for i, col in enumerate(header)}
        for line in tqdm.tqdm(f, unit=" lines", unit_scale=True):
            fields = line.strip().split('\t')
            locus_id = fields[col_indices['locus_id']]
            lookup[locus_id] = {
                'allele_size_histogram': fields[col_indices['allele_size_histogram']],
                'mode_allele': int(fields[col_indices['mode_allele']]),
                'stdev': float(fields[col_indices['stdev']]),
                #'mean': float(fields[col_indices['mean']]),
                'median': float(fields[col_indices['median']]),
                '99th_percentile': float(fields[col_indices['99th_percentile']]),
            }
    
    return lookup

print(f"Parsing tenk10k data from {args.tenk10k_tsv}")
tenk10k_lookup = parse_allele_histograms_from_tsv(args.tenk10k_tsv)
print(f"Parsed {len(tenk10k_lookup):,d} records from the tenk10k table: {args.tenk10k_tsv}")

print(f"Parsing hprc100 data from {args.hprc100_tsv}")
hprc100_lookup = parse_allele_histograms_from_tsv(args.hprc100_tsv)
print(f"Parsed {len(hprc100_lookup):,d} records from the hprc100 table: {args.hprc100_tsv}")



print(f"Parsing known disease-associated loci from {args.known_disease_associated_loci}")
if os.path.isfile(args.known_disease_associated_loci):
    fopen = open if args.known_disease_associated_loci.endswith("gz") else gzip.open 
    with fopen(args.known_disease_associated_loci) as f:
        known_disease_associated_loci = ijson.items(f, "item", use_float=True)
elif args.known_disease_associated_loci.startswith("http"):
    response = requests.get(args.known_disease_associated_loci)
    known_disease_associated_loci = ijson.items(response.content, "item", use_float=True)
else:
    parser.error(f"Invalid catalog path: {args.known_disease_associated_loci}")

known_disease_associated_loci = {
    x["MainReferenceRegion"].replace("chr", ""): x for x in known_disease_associated_loci if x.get("Diseases")
}
known_disease_associated_locus_ids = {x["LocusId"] for x in known_disease_associated_loci.values()}
print(f"Parsed {len(known_disease_associated_loci)} known disease-associated loci:")
print(", ".join(sorted(known_disease_associated_locus_ids)))

strchive_data_json = requests.get("https://raw.githubusercontent.com/dashnowlab/STRchive/refs/heads/main/data/STRchive-loci.json").json()
strchive_id_lookup = {}
for locus in strchive_data_json:
    strchive_gene = locus["gene"].upper()
    strchive_locus_id = locus["id"].upper()
    if strchive_gene in strchive_id_lookup:
        print(f"WARNING: {strchive_gene} already has an ID in strchive_gene_to_id: {strchive_id_lookup[strchive_gene]}")
    strchive_id_lookup[strchive_gene] = strchive_locus_id

strchive_id_lookup["ARX_1"] = strchive_id_lookup["ARX"]
strchive_id_lookup["ARX_2"] = strchive_id_lookup["ARX"]
strchive_id_lookup["HOXA13_1"] = strchive_id_lookup["HOXA13"]
strchive_id_lookup["HOXA13_2"] = strchive_id_lookup["HOXA13"]
strchive_id_lookup["HOXA13_3"] = strchive_id_lookup["HOXA13"]
del strchive_id_lookup["ARX"]
del strchive_id_lookup["HOXA13"]

stripy_id_lookup = {}
# Compute STRipy urls
for locus_id in known_disease_associated_locus_ids:
    stripy_name = locus_id
    stripy_url = f"https://stripy.org/database/{stripy_name}"
    r = requests.get(stripy_url)
    if r.ok and "invalid locus" not in r.content.decode("UTF-8").lower():
        stripy_id_lookup[locus_id] = stripy_name
    else:
        print(f"WARNING: STRipy page not found for {locus_id}")

for locus_id in set(strchive_id_lookup.keys()) - {x["LocusId"] for x in known_disease_associated_loci.values()}:
    print(f"WARNING: {locus_id} is in STRchive but not in {args.known_disease_associated_loci}")

known_disease_associated_loci_interval_trees = collections.defaultdict(IntervalTree)
for reference_region, locus_info in known_disease_associated_loci.items():
    if locus_info["LocusId"] in strchive_id_lookup:
        locus_info["STRchiveId"] = strchive_id_lookup[locus_info["LocusId"]]
    else:
        print(f"WARNING: {locus_info['LocusId']} is not in STRchive")

    if locus_info["LocusId"] in stripy_id_lookup:
        locus_info["STRipyId"] = stripy_id_lookup[locus_info["LocusId"]]
    else:
        print(f"WARNING: {locus_info['LocusId']} is not in STRipy")

    chrom, start_0based, end_1based = parse_interval(reference_region)
    known_disease_associated_loci_interval_trees[chrom].addi(start_0based, end_1based, data=locus_info)

    # drop all keys from  locus_info except LocusId, STRchiveId, and Diseases
    known_disease_associated_loci[reference_region] = {
        key: locus_info[key] for key in ["LocusId", "STRchiveId", "STRipyId", "Diseases", "RepeatUnit"] if key in locus_info
    }

# Define the schema
schema = [
    bigquery.SchemaField("id", "INTEGER", mode="REQUIRED"),
    bigquery.SchemaField("chrom_index", "INTEGER", mode="REQUIRED"),
    bigquery.SchemaField("chrom", "STRING", mode="REQUIRED"),
    bigquery.SchemaField("start_0based", "INTEGER", mode="REQUIRED"),
    bigquery.SchemaField("end_1based", "INTEGER", mode="REQUIRED"),
    bigquery.SchemaField("ReferenceRegion", "STRING", mode="REQUIRED"),
    bigquery.SchemaField("LocusId", "STRING"),
    bigquery.SchemaField("MotifSize", "INTEGER"),
    bigquery.SchemaField("ReferenceMotif", "STRING"),
    bigquery.SchemaField("CanonicalMotif", "STRING"),
    bigquery.SchemaField("NumRepeatsInReference", "INTEGER"),
    bigquery.SchemaField("ReferenceRepeatPurity", "FLOAT"),
    bigquery.SchemaField("NsInFlanks", "INTEGER"),
    bigquery.SchemaField("TRsInRegion", "INTEGER"),
    bigquery.SchemaField("Source", "STRING"),

    bigquery.SchemaField("FoundInKnownDiseaseAssociatedLoci", "STRING"),
    bigquery.SchemaField("FoundInIllumina174kPolymorphicTRs", "STRING"),
    bigquery.SchemaField("FoundInPerfectRepeatsInReference", "STRING"),
    bigquery.SchemaField("FoundInPolymorphicTRsInT2TAssemblies", "STRING"),

    bigquery.SchemaField("LeftFlankMappability", "FLOAT"),
    bigquery.SchemaField("RightFlankMappability", "FLOAT"),
    bigquery.SchemaField("FlanksAndLocusMappability", "FLOAT"),
    bigquery.SchemaField("VariationCluster", "STRING"),
    bigquery.SchemaField("VariationClusterSizeDiff", "INTEGER"),

    bigquery.SchemaField("KnownDiseaseAssociatedLocus", "STRING"),
    bigquery.SchemaField("KnownDiseaseAssociatedMotif", "STRING"),
    bigquery.SchemaField("DiseaseInfo", "STRING"),

    bigquery.SchemaField("GencodeGeneRegion", "STRING"),
    bigquery.SchemaField("GencodeGeneName", "STRING"),
    bigquery.SchemaField("GencodeGeneId", "STRING"),
    bigquery.SchemaField("GencodeTranscriptId", "STRING"),
    bigquery.SchemaField("RefseqGeneRegion", "STRING"),
    bigquery.SchemaField("RefseqGeneName", "STRING"),
    bigquery.SchemaField("RefseqGeneId", "STRING"),
    bigquery.SchemaField("RefseqTranscriptId", "STRING"),
    bigquery.SchemaField("ManeGeneRegion", "STRING"),
    bigquery.SchemaField("ManeGeneName", "STRING"),
    bigquery.SchemaField("ManeGeneId", "STRING"),
    bigquery.SchemaField("ManeTranscriptId", "STRING"),
    bigquery.SchemaField("LPSLengthStdevFromHPRC100", "FLOAT"),
    bigquery.SchemaField("LPSMotifFromHPRC100", "STRING"),
    bigquery.SchemaField("LPSMotifFractionFromHPRC100", "FLOAT"),
    bigquery.SchemaField("LPSMotifDenominatorFromHPRC100", "INTEGER"),
    bigquery.SchemaField("AlleleFrequenciesFromIllumina174k", "STRING"),
    bigquery.SchemaField("StdevFromIllumina174k", "FLOAT"),
    bigquery.SchemaField("AlleleFrequenciesFromT2TAssemblies", "STRING"),
    bigquery.SchemaField("StdevFromT2TAssemblies", "FLOAT"),

    bigquery.SchemaField("TenK10K_AlleleHistogram", "STRING"),
    bigquery.SchemaField("TenK10K_ModeAllele", "INTEGER"),
    bigquery.SchemaField("TenK10K_Stdev", "FLOAT"),
    bigquery.SchemaField("TenK10K_Median", "INTEGER"),
    bigquery.SchemaField("TenK10K_99thPercentile", "INTEGER"),

    bigquery.SchemaField("HPRC100_AlleleHistogram", "STRING"),
    bigquery.SchemaField("HPRC100_ModeAllele", "INTEGER"),
    bigquery.SchemaField("HPRC100_Stdev", "FLOAT"),
    bigquery.SchemaField("HPRC100_Median", "INTEGER"),
    bigquery.SchemaField("HPRC100_99thPercentile", "INTEGER"),

]

field_names = {field.name for field in schema}

# Initialize BigQuery client
client = bigquery.Client(project=PROJECT_ID)

def does_table_exist(table_ref):
    try:
        client.get_table(table_ref)
        return True
    except Exception:
        return False

def insert_with_retries(table_ref, rows_to_insert, batch_size=1000, max_retries=5):
    for i in range(0, len(rows_to_insert), batch_size):
        batch = rows_to_insert[i:i+batch_size]
        retries = 0
        while retries < max_retries:
            try:
                errors = client.insert_rows_json(table_ref, batch)
                if errors:
                    print(f"Encountered errors while inserting rows: {errors}")
                    raise Exception(f"Encountered errors while inserting rows: {errors}")
                break
            except Exception as e:
                print(f"Error inserting batch: {e}. Retrying...")
                retries += 1
                time.sleep(5)
        if retries == max_retries:
            raise Exception(f"Failed to insert batch after {max_retries} retries")

# Create dataset if it doesn't exist
dataset_ref = client.dataset(DATASET_ID)
try:
    client.get_dataset(dataset_ref)
    print(f"BigQuery dataset {DATASET_ID} already exists")
except Exception:
    dataset = bigquery.Dataset(dataset_ref)
    dataset.location = "US-CENTRAL1"
    client.create_dataset(dataset)
    print(f"Created dataset {DATASET_ID}")


new_table_id = f"{TABLE_ID}_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}"
# Check if table exists
new_table_ref = dataset_ref.table(new_table_id)

# Delete table if it exists
if does_table_exist(new_table_ref):
    print(f"ERROR: Table {new_table_id} already exists")
    sys.exit(1)

    
# Create table if it doesn't exist
new_table = bigquery.Table(new_table_ref, schema=schema)
client.create_table(new_table)
print(f"Created table {new_table_id}")

# Read catalog data
print(f"Reading catalog from {args.catalog_path}")
is_gzipped = args.catalog_path.endswith('gz')
if args.catalog_path.startswith("http"):
    response = requests.get(args.catalog_path)
    catalog = get_json_iterator(response.content, is_gzipped)
elif os.path.isfile(args.catalog_path):
    catalog = get_json_iterator(args.catalog_path, is_gzipped)
else:
    parser.error(f"Invalid catalog path: {args.catalog_path}")


chrom_indices = {str(i): i for i in range(1, 23)}
chrom_indices.update({"X": 23, "Y": 24, "M": 25, "MT": 25})

# Prepare data for loading
rows_to_insert = []
batch_size = 1000
locus_ids_with_added_disease_info = set()

counters = collections.Counter()
for i, record in tqdm.tqdm(enumerate(catalog), unit=" records", unit_scale=True):
    if args.n and i >= args.n:
        break

    chrom, start_0based, end_1based = parse_interval(record["ReferenceRegion"])
    record["chrom"] = chrom.replace("chr", "").upper()
    record["chrom_index"] = chrom_indices[record["chrom"]]
    record["start_0based"] = start_0based
    record["end_1based"] = end_1based
    record["MotifSize"] = len(record["CanonicalMotif"])
    if record.get("LPSMotifFractionFromHPRC100"):
        counters["rows_with_original_LPS_data_from_hprc100"] += 1
        tokens = record["LPSMotifFractionFromHPRC100"].split(": ")  # example value: "AACCCT: 112/187"
        record["LPSMotifFromHPRC100"] = tokens[0]
        numerator, denominator = map(int, tokens[1].split("/"))
        record["LPSMotifFractionFromHPRC100"] = numerator/denominator if denominator > 0 else None
        record["LPSMotifDenominatorFromHPRC100"] = denominator

    if record.get("StdevFromT2TAssemblies"):
        record["StdevFromT2TAssemblies"] = round(record["StdevFromT2TAssemblies"], 3)

    if record.get("StdevFromIllumina174k"):
        record["StdevFromIllumina174k"] = round(record["StdevFromIllumina174k"], 3)

    motif_match = re.match(r"^[(]([A-Z]+)[)][+*]", record["LocusStructure"])
    record["ReferenceMotif"] = motif_match.group(1) if motif_match else None

    catalog_reference_region = record["ReferenceRegion"].replace("chr", "")
    if catalog_reference_region in known_disease_associated_loci:
        known_locus_info = known_disease_associated_loci[catalog_reference_region]
        locus_ids_with_added_disease_info.add(known_locus_info["LocusId"])
        record["DiseaseInfo"] = json.dumps(known_locus_info)
        #print(f"Added disease info for locus #{len(locus_ids_with_added_disease_info)}:", known_locus_info["LocusId"])
    else:
        chrom, start_0based, end_1based = parse_interval(catalog_reference_region)
        for overlapping_interval in known_disease_associated_loci_interval_trees[chrom].overlap(start_0based, end_1based):
            known_locus_info = overlapping_interval.data
            known_locus_canonical_motif = compute_canonical_motif(known_locus_info["RepeatUnit"])
            catalog_record_canonical_motif = compute_canonical_motif(record["ReferenceMotif"])

            if known_locus_canonical_motif == catalog_record_canonical_motif and overlapping_interval.overlap_size(start_0based, end_1based) > len(catalog_record_canonical_motif):
                record["DiseaseInfo"] = json.dumps(known_locus_info)
                locus_ids_with_added_disease_info.add(known_locus_info["LocusId"])
                print()
                print(f"Added disease info for overlapping locus #{len(locus_ids_with_added_disease_info)}:", known_locus_info["LocusId"])
                print(f"    Known locus definition: {known_locus_info['ReferenceRegion']} {known_locus_info['RepeatUnit']}")
                print(f"  Catalog locus definition: {record['ReferenceRegion']} {record['ReferenceMotif']}")
                break

    if record.get("StdevFromIllumina174k") and not record.get("AlleleFrequenciesFromIllumina174k"):
        print("WARNING: StdevFromIllumina174k is present but AlleleFrequenciesFromIllumina174k is missing in record:", pformat(record, indent=4))
    if record.get("LPSLengthStdevFromHPRC100") and not record.get("LPSMotifFromHPRC100"):
        print("WARNING: LPSLengthStdevFromHPRC100 is present but LPSMotifFromHPRC100 is missing in record:", pformat(record, indent=4))

    #if record.get("StdevFromT2TAssemblies") and not record.get("AlleleFrequenciesFromT2TAssemblies"):
    #    This happens where loci in source catalogs overlapped
    #    print("WARNING: StdevFromT2TAssemblies is present but AlleleFrequenciesFromT2TAssemblies is missing in record:", pformat(record, indent=4))

    # Add id field
    record["id"] = i + 1

    if record["LocusId"] in tenk10k_lookup:
        counters["rows_with_tenk10k_data"] += 1
        record["TenK10K_AlleleHistogram"] = tenk10k_lookup[record["LocusId"]]["allele_size_histogram"]
        record["TenK10K_ModeAllele"] = tenk10k_lookup[record["LocusId"]]["mode_allele"]
        record["TenK10K_Stdev"] = tenk10k_lookup[record["LocusId"]]["stdev"]
        record["TenK10K_Median"] = tenk10k_lookup[record["LocusId"]]["median"]
        record["TenK10K_99thPercentile"] = tenk10k_lookup[record["LocusId"]]["99th_percentile"]

    if record["LocusId"] in hprc100_lookup:
        counters["rows_with_hprc100_data"] += 1
        record["HPRC100_AlleleHistogram"] = hprc100_lookup[record["LocusId"]]["allele_size_histogram"]
        record["HPRC100_ModeAllele"] = hprc100_lookup[record["LocusId"]]["mode_allele"]
        record["HPRC100_Stdev"] = hprc100_lookup[record["LocusId"]]["stdev"]
        record["HPRC100_Median"] = hprc100_lookup[record["LocusId"]]["median"]
        record["HPRC100_99thPercentile"] = hprc100_lookup[record["LocusId"]]["99th_percentile"]

    counters["total_rows"] += 1
    # Convert any None values to None (BigQuery will handle NULL)
    row = {k: v for k, v in record.items() if k in field_names}
    rows_to_insert.append(row)
    if len(rows_to_insert) >= batch_size:
        insert_with_retries(new_table_ref, rows_to_insert)
        rows_to_insert = []

if rows_to_insert:
    insert_with_retries(new_table_ref, rows_to_insert)

if len(locus_ids_with_added_disease_info) != len(known_disease_associated_locus_ids):
    print(f"WARNING: {len(known_disease_associated_loci) - len(locus_ids_with_added_disease_info)} out of {len(known_disease_associated_loci)} known disease-associated loci were not found in the catalog. "
           "Missing LocusIds:", ", ".join(known_disease_associated_locus_ids - locus_ids_with_added_disease_info))

with open("../index.html", "r") as f:
    html_content = f.read()

html_content = re.sub(
    r"const TABLE_ID[\s]*=[\s]*'catalog[^']*'",
    f"const TABLE_ID = '{new_table_id}'", 
    html_content)

with open("../index.html", "wt") as f:
    f.write(html_content)

print("Done!")

print("\nCounters:")
for key, count in sorted(counters.items(), key=lambda x: x[1], reverse=True):
    print(f"{count:10,d} {key}")
