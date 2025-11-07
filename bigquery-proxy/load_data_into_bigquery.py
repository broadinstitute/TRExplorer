"""This script takes the annotated catalog JSON and loads it into a BigQuery table."""

import argparse
import collections
import datetime
import gzip
import json
from typing import Any
import ijson
import os
from pprint import pformat
import re
import requests
from google.cloud import bigquery
import time
import tqdm
import scipy.stats
import sys


from str_analysis.utils.misc_utils import parse_interval
from str_analysis.utils.canonical_repeat_unit import compute_canonical_motif
from intervaltree import IntervalTree, Interval

PROJECT_ID = "cmg-analysis"
DATASET_ID = "tandem_repeat_explorer"
TABLE_ID = "catalog"

parser = argparse.ArgumentParser(description="Load data into BigQuery from the annotated catalog JSON file.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-c", "--catalog-path", help="Path to the annotated catalog JSON file",
                    #default="https://github.com/broadinstitute/tandem-repeat-catalog/releases/download/v1.0/repeat_catalog_v1.hg38.1_to_1000bp_motifs.EH.with_annotations.json.gz")
                    default="~/code/tandem-repeat-catalogs/results__2025-11-03/release_draft_2025-11-03/repeat_catalog_v2.hg38.1_to_1000bp_motifs.EH.with_annotations.json.gz")
parser.add_argument("-n", type=int, help="Number of records to read from the catalog")
parser.add_argument("-d", "--known-disease-associated-loci",
                    help="ExpansionHunter catalog with the latest known disease-associated loci",
                    default="https://raw.githubusercontent.com/broadinstitute/str-analysis/refs/heads/main/str_analysis/variant_catalogs/variant_catalog_without_offtargets.GRCh38.json")
parser.add_argument("--tenk10k-tsv", default="../data-prep/tenk10k_str_mt_rows.reformatted.tsv.gz")
parser.add_argument("--hprc100-tsv", default="../data-prep/hprc_lps.2025_05.grouped_by_locus_and_motif.with_biallelic_histogram.tsv.gz")
parser.add_argument("--aou1027-tsv", default="../data-prep/AoULR_phase1_TRGT_Weisburd_v1.0.1_combined.txt.gz")
parser.add_argument("--vamos-ori-motifs-tsv", default="../data-prep/oriMotifFinal.all.tsv.gz")
parser.add_argument("--vamos-eff-motifs-tsv", default="../data-prep/effMotifFinal.all.tsv.gz")
parser.add_argument("--repeat-masker-lookup-json", default="../data-prep/hg38.RepeatMasker.lookup.json.gz")
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
            if 'biallelic_histogram' in col_indices:
                lookup[locus_id]['biallelic_histogram'] = fields[col_indices['biallelic_histogram']]
    
    return lookup


def parse_AoU1027_data_from_tsv(tsv_path):
    """Parse the AoU1027 TSV file line by line and create a lookup dictionary."""
    lookup = {}
    with gzip.open(tsv_path, 'rt') as f:
        header = f.readline().strip().split('\t')
        col_indices = {col: i for i, col in enumerate(header)}
        for line in tqdm.tqdm(f, unit=" lines", unit_scale=True):
            fields = line.strip().split('\t')
            locus_id = fields[col_indices['TRID2']]
            motif_size = len(fields[col_indices['longestPureSegmentMotif']])
            lookup[locus_id] = {
                'min_allele': int(float(fields[col_indices['0thPercentile']]))//motif_size,
                'mode_allele': int(float(fields[col_indices['Mode']]))//motif_size,  # mode allele is in the table as number of repeats
                'stdev': float(fields[col_indices['Stdev']])/motif_size,
                'median': int(float(fields[col_indices['50thPercentile']]))//motif_size,
                '99th_percentile': int(float(fields[col_indices['99thPercentile']]))//motif_size,

                'max_allele': int(float(fields[col_indices['100thPercentile']]))//motif_size,
                'unique_alleles': None, # int(fields[col_indices['numAlleles']]),  this field currently has a bug which needs to be fixed
                'num_called_alleles': int(fields[col_indices['numCalledAlleles']]),

                #"combined_lps_stdev": float(fields[col_indices['combinedLPSStdev']]),
                #"expected_lps_stdev": float(fields[col_indices['expectedCombinedLPSStdev']]),
                "oe_length": float(fields[col_indices['OE_len']]),
                "oe_length_percentile": float(fields[col_indices['OE_len_percentile']]),
            }

    return lookup


def parse_vamos_unique_motif_frequencies_from_tsv(unique_motif_tsv_path, efficient_motif_tsv_path=None):
    """Parse the Vamos motif frequencies TSV file line by line and create a lookup dictionary."""

    if efficient_motif_tsv_path is not None:
        vamos_ori_to_efficient_motif_map = parse_vamos_efficient_motif_lookup(args.vamos_eff_motifs_tsv)
    else:
        vamos_ori_to_efficient_motif_map = {}

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
            motif_list = []
            for motif, count in zip(ori_motifs.split(","), ori_motif_counts.split(",")):
                if motif in motif_list:
                    raise ValueError(f"ERROR: Motif {motif} specified more than once in line: {line}")

                count = int(count)
                total_counts += count
                canonical_motif = compute_canonical_motif(motif)
                if canonical_motif not in canonical_motif_to_motif:
                    canonical_motif_to_motif[canonical_motif] = motif
                    canonical_motif_to_count[canonical_motif] = count
                    motif_list.append(motif)
                else:
                    canonical_motif_to_count[canonical_motif] += count

            ori_to_efficient_motif_map = vamos_ori_to_efficient_motif_map.get(reference_region, {})
            total_frequent_motifs = 0
            for motif in motif_list:
                canonical_motif = compute_canonical_motif(motif)
                count = canonical_motif_to_count[canonical_motif]
                if vamos_ori_to_efficient_motif_map:
                    efficient_motif = ori_to_efficient_motif_map.get(motif, "")
                    value = f"{motif}:{efficient_motif}:{count}"
                else:
                    value = f"{motif}:{count}"

                motif_and_count_list.append(value)


                if int(count) / total_counts >= 0.01:
                    total_frequent_motifs += 1

            result[reference_region] = {
                "unique_motifs": ",".join(motif_list),
                "motif_frequencies": ",".join(motif_and_count_list),
                "total_frequent_motif_count": total_frequent_motifs,
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

tenk10k_lookup = {}
if args.tenk10k_tsv:
    print(f"Parsing Tenk10k data from {args.tenk10k_tsv}")
    tenk10k_lookup = parse_allele_histograms_from_tsv(args.tenk10k_tsv)
    print(f"Parsed {len(tenk10k_lookup):,d} records from the tenk10k table: {args.tenk10k_tsv}")
    
hprc100_lookup = {}
if args.hprc100_tsv:
    print(f"Parsing HPRC100 data from {args.hprc100_tsv}")
    hprc100_lookup = parse_allele_histograms_from_tsv(args.hprc100_tsv)
    print(f"Parsed {len(hprc100_lookup):,d} records from the hprc100 table: {args.hprc100_tsv}")

aou1027_lookup = {}
if args.aou1027_tsv:
    print(f"Parsing AoU1027 data from {args.aou1027_tsv}")
    aou1027_lookup = parse_AoU1027_data_from_tsv(args.aou1027_tsv)
    print(f"Parsed {len(aou1027_lookup):,d} records from the AoU1027 table: {args.aou1027_tsv}")

if hprc100_lookup and aou1027_lookup:
    # compute correlation between medians, stdevs, 
    shared_keys = set(hprc100_lookup.keys()) & set(aou1027_lookup.keys())
    for column in "mode_allele", "stdev", "median", "99th_percentile":
        pearsonr, _ = scipy.stats.pearsonr(
            [hprc100_lookup[k][column] for k in shared_keys], 
            [aou1027_lookup[k][column] for k in shared_keys])
        if pearsonr < 0.66:
            print(f"WARNING: pearson correlation coefficient (R) between the {column} values of HPRC100 and AoU1027 is only: {pearsonr}")
        else:
            print(f"Correlation between the {column} at {len(shared_keys):,d} TR loci in HPRC100 and AoU1027 is {pearsonr:0.2f}")

vamos_ori_motif_frequencies_lookup = parse_vamos_unique_motif_frequencies_from_tsv(
    args.vamos_ori_motifs_tsv,
    args.vamos_eff_motifs_tsv)



print(f"Parsing known disease-associated loci from {args.known_disease_associated_loci}")
if os.path.isfile(args.known_disease_associated_loci):
    fopen = gzip.open if args.known_disease_associated_loci.endswith(".gz") else open 
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

repeat_masker_lookup = {}
if args.repeat_masker_lookup_json:
    print(f"Parsing repeat masker lookup from {args.repeat_masker_lookup_json}")
    with gzip.open(args.repeat_masker_lookup_json, 'rt') as f:
        repeat_masker_lookup = json.load(f)
    print(f"Parsed {len(repeat_masker_lookup):,d} records from the repeat masker lookup: {args.repeat_masker_lookup_json}")
    # convert to repeat masker lookup entries to string format
    print(f"Converting {len(repeat_masker_lookup):,d} repeat masker lookup entries to string format")
    for locus_id, repeat_masker_intervals in tqdm.tqdm(repeat_masker_lookup.items(), unit=" loci", unit_scale=True, total=len(repeat_masker_lookup)):
        repeat_masker_lookup[locus_id] = ", ".join([f"{interval['repFamily']}:{interval['repName']}" for interval in repeat_masker_intervals])  # {interval['repClass']}:
    print(f"Done converting {len(repeat_masker_lookup):,d} repeat masker lookup entries to string format")

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

    # Illumina174k data
    bigquery.SchemaField("AlleleFrequenciesFromIllumina174k", "STRING"),
    bigquery.SchemaField("StdevFromIllumina174k", "FLOAT"),
    bigquery.SchemaField("AlleleFrequenciesFromT2TAssemblies", "STRING"),
    bigquery.SchemaField("StdevFromT2TAssemblies", "FLOAT"),

    # TenK10K data
    bigquery.SchemaField("TenK10K_AlleleHistogram", "STRING"),
    bigquery.SchemaField("TenK10K_BiallelicHistogram", "STRING"),
    bigquery.SchemaField("TenK10K_ModeAllele", "INTEGER"),
    bigquery.SchemaField("TenK10K_Stdev", "FLOAT"),
    bigquery.SchemaField("TenK10K_Median", "INTEGER"),
    bigquery.SchemaField("TenK10K_99thPercentile", "INTEGER"),

    # HPRC100 data
    bigquery.SchemaField("HPRC100_AlleleHistogram", "STRING"),
    bigquery.SchemaField("HPRC100_BiallelicHistogram", "STRING"),
    bigquery.SchemaField("HPRC100_ModeAllele", "INTEGER"),
    bigquery.SchemaField("HPRC100_Stdev", "FLOAT"),
    bigquery.SchemaField("HPRC100_Median", "INTEGER"),
    bigquery.SchemaField("HPRC100_99thPercentile", "INTEGER"),

    # AoU1027 data
    bigquery.SchemaField("AoU1027_MinAllele", "INTEGER"),
    bigquery.SchemaField("AoU1027_ModeAllele", "INTEGER"),
    bigquery.SchemaField("AoU1027_Stdev", "FLOAT"),
    bigquery.SchemaField("AoU1027_Median", "INTEGER"),
    bigquery.SchemaField("AoU1027_99thPercentile", "INTEGER"),
    bigquery.SchemaField("AoU1027_MaxAllele", "INTEGER"),
    bigquery.SchemaField("AoU1027_UniqueAlleles", "INTEGER"),
    bigquery.SchemaField("AoU1027_NumCalledAlleles", "INTEGER"),

    #bigquery.SchemaField("AoU1027_CombinedLPSStdev", "FLOAT"),
    #bigquery.SchemaField("AoU1027_ExpectedLPSStdev", "FLOAT"),
    bigquery.SchemaField("AoU1027_OE_Length", "FLOAT"),
    bigquery.SchemaField("AoU1027_OE_LengthPercentile", "FLOAT"),

    bigquery.SchemaField("RepeatMaskerIntervals", "STRING"),
    bigquery.SchemaField("VamosOriUniqueMotifs", "STRING"),
    bigquery.SchemaField("VamosOriMotifFrequencies", "STRING"),
    bigquery.SchemaField("VamosOriMotifCount", "INTEGER"),
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
elif os.path.isfile(os.path.expanduser(args.catalog_path)):
    catalog = get_json_iterator(os.path.expanduser(args.catalog_path), is_gzipped)
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

    #if record.get("StdevFromT2TAssemblies") and not record.get("AlleleFrequenciesFromT2TAssemblies"):
    #    This happens where loci in source catalogs overlapped
    #    print("WARNING: StdevFromT2TAssemblies is present but AlleleFrequenciesFromT2TAssemblies is missing in record:", pformat(record, indent=4))

    # Add id field
    record["id"] = i + 1

    if record["LocusId"] in tenk10k_lookup:
        counters["rows_with_tenk10k_data"] += 1
        tenk10k_record = tenk10k_lookup[record["LocusId"]]
        record["TenK10K_AlleleHistogram"] = tenk10k_record["allele_size_histogram"]
        record["TenK10K_BiallelicHistogram"] = tenk10k_record.get("biallelic_histogram")
        record["TenK10K_ModeAllele"] = tenk10k_record["mode_allele"]
        record["TenK10K_Stdev"] = tenk10k_record["stdev"]
        record["TenK10K_Median"] = tenk10k_record["median"]
        record["TenK10K_99thPercentile"] = tenk10k_record["99th_percentile"]

        # sanity checks: 
        if record["TenK10K_Median"] > record["TenK10K_99thPercentile"]:
            print(f"WARNING: TenK10K_Median > TenK10K_99thPercentile for {record['LocusId']}: Median == {record['TenK10K_Median']} and 99thPercentile == {record['TenK10K_99thPercentile']}")
        allele_set = set([record["TenK10K_ModeAllele"], record["TenK10K_Median"], record["TenK10K_99thPercentile"]])
        if len(allele_set) > 1 and record["TenK10K_Stdev"] == 0:
            print(f"WARNING: len({allele_set}) > 1 and TenK10K_Stdev == 0 for {record['LocusId']}: len({allele_set}) == {len(allele_set)}")

    if record["LocusId"] in hprc100_lookup:
        counters["rows_with_hprc100_data"] += 1
        hprc100_record = hprc100_lookup[record["LocusId"]]
        record["HPRC100_AlleleHistogram"] = hprc100_record["allele_size_histogram"]
        record["HPRC100_BiallelicHistogram"] = hprc100_record.get("biallelic_histogram")
        record["HPRC100_ModeAllele"] = hprc100_record["mode_allele"]
        record["HPRC100_Stdev"] = hprc100_record["stdev"]
        record["HPRC100_Median"] = hprc100_record["median"]
        record["HPRC100_99thPercentile"] = hprc100_record["99th_percentile"]

        # sanity checks: 
        if record["HPRC100_Median"] > record["HPRC100_99thPercentile"]:
            print(f"WARNING: HPRC100_Median > HPRC100_99thPercentile for {record['LocusId']}: Median == {record['HPRC100_Median']} and 99thPercentile == {record['HPRC100_99thPercentile']}")
        allele_set = set([record["HPRC100_ModeAllele"], record["HPRC100_Median"], record["HPRC100_99thPercentile"]])
        if len(allele_set) > 1 and record["HPRC100_Stdev"] == 0:
            print(f"WARNING: len({allele_set}) > 1 and HPRC100_Stdev == 0 for {record['LocusId']}: len({allele_set}) == {len(allele_set)}")

    if record["LocusId"] in aou1027_lookup:
        counters["rows_with_aou1027_data"] += 1
        aou1027_record = aou1027_lookup[record["LocusId"]]
        record["AoU1027_MinAllele"] = aou1027_record["min_allele"]
        record["AoU1027_ModeAllele"] = aou1027_record["mode_allele"]
        record["AoU1027_Stdev"] = aou1027_record["stdev"]
        record["AoU1027_Median"] = aou1027_record["median"]
        record["AoU1027_99thPercentile"] = aou1027_record["99th_percentile"]
        record["AoU1027_MaxAllele"] = aou1027_record["max_allele"]
        record["AoU1027_UniqueAlleles"] = aou1027_record["unique_alleles"]
        record["AoU1027_NumCalledAlleles"] = aou1027_record["num_called_alleles"]
        #record["AoU1027_CombinedLPSStdev"] = aou1027_record["combined_lps_stdev"]
        #record["AoU1027_ExpectedLPSStdev"] = aou1027_record["expected_lps_stdev"]
        record["AoU1027_OE_Length"] = aou1027_record["oe_length"]
        record["AoU1027_OE_LengthPercentile"] = aou1027_record["oe_length_percentile"]

        # sanity checks: 
        if record["AoU1027_Median"] > record["AoU1027_MaxAllele"]:
            print(f"WARNING: For {record['LocusId']}, AoU1027 Median > 100thPercentile: Median == {record['AoU1027_Median']} and MaxAllele == {record['AoU1027_MaxAllele']}")
        if record["AoU1027_Median"] > record["AoU1027_99thPercentile"]:
            print(f"WARNING: For {record['LocusId']}, AoU1027 Median > 99thPercentile: Median == {record['AoU1027_Median']} and 99thPercentile == {record['AoU1027_99thPercentile']}")
        if record["AoU1027_99thPercentile"] > record["AoU1027_MaxAllele"]:
            print(f"WARNING: For {record['LocusId']}, AoU1027 99thPercentile > 100thPercentile: 99thPercentile == {record['AoU1027_99thPercentile']} and MaxAllele == {record['AoU1027_MaxAllele']}")
        if record["AoU1027_ModeAllele"] > record["AoU1027_MaxAllele"]:
            print(f"WARNING: For {record['LocusId']}, AoU1027 Mode > 100thPercentile: Mode == {record['AoU1027_ModeAllele']} and MaxAllele == {record['AoU1027_MaxAllele']}")
        allele_set = set([record["AoU1027_MinAllele"], record["AoU1027_ModeAllele"], record["AoU1027_Median"], record["AoU1027_99thPercentile"], record["AoU1027_MaxAllele"]])
        allele_set_string = f"0thPercentile == {record['AoU1027_MinAllele']}, Mode == {record['AoU1027_ModeAllele']}, Median == {record['AoU1027_Median']}, 99thPercentile == {record['AoU1027_99thPercentile']}, 100thPercentile == {record['AoU1027_MaxAllele']}"
        #if record["AoU1027_UniqueAlleles"] < len(allele_set):
        #    print(f"WARNING: For {record['LocusId']}, AoU1027 numAlleles == {record['AoU1027_UniqueAlleles']} while {allele_set_string}")
        #if record["AoU1027_Stdev"] > 0 and record["AoU1027_UniqueAlleles"] == 1:
        #    print(f"WARNING: For {record['LocusId']}, AoU1027 Stdev > 0 while numAlleles == 1: Stdev == {record['AoU1027_Stdev']}")
        if len(allele_set) > 1 and record["AoU1027_Stdev"] == 0:
            print(f"WARNING: For {record['LocusId']}, AoU1027 Stdev == 0 while {allele_set_string}")


    if record["ReferenceRegion"] in vamos_ori_motif_frequencies_lookup:
        counters["rows_with_vamos_ori_motif_frequencies"] += 1
        vamos_ori_data = vamos_ori_motif_frequencies_lookup[record["ReferenceRegion"]]
        record["VamosOriUniqueMotifs"] = vamos_ori_data["unique_motifs"]
        record["VamosOriMotifFrequencies"] = vamos_ori_data["motif_frequencies"]
        record["VamosOriMotifCount"] = vamos_ori_data["total_frequent_motif_count"]

    if record["LocusId"] in repeat_masker_lookup:
        counters["rows_with_repeat_masker_intervals"] += 1
        record["RepeatMaskerIntervals"] = repeat_masker_lookup[record["LocusId"]]
        
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

for html_path in "../website/header_template.html", "../index.html", "../locus.html":
    print(f"Update TABLE_ID to '{new_table_id}' in {html_path}")
    with open(html_path, "r") as f:
        html_content = f.read()

    html_content = re.sub(
        r"const TABLE_ID[\s]*=[\s]*'catalog[^']*'",
        f"const TABLE_ID = '{new_table_id}'",
        html_content)

    with open(html_path, "wt") as f:
        f.write(html_content)

print("Done!")

print("\nCounters:")
total_rows = counters["total_rows"]
for key, count in sorted(counters.items(), key=lambda x: x[1], reverse=True):
    if key == "total_rows":
        print(f"{count:10,d}  {key}")
    else:
        print(f"{count:10,d} ( {count/total_rows*100:.1f}%)  {key}")
