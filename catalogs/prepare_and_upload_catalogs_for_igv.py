import argparse
import datetime
import gzip
import json
import os
import pandas as pd
import re
import subprocess
import sys
import time


parser = argparse.ArgumentParser()
parser.add_argument("--hg38-reference-fasta", default="hg38.fa", help="Path of hg38 reference genome FASTA file")
parser.add_argument("--dry-run", action="store_true", help="Print commands without running them")
parser.add_argument("-k", "--keyword", help="Only process catalogs that contain this keyword")
args = parser.parse_args()


def run(command, verbose=True):
	if verbose:
		print("-"*100)
	command = re.sub("[ \\t]{2,}", "  ", command)  # remove extra spaces
	if verbose:
		print(command)
	# run command and return output
	if not args.dry_run or command.startswith("mkdir"):
		output = subprocess.check_output(command, shell=True, encoding="utf-8")
		return output
	else:
		return ""

def chdir(d):
	print(f"cd {d}")
	os.chdir(d)

scripts_dir = os.path.expanduser("~/code/tandem-repeat-catalogs/scripts")
working_dir = os.path.join(os.path.abspath("."), "data")
run(f"mkdir -p {working_dir}")
chdir(working_dir)


start_time = time.time()

# list of source catalogs, in order. If a locus is defined in more than one catalog (ie. overlapping boundaries,
# same motif), then the definition in the catalog that's earlier in the list will take precedence over definitions in
# subsequent catalogs.
catalogs_in_order = [
	("KnownDiseaseAssociatedLoci", "https://raw.githubusercontent.com/broadinstitute/str-analysis/main/str_analysis/variant_catalogs/variant_catalog_without_offtargets.GRCh38.json"),
	("Illumina174kPolymorphicTRs", "https://storage.googleapis.com/str-truth-set/hg38/ref/other/illumina_variant_catalog.sorted.bed.gz"),
	("UCSC_SimpleRepeatTrack", "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/simpleRepeat.txt.gz"),
	("VamosCatalog_v2.1", "https://storage.googleapis.com/str-truth-set/hg38/ref/other/vamos_catalog.v2.1.bed.gz"),
	("GangSTR_v17", "https://s3.amazonaws.com/gangstr/hg38/genomewide/hg38_ver17.bed.gz"),
	("HipSTR", "https://github.com/HipSTR-Tool/HipSTR-references/raw/master/human/hg38.hipstr_reference.bed.gz"),
	("TRGT", "https://raw.githubusercontent.com/PacificBiosciences/trgt/refs/heads/main/repeats/repeat_catalog.hg38.bed"),
	("Adotto_v1.2", "https://storage.googleapis.com/str-truth-set/hg38/ref/other/adotto_tr_catalog_v1.2.bed.gz"),
	("PopSTR_v2", "https://storage.googleapis.com/str-truth-set/hg38/ref/other/popstr_catalog_v2.bed.gz"),
	("PlatinumTRs", "https://zenodo.org/records/13178746/files/human_GRCh38_no_alt_analysis_set.platinumTRs-v1.0.trgt.bed.gz"),
	("Chiu_et_al", "https://zenodo.org/records/11522276/files/hg38.v1.bed.gz"),
	("PerfectRepeatsInGRCh38", "https://storage.googleapis.com/str-truth-set/hg38/ref/other/colab-repeat-finder/hg38_repeats.motifs_1_to_1000bp.repeats_3x_and_spans_9bp/hg38_repeats.motifs_1_to_1000bp.repeats_3x_and_spans_9bp.bed.gz"),
	("PolymorphicTRsInT2TAssemblies", "https://storage.googleapis.com/str-truth-set-v2/filter_vcf/all_repeats_including_homopolymers_keeping_loci_that_have_overlapping_variants/combined/merged_expansion_hunter_catalog.78_samples.json.gz"),
	("KnownDiseaseAssociatedLoci_July2024", "https://raw.githubusercontent.com/broadinstitute/str-analysis/69dd90ecbc1dcbb23d5ca84ab4022850a283114f/str_analysis/variant_catalogs/variant_catalog_without_offtargets.GRCh38.json"),
	#("MukamelVNTRs", "https://storage.googleapis.com/str-truth-set/hg38/ref/other/mukamel_VNTR_catalog.bed.gz"),
]


catalog_paths = {}
for catalog_name, url in catalogs_in_order:
	if url.endswith(".bed.gz") or url.endswith(".bed"):
		suffix = ".bed.gz"
	elif url.endswith(".txt.gz"):
		suffix = ".txt.gz"
	elif url.endswith(".json"):
		suffix = ".json"
	elif url.endswith(".json.gz"):
		suffix = ".json.gz"
	else:
		raise ValueError(f"Unexpected file extension in: {url}")
	output_path = f"{catalog_name}{suffix}"
	if not os.path.isfile(output_path):
		bgzip_command = "| bgzip" if url.endswith(".bed") else ""
		run(f"curl --silent {url} {bgzip_command} > {output_path}.tmp && mv {output_path}.tmp {output_path}")
	catalog_paths[catalog_name] = output_path

# preprocess catalog of known disease-associated loci: split compound definitions
for key in "KnownDiseaseAssociatedLoci_July2024", "KnownDiseaseAssociatedLoci":
	output_path = catalog_paths[key].replace(".json", ".split.json")
	if not os.path.isfile(output_path):
		run(f"python3 -m str_analysis.split_adjacent_loci_in_expansion_hunter_catalog {catalog_paths[key]}")
	catalog_paths[key] = output_path

	# change motif definition for the RFC1 locus from AARRG => AAAAG since our catalog doesn't currently support IUPAC codes
	run(f"sed -i 's/AARRG/AAAAG/g' {catalog_paths[key]}")

	# compute stats for primary disease-associated loci
	with open(catalog_paths[key]) as f:
		known_disease_associated_loci = json.load(f)

	primary_disease_associated_loci = [
		x for x in known_disease_associated_loci if x["Diseases"] and (
			x["LocusId"].startswith("HOXA") or x["LocusId"].startswith("ARX") or "_" not in x["LocusId"]
		)
	]

	primary_disease_associated_loci_path = catalog_paths[key].replace(
		".json", ".primary_disease_associated_loci.json")

	with open(primary_disease_associated_loci_path, "wt") as f:
		json.dump(primary_disease_associated_loci, f, indent=4)

	run(f"python3 -u -m str_analysis.convert_expansion_hunter_catalog_to_bed {primary_disease_associated_loci_path} -o {key}.bed.gz")
	catalog_paths[key] = f"{key}.bed.gz"

# convert TRGT catalog to bgzipped bed
if "TRGT" in catalog_paths and (
	not args.keyword or args.keyword.lower() in "TRGT".lower() or args.keyword.lower() in catalog_paths["TRGT"].lower()):
	path_after_conversion = catalog_paths["TRGT"].replace(".bed.gz", ".catalog.json.gz")
	if not os.path.isfile(path_after_conversion):
		run(f"python3 -u -m str_analysis.convert_trgt_catalog_to_expansion_hunter_catalog -r {args.hg38_reference_fasta} {catalog_paths['TRGT']} -o {path_after_conversion}")
		run(f"python3 -u -m str_analysis.convert_expansion_hunter_catalog_to_bed {path_after_conversion} -o {catalog_paths['TRGT']}")


# convert PlatinumTR catalog to bgzipped bed
if "PlatinumTRs" in catalog_paths and (
	not args.keyword or args.keyword.lower() in "PlatinumTRs".lower() or args.keyword.lower() in catalog_paths["PlatinumTRs"].lower()):

	#platinum_trs_original_trgt_bed =  catalog_paths["PlatinumTRs"].replace(".bed", ".trgt.bed")
	#run(f"mv {catalog_paths['PlatinumTRs']} {platinum_trs_original_trgt_bed}")

	path_after_conversion = catalog_paths["PlatinumTRs"].replace(".bed.gz", ".catalog.json.gz")
	if not os.path.isfile(path_after_conversion):
		run(f"python3 -u -m str_analysis.convert_trgt_catalog_to_expansion_hunter_catalog -r {args.hg38_reference_fasta} {catalog_paths['PlatinumTRs']} -o {path_after_conversion}")
		run(f"python3 -u -m str_analysis.convert_expansion_hunter_catalog_to_bed {path_after_conversion} -o {catalog_paths['PlatinumTRs']}")

else:
	print("WARNING: PlatinumTRs catalog not included in list of catalogs")

# convert GangSTR catalog to bgzipped bed
if "GangSTR_v17" in catalog_paths and (
	not args.keyword or args.keyword.lower() in "GangSTR_v17".lower() or args.keyword.lower() in catalog_paths["GangSTR_v17"].lower()):
	path_after_conversion = catalog_paths["GangSTR_v17"].replace(".bed.gz", ".catalog.json.gz")
	if not os.path.isfile(path_after_conversion):
		run(f"python3 -u -m str_analysis.convert_gangstr_spec_to_expansion_hunter_catalog --verbose {catalog_paths['GangSTR_v17']} -o {path_after_conversion}")
		run(f"python3 -u -m str_analysis.convert_expansion_hunter_catalog_to_bed {path_after_conversion} -o GangSTR_v17.bed.gz")

	catalog_paths["GangSTR_v17"] = "GangSTR_v17.bed.gz"
else:
	print("WARNING: GangSTR catalog not included in list of catalogs")

# convert HipSTR catalog to bgzipped bed
if "HipSTR" in catalog_paths and (
	not args.keyword or args.keyword.lower() in "HipSTR".lower() or args.keyword.lower() in catalog_paths["HipSTR"].lower()):
	original_catalog_path = catalog_paths["HipSTR"].replace(".bed.gz", ".original.bed.gz")
	if not os.path.isfile(original_catalog_path):
		run(f"mv {catalog_paths['HipSTR']} {original_catalog_path}")
		run(f"python3 -u {scripts_dir}/convert_hipstr_catalog_to_regular_bed_file.py {original_catalog_path} -o {catalog_paths['HipSTR']}")
else:
	print("WARNING: HipSTR catalog not included in list of catalogs")

# convert UCSC track to BED format
if "UCSC_SimpleRepeatTrack" in catalog_paths and (
	not args.keyword or args.keyword.lower() in "UCSC_SimpleRepeatTrack".lower() or args.keyword.lower() in catalog_paths["UCSC_SimpleRepeatTrack"].lower()):
	path_after_conversion = catalog_paths["UCSC_SimpleRepeatTrack"].replace(".txt.gz", "") + ".bed.gz"
	if not os.path.isfile(path_after_conversion):
		run(f"python3 -u {scripts_dir}/convert_ucsc_simple_repeat_track_to_bed.py {catalog_paths['UCSC_SimpleRepeatTrack']} -o {path_after_conversion}")
	catalog_paths["UCSC_SimpleRepeatTrack"] = path_after_conversion

if "PolymorphicTRsInT2TAssemblies" in catalog_paths and (
	not args.keyword or args.keyword.lower() in "PolymorphicTRsInT2TAssemblies".lower() or args.keyword.lower() in catalog_paths["PolymorphicTRsInT2TAssemblies"].lower()):

	path_after_conversion = "PolymorphicTRsInT2TAssemblies.bed.gz"
	if not os.path.isfile(path_after_conversion):
		run(f"python3 -u -m str_analysis.convert_expansion_hunter_catalog_to_bed {catalog_paths['PolymorphicTRsInT2TAssemblies']} -o {path_after_conversion}")

	catalog_paths["PolymorphicTRsInT2TAssemblies"] = path_after_conversion



print("Catalogs: ")
for label, path in catalog_paths.items():
	print(f"{label:>30s}: {path}")

# process Chiu et al catalog
if "Chiu_et_al" in catalog_paths and (
	not args.keyword or args.keyword.lower() in "Chiu_et_al".lower() or args.keyword.lower() in catalog_paths["Chiu_et_al"].lower()):
	original_catalog_path = catalog_paths["Chiu_et_al"].replace(".bed.gz", ".original.bed.gz")
	if not os.path.isfile(original_catalog_path):
		run(f"mv {catalog_paths['Chiu_et_al']} {original_catalog_path}")
		# convert the comprehensive catalog from Chiu et al. to a regular bed file (removing the header and converting to 0-basd coords)
		run(f"gunzip -c {original_catalog_path} | tail -n +3 | cut -f 1-4 | awk 'BEGIN {{OFS=\"\\t\"}} {{ print( $1, $2 - 1, $3, $4 ) }}' | bgzip > {catalog_paths['Chiu_et_al']}")

all_stats_tsv_paths = {}
for catalog_name, _ in catalogs_in_order:
	path = catalog_paths[catalog_name]

	if args.keyword and args.keyword.lower() not in catalog_name.lower() and args.keyword.lower() not in path.lower():
		continue

	if not os.path.isfile(f"{path}.gbi"):
		run(f"grabix index {path}")

	# just compute stats
	stats_tsv_path = re.sub("(.json|.bed)(.gz)?$", "", path) + ".catalog_stats.tsv"
	if not os.path.isfile(stats_tsv_path):
		print(f"Generating {stats_tsv_path}")
		run(f"python3 -m str_analysis.compute_catalog_stats --verbose {path}")

	all_stats_tsv_paths[catalog_name] = stats_tsv_path

# upload catalogs to Google Cloud Storage
google_storage_bucket = "gs://tandem-repeat-catalog"
run(f"gsutil -m cp -n {' '.join([path+'*' for path in catalog_paths.values() ])}  {google_storage_bucket}/other_catalogs/")
run(f"gsutil -m cp -n {' '.join([path for path in all_stats_tsv_paths.values()])} {google_storage_bucket}/other_catalogs/")

for catalog_name, path in catalog_paths.items():
	num_records = run(f"grabix size {path}", verbose=False).strip()
	print(f"{num_records:10s}    {catalog_name:<30s}   {os.path.join(google_storage_bucket, 'other_catalogs', os.path.basename(path))}")

	#with gzip.open(path, "rt") as f:
	#	for i, line in enumerate(f):
	#		print(f"\t#{i+1}:  {line.strip()}")
	#		if i >= 2:
	#			break


# process TR_Explorer catalog files
tr_explorer_catalog_paths = {}
tr_explorer_catalog_paths["TR_Explorer_catalog"] = "TR_Explorer.repeat_catalog_v1.hg38.1_to_1000bp_motifs.bed.gz"

variant_clusters_path = "variation_clusters_v1.hg38.TRGT.bed.gz"
variant_clusters_path_after_conversion = "variation_clusters_v1.hg38.bed.gz"
if not os.path.isfile(variant_clusters_path_after_conversion):
	run(f"python3 -u -m str_analysis.convert_trgt_catalog_to_expansion_hunter_catalog -r {args.hg38_reference_fasta}  -o temp.json.gz {variant_clusters_path}")
	run(f"python3 -u -m str_analysis.convert_expansion_hunter_catalog_to_bed  temp.json.gz -o {variant_clusters_path_after_conversion}")
	run("rm temp.json.gz")

tr_explorer_catalog_paths["variation_clusters"] = variant_clusters_path_after_conversion

for catalog_name, path in tr_explorer_catalog_paths.items():
	if not os.path.isfile(f"{path}.gbi"):
		run(f"grabix index {path}")

	num_records = run(f"grabix size {path}", verbose=False).strip()
	print(f"{num_records:10s}    {catalog_name:<30s}   {os.path.join(google_storage_bucket, 'other_catalogs', os.path.basename(path))}")

run(f"gsutil -m cp -n {' '.join([path+'*' for path in tr_explorer_catalog_paths.values() ])}  {google_storage_bucket}/v1.0/")

diff = time.time() - start_time
print(f"Done. Took {diff//3600:.0f}h, {(diff%3600)//60:.0f}m, {diff%60:.0f}s")

