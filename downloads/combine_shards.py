import collections
import glob
import gzip
import json
import re
import tqdm

# Find all sharded JSON files matching the pattern
pattern = "TR_catalog.shard_*_of_*.4863041_loci.20251030*.json.gz"
input_file_paths = glob.glob(pattern)

if not input_file_paths:
    print(f"{pattern} files not found in the current directory")
    exit()

print(f"Input paths:\n   ", "\n    ".join(sorted(input_file_paths)))


# Construct output filename
first_file = input_file_paths[0]
# Match pattern: TR_catalog.shard_XX_of_YY.ZZZZZ_loci.YYYYMMDD_HHMMSS.json.gz
match = re.match(r"TR_catalog\.shard_\d+_of_\d+\.(.+)\.json\.gz", first_file)
if not match:
    print(f"Could not parse filename pattern from {first_file}")

suffix = match.group(1)  # e.g., "4863041_loci.20250611_144717"
output_json_path = f"TR_catalog.{suffix}.json.gz"
output_tsv_path = f"TR_catalog.{suffix}.tsv.gz"

tsv_header = [
    "LocusId",
    "ReferenceRegion",
    "ReferenceMotif",
    "CanonicalMotif",
    "MotifSize",
    "NumRepeatsInReference",
    "ReferenceRepeatPurity",
    #"NsInFlanks",
    "TRsInRegion",
    "Source",
    "FlanksAndLocusMappability",
    "GencodeGeneRegion",
    "RefseqGeneRegion",
    "ManeGeneRegion",
    #"LocusStructure",
    #"HPRC256_AlleleHistogram",
    "HPRC256_ModeAllele",
    "HPRC256_Stdev",
    "HPRC256_Median",
    "HPRC256_99thPercentile",
    "AoU1027_ModeAllele",
    "AoU1027_Stdev",
    "AoU1027_Median",
    "AoU1027_99thPercentile",
    "AoU1027_MaxAllele",
    "AoU1027_OE_Length",
    "AoU1027_OE_LengthPercentile",
    "GencodeGeneName",
    #"GencodeGeneId",
    #"GencodeTranscriptId",
    "RefseqGeneName",
    #"RefseqGeneId",
    #"RefseqTranscriptId",
    #"TenK10K_AlleleHistogram",
    "TenK10K_ModeAllele",
    "TenK10K_Stdev",
    "TenK10K_Median",
    "TenK10K_99thPercentile",
    "ManeGeneName",
    #"ManeGeneId",
    #"ManeTranscriptId",
    #"StdevFromT2TAssemblies",
    #"AlleleFrequenciesFromT2TAssemblies",
    "VariationCluster",
    "VariationClusterSizeDiff",
    #"KnownDiseaseAssociatedMotif",
    #"AlleleFrequenciesFromIllumina174k",
    "StdevFromIllumina174k",
    #"KnownDiseaseAssociatedLocus",
]

# Combine all shards into a single output file
total_rows = 0
key_counters = collections.Counter()
with gzip.open(output_json_path, "wt") as out_json, gzip.open(output_tsv_path, "wt") as out_tsv:
    out_json.write("[")
    out_tsv.write("\t".join(tsv_header) + "\n")
    first_object = True
    for path in sorted(input_file_paths):
        print(f"Processing {path}")
        with gzip.open(path, "rt") as f:
            for line in tqdm.tqdm(f, unit=" lines", unit_scale=True):
                if first_object:
                    first_object = False
                else:
                    out_json.write(", ")
                
                # process json
                record = json.loads(line)
                del record["chrom"]
                del record["start_0based"]
                del record["end_1based"]
                record["LocusStructure"] = f'({record["ReferenceMotif"]})*'

                json.dump(record, out_json, indent=4)

                # update counters
                total_rows += 1
                for k, v in record.items():
                    if v is not None:
                        key_counters[k] += 1

                # write to tsv
                out_tsv.write("\t".join([str(record.get(k, "")) for k in tsv_header]) + "\n")

    out_json.write("\n]")


print(f"Wrote {total_rows:,d} records from {len(input_file_paths)} shards to {output_json_path} and {output_tsv_path}")

for k, c in sorted(key_counters.items(), key=lambda i: -i[1]):
    print(f"{c:10,d} out of {total_rows:,d} ({c/total_rows:.1%}) rows have {k}")
print("-" * 80)


"""
 4,863,041 out of 4,863,041 (100.0%) rows have LocusId
 4,863,041 out of 4,863,041 (100.0%) rows have ReferenceRegion
 4,863,041 out of 4,863,041 (100.0%) rows have ReferenceMotif
 4,863,041 out of 4,863,041 (100.0%) rows have CanonicalMotif
 4,863,041 out of 4,863,041 (100.0%) rows have MotifSize
 4,863,041 out of 4,863,041 (100.0%) rows have NumRepeatsInReference
 4,863,041 out of 4,863,041 (100.0%) rows have ReferenceRepeatPurity
 4,863,041 out of 4,863,041 (100.0%) rows have NsInFlanks
 4,863,041 out of 4,863,041 (100.0%) rows have TRsInRegion
 4,863,041 out of 4,863,041 (100.0%) rows have Source
 4,863,041 out of 4,863,041 (100.0%) rows have FlanksAndLocusMappability
 4,863,041 out of 4,863,041 (100.0%) rows have GencodeGeneRegion
 4,863,041 out of 4,863,041 (100.0%) rows have RefseqGeneRegion
 4,863,041 out of 4,863,041 (100.0%) rows have ManeGeneRegion
 4,863,041 out of 4,863,041 (100.0%) rows have LocusStructure
 4,609,399 out of 4,863,041 (94.8%) rows have HPRC256_AlleleHistogram
 4,609,399 out of 4,863,041 (94.8%) rows have HPRC256_ModeAllele
 4,609,399 out of 4,863,041 (94.8%) rows have HPRC256_Stdev
 4,609,399 out of 4,863,041 (94.8%) rows have HPRC256_Median
 4,609,399 out of 4,863,041 (94.8%) rows have HPRC256_99thPercentile
 4,441,086 out of 4,863,041 (91.3%) rows have AoU1027_ModeAllele
 4,441,086 out of 4,863,041 (91.3%) rows have AoU1027_Stdev
 4,441,086 out of 4,863,041 (91.3%) rows have AoU1027_Median
 4,441,086 out of 4,863,041 (91.3%) rows have AoU1027_99thPercentile
 4,441,086 out of 4,863,041 (91.3%) rows have AoU1027_MaxAllele
 4,441,086 out of 4,863,041 (91.3%) rows have AoU1027_OE_Length
 4,441,086 out of 4,863,041 (91.3%) rows have AoU1027_OE_LengthPercentile
 2,912,585 out of 4,863,041 (59.9%) rows have GencodeGeneName
 2,912,585 out of 4,863,041 (59.9%) rows have GencodeGeneId
 2,912,585 out of 4,863,041 (59.9%) rows have GencodeTranscriptId
 2,848,288 out of 4,863,041 (58.6%) rows have RefseqGeneName
 2,848,288 out of 4,863,041 (58.6%) rows have RefseqGeneId
 2,848,288 out of 4,863,041 (58.6%) rows have RefseqTranscriptId
 2,459,897 out of 4,863,041 (50.6%) rows have TenK10K_AlleleHistogram
 2,459,897 out of 4,863,041 (50.6%) rows have TenK10K_ModeAllele
 2,459,897 out of 4,863,041 (50.6%) rows have TenK10K_Stdev
 2,459,897 out of 4,863,041 (50.6%) rows have TenK10K_Median
 2,459,897 out of 4,863,041 (50.6%) rows have TenK10K_99thPercentile
 1,956,125 out of 4,863,041 (40.2%) rows have ManeGeneName
 1,956,125 out of 4,863,041 (40.2%) rows have ManeGeneId
 1,956,125 out of 4,863,041 (40.2%) rows have ManeTranscriptId
 1,943,676 out of 4,863,041 (40.0%) rows have StdevFromT2TAssemblies
 1,928,301 out of 4,863,041 (39.7%) rows have AlleleFrequenciesFromT2TAssemblies
   593,325 out of 4,863,041 (12.2%) rows have VariationCluster
   593,325 out of 4,863,041 (12.2%) rows have VariationClusterSizeDiff
   419,191 out of 4,863,041 (8.6%) rows have KnownDiseaseAssociatedMotif
   173,885 out of 4,863,041 (3.6%) rows have AlleleFrequenciesFromIllumina174k
   173,885 out of 4,863,041 (3.6%) rows have StdevFromIllumina174k
        63 out of 4,863,041 (0.0%) rows have KnownDiseaseAssociatedLocus
"""