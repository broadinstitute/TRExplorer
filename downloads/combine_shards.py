import collections
import glob
import gzip
import json
import re
import tqdm

# Find all sharded JSON files matching the pattern
pattern = "TR_catalog.shard_*_of_*.4863041_loci.20251021*.json.gz"
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
output_path = f"TR_catalog.{suffix}.json.gz"

# Combine all shards into a single output file
total_rows = 0
key_counters = collections.Counter()
with gzip.open(output_path, "wt") as out_f:
    out_f.write("[")
    first_object = True
    for path in sorted(input_file_paths):
        print(f"Processing {path}")
        with gzip.open(path, "rt") as f:
            for line in tqdm.tqdm(f, unit=" lines", unit_scale=True):
                if first_object:
                    first_object = False
                else:
                    out_f.write(", ")

                record = json.loads(line)
                del record["chrom"]
                del record["start_0based"]
                del record["end_1based"]
                record["LocusStructure"] = f'({record["ReferenceMotif"]})*'

                json.dump(record, out_f, indent=4)

                total_rows += 1
                for k, v in record.items():
                    if v is not None:
                        key_counters[k] += 1

    out_f.write("\n]")

print(f"Wrote {total_rows:,d} records from {len(input_file_paths)} shards to {output_path}")

for k, c in sorted(key_counters.items(), key=lambda i: -i[1]):
    print(f"{c:10,d} out of {total_rows:,d} ({c/total_rows:.1%}) rows have {k}")

