"""Merge the 5 TRExplorer-v2 ExpansionHunter ndjson.gz shards into one valid bracketed JSON array,
and write a genome-wide-representative pilot subset (every Nth locus) for a quick cost dry-run.

Streams line-by-line (no full-file buffering) since the merged catalog is ~5.6M records.
"""
import glob
import gzip
import os

DOWNLOADS_DIR = os.path.expanduser("~/Downloads")
OUT_DIR = os.path.dirname(os.path.abspath(__file__))
PILOT_EVERY_NTH = 50  # ~112K loci pilot subset (2%) for a fast Hail Batch dry-run

shard_paths = sorted(glob.glob(os.path.join(
    DOWNLOADS_DIR, "TR_catalog.shard_*_of_05.5594988_loci.*.ExpansionHunter.ndjson.gz")))
assert len(shard_paths) == 5, f"expected 5 shards, found {len(shard_paths)}: {shard_paths}"

full_out_path = os.path.join(OUT_DIR, "TR_catalog.TRExplorer-v2.5594988_loci.ExpansionHunter.json.gz")
pilot_out_path = os.path.join(OUT_DIR, f"TR_catalog.TRExplorer-v2.pilot_every_{PILOT_EVERY_NTH}th.ExpansionHunter.json.gz")

total_written = 0
pilot_written = 0
with gzip.open(full_out_path, "wt") as full_out, gzip.open(pilot_out_path, "wt") as pilot_out:
    full_out.write("[\n")
    pilot_out.write("[\n")
    for shard_path in shard_paths:
        print(f"Reading {shard_path}...")
        with gzip.open(shard_path, "rt") as f:
            for line in f:
                line = line.rstrip("\n")
                if not line:
                    continue
                if total_written > 0:
                    full_out.write(",\n")
                full_out.write(line)
                if total_written % PILOT_EVERY_NTH == 0:
                    if pilot_written > 0:
                        pilot_out.write(",\n")
                    pilot_out.write(line)
                    pilot_written += 1
                total_written += 1
    full_out.write("\n]\n")
    pilot_out.write("\n]\n")

print(f"Wrote {total_written:,} loci to {full_out_path}")
print(f"Wrote {pilot_written:,} loci to {pilot_out_path}")
