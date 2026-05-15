"""Decompose every TRGT allele in the HPRC256 multisample VCF into motif compositions.

Hail Batch orchestrator that:

  1. Streams just the ``.gbi`` grabix index from GCS (~8 KB) and reads the
     total VCF record count via ``grabix size`` — no local copy of the
     9.8 GB VCF is required on the orchestrator host.
  2. Splits records into equal-sized row-range shards (``--shard-size``).
  3. Submits one shard worker per range. Each worker runs
     ``grabix grab vcf $START $END | decompose_hprc_alleles.py --output …``.
  4. Submits a final concat job that merges all shard parquets into one
     ``decomposed_alleles.HPRC256.<timestamp>.parquet`` file.

Mirrors the structure of ``run_trgt_lps.py``. Differs in one way: the worker
script ``decompose_hprc_alleles.py`` is a standalone file (not a HEREDOC), so
the orchestrator uploads it to GCS once and each shard job localizes it via
``s.input()``.
"""

import math
import os
import subprocess
import sys
import tempfile
from datetime import datetime

import hailtop.fs as hfs
from step_pipeline import pipeline, Backend, Localize, Delocalize


DOCKER_IMAGE = "us-central1-docker.pkg.dev/cmg-analysis/docker-repo/str-analysis-with-trgt@sha256:aa267fe4b209dd07dce8ff1cdd3bc98d79f4e9b16a9038012a5708dd6b96152e"

# Same base dir used by run_trgt_lps.py — keeps all HPRC256 artifacts colocated.
GCS_BASE_DIR = "gs://tandem-repeat-catalog/v2.0/HPRC256_TRGT_vcfs_and_lps.2025_12"

# VCF + grabix index (assumed already uploaded; see plan Stage 1 — Critical findings #1).
VCF_GCS_PATH = f"{GCS_BASE_DIR}/trgt-hprc.vcf.gz"
GBI_GCS_PATH = f"{GCS_BASE_DIR}/trgt-hprc.vcf.gz.gbi"

# Shard parquets and concatenated final parquet live under decomposed_alleles/.
SHARDS_GCS_DIR = f"{GCS_BASE_DIR}/decomposed_alleles/shards"
FINAL_GCS_DIR = f"{GCS_BASE_DIR}/decomposed_alleles"

# Worker script staging path — uploaded once by the orchestrator, then localized
# by each shard job. Lives alongside the outputs so it ships with the artifacts.
WORKER_SCRIPT_NAME = "decompose_hprc_alleles.py"
WORKER_SCRIPT_GCS_PATH = f"{FINAL_GCS_DIR}/{WORKER_SCRIPT_NAME}"


# Tiny merge utility for the concat step. Embedded as a HEREDOC because the
# logic is trivial (concat + dedup-check) and doesn't warrant a separate file.
#
# Implementation note: shards are streamed one-at-a-time through the writer
# rather than concat-ed in memory. Holding all 35 shards' parquet data + casting
# allele_data to large_string in memory uses ~10-15 GB and OOM-kills even
# highmem workers. Streaming bounds peak memory to ~one shard worth.
MERGE_SCRIPT = r'''#!/usr/bin/env python3
"""Stream-concatenate shard parquets into a single parquet and run sanity checks."""

import argparse
import sys

import pyarrow as pa
import pyarrow.parquet as pq


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True, help="Output parquet path.")
    parser.add_argument("--expected-rows", type=int, default=None,
                        help="If set, log a diff against the merged row count "
                             "(typically the VCF record count from `grabix size`).")
    parser.add_argument("shards", nargs="+", help="Input shard parquet paths.")
    args = parser.parse_args()

    if not args.shards:
        print("ERROR: no shards passed", file=sys.stderr)
        sys.exit(1)

    # Build the output schema from the first shard, but widen allele_data to
    # large_string (int64 offsets) so the cumulative column can exceed 2GB at
    # full scale without overflowing the int32 offsets of the default `string`.
    sample_schema = pq.read_schema(args.shards[0])
    out_schema = pa.schema([
        pa.field(f.name, pa.large_string() if f.name == "allele_data" else f.type)
        for f in sample_schema
    ])
    ad_idx = sample_schema.get_field_index("allele_data")

    total_rows = 0
    total_haps = 0
    max_haps = 0
    pool = pa.default_memory_pool()
    print(f"Stream-merging {len(args.shards):,d} shard parquets", flush=True)
    with pq.ParquetWriter(args.output, out_schema) as writer:
        for path in sorted(args.shards):
            t = pq.read_table(path)
            n = t.num_rows
            haps_col = t.column("total_allele_number")
            total_haps += pa.compute.sum(haps_col).as_py() or 0
            shard_max = pa.compute.max(haps_col).as_py() or 0
            if shard_max > max_haps:
                max_haps = shard_max
            t = t.set_column(ad_idx, "allele_data",
                             t.column("allele_data").cast(pa.large_string()))
            writer.write_table(t)
            total_rows += n
            # Release the per-shard pyarrow buffers back to the OS before the next
            # iteration; without this, the memory pool accumulates ~one shard's
            # worth of data per iteration even though the Python reference is gone.
            del t, haps_col
            pool.release_unused()
            print(f"  loaded {path}: {n:,d} rows  (running total: {total_rows:,d}, "
                  f"pool_bytes={pool.bytes_allocated():,d})", flush=True)

    print(f"Total rows: {total_rows:,d}")
    print(f"Sum of total_allele_number: {total_haps:,d}")
    print(f"max(total_allele_number) = {max_haps}")

    # Informational row-count comparison. The worker emits one row per sub-TRID,
    # so compound TRGT records (TRID with comma) yield more output rows than
    # input VCF records. The diff is informational.
    if args.expected_rows is not None:
        diff = total_rows - args.expected_rows
        print(f"Final {total_rows:,d} rows; source VCF had {args.expected_rows:,d} records "
              f"(diff {diff:+,d} = sub-TRID rows added by compound records)")

    # Sanity: total_allele_number never exceeds 256 samples * 2 = 512.
    if max_haps > 512:
        print(f"ERROR: total_allele_number > 512", file=sys.stderr)
        sys.exit(1)

    print(f"Wrote {total_rows:,d} rows to {args.output}")


if __name__ == "__main__":
    main()
'''


def get_total_record_count_from_gcs(gbi_gcs_path):
    """Returns the VCF record count reported by ``grabix size``, with no local VCF.

    ``grabix size <path.vcf.gz>`` only ever reads ``<path.vcf.gz>.gbi`` — it
    never opens the ``.vcf.gz`` body. So we stream just the index (~8 KB) into
    a temporary directory, point grabix at the matching ``.vcf.gz`` placeholder
    path next to it (the placeholder file doesn't need to exist), and parse
    the integer count from stdout.

    Args:
        gbi_gcs_path: ``gs://`` URI of the grabix index file
            (``<vcf>.vcf.gz.gbi``).

    Returns:
        Integer record count.
    """
    if not gbi_gcs_path.endswith(".gbi"):
        raise ValueError(f"Expected a .gbi path, got: {gbi_gcs_path}")

    with tempfile.TemporaryDirectory() as tmpdir:
        local_gbi = os.path.join(tmpdir, os.path.basename(gbi_gcs_path))
        local_vcf = local_gbi[: -len(".gbi")]  # path grabix derives `.gbi` from
        print(f"Streaming {gbi_gcs_path} -> {local_gbi}")
        hfs.copy(gbi_gcs_path, local_gbi)

        print(f"Running: grabix size {local_vcf}")
        result = subprocess.run(
            ["grabix", "size", local_vcf],
            check=True,
            capture_output=True,
            text=True,
        )
        return int(result.stdout.strip())


def upload_worker_script(force):
    """Uploads the local ``decompose_hprc_alleles.py`` to GCS if needed.

    Idempotent: skips the upload when the GCS object already exists and
    ``force`` is False. Always overwrites when ``force`` is True so that
    in-flight orchestrator changes to the worker reach the cluster.

    Args:
        force: When True, re-upload even if the GCS object already exists.
    """
    local_worker = os.path.join(os.path.dirname(__file__), WORKER_SCRIPT_NAME)
    if not os.path.exists(local_worker):
        raise FileNotFoundError(f"Worker script not found at {local_worker}")

    if not force and hfs.exists(WORKER_SCRIPT_GCS_PATH):
        print(f"Worker script already at {WORKER_SCRIPT_GCS_PATH} — skipping upload (use --force to overwrite)")
        return

    print(f"Uploading {local_worker} -> {WORKER_SCRIPT_GCS_PATH}")
    hfs.copy(os.path.abspath(local_worker), WORKER_SCRIPT_GCS_PATH)


def shard_ranges(total_records, shard_size):
    """Computes 1-based inclusive grabix row ranges for each shard.

    grabix uses 1-based inclusive line numbering (``grabix grab vcf 1 100``
    returns lines 1..100). The first shard starts at row 1; each subsequent
    shard's start row is the previous shard's end + 1. The final shard's
    end row is clamped to ``total_records``.

    Args:
        total_records: Total number of data records in the VCF (excluding
            header lines — this is what ``grabix size`` returns).
        shard_size: Number of records per shard. Must be >= 1.

    Returns:
        List of ``(start_row, end_row)`` tuples, both 1-based inclusive.
    """
    if shard_size < 1:
        raise ValueError(f"--shard-size must be >= 1, got {shard_size}")
    num_shards = math.ceil(total_records / shard_size)
    ranges = []
    for k in range(num_shards):
        start = k * shard_size + 1
        end = min((k + 1) * shard_size, total_records)
        ranges.append((start, end))
    return ranges


def main():
    bp = pipeline(backend=Backend.HAIL_BATCH_SERVICE, config_file_path="~/.step_pipeline")

    parser = bp.get_config_arg_parser()
    # Calibrated against the patched bw2/trviz image (May 2026 runs):
    # rows 1-60000 averaged ~6ms/record of decompose work, plus ~95s fixed overhead
    # per shard (Docker pull, gcloud auth, 9.2 GB VCF localization). The sample was
    # entirely on chr1, so genome-wide cost may differ. 150k records targets ~15 min
    # mean wall time per shard — well under the 30-min budget, giving room for slow
    # chromosomes and cluster-load variance.
    parser.add_argument("--shard-size", type=int, default=150_000,
                        help="Number of VCF records per shard (default: 150_000, ~15 min/shard).")
    parser.add_argument("-n", "--max-shards", type=int, default=None,
                        help="Only submit the first N shards (for dry-runs).")

    # step_pipeline's parse_known_args swallows -h/--help so it can be re-parsed
    # later by bp.run(). Short-circuit here so `--help` doesn't trigger the
    # side-effecting `grabix size` + worker-script GCS upload below.
    if "-h" in sys.argv[1:] or "--help" in sys.argv[1:]:
        parser.parse_args()
        return

    args = bp.parse_known_args()

    # Stream just the .gbi from GCS and run `grabix size` on it — no need for a
    # local copy of the 9.8 GB VCF on the orchestrator host.
    total_records = get_total_record_count_from_gcs(GBI_GCS_PATH)
    print(f"Total VCF records: {total_records:,d}")

    ranges = shard_ranges(total_records, args.shard_size)
    print(f"Computed {len(ranges):,d} shards of size {args.shard_size:,d}")
    if args.max_shards is not None:
        ranges = ranges[: args.max_shards]
        print(f"--max-shards={args.max_shards} -> submitting only {len(ranges):,d} shards")

    # Upload the standalone worker script to GCS so each shard job can localize it.
    upload_worker_script(force=args.force)

    # Tell step_pipeline which shard outputs already exist in GCS so step 1
    # workers can be skipped when their output parquet is already there.
    # The Step constructors below set `all_outputs_precached=True` which makes
    # them consult this cache; without the precache call the cache is empty
    # and step_pipeline would re-run every shard from scratch.
    bp.precache_file_paths(f"{SHARDS_GCS_DIR}/shard_*.parquet")

    bp.set_name(f"Decompose HPRC alleles: {len(ranges)} shards x {args.shard_size:,d}")

    # Step 1: per-shard worker jobs.
    shard_jobs = []
    shard_output_paths = []
    for k, (start_row, end_row) in enumerate(ranges):
        shard_filename = f"shard_{start_row}_{end_row}.parquet"
        shard_output_paths.append(f"{SHARDS_GCS_DIR}/{shard_filename}")

        s1 = bp.new_step(
            f"decompose shard #{k + 1}: rows {start_row}-{end_row}",
            arg_suffix="decompose",
            step_number=1,
            image=DOCKER_IMAGE,
            cpu=1,
            # highmem (~6.5GB) is needed because deeply-compound records (some with
            # 100+ sub-TRIDs and 100+ unique alleles) produce megabyte-scale rows;
            # n1-standard-1 (~3.75GB) OOM-killed shard 33 on the 4M-record VCF.
            memory="highmem",
            storage="20Gi",
            all_outputs_precached=True,
            localize_by=Localize.GSUTIL_COPY,
            delocalize_by=Delocalize.GSUTIL_COPY,
            output_dir=SHARDS_GCS_DIR,
        )

        s1.command("set -euxo pipefail")
        s1.switch_gcloud_auth_to_user_account()
        s1.command("cd /io/")

        # Localize the VCF, its grabix index, and the standalone worker script.
        # NOTE: input() returns a localized path. grabix expects the .gbi next to
        # the .vcf.gz, so we localize the .gbi into the same dir (Localize.GSUTIL_COPY
        # preserves filenames into a common /io/ working dir).
        local_vcf_path = s1.input(VCF_GCS_PATH)
        s1.input(GBI_GCS_PATH)
        local_worker_path = s1.input(WORKER_SCRIPT_GCS_PATH)

        # grabix uses 1-based inclusive line ranges. start_row/end_row are baked
        # into the command at orchestrator-submit time, so the worker doesn't need
        # to know about shard size or total count.
        s1.command(
            f"grabix grab {local_vcf_path} {start_row} {end_row} "
            f"| python3 {local_worker_path} --output {shard_filename}"
        )
        s1.command(f"ls -lh {shard_filename}")
        s1.output(f"/io/{shard_filename}")
        shard_jobs.append(s1)

    # Step 2: concat all shards into a single timestamped parquet + sanity-check.
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    final_filename = f"decomposed_alleles.HPRC256.{timestamp}.parquet"

    s2 = bp.new_step(
        f"Concat {len(ranges)} shards -> {final_filename}",
        arg_suffix="concat",
        step_number=2,
        image=DOCKER_IMAGE,
        # cpu=2 + highmem (~13 GB) OOM-killed the streaming merge twice in a row.
        # pyarrow's memory pool seems not to aggressively release between
        # `pq.read_table` iterations even with explicit `del`, so the merge
        # accumulates memory across 35 shards. cpu=8 + highmem (~52 GB) gives
        # the merge enough headroom to absorb that growth.
        cpu=8,
        memory="highmem",
        storage="50Gi",
        all_outputs_precached=True,
        localize_by=Localize.GSUTIL_COPY,
        delocalize_by=Delocalize.GSUTIL_COPY,
        output_dir=FINAL_GCS_DIR,
    )
    for j in shard_jobs:
        s2.depends_on(j)

    s2.command("set -euxo pipefail")
    s2.switch_gcloud_auth_to_user_account()
    s2.command("cd /io/")

    local_shard_paths = [s2.input(p) for p in shard_output_paths]

    s2.command(f"cat > merge_decomposed_alleles_shards.py << 'PYEOF'\n{MERGE_SCRIPT}\nPYEOF")
    # `python3 -u` disables stdout/stderr buffering so progress prints reach the
    # Hail Batch log even if the process gets SIGKILL'd before clean shutdown.
    s2.command(
        "python3 -u merge_decomposed_alleles_shards.py "
        f"--output {final_filename} "
        f"--expected-rows {total_records} "
        + " ".join(str(p) for p in local_shard_paths)
    )
    s2.command(f"ls -lh {final_filename}")
    s2.output(f"/io/{final_filename}")

    bp.run()


if __name__ == "__main__":
    main()
