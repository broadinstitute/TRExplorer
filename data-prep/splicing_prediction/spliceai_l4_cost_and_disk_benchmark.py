"""Estimate the Modal cost and disk usage of scoring simulated TR alleles with SpliceAI + SAI-10k on L4 GPUs.

A focused follow-up to an earlier benchmark of Modal GPU types, which found the L4 to be the cheapest
option. It runs only the configuration chosen there: an L4 with TensorFloat-32 off (so scores match the
production CPU code exactly), 4 worker processes per GPU, one compiled TensorFlow function for the
5 SpliceAI models instead of 10 keras predict() calls per transcript window, and the REF prediction
reused across the alleles of a locus.

It scores the alleles in benchmark/tr_benchmark_alleles.<gene_set>.json (see make_tr_expansion_benchmark_alleles.py)
through the production SpliceAI-lookup code path (server.py get_spliceai_scores: SpliceAI, then
SAI-10k), at distance 10,000 with no masking, split across several L4 containers. For every allele
it measures the scoring time and builds the record that would be stored (see make_storage_records),
then scales both to every locus the full run would score. It also checks, in one process, that the
optimized code gives the same full output as the original on a subset of loci.

Usage:
    modal run spliceai_l4_cost_and_disk_benchmark.py --gene-sets basic,comprehensive

Writes benchmark/l4_benchmark_results.<gene_set>.json (per-allele measurements) and
benchmark/l4_benchmark_estimate.json.
"""
import collections
import gzip
import hashlib
import json
import math
import os
import pathlib
import random
import time

import modal

HERE = pathlib.Path(__file__).parent
# server.py, sai10k_predictions.py and the GENCODE v50 annotation files, copied from the SpliceAI-lookup
# repo at commit 78030ec (google_cloud_run_services/ and its docker/spliceai and docker/ref folders).
PROD_FILES = HERE / "spliceai_lookup_files_gencode_v50"
# The benchmark's sample (from make_tr_expansion_benchmark_alleles.py), results and estimate
BENCHMARK_DIR = HERE / "benchmark"
# The transcripts_hg38 table production loads for GENCODE v50, built by make_transcript_structures_table.py
TRANSCRIPT_STRUCTURES_TSV = "/bench/transcripts_hg38.gencode_v50_comprehensive.tsv"
REFERENCE_DIR = pathlib.Path("/Users/weisburd/code/SpliceAI-lookup/google_cloud_run_services/docker/ref/GRCh38")
SPLICEAI_COMMIT = "5854a4ef2662966c6b6d48f4049fad2b7cb150cc"
DISTANCE = 10000
N_WORKER_PROCESSES = 4  # 4 processes already saturate an L4; 8 gave the same throughput
RESERVED_CORES = 4.0
RESERVED_MEMORY_MIB = 32768
# Loci per full-run input chunk (make_full_run_input_chunks.py); the full run starts one container per chunk
FULL_RUN_LOCI_PER_CHUNK = 1000

# Modal list prices: $ per physical core-second (1 core = 2 vCPUs), per GiB-second, per L4-second
CORE_PRICE_PER_SECOND = 0.0000131
GIB_PRICE_PER_SECOND = 0.00000222
L4_PRICE_PER_SECOND = 0.000222

DELTA_SCORES = ("DS_AG", "DS_AL", "DS_DG", "DS_DL")
DELTA_POSITIONS = ("DP_AG", "DP_AL", "DP_DG", "DP_DL")
MIN_DELTA_SCORE_TO_STORE = 0.01
STORAGE_FORMATS = ("deltas_only", "deltas_with_ref_alt")

# tensorflow[and-cuda]==2.16.1 installs the CUDA libraries as pip packages but does not put them on
# the loader path; without this TensorFlow silently runs on the CPU ("Skipping registering GPU devices").
NVIDIA_PIP_LIB_DIRS = ":".join(
    f"/usr/local/lib/python3.10/site-packages/nvidia/{name}/lib"
    for name in ["cublas", "cuda_cupti", "cuda_nvcc", "cuda_nvrtc", "cuda_runtime", "cudnn", "cufft",
                 "curand", "cusolver", "cusparse", "nccl", "nvjitlink"])

# Same steps as google_cloud_run_services/docker/spliceai/Dockerfile (GRCh38) with the CUDA build of
# the same TensorFlow version. Python 3.10 rather than production's 3.9 because the Modal container
# runtime requires >= 3.10; every package is pinned to the version `pip freeze` reports in production.
image = (
    modal.Image.from_registry("python:3.10-slim-bullseye")
    .run_commands(
        "sed -i 's|^deb .*debian-security.*|deb [check-valid-until=no] http://snapshot.debian.org/archive/debian-security/20260901T000000Z bullseye-security main|' /etc/apt/sources.list",
        "apt update && apt-get install --no-install-recommends -y ca-certificates wget bzip2 unzip git "
        "libcurl4-openssl-dev libbz2-dev liblzma-dev zlib1g-dev build-essential libpq-dev",
        "python3 -m pip install 'tensorflow[and-cuda]==2.16.1'",
        "python3 -m pip install flask flask_cors flask-talisman==1.1.0 gunicorn pandas==2.2.2 pyfastx==2.1.0 "
        "psycopg2==2.9.9 psutil keras==3.10.0 numpy==1.26.4 h5py==3.14.0",
        f"python3 -m pip install https://github.com/bw2/SpliceAI/archive/{SPLICEAI_COMMIT}.zip",
    )
    .add_local_file(REFERENCE_DIR / "hg38.fa.gz", "/hg38.fa.gz", copy=True)
    .add_local_file(REFERENCE_DIR / "hg38.fa.gz.fai", "/hg38.fa.gz.fai", copy=True)
    .add_local_file(REFERENCE_DIR / "hg38.fa.gz.fxi", "/hg38.fa.gz.fxi", copy=True)
    .add_local_file(REFERENCE_DIR / "hg38.fa.gz.gzi", "/hg38.fa.gz.gzi", copy=True)
    .add_local_file(PROD_FILES / "gencode.v50.basic.annotation.txt.gz", "/gencode.v50.basic.annotation.txt.gz", copy=True)
    .add_local_file(PROD_FILES / "gencode.v50.annotation.txt.gz", "/gencode.v50.annotation.txt.gz", copy=True)
    .add_local_file(PROD_FILES / "gencode.v50.basic.annotation.transcript_annotations.json.gz",
                    "/gencode.v50.basic.annotation.transcript_annotations.json.gz", copy=True)
    .add_local_file(PROD_FILES / "gencode.v50.annotation.transcript_annotations.json.gz",
                    "/gencode.v50.annotation.transcript_annotations.json.gz", copy=True)
    # The transcripts_hg38 table, which server.py reads from Cloud SQL
    .add_local_file(HERE / "transcripts_hg38.gencode_v50_comprehensive.tsv", TRANSCRIPT_STRUCTURES_TSV, copy=True)
    .env({
        "MODEL_COMMIT": SPLICEAI_COMMIT,
        "TF_CPP_MIN_LOG_LEVEL": "3",
        "OMP_NUM_THREADS": "1",
        "TF_NUM_INTRAOP_THREADS": "1",
        "TF_NUM_INTEROP_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
        "TOOL": "spliceai",
        "GENOME_VERSION": "38",
        "RUNNING_ON_GOOGLE_CLOUD_RUN": "1",
        "LD_LIBRARY_PATH": NVIDIA_PIP_LIB_DIRS,
        # Lets the worker processes share one GPU instead of the first one reserving all of it.
        "TF_FORCE_GPU_ALLOW_GROWTH": "true",
    })
    .add_local_file(PROD_FILES / "server.py", "/root/server.py")
    .add_local_file(PROD_FILES / "sai10k_predictions.py", "/root/sai10k_predictions.py")
)

app = modal.App("spliceai-l4-cost-and-disk-benchmark")

WORKER_STATE = {}


def load_transcript_structures():
    """Parse the transcript-structure table the same way server.get_transcript_structures does."""
    structures = {}
    for line in open(TRANSCRIPT_STRUCTURES_TSV):
        transcript_id, strand, cds_start, cds_end, exon_starts, exon_ends = line.rstrip("\n").split("\t")
        structures[transcript_id] = {
            "EXON_STARTS": [int(s) + 1 for s in exon_starts.rstrip(",").split(",") if s],
            "EXON_ENDS": [int(s) for s in exon_ends.rstrip(",").split(",") if s],
            "CDS_START": int(cds_start) + 1 if cds_start != "\\N" else None,
            "CDS_END": int(cds_end) if cds_end != "\\N" else None,
            "STRAND": strand,
        }
    return structures


def load_production_scoring_code(gene_set):
    """Import production server.py for one gene set, with the Cloud SQL lookup served from memory.

    Each worker process calls this once, for a single gene set: server.py reads GENE_SET when it is
    first imported, so a process must never be reused for another gene set.
    """
    import sys
    from contextlib import contextmanager
    import tensorflow as tf
    # Ampere and newer GPUs default to TensorFloat-32 for convolutions, which shifts about 1 in 6
    # variants' delta scores by 0.001 relative to the production CPU scores.
    tf.config.experimental.enable_tensor_float_32_execution(False)
    os.environ["GENE_SET"] = gene_set
    sys.path.insert(0, "/root")
    import server
    if server.GENE_SET != gene_set:
        raise RuntimeError(f"server.py was already imported for {server.GENE_SET}, not {gene_set}")
    structures = load_transcript_structures()

    @contextmanager
    def in_memory_db_connection():
        yield "in-memory"

    server.get_db_connection = in_memory_db_connection
    server.get_transcript_structures = lambda conn, ids, genome_version, **kwargs: {
        i: structures[i] for i in ids if i in structures}
    return server


def use_compiled_ensemble_and_ref_cache():
    """Replace spliceai.utils.get_delta_scores_for_transcript with a faster version of the same computation.

    The original makes 10 keras predict() calls per transcript window (5 models x REF/ALT), each on a
    batch of one; in Keras 3 every predict() call builds a data pipeline and syncs with the device,
    which dominates the runtime of one sequence on a GPU. This version runs all 5 models in one
    tf.function and averages them with np.mean as the original does. It also reuses the REF
    prediction from WORKER_STATE["ref_prediction_cache"], which score_locus empties at the start of
    each locus: all alleles at a locus share the same REF window, and the original recomputes it
    for every allele. Everything after the model calls is the original code.

    Output is not guaranteed to be bit-identical on a GPU. In same-process comparisons with the
    original (check_optimized_matches_original) on CPU it matched for every allele tested; on an L4
    one allele per gene set differed (90 of 91 basic and 94 of 95 comprehensive alleles identical,
    both differing alleles long insertions), and a similar one-off did not reproduce on rerun. The
    size of those differences was not recorded; the check now keeps both records for differing alleles.
    """
    import numpy as np
    import tensorflow as tf
    import spliceai.utils as spliceai_utils

    compiled_ensembles = {}
    WORKER_STATE["ref_prediction_cache"] = {}

    def predict_with_all_models(models, one_hot_sequence):
        if id(models) not in compiled_ensembles:
            @tf.function(input_signature=[tf.TensorSpec([None, None, 4], tf.float32)], reduce_retracing=True)
            def run_all_models(batch):
                return tf.stack([model(batch, training=False) for model in models])
            compiled_ensembles[id(models)] = run_all_models
        per_model_predictions = compiled_ensembles[id(models)](tf.constant(one_hot_sequence, dtype=tf.float32)).numpy()
        return np.mean(list(per_model_predictions), axis=0)

    def get_delta_scores_for_transcript(x_ref, x_alt, ref, alt, strand, cov, ann):
        ref_cache_key = (x_ref, strand)
        x_ref = spliceai_utils.one_hot_encode(x_ref)[None, :]
        x_alt = spliceai_utils.one_hot_encode(x_alt)[None, :]

        if strand == '-':
            x_ref = x_ref[:, ::-1, ::-1]
            x_alt = x_alt[:, ::-1, ::-1]

        ref_prediction_cache = WORKER_STATE["ref_prediction_cache"]
        if ref_cache_key not in ref_prediction_cache:
            ref_prediction_cache[ref_cache_key] = predict_with_all_models(ann.models, x_ref)
        y_ref = ref_prediction_cache[ref_cache_key].copy()
        y_alt = predict_with_all_models(ann.models, x_alt)

        if strand == '-':
            y_ref = y_ref[:, ::-1]
            y_alt = y_alt[:, ::-1]

        _, trimmed_ref, trimmed_alt = spliceai_utils.trim_shared_bases(ref, alt)
        y_alt_with_inserted_bases = y_alt if len(trimmed_ref) == 1 and len(trimmed_alt) > 1 else None
        y_ref, y_alt = spliceai_utils.align_ref_and_alt_scores(y_ref, y_alt, ref, alt, cov)

        return y_ref, y_alt, y_alt_with_inserted_bases

    spliceai_utils.get_delta_scores_for_transcript = get_delta_scores_for_transcript


def thousandths(score_string):
    """Returns a 3-decimal score string such as "0.026" as an exact whole number of thousandths (26)."""
    return round(float(score_string) * 1000)


def find_sai10k_problem(result):
    """Returns why this allele's SAI-10k prediction is missing or degraded, or None if it is complete.

    The server does not report either case as an error: when SAI-10k raises, it returns the message
    in sai10kPredictionsError, and when the selected transcript is missing from the transcript-
    structure table, SAI-10k runs with zero exons and silently returns no aberrations.
    """
    if result.get("sai10kPredictionsError"):
        return result["sai10kPredictionsError"]
    selected = [t for t in result.get("scores") or [] if t.get("t_id") == result.get("allNonZeroScoresTranscriptId")]
    if selected and "EXON_STARTS" not in selected[0]:
        return f"{selected[0].get('t_id')} is not in the transcript-structure table"
    return None


def make_storage_records(variant, result):
    """Returns the JSON line that would be stored for this allele, in each storage format, or None.

    The stored transcript is the one the server selects for its per-position table and for SAI-10k
    (allNonZeroScoresTranscriptId: MANE Select, then MANE Plus Clinical, then canonical, then the
    largest sum of delta scores), so the SpliceAI scores and the SAI-10k predictions describe the
    same transcript. An allele is stored only if that transcript has a delta score >= 0.01; an
    effect confined to another transcript is not stored. A record has the transcript's ID, four
    delta scores and their positions, the SAI-10k predictions (plus "sai10k_problem" when they are
    missing or degraded; see find_sai10k_problem), and every
    position where the acceptor or donor delta is at least 0.01 in absolute value, as
    [pos, delta acceptor, delta donor] ("deltas_only") or [pos, REF acceptor, ALT acceptor,
    REF donor, ALT donor] ("deltas_with_ref_alt"; the deltas follow from these). Scores at inserted
    bases have no REF counterpart and are not stored.

    Every position with a delta of at least 0.01 is in the transcript's ALL_NON_ZERO_SCORES rows,
    since those include every position where a REF or ALT probability is at least 0.01.

    Returns:
        dict: {storage format: JSON string or None}
    """
    selected = [t for t in result.get("scores") or [] if t.get("t_id") == result.get("allNonZeroScoresTranscriptId")]
    if not selected or max(float(selected[0][k]) for k in DELTA_SCORES) < MIN_DELTA_SCORE_TO_STORE:
        return {record_format: None for record_format in STORAGE_FORMATS}
    selected_transcript = selected[0]

    # The scores are 3-decimal strings. Compare them as whole thousandths: subtracting them as floats
    # turns an exact 0.010 difference (e.g. "0.026" - "0.016") into 0.00999..., dropping it.
    min_delta_in_thousandths = round(MIN_DELTA_SCORE_TO_STORE * 1000)
    rows = [
        row for row in selected_transcript.get("ALL_NON_ZERO_SCORES", [])
        if abs(thousandths(row["AA"]) - thousandths(row["RA"])) >= min_delta_in_thousandths
        or abs(thousandths(row["AD"]) - thousandths(row["RD"])) >= min_delta_in_thousandths]
    shared_fields = {
        "variant": variant,
        "transcript": selected_transcript.get("t_id"),
        **{k: selected_transcript[k] for k in DELTA_SCORES + DELTA_POSITIONS},
        "sai10k": result.get("sai10kPredictions"),
    }
    sai10k_problem = find_sai10k_problem(result)
    if sai10k_problem:
        shared_fields["sai10k_problem"] = sai10k_problem
    positions = {
        "deltas_only": [
            [row["pos"], f"{(thousandths(row['AA']) - thousandths(row['RA'])) / 1000:.3f}",
             f"{(thousandths(row['AD']) - thousandths(row['RD'])) / 1000:.3f}"]
            for row in rows],
        "deltas_with_ref_alt": [[row["pos"], row["RA"], row["AA"], row["RD"], row["AD"]] for row in rows],
    }
    return {
        record_format: json.dumps({**shared_fields, "positions": positions[record_format]}, separators=(",", ":"))
        for record_format in STORAGE_FORMATS
    }


def full_run_record(locus_id, target_labels, record):
    """Returns the JSON line the full run stores for one allele: its locus ID and target labels, then the record."""
    return json.dumps({"locus": locus_id, "target_labels": target_labels, **json.loads(record)}, separators=(",", ":"))


def add_full_run_record_metadata(shards, sample):
    """Rewrites the shards' storage records in the full run's format (see full_run_record), so the disk estimate
    measures what the full run writes. Uses the sample's locus_ids and target_labels_by_locus; a sample made
    before those existed is left as it is."""
    if "locus_ids" not in sample or "target_labels_by_locus" not in sample:
        return
    # Overlapping loci can have identical allele lists, so each list maps to all of its loci, used one at a time
    metadata_by_alleles = collections.defaultdict(list)
    for locus_id, alleles, labels in zip(sample["locus_ids"], sample["alleles_by_locus"], sample["target_labels_by_locus"]):
        metadata_by_alleles[tuple(alleles)].append((locus_id, labels))
    for shard in shards:
        for results in shard["results_by_locus"]:
            locus_id, labels = metadata_by_alleles[tuple(r["variant"] for r in results)].pop()
            for result, allele_labels in zip(results, labels):
                result["storage_records"] = {
                    record_format: full_run_record(locus_id, allele_labels, record) if record else None
                    for record_format, record in result["storage_records"].items()}


def score_allele(variant):
    """Score one allele with the production code path. Returns its timing, output fingerprint and storage records."""
    server = WORKER_STATE["server"]
    t0 = time.perf_counter()
    result = server.get_spliceai_scores(variant, "38", DISTANCE, 0, WORKER_STATE["gene_set"])
    seconds = time.perf_counter() - t0
    # Everything the model output feeds: per-transcript scores with their per-position rows, the
    # per-position rows returned to the page, and the SAI-10k predictions.
    full_output = {k: result.get(k) for k in ("scores", "allNonZeroScores", "sai10kPredictions", "error")}
    return {
        "variant": variant,
        "seconds": seconds,
        "error": result.get("error"),
        "sai10k_problem": find_sai10k_problem(result),
        # Largest of the four delta scores on any transcript, and on the one the server selects (the one stored)
        "max_delta_score_any_transcript": max(
            (float(t[k]) for t in result.get("scores") or [] for k in DELTA_SCORES), default=0.0),
        "max_delta_score_selected_transcript": max(
            (float(t[k]) for t in result.get("scores") or [] if t.get("t_id") == result.get("allNonZeroScoresTranscriptId")
             for k in DELTA_SCORES), default=0.0),
        # The server turns exceptions into results rather than raising: a SpliceAI exception becomes an
        # error of the form "<class '...'>: message", and a SAI-10k exception sets sai10kPredictionsError.
        # Its deliberate errors (no transcript, REF checks) are plain sentences and are not flagged.
        "failed_unexpectedly": (result.get("error") or "").startswith("<class ") or bool(result.get("sai10kPredictionsError")),
        "full_output_sha256": hashlib.sha256(json.dumps(full_output, sort_keys=True, default=str).encode()).hexdigest(),
        "storage_records": make_storage_records(variant, result),
    }


def score_locus(alleles):
    """Score the alleles of one locus in order, sharing REF predictions only within the locus."""
    WORKER_STATE["ref_prediction_cache"].clear()
    return [score_allele(variant) for variant in alleles]


def worker_init(gene_set, warmup_alleles, ready_barrier):
    """Load the scoring code in this worker process, warm it up, then wait until every worker is ready.

    Warming up in the initializer (rather than with warm-up tasks, which a pool hands to whichever
    worker is free) guarantees each worker has built its TF graph and initialized CUDA before
    timing starts. A new ALT length can still pay a one-time setup cost inside the timed region.
    """
    WORKER_STATE["gene_set"] = gene_set
    WORKER_STATE["server"] = load_production_scoring_code(gene_set)
    use_compiled_ensemble_and_ref_cache()
    score_locus(warmup_alleles)
    ready_barrier.wait()


def get_cpu_seconds_and_memory_gib():
    """CPU seconds (user + system) used so far by this process and its workers, and their memory.

    Modal's cpu= setting is a reservation that a container may exceed, and Modal bills the higher of
    reserved and used cores. Memory is the summed proportional set size (PSS), which splits shared
    library pages among the processes that map them; the Modal sandbox has no cgroup memory counter.
    """
    import psutil
    this_process = psutil.Process()
    cpu_seconds = memory_bytes = 0
    for process in [this_process] + this_process.children(recursive=True):
        try:
            times = process.cpu_times()
            cpu_seconds += times.user + times.system
            memory_bytes += process.memory_full_info().pss
        except psutil.NoSuchProcess:
            pass
    return cpu_seconds, memory_bytes / 2**30


@app.function(image=image, gpu="L4", cpu=RESERVED_CORES, memory=RESERVED_MEMORY_MIB, timeout=3 * 3600, retries=0)
def benchmark_shard(gene_set, alleles_by_locus):
    """Score a shard of loci with N_WORKER_PROCESSES workers on one L4 and measure time, CPU and memory.

    Each task is one locus, so all its alleles run in one worker and share the REF prediction.
    """
    import multiprocessing
    import tensorflow as tf
    if not tf.config.list_physical_devices("GPU"):
        raise RuntimeError("TensorFlow does not see the GPU; check LD_LIBRARY_PATH")
    # Warm up on one contraction and one expansion, so both allele shapes have run once.
    all_alleles = [variant for alleles in alleles_by_locus for variant in alleles]
    warmup_alleles = [next(v for v in all_alleles if len(v.split("-")[2]) > 1),
                      next(v for v in all_alleles if len(v.split("-")[3]) > 1)]

    context = multiprocessing.get_context("spawn")
    ready_barrier = context.Barrier(N_WORKER_PROCESSES + 1)
    t_load = time.time()
    pool = context.Pool(N_WORKER_PROCESSES, initializer=worker_init, initargs=(gene_set, warmup_alleles, ready_barrier))
    ready_barrier.wait(timeout=1800)
    load_seconds = time.time() - t_load

    cpu_seconds_at_start, _ = get_cpu_seconds_and_memory_gib()
    t0 = time.time()
    results_by_locus = list(pool.imap_unordered(score_locus, alleles_by_locus, chunksize=1))
    wall_seconds = time.time() - t0
    # Measured while the workers are still alive; their CPU time leaves the process tree when they exit.
    cpu_seconds_at_end, memory_gib = get_cpu_seconds_and_memory_gib()
    pool.close()

    busy_seconds = sum(r["seconds"] for results in results_by_locus for r in results)
    return {
        "gene_set": gene_set,
        "n_loci": len(alleles_by_locus),
        "n_alleles": len(all_alleles),
        "load_and_warmup_seconds": load_seconds,
        "wall_seconds": wall_seconds,
        "busy_seconds": busy_seconds,
        "fraction_of_worker_time_busy": busy_seconds / (N_WORKER_PROCESSES * wall_seconds),
        "cpu_seconds": cpu_seconds_at_end - cpu_seconds_at_start,
        # Physical cores busy on average while the workers were scoring (1 core = 2 vCPUs)
        "cores_used_while_busy": (cpu_seconds_at_end - cpu_seconds_at_start) / (busy_seconds / N_WORKER_PROCESSES) / 2,
        "memory_gib": memory_gib,
        "results_by_locus": results_by_locus,
    }


@app.function(image=image, gpu="L4", cpu=2.0, memory=16384, timeout=3 * 3600, retries=0)
def check_optimized_matches_original(gene_set, alleles_by_locus):
    """Score the loci with the original code, then with the optimized code, in this one process.

    Running both in one process is the only valid equivalence test: across containers even the
    original code's output can differ in the last bit of a float. Each locus runs through
    score_locus, so the REF cache is exercised as in the benchmark.
    """
    WORKER_STATE["gene_set"] = gene_set
    WORKER_STATE["server"] = load_production_scoring_code(gene_set)
    all_alleles = [variant for alleles in alleles_by_locus for variant in alleles]
    original = {variant: score_allele(variant) for variant in all_alleles}
    use_compiled_ensemble_and_ref_cache()
    optimized = {r["variant"]: r for alleles in alleles_by_locus for r in score_locus(alleles)}
    differing = [v for v in all_alleles if original[v]["full_output_sha256"] != optimized[v]["full_output_sha256"]]
    return {
        "gene_set": gene_set,
        "n_alleles": len(all_alleles),
        # Alleles that error out match trivially; these are the ones that actually ran the model.
        "n_alleles_scored": sum(1 for v in all_alleles if not original[v]["error"]),
        "n_identical_full_output": len(all_alleles) - len(differing),
        "differing_alleles": differing,
        # What would be stored for each differing allele, from each code path, to show the size of the difference
        "stored_records_of_differing_alleles": {
            v: {"original": original[v]["storage_records"]["deltas_with_ref_alt"],
                "optimized": optimized[v]["storage_records"]["deltas_with_ref_alt"]} for v in differing},
    }


def bootstrap_range_of_mean(values, n_iterations=5000, seed=1):
    """Returns the 5th and 95th percentile of the mean of values resampled with replacement."""
    rng = random.Random(seed)
    means = sorted(sum(rng.choice(values) for _ in values) / len(values) for _ in range(n_iterations))
    return means[int(0.05 * n_iterations)], means[int(0.95 * n_iterations)]


def estimate_full_run(shards, n_loci_in_full_run):
    """Scale the shards' per-locus time and storage to the full run.

    Cost uses the steady-state time (busy seconds divided by the worker count), since in a long run
    the workers are never idle waiting for the last loci of a shard, plus the measured model loading
    and warm-up (mean load_and_warmup_seconds) for each of the full run's containers, one per chunk of
    FULL_RUN_LOCI_PER_CHUNK loci. Each container is billed at the L4 rate plus the higher of reserved
    and used cores, plus the higher of reserved memory and measured memory with 25% headroom.

    Returns:
        dict: per-allele and per-locus figures and the full-run totals, with 90% ranges from resampling loci
    """
    loci = [results for shard in shards for results in shard["results_by_locus"]]
    alleles = [r for results in loci for r in results]
    busy_seconds = sum(shard["busy_seconds"] for shard in shards)
    dollars_per_container_second = sum(
        shard["busy_seconds"] * (L4_PRICE_PER_SECOND
                                 + max(RESERVED_CORES, shard["cores_used_while_busy"]) * CORE_PRICE_PER_SECOND
                                 + max(RESERVED_MEMORY_MIB / 1024, shard["memory_gib"] * 1.25) * GIB_PRICE_PER_SECOND)
        for shard in shards) / busy_seconds
    container_seconds_per_locus = [sum(r["seconds"] for r in results) / N_WORKER_PROCESSES for results in loci]
    mean_container_seconds_per_locus = sum(container_seconds_per_locus) / len(loci)
    low, high = bootstrap_range_of_mean(container_seconds_per_locus)
    n_full_run_containers = math.ceil(n_loci_in_full_run / FULL_RUN_LOCI_PER_CHUNK)
    load_and_warmup_seconds = n_full_run_containers * sum(shard["load_and_warmup_seconds"] for shard in shards) / len(shards)

    estimate = {
        "n_loci_in_full_run": n_loci_in_full_run,
        "n_sampled_loci": len(loci),
        "n_sampled_alleles": len(alleles),
        "n_full_run_alleles": round(len(alleles) / len(loci) * n_loci_in_full_run),
        "fraction_of_alleles_with_errors": sum(1 for r in alleles if r["error"]) / len(alleles),
        "container_seconds_per_allele": busy_seconds / N_WORKER_PROCESSES / len(alleles),
        "dollars_per_l4_container_hour": dollars_per_container_second * 3600,
        "n_full_run_containers": n_full_run_containers,
        "full_run_load_and_warmup_l4_container_hours": load_and_warmup_seconds / 3600,
        "full_run_l4_container_hours": (mean_container_seconds_per_locus * n_loci_in_full_run + load_and_warmup_seconds) / 3600,
        "full_run_dollars": (mean_container_seconds_per_locus * n_loci_in_full_run + load_and_warmup_seconds) * dollars_per_container_second,
        "full_run_dollars_90pct_range": [(x * n_loci_in_full_run + load_and_warmup_seconds) * dollars_per_container_second
                                         for x in (low, high)],
    }

    stored = [r for r in alleles if r["storage_records"]["deltas_only"]]
    estimate["fraction_of_alleles_stored"] = len(stored) / len(alleles)
    estimate["stored_positions_per_stored_allele"] = sum(
        len(json.loads(r["storage_records"]["deltas_only"])["positions"]) for r in stored) / len(stored)
    for record_format in STORAGE_FORMATS:
        records = [r["storage_records"][record_format] for r in stored]
        raw_bytes = sum(len(record) + 1 for record in records)
        gzip_ratio = len(gzip.compress(("\n".join(records) + "\n").encode(), compresslevel=6)) / raw_bytes
        bytes_per_locus = [sum(len(r["storage_records"][record_format]) + 1 for r in results
                               if r["storage_records"][record_format]) for results in loci]
        mean_bytes_per_locus = sum(bytes_per_locus) / len(loci)
        low, high = bootstrap_range_of_mean(bytes_per_locus)
        estimate[record_format] = {
            "bytes_per_stored_allele": raw_bytes / len(records),
            "gzip_compression_ratio": gzip_ratio,
            "full_run_gb_uncompressed": mean_bytes_per_locus * n_loci_in_full_run / 1e9,
            "full_run_gb_uncompressed_90pct_range": [x * n_loci_in_full_run / 1e9 for x in (low, high)],
            "full_run_gb_gzipped": mean_bytes_per_locus * n_loci_in_full_run * gzip_ratio / 1e9,
        }
    estimate["sai10k_bytes_per_stored_allele"] = sum(
        len(json.dumps(json.loads(r["storage_records"]["deltas_only"])["sai10k"], separators=(",", ":")))
        for r in stored) / len(stored)
    return estimate


@app.local_entrypoint()
def main(gene_sets: str = "basic,comprehensive", n_shards: int = 4, n_equivalence_loci: int = 10, distance_suffix: str = "",
         alleles_from_gene_set: str = "", allow_stale_sample: bool = False):
    """Benchmark each gene set's allele sample on n_shards L4 containers in parallel, then estimate the full run.

    distance_suffix selects an allele sample made with a non-default distance range, e.g. ".501-600bp"
    for tr_benchmark_alleles.basic.501-600bp.json; the output files get the same suffix. With
    n_equivalence_loci 0 the equivalence check is skipped. alleles_from_gene_set scores every gene
    set on the sample made for that gene set (e.g. "comprehensive"), so the same alleles are compared;
    the output files then say which sample was used, and the estimate scales to that gene set's loci.

    Refuses, before any GPU work, a sample whose recorded allele_design_fingerprint differs from the
    current allele design's (make_full_run_input_chunks.compute_allele_design_fingerprint), since its
    cost would not describe the current full run; pass --allow-stale-sample to benchmark an older
    sample on purpose (e.g. the size ladder or the distance bands, which predate the fingerprint).
    """
    from make_full_run_input_chunks import compute_allele_design_fingerprint
    samples, shard_calls, equivalence_calls = {}, {}, {}
    for gene_set in gene_sets.split(","):
        sample_gene_set = alleles_from_gene_set or gene_set
        sample_path = BENCHMARK_DIR / f"tr_benchmark_alleles.{sample_gene_set}{distance_suffix}.json"
        samples[gene_set] = json.load(open(sample_path))
        if not allow_stale_sample and samples[gene_set].get("allele_design_fingerprint") != compute_allele_design_fingerprint(sample_gene_set):
            raise SystemExit(f"{sample_path.name} was made under other allele-design rules or data than the current ones; "
                             f"regenerate it with make_tr_expansion_benchmark_alleles.py, or pass --allow-stale-sample.")
    for gene_set, sample in samples.items():
        alleles_by_locus = sample["alleles_by_locus"]
        shard_calls[gene_set] = [benchmark_shard.spawn(gene_set, alleles_by_locus[i::n_shards]) for i in range(n_shards)]
        if n_equivalence_loci:
            equivalence_calls[gene_set] = check_optimized_matches_original.spawn(gene_set, alleles_by_locus[:n_equivalence_loci])

    estimates = {}
    for gene_set, calls in shard_calls.items():
        shards = [call.get() for call in calls]
        add_full_run_record_metadata(shards, samples[gene_set])
        sample_label = f".alleles_from_{alleles_from_gene_set}" if alleles_from_gene_set else ""
        (BENCHMARK_DIR / f"l4_benchmark_results.{gene_set}{sample_label}{distance_suffix}.json").write_text(json.dumps(shards))
        estimates[gene_set] = {
            "equivalence_check": equivalence_calls[gene_set].get() if gene_set in equivalence_calls else None,
            "shards": [{k: v for k, v in shard.items() if k != "results_by_locus"} for shard in shards],
            # Sample files made before the polymorphism filter have no n_loci_in_full_run; they scored every locus inside a transcript.
            "estimate": estimate_full_run(shards, samples[gene_set].get("n_loci_in_full_run", samples[gene_set]["n_loci_inside_a_transcript"])),
        }
        print(json.dumps({gene_set: {k: estimates[gene_set][k] for k in ("equivalence_check", "estimate")}}, indent=1))
    sample_label = f".alleles_from_{alleles_from_gene_set}" if alleles_from_gene_set else ""
    (BENCHMARK_DIR / f"l4_benchmark_estimate{sample_label}{distance_suffix}.json").write_text(json.dumps(estimates, indent=1))
