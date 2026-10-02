"""Score simulated TR contractions and expansions at every TRExplorer locus near a splice site on Modal L4 GPUs.

Runs the production SpliceAI-lookup code path (SpliceAI, then SAI-10k; distance 10,000, no masking,
GENCODE v50) on the input chunks written by make_full_run_input_chunks.py, with the configuration,
image and code measured by spliceai_l4_cost_and_disk_benchmark.py, which this module imports: an L4
with TensorFloat-32 off, 4 worker processes, one compiled TensorFlow function for the 5 SpliceAI
models, and the REF prediction reused across the alleles of a locus. SAI-10k reads transcript
structures from the GENCODE v50 comprehensive table built by make_transcript_structures_table.py.
Any allele whose SAI-10k prediction is missing or degraded is flagged ("sai10k_problem" in its
record, and listed in the chunk summary).

For each chunk, one L4 container writes to the Modal Volume VOLUME_NAME:
    /outputs/<gene_set>/<run_id>/<chunk>.jsonl.gz      one JSON line per allele whose selected transcript
                                                       has a delta score >= 0.01 (see make_storage_records
                                                       in the benchmark; the "deltas_with_ref_alt" format,
                                                       plus the locus and the allele's target labels)
    /outputs/<gene_set>/<run_id>/<chunk>.summary.json  counts, errors, timing, the input chunk's SHA-256
                                                       and the scoring fingerprint (see SCORING_FILES)
The run ID is a hash of the scoring fingerprint (production code, annotation, transcript table, this
pipeline's scoring code and settings), so scoring with changed code or data writes to a new folder
and a folder never mixes scores from different code; earlier runs' folders are left as they are.
Within a folder each chunk's summary records its input SHA-256, so after the inputs are regenerated
only the chunks whose contents changed are scored again. The outputs of a chunk that is no longer in
the manifest (e.g. its first locus changed) stay in the folder, so each launch also writes
    /outputs/<gene_set>/<run_id>/current_chunks.json    {chunk: input SHA-256} of the manifest launched
as a record. Only the outputs of chunks in the local manifest with the same SHA-256 in their summary are
current; report lists the others. main and report print the run's folder.
An allele whose scoring raised an exception (which the server returns as a result, not an error) is
rescored once in the same container; one that fails again is listed in the summary as failed_twice.
Both are written to a temporary name and renamed, the summary last, so a chunk counts as done only
when its summary exists and names the same input SHA-256 and scoring fingerprint. Rerunning with
unchanged inputs and code skips done chunks and redoes the rest, so a crashed, preempted or
interrupted run can simply be launched again, once the earlier app has stopped (main refuses to
launch while another app of this pipeline is running, and when another run folder already holds
these input chunks; see main).

Usage (from this directory):
    python3 make_transcript_structures_table.py                              # once, locally
    python3 make_full_run_input_chunks.py --gene-set basic                   # once per gene set, locally
    modal run spliceai_full_run_pipeline.py::main --gene-set basic           # dry run: shows what would run
    modal run spliceai_full_run_pipeline.py::main --gene-set basic --max-chunks 2 --launch      # small test
    modal run --detach spliceai_full_run_pipeline.py::main --gene-set basic --launch            # full run
    modal run spliceai_full_run_pipeline.py::report --gene-set basic         # progress, totals, cost, folder
    modal volume get spliceai-tr-full-run /outputs/basic/<run_id> ./full_run_outputs/basic      # download

The ::main or ::report is required because the module has two entry points. --detach keeps the run
going after this terminal disconnects (e.g. the laptop sleeps); all chunks are submitted before
waiting, so Modal holds the queue.
"""
import gzip
import hashlib
import io
import json
import os
import pathlib
import time

import modal

from spliceai_l4_cost_and_disk_benchmark import (
    CORE_PRICE_PER_SECOND, DISTANCE, GIB_PRICE_PER_SECOND, L4_PRICE_PER_SECOND, MIN_DELTA_SCORE_TO_STORE,
    N_WORKER_PROCESSES, PROD_FILES, REFERENCE_DIR, RESERVED_CORES, RESERVED_MEMORY_MIB, SPLICEAI_COMMIT,
    full_run_record, get_cpu_seconds_and_memory_gib, image, score_locus, worker_init)

HERE = pathlib.Path(__file__).parent
VOLUME_NAME = "spliceai-tr-full-run"
VOLUME_MOUNT = "/data"
STORAGE_FORMAT = "deltas_with_ref_alt"
# The local files that decide the scores and the stored records: the production code, the GENCODE
# annotation, the transcript structures SAI-10k reads, the reference's index files (they change with
# the reference), the scoring and storage code in the benchmark module, and this file. Their contents
# and the scoring settings make up a chunk's scoring fingerprint, so changing any of them, even a
# comment in one of the two .py files here, makes every chunk pending again.
SCORING_FILES = [PROD_FILES / name for name in (
    "server.py", "sai10k_predictions.py", "gencode.v50.basic.annotation.txt.gz", "gencode.v50.annotation.txt.gz",
    "gencode.v50.basic.annotation.transcript_annotations.json.gz", "gencode.v50.annotation.transcript_annotations.json.gz")
] + [REFERENCE_DIR / "hg38.fa.gz.fai", REFERENCE_DIR / "hg38.fa.gz.gzi", HERE / "transcripts_hg38.gencode_v50_comprehensive.tsv",
     HERE / "spliceai_l4_cost_and_disk_benchmark.py", HERE / "spliceai_full_run_pipeline.py"]
# Chunks run at the same time; each is one L4. Modal queues the rest.
# The ben-weisburd workspace's plan allows 10 GPUs at a time (Modal email, 2026-10-01); more just queue.
MAX_CONCURRENT_CONTAINERS = 10

volume = modal.Volume.from_name(VOLUME_NAME, create_if_missing=True)
app = modal.App("spliceai-tr-full-run")
pipeline_image = image.add_local_python_source("spliceai_l4_cost_and_disk_benchmark")


def input_path(gene_set, chunk):
    return f"/inputs/{gene_set}/{chunk}.json.gz"


def output_dir(gene_set, run_id):
    return f"/outputs/{gene_set}/{run_id}"


def output_path(gene_set, run_id, chunk):
    return f"{output_dir(gene_set, run_id)}/{chunk}.jsonl.gz"


def summary_path(gene_set, run_id, chunk):
    return f"{output_dir(gene_set, run_id)}/{chunk}.summary.json"


def compute_run_id(scoring_fingerprint):
    """Returns 16 hex digits identifying the run's output folder: the start of the scoring fingerprint.

    The inputs are deliberately not part of it: each chunk's summary records its own input SHA-256, so
    regenerating the inputs reruns only the chunks whose contents changed.
    """
    return scoring_fingerprint[:16]


def compute_scoring_fingerprint(paths, settings):
    """Returns the SHA-256 of the files' contents, in order, followed by the settings dict as JSON."""
    sha256 = hashlib.sha256()
    for path in paths:
        with open(path, "rb") as f:
            for block in iter(lambda: f.read(1 << 20), b""):
                sha256.update(block)
    sha256.update(json.dumps(settings, sort_keys=True).encode())
    return sha256.hexdigest()


def is_chunk_done(summary, input_sha256, scoring_fingerprint):
    """True when a chunk's summary was written for the same input file and the same scoring fingerprint."""
    return summary.get("input_sha256") == input_sha256 and summary.get("scoring_fingerprint") == scoring_fingerprint


def score_indexed_locus(index_and_alleles):
    """Pool task: score one locus's alleles and return them with the locus's position in the chunk."""
    index, alleles = index_and_alleles
    return index, score_locus(alleles)


def make_output_lines(loci, results_by_locus_index):
    """Returns the JSON lines to store, in the chunk's locus and allele order.

    Each record gets its locus ID first, then the allele's target labels (the design's target lengths it
    stands for, e.g. ["99.5th + 3x motif range"]; from the chunk's target_labels), then the scores.
    """
    lines = []
    for index, locus in enumerate(loci):
        for position, result in enumerate(results_by_locus_index[index]):
            record = result["storage_records"][STORAGE_FORMAT]
            if record:
                lines.append(full_run_record(locus["locus"], locus["target_labels"][position], record))
    return lines


def write_atomically(path, data):
    """Writes bytes to path via a temporary file and a rename, so readers never see a partial file."""
    temporary_path = f"{path}.tmp"
    with open(temporary_path, "wb") as f:
        f.write(data)
    os.replace(temporary_path, path)


@app.function(image=pipeline_image, gpu="L4", cpu=RESERVED_CORES, memory=RESERVED_MEMORY_MIB, timeout=4 * 3600,
              retries=modal.Retries(max_retries=2, initial_delay=10.0), volumes={VOLUME_MOUNT: volume},
              max_containers=MAX_CONCURRENT_CONTAINERS)
def score_chunk(gene_set, run_id, chunk, input_sha256, scoring_fingerprint):
    """Score one input chunk and write its records and summary to the run's folder on the volume. Returns the summary."""
    import multiprocessing
    import tensorflow as tf
    t_start = time.time()
    volume.reload()
    existing_summary = pathlib.Path(VOLUME_MOUNT + summary_path(gene_set, run_id, chunk))
    if existing_summary.exists() and is_chunk_done(json.loads(existing_summary.read_text()), input_sha256, scoring_fingerprint):
        return {**json.loads(existing_summary.read_text()), "status": "already done"}

    data = pathlib.Path(VOLUME_MOUNT + input_path(gene_set, chunk)).read_bytes()
    if hashlib.sha256(data).hexdigest() != input_sha256:
        raise ValueError(f"{chunk}: the input file on the volume does not match the manifest's SHA-256")
    loci = json.loads(gzip.decompress(data))["loci"]
    if not tf.config.list_physical_devices("GPU"):
        raise RuntimeError("TensorFlow does not see the GPU; check LD_LIBRARY_PATH")

    # Warm up every worker on one contraction and one expansion (when the chunk has them).
    all_alleles = [variant for locus in loci for variant in locus["alleles"]]
    contractions = [v for v in all_alleles if len(v.split("-")[2]) > 1]
    expansions = [v for v in all_alleles if len(v.split("-")[3]) > 1]
    warmup_alleles = contractions[:1] + expansions[:1]
    context = multiprocessing.get_context("spawn")
    ready_barrier = context.Barrier(N_WORKER_PROCESSES + 1)
    pool = context.Pool(N_WORKER_PROCESSES, initializer=worker_init, initargs=(gene_set, warmup_alleles, ready_barrier))
    ready_barrier.wait(timeout=1800)
    load_seconds = time.time() - t_start

    cpu_seconds_at_start, _ = get_cpu_seconds_and_memory_gib()
    t0 = time.time()
    results_by_locus_index = dict(pool.imap_unordered(
        score_indexed_locus, [(index, locus["alleles"]) for index, locus in enumerate(loci)], chunksize=1))
    # Rescore once each allele whose scoring raised (failed_unexpectedly), e.g. a transient GPU error.
    # An allele that fails again is taken to fail deterministically: it is listed in the summary as
    # failed_twice rather than blocking the chunk, which would otherwise never complete.
    retry_tasks = [((index, position), [result["variant"]]) for index, results in results_by_locus_index.items()
                   for position, result in enumerate(results) if result["failed_unexpectedly"]]
    for (index, position), rescored in pool.imap_unordered(score_indexed_locus, retry_tasks, chunksize=1):
        results_by_locus_index[index][position] = rescored[0]
    wall_seconds = time.time() - t0
    cpu_seconds_at_end, memory_gib = get_cpu_seconds_and_memory_gib()
    pool.close()
    pool.join()

    lines = make_output_lines(loci, results_by_locus_index)
    output_data = gzip.compress(("".join(line + "\n" for line in lines)).encode(), mtime=0)
    os.makedirs(VOLUME_MOUNT + output_dir(gene_set, run_id), exist_ok=True)
    write_atomically(VOLUME_MOUNT + output_path(gene_set, run_id, chunk), output_data)

    results = [r for index in range(len(loci)) for r in results_by_locus_index[index]]
    busy_seconds = sum(r["seconds"] for r in results)
    summary = {
        "gene_set": gene_set,
        "chunk": chunk,
        "input_sha256": input_sha256,
        "scoring_fingerprint": scoring_fingerprint,
        "n_loci": len(loci),
        "n_alleles": len(results),
        "n_alleles_stored": len(lines),
        "errors": [[r["variant"], r["error"]] for r in results if r["error"]],
        "sai10k_problems": [[r["variant"], r["sai10k_problem"]] for r in results if r["sai10k_problem"]],
        "n_alleles_rescored": len(retry_tasks),
        "failed_twice": [r["variant"] for r in results if r["failed_unexpectedly"]],
        "output_bytes": len(output_data),
        "load_and_warmup_seconds": load_seconds,
        "wall_seconds": wall_seconds,
        "busy_seconds": busy_seconds,
        "cores_used_while_busy": (cpu_seconds_at_end - cpu_seconds_at_start) / (busy_seconds / N_WORKER_PROCESSES) / 2,
        "memory_gib": memory_gib,
        "container_seconds": time.time() - t_start,
    }
    write_atomically(VOLUME_MOUNT + summary_path(gene_set, run_id, chunk), json.dumps(summary, indent=1).encode())
    volume.commit()
    return {**summary, "status": "scored"}


def read_manifest(gene_set, require_current_allele_design=True):
    """Returns the input manifest, refusing (by default) inputs built under other allele-design rules.

    The manifest records the fingerprint of the files that decided its alleles
    (make_full_run_input_chunks.compute_allele_design_fingerprint); one that differs from the current
    fingerprint, or is missing, means the chunks may hold obsolete alleles. Only runs locally: the
    allele design modules are not in the Modal image.
    """
    from make_full_run_input_chunks import compute_allele_design_fingerprint
    manifest = json.loads((HERE / "full_run_inputs" / gene_set / "manifest.json").read_text())
    if require_current_allele_design and manifest.get("allele_design_fingerprint") != compute_allele_design_fingerprint(gene_set):
        raise SystemExit(f"full_run_inputs/{gene_set} was built under other allele-design rules or data than the current "
                         f"ones (the manifest's allele_design_fingerprint differs); regenerate the inputs with "
                         f"`python3 make_full_run_input_chunks.py --gene-set {gene_set}`. Unchanged alleles give identical "
                         f"chunks, so nothing already scored is redone.")
    return manifest


def find_other_runs_with_these_inputs(summaries_by_run_id, this_run_id, manifest_chunks):
    """Returns {run_id: number of chunks done} for other runs whose summaries cover some of these exact input chunks."""
    input_sha256s = {c["sha256"] for c in manifest_chunks}
    counts = {run_id: sum(s.get("input_sha256") in input_sha256s for s in summaries.values())
              for run_id, summaries in summaries_by_run_id.items() if run_id != this_run_id}
    return {run_id: n for run_id, n in counts.items() if n}


def read_all_runs_summaries_from_volume(gene_set):
    """Returns {run_id: {chunk: summary}} for every run folder of the gene set on the volume."""
    try:
        entries = volume.listdir(f"/outputs/{gene_set}")
    except (FileNotFoundError, modal.exception.NotFoundError):  # nothing written for this gene set yet
        return {}
    return {os.path.basename(entry.path): read_summaries_from_volume(gene_set, os.path.basename(entry.path))
            for entry in entries if entry.type == modal.volume.FileEntryType.DIRECTORY}


def compute_this_runs_scoring_fingerprint():
    return compute_scoring_fingerprint(SCORING_FILES, {
        "spliceai_commit": SPLICEAI_COMMIT, "distance": DISTANCE, "storage_format": STORAGE_FORMAT,
        "min_delta_score_to_store": MIN_DELTA_SCORE_TO_STORE})


def read_summaries_from_volume(gene_set, run_id):
    """Returns {chunk: summary} for every chunk summary in the run's folder on the volume."""
    summaries = {}
    try:
        entries = volume.listdir(output_dir(gene_set, run_id))
    except (FileNotFoundError, modal.exception.NotFoundError):  # nothing written for this run yet
        return summaries
    for entry in entries:
        if entry.path.endswith(".summary.json"):
            summary = json.loads(b"".join(volume.read_file(entry.path)))
            summaries[summary["chunk"]] = summary
    return summaries


def find_pending_chunks(manifest_chunks, summaries, scoring_fingerprint):
    """Returns the manifest entries with no summary for the same input SHA-256 and scoring fingerprint."""
    return [c for c in manifest_chunks
            if c["chunk"] not in summaries or not is_chunk_done(summaries[c["chunk"]], c["sha256"], scoring_fingerprint)]


def find_other_running_apps_of_this_pipeline(app_list, this_app_id):
    """Returns the entries of `modal app list --json` for other apps of this pipeline that have not stopped."""
    return [a for a in app_list
            if a["description"] == app.name and a["app_id"] != this_app_id and a["state"] != "stopped"]


def dollars_for_chunk(summary):
    """Modal cost of one chunk's container: L4, plus the higher of reserved and used cores and memory."""
    dollars_per_second = (L4_PRICE_PER_SECOND
                          + max(RESERVED_CORES, summary["cores_used_while_busy"]) * CORE_PRICE_PER_SECOND
                          + max(RESERVED_MEMORY_MIB / 1024, summary["memory_gib"]) * GIB_PRICE_PER_SECOND)
    return summary["container_seconds"] * dollars_per_second


@app.local_entrypoint()
def main(gene_set: str, max_chunks: int = 0, launch: bool = False, force: bool = False,
         rescore_already_scored_chunks: bool = False):
    """Upload the pending input chunks and score them; without --launch, only report what would run.

    Refuses to launch:
    - while another app of this pipeline is still running (e.g. an earlier `--detach` run), since its
      queued chunks would all be scored a second time. Stop the old app with `modal app stop <app_id>`,
      or pass --force.
    - when another run folder already holds summaries for these exact input chunks: the scoring
      fingerprint covers both .py files in full, so any edit to them, even one that cannot change the
      scores, starts a new run folder, and relaunching would rescore and pay for every chunk again.
      Revert the edit, or pass --rescore-already-scored-chunks if the scores really should be redone.
    The two are separate flags, so asking to rescore never also lets an earlier run keep running.
    """
    import subprocess
    running = find_other_running_apps_of_this_pipeline(
        json.loads(subprocess.run(["modal", "app", "list", "--json"], check=True, capture_output=True, text=True).stdout),
        app.app_id)
    if running and launch and not force:
        raise SystemExit(f"Another {app.name} app is still running ({', '.join(a['app_id'] for a in running)}); "
                         f"its chunks would be scored twice. Stop it with `modal app stop <app_id>`, or pass --force.")
    manifest = read_manifest(gene_set)
    scoring_fingerprint = compute_this_runs_scoring_fingerprint()
    run_id = compute_run_id(scoring_fingerprint)
    summaries = read_summaries_from_volume(gene_set, run_id)
    all_pending = find_pending_chunks(manifest["chunks"], summaries, scoring_fingerprint)
    pending = all_pending[:max_chunks] if max_chunks else all_pending
    print(f"{gene_set}: outputs in {VOLUME_NAME}:{output_dir(gene_set, run_id)}; {len(manifest['chunks'])} chunks, "
          f"{len(manifest['chunks']) - len(all_pending)} done; this run: {len(pending)} chunks, "
          f"{sum(c['n_loci'] for c in pending):,} loci, {sum(c['n_alleles'] for c in pending):,} alleles")
    other_runs = find_other_runs_with_these_inputs(read_all_runs_summaries_from_volume(gene_set), run_id, manifest["chunks"])
    for other_run_id, n_done in other_runs.items():
        print(f"Run folder {output_dir(gene_set, other_run_id)} already holds {n_done:,} of these input chunks, scored "
              f"with a different scoring fingerprint (code, data or settings).")
    if not launch:
        print("Dry run: pass --launch to upload the inputs and score them.")
        return
    if other_runs and not rescore_already_scored_chunks:
        raise SystemExit("Not launching: these chunks were already scored in another run folder (see above), and this "
                         "launch would score and pay for them again. If the change since then cannot affect the scores "
                         "(e.g. a comment or MAX_CONCURRENT_CONTAINERS), revert it; to rescore anyway, pass "
                         "--rescore-already-scored-chunks.")

    input_dir = HERE / "full_run_inputs" / gene_set
    with volume.batch_upload(force=True) as batch:
        # Which chunk outputs in the run folder belong to this manifest (others are from earlier inputs)
        batch.put_file(io.BytesIO(json.dumps({c["chunk"]: c["sha256"] for c in manifest["chunks"]}, indent=1).encode()),
                       f"{output_dir(gene_set, run_id)}/current_chunks.json")
        for c in pending:
            batch.put_file(str(input_dir / f"{c['chunk']}.json.gz"), input_path(gene_set, c["chunk"]))
    calls = [score_chunk.spawn(gene_set, run_id, c["chunk"], c["sha256"], scoring_fingerprint) for c in pending]
    print(f"Submitted {len(calls)} chunks (at most {MAX_CONCURRENT_CONTAINERS} at a time)", flush=True)
    for c, call in zip(pending, calls):
        try:
            summary = call.get()
            print(f"{c['chunk']}: {summary['status']}, {summary['n_alleles_stored']:,} of {summary['n_alleles']:,} alleles stored, "
                  f"{len(summary['errors'])} errors, {len(summary['sai10k_problems'])} SAI-10k problems, "
                  f"{len(summary['failed_twice'])} failed twice, "
                  f"{summary['container_seconds'] / 60:.0f} min", flush=True)
        except Exception as e:
            print(f"{c['chunk']}: FAILED {type(e).__name__}: {e}", flush=True)


@app.local_entrypoint()
def report(gene_set: str):
    """Print progress, totals and the Modal cost of the chunks done so far.

    The cost counts only each chunk's completed attempt: attempts that were preempted, crashed,
    timed out or failed leave no summary. For the actual bill, run `modal billing report --for today`
    (or over the run's dates) and look at the spliceai-tr-full-run app.
    """
    # Progress is reported even if the allele design files changed after launch (e.g. a comment edit).
    manifest = read_manifest(gene_set, require_current_allele_design=False)
    scoring_fingerprint = compute_this_runs_scoring_fingerprint()
    run_id = compute_run_id(scoring_fingerprint)
    summaries = read_summaries_from_volume(gene_set, run_id)
    pending_chunks = {c["chunk"] for c in find_pending_chunks(manifest["chunks"], summaries, scoring_fingerprint)}
    done = [summaries[c["chunk"]] for c in manifest["chunks"] if c["chunk"] not in pending_chunks]
    n_alleles = sum(s["n_alleles"] for s in done)
    print(f"{gene_set}: outputs in {VOLUME_NAME}:{output_dir(gene_set, run_id)} "
          f"(download: modal volume get {VOLUME_NAME} {output_dir(gene_set, run_id)} ./full_run_outputs/{gene_set})")
    # After an edit to either .py file the fingerprint, and so the folder, changes; point to the folders that
    # hold these inputs' scores so an in-progress or finished run is not mistaken for an empty one.
    for other_run_id, n_done in find_other_runs_with_these_inputs(
            read_all_runs_summaries_from_volume(gene_set), run_id, manifest["chunks"]).items():
        print(f"Run folder {output_dir(gene_set, other_run_id)} holds {n_done:,} of these input chunks, scored with a "
              f"different scoring fingerprint (e.g. before a later edit to spliceai_full_run_pipeline.py or the benchmark "
              f"module); this report covers only {output_dir(gene_set, run_id)}.")
    current_inputs = {(c["chunk"], c["sha256"]) for c in manifest["chunks"]}
    stale = sorted(chunk for chunk, s in summaries.items() if (chunk, s.get("input_sha256")) not in current_inputs)
    if stale:
        print(f"{len(stale):,} chunk outputs in this folder are from inputs no longer in the manifest; they are excluded "
              f"from the totals below and should be ignored after downloading (current_chunks.json in the folder lists "
              f"the current ones): {', '.join(stale)}")
    print(f"{gene_set}: {len(done)} of {len(manifest['chunks'])} chunks done; {n_alleles:,} of {manifest['n_alleles']:,} alleles scored")
    if not done:
        return
    dollars = sum(dollars_for_chunk(s) for s in done)
    print(f"stored {sum(s['n_alleles_stored'] for s in done):,} alleles in {sum(s['output_bytes'] for s in done) / 1e9:.2f} GB; "
          f"{sum(len(s['errors']) for s in done):,} alleles with errors "
          f"({sum(len(s['failed_twice']) for s in done):,} of them raised an exception twice), "
          f"{sum(len(s['sai10k_problems']) for s in done):,} with a missing or degraded SAI-10k prediction")
    print(f"completed attempts: {sum(s['container_seconds'] for s in done) / 3600:,.1f} L4-hours, ${dollars:,.2f}; "
          f"projected total ${dollars / n_alleles * manifest['n_alleles']:,.0f}. Failed or preempted attempts are not "
          f"included; see `modal billing report` for the actual bill.")
