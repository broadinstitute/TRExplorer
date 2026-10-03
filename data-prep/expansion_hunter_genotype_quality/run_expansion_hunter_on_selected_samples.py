"""Genotype the TRExplorer-v2 catalog, with the gap-purity rule's extended locus definitions added
alongside the original ones, using ExpansionHunter-bw2 (EHv5-bw2-optimized, optimized-streaming) on
the samples in selected_1kGP_samples.tsv (see select_1kGP_samples.py).

The catalog holds 5,657,854 locus definitions: all 5,599,658 v2 loci plus the 58,196 extended
definitions that survived deduplication, built by run_extend_v2_catalog_boundaries.sh and then
run_convert_extended_catalog_to_eh_catalog.sh. Both definitions of an extended locus are genotyped in
the same run, against the same reads, so comparing them measures the boundary and nothing else. Each
LocusId is the definition's own "{chrom}-{start_0based}-{end_1based}-{motif}", so an extended
definition and the locus it grew from are linked by overlap rather than by a shared id.

Config is cpu=2/threads=4/highmem, one job per sample. Measured by the 2026-08-01 cost benchmark on 3
replicate WGS samples against the 5,594,988-locus catalog:

    cpu=1/threads=2   OOM before genotyping starts (highmem at 1 core is only 6.5GB)
    cpu=2/threads=4   5.77, 5.79, 5.87 hours   $0.482, $0.484, $0.490
    cpu=4/threads=8   1.96, 4.04, 4.44 hours   $0.365, $0.637, $0.716

These jobs run on preemptible VMs, and Hail Batch reruns a preempted job from the start on a new VM. A
10-sample pilot at cpu=4 without checkpointing was preempted at least 13 times in 3 hours and finished
nothing. So each job runs ExpansionHunter with --resume and copies its resume files to a per-job folder
under --checkpoint-dir every --checkpoint-interval-seconds; an attempt that follows a preemption restores
them and continues from the last finished locus (see resume_checkpoint_dir in
create_expansion_hunter_steps). A preemption then costs the CRAM copy, the decompression scan and at most
one checkpoint interval of genotyping.

--sample-table-path also accepts short_read_samples_with_truth_data.tsv, the table of the 135 short-read samples
for the genotype-quality model (keyed by sample_label, so HG002 at 10x/20x/31x are separate samples), copied from
eh_on_modal/genotype_quality_model_runs.tsv on 2026-10-01 and given a truth VCF column. Results for such a table
land in {OUTPUT_DIR}/{sample_label}/json/. For the genotype-quality model samples on the TRExplorer v2.1 catalog:

    python3 run_expansion_hunter_on_selected_samples.py --no-wait \\
        --sample-table-path short_read_samples_with_truth_data.tsv \\
        --catalog-path gs://tandem-repeat-catalog/v2.1/TRExplorer.repeat_catalog_v2.1.hg38.1_to_1000bp_motifs.EH.json.gz \\
        --output-dir gs://tandem-repeat-explorer/tool_genotype_quality/expansion_hunter_v2.1

Memory is highmem (13GB at cpu=2) because peak RSS on the whole catalog was ~10.7GB in the
2026-08-01 benchmark.

Sample sex is read from the sample table's Gender column ("male"/"female") and passed through to
--sex, per-sample -- required for correct hemizygous calling on chrX/chrY loci. select_1kGP_samples.py
takes it from the assembly rather than from the 1kGP metadata sheet, which has it wrong for HG02300.

Results land at {OUTPUT_DIR}/{sample_id}/json/. The file name comes from the CATALOG, not the sample:
create_expansion_hunter_steps names step 1's output after the catalog file, and its output_prefix
argument only names the combined TSV/BED that step 2 would have written. Its prefix strips a trailing
".json" but not ".json.gz", so each sample directory holds one "<catalog file name>.json.gz".

Since that name records nothing about what produced the file, the run writes {OUTPUT_DIR}/eh_run.json
holding the ExpansionHunter image, analysis mode, motif composition and --max-depth settings, catalog path and
content hash, and reference, and
refuses to add results to an output directory whose record disagrees with the current run.

create_expansion_hunter_steps' second step, which flattens the genotyping JSON into TSV/BED, is
cancelled rather than run. It is hardcoded to cpu=4/highmem (26GB), sized for a ~1.6M-locus catalog,
and OOM-kills on a catalog this size (confirmed in the 2026-08-01 cost benchmark). Cancelling it costs
nothing here: step 1 writes the per-sample JSON to GCS on its own, and the analysis downstream of this
reads that JSON directly. Delete the .skip() call if the TSV/BED conversion is wanted, but raise that
step's memory first.

Run from anywhere: the hail_batch_pipelines directory is added to sys.path so expansion_hunter_pipeline
can be imported, and every other path this script uses is absolute or resolved against SCRIPT_DIR.
"""
import json
import os
import re
import subprocess
import sys
import hailtop.fs as hfs
import pandas as pd
from step_pipeline import pipeline, Backend

HAIL_BATCH_PIPELINES_DIR = os.path.expanduser(
    "~/code/str-truth-set-v2/str-truth-set/tool_comparison/hail_batch_pipelines")
sys.path.append(HAIL_BATCH_PIPELINES_DIR)
from expansion_hunter_pipeline import create_expansion_hunter_steps, DOCKER_IMAGE, REFERENCE_FASTA_PATH, \
    REFERENCE_FASTA_FAI_PATH

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SAMPLE_TABLE_PATH = os.path.join(SCRIPT_DIR, "selected_1kGP_samples.tsv")
CATALOG_PATH = ("gs://tandem-repeat-catalog/v2.0/"
                "TR_catalog.TRExplorer-v2.with_extended_definitions.5657854_loci.ExpansionHunter.json.gz")
OUTPUT_DIR = "gs://tandem-repeat-explorer/tool_genotype_quality/expansion_hunter"
CHECKPOINT_DIR = "gs://tandem-repeat-explorer/tool_genotype_quality/expansion_hunter_checkpoints"

ANALYSIS_MODE = "optimized-streaming"
CPU = 2
THREADS = 4
MEMORY = "highmem"
CHECKPOINT_INTERVAL_SECONDS = 600
# Motif composition counts are written for every locus that has a record. Turning them on raises EH's default
# --max-depth from 150 to 500, which changes the calls at high-coverage loci; MAX_DEPTH keeps it at 150 so the calls
# match a run without motif composition, which is what the genotype quality model scores in normal use.
OUTPUT_MOTIF_COMPOSITION = "all-loci"
MAX_DEPTH = 150

bp = pipeline("run_expansion_hunter_on_selected_samples", backend=Backend.HAIL_BATCH_SERVICE, config_file_path="~/.step_pipeline")
parser = bp.get_config_arg_parser()
parser.add_argument("--sample-table-path", default=SAMPLE_TABLE_PATH)
parser.add_argument("--catalog-path", default=CATALOG_PATH)
parser.add_argument("--output-dir", default=OUTPUT_DIR)
parser.add_argument("-s", "--sample-id", action="append", help="Process only this sample. Can be specified more than once.")
parser.add_argument("--checkpoint-dir", default=CHECKPOINT_DIR,
                    help="Where each job keeps its ExpansionHunter --resume files, so an attempt that follows a "
                         "preemption continues instead of starting over.")
parser.add_argument("--checkpoint-interval-seconds", type=int, default=CHECKPOINT_INTERVAL_SECONDS,
                    help="How often each job copies its --resume files to --checkpoint-dir.")
parser.add_argument("--cpu", type=int, default=CPU,
                    help="Cores per job. Memory is highmem, 6.5GB per core, so pass 4 (26GB) to rerun a sample that "
                         "ran out of memory at the default. Not part of the checkpoint key, so the rerun resumes.")
parser.add_argument("--no-resume", action="store_true",
                    help="Do not checkpoint, and do not pass --resume to ExpansionHunter.")
args = bp.parse_known_args()

df = pd.read_table(args.sample_table_path)
if "sample_label" in df.columns:
    # short_read_samples_with_truth_data.tsv: one row per short-read sample, keyed by sample_label, so a genome
    # sequenced at several depths (HG002 at 10x/20x/31x) has one row per depth.
    df = df.rename(columns={"sex": "Gender", "reads_path": "cram_path", "reads_index_path": "crai_path"})
else:
    df["sample_label"] = df["sample_id"]
assert df.sample_label.is_unique, f"{args.sample_table_path} has duplicate sample labels"
if args.sample_id:
    # A mistyped -s would otherwise just drop out of the filter, and the run would submit a smaller
    # set of samples than was asked for without saying so.
    missing = set(args.sample_id) - set(df.sample_label)
    assert not missing, f"sample id(s) not in {args.sample_table_path}: {sorted(missing)}"
    df = df[df.sample_label.isin(args.sample_id)]

if len(df) == 0:
    parser.error(f"{args.sample_table_path} lists no samples, so there is nothing to genotype. "
                 f"Regenerate it with select_1kGP_samples.py, or pass --sample-table-path.")
assert set(df["Gender"].unique()) <= {"male", "female"}, f"unexpected Gender values: {df['Gender'].unique()}"

# Result files are named after the catalog, and step_pipeline skips a sample whose output already
# exists, so nothing in the output path records which ExpansionHunter build, analysis mode, catalog or
# reference produced a file. Recording that once per output directory turns "re-ran after the image
# changed" from a silent mix of old and new genotypes into an error.
#
# The catalog is recorded by content, not just by path: run_convert_extended_catalog_to_eh_catalog.sh
# overwrites one fixed object name whose only varying part is the locus count, so a re-converted
# catalog with the same count would otherwise pass this check and reuse the earlier results. crc32c
# rather than md5, because gsutil uploads a catalog this size as a composite object, and GCS reports
# no md5 for those.
catalog_stat = subprocess.run(["gsutil", "stat", args.catalog_path], capture_output=True, text=True, check=True).stdout
catalog_crc32c = re.search(r"Hash \(crc32c\):\s*(\S+)", catalog_stat)
assert catalog_crc32c, f"gsutil stat {args.catalog_path} reported no crc32c:\n{catalog_stat}"
run_record = {
    "docker_image": DOCKER_IMAGE,
    "analysis_mode": ANALYSIS_MODE,
    "output_motif_composition": OUTPUT_MOTIF_COMPOSITION,
    "max_depth": MAX_DEPTH,
    "catalog_path": args.catalog_path,
    "catalog_crc32c": catalog_crc32c.group(1),
    "reference_fasta": REFERENCE_FASTA_PATH,
}
run_record_path = os.path.join(args.output_dir, "eh_run.json")
if hfs.exists(run_record_path):
    with hfs.open(run_record_path, "r") as f:
        previous_run_record = json.load(f)
    assert previous_run_record == run_record, (
        f"{run_record_path} says the results already in {args.output_dir} were produced with "
        f"{previous_run_record}, which differs from this run's {run_record}. Point --output-dir "
        f"somewhere else rather than mixing the two.")

for _, row in df.iterrows():
    combine_step, _ = create_expansion_hunter_steps(
        bp,
        reference_fasta=REFERENCE_FASTA_PATH,
        reference_fasta_fai=REFERENCE_FASTA_FAI_PATH,
        input_bam=row.cram_path,
        input_bai=row.crai_path,
        male_or_female=row.Gender,
        variant_catalog_file_paths=[args.catalog_path],
        output_dir=os.path.join(args.output_dir, row.sample_label),
        output_prefix=f"{row.sample_label}.EHv5-bw2-optimized",
        analysis_mode=ANALYSIS_MODE,
        loci_to_exclude=None,
        min_locus_coverage=None,
        use_illumina_expansion_hunter=False,
        catalog_prefilter_step=None,
        num_shards=1,
        streaming_cpu=args.cpu,
        streaming_threads=THREADS,
        streaming_memory=MEMORY,
        resume_checkpoint_dir=None if args.no_resume else args.checkpoint_dir,
        checkpoint_interval_seconds=args.checkpoint_interval_seconds,
        output_motif_composition=OUTPUT_MOTIF_COMPOSITION,
        max_depth=MAX_DEPTH)
    combine_step.skip()

result = bp.run()
batch_id = getattr(result, "id", None)
# No batch id means step_pipeline transferred no steps: --dry-run, every sample's output already
# exists, or one of its own --skip-* flags was passed. Printing a URL in that case would name a batch
# that does not exist (.../batches/None), and naming one of those causes would be a guess. len(df) is
# what was SELECTED, never what was submitted, since a sample that already has output never becomes a
# job -- the per-step lines above are what say which ones ran.
if args.dry_run:
    print(f"Dry run, nothing submitted. {len(df)} samples selected.")
elif batch_id is None:
    print(f"No jobs submitted, out of {len(df)} samples selected. See the per-step lines above for why.")
else:
    print(f"Submitted batch: https://batch.hail.is/batches/{batch_id}")
    print(f"{len(df)} samples selected -> {args.output_dir}")

# Written only now, so a run that dies before submitting anything (bp.run() does a full parse_args and
# exits on an unrecognized flag, which bp.parse_known_args() above lets through) does not leave a
# record claiming results that were never produced.
if batch_id is not None and not hfs.exists(run_record_path):
    with hfs.open(run_record_path, "w") as f:
        json.dump(run_record, f, indent=4)
