"""Genotype the TRExplorer-v2 catalog, with the gap-purity rule's extended locus definitions added
alongside the original ones, using ExpansionHunter-bw2 (EHv5-bw2-optimized, optimized-streaming) on
the samples in selected_1kGP_samples.tsv (see select_1kGP_samples.py).

The catalog holds 5,657,854 locus definitions: all 5,599,658 v2 loci plus the 58,196 extended
definitions that survived deduplication, built by run_extend_v2_catalog_boundaries.sh and then
run_convert_extended_catalog_to_eh_catalog.sh. Both definitions of an extended locus are genotyped in
the same run, against the same reads, so comparing them measures the boundary and nothing else. Each
LocusId is the definition's own "{chrom}-{start_0based}-{end_1based}-{motif}", so an extended
definition and the locus it grew from are linked by overlap rather than by a shared id.

Config is cpu=4/threads=8/highmem (26GiB). Measured by the 2026-08-01 cost benchmark on 3 replicate
WGS samples against the 5,594,988-locus catalog, which this one exceeds by 1.1%:

    cpu=1/threads=2   OOM before genotyping starts (highmem at 1 core is only 6.5GB)
    cpu=2/threads=4   5.77, 5.79, 5.87 hours   $0.482, $0.484, $0.490
    cpu=4/threads=8   1.96, 4.04, 4.44 hours   $0.365, $0.637, $0.716

cpu=2 is the cheaper config per sample, but these jobs run on preemptible VMs, and a ~5.8 hour job is
a long window in which to be preempted; a preempted job restarts from scratch. cpu=4 cuts the typical
runtime to ~3.5 hours for about 17% more money on average, which is the trade this run takes.

Memory is left at highmem out of caution rather than measurement: peak RSS was ~10.7GB, so standard
(4GiB/core = 16GiB at cpu=4) would also fit and would cost less. Change it only with a pilot sample.

Sample sex is read from the sample table's Gender column ("male"/"female") and passed through to
--sex, per-sample -- required for correct hemizygous calling on chrX/chrY loci. select_1kGP_samples.py
takes it from the assembly rather than from the 1kGP metadata sheet, which has it wrong for HG02300.

Results land at {OUTPUT_DIR}/{sample_id}/json/. The file name comes from the CATALOG, not the sample:
create_expansion_hunter_steps names step 1's outputs after the catalog file, and its output_prefix
argument only names the combined TSV/BED that step 2 would have written. Its prefix strips a trailing
".json" but not ".json.gz", so with this catalog every sample directory holds one
"TR_catalog.TRExplorer-v2.with_extended_definitions.5657854_loci.ExpansionHunter.json.gz.json.gz".

Since that name records nothing about what produced the file, the run writes {OUTPUT_DIR}/eh_run.json
holding the ExpansionHunter image, analysis mode, catalog path and content hash, and reference, and
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
OUTPUT_DIR = "gs://str-truth-set-v2/tool_genotype_quality/expansion_hunter"

ANALYSIS_MODE = "optimized-streaming"
CPU = 4
THREADS = 8
MEMORY = "highmem"

bp = pipeline("run_expansion_hunter_on_selected_samples", backend=Backend.HAIL_BATCH_SERVICE, config_file_path="~/.step_pipeline")
parser = bp.get_config_arg_parser()
parser.add_argument("--sample-table-path", default=SAMPLE_TABLE_PATH)
parser.add_argument("--catalog-path", default=CATALOG_PATH)
parser.add_argument("--output-dir", default=OUTPUT_DIR)
parser.add_argument("-s", "--sample-id", action="append", help="Process only this sample. Can be specified more than once.")
args = bp.parse_known_args()

df = pd.read_table(args.sample_table_path)
if args.sample_id:
    # A mistyped -s would otherwise just drop out of the filter, and the run would submit a smaller
    # set of samples than was asked for without saying so.
    missing = set(args.sample_id) - set(df.sample_id)
    assert not missing, f"sample id(s) not in {args.sample_table_path}: {sorted(missing)}"
    df = df[df.sample_id.isin(args.sample_id)]

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
        output_dir=os.path.join(args.output_dir, row.sample_id),
        output_prefix=f"{row.sample_id}.EHv5-bw2-optimized",
        analysis_mode=ANALYSIS_MODE,
        loci_to_exclude=None,
        min_locus_coverage=None,
        use_illumina_expansion_hunter=False,
        catalog_prefilter_step=None,
        num_shards=1,
        streaming_cpu=CPU,
        streaming_threads=THREADS,
        streaming_memory=MEMORY)
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
