"""Score each sample's ExpansionHunter genotypes against its assembly-derived truth genotypes, one
Hail Batch job per sample, by running compare_eh_to_truth.py in the cloud next to the data.

The comparison runs here rather than locally because each sample's ExpansionHunter JSON holds a record
per catalog locus and is hundreds of megabytes; pulling 100 of them down would move tens of gigabytes
to do arithmetic that produces a few megabytes per sample.

compare_eh_to_truth.py is uploaded to the cloud under a name carrying its own content hash, so a job
can never run a stale copy of it. Results are written flat, at {OUTPUT_DIR}/{sample}/, one directory
per sample.

Note what that gives up: step_pipeline skips a step whose output already exists, so a staged run made
after the comparison logic changed leaves the already-finished samples on the old logic and scores
them alongside samples computed with the new one. The same applies when the truth genotypes change,
since nothing in the output path reflects them. Pass --force-comparison-step, or delete the output
directory, whenever either side has changed.

Inputs per sample, all produced by earlier steps:
    {EH_OUTPUT_DIR}/{sample}/json/*.json.gz                  run_expansion_hunter_on_selected_samples.py,
                                                             which writes exactly one JSON per sample
    {TRUTH_DIR}/{sample}/{sample}.tandem_repeat_genotypes.tsv.gz
                                                             run_truth_genotyping_on_selected_samples.py
    the sample's high_confidence_bed_path                     DipCall, run previously

Usage:
    python3 run_comparison_on_selected_samples.py --no-wait
    python3 run_comparison_on_selected_samples.py -s HG01993   # subset of samples
"""
import hashlib
import os

import hailtop.fs as hfs
import pandas as pd
from step_pipeline import pipeline, Backend, Localize

DOCKER_IMAGE = "weisburd/str-analysis-with-expansion-hunter@sha256:0f6cd8efbae6b2c35837c856267347e94d0d86cdce60cd80881153fe2d0e57f7"

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
COMPARISON_SCRIPT_PATH = os.path.join(SCRIPT_DIR, "compare_eh_to_truth.py")
SAMPLE_TABLE_PATH = os.path.join(SCRIPT_DIR, "selected_1kGP_samples.tsv")
CATALOG_BED_PATH = ("gs://tandem-repeat-catalog/v2.0/"
                    "TRExplorer.repeat_catalog_v2.hg38.1_to_1000bp_motifs.with_extended_definitions.bed.gz")
EH_OUTPUT_DIR = "gs://str-truth-set-v2/tool_genotype_quality/expansion_hunter"
TRUTH_DIR = "gs://str-truth-set-v2/tool_genotype_quality/truth_genotypes"
SCRIPTS_DIR = "gs://str-truth-set-v2/tool_genotype_quality/scripts"
OUTPUT_DIR = "gs://str-truth-set-v2/tool_genotype_quality/comparisons"

CPU = 2
MEMORY = "highmem"
STORAGE = "20Gi"


def upload_comparison_script(dry_run):
    """Copy compare_eh_to_truth.py to the cloud under a hash-suffixed name, and return that path."""
    script_bytes = open(COMPARISON_SCRIPT_PATH, "rb").read()
    cloud_path = os.path.join(
        SCRIPTS_DIR, f"compare_eh_to_truth.{hashlib.sha256(script_bytes).hexdigest()[:12]}.py")
    if hfs.exists(cloud_path):
        print(f"{cloud_path} is already up to date")
    elif dry_run:
        print(f"would upload {COMPARISON_SCRIPT_PATH} to {cloud_path}")
    else:
        print(f"uploading {COMPARISON_SCRIPT_PATH} to {cloud_path}")
        with hfs.open(cloud_path, "wb") as f:
            f.write(script_bytes)
    return cloud_path


def find_eh_json_path(eh_output_dir, sample_id):
    """Return the sample's ExpansionHunter JSON, or None if the run has not written it yet.

    run_expansion_hunter_on_selected_samples.py names the file after the catalog, not the sample, so it is
    found by listing the sample's json/ directory. More than one JSON there means results from different
    runs were mixed in one directory, which is an error rather than something to pick from.
    """
    json_dir = os.path.join(eh_output_dir, sample_id, "json")
    if not hfs.exists(json_dir):
        return None
    paths = [entry.path for entry in hfs.ls(json_dir) if entry.path.endswith(".json.gz")]
    assert len(paths) <= 1, f"expected at most 1 ExpansionHunter JSON in {json_dir}, found {len(paths)}: {paths}"
    return paths[0] if paths else None


bp = pipeline("run_comparison_on_selected_samples", backend=Backend.HAIL_BATCH_SERVICE,
              config_file_path="~/.step_pipeline")
parser = bp.get_config_arg_parser()
parser.add_argument("--sample-table-path", default=SAMPLE_TABLE_PATH)
parser.add_argument("--catalog-bed-path", default=CATALOG_BED_PATH)
parser.add_argument("--eh-output-dir", default=EH_OUTPUT_DIR)
parser.add_argument("--truth-dir", default=TRUTH_DIR)
parser.add_argument("--output-dir", default=OUTPUT_DIR)
parser.add_argument("-s", "--sample-id", action="append",
                    help="Process only this sample. Can be specified more than once.")
args = bp.parse_known_args()

df = pd.read_table(args.sample_table_path)
if args.sample_id:
    missing = set(args.sample_id) - set(df.sample_id)
    assert not missing, f"sample id(s) not in {args.sample_table_path}: {sorted(missing)}"
    df = df[df.sample_id.isin(args.sample_id)]

if len(df) == 0:
    parser.error(f"{args.sample_table_path} lists no samples, so there is nothing to compare.")

comparison_script_cloud_path = upload_comparison_script(args.dry_run)

# A sample whose ExpansionHunter or truth output is not there yet is skipped rather than submitted:
# the job would only fail on a missing input several minutes in, and both of those steps are expected
# to be run in stages.
missing_inputs = []
for _, row in df.iterrows():
    eh_json_path = find_eh_json_path(args.eh_output_dir, row.sample_id)
    truth_tsv_path = os.path.join(args.truth_dir, row.sample_id,
                                  f"{row.sample_id}.tandem_repeat_genotypes.tsv.gz")
    if eh_json_path is None or not hfs.exists(truth_tsv_path):
        missing_inputs.append(row.sample_id)
        continue

    s1 = bp.new_step(
        f"Compare ExpansionHunter to truth for {row.sample_id}",
        arg_suffix="comparison-step",
        step_number=1,
        image=DOCKER_IMAGE,
        cpu=CPU,
        memory=MEMORY,
        storage=STORAGE,
        localize_by=Localize.GSUTIL_COPY,
        output_dir=os.path.join(args.output_dir, row.sample_id),
    )

    local_script = s1.input(comparison_script_cloud_path)
    local_catalog = s1.input(args.catalog_bed_path)
    local_eh_json = s1.input(eh_json_path)
    local_truth_tsv = s1.input(truth_tsv_path)
    local_high_confidence_bed = s1.input(row.high_confidence_bed_path)

    s1.command("set -euxo pipefail")
    s1.command(f"""/usr/bin/time --verbose python3 {local_script} \\
        --catalog-bed {local_catalog} \\
        --eh-json {local_eh_json} \\
        --truth-tsv {local_truth_tsv} \\
        --high-confidence-bed {local_high_confidence_bed} \\
        --output-npz {row.sample_id}.comparison.npz""")
    s1.command("ls -lhrt")

    s1.output(f"{row.sample_id}.comparison.npz")

if missing_inputs:
    print(f"skipping {len(missing_inputs)} samples that do not have both an ExpansionHunter and a "
          f"truth output yet: {', '.join(missing_inputs)}")

result = bp.run()
batch_id = getattr(result, "id", None)
if args.dry_run:
    print(f"Dry run, nothing submitted. {len(df) - len(missing_inputs)} samples ready to compare.")
elif batch_id is None:
    print(f"No jobs submitted, out of {len(df) - len(missing_inputs)} samples ready to compare. "
          f"See the per-step lines above for why.")
else:
    print(f"Submitted batch: https://batch.hail.is/batches/{batch_id}")
    print(f"{len(df) - len(missing_inputs)} samples -> {args.output_dir}/")
