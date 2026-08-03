"""Genotype the gap-purity-rule catalog (16,264 loci x original+extended = 32,528 entries,
gs://tandem-repeat-catalog/v2.0/boundary_optimization/boundary_optimization_gap_purity_extended_loci.EH_catalog.json,
built by build_eh_catalog_gap_purity_extended_loci.py) with ExpansionHunter-bw2 on the 1kGP samples
listed in data-prep/expansion_hunter_genotype_quality/selected_1kGP_samples.tsv.

--analysis-mode optimized-streaming, cpu=2/threads=4/highmem -- same config as the sibling
run_expansion_hunter_on_selected_samples.py, per instructions.

Sample sex is read from the sample table's Gender column, passed through to --sex per-sample.

Usage:
    python3 run_expansion_hunter_gap_purity_extended_loci.py --output-dir gs://bw2-delete-after-5-days/boundary_optimization_validation/gap_purity_eh_run
    python3 run_expansion_hunter_gap_purity_extended_loci.py --output-dir ... -s HG01884 -s HG01891   # pilot on specific samples
"""

import os
import sys

import pandas as pd
from step_pipeline import pipeline, Backend

HAIL_BATCH_PIPELINES_DIR = os.path.expanduser(
    "~/code/str-truth-set-v2/str-truth-set/tool_comparison/hail_batch_pipelines")
sys.path.append(HAIL_BATCH_PIPELINES_DIR)
from expansion_hunter_pipeline import create_expansion_hunter_steps, REFERENCE_FASTA_PATH, REFERENCE_FASTA_FAI_PATH

SAMPLE_TABLE_PATH = os.path.expanduser(
    "~/code/tandem-repeat-explorer/data-prep/expansion_hunter_genotype_quality/selected_1kGP_samples.tsv")
CATALOG_PATH = "gs://tandem-repeat-catalog/v2.0/boundary_optimization/boundary_optimization_gap_purity_extended_loci.EH_catalog.json"
OUTPUT_DIR = "gs://bw2-delete-after-5-days/boundary_optimization_validation/gap_purity_eh_run"

ANALYSIS_MODE = "optimized-streaming"
CPU = 2
THREADS = 4
MEMORY = "highmem"

bp = pipeline("run_expansion_hunter_gap_purity_extended_loci", backend=Backend.HAIL_BATCH_SERVICE, config_file_path="~/.step_pipeline")
parser = bp.get_config_arg_parser()
parser.add_argument("--sample-table-path", default=SAMPLE_TABLE_PATH)
parser.add_argument("--catalog-path", default=CATALOG_PATH)
parser.add_argument("--output-dir", default=OUTPUT_DIR)
parser.add_argument("-s", "--sample-id", action="append", help="Process only this sample. Can be specified more than once.")
args = bp.parse_known_args()

df = pd.read_table(args.sample_table_path)
if args.sample_id:
    df = df[df.sample_id.isin(args.sample_id)]

assert set(df["Gender"].unique()) <= {"male", "female"}, f"unexpected Gender values: {df['Gender'].unique()}"
assert len(df) > 0, "no samples selected"

bp.set_name(f"Boundary optimization: EHv5-bw2 ({ANALYSIS_MODE}) on {len(df)} sample(s), gap-purity-extended catalog")

for _, row in df.iterrows():
    create_expansion_hunter_steps(
        bp,
        reference_fasta=REFERENCE_FASTA_PATH,
        reference_fasta_fai=REFERENCE_FASTA_FAI_PATH,
        input_bam=row.cram_path,
        input_bai=row.crai_path,
        male_or_female=row.Gender,
        variant_catalog_file_paths=[args.catalog_path],
        output_dir=os.path.join(args.output_dir, row.sample_id),
        output_prefix=f"{row.sample_id}.gap_purity_extended_loci",
        analysis_mode=ANALYSIS_MODE,
        loci_to_exclude=None,
        min_locus_coverage=None,
        use_illumina_expansion_hunter=False,
        catalog_prefilter_step=None,
        num_shards=1,
        streaming_cpu=CPU,
        streaming_threads=THREADS,
        streaming_memory=MEMORY)

result = bp.run()
batch_id = getattr(result, "id", None)
print(f"Submitted batch: https://batch.hail.is/batches/{batch_id}")
print(f"Genotyped {len(df)} sample(s) -> {args.output_dir}")
