"""Cost benchmark: EHv5-bw2-optimized genotyping the full TRExplorer-v2 catalog (or a pilot subset) at
several cpu/thread configs, on the 3 standard WGS replicate samples (HG00438, NA12878, HG00733).

Bypasses run_genotyping_tools.py's main() on purpose: that function hardcodes streaming_cpu=2,
streaming_threads=4 when calling create_expansion_hunter_steps ("the June cost benchmark's balanced
optimum ... Baked in code, not env vars" -- see its comment above that call), so the EHV5_STREAMING_CPU/
THREADS/MEMORY env vars the older 2026-06-25 benchmark scripts relied on no longer have any effect there.
This script calls create_expansion_hunter_steps directly so streaming_cpu/threads/memory can actually vary,
without touching the shared production script.

Run from anywhere (cd's into the hail_batch_pipelines dir itself so its relative imports resolve).
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
    "~/code/str-truth-set-v2/run_tools/HPRC_all_aligned_short_read_and_long_read_samples.tsv")
OUTPUT_DIR = "gs://str-truth-set-v2/benchmark_cost__2026_08_01/tool_results"
SAMPLES = ["HG00438", "NA12878", "HG00733"]

bp = pipeline("benchmark_full_catalog_configs", backend=Backend.HAIL_BATCH_SERVICE, config_file_path="~/.step_pipeline")
parser = bp.get_config_arg_parser()
parser.add_argument("--catalog-path", required=True, help="gs:// path to the ExpansionHunter variant catalog")
parser.add_argument("--label", required=True, help="short label for this catalog (used in the output subdir)")
parser.add_argument("--cpu", type=int, required=True)
parser.add_argument("--threads", type=int, required=True)
parser.add_argument("--memory", default="highmem")
args = bp.parse_known_args()

df = pd.read_table(SAMPLE_TABLE_PATH)
df = df[(df.sample_id.isin(SAMPLES)) & (df.sequencing_data_type == "illumina")]
assert len(df) == len(SAMPLES), f"expected {len(SAMPLES)} illumina rows, got {len(df)}:\n{df}"

subdir = f"{args.label}/cpu{args.cpu}_t{args.threads}"
for _, row in df.iterrows():
    coverage = int(round(float(row.depth_of_coverage)))
    output_dir = os.path.join(OUTPUT_DIR, row.sample_id, "illumina", "EHv5-bw2-optimized", f"{coverage}x_coverage", subdir)
    create_expansion_hunter_steps(
        bp,
        reference_fasta=REFERENCE_FASTA_PATH,
        reference_fasta_fai=REFERENCE_FASTA_FAI_PATH,
        input_bam=row.read_data_path,
        input_bai=row.read_data_index_path,
        male_or_female=row.male_or_female,
        variant_catalog_file_paths=[args.catalog_path],
        output_dir=output_dir,
        output_prefix=f"{row.sample_id}.EHv5-bw2-optimized",
        analysis_mode="optimized-streaming",
        loci_to_exclude=None,
        min_locus_coverage=None,
        use_illumina_expansion_hunter=False,
        catalog_prefilter_step=None,
        num_shards=1,
        streaming_cpu=args.cpu,
        streaming_threads=args.threads,
        streaming_memory=args.memory)

result = bp.run()
batch_id = getattr(result, "id", None)
print(f"Submitted batch: https://batch.hail.is/batches/{batch_id}")
