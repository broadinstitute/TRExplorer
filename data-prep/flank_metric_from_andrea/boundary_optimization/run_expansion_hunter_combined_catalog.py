"""Genotype the combined catalog (364,030 entries: all Andrea v2-exact-match loci -- extended and
not-extended -- plus all non-GCN known disease loci, built by build_combined_eh_catalog.py) with
ExpansionHunter-bw2 on the 1kGP samples in selected_1kGP_samples.tsv.

Single step per sample: runs ExpansionHunter (--analysis-mode optimized-streaming, cpu=4/threads=4/
highmem) AND the combine_str_json_to_tsv step in the same container/step (folded together, rather
than the sibling expansion_hunter_pipeline.py's separate step1/step2), so there's one job per sample
instead of two. All json/tsv outputs are gzip-compressed (EH via --compress-output-files, the combine
script's tsv outputs are gzip natively, the bed output is bgzipped explicitly).

Usage:
    python3 run_expansion_hunter_combined_catalog.py --output-dir gs://bw2-delete-after-5-days/boundary_optimization_validation/combined_catalog_eh_run
    python3 run_expansion_hunter_combined_catalog.py --output-dir ... -s HG01884 -s HG01891   # subset of samples
"""

import hashlib
import os

import hailtop.fs as hfs
import pandas as pd
from step_pipeline import pipeline, Backend, Localize

DOCKER_IMAGE = "weisburd/str-analysis-with-expansion-hunter@sha256:7a006e45a297d42a914fb5b36ee2866308d0a4c3f0a60360f56897a53a76f1ec"
REFERENCE_FASTA_PATH = "gs://str-truth-set/hg38/ref/hg38.fa"
REFERENCE_FASTA_FAI_PATH = "gs://str-truth-set/hg38/ref/hg38.fa.fai"

SAMPLE_TABLE_PATH = os.path.expanduser(
    "~/code/tandem-repeat-explorer/data-prep/expansion_hunter_genotype_quality/selected_1kGP_samples.tsv")
CATALOG_PATH = "gs://tandem-repeat-catalog/v2.0/boundary_optimization/combined_eh_catalog.json"
OUTPUT_DIR = "gs://bw2-delete-after-5-days/boundary_optimization_validation/combined_catalog_eh_run"

ANALYSIS_MODE = "optimized-streaming"
CPU = 4
THREADS = 4
MEMORY = "highmem"

bp = pipeline("run_expansion_hunter_combined_catalog", backend=Backend.HAIL_BATCH_SERVICE, config_file_path="~/.step_pipeline")
parser = bp.get_config_arg_parser()
parser.add_argument("--sample-table-path", default=SAMPLE_TABLE_PATH)
parser.add_argument("--catalog-path", default=CATALOG_PATH)
parser.add_argument("--output-dir", default=OUTPUT_DIR)
parser.add_argument("-s", "--sample-id", action="append", help="Process only this sample. Can be specified more than once.")
args = bp.parse_known_args()
bp.gcloud_project(args.gcloud_project)  # requester-pays project for the 1kGP CRAM bucket

df = pd.read_table(args.sample_table_path)
if args.sample_id:
    df = df[df.sample_id.isin(args.sample_id)]

assert set(df["Gender"].unique()) <= {"male", "female"}, f"unexpected Gender values: {df['Gender'].unique()}"
assert len(df) > 0, "no samples selected"

with hfs.open(args.catalog_path, "rb") as f:
    catalog_hash = hashlib.sha256(f.read()).hexdigest()[:12]
print(f"Catalog: {args.catalog_path} (sha256:{catalog_hash})")

bp.set_name(f"Combined catalog: EHv5-bw2 ({ANALYSIS_MODE}, single-step) on {len(df)} sample(s) [catalog {catalog_hash}]")

if not args.force:
    # precache_file_paths skips re-running any sample whose output already exists in --output-dir --
    # keyed only on sample_id, with no record of which catalog version produced that cached output.
    # If the catalog changed since a cached output was produced, that stale output is silently reused.
    print(f"NOTE: --force not passed -- samples with existing output in {args.output_dir} will be "
          f"skipped even if they were genotyped against a different catalog than the current one "
          f"(sha256:{catalog_hash}). Pass --force to regenotype everyone against the current catalog.")
    bp.precache_file_paths(os.path.join(args.output_dir, "**/*.*"))

for _, row in df.iterrows():
    output_prefix = f"{row.sample_id}.combined_catalog"

    cram_size_bytes = hfs.ls(row.cram_path)[0].size
    storage = f"{int(cram_size_bytes / 10**9) + 14}Gi"  # CRAM + hg38.fa (~3GB) + catalog + output, w/ margin

    s1 = bp.new_step(
        name=f"EHv5:{ANALYSIS_MODE}+combine on {row.sample_id}",
        step_number=1,
        arg_suffix="run-and-combine",
        image=DOCKER_IMAGE,
        cpu=CPU,
        memory=MEMORY,
        storage=storage,
        localize_by=Localize.GSUTIL_COPY,
        output_dir=args.output_dir,
    )

    local_fasta = s1.input(REFERENCE_FASTA_PATH)
    s1.input(REFERENCE_FASTA_FAI_PATH)
    local_bam = s1.input(row.cram_path)
    s1.input(row.crai_path)
    local_catalog = s1.input(args.catalog_path)

    s1.command("set -euxo pipefail")
    s1.command(f"""/usr/bin/time --verbose ExpansionHunter --threads {THREADS} --cache-mates \\
        --dont-output-consensus-sequences --output-genotype-timing --compress-output-files \\
        --reference {local_fasta} \\
        --reads {local_bam} \\
        --analysis-mode {ANALYSIS_MODE} \\
        --sex {row.Gender} \\
        --variant-catalog {local_catalog} \\
        --output-prefix {output_prefix}""")
    s1.command("ls -lhrt")

    # fold the combine step in here. combine_str_json_to_tsv accepts an explicit json_paths
    # positional arg -- pass the exact file rather than relying on its default cwd-glob discovery,
    # which would otherwise also pick up whatever other .json(.gz) files happen to be in the step's
    # working directory (reference/CRAM/catalog inputs, docker-image package data, etc.) and inflate
    # its "{prefix}.{n}_json_files...." output naming past n=1.
    s1.command(f"python3.9 -m str_analysis.combine_str_json_to_tsv {output_prefix}.json.gz "
               f"--include-extra-expansion-hunter-fields --output-prefix {output_prefix}")
    s1.command(f"bgzip {output_prefix}.1_json_files.bed")
    s1.command(f"tabix {output_prefix}.1_json_files.bed.gz")
    s1.command("ls -lhrt")

    s1.output(f"{output_prefix}.json.gz", output_dir=os.path.join(args.output_dir, "json"))
    s1.output(f"{output_prefix}.1_json_files.variants.tsv.gz")
    s1.output(f"{output_prefix}.1_json_files.alleles.tsv.gz")
    s1.output(f"{output_prefix}.1_json_files.bed.gz")
    s1.output(f"{output_prefix}.1_json_files.bed.gz.tbi")

result = bp.run()
batch_id = getattr(result, "id", None)
print(f"Submitted batch: https://batch.hail.is/batches/{batch_id}")
print(f"Genotyped {len(df)} sample(s) -> {args.output_dir}")
