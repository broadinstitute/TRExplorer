"""Derive an assembly-based truth genotype for every locus definition in the catalog that
run_expansion_hunter_on_selected_samples.py genotyped, for each of the samples in
selected_1kGP_samples.tsv.

Truth comes from each sample's telomere-to-telomere assembly: the assembly was aligned to hg38 and its
variants called by DipCall, then restricted to DipCall's high-confidence regions (that step ran
previously; its output is the high_confidence_vcf_path column of the sample table). This step reads
those variant calls and works out how many repeat copies each haplotype carries at each catalog
interval, so "truth" here is independent of any short-read caller, which is what makes it a fair
reference to score ExpansionHunter against.

Both sides identify a locus the same way. filter_vcf_to_tandem_repeats writes
LocusId = "{chrom}-{start_0based}-{end_1based}-{motif}", which is exactly what
convert_bed_to_expansion_hunter_catalog assigned when the ExpansionHunter catalog was built from the
same BED, so the two tables join on LocusId with nothing to reconstruct.

One Hail Batch job per sample, mirroring run_expansion_hunter_on_selected_samples.py. The sibling
boundary_optimization/run_truth_genotyping.py runs this same tool locally, which was reasonable for
26 samples against 407,565 intervals; this catalog is 14x larger and there are 100 samples.

Measured locally on HG00438 against this catalog: 5.5 minutes, peak RSS 4.8GB, and a 147MB output
holding one row per catalog definition. cpu=2/highmem leaves room above that peak; the job is
single-threaded, so the second core is only there because highmem at 1 core is 6.5GB.

Note the output includes every locus, including those with no overlapping variant, where the tool
reports the reference length. Those are indistinguishable from loci DipCall could not call at all, so
the comparison step, not this one, drops loci outside each sample's high-confidence regions.

Usage:
    python3 run_truth_genotyping_on_selected_samples.py --no-wait
    python3 run_truth_genotyping_on_selected_samples.py -s HG01993 -s HG02293   # subset of samples

    # Truth for the genotype-quality model samples on the TRExplorer v2.1 catalog:
    python3 run_truth_genotyping_on_selected_samples.py --no-wait \\
        --sample-table-path short_read_samples_with_truth_data.tsv \\
        --catalog-bed-path gs://tandem-repeat-catalog/v2.1/TRExplorer.repeat_catalog_v2.1.hg38.1_to_1000bp_motifs.bed.gz \\
        --output-dir gs://str-truth-set-v2/tool_genotype_quality/truth_genotypes_v2.1
"""
import os

import pandas as pd
from step_pipeline import pipeline, Backend, Localize

# Built from str-analysis 2d6c4fd, which includes 54fd419: a non-repeat insertion on either allele now sets the
# whole locus to no call rather than dropping just that allele. The 10 samples genotyped on 2026-08-25 used
# sha256:0f6cd8ef, which predates the insertion filter entirely.
DOCKER_IMAGE = "weisburd/str-analysis-with-expansion-hunter@sha256:5990e80cd34ebf69e624c824b530504a03476d23ed0a421017f281587555a162"
REFERENCE_FASTA_PATH = "gs://str-truth-set/hg38/ref/hg38.fa"
REFERENCE_FASTA_FAI_PATH = "gs://str-truth-set/hg38/ref/hg38.fa.fai"

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SAMPLE_TABLE_PATH = os.path.join(SCRIPT_DIR, "selected_1kGP_samples.tsv")
CATALOG_BED_PATH = ("gs://tandem-repeat-catalog/v2.0/"
                    "TRExplorer.repeat_catalog_v2.hg38.1_to_1000bp_motifs.with_extended_definitions.bed.gz")
OUTPUT_DIR = "gs://str-truth-set-v2/tool_genotype_quality/truth_genotypes"

CPU = 2
MEMORY = "highmem"
STORAGE = "30Gi"

bp = pipeline("run_truth_genotyping_on_selected_samples", backend=Backend.HAIL_BATCH_SERVICE,
              config_file_path="~/.step_pipeline")
parser = bp.get_config_arg_parser()
parser.add_argument("--sample-table-path", default=SAMPLE_TABLE_PATH)
parser.add_argument("--catalog-bed-path", default=CATALOG_BED_PATH)
parser.add_argument("--output-dir", default=OUTPUT_DIR)
parser.add_argument("-s", "--sample-id", action="append",
                    help="Process only this sample. Can be specified more than once.")
args = bp.parse_known_args()

df = pd.read_table(args.sample_table_path)
# short_read_samples_with_truth_data.tsv has one row per short-read sample, so a genome sequenced at several
# depths (HG002 at 10x/20x/31x) appears more than once. Its truth comes from the assembly, not the reads, so it
# is genotyped once.
assert df.groupby("sample_id").high_confidence_vcf_path.nunique().max() == 1, (
    f"{args.sample_table_path} lists different truth VCFs for the same sample_id")
df = df.drop_duplicates("sample_id")
if args.sample_id:
    missing = set(args.sample_id) - set(df.sample_id)
    assert not missing, f"sample id(s) not in {args.sample_table_path}: {sorted(missing)}"
    df = df[df.sample_id.isin(args.sample_id)]

if len(df) == 0:
    parser.error(f"{args.sample_table_path} lists no samples, so there is nothing to genotype.")

for _, row in df.iterrows():
    s1 = bp.new_step(
        f"Truth genotypes for {row.sample_id}",
        arg_suffix="truth-genotyping-step",
        step_number=1,
        image=DOCKER_IMAGE,
        cpu=CPU,
        memory=MEMORY,
        storage=STORAGE,
        localize_by=Localize.GSUTIL_COPY,
        output_dir=os.path.join(args.output_dir, row.sample_id),
    )

    local_fasta = s1.input(REFERENCE_FASTA_PATH)
    s1.input(REFERENCE_FASTA_FAI_PATH)
    local_catalog = s1.input(args.catalog_bed_path)
    s1.input(args.catalog_bed_path + ".tbi")
    local_vcf = s1.input(row.high_confidence_vcf_path)
    s1.input(row.high_confidence_vcf_path + ".tbi")

    s1.command("set -euxo pipefail")
    s1.command(f"""/usr/bin/time --verbose python3 -m str_analysis.filter_vcf_to_tandem_repeats genotype \\
        --reference-fasta-path {local_fasta} \\
        --catalog-bed {local_catalog} \\
        --output-prefix {row.sample_id} \\
        {local_vcf}""")
    s1.command("ls -lhrt")

    s1.output(f"{row.sample_id}.tandem_repeat_genotypes.tsv.gz")

result = bp.run()
batch_id = getattr(result, "id", None)
if args.dry_run:
    print(f"Dry run, nothing submitted. {len(df)} samples selected.")
elif batch_id is None:
    print(f"No jobs submitted, out of {len(df)} samples selected. See the per-step lines above for why.")
else:
    print(f"Submitted batch: https://batch.hail.is/batches/{batch_id}")
    print(f"{len(df)} samples selected -> {args.output_dir}")
