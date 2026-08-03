"""Round 2: run ExpansionHunter-bw2 on the motif-reanchoring-rule validation catalog (12 loci x
baseline/reanchored, built by build_reanchor_validation_catalog_2.py) against real HG002 short reads
via Hail Batch. Same pattern as run_hail_validation.py (round 1).

Usage: python3 run_hail_validation_2.py --output-dir gs://bw2-delete-after-5-days/boundary_optimization_validation
"""

import os

from step_pipeline import pipeline, Backend, Localize

CATALOG_PATH = "gs://bw2-delete-after-5-days/boundary_optimization_validation/reanchor_validation_catalog_2.json"
REFERENCE_FASTA = "gs://str-truth-set/hg38/ref/hg38.fa"
REFERENCE_FAI = "gs://str-truth-set/hg38/ref/hg38.fa.fai"
HG002_CRAM = "gs://str-truth-set-v2/raw_data/HG002/illumina/HG002.pcr_free.cram"
HG002_CRAI = "gs://str-truth-set-v2/raw_data/HG002/illumina/HG002.pcr_free.cram.crai"

DOCKER_IMAGE = "docker.io/weisburd/expansion-hunter:latest"


def main():
    bp = pipeline(backend=Backend.HAIL_BATCH_SERVICE, config_file_path="~/.step_pipeline")

    parser = bp.get_config_arg_parser()
    parser.add_argument("--output-dir", required=True)
    args = bp.parse_known_args()

    bp.set_name("Motif-reanchoring rule validation round 2: EH on 24 catalog variants (12 loci) vs HG002")

    if not args.force:
        bp.precache_file_paths(os.path.join(args.output_dir, "**/*.*"))

    s1 = bp.new_step(
        name="ExpansionHunter on motif-reanchoring validation catalog (HG002)",
        step_number=1,
        arg_suffix="eh_hg002_round2",
        image=DOCKER_IMAGE,
        cpu=2,
        memory="standard",
        storage="20Gi",
        output_dir=args.output_dir,
    )

    local_ref = s1.input(REFERENCE_FASTA, localize_by=Localize.COPY)
    s1.input(REFERENCE_FAI, localize_by=Localize.COPY)
    local_cram = s1.input(HG002_CRAM, localize_by=Localize.HAIL_BATCH_CLOUDFUSE)
    s1.input(HG002_CRAI, localize_by=Localize.HAIL_BATCH_CLOUDFUSE)
    local_catalog = s1.input(CATALOG_PATH, localize_by=Localize.COPY)

    s1.command("set -euxo pipefail")
    s1.command("cd /io/")
    s1.command(
        f"""/usr/bin/time --verbose /usr/bin/ExpansionHunter \\
    --reads {local_cram} \\
    --reference {local_ref} \\
    --catalog {local_catalog} \\
    --output-prefix HG002.reanchor_validation_round2 \\
    --sex male
"""
    )
    s1.command("ls -lhrt")
    s1.output("HG002.reanchor_validation_round2.json")
    s1.output("HG002.reanchor_validation_round2.vcf")

    bp.run()


if __name__ == "__main__":
    main()
