"""Run trgt-lps on 256 HPRC TRGT VCF files, merge results into wide format,
and compute per-locus statistics for BigQuery loading.

Step 1: Per-sample trgt-lps (256 parallel jobs)
Step 2: Merge per-sample files into wide format (1 job)
Step 3: Compute per-locus statistics with population/sex stratification (1 job)
"""

import logging
import os
import re
import hailtop.fs as hfs
from step_pipeline import pipeline, Backend, Localize, Delocalize


DOCKER_IMAGE = "us-central1-docker.pkg.dev/cmg-analysis/docker-repo/str-analysis-with-trgt@sha256:75798915d2ae2f235ab9da470f718a7a7e849c2807f5837af99e4e4658f0334d"

GCS_BASE_DIR = "gs://tandem-repeat-catalog/v2.0/HPRC256_TRGT_vcfs_and_lps.2026_02"


CONVERT_SCRIPT = r'''#!/usr/bin/env python3
"""Compute per-locus LPS statistics from wide-format LPS table."""

import argparse
import collections
import gzip
import numpy as np
import pandas as pd

HEADER_FIELDS = [
    "locus_id", "motif", "allele_size_histogram", "biallelic_histogram",
    "min_allele", "mode_allele", "mean", "stdev", "median",
    "99th_percentile", "max_allele", "unique_allele_lengths", "num_called_alleles",
]

def write_line(locus_id, motif, allele_sizes, alleles_by_sample_id, output_file):
    if not allele_sizes:
        return
    if "," in locus_id:
        for specific_locus_id in locus_id.split(","):
            if specific_locus_id.endswith(f"-{motif}"):
                locus_id = specific_locus_id
                break
        else:
            raise ValueError(f"Couldn't resolve locus id for motif {motif} in {locus_id}")

    allele_sizes = sorted(allele_sizes)
    allele_counts = collections.Counter(allele_sizes)
    mode_allele, _ = allele_counts.most_common(1)[0]

    genotype_counts = collections.defaultdict(int)
    for sample_id, allele_list in alleles_by_sample_id.items():
        if len(allele_list) == 1:
            allele_list = allele_list * 2
        elif len(allele_list) != 2:
            raise ValueError(f"Found {len(allele_list)} alleles for {sample_id} in {locus_id} {motif}")
        genotype_counts[tuple(sorted(allele_list))] += 1

    allele_size_histogram = ",".join(f"{a}x:{c}" for a, c in sorted(allele_counts.items()))
    biallelic_histogram = ",".join(
        f"{g[0]}/{g[1]}:{c}" for g, c in sorted(genotype_counts.items(), key=lambda x: (x[0][1], x[0][0]))
    )
    output_file.write("\t".join(map(str, [
        locus_id, motif, allele_size_histogram, biallelic_histogram,
        int(min(allele_sizes)), mode_allele, f"{np.mean(allele_sizes):.3f}",
        f"{np.std(allele_sizes):.3f}", int(np.median(allele_sizes)),
        int(np.percentile(allele_sizes, 99)), int(max(allele_sizes)),
        len(set(allele_sizes)), len(allele_sizes),
    ])) + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-table", required=True)
    parser.add_argument("--sample-metadata-tsv", required=True)
    parser.add_argument("--population", choices=["AFR", "AMR", "EAS", "EUR", "SAS"])
    parser.add_argument("--sex", choices=["male", "female"])
    parser.add_argument("--output-prefix", default="hprc_lps.2026_02.per_locus_and_motif")
    args = parser.parse_args()

    fopen = gzip.open if args.input_table.endswith("gz") else open
    with fopen(args.input_table, "rt") as infile:
        header_fields = next(infile).strip().split("\t")
    sample_ids_in_input = header_fields[2:]
    print(f"Parsed {len(sample_ids_in_input):,d} sample ids from {args.input_table}")

    df_metadata = pd.read_table(args.sample_metadata_tsv)
    df_metadata = df_metadata[df_metadata.SampleId.isin(set(sample_ids_in_input))]
    if args.population:
        df_metadata = df_metadata[df_metadata.Population == args.population]
        print(f"Kept {len(df_metadata):,d} samples from population {args.population}")
    if args.sex:
        df_metadata = df_metadata[df_metadata.Sex == args.sex]
        print(f"Kept {len(df_metadata):,d} samples from sex {args.sex}")

    sample_ids_to_include = set(df_metadata.SampleId)

    output_path = args.output_prefix
    if args.population:
        output_path += f".only_{args.population}"
    if args.sex:
        output_path += f".only_{args.sex}"
    output_path += f".{len(sample_ids_to_include)}_samples.tsv.gz"

    print(f"Writing data from {len(sample_ids_to_include):,d} samples to {output_path}")
    with fopen(args.input_table, "rt") as infile, gzip.open(output_path, "wt") as outfile:
        next(infile)
        outfile.write("\t".join(HEADER_FIELDS) + "\n")
        for line_number, line in enumerate(infile):
            fields = line.strip().split("\t")
            locus_id, motif = fields[0], fields[1]
            alleles = []
            alleles_by_sample_id = collections.defaultdict(list)
            for sample_id, allele_sizes in zip(header_fields[2:], fields[2:]):
                if sample_id not in sample_ids_to_include or not allele_sizes:
                    continue
                for a in allele_sizes.split(","):
                    alleles.append(int(a))
                    alleles_by_sample_id[sample_id].append(int(a))
            write_line(locus_id, motif, alleles, alleles_by_sample_id, outfile)
    print(f"Wrote {line_number + 1:9,d} lines to {output_path}")

if __name__ == "__main__":
    main()
'''

MERGE_SCRIPT = r'''#!/usr/bin/env python3
"""Merge per-sample 3-col LPS files into wide format."""

import gzip
import collections
import sys

files = sys.argv[1:]
sample_ids = []
data = collections.OrderedDict()

for filepath in files:
    with gzip.open(filepath, "rt") as f:
        header = next(f).strip().split("\t")
        sample_id = header[2]
        sample_ids.append(sample_id)
        sample_idx = len(sample_ids) - 1
        for line in f:
            fields = line.strip().split("\t")
            key = (fields[0], fields[1])
            if key not in data:
                data[key] = [""] * len(files)
            data[key][sample_idx] = fields[2]
    print(f"Loaded {filepath}: {sample_id}", file=sys.stderr)

print(f"Writing {len(data)} rows for {len(sample_ids)} samples", file=sys.stderr)
with gzip.open("hprc-lps.txt.gz", "wt") as f:
    f.write("trid\tmotif\t" + "\t".join(sample_ids) + "\n")
    for (trid, motif), values in data.items():
        f.write(trid + "\t" + motif + "\t" + "\t".join(values) + "\n")
print("Done", file=sys.stderr)
'''


def main():
    bp = pipeline(backend=Backend.HAIL_BATCH_SERVICE, config_file_path="~/.step_pipeline")

    parser = bp.get_config_arg_parser()
    parser.add_argument("-n", "--num-samples", type=int, help="Only process first n samples (for testing)")
    args = bp.parse_known_args()

    # Upload metadata and helper scripts to GCS
    #metadata_local = os.path.join(os.path.dirname(__file__), "..", "hprc-lps_2025-12-06", "1kGP_metadata.tsv")
    #metadata_gcs = f"{GCS_BASE_DIR}/1kGP_metadata.tsv"
    #if not hfs.exists(metadata_gcs):
    #    print(f"Uploading {metadata_local} to {metadata_gcs}")
    #    hfs.copy(os.path.abspath(metadata_local), metadata_gcs)

    vcf_paths_precached = bp.precache_file_paths(os.path.join(GCS_BASE_DIR, "vcfs/*.vcf.gz"))

    sample_ids = []
    vcf_paths = []
    for vcf_path_dict in vcf_paths_precached:
        vcf_path = vcf_path_dict["path"]
        sample_id = re.sub(".vcf.gz$", "", os.path.basename(vcf_path))
        sample_ids.append(sample_id)
        vcf_paths.append(vcf_path)
        if args.num_samples and len(sample_ids) >= args.num_samples:
            break

    bp.set_name(f"TRGT-LPS: {len(vcf_paths)} HPRC samples")

    lps_output_dir = os.path.join(GCS_BASE_DIR, "lps")
    if not args.force:
        precached_paths = bp.precache_file_paths(os.path.join(lps_output_dir, "*.*"))
        logging.info(f"Precached {len(precached_paths):,d} output files")

    # Step 1: Per-sample trgt-lps

    step1_jobs = []
    for i, (sample_id, vcf_path) in enumerate(zip(sample_ids, vcf_paths)):

        s1 = bp.new_step(
            f"trgt-lps sample #{i+1}: {sample_id}",
            arg_suffix="lps",
            step_number=1,
            image=DOCKER_IMAGE,
            cpu=0.5,
            storage="30Gi",
            all_outputs_precached=True,
            localize_by=Localize.GSUTIL_COPY,
            delocalize_by=Delocalize.GSUTIL_COPY,
            output_dir=lps_output_dir,
        )

        s1.command("set -euxo pipefail")
        s1.switch_gcloud_auth_to_user_account()
        s1.command("cd /io/")

        local_vcf = s1.input(vcf_path)
        s1.command(f"trgt-lps --vcf {local_vcf} | bgzip > {sample_id}.lps.txt.gz")

        s1.output(f"/io/{sample_id}.lps.txt.gz")
        step1_jobs.append(s1)

    # Step 2: Merge per-sample files into wide format
    s2 = bp.new_step(
        f"Merge {len(sample_ids)} LPS files",
        arg_suffix="merge",
        step_number=2,
        image=DOCKER_IMAGE,
        cpu=2,
        memory="highmem",
        storage="100Gi",
        all_outputs_precached=True,
        localize_by=Localize.GSUTIL_COPY,
        delocalize_by=Delocalize.GSUTIL_COPY,
        output_dir=GCS_BASE_DIR,
    )
    for s1_job in step1_jobs:
        s2.depends_on(s1_job)

    s2.command("set -euxo pipefail")
    s2.switch_gcloud_auth_to_user_account()
    s2.command("cd /io/")

    lps_file_paths = []
    for sample_id in sample_ids:
        local_lps = s2.input(f"{lps_output_dir}/{sample_id}.lps.txt.gz")
        lps_file_paths.append(str(local_lps))

    s2.command(f"cat > merge_lps.py << 'PYEOF'\n{MERGE_SCRIPT}\nPYEOF")
    s2.command(f"python3 merge_lps.py {' '.join(lps_file_paths)}")
    s2.command("ls -lh hprc-lps.txt.gz")
    s2.output("/io/hprc-lps.txt.gz")

    bp.run()


if __name__ == "__main__":
    main()
