"""Select the 50 male and 50 female 1kGP samples with the most autosomal sequence inside their
assembly's DipCall high-confidence regions.

A sample qualifies only if all three pieces of the comparison exist:
  - a short-read CRAM (the 698-sample Broad 1kGP list), so ExpansionHunter can genotype it
  - DipCall output from a de novo assembly, which defines the high-confidence regions
  - filter_vcf_to_tandem_repeats output, i.e. the variant calls a truth genotype is derived from

DipCall was run in three batches that write to different directories under
gs://str-truth-set-v2/dipcall_pipeline/ (the original HPRC-release-1 + HGSVC2 assemblies at the top
level, HPRC release 2 under HPRC_release2/, and the OpenHGL assemblies under human579_assemblies/),
so which samples have truth data is read from the outputs themselves rather than from any one
assembly table. As of 2026-08-24 that is 138 samples: 70 male, 68 female.

Samples are ranked by total autosomal high-confidence bases, so the ones picked are those where the
assembly can adjudicate the most of the genome. chrX and chrY are left out of the ranking because
their coverage is a function of the sample's sex and of which assembly batch produced it, not of
assembly quality: a male's chrX is haploid, and three of the HGSVC2 male assemblies have essentially
no high-confidence chrX at all.

Inputs:
  ~/code/str-truth-set-v2/20130606_sample_info_1kGP.tsv                            -- 1kGP phase-3 sample metadata
  ~/code/str-truth-set-v2/run_tools/broad_short_read_cram_paths_for_all_1kGP.txt   -- the 698-sample short-read CRAM list
  gs://str-truth-set-v2/dipcall_pipeline/**/{sample}.dip.bed.gz                    -- high-confidence regions
  gs://str-truth-set-v2/filter_vcf_v2__2025_12_29/{sample}/                        -- truth variant calls

Output:
  selected_1kGP_samples.tsv (in this script's own directory)

Usage:
    python3 select_1kGP_samples.py
    python3 select_1kGP_samples.py --refresh   # re-list the cloud directories instead of using the cache
"""
import argparse
import gzip
import subprocess
from pathlib import Path

import pandas as pd

SCRIPT_DIR = Path(__file__).parent
SAMPLE_INFO_PATH = Path("~/code/str-truth-set-v2/20130606_sample_info_1kGP.tsv").expanduser()
CRAM_LIST_PATH = Path("~/code/str-truth-set-v2/run_tools/broad_short_read_cram_paths_for_all_1kGP.txt").expanduser()
DIPCALL_DIR = "gs://str-truth-set-v2/dipcall_pipeline"
DIPCALL_SUBDIRS = ("", "HPRC_release2/", "human579_assemblies/")
TRUTH_VCF_DIR = "gs://str-truth-set-v2/filter_vcf_v2__2025_12_29"
CACHE_DIR = SCRIPT_DIR / "data"
BED_CACHE_DIR = CACHE_DIR / "high_confidence_beds"
OUTPUT_PATH = SCRIPT_DIR / "selected_1kGP_samples.tsv"

SAMPLES_PER_SEX = 50
AUTOSOMES = frozenset(f"chr{i}" for i in range(1, 23))
# Every female assembly has exactly 0 bases of high-confidence chrY and the lowest male has 889,749,
# so any threshold between the two separates them; this one is far from both.
MIN_MALE_CHR_Y_BASES = 500_000


def gsutil_ls(patterns, cache_path, refresh):
    """Return the paths matching one or more gs:// wildcard patterns, caching the listing on disk."""
    if cache_path.exists() and not refresh:
        return cache_path.read_text().split()

    paths = []
    for pattern in patterns:
        print(f"listing {pattern}")
        result = subprocess.run(["gsutil", "ls", pattern], capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"gsutil ls {pattern} failed (exit {result.returncode}): {result.stderr[-2000:]}")
        paths.extend(result.stdout.split())

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache_path.write_text("\n".join(paths) + "\n")
    return paths


def high_confidence_bed_paths(refresh):
    """Return {sample_id: gs:// path of its DipCall high-confidence BED}.

    A sample assembled in more than one batch (HPRC release 1 samples were also assembled for
    release 2) has a BED under each batch's directory. The first listed wins, which makes the choice
    deterministic and keeps it stable as later batches are added.
    """
    paths = gsutil_ls([f"{DIPCALL_DIR}/{subdir}*/*.dip.bed.gz" for subdir in DIPCALL_SUBDIRS],
                      CACHE_DIR / "dipcall_bed_paths.txt", refresh)
    bed_path_by_sample = {}
    for path in paths:
        bed_path_by_sample.setdefault(path.split("/")[-2], path)
    return bed_path_by_sample


def truth_vcf_paths(refresh):
    """Return {sample_id: gs:// path of its filter_vcf_to_tandem_repeats high-confidence VCF}."""
    return {path.split("/")[-2]: path for path in gsutil_ls(
        [f"{TRUTH_VCF_DIR}/*/*.high_confidence_regions.vcf.gz"],
        CACHE_DIR / "truth_vcf_paths.txt", refresh)}


def download_beds(bed_path_by_sample):
    """Copy any high-confidence BEDs that aren't cached locally yet. Returns {sample_id: local path}."""
    BED_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    local_path_by_sample = {sample_id: BED_CACHE_DIR / f"{sample_id}.dip.bed.gz"
                            for sample_id in bed_path_by_sample}
    missing = [bed_path_by_sample[sample_id] for sample_id, local_path in local_path_by_sample.items()
               if not local_path.exists()]
    if missing:
        print(f"downloading {len(missing)} high-confidence BED files to {BED_CACHE_DIR}")
        subprocess.run(["gsutil", "-m", "cp", *missing, str(BED_CACHE_DIR)], check=True)
    else:
        print(f"all {len(local_path_by_sample)} high-confidence BED files already cached in {BED_CACHE_DIR}")
    return local_path_by_sample


def high_confidence_bases_per_chrom(bed_path):
    """Return {chrom: total bases} covered by one sample's high-confidence regions."""
    totals = {}
    with gzip.open(bed_path, "rt") as f:
        for line in f:
            chrom, start, end = line.split("\t")[:3]
            totals[chrom] = totals.get(chrom, 0) + int(end) - int(start)
    return totals


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--refresh", action="store_true",
                        help="re-list the cloud directories instead of reusing the cached listings")
    args = parser.parse_args()

    crams = pd.read_table(CRAM_LIST_PATH).dropna(subset=["sample_id"])
    df = crams.merge(pd.read_table(SAMPLE_INFO_PATH), left_on="sample_id", right_on="Sample", how="left")
    df = df.dropna(subset=["Population", "Gender"])  # drops samples absent from the 1kGP metadata sheet
    df = df.rename(columns={"Gender": "Gender_1kGP_metadata"})
    print(f"{len(df):,} 1kGP samples have a short-read CRAM and 1kGP metadata")

    bed_path_by_sample = high_confidence_bed_paths(args.refresh)
    vcf_path_by_sample = truth_vcf_paths(args.refresh)
    df = df[df["sample_id"].isin(bed_path_by_sample) & df["sample_id"].isin(vcf_path_by_sample)]
    print(f"{len(df):,} of them also have DipCall and filter_vcf_to_tandem_repeats output")

    df = df.assign(
        high_confidence_bed_path=df["sample_id"].map(bed_path_by_sample),
        high_confidence_vcf_path=df["sample_id"].map(vcf_path_by_sample))

    local_bed_path_by_sample = download_beds({s: bed_path_by_sample[s] for s in df["sample_id"]})
    bases_per_chrom = {sample_id: high_confidence_bases_per_chrom(local_path)
                       for sample_id, local_path in local_bed_path_by_sample.items()}
    df = df.assign(
        autosomal_high_confidence_bp=df["sample_id"].map(
            lambda s: sum(bases for chrom, bases in bases_per_chrom[s].items() if chrom in AUTOSOMES)),
        chrY_high_confidence_bp=df["sample_id"].map(lambda s: bases_per_chrom[s].get("chrY", 0)))

    # Sex decides the 50/50 split and is passed to ExpansionHunter as --sex, so it is taken from the
    # assembly rather than from the 1kGP metadata sheet, which is wrong for at least one sample
    # (HG02300, listed male in 20130606_sample_info_1kGP.tsv but female in the official 1kGP metadata,
    # in both assembly tables, and by its own assembly). The two groups are cleanly separated: every
    # female assembly has exactly zero high-confidence chrY, and the lowest male has 889,749 bp.
    df = df.assign(Gender=df["chrY_high_confidence_bp"].map(
        lambda bp: "male" if bp >= MIN_MALE_CHR_Y_BASES else "female"))
    disagrees = df["Gender"] != df["Gender_1kGP_metadata"]
    if disagrees.any():
        print(f"\nNOTE: {disagrees.sum()} samples' sex was taken from the assembly rather than from "
              f"{SAMPLE_INFO_PATH.name}:")
        for row in df[disagrees].itertuples(index=False):
            print(f"  {row.sample_id}: metadata says {row.Gender_1kGP_metadata}, assembly has "
                  f"{row.chrY_high_confidence_bp:,} bp of high-confidence chrY, so it is treated as {row.Gender}")

    df = df.sort_values(["autosomal_high_confidence_bp", "sample_id"], ascending=[False, True])
    selected = df.groupby("Gender", group_keys=False).head(SAMPLES_PER_SEX)
    for sex in ("male", "female"):
        in_sex = selected[selected["Gender"] == sex]["autosomal_high_confidence_bp"]
        if len(in_sex) < SAMPLES_PER_SEX:
            print(f"\nWARNING: only {len(in_sex)} {sex} samples qualify, fewer than the {SAMPLES_PER_SEX} requested")
        print(f"\n{len(in_sex)} {sex} samples selected, autosomal high-confidence bases "
              f"{in_sex.min():,} to {in_sex.max():,} (cutoff excluded "
              f"{(df[df['Gender'] == sex]['autosomal_high_confidence_bp'] < in_sex.min()).sum()} samples)")

    selected = selected.sort_values(["Gender", "autosomal_high_confidence_bp"], ascending=[True, False])
    selected[["sample_id", "Population", "Population Description", "Gender", "cram_path", "crai_path",
              "high_confidence_vcf_path", "high_confidence_bed_path",
              "autosomal_high_confidence_bp"]].to_csv(OUTPUT_PATH, sep="\t", index=False)
    print(f"\nwrote {len(selected)} samples to {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
