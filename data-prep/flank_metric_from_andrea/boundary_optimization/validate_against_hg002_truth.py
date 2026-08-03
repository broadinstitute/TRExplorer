"""Test whether reference-boundary purity (at the CURRENT annotated boundary, computed independently
from the genome, not from any sample) predicts ExpansionHunter genotype accuracy against real HG002
truth, using the existing genome-wide comparison table from str-truth-set-v2 (578k+ loci) -- no new
runs needed, this is real long-read-truth-set-derived data at full-genome scale.

This answers a correlational version of the boundary-optimization question directly: loci whose
CURRENT reference definition already has low purity are proxies for "badly-bounded" loci (either the
true repeat wasn't a clean match to begin with, or the boundary should have been extended further to
capture more of a genuinely repetitive region but wasn't) -- if EH accuracy tracks this purity, and
if errors at low purity skew toward overestimation (the EP400 symptom), that's strong independent
support for a purity-based boundary rule. It can't test "does extending the boundary fix it" though
-- that's what the small Hail Batch run is for.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

TRUTH_TABLE = "/Users/weisburd/code/str-truth-set-v2/results/HG002.tandem_repeat_genotypes.for_comparison.with_EHv5-bw2-optimized_vs_Truth_columns.tsv.gz"
FASTA_PATH = "/Users/weisburd/hg38.fa"
OUTPUT_DIR = Path(__file__).parent / "data"

PURITY_BINS = [0.0, 0.65, 0.75, 0.85, 0.90, 0.95, 1.0001]
PURITY_BIN_LABELS = ["0.00-0.65", "0.65-0.75", "0.75-0.85", "0.85-0.90", "0.90-0.95", "0.95-1.00"]


def fast_core_purity(core_seq, motif):
    """Vectorized best-phase purity: max over motif rotations of (matched bases / length), using
    numpy broadcasting across phases instead of a python double loop -- needed here because this
    runs across ~578k loci, some with long (VNTR-scale) reference regions.
    """
    core_seq = core_seq.upper()
    motif = motif.upper()
    motif_length = len(motif)
    core_length = len(core_seq)
    if motif_length == 0 or core_length == 0:
        return 0.0, motif

    core_array = np.frombuffer(core_seq.encode(), dtype=np.uint8)
    motif_array = np.frombuffer(motif.encode(), dtype=np.uint8)
    positions = np.arange(core_length)

    best_purity, best_phase = -1.0, 0
    for phase in range(motif_length):
        phase_index = (positions + phase) % motif_length
        matches = np.count_nonzero(core_array == motif_array[phase_index])
        purity = matches / core_length
        if purity > best_purity:
            best_purity, best_phase = purity, phase

    rotated_motif = motif[best_phase:] + motif[:best_phase]
    return best_purity, rotated_motif


def main():
    print(f"loading {TRUTH_TABLE} ...")
    df = pd.read_csv(
        TRUTH_TABLE,
        sep="\t",
        usecols=[
            "LocusId",
            "Motif",
            "MotifSize",
            "NumRepeatsInReference",
            "ReferenceLocusSize (bp)",
            "TruthSetOrNegativeLocus",
            "IsPureRepeat",
            "Coverage: EHv5-bw2-optimized",
            "Variant: Concordance: EHv5-bw2-optimized vs Truth",
            "DiffRepeats: Allele 1: EHv5-bw2-optimized - Truth",
            "DiffRepeats: Allele 2: EHv5-bw2-optimized - Truth",
        ],
    )
    print(f"loaded {len(df)} rows")

    df = df[df["TruthSetOrNegativeLocus"] == "TruthSet"].copy()
    df = df.dropna(subset=["Variant: Concordance: EHv5-bw2-optimized vs Truth"])
    print(f"{len(df)} rows after restricting to TruthSet loci with a concordance call")

    parts = df["LocusId"].str.split("-", n=3, expand=True)
    df["chrom"] = "chr" + parts[0]
    df["start_0based"] = parts[1].astype(int)
    df["end"] = parts[2].astype(int)

    fasta = pysam.FastaFile(FASTA_PATH)
    purities = np.empty(len(df))
    rotated_motifs = [None] * len(df)
    print("computing reference-boundary purity for each locus (this may take a few minutes)...")
    for i, row in enumerate(df.itertuples()):
        if i % 50000 == 0:
            print(f"  {i}/{len(df)}")
        core_seq = fasta.fetch(row.chrom, row.start_0based, row.end)
        purity, rotated_motif = fast_core_purity(core_seq, row.Motif)
        purities[i] = purity
        rotated_motifs[i] = rotated_motif

    df["reference_boundary_purity"] = purities
    df["purity_bin"] = pd.cut(df["reference_boundary_purity"], PURITY_BINS, labels=PURITY_BIN_LABELS, right=False)

    df["is_exact_match"] = df["Variant: Concordance: EHv5-bw2-optimized vs Truth"] == "ExactlyTheSame"
    df["is_discordant"] = df["Variant: Concordance: EHv5-bw2-optimized vs Truth"] == "Discordant"
    df["mean_diff_repeats"] = df[
        ["DiffRepeats: Allele 1: EHv5-bw2-optimized - Truth", "DiffRepeats: Allele 2: EHv5-bw2-optimized - Truth"]
    ].mean(axis=1)

    summary = (
        df.groupby("purity_bin", observed=True)
        .agg(
            n_loci=("LocusId", "count"),
            exact_match_rate=("is_exact_match", "mean"),
            discordant_rate=("is_discordant", "mean"),
            median_diff_repeats=("mean_diff_repeats", "median"),
            mean_diff_repeats=("mean_diff_repeats", "mean"),
            frac_overestimated=("mean_diff_repeats", lambda x: (x.dropna() > 0).mean()),
            frac_underestimated=("mean_diff_repeats", lambda x: (x.dropna() < 0).mean()),
        )
        .reset_index()
    )
    print("\n=== EH accuracy vs. current reference-boundary purity (all loci) ===")
    print(summary.to_string(index=False))

    pure_only = df[df["IsPureRepeat"]]
    summary_pure = (
        pure_only.groupby("purity_bin", observed=True)
        .agg(
            n_loci=("LocusId", "count"),
            exact_match_rate=("is_exact_match", "mean"),
            discordant_rate=("is_discordant", "mean"),
            median_diff_repeats=("mean_diff_repeats", "median"),
            frac_overestimated=("mean_diff_repeats", lambda x: (x > 0).mean()),
        )
        .reset_index()
    )
    print("\n=== Same, restricted to IsPureRepeat==True loci (removes interruption confound) ===")
    print(summary_pure.to_string(index=False))

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    df_out_columns = [
        "LocusId",
        "chrom",
        "start_0based",
        "end",
        "Motif",
        "MotifSize",
        "IsPureRepeat",
        "reference_boundary_purity",
        "purity_bin",
        "Variant: Concordance: EHv5-bw2-optimized vs Truth",
        "mean_diff_repeats",
        "Coverage: EHv5-bw2-optimized",
    ]
    df[df_out_columns].to_csv(OUTPUT_DIR / "hg002_purity_vs_accuracy.tsv", sep="\t", index=False)
    summary.to_csv(OUTPUT_DIR / "hg002_purity_vs_accuracy.summary_all.tsv", sep="\t", index=False)
    summary_pure.to_csv(OUTPUT_DIR / "hg002_purity_vs_accuracy.summary_pure_only.tsv", sep="\t", index=False)
    print(f"\nwrote outputs to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
