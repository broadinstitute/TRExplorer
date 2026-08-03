"""Relate the EXISTING reference-repeat-purity value of each locus (as currently annotated, no
boundary changes) to REAL empirical genotype quality, by joining Andrea's flank-metrics catalog
(415k loci, 2-6bp motifs, has ReferenceRepeatPurity) against ExpansionHunter-bw2's own truth-set
comparison tables for HG002 and HG00438 (real short reads vs. assembly-derived truth genotypes).

This answers a different question than compute_profiles.py: not "how far COULD we extend the
boundary before purity drops" but "does the CURRENT purity value actually predict genotyping
quality on real data" -- i.e. does the signal we're proposing to use for boundary decisions have
real predictive power. Purity bins match benchmark_on_truth_set_v2/make_purity_plots.py for
consistency with prior analysis in this codebase.

Usage:
    python3 join_with_truth_quality.py --output data/purity_vs_quality.tsv
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

CATALOG_PATH = "/Users/weisburd/code/tandem-repeat-explorer/data-prep/flank_metric_from_andrea/TR_catalog_2_6bp_with_flank_metrics.tsv.gz"

TRUTH_TABLES = {
    "HG002": "/Users/weisburd/code/str-truth-set-v2/results/HG002.tandem_repeat_genotypes.for_comparison.with_EHv5-bw2-optimized_vs_Truth_columns.tsv.gz",
    "HG00438": "/Users/weisburd/code/str-truth-set-v2/results/HG00438.tandem_repeat_genotypes.for_comparison.with_EHv5-bw2-optimized_vs_Truth_columns.tsv.gz",
}

# Same bin edges as benchmark_on_truth_set_v2/make_purity_plots.py, for consistency with prior art.
PURITY_BINS = [
    ("<0.70", -1e9, 0.70),
    ("0.70-0.75", 0.70, 0.75),
    ("0.75-0.80", 0.75, 0.80),
    ("0.80-0.85", 0.80, 0.85),
    ("0.85-0.90", 0.85, 0.90),
    ("0.90-0.95", 0.90, 0.95),
    ("0.95-1.00", 0.95, 1e9),
]

TRUTH_COLUMNS = [
    "LocusId",
    "TruthSetOrNegativeLocus",
    "NumReadsTotal: EHv5-bw2-optimized",
    "Genotype: EHv5-bw2-optimized",
    "Variant: Concordance: EHv5-bw2-optimized vs Truth",
    "Q: Allele 1: EHv5-bw2-optimized",
    "Q: Allele 2: EHv5-bw2-optimized",
    "CI size: Allele 1: EHv5-bw2-optimized",
    "CI size: Allele 2: EHv5-bw2-optimized",
    "DiffRepeats: Allele 1: EHv5-bw2-optimized - Truth",
    "DiffRepeats: Allele 2: EHv5-bw2-optimized - Truth",
]

CATALOG_COLUMNS = [
    "LocusId",
    "MotifSize",
    "ReferenceMotif",
    "NumRepeatsInReference",
    "ReferenceRepeatPurity",
    "flank_complexity_left",
    "flank_complexity_right",
    "flank_motif_similarity_left",
    "flank_motif_similarity_right",
]


def purity_bin(value):
    for label, lo, hi in PURITY_BINS:
        if lo <= value < hi:
            return label
    return None


def load_catalog():
    df = pd.read_table(CATALOG_PATH, usecols=CATALOG_COLUMNS)
    return df


def load_and_join_one_sample(sample_name, truth_path, catalog):
    print(f"Loading {sample_name} truth-comparison table...")
    truth = pd.read_table(truth_path, usecols=TRUTH_COLUMNS)
    n_truth = len(truth)

    truth = truth[truth["TruthSetOrNegativeLocus"] == "TruthSet"]
    truth = truth[truth["Genotype: EHv5-bw2-optimized"].notna()]  # exclude no-calls
    print(f"  {len(truth):,} / {n_truth:,} rows are TruthSet loci with an EH call")

    n_before_dedup = len(truth)
    truth = truth.drop_duplicates(subset=["LocusId"])
    if len(truth) < n_before_dedup:
        print(f"  dropped {n_before_dedup - len(truth):,} duplicate-LocusId truth rows")

    merged = truth.merge(catalog, on="LocusId", how="inner", validate="one_to_one")
    print(f"  {len(merged):,} rows after joining to the flank-metrics catalog (2-6bp motifs only)")

    merged["sample"] = sample_name
    diff1 = merged["DiffRepeats: Allele 1: EHv5-bw2-optimized - Truth"].abs()
    diff2 = merged["DiffRepeats: Allele 2: EHv5-bw2-optimized - Truth"].abs()
    merged["max_abs_diff_repeats"] = np.fmax(diff1, diff2)
    merged["is_concordant"] = merged["Variant: Concordance: EHv5-bw2-optimized vs Truth"] == "ExactlyTheSame"
    merged["max_ci_size"] = np.fmax(
        merged["CI size: Allele 1: EHv5-bw2-optimized"], merged["CI size: Allele 2: EHv5-bw2-optimized"]
    )
    merged["min_q"] = np.fmin(merged["Q: Allele 1: EHv5-bw2-optimized"], merged["Q: Allele 2: EHv5-bw2-optimized"])
    return merged


def summarize_by_purity_bin(df):
    df = df.copy()
    df["purity_bin"] = df["ReferenceRepeatPurity"].apply(purity_bin)
    summary = (
        df.groupby("purity_bin", observed=True)
        .agg(
            n_loci=("LocusId", "count"),
            pct_concordant=("is_concordant", lambda s: 100 * s.mean()),
            median_abs_diff_repeats=("max_abs_diff_repeats", "median"),
            mean_abs_diff_repeats=("max_abs_diff_repeats", "mean"),
            median_ci_size=("max_ci_size", "median"),
            median_min_q=("min_q", "median"),
        )
        .reindex([label for label, _, _ in PURITY_BINS])
    )
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, help="output .tsv path for the per-locus joined table")
    parser.add_argument("--summary-output", help="output .tsv path for the purity-bin summary (default: derived from --output)")
    args = parser.parse_args()

    catalog = load_catalog()
    print(f"Loaded catalog: {len(catalog):,} loci (2-6bp motifs)")

    joined_frames = [load_and_join_one_sample(name, path, catalog) for name, path in TRUTH_TABLES.items()]
    joined = pd.concat(joined_frames, ignore_index=True)

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    joined.to_csv(output_path, sep="\t", index=False)
    print(f"\nWrote {len(joined):,} joined rows to {output_path}")

    summary = summarize_by_purity_bin(joined)
    summary_path = Path(args.summary_output) if args.summary_output else output_path.with_name(output_path.stem + "_summary_by_purity_bin.tsv")
    summary.to_csv(summary_path, sep="\t")
    print(f"\n=== Genotype quality by ReferenceRepeatPurity bin (pooled across {list(TRUTH_TABLES)}) ===")
    print(summary.to_string())
    print(f"\nWrote summary to {summary_path}")

    correlation = joined[["ReferenceRepeatPurity", "max_abs_diff_repeats"]].corr(method="spearman").iloc[0, 1]
    print(f"\nSpearman correlation(ReferenceRepeatPurity, |diff repeats from truth|) = {correlation:.3f}")


if __name__ == "__main__":
    main()
