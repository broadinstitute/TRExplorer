"""Aggregate boundary-extension profiles across loci within each tier, to find where each metric's
drop-off "knee" sits and how far a purity threshold rule would extend the boundary.

Usage: python3 find_thresholds.py <tier1_prefix> <tier2_prefix> <tier3_prefix> <output_dir>
    (prefixes are the --output-prefix values used with compute_boundary_profiles.py, i.e. this
    script expects <prefix>.profiles.pkl and <prefix>.purity_crossings.tsv to exist for each)

Writes to <output_dir>/:
  aggregate_curves.pkl          -- per-tier, per-metric median/p25/p75 curves vs offset (both sides pooled)
  threshold_summary.tsv         -- per tier x threshold, median/p25/p75 safe-extension-offset
  figures/*.png                 -- one plot per metric x tier-comparison, for embedding in the report
"""

import argparse
import pickle
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams.update(
    {
        "figure.facecolor": "#fcfcfb",
        "axes.facecolor": "#fcfcfb",
        "axes.edgecolor": "#c3c2b7",
        "axes.labelcolor": "#0b0b0b",
        "text.color": "#0b0b0b",
        "xtick.color": "#52514e",
        "ytick.color": "#52514e",
        "grid.color": "#e1e0d9",
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "sans-serif"],
        "axes.grid": True,
        "grid.linewidth": 0.6,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "lines.linewidth": 1.8,
    }
)

TIER_LABELS = {
    "known_disease": "Known disease-associated loci",
    "same_canonical_motif": "Loci sharing a disease motif",
    "trexplorer_v1_sample": "TRExplorer v1 (sampled)",
}
TIER_COLORS = {
    "known_disease": "#2a78d6",
    "same_canonical_motif": "#eb6834",
    "trexplorer_v1_sample": "#1baf7a",
}
PURITY_THRESHOLDS = [1.00, 0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.60, 0.50]


def load_tier(prefix):
    with open(f"{prefix}.profiles.pkl", "rb") as pickle_file:
        profiles = pickle.load(pickle_file)
    crossings = pd.read_csv(f"{prefix}.purity_crossings.tsv", sep="\t")
    return profiles, crossings


def stack_curves(profiles, metric, max_offset):
    """Pool left+right curves for `metric` across all loci in `profiles` onto a common 1..max_offset
    grid (loci with shorter curves contribute NaN past their own max_offset_used, and are excluded
    from that offset's percentile calculation via nanpercentile).
    """
    stacked = np.full((len(profiles) * 2, max_offset), np.nan)
    row = 0
    for locus in profiles.values():
        for side in ("left", "right"):
            curve = locus["profile"][side].get(metric)
            if curve is None or len(curve) == 0:
                row += 1
                continue
            n = min(len(curve), max_offset)
            stacked[row, :n] = curve[:n]
            row += 1
    return stacked


def summarize_curve(stacked):
    with np.errstate(all="ignore"):
        median = np.nanmedian(stacked, axis=0)
        p25 = np.nanpercentile(stacked, 25, axis=0)
        p75 = np.nanpercentile(stacked, 75, axis=0)
        n_loci = np.sum(~np.isnan(stacked), axis=0)
    return {"median": median, "p25": p25, "p75": p75, "n_loci": n_loci}


def plot_metric_comparison(aggregate_by_tier, metric, ylabel, output_path, hline=None, ylim=None):
    fig, ax = plt.subplots(figsize=(8, 5))
    for tier, curves in aggregate_by_tier.items():
        if metric not in curves:
            continue
        summary = curves[metric]
        offsets = np.arange(1, len(summary["median"]) + 1)
        ax.plot(offsets, summary["median"], label=TIER_LABELS.get(tier, tier), color=TIER_COLORS.get(tier))
        ax.fill_between(offsets, summary["p25"], summary["p75"], color=TIER_COLORS.get(tier), alpha=0.15)
    if hline is not None:
        ax.axhline(hline, color="black", linestyle="--", linewidth=1, label=f"reference = {hline}")
    ax.set_xlabel("Distance from current annotated boundary (bp)")
    ax.set_ylabel(ylabel)
    if ylim:
        ax.set_ylim(*ylim)
    ax.legend(fontsize=9)
    ax.set_title(f"{ylabel} vs. extension distance (median, IQR band; left+right pooled)")
    fig.tight_layout()
    fig.savefig(output_path, dpi=130)
    plt.close(fig)


def threshold_summary_table(all_crossings):
    """median/p25/p75_offset_among_crossed: percentiles of the safe-extension offset among locus-sides
    whose purity actually crossed below the threshold within the scanned window -- excludes never-
    crossed (right-censored) sides entirely, so these are NOT tier-wide safe-extension percentiles.
    median/p25/p75_offset_censored_lower_bound: the same percentiles but with never-crossed sides
    included at their max_offset_used (the offset actually scanned to) -- a conservative LOWER BOUND
    on the true tier-wide percentiles, since a never-crossed side's real safe extension could be
    longer than what was scanned.
    """
    rows = []
    for tier, group in all_crossings.groupby("tier"):
        for threshold in PURITY_THRESHOLDS:
            sub = group[group["threshold"] == threshold]
            for column, label in [
                ("full_region_first_offset_below", "full_region"),
                ("local_flank_first_offset_below", "local_flank"),
            ]:
                # column is the first 1-based offset where purity DROPS BELOW threshold (i.e. that
                # offset already fails) -- the actual safe extension is one base less, floored at 0.
                crossed_mask = sub[column].notna()
                if crossed_mask.sum() == 0:
                    continue
                values_among_crossed = (sub.loc[crossed_mask, column] - 1).clip(lower=0)
                censored_lower_bound = sub[column].sub(1).clip(lower=0)
                censored_lower_bound[~crossed_mask] = sub.loc[~crossed_mask, "max_offset_used"]
                rows.append(
                    {
                        "tier": tier,
                        "threshold": threshold,
                        "curve": label,
                        "n_loci_sides": len(values_among_crossed),
                        "median_offset_among_crossed": values_among_crossed.median(),
                        "p25_offset_among_crossed": values_among_crossed.quantile(0.25),
                        "p75_offset_among_crossed": values_among_crossed.quantile(0.75),
                        "median_offset_censored_lower_bound": censored_lower_bound.median(),
                        "p25_offset_censored_lower_bound": censored_lower_bound.quantile(0.25),
                        "p75_offset_censored_lower_bound": censored_lower_bound.quantile(0.75),
                        "fraction_never_crossed": 1 - crossed_mask.mean(),
                    }
                )
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("tier1_prefix")
    parser.add_argument("tier2_prefix")
    parser.add_argument("tier3_prefix")
    parser.add_argument("output_dir")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    figures_dir = output_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    tier_prefixes = {
        "known_disease": args.tier1_prefix,
        "same_canonical_motif": args.tier2_prefix,
        "trexplorer_v1_sample": args.tier3_prefix,
    }

    all_crossings = []
    aggregate_by_tier = {}
    for tier, prefix in tier_prefixes.items():
        print(f"loading {tier} ({prefix})...")
        profiles, crossings = load_tier(prefix)
        all_crossings.append(crossings)

        max_offset = max(locus["max_offset_used"] for locus in profiles.values())
        metrics_present = next(iter(profiles.values()))["profile"]["left"].keys()

        curves = {}
        for metric in ["full_region_purity", "purity", "linguistic_complexity", "kmer_jaccard_similarity"]:
            if metric not in metrics_present:
                continue
            print(f"  stacking {metric}...")
            curves[metric] = summarize_curve(stack_curves(profiles, metric, max_offset))
        if "normalized_edit_distance" in metrics_present:
            curves["normalized_edit_distance"] = summarize_curve(
                1 - stack_curves(profiles, "normalized_edit_distance", max_offset)
            )
        aggregate_by_tier[tier] = curves

    all_crossings_df = pd.concat(all_crossings, ignore_index=True)

    with open(output_dir / "aggregate_curves.pkl", "wb") as pickle_file:
        pickle.dump(aggregate_by_tier, pickle_file)

    summary_df = threshold_summary_table(all_crossings_df)
    summary_df.to_csv(output_dir / "threshold_summary.tsv", sep="\t", index=False)
    print(summary_df.to_string(index=False))

    plot_metric_comparison(
        aggregate_by_tier,
        "full_region_purity",
        "Full-region reference-repeat purity\n(EH ReferenceRepeatPurity convention: core + extension)",
        figures_dir / "full_region_purity.png",
        hline=0.90,
        ylim=(0, 1.02),
    )
    plot_metric_comparison(
        aggregate_by_tier,
        "purity",
        "Local flank purity\n(added flank segment only vs. perfect motif repeat)",
        figures_dir / "local_flank_purity.png",
        hline=0.90,
        ylim=(0, 1.02),
    )
    plot_metric_comparison(
        aggregate_by_tier,
        "linguistic_complexity",
        "Linguistic complexity of flank\n(0 = perfectly repetitive, 1 = maximally complex; 24bp sliding window)",
        figures_dir / "linguistic_complexity.png",
        ylim=(0, 1.02),
    )
    plot_metric_comparison(
        aggregate_by_tier,
        "kmer_jaccard_similarity",
        "K-mer Jaccard similarity of flank to motif\n(24bp sliding window)",
        figures_dir / "kmer_jaccard_similarity.png",
        ylim=(0, 1.02),
    )
    if "normalized_edit_distance" in aggregate_by_tier.get("known_disease", {}):
        plot_metric_comparison(
            aggregate_by_tier,
            "normalized_edit_distance",
            "1 - normalized DP edit distance of flank to perfect motif repeat\n(1 = identical, 0 = maximally different)",
            figures_dir / "edit_distance_identity.png",
            ylim=(0, 1.02),
        )

    print(f"\nwrote aggregate_curves.pkl, threshold_summary.tsv, and figures/ to {output_dir}")


if __name__ == "__main__":
    main()
