"""Analyze the round-3 Hail Batch validation: compares gap-mismatch thresholds 0.15 vs 0.20 for the
motif-reanchoring rule (see motif_reanchor_extension.py, report_v2_motif_reanchoring.html) across the
98 TRExplorer v1 loci that get nonzero extension under either threshold (threshold_sweep_trexplorer_v1.tsv),
plus 4 ARX_1/ARX_2 entries re-run with the corrected --sex male flag (the original round-2 validation
used --sex female for HG002, which is actually male, giving wrong ploidy for all chrX/chrY loci).

Usage:
    python3 analyze_round3_threshold_sweep.py
"""

import json
from pathlib import Path

import pandas as pd

DATA_DIR = Path(__file__).parent / "data"
RESULTS_PATH = DATA_DIR / "HG002.reanchor_validation_round3_male.json"
SWEEP_TSV = DATA_DIR / "threshold_sweep_trexplorer_v1.tsv"
CATALOG_PATH = DATA_DIR / "reanchor_validation_catalog_3_male.json"


def mean_pok(locus_result):
    variant = next(iter(locus_result["Variants"].values()))
    if "AlleleQualityMetrics" not in variant:
        return None
    alleles = variant["AlleleQualityMetrics"]["Alleles"]
    return sum(a["pOk"] for a in alleles) / len(alleles)


def max_ci_size(locus_result):
    variant = next(iter(locus_result["Variants"].values()))
    if "AlleleQualityMetrics" not in variant:
        return None
    alleles = variant["AlleleQualityMetrics"]["Alleles"]
    return max(a["ConfidenceIntervalDividedByAlleleSize"] for a in alleles)


def genotype_str(locus_result):
    variant = next(iter(locus_result["Variants"].values()))
    return variant.get("Genotype", "NO_CALL")


def main():
    results = json.load(open(RESULTS_PATH))["LocusResults"]
    catalog = json.load(open(CATALOG_PATH))
    catalog_ids = {e["LocusId"] for e in catalog}
    print(f"catalog: {len(catalog)} entries, results: {len(results)} loci ({len(catalog_ids) - len(results)} skipped, e.g. extended region > 2x read length)")

    n_no_call = 0
    rows = []
    for locus_id, result in results.items():
        if "__" in locus_id:
            base_id, labels_str = locus_id.split("__", 1)
            labels = labels_str.split("+")
        else:
            # ARX_*_MALE_FIX entries: "{base}_{baseline|reanchored}_MALE_FIX"
            label = "baseline" if "_baseline_" in locus_id else "reanchored"
            base_id = locus_id.replace(f"_{label}_MALE_FIX", "_MALE_FIX")
            labels = [label]
        pok = mean_pok(result)
        if pok is None:
            n_no_call += 1
            continue
        ci = max_ci_size(result)
        gt = genotype_str(result)
        for label in labels:
            rows.append({"base_id": base_id, "label": label, "pOk": pok, "ci_size": ci, "genotype": gt})

    print(f"{n_no_call} locus entries excluded (no genotype call)")

    long_df = pd.DataFrame(rows)
    wide = long_df.pivot_table(index="base_id", columns="label", values="pOk", aggfunc="first").reset_index()
    wide_ci = long_df.pivot_table(index="base_id", columns="label", values="ci_size", aggfunc="first").reset_index()
    wide_gt = long_df.pivot_table(index="base_id", columns="label", values="genotype", aggfunc="first").reset_index()

    print("\n=== ARX_1/ARX_2 with corrected --sex male ===")
    arx = wide[wide["base_id"].str.contains("MALE_FIX")]
    print(arx.to_string(index=False))
    arx_gt = wide_gt[wide_gt["base_id"].str.contains("MALE_FIX")]
    print(arx_gt.to_string(index=False))

    sweep = pd.read_csv(SWEEP_TSV, sep="\t")
    disagree = sweep[~sweep["agree"]].copy()
    print(f"\n{len(disagree)} disagreement loci (0.15 != 0.20) out of {len(sweep)} candidates")

    merged = disagree.merge(wide, left_on="locus_id", right_on="base_id", how="left")
    n_missing = merged["baseline"].isna().sum()
    print(f"{n_missing} disagreement loci missing EH results (skipped, e.g. region > 2x read length)")
    merged = merged.dropna(subset=["baseline"])

    merged["delta_015"] = merged["reanchored_015"] - merged["baseline"]
    merged["delta_020"] = merged["reanchored_020"] - merged["baseline"]
    merged["020_worse_than_015"] = merged["delta_020"] < merged["delta_015"]

    out_path = DATA_DIR / "round3_threshold_comparison.tsv"
    merged.to_csv(out_path, sep="\t", index=False)
    print(f"\nWrote {len(merged)} rows to {out_path}")

    print("\n=== Summary: 0.15 vs 0.20 on the disagreement set ===")
    print(f"mean delta_pOk @0.15: {merged['delta_015'].mean():+.4f}")
    print(f"mean delta_pOk @0.20: {merged['delta_020'].mean():+.4f}")
    print(f"median delta_pOk @0.15: {merged['delta_015'].median():+.4f}")
    print(f"median delta_pOk @0.20: {merged['delta_020'].median():+.4f}")
    print(f"n loci worse at 0.15 (delta<-0.05): {(merged['delta_015'] < -0.05).sum()}")
    print(f"n loci worse at 0.20 (delta<-0.05): {(merged['delta_020'] < -0.05).sum()}")
    print(f"n loci better at 0.15 (delta>0.05): {(merged['delta_015'] > 0.05).sum()}")
    print(f"n loci better at 0.20 (delta>0.05): {(merged['delta_020'] > 0.05).sum()}")
    print(f"n loci where 0.20 is strictly worse than 0.15: {merged['020_worse_than_015'].sum()} / {len(merged)}")

    # "actually extended" must come from the extension columns, not from whether pOk happened to
    # change -- an extension can land on an identical pOk, which a delta!=0 check would wrongly
    # exclude.
    extended_015 = (merged["left_ext_015"] > 0) | (merged["right_ext_015"] > 0)
    extended_020 = (merged["left_ext_020"] > 0) | (merged["right_ext_020"] > 0)
    print(f"n loci actually extended @0.15: {extended_015.sum()} / {len(merged)}")
    print(f"n loci actually extended @0.20: {extended_020.sum()} / {len(merged)}")
    print(f"sum(delta_015): {merged['delta_015'].sum():+.3f}   sum(delta_020): {merged['delta_020'].sum():+.3f}")

    print("\n=== Worst regressions at 0.20 (sorted by delta_020) ===")
    print(merged.sort_values("delta_020")[["locus_id", "motif", "baseline", "reanchored_015", "reanchored_020", "delta_015", "delta_020"]].head(15).to_string(index=False))

    print("\n=== Biggest improvements at 0.20 (sorted by -delta_020) ===")
    print(merged.sort_values("delta_020", ascending=False)[["locus_id", "motif", "baseline", "reanchored_015", "reanchored_020", "delta_015", "delta_020"]].head(10).to_string(index=False))

    print("\n=== Any nonzero movement at 0.15 ===")
    moved_015 = merged[merged["delta_015"] != 0]
    print(moved_015[["locus_id", "motif", "baseline", "reanchored_015", "delta_015"]].to_string(index=False))


if __name__ == "__main__":
    main()
