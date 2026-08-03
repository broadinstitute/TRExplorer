"""Analyze round 4: narrow-vs-wide pOk comparison for the 534 loci (of the top 20,000 ranked by
flank repetitiveness) that get extended by the motif-reanchoring rule at threshold 0.15 or 0.20.

Usage:
    python3 analyze_round4_flank_repetitive.py
"""

import json
from pathlib import Path

import pandas as pd

DATA_DIR = Path(__file__).parent / "data"
RESULTS_PATH = DATA_DIR / "HG002.flank_repetitive_validation.json"
EXTENSIONS_TSV = DATA_DIR / "top20k_flank_repetitive_extensions.tsv"
CATALOG_PATH = DATA_DIR / "flank_repetitive_validation_catalog.json"


def mean_pok(locus_result):
    variant = next(iter(locus_result["Variants"].values()))
    if "AlleleQualityMetrics" not in variant:
        return None
    alleles = variant["AlleleQualityMetrics"]["Alleles"]
    return sum(a["pOk"] for a in alleles) / len(alleles)


def genotype_str(locus_result):
    variant = next(iter(locus_result["Variants"].values()))
    return variant.get("Genotype", "NO_CALL")


def main():
    results = json.load(open(RESULTS_PATH))["LocusResults"]
    catalog = json.load(open(CATALOG_PATH))
    print(f"catalog: {len(catalog)} entries, results: {len(results)} loci ({len(catalog) - len(results)} skipped, e.g. extended region > 2x read length)")

    n_no_call = 0
    rows = []
    for locus_id, result in results.items():
        base_id, labels_str = locus_id.split("__", 1)
        labels = labels_str.split("+")
        pok = mean_pok(result)
        if pok is None:
            n_no_call += 1
            continue
        gt = genotype_str(result)
        for label in labels:
            rows.append({"base_id": base_id, "label": label, "pOk": pok, "genotype": gt})
    print(f"{n_no_call} locus entries excluded (no genotype call)")

    long_df = pd.DataFrame(rows)
    wide = long_df.pivot_table(index="base_id", columns="label", values="pOk", aggfunc="first").reset_index()

    ext = pd.read_csv(EXTENSIONS_TSV, sep="\t")
    merged = ext.merge(wide, left_on="LocusId", right_on="base_id", how="inner")
    print(f"\n{len(merged)} / {len(ext)} extended loci have scored baseline pOk")

    merged["delta_015"] = merged["reanchored_015"] - merged["baseline"]
    merged["delta_020"] = merged["reanchored_020"] - merged["baseline"]

    out_path = DATA_DIR / "round4_flank_repetitive_comparison.tsv"
    merged.to_csv(out_path, sep="\t", index=False)
    print(f"Wrote {len(merged)} rows to {out_path}")

    # has_015/has_020 must check the ACTUAL extension amount, not just whether a reanchored_XXX pOk
    # value exists -- build_flank_repetitive_validation_catalog.py combines labels (e.g.
    # "baseline+reanchored_015") when a locus's region is unchanged at that threshold, so an
    # unextended locus still gets a (trivially-equal) reanchored_015 value and would otherwise be
    # miscounted as genuinely extended.
    has_015 = (merged["left_ext_015"] > 0) | (merged["right_ext_015"] > 0)
    has_020 = (merged["left_ext_020"] > 0) | (merged["right_ext_020"] > 0)
    print(f"\n{has_015.sum()} loci genuinely extended at 0.15")
    print(f"{has_020.sum()} loci genuinely extended at 0.20")

    for label, mask, col in [("0.15", has_015, "delta_015"), ("0.20", has_020, "delta_020")]:
        sub = merged.loc[mask, col].dropna()
        n_pos = (sub > 0.05).sum()
        n_neg = (sub < -0.05).sum()
        n_flat = len(sub) - n_pos - n_neg
        print(
            f"\n=== threshold {label}: n={len(sub)} ===\n"
            f"  mean delta_pOk: {sub.mean():+.4f}   median: {sub.median():+.4f}\n"
            f"  better (>+0.05): {n_pos} ({100*n_pos/len(sub):.1f}%)   "
            f"worse (<-0.05): {n_neg} ({100*n_neg/len(sub):.1f}%)   flat: {n_flat} ({100*n_flat/len(sub):.1f}%)\n"
            f"  sum(delta): {sub.sum():+.2f}   sum(positive): {sub[sub>0].sum():+.2f}   sum(negative): {sub[sub<0].sum():+.2f}"
        )

    print("\n=== Correlation: does max_flank_motif_similarity predict delta_pOk @0.20? ===")
    sub020 = merged[has_020].dropna(subset=["delta_020"])
    corr = sub020[["max_flank_motif_similarity", "delta_020"]].corr(method="spearman").iloc[0, 1]
    print(f"Spearman correlation(max_flank_motif_similarity, delta_pOk@0.20) = {corr:.3f}  (n={len(sub020)})")

    print("\n=== Worst regressions @0.20 ===")
    print(sub020.sort_values("delta_020")[["LocusId", "motif", "max_flank_motif_similarity", "baseline", "reanchored_020", "delta_020"]].head(10).to_string(index=False))

    print("\n=== Biggest improvements @0.20 ===")
    print(sub020.sort_values("delta_020", ascending=False)[["LocusId", "motif", "max_flank_motif_similarity", "baseline", "reanchored_020", "delta_020"]].head(10).to_string(index=False))

    print("\n=== All loci extended @0.15 (n small enough to show all) ===")
    sub015 = merged[has_015].dropna(subset=["delta_015"])
    print(sub015.sort_values("delta_015")[["LocusId", "motif", "max_flank_motif_similarity", "baseline", "reanchored_015", "delta_015"]].to_string(index=False))


if __name__ == "__main__":
    main()
