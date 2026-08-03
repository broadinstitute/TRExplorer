"""Compare ExpansionHunter's per-allele pOk under the gap-purity rule's original vs. extended locus
definitions, using the corrected combined-catalog run (samples listed in
../expansion_hunter_genotype_quality/selected_1kGP_samples.tsv, optimized-streaming -- built from
the corrected combined_eh_catalog.json, which covers ALL 347,704 Andrea v2-exact-match loci as
__original entries, plus __extended entries for the 16,264 that the gap-purity rule extends, plus
the non-GCN/NGC known disease loci).

LocusId convention (post code-review fix): "{stable source LocusId}__{original|extended}" -- the
base id is now trivially recoverable via a single rsplit, unlike the earlier catalog where a locus's
own extended-region coordinates were embedded in the id and collided across different loci.

Produces build_pOk_delta_report.html with three sections:
  1. How many of Andrea's v2-exact-match loci were extended, by motif size and by quintile bin of
     max(flank_motif_similarity_left, flank_motif_similarity_right).
  2. Absolute pOk, 3 groups per bin (original pre-extension / extended post-extension / not
     extended), same x-axis binnings.
  3. Delta pOk (extended - original) for the 16,264 genuinely-extended loci, violin plots by motif
     size, and by quintile bin of max(flank_motif_similarity_left, flank_motif_similarity_right) --
     plus EP400 shown as its own violin on the motif-size plot (EP400 isn't in Andrea's catalog, so
     it has no flank_motif_similarity data and can't go on the other plot in that section).

Usage:
    python3 build_pOk_delta_report.py
    python3 build_pOk_delta_report.py --exclude-near-reference   # exclude alleles within +/-2
                                                                   # repeats of the reference allele
                                                                   # size (per locus definition). For
                                                                   # a locus with no extended
                                                                   # definition this filters per
                                                                   # allele; for a genuinely-extended
                                                                   # locus it's all-or-nothing across
                                                                   # BOTH definitions' alleles, so the
                                                                   # original/extended arms stay
                                                                   # matched populations
"""

import argparse
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DATA_DIR = Path(__file__).parent / "data"
JSON_DIR = DATA_DIR / "combined_catalog_eh_output" / "json"
ANDREA_CATALOG_PATH = Path(__file__).parent.parent / "TR_catalog_2_6bp_with_flank_metrics.tsv.gz"
ANDREA_SCAN_TSV = DATA_DIR / "gap_purity_scan_v2_exact_matches.tsv"
SAMPLE_TABLE_PATH = Path(__file__).parent.parent.parent / "expansion_hunter_genotype_quality" / "selected_1kGP_samples.tsv"

SAMPLES = pd.read_table(SAMPLE_TABLE_PATH)["sample_id"].tolist()

BLUE_LIGHT, BLUE_DARK = "#2a78d6", "#3987e5"
ORANGE_LIGHT, ORANGE_DARK = "#eb6834", "#d95926"

# Andrea-derived LocusIds are always "{chrom}-{start}-{end}-{motif}"; known-disease-locus LocusIds are
# gene names (e.g. "EP400", "RFC1") and don't fit this shape -- used to keep the two populations apart.
ANDREA_ID_PATTERN = r"^(chr)?[0-9XYM]+-\d+-\d+-[ACGTN]+$"

NEAR_REFERENCE_REPEAT_TOLERANCE = 2


def load_sample_pok(sample_id, exclude_near_reference=False):
    """Return list of (base_locus_id, definition, allele_number, pOk, allele_size, ref_repeat_count,
    is_near_reference) for one sample's JSON. allele_size and ref_repeat_count (both in repeat
    units, from AlleleSize and ReferenceRegion/RepeatUnit) are always populated; is_near_reference
    is None unless exclude_near_reference is set.

    exclude_near_reference: alleles within NEAR_REFERENCE_REPEAT_TOLERANCE repeats of the reference
    allele size (recomputed per definition, since the reference region grows with the definition)
    are near-reference -- but how that's turned into an exclusion depends on whether the locus has
    an extended definition:
      - No extended definition (just "original"): no pairing is possible for this locus either way,
        so near-reference alleles are dropped directly here, per allele.
      - Has an extended definition: alleles are NOT dropped here. Every allele is returned, tagged
        with is_near_reference, and it's on the caller to drop a PAIRED (same allele_number, both
        "original" and "extended" present) row only when BOTH sides are near-reference. Dropping
        near-reference alleles independently per definition at this stage -- before pairing -- can
        let an allele that's near-reference under one definition but unpaired (no counterpart under
        the other definition, e.g. a no-call) admit the whole locus/sample past the filter, only for
        that same allele to then be dropped by the caller's pairing step anyway, silently leaving
        only near-reference pairs in the plotted output.
    """
    path = JSON_DIR / f"{sample_id}.combined_catalog.json.gz"
    with gzip.open(path) as f:
        d = json.load(f)

    by_base_id = {}
    for locus_id, result in d["LocusResults"].items():
        base_id, _, definition = locus_id.rpartition("__")
        variant = next(iter(result["Variants"].values()))
        if "AlleleQualityMetrics" not in variant:
            continue
        alleles = variant["AlleleQualityMetrics"]["Alleles"]
        _, region = variant["ReferenceRegion"].split(":")
        region_start, region_end = (int(x) for x in region.split("-"))
        ref_repeat_count = (region_end - region_start) / len(variant["RepeatUnit"])
        by_base_id.setdefault(base_id, {})[definition] = (alleles, ref_repeat_count)

    rows = []
    for base_id, by_definition in by_base_id.items():
        multi_definition = len(by_definition) > 1
        for definition, (alleles, ref_repeat_count) in by_definition.items():
            for allele in alleles:
                is_near_reference = abs(allele["AlleleSize"] - ref_repeat_count) <= NEAR_REFERENCE_REPEAT_TOLERANCE
                if exclude_near_reference and not multi_definition and is_near_reference:
                    continue  # single-definition locus: no pairing possible, filter directly
                rows.append((
                    base_id, definition, allele["AlleleNumber"], allele["pOk"],
                    allele["AlleleSize"], ref_repeat_count,
                    is_near_reference if exclude_near_reference else None,
                ))
    return rows


def pivot_paired(df, index_cols):
    """Pivot df on 'definition' to wide columns 'original'/'extended' for pOk, keeping only rows
    with a pOk under both definitions. When df's 'is_near_reference' column is populated
    (--exclude-near-reference mode), also drop pairs where BOTH sides are near-reference -- this
    applies the near-reference exclusion AFTER pairing, so a locus/sample can't be admitted by an
    allele that never itself survives pairing (see load_sample_pok's docstring)."""
    pok_wide = df.pivot_table(index=index_cols, columns="definition", values="pOk", aggfunc="first")
    wide = pok_wide.dropna(subset=["original", "extended"])
    if df["is_near_reference"].notna().any():
        near_wide = df.pivot_table(index=index_cols, columns="definition", values="is_near_reference", aggfunc="first")
        both_near_reference = near_wide.reindex(wide.index)[["original", "extended"]].all(axis=1)
        wide = wide[~both_near_reference]
    return wide.reset_index()


# Same width-2 (widening to width-5 beyond +/-20, open tails beyond +/-31) binning as
# str-truth-set-v2's tool_comparison_viewer.html plots of allele size minus reference repeat count
# (plot_tool_accuracy_by_allele_size.bin_num_repeats_wrapper(bin_size=2), non-exome case), reimplemented
# here rather than imported since it lives in a separate sibling repo.
DIFF_FROM_REF_BIN_ORDER = [
    "-31 or more", "-26 to -30", "-21 to -25",
    "-19 to -20", "-17 to -18", "-15 to -16", "-13 to -14", "-11 to -12", "-9 to -10", "-7 to -8", "-5 to -6", "-3 to -4", "-1 to -2",
    "0",
    "+1 to +2", "+3 to +4", "+5 to +6", "+7 to +8", "+9 to +10", "+11 to +12", "+13 to +14", "+15 to +16", "+17 to +18", "+19 to +20",
    "+21 to +25", "+26 to +30", "+31 or more",
]


def bin_diff_from_ref(diff):
    """Bin (AlleleSize - ref_repeat_count) into a DIFF_FROM_REF_BIN_ORDER label, or None if NaN."""
    if pd.isna(diff):
        return None
    if diff == 0:
        return "0"
    sign = "-" if diff < 0 else "+"
    diff = abs(diff)
    if 21 <= diff <= 25:
        return f"{sign}21 to {sign}25"
    if 26 <= diff <= 30:
        return f"{sign}26 to {sign}30"
    if diff >= 31:
        return f"{sign}31 or more"
    start = int((diff - 1) / 2) * 2
    return f"{sign}{int(start + 1)} to {sign}{int(start + 2)}"


def quantile_bin(df, value_col, n_bins=5):
    """Bin df[value_col] into n_bins quantiles, with edges computed from one value per
    base_locus_id rather than from df[value_col] directly -- df[value_col] is a per-row Series
    where a locus's value repeats once per surviving (sample, allele[, definition]) row, so
    quantiling it directly would give equal-ROW-count bins (skewed toward loci with more surviving
    rows), not the equal-LOCUS-count bins the "(quintile)" plot labels imply.
    """
    locus_values = df.drop_duplicates("base_locus_id")[value_col].dropna()
    try:
        bins = pd.qcut(locus_values, n_bins, duplicates="drop", retbins=True)[1]
        return pd.cut(df[value_col], bins=bins, include_lowest=True)
    except ValueError:
        bins = pd.qcut(locus_values.rank(method="first"), n_bins, retbins=True)[1]
        return pd.cut(df[value_col].rank(method="first"), bins=bins, include_lowest=True)


def make_violin_figure(groups_by_label, filename, xlabel, legend_labels=None, colors=(BLUE_LIGHT,), big=False):
    """groups_by_label: list of (x_label, [group_values...], [n_loci...]) -- 1 or more parallel groups
    (arrays) per bin, plus the number of distinct loci backing each series (for the "(n=...)" label --
    a series' value array has one entry per (locus, sample, allele), so len() would overcount loci).
    big: 1.5x taller figure and 2x label fonts (used for the Section 2 absolute-pOk plots).
    """
    n_series = len(groups_by_label[0][1])
    fig_height = 6.75 if big else 4.5
    tick_fontsize = 18 if big else 9
    label_fontsize = 22 if big else 11
    xticklabels = []
    for lab, group_vals, n_loci in groups_by_label:
        ns = "/".join(f"{n:,}" for n in n_loci)
        xticklabels.append(f"{lab}\n(n={ns})")
    # width needed so neighboring bins' (possibly two-line) tick labels don't collide, on top of the
    # existing per-series-count heuristic -- whichever demands more room wins.
    max_label_chars = max(max(len(line) for line in t.split("\n")) for t in xticklabels)
    label_width_in = max_label_chars * tick_fontsize * 0.6 / 72
    fig_width = max(6, len(groups_by_label) * max(1.4 * n_series, label_width_in + 1.2))
    # a plot with many bins (e.g. the 27-bin diff-from-ref plot) gets very wide at a fixed height,
    # which -- once scaled to the report's fixed-width .fig container -- makes it look squashed next
    # to narrower plots. Scale height up with width so all plots keep a similar on-screen aspect ratio.
    fig_height = max(fig_height, fig_width * 0.28)

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    fig.patch.set_alpha(0)
    ax.set_facecolor("none")

    width = 0.8 / n_series
    for series_i in range(n_series):
        positions = [i + 1 + (series_i - (n_series - 1) / 2) * width * 1.15 for i in range(len(groups_by_label))]
        data = [g[1][series_i] for g in groups_by_label]
        color = colors[series_i % len(colors)]
        parts = ax.violinplot(data, positions=positions, widths=width, showmedians=True, showextrema=True)
        for body in parts["bodies"]:
            body.set_facecolor(color)
            body.set_edgecolor(color)
            body.set_alpha(0.55)
        for key in ("cmedians", "cbars", "cmins", "cmaxes"):
            parts[key].set_color(color)
            parts[key].set_linewidth(1.2)

    ax.axhline(0, color="#898781", linewidth=1, linestyle="--", zorder=0) if n_series == 1 else None

    ax.set_xticks(range(1, len(groups_by_label) + 1))
    ax.set_xticklabels(xticklabels, fontsize=tick_fontsize)
    ax.set_xlabel(xlabel, fontsize=label_fontsize, color="#52514e", labelpad=16)
    ax.set_ylabel("pOk" if n_series > 1 else "Δ pOk  (extended − original)", fontsize=label_fontsize, color="#52514e")
    ax.tick_params(axis="x", colors="#52514e", labelsize=tick_fontsize, pad=10)
    ax.tick_params(axis="y", colors="#52514e", labelsize=tick_fontsize)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color("#e1e0d9")
    ax.grid(axis="y", color="#e1e0d9", linewidth=0.6, zorder=-1)
    if legend_labels:
        handles = [plt.Rectangle((0, 0), 1, 1, color=colors[i], alpha=0.55) for i in range(n_series)]
        # anchor the legend's bottom-right corner just outside the axes, bottom-right of the figure
        ax.legend(handles, legend_labels, fontsize=tick_fontsize, frameon=False, loc="lower right", bbox_to_anchor=(1.22, 0))
    fig.tight_layout()

    out_path = DATA_DIR / filename
    fig.savefig(out_path, format="svg", transparent=True, bbox_inches="tight")
    plt.close(fig)
    return out_path


def make_stacked_bar(labels, extended_counts, not_extended_counts, xlabel, filename):
    """Stacked bar chart: locus count per bin, split into extended (orange) / not extended (aqua)."""
    fig, ax = plt.subplots(figsize=(max(6, 1.3 * len(labels)), 4.5))
    fig.patch.set_alpha(0)
    ax.set_facecolor("none")
    x = range(1, len(labels) + 1)
    ax.bar(x, extended_counts, color=ORANGE_LIGHT, label="extended", zorder=3)
    ax.bar(x, not_extended_counts, bottom=extended_counts, color="#1baf7a", label="not extended", zorder=3)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_xlabel(xlabel, fontsize=11, color="#52514e", labelpad=16)
    ax.set_ylabel("Locus count", fontsize=11, color="#52514e")
    ax.tick_params(axis="x", colors="#52514e", pad=10)
    ax.tick_params(axis="y", colors="#52514e")
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color("#e1e0d9")
    ax.grid(axis="y", color="#e1e0d9", linewidth=0.6, zorder=0)
    ax.legend(fontsize=9, frameon=False, loc="upper right")
    fig.tight_layout()

    out_path = DATA_DIR / filename
    fig.savefig(out_path, format="svg", transparent=True)
    plt.close(fig)
    return out_path


def make_stacked_bar_by_group(df, group_col, group_values, label_fn, xlabel, filename):
    labels = [label_fn(v) for v in group_values]
    extended_counts = [int(((df[group_col] == v) & df["is_extended"]).sum()) for v in group_values]
    not_extended_counts = [int(((df[group_col] == v) & ~df["is_extended"]).sum()) for v in group_values]
    return make_stacked_bar(labels, extended_counts, not_extended_counts, xlabel, filename)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--exclude-near-reference", action="store_true",
                         help=f"drop alleles within +/-{NEAR_REFERENCE_REPEAT_TOLERANCE} repeats of the reference allele size")
    args = parser.parse_args()

    suffix = "_exclude_near_reference" if args.exclude_near_reference else ""
    output_tsv = DATA_DIR / f"pOk_by_locus_sample_allele{suffix}.tsv"
    output_html = Path(__file__).parent / f"build_pOk_delta_report{suffix}.html"

    print(f"Loading pOk from {len(SAMPLES)} sample JSONs...")
    all_rows = []
    for sample_id in SAMPLES:
        rows = load_sample_pok(sample_id, exclude_near_reference=args.exclude_near_reference)
        print(f"  {sample_id}: {len(rows):,} allele records")
        all_rows.extend((sample_id, *r) for r in rows)

    long_df = pd.DataFrame(all_rows, columns=["sample_id", "base_locus_id", "definition", "allele_number", "pOk",
                                               "allele_size", "ref_repeat_count", "is_near_reference"])
    long_df.to_csv(output_tsv, sep="\t", index=False)
    print(f"Wrote {len(long_df):,} rows to {output_tsv}")

    is_andrea = long_df["base_locus_id"].str.match(ANDREA_ID_PATTERN)
    andrea_long = long_df[is_andrea].copy()
    andrea_long["motif"] = andrea_long["base_locus_id"].str.rsplit("-", n=1).str[-1]
    andrea_long["motif_size"] = andrea_long["motif"].str.len()

    flank = pd.read_csv(
        ANDREA_CATALOG_PATH, sep="\t", usecols=["LocusId", "flank_motif_similarity_left", "flank_motif_similarity_right"]
    )
    andrea_long = andrea_long.merge(flank, left_on="base_locus_id", right_on="LocusId", how="left", validate="many_to_one")

    scan = pd.read_csv(ANDREA_SCAN_TSV, sep="\t")
    extended_ids = set(scan.loc[(scan["left_ext"] > 0) | (scan["right_ext"] > 0), "LocusId"])
    print(f"{len(extended_ids):,} genuinely-extended Andrea loci")

    # --- Section 1: how many of Andrea's 347,704 v2-exact-match loci were extended, by motif size
    # and by quintile of max(flank_motif_similarity_left, flank_motif_similarity_right) -- one row
    # per locus (not per sample). ---
    locus_extension = scan.merge(flank, on="LocusId", how="left", validate="one_to_one")
    locus_extension["motif_size"] = locus_extension["motif"].str.len()
    locus_extension["is_extended"] = locus_extension["LocusId"].isin(extended_ids)
    # one row per locus already, so no repeated-row weighting concern here -- qcut directly.
    locus_extension["fms_max"] = locus_extension[["flank_motif_similarity_left", "flank_motif_similarity_right"]].max(axis=1)
    locus_extension["fms_max_bin"] = pd.qcut(locus_extension["fms_max"], 5, duplicates="drop")

    fig_n_extended_motif = make_stacked_bar_by_group(
        locus_extension, "motif_size", sorted(locus_extension["motif_size"].dropna().unique()),
        lambda v: f"{int(v)}bp", "Motif size", "n_extended_by_motif_size.svg",
    )
    fig_n_extended_fms_max = make_stacked_bar_by_group(
        locus_extension, "fms_max_bin", sorted(locus_extension["fms_max_bin"].dropna().unique(), key=lambda c: c.left),
        lambda c: f"{c.left:.3f}–{c.right:.3f}", "max(flank_motif_similarity_left, flank_motif_similarity_right) (quintile)",
        "n_extended_by_flank_motif_similarity_max.svg",
    )

    # --- Section 3: delta pOk for genuinely-extended loci ---
    ext_wide = pivot_paired(
        andrea_long[andrea_long["base_locus_id"].isin(extended_ids)],
        index_cols=["base_locus_id", "sample_id", "allele_number", "motif_size", "flank_motif_similarity_left", "flank_motif_similarity_right"],
    )
    ext_wide["delta_pOk"] = ext_wide["extended"] - ext_wide["original"]
    print(f"\n{len(ext_wide):,} (locus, sample, allele) triples with delta pOk, across {ext_wide['base_locus_id'].nunique():,} loci")

    # x-axis for the diff-from-ref plot below: the ORIGINAL (pre-extension) definition's own
    # AlleleSize minus its own reference repeat count -- i.e. how repeat-expanded the allele already
    # is relative to the reference genome, same quantity/binning as str-truth-set-v2's
    # tool_comparison_viewer.html plots (which use the true/original reference, not an extended one).
    original_allele_size = andrea_long[andrea_long["definition"] == "original"][
        ["base_locus_id", "sample_id", "allele_number", "allele_size", "ref_repeat_count"]
    ]
    ext_wide = ext_wide.merge(original_allele_size, on=["base_locus_id", "sample_id", "allele_number"],
                               how="left", validate="one_to_one")
    ext_wide["diff_from_ref"] = ext_wide["allele_size"] - ext_wide["ref_repeat_count"]
    ext_wide["diff_from_ref_bin"] = ext_wide["diff_from_ref"].map(bin_diff_from_ref)

    ep400 = load_ep400(long_df)

    fig_motif_size = make_delta_violin_by_motif_size(ext_wide, ep400, f"delta_pOk_by_motif_size{suffix}.svg")
    ext_wide["fms_max"] = ext_wide[["flank_motif_similarity_left", "flank_motif_similarity_right"]].max(axis=1)
    ext_wide["fms_max_bin"] = quantile_bin(ext_wide, "fms_max")
    fig_fms_max = make_violin_figure(
        [(f"{c.left:.3f}–{c.right:.3f}",
          [ext_wide.loc[ext_wide["fms_max_bin"] == c, "delta_pOk"].dropna().values],
          [ext_wide.loc[ext_wide["fms_max_bin"] == c, "base_locus_id"].nunique()])
         for c in sorted(ext_wide["fms_max_bin"].dropna().unique(), key=lambda c: c.left)],
        f"delta_pOk_by_flank_motif_similarity_max{suffix}.svg", "max(flank_motif_similarity_left, flank_motif_similarity_right) (quintile)",
    )
    present_diff_bins = [b for b in DIFF_FROM_REF_BIN_ORDER if b in set(ext_wide["diff_from_ref_bin"].dropna())]
    fig_diff_from_ref = make_violin_figure(
        [(b,
          [ext_wide.loc[ext_wide["diff_from_ref_bin"] == b, "delta_pOk"].dropna().values],
          [ext_wide.loc[ext_wide["diff_from_ref_bin"] == b, "base_locus_id"].nunique()])
         for b in present_diff_bins],
        f"delta_pOk_by_diff_from_ref{suffix}.svg", "Original allele size minus reference repeat count",
        big=True,
    )

    # --- Section 2: absolute pOk, 3 groups -- original (pre-extension pOk of loci the rule DOES
    # extend), extended (post-extension pOk of those same loci), not_extended (the only pOk a locus
    # the rule never touches has). "original" and "not_extended" are both definition=="original" rows
    # but describe different populations, so they're kept as distinct groups, not merged.
    andrea_long["is_extended"] = andrea_long["base_locus_id"].isin(extended_ids)
    # for a genuinely-extended locus, only keep (locus, sample, allele) triples with a pOk under BOTH
    # definitions -- a sample where only one definition produced a call would otherwise contribute to
    # just one of "original"/"extended", making them unmatched populations rather than the same loci
    # compared before/after extension. Loci the rule never touches have only "original" and need no
    # such pairing.
    key_cols = ["base_locus_id", "sample_id", "allele_number"]
    paired_index = pd.MultiIndex.from_frame(
        pivot_paired(andrea_long[andrea_long["is_extended"]], index_cols=key_cols)[key_cols]
    )
    is_paired = andrea_long.set_index(key_cols).index.isin(paired_index)
    abs_df = andrea_long[~andrea_long["is_extended"] | is_paired].copy()
    abs_df["group"] = np.where(abs_df["is_extended"], abs_df["definition"], "not_extended")
    group_counts = abs_df["group"].value_counts()
    print(f"rows for the 3-group absolute-pOk comparison: {dict(group_counts)}")

    fig_abs_motif = make_abs_violin(abs_df, "motif_size", sorted(abs_df["motif_size"].dropna().unique()),
                                     lambda v: f"{int(v)}bp", "Motif size", f"abs_pOk_by_motif_size{suffix}.svg")
    abs_df["fms_max"] = abs_df[["flank_motif_similarity_left", "flank_motif_similarity_right"]].max(axis=1)
    abs_df["fms_max_bin"] = quantile_bin(abs_df, "fms_max")
    fig_abs_fms_max = make_abs_violin(
        abs_df, "fms_max_bin", sorted(abs_df["fms_max_bin"].dropna().unique(), key=lambda c: c.left),
        lambda c: f"{c.left:.3f}–{c.right:.3f}", "max(flank_motif_similarity_left, flank_motif_similarity_right) (quintile)",
        f"abs_pOk_by_flank_motif_similarity_max{suffix}.svg",
    )

    write_html_report(ext_wide, ep400, abs_df, fig_motif_size, fig_fms_max, fig_diff_from_ref, fig_abs_motif, fig_abs_fms_max,
                       fig_n_extended_motif, fig_n_extended_fms_max,
                       len(scan), output_html, args.exclude_near_reference)
    print(f"\nWrote {output_html}")


def load_ep400(long_df):
    ep400 = long_df[long_df["base_locus_id"] == "EP400"]
    wide = pivot_paired(ep400, index_cols=["sample_id", "allele_number"])
    wide["delta_pOk"] = wide["extended"] - wide["original"]
    return wide


def make_delta_violin_by_motif_size(ext_wide, ep400, filename):
    sizes = sorted(ext_wide["motif_size"].dropna().unique())
    groups = [ext_wide.loc[ext_wide["motif_size"] == s, "delta_pOk"].dropna().values for s in sizes]
    n_loci = [ext_wide.loc[ext_wide["motif_size"] == s, "base_locus_id"].nunique() for s in sizes]
    labels = [f"{s}bp\n(n={n:,})" for s, n in zip(sizes, n_loci)]
    colors = [BLUE_LIGHT] * len(groups)

    # EP400 (motif CAG = 3bp): its own violin, inserted as a category right after the 3bp bin --
    # since it's a single locus and not part of Andrea's catalog. Slotted into the same uniform
    # position/width grid as the bp bins (not squeezed between them), so it can't overlap its
    # neighbors any more than the bp bins overlap each other.
    ep400_vals = ep400["delta_pOk"].dropna().values
    if 3 in sizes and len(ep400_vals) >= 2:
        insert_at = sizes.index(3) + 1
        groups.insert(insert_at, ep400_vals)
        labels.insert(insert_at, "EP400\n(n=1)")
        colors.insert(insert_at, ORANGE_LIGHT)

    fig, ax = plt.subplots(figsize=(max(6, 1.3 * len(groups)), 4.5))
    fig.patch.set_alpha(0)
    ax.set_facecolor("none")
    positions = list(range(1, len(groups) + 1))
    for pos, vals, color in zip(positions, groups, colors):
        parts = ax.violinplot([vals], positions=[pos], showmedians=True, showextrema=True, widths=0.8)
        for body in parts["bodies"]:
            body.set_facecolor(color); body.set_edgecolor(color); body.set_alpha(0.55)
        for key in ("cmedians", "cbars", "cmins", "cmaxes"):
            parts[key].set_color(color); parts[key].set_linewidth(1.2)
    ax.axhline(0, color="#898781", linewidth=1, linestyle="--", zorder=0)

    ax.set_xticks(positions)
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_xlabel("Motif size", fontsize=11, color="#52514e", labelpad=16)
    ax.set_ylabel("Δ pOk  (extended − original)", fontsize=11, color="#52514e")
    ax.tick_params(axis="x", colors="#52514e", pad=10)
    ax.tick_params(axis="y", colors="#52514e")
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color("#e1e0d9")
    ax.grid(axis="y", color="#e1e0d9", linewidth=0.6, zorder=-1)
    fig.tight_layout()

    out_path = DATA_DIR / filename
    fig.savefig(out_path, format="svg", transparent=True)
    plt.close(fig)
    return out_path


ABS_GROUPS = ["original", "extended", "not_extended"]
ABS_GROUP_LABELS = ["original (pre-extension)", "extended (post-extension)", "not extended"]
ABS_GROUP_COLORS = (BLUE_LIGHT, ORANGE_LIGHT, "#1baf7a")  # fixed categorical order: blue, orange, aqua


def make_abs_violin(df, group_col, group_values, label_fn, xlabel, filename):
    groups_by_label = []
    for v in group_values:
        sub = df[df[group_col] == v]
        vals = [sub.loc[sub["group"] == g, "pOk"].dropna().values for g in ABS_GROUPS]
        n_loci = [sub.loc[sub["group"] == g, "base_locus_id"].nunique() for g in ABS_GROUPS]
        groups_by_label.append((label_fn(v), vals, n_loci))

    return make_violin_figure(groups_by_label, filename, xlabel, legend_labels=ABS_GROUP_LABELS, colors=ABS_GROUP_COLORS, big=True)


def write_html_report(ext_wide, ep400, abs_df, fig_motif_size, fig_fms_max, fig_diff_from_ref, fig_abs_motif, fig_abs_fms_max,
                       fig_n_extended_motif, fig_n_extended_fms_max,
                       n_scanned_loci, output_html, exclude_near_reference):
    def svg_inline(path):
        return path.read_text()

    title_suffix = ", near-reference alleles excluded" if exclude_near_reference else ""
    subtitle_suffix = (f" &middot; <strong>near-reference alleles excluded</strong> -- alleles within "
                        f"&plusmn;{NEAR_REFERENCE_REPEAT_TOLERANCE} repeats of the reference allele size "
                        f"(recomputed per locus definition) are dropped; both short and long alleles of a "
                        f"genotype are otherwise kept") if exclude_near_reference else ""

    overall_mean = ext_wide["delta_pOk"].mean()
    overall_median = ext_wide["delta_pOk"].median()
    locus_mean_delta = ext_wide.groupby("base_locus_id")["delta_pOk"].mean()
    n_loci = len(locus_mean_delta)
    n_better = (locus_mean_delta > 0.05).sum()
    n_worse = (locus_mean_delta < -0.05).sum()
    n_flat = n_loci - n_better - n_worse

    motif_size_table_rows = ""
    for size, sub in sorted(ext_wide.groupby("motif_size")):
        motif_size_table_rows += (
            f"<tr><td>{int(size)}bp</td><td class='mono'>{sub['base_locus_id'].nunique():,}</td>"
            f"<td class='mono'>{sub['delta_pOk'].mean():+.4f}</td>"
            f"<td class='mono'>{sub['delta_pOk'].median():+.4f}</td></tr>\n"
        )
        if size == 3:
            motif_size_table_rows += (
                f"<tr><td>EP400</td><td class='mono'>1</td>"
                f"<td class='mono'>{ep400['delta_pOk'].mean():+.4f}</td>"
                f"<td class='mono'>{ep400['delta_pOk'].median():+.4f}</td></tr>\n"
            )

    fms_max_table_rows = ""
    for c, sub in sorted(ext_wide.groupby("fms_max_bin", observed=True), key=lambda item: item[0].left):
        fms_max_table_rows += (
            f"<tr><td>{c.left:.3f}–{c.right:.3f}</td><td class='mono'>{sub['base_locus_id'].nunique():,}</td>"
            f"<td class='mono'>{sub['delta_pOk'].mean():+.4f}</td>"
            f"<td class='mono'>{sub['delta_pOk'].median():+.4f}</td></tr>\n"
        )

    diff_from_ref_table_rows = ""
    for b in [b for b in DIFF_FROM_REF_BIN_ORDER if b in set(ext_wide["diff_from_ref_bin"].dropna())]:
        sub = ext_wide[ext_wide["diff_from_ref_bin"] == b]
        diff_from_ref_table_rows += (
            f"<tr><td>{b}</td><td class='mono'>{sub['base_locus_id'].nunique():,}</td>"
            f"<td class='mono'>{sub['delta_pOk'].mean():+.4f}</td>"
            f"<td class='mono'>{sub['delta_pOk'].median():+.4f}</td></tr>\n"
        )

    abs_group_stats = abs_df.groupby("group").agg(n_loci=("base_locus_id", "nunique"), mean=("pOk", "mean")).reindex(ABS_GROUPS)

    html = f"""<title>pOk: Original vs. Extended Definitions ({len(SAMPLES)}-sample, corrected{title_suffix})</title>
<style>
  :root {{
    --surface: #fcfcfb; --page: #f9f9f7; --ink: #0b0b0b; --ink2: #52514e; --muted: #898781;
    --grid: #e1e0d9; --border: rgba(11,11,11,0.10); --blue: #2a78d6; --good: #0ca30c; --crit: #d03b3b;
  }}
  @media (prefers-color-scheme: dark) {{
    :root:not([data-theme="light"]) {{
      --surface: #1a1a19; --page: #0d0d0d; --ink: #ffffff; --ink2: #c3c2b7; --muted: #898781;
      --grid: #2c2c2a; --border: rgba(255,255,255,0.10); --blue: #3987e5;
    }}
  }}
  :root[data-theme="dark"] {{
    --surface: #1a1a19; --page: #0d0d0d; --ink: #ffffff; --ink2: #c3c2b7; --muted: #898781;
    --grid: #2c2c2a; --border: rgba(255,255,255,0.10); --blue: #3987e5;
  }}
  * {{ box-sizing: border-box; }}
  body {{
    background: var(--page); color: var(--ink); font-family: system-ui, -apple-system, "Segoe UI", sans-serif;
    line-height: 1.55; margin: 0; padding: 0;
  }}
  .wrap {{ max-width: 980px; margin: 0 auto; padding: 40px 24px 120px; }}
  h1 {{ font-size: 1.7rem; margin-bottom: 4px; }}
  h2 {{ font-size: 1.25rem; margin-top: 40px; border-bottom: 1px solid var(--grid); padding-bottom: 8px; }}
  .subtitle {{ color: var(--ink2); font-size: 1rem; margin-top: 0; }}
  .meta {{ color: var(--muted); font-size: 0.85rem; }}
  .card {{ background: var(--surface); border: 1px solid var(--border); border-radius: 10px; padding: 20px 24px; margin: 18px 0; }}
  .callout {{ background: var(--surface); border-left: 4px solid var(--blue); border-radius: 6px; padding: 14px 18px; margin: 16px 0; }}
  .callout.warn {{ border-left-color: #eda100; }}
  table {{ border-collapse: collapse; width: 100%; margin: 14px 0; font-size: 0.92rem; }}
  th, td {{ text-align: left; padding: 7px 10px; border-bottom: 1px solid var(--grid); font-variant-numeric: tabular-nums; }}
  th {{ color: var(--ink2); font-weight: 600; border-bottom: 1px solid var(--muted); }}
  code, .mono {{ font-family: ui-monospace, "SF Mono", Menlo, monospace; font-size: 0.88em; }}
  .fig {{ background: var(--surface); border: 1px solid var(--border); border-radius: 10px; padding: 16px; margin: 16px 0; overflow-x: auto; cursor: zoom-in; }}
  .fig svg {{ max-width: 100%; height: auto; }}
  footer {{ color: var(--muted); font-size: 0.82rem; margin-top: 60px; }}
  .modal-overlay {{
    display: none; position: fixed; inset: 0; background: rgba(11,11,11,0.75); z-index: 1000;
    align-items: center; justify-content: center; padding: 24px; cursor: zoom-out;
  }}
  .modal-overlay.open {{ display: flex; }}
  .modal-content {{
    background: var(--surface); border-radius: 10px; padding: 24px; max-width: 95vw; max-height: 95vh;
    overflow: auto; cursor: default;
  }}
  .modal-content svg {{ max-width: 88vw; max-height: 85vh; width: auto; height: auto; }}
  .modal-close {{
    position: fixed; top: 20px; right: 28px; color: white; font-size: 2rem; line-height: 1;
    cursor: pointer; background: none; border: none; z-index: 1001;
  }}
</style>

<div class="wrap">
<h1>pOk: Original vs. Extended Definitions</h1>

<div class="callout">
Evaluating <strong>{ext_wide['base_locus_id'].nunique():,}</strong> extended locus definitions, genotyping
them (and the original definitions) in {len(SAMPLES)} genome samples from 1kGP.
Mean &Delta;pOk = <strong>{overall_mean:+.4f}</strong>, median = <strong>{overall_median:+.4f}</strong>.
By per-locus mean &Delta;pOk (averaged across its samples/alleles) -- Better (&gt;+0.05): {n_better:,} loci ({100*n_better/n_loci:.1f}%).
Worse (&lt;&minus;0.05): {n_worse:,} loci ({100*n_worse/n_loci:.1f}%). Flat: {n_flat:,} loci ({100*n_flat/n_loci:.1f}%).
</div>
<p class="subtitle">{len(SAMPLES)} 1kGP samples (up to 3 per Population/Gender group){subtitle_suffix}</p>

<h2 id="counts">1. How many loci were extended</h2>
<p class="meta">Locus counts (not sample-dependent) across all {n_scanned_loci:,} Andrea v2-exact-match loci, split into extended vs. not extended by the gap-purity rule.</p>

<h3>By motif size</h3>
<div class="fig" onclick="openModal(this)">{svg_inline(fig_n_extended_motif)}</div>

<h3>By max(flank_motif_similarity_left, flank_motif_similarity_right) (quintiles)</h3>
<div class="fig" onclick="openModal(this)">{svg_inline(fig_n_extended_fms_max)}</div>

<h2 id="absolute">2. Absolute pOk: original vs. extended vs. not-extended</h2>
<div class="callout">
Three groups: <strong>original</strong> = pOk of locus definitions pre-extension;
<strong>extended</strong> = pOk of the same loci post-extension; <strong>not extended</strong> =
pOk of loci whose boundaries were not extended.
<table>
<tr><th>Group</th><th>n</th><th>mean pOk</th></tr>
{"".join(f"<tr><td>{label}</td><td class='mono'>{int(abs_group_stats.loc[g,'n_loci']):,}</td><td class='mono'>{abs_group_stats.loc[g,'mean']:.4f}</td></tr>" for g, label in zip(ABS_GROUPS, ABS_GROUP_LABELS))}
</table>
</div>

<h3>By motif size</h3>
<div class="fig" onclick="openModal(this)">{svg_inline(fig_abs_motif)}</div>

<h3>By max(flank_motif_similarity_left, flank_motif_similarity_right) (quintiles)</h3>
<div class="fig" onclick="openModal(this)">{svg_inline(fig_abs_fms_max)}</div>

<h2 id="delta">3. Plots comparing ExpansionHunter genotype quality for extended loci vs. their original definitions (&Delta;pOk)</h2>

<h3>By motif size</h3>
<div class="fig" onclick="openModal(this)">{svg_inline(fig_motif_size)}</div>
<table>
<tr><th>Motif size</th><th>n</th><th>mean &Delta;pOk</th><th>median &Delta;pOk</th></tr>
{motif_size_table_rows}
</table>

<h3>By max(flank_motif_similarity_left, flank_motif_similarity_right) (original locus definition, quintiles)</h3>
<p class="meta">How closely the ORIGINAL locus's flanks already match a perfect repeat of its own motif, before extension (max of left/right).</p>
<div class="fig" onclick="openModal(this)">{svg_inline(fig_fms_max)}</div>
<table>
<tr><th>max(flank_motif_similarity_left, flank_motif_similarity_right)</th><th>n</th><th>mean &Delta;pOk</th><th>median &Delta;pOk</th></tr>
{fms_max_table_rows}
</table>

<h3>By original allele size minus reference repeat count</h3>
<p class="meta">Computed from the ORIGINAL (pre-extension) locus definition's own AlleleSize and reference repeat count.</p>
<div class="fig" onclick="openModal(this)">{svg_inline(fig_diff_from_ref)}</div>
<table>
<tr><th>Original allele size &minus; reference repeat count</th><th>n</th><th>mean &Delta;pOk</th><th>median &Delta;pOk</th></tr>
{diff_from_ref_table_rows}
</table>

</div>

<div class="modal-overlay" id="modal" onclick="closeModal()">
<button class="modal-close" onclick="closeModal()">&times;</button>
<div class="modal-content" id="modal-content" onclick="event.stopPropagation()"><div id="modal-svg-container"></div></div>
</div>

<script>
function openModal(figElement) {{
  const svg = figElement.querySelector("svg");
  if (!svg) return;
  document.getElementById("modal-svg-container").innerHTML = "";
  document.getElementById("modal-svg-container").appendChild(svg.cloneNode(true));
  document.getElementById("modal").classList.add("open");
}}
function closeModal() {{
  document.getElementById("modal").classList.remove("open");
}}
document.addEventListener("keydown", function(event) {{
  if (event.key === "Escape") closeModal();
}});
</script>
"""
    output_html.write_text(html)


if __name__ == "__main__":
    main()
