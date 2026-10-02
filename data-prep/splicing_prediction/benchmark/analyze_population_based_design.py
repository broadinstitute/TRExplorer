"""Which simulated allele sizes are worth keeping: analysis of the population-based design benchmark.

Reads the benchmark sample (make_tr_expansion_benchmark_alleles.py: random polymorphic GENCODE v50 basic
loci, each with an allele for every target length the sample was made with, which can be more than the
current design's TARGET_LABELS) and its L4 scores (spliceai_l4_cost_and_disk_benchmark.py), and reports:
- for each target length, the share of loci where its allele changes splicing by >= 0.2 (largest
  SpliceAI delta score on the selected, stored transcript), the share found by no other target, and
  what dropping it would save and lose
- the best set of targets of each size, by loci found, and its full-run cost and disk use
- the full-run cost, time and disk estimate for every target in the sample, and for the current design
- plots: the share of loci with a change by distance to the nearest splice site, the cost and loci
  found of every set of targets, and the cumulative cost and loci found as the distance limit grows
and writes them to population_based_design_report.html.

Cost of a set of targets is the benchmark's full-run estimate times the share of scoring seconds the
set keeps. Each locus's first scored allele also paid for the shared REF prediction, which every set
that keeps an allele at the locus needs; when such a set drops that allele, the REF part (its time
minus the median time of the locus's other alleles) is still counted. Loci with no allele for a target (e.g. its length rounds to the reference) count as not
found by it.

Usage:
    python3 analyze_population_based_design.py
"""
import collections
import csv
import gzip
import itertools
import json
import math
import os
import sys

import markdown

# The pipeline modules are in the parent folder
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from make_tr_expansion_benchmark_alleles import (  # noqa: E402
    DISTANCE_COLUMN, DISTANCES_TSV, TARGET_LABELS, load_polymorphic_scored_loci, short_target_label, target_sort_key)

HERE = os.path.dirname(os.path.abspath(__file__))
ALLELES_JSON = os.path.join(HERE, "tr_benchmark_alleles.basic.json")
RESULTS_JSON = os.path.join(HERE, "l4_benchmark_results.basic.json")
ESTIMATE_JSON = os.path.join(HERE, "l4_benchmark_estimate.json")
REPORT_HTML = os.path.join(HERE, "population_based_design_report.html")
STORAGE_FORMAT = "deltas_with_ref_alt"
CUTOFFS = (0.2, 0.5)
MAX_CONCURRENT_GPUS = 10
# Coarser at both ends than in summarize_splicing_effects_by_distance.py: a random sample of polymorphic
# loci has only a few loci at distance 0 or more than 200 kb away.
DISTANCE_BINS = [(0, 50), (51, 188), (189, 500), (501, 2000), (2001, 10000), (10001, 50000), (50001, 10 ** 9)]
COLORS = ("#2a78d6", "#eb6834", "#1baf7a")  # reference palette, categorical slots 1-3, in fixed order


def load_scored_loci():
    """Returns [{"locus_id", "alleles": [(set of target labels, largest delta, seconds, stored bytes, variant)]}]."""
    sample = json.load(open(ALLELES_JSON))
    # Overlapping loci can have identical allele lists, so each list maps to all of its loci, used one at a time
    by_variants = collections.defaultdict(list)
    for locus_id, alleles, labels in zip(sample["locus_ids"], sample["alleles_by_locus"], sample["target_labels_by_locus"]):
        by_variants[tuple(alleles)].append((locus_id, labels))
    loci = []
    for shard in json.load(open(RESULTS_JSON)):
        for results in shard["results_by_locus"]:
            locus_id, labels = by_variants[tuple(r["variant"] for r in results)].pop()
            loci.append({"locus_id": locus_id, "alleles": [
                (set(allele_labels), largest_delta_score(r), r["seconds"],
                 len(r["storage_records"][STORAGE_FORMAT]) + 1 if r["storage_records"][STORAGE_FORMAT] else 0, r["variant"])
                for allele_labels, r in zip(labels, results)]})
    return loci


def largest_delta_score(result):
    """The largest of the four delta scores on the transcript the server selects (the one stored)."""
    return result["max_delta_score_selected_transcript"]


def wilson_interval(successes, n, z=1.96):
    """95% Wilson score interval for a binomial proportion."""
    p = successes / n
    center = (p + z * z / (2 * n)) / (1 + z * z / n)
    half_width = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return center - half_width, center + half_width


def short_label(target):
    """e.g. "2.5th pct", "+1x"."""
    return short_target_label(target) + ("" if " + " in target else " pct")


def is_observed_length(target):
    """True for the HPRC256 percentile lengths, False for the 99.5th + k x motif range sizes."""
    return " + " not in target


def is_found(locus, kept_targets, cutoff):
    """True if an allele for one of kept_targets changes splicing by at least cutoff at this locus."""
    return any(delta >= cutoff for labels, delta, _, _, _ in locus["alleles"] if labels & kept_targets)


def evaluate(loci, kept_targets, cutoff=CUTOFFS[0]):
    """Returns (loci found, scoring seconds, stored bytes) for the alleles of kept_targets."""
    kept_targets = set(kept_targets)
    found = seconds = stored_bytes = 0
    for locus in loci:
        alleles = locus["alleles"]
        found += is_found(locus, kept_targets, cutoff)
        kept = [a for a in alleles if a[0] & kept_targets]
        seconds += sum(a[2] for a in kept)
        if kept and not alleles[0][0] & kept_targets:
            # The first allele scored also paid for the REF prediction, which every set of targets that keeps
            # an allele at this locus needs: keep that part, estimated as its time minus the median time of
            # the locus's other alleles.
            others = sorted(a[2] for a in alleles[1:])
            seconds += max(0.0, alleles[0][2] - others[len(others) // 2]) if others else alleles[0][2]
        stored_bytes += sum(a[3] for a in kept)
    return found, seconds, stored_bytes


def load_distances(locus_ids):
    """Returns {locus_id: distance to the nearest GENCODE v50 basic splice site} for the given loci."""
    wanted = set(locus_ids)
    distances = {}
    with gzip.open(DISTANCES_TSV, "rt") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            locus_id = f"{row['chrom']}-{row['start_0based']}-{row['end']}-{row['motif']}"
            distance = row[DISTANCE_COLUMN["basic"]]
            if locus_id in wanted and distance != ".":
                distances[locus_id] = int(distance)
    return distances


def distance_bin(distance):
    return next(b for b in DISTANCE_BINS if b[0] <= distance <= b[1])


def bin_label(b):
    return f"{b[0]:,}" if b[0] == b[1] else (f"> {b[0] - 1:,}" if b[1] == 10 ** 9 else f"{b[0]:,}-{b[1]:,}")


def percent(x):
    return f"{100 * x:.1f}%"


def markdown_table(header, rows):
    return "\n".join(["| " + " | ".join(header) + " |", "|" + "---|" * len(header)] +
                     ["| " + " | ".join(str(c) for c in row) + " |" for row in rows])


def style_axes(ax):
    ax.grid(axis="y", color="#e8e7e3", linewidth=0.8)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)


def plot_fraction_by_distance(loci_by_bin, targets, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    observed = [t for t in targets if is_observed_length(t)]
    motif_range_sizes = [t for t in targets if not is_observed_length(t)]
    series = [("any target", set(targets)),
              (f"observed lengths ({', '.join(short_target_label(t) for t in observed)} pct)", set(observed)),
              (f"motif-range sizes ({', '.join(short_target_label(t) for t in motif_range_sizes)})", set(motif_range_sizes))]
    bins = [b for b in DISTANCE_BINS if loci_by_bin.get(b)]
    fig, axes = plt.subplots(len(CUTOFFS), 1, figsize=(10, 8), dpi=150, sharex=True)
    for ax, cutoff in zip(axes, CUTOFFS):
        for (name, targets), color, offset in zip(series, COLORS, (-0.2, 0, 0.2)):
            x, y, lower, upper = [], [], [], []
            for i, b in enumerate(bins):
                n_found = sum(is_found(locus, targets, cutoff) for locus in loci_by_bin[b])
                fraction = n_found / len(loci_by_bin[b])
                low, high = wilson_interval(n_found, len(loci_by_bin[b]))
                x.append(i + offset)
                y.append(fraction)
                # max(0, ...): at a fraction of 0 or 1 the interval's edge can differ from it by a float rounding error
                lower.append(max(0, fraction - low))
                upper.append(max(0, high - fraction))
            ax.errorbar(x, y, yerr=[lower, upper], color=color, fmt="o", markersize=7, markeredgecolor="white",
                        markeredgewidth=1.5, elinewidth=1.2, capsize=2.5, label=name)
        ax.set_title(f"change >= {cutoff}", fontsize=10, color="#0b0b0b", loc="left")
        ax.set_ylabel(f"Loci with a change >= {cutoff}", color="#0b0b0b", fontsize=9)
        ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0))
        ax.set_ylim(0, None)
        style_axes(ax)
        ax.legend(frameon=False, fontsize=9, loc="upper right")
    axes[-1].set_xticks(range(len(bins)))
    axes[-1].set_xticklabels([f"{bin_label(b)}\n(n={len(loci_by_bin[b])})" for b in bins], fontsize=8)
    axes[-1].set_xlabel("Distance from the TR locus to the nearest GENCODE v50 basic splice site (bp)", color="#0b0b0b")
    fig.suptitle("Polymorphic loci where a simulated allele changes splicing, by distance to a splice site\n"
                 "(largest SpliceAI delta on the selected transcript; error bars: 95% interval)", fontsize=11, color="#0b0b0b")
    fig.tight_layout()
    fig.savefig(path)


def plot_target_sets(points, best_by_size, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(10, 6), dpi=150)
    ax.scatter([p[1] for p in points], [p[0] for p in points], s=18, color="#b9b8b3", label="every set of targets")
    ax.plot([b[1] for b in best_by_size], [b[0] for b in best_by_size], color=COLORS[0], linewidth=2, marker="o",
            markersize=8, markeredgecolor="white", label="best set with N targets (N above each point)")
    # The sets themselves are listed in the report's table; points are labeled by their number of targets.
    for found, dollars, targets in best_by_size:
        ax.annotate(f"{len(targets)}", (dollars, found), textcoords="offset points", xytext=(0, 9), ha="center",
                    fontsize=9, color="#0b0b0b")
    ax.set_xlabel("Full-run cost on Modal L4s ($)", color="#0b0b0b")
    ax.set_ylabel("Polymorphic loci with a change >= 0.2 (full run)", color="#0b0b0b")
    ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:,.0f}"))
    style_axes(ax)
    ax.legend(frameon=False, fontsize=9, loc="lower right")
    ax.set_title("Loci found vs cost for every set of target sizes", fontsize=11, color="#0b0b0b", loc="left")
    fig.tight_layout()
    fig.savefig(path)


def plot_cumulative_by_distance(rows, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(10, 6), dpi=150)
    for (name, key_found, key_dollars, offset), color in zip((("every target", "found", "dollars", (6, -12)),
                                                               ("best 4 targets", "found_best4", "dollars_best4", (-6, 6))), COLORS):
        ax.plot([r[key_dollars] for r in rows], [r[key_found] for r in rows], color=color, linewidth=2, marker="o",
                markersize=7, markeredgecolor="white", label=name)
        for r in rows:
            ax.annotate(f"≤ {r['threshold_label']}", (r[key_dollars], r[key_found]), textcoords="offset points",
                        xytext=offset, ha="left" if offset[0] > 0 else "right", fontsize=7, color="#3d3d3a")
    ax.set_xlabel("Cumulative full-run cost ($)", color="#0b0b0b")
    ax.set_ylabel("Cumulative polymorphic loci with a change >= 0.2", color="#0b0b0b")
    ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:,.0f}"))
    style_axes(ax)
    ax.legend(frameon=False, fontsize=9, loc="lower right")
    ax.set_title("Loci found vs cost as the distance-to-splice-site limit grows (bp)", fontsize=11, color="#0b0b0b", loc="left")
    fig.tight_layout()
    fig.savefig(path)


def count_polymorphic_loci_by_distance_bin():
    """Returns {distance bin: number of polymorphic loci in the full run}."""
    _, polymorphic_loci = load_polymorphic_scored_loci("basic")
    distances = load_distances(f"{c}-{s}-{e}-{m}" for c, s, e, m, _ in polymorphic_loci)
    return collections.Counter(distance_bin(d) for d in distances.values())


def main():
    loci = load_scored_loci()
    estimate = json.load(open(ESTIMATE_JSON))["basic"]
    full = estimate["estimate"]
    n_full_run = full["n_loci_in_full_run"]
    targets = sorted({t for locus in loci for allele in locus["alleles"] for t in allele[0]}, key=target_sort_key)
    all_targets = set(targets)
    found_all, seconds_all, bytes_all = evaluate(loci, all_targets)
    dollars_per_second = full["full_run_dollars"] / seconds_all
    gb_per_byte = full[STORAGE_FORMAT]["full_run_gb_gzipped"] / bytes_all
    scale = n_full_run / len(loci)
    sections = []

    # 1. The whole design
    low, high = wilson_interval(found_all, len(loci))
    strong = evaluate(loci, all_targets, CUTOFFS[1])[0]
    equivalence = estimate["equivalence_check"]
    sections.append(f"## Full-run estimate (all {len(targets)} target lengths in the benchmark)\n\n" + markdown_table(
        ["", "Value"], [
            ["Polymorphic GENCODE v50 basic loci in the full run", f"{n_full_run:,}"],
            ["Sampled loci / alleles", f"{len(loci):,} / {sum(len(l['alleles']) for l in loci):,}"],
            ["Full-run alleles", f"{full['n_full_run_alleles']:,}"],
            ["Loci with a change >= 0.2", f"{percent(found_all / len(loci))} (95% interval {percent(low)}-{percent(high)}), "
                                          f"about {found_all * scale:,.0f} loci"],
            ["Loci with a change >= 0.5", f"{percent(strong / len(loci))}, about {strong * scale:,.0f} loci"],
            ["Cost on Modal L4s", f"${full['full_run_dollars']:,.0f} (90% range ${full['full_run_dollars_90pct_range'][0]:,.0f}-"
                                  f"${full['full_run_dollars_90pct_range'][1]:,.0f})"],
            ["L4 container-hours / wall time at 10 GPUs", f"{full['full_run_l4_container_hours']:,.0f} h / "
                                                          f"{full['full_run_l4_container_hours'] / MAX_CONCURRENT_GPUS:,.0f} h"],
            ["Disk, gzipped (deltas + REF/ALT)", f"{full[STORAGE_FORMAT]['full_run_gb_gzipped']:.2f} GB"],
            ["Errors / optimized code identical to original",
             f"{percent(full['fraction_of_alleles_with_errors'])} / {equivalence['n_identical_full_output']} of {equivalence['n_alleles']} alleles"],
        ]))

    # 2. Each target length's contribution
    rows = []
    for target in targets:
        has_allele = sum(any(target in a[0] for a in locus["alleles"]) for locus in loci)
        found = sum(is_found(locus, {target}, CUTOFFS[0]) for locus in loci)
        only = sum(is_found(locus, {target}, CUTOFFS[0]) and not is_found(locus, all_targets - {target}, CUTOFFS[0]) for locus in loci)
        found_without, seconds_without, _ = evaluate(loci, all_targets - {target})
        rows.append([short_label(target), f"{has_allele:,}", percent(found / len(loci)), f"{found * scale:,.0f}",
                     percent(only / len(loci)), f"{only * scale:,.0f}", f"{percent(only / found_all)}",
                     f"${(seconds_all - seconds_without) * dollars_per_second:,.0f}"])
    sections.append(
        "## What each target length contributes\n\n"
        "**Finds**: loci where this target's allele changes splicing by >= 0.2. **Only this target**: loci no other "
        "target finds, i.e. what dropping it loses. **Saves if dropped**: full-run cost of its alleles. Loci counts "
        f"are scaled to the full run; they rest on {found_all} found loci in the sample, so a count from a handful of "
        "loci (e.g. under 1,000 in the full run, about 3 sampled) is very uncertain.\n\n" + markdown_table(
            ["Target", "Sampled loci with this allele", "Finds", "Finds (full run)", "Only this target",
             "Only this target (full run)", "Share of all found loci", "Saves if dropped"], rows))

    # 3. Best set of each size
    points, best_by_size, best_rows = [], [], []
    for size in range(1, len(targets) + 1):
        candidates = []
        for target_set in itertools.combinations(targets, size):
            found, seconds, stored_bytes = evaluate(loci, target_set)
            points.append((found * scale, seconds * dollars_per_second))
            candidates.append((found, -seconds, target_set, stored_bytes))
        found, neg_seconds, target_set, stored_bytes = max(candidates)
        best_by_size.append((found * scale, -neg_seconds * dollars_per_second, set(target_set)))
        best_rows.append([size, ", ".join(short_label(t) for t in target_set), percent(found / len(loci)),
                          f"{found * scale:,.0f}", percent(found / found_all), f"${-neg_seconds * dollars_per_second:,.0f}",
                          f"{stored_bytes * gb_per_byte:.2f}"])
    # The design make_tr_expansion_benchmark_alleles.py now uses, if this benchmark scored all its targets
    if set(TARGET_LABELS) <= all_targets:
        found, seconds, stored_bytes = evaluate(loci, TARGET_LABELS)
        best_rows.append([f"current design ({len(TARGET_LABELS)})", ", ".join(short_label(t) for t in TARGET_LABELS),
                          percent(found / len(loci)), f"{found * scale:,.0f}", percent(found / found_all),
                          f"${seconds * dollars_per_second:,.0f}", f"{stored_bytes * gb_per_byte:.2f}"])
    sections.append(
        "## Best set of targets of each size\n\nChosen by loci found (ties: cheaper). Ranked on the same loci it is "
        "evaluated on, so small differences between sets are within sampling noise. The last row is the design "
        "make_tr_expansion_benchmark_alleles.py now uses.\n\n" + markdown_table(
            ["Targets kept", "Set", "Loci found", "Loci found (full run)", "Share of all found loci", "Cost",
             "Disk (GB gzipped)"], best_rows))
    plot_target_sets(points, best_by_size, os.path.join(HERE, "population_based_design_target_sets.png"))

    # 4. By distance to the nearest splice site
    distances = load_distances(locus["locus_id"] for locus in loci)
    loci_by_bin = collections.defaultdict(list)
    for locus in loci:
        if locus["locus_id"] in distances:
            loci_by_bin[distance_bin(distances[locus["locus_id"]])].append(locus)
    plot_fraction_by_distance(loci_by_bin, targets, os.path.join(HERE, "population_based_design_fraction_by_distance.png"))
    best4 = best_by_size[3][2]
    n_full_by_bin = count_polymorphic_loci_by_distance_bin()
    cumulative, totals = [], collections.Counter()
    for b in DISTANCE_BINS:
        sampled = loci_by_bin.get(b, [])
        if sampled:
            for key, target_set in (("", all_targets), ("_best4", best4)):
                found, seconds, _ = evaluate(sampled, target_set)
                totals["found" + key] += n_full_by_bin[b] * found / len(sampled)
                # dollars_per_second turns the sample's seconds into full-run dollars; per locus, that is
                # its seconds x dollars_per_second x (sampled loci / full-run loci)
                totals["dollars" + key] += n_full_by_bin[b] * seconds / len(sampled) * dollars_per_second / scale
        totals["loci"] += n_full_by_bin[b]
        cumulative.append({"threshold_label": f"{b[1]:,}" if b[1] < 10 ** 9 else "any", "loci": totals["loci"], **totals})
    plot_cumulative_by_distance(cumulative, os.path.join(HERE, "population_based_design_cumulative_by_distance.png"))
    sections.append(
        "## By distance to the nearest splice site\n\n![](population_based_design_fraction_by_distance.png)\n\n" + markdown_table(
            ["Distance (bp)", "Sampled loci", "Loci with a change >= 0.2", ">= 0.5", "Polymorphic loci in the full run"],
            [[bin_label(b), len(loci_by_bin[b]),
              percent(sum(is_found(l, all_targets, 0.2) for l in loci_by_bin[b]) / len(loci_by_bin[b])),
              percent(sum(is_found(l, all_targets, 0.5) for l in loci_by_bin[b]) / len(loci_by_bin[b])),
              f"{n_full_by_bin[b]:,}"] for b in DISTANCE_BINS if loci_by_bin.get(b)]) +
        "\n\n![](population_based_design_cumulative_by_distance.png)\n\n" + markdown_table(
            ["Distance limit (bp)", "Polymorphic loci", "Found, every target", "Cost, every target",
             f"Found, best 4 ({', '.join(short_label(t) for t in targets if t in best4)})", "Cost, best 4"],
            [[r["threshold_label"], f"{r['loci']:,}", f"{r['found']:,.0f}", f"${r['dollars']:,.0f}",
              f"{r['found_best4']:,.0f}", f"${r['dollars_best4']:,.0f}"] for r in cumulative]))
    sections.insert(3, "![](population_based_design_target_sets.png)")

    body = markdown.markdown(
        "# Population-based allele design: which target sizes are worth keeping\n\n"
        "Random polymorphic GENCODE v50 basic loci (at least 2 distinct HPRC256 allele sizes), scored on Modal "
        "L4s with SpliceAI + SAI-10k (production code path, distance 10,000). Targets: the HPRC256 "
        f"{', '.join(short_target_label(t) for t in targets if is_observed_length(t))} percentile total lengths, and "
        f"the 99.5th plus {', '.join(short_target_label(t).lstrip('+') for t in targets if not is_observed_length(t))} "
        "the motif's range (75th percentile, "
        "across loci with that motif and 5+ distinct lengths, of the 97.5th-2.5th percentile range). A locus is found "
        "when a kept allele's largest delta score on the selected transcript is >= 0.2.\n\n" + "\n\n".join(sections),
        extensions=["tables"])
    with open(REPORT_HTML, "w") as f:
        f.write("<!doctype html><html><head><meta charset='utf-8'><title>Population-based design</title><style>"
                "body{font:15px/1.5 -apple-system,Segoe UI,Roboto,sans-serif;max-width:1100px;margin:24px auto;padding:0 16px;"
                "color:#1d1d1f;background:#fff}table{border-collapse:collapse;margin:12px 0;font-size:.9em}"
                "th,td{border-bottom:1px solid #dde1e6;padding:4px 9px;text-align:right}th{background:#f6f7f9}"
                "td:first-child,th:first-child{text-align:left}img{max-width:100%}</style></head><body>" + body + "</body></html>")
    print("\n\n".join(sections))
    print(f"\nWrote {REPORT_HTML}")


if __name__ == "__main__":
    main()
