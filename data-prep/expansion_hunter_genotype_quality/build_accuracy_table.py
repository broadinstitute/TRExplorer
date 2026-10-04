"""Sum the per-sample comparison arrays into one table with a row per locus definition.

Reads the .npz files that compare_eh_to_truth.py wrote, one per sample, plus the catalog BED they are
aligned to.

Each genome counts once. The sample table can list one genome more than once (HG002 at 10x, 20x and
full depth share sample_id HG002), and counting all of those would weight that genome three times, so
only the sample whose sample_label equals its sample_id (the full-depth one) is scored; the others are
skipped and named in the output.

Per locus, over the samples that had something to contribute:
    n_samples_with_call   samples where ExpansionHunter produced a genotype
    mean_pOk              its per-allele confidence, averaged over those calls
    n_compared            samples where the genotype could be scored against truth
    fraction_exact        of those, how often every allele matched truth
    fraction_within1      how often every allele was within one repeat
    mean_abs_error        mean |ExpansionHunter - truth| in repeat copies
    mean_pOk_on_compared  pOk over exactly the scored samples, so the rate and the confidence for one
                          locus come from the same observations
    n_alleles_compared, n_alleles_within_10_percent, fraction_alleles_within_10_percent
                          the same samples counted by allele: how many alleles were scored, and how many
                          of those were within 10% of the truth allele size
    n_longest_truth_alleles_scored, shortest_of_the_longest_truth_alleles,
    fraction_of_longest_truth_alleles_within_10_percent, mean_abs_error_of_longest_truth_alleles
                          the per-locus accuracy score: across all scored genomes, the
                          N_LONGEST_TRUTH_ALLELES longest truth alleles at the locus (fewer if fewer
                          were scored) and how well ExpansionHunter called each of them. Long alleles
                          are where a short-read caller is most likely to fail, so this says how far its
                          calls at this locus can be trusted when an allele is long. Ties in truth
                          length go to the sample that sorts first by sample label.

Rates are over whatever samples had a call and a truth genotype at that locus, which differs from
locus to locus, so n_compared is carried alongside them to weight or filter by.

Usage:
    python3 build_accuracy_table.py --catalog-bed TRExplorer.repeat_catalog_v2.1.hg38.1_to_1000bp_motifs.bed.gz \\
        --comparison-dir data/comparisons_v2.1
"""

import argparse
import glob
import gzip
import os

import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_COMPARISON_DIR = os.path.join(SCRIPT_DIR, "data", "comparisons")
DEFAULT_BY_LOCUS_PATH = os.path.join(SCRIPT_DIR, "data", "accuracy_by_locus.tsv.gz")
DEFAULT_SAMPLE_TABLE_PATH = os.path.join(SCRIPT_DIR, "short_read_samples_with_truth_data.tsv")

N_LONGEST_TRUTH_ALLELES = 10
# Same rule as compare_eh_to_truth.WITHIN_PERCENT_TOLERANCE: |ExpansionHunter - truth| <= 10% of truth.
WITHIN_PERCENT_TOLERANCE = 10

ARRAY_NAMES = ("has_call", "pok_sum", "pok_sum_compared", "has_truth", "is_exact", "is_within1",
               "abs_error", "n_alleles_compared", "n_alleles_within_10_percent")
FLOAT_ARRAYS = ("pok_sum", "pok_sum_compared", "abs_error")
PER_ALLELE_ARRAYS = ("truth_short_allele", "truth_long_allele", "eh_short_allele", "eh_long_allele")


def read_cohort(sample_table_path):
    """Return (sample labels to score, sample labels to skip): one sample label per genome.

    A table without a sample_label column uses the sample_id. When several labels share a sample_id,
    the one equal to the sample_id is scored and the rest are skipped.
    """
    with open(sample_table_path) as f:
        column = {name: i for i, name in enumerate(next(f).rstrip("\n").split("\t"))}
        rows = [line.rstrip("\n").split("\t") for line in f if line.strip()]
    label_column = column.get("sample_label", column["sample_id"])
    labels_by_sample_id = {}
    for row in rows:
        labels_by_sample_id.setdefault(row[column["sample_id"]], []).append(row[label_column])
    scored, skipped = [], []
    for sample_id, labels in labels_by_sample_id.items():
        if len(labels) == 1:
            scored.append(labels[0])
            continue
        if sample_id not in labels:
            raise ValueError(f"{sample_table_path} lists {sample_id} as {labels}, none of them labeled "
                             f"{sample_id}, so there is no way to pick which one represents the genome")
        scored.append(sample_id)
        skipped.extend(label for label in labels if label != sample_id)
    return sorted(scored), sorted(skipped)


def read_catalog(bed_path):
    """Return (locus_ids, chroms, starts, ends, motifs) in file order, which is the array order."""
    locus_ids, chroms, starts, ends, motifs = [], [], [], [], []
    with gzip.open(bed_path, "rt") as f:
        for line in f:
            chrom, start, end, motif = line.split("\t")[:4]
            locus_ids.append(f"{chrom}-{start}-{end}-{motif}")
            chroms.append(chrom)
            starts.append(int(start))
            ends.append(int(end))
            motifs.append(motif)
    return locus_ids, chroms, starts, ends, motifs


def keep_longest_truth_alleles(longest_truth, longest_eh, truth_alleles, eh_alleles):
    """Merge one sample's alleles into the running per-locus N_LONGEST_TRUTH_ALLELES longest.

    All arrays are (rows, n_loci) with nan for an empty slot. Returns the updated (truth, eh) pair. The
    sort is stable and the running set comes first, so a tie in truth length keeps the allele seen
    earlier.
    """
    truth = np.vstack([longest_truth, truth_alleles])
    eh = np.vstack([longest_eh, eh_alleles])
    order = np.argsort(np.where(np.isnan(truth), np.inf, -truth), axis=0, kind="stable")[
        :N_LONGEST_TRUTH_ALLELES]
    return np.take_along_axis(truth, order, axis=0), np.take_along_axis(eh, order, axis=0)


def sample_alleles(arrays):
    """Return (truth, eh) arrays of shape (2, n_loci) for one sample's scored alleles, nan elsewhere.

    A haploid call is stored in both the short and the long slot, so its long slot is dropped to count
    it once.
    """
    scored = arrays["has_truth"] == 1
    haploid = arrays["n_alleles_compared"] == 1
    truth = np.vstack([np.where(scored, arrays["truth_short_allele"], np.nan),
                       np.where(scored & ~haploid, arrays["truth_long_allele"], np.nan)])
    eh = np.vstack([np.where(scored, arrays["eh_short_allele"], np.nan),
                    np.where(scored & ~haploid, arrays["eh_long_allele"], np.nan)])
    return truth.astype(np.float32), eh.astype(np.float32)


def sum_comparisons(comparison_paths, n_loci):
    """Add up the per-sample arrays. Returns ({array name: per-locus total}, longest truth alleles,
    the ExpansionHunter calls paired with them)."""
    totals = {name: np.zeros(n_loci, dtype=np.float64 if name in FLOAT_ARRAYS else np.int64)
              for name in ARRAY_NAMES}
    longest_truth = np.full((0, n_loci), np.nan, dtype=np.float32)
    longest_eh = np.full((0, n_loci), np.nan, dtype=np.float32)
    for path in comparison_paths:
        with np.load(path) as arrays:
            missing = [name for name in ARRAY_NAMES + PER_ALLELE_ARRAYS if name not in arrays.files]
            if missing:
                raise ValueError(f"{path} lacks {missing}: it was written by an older "
                                 f"compare_eh_to_truth.py. Rerun run_comparison_on_selected_samples.py "
                                 f"with --force-comparison-step.")
            if len(arrays["has_call"]) != n_loci:
                raise ValueError(f"{path} has {len(arrays['has_call']):,} loci but the catalog has "
                                 f"{n_loci:,}, so they were not built from the same catalog")
            for name in ARRAY_NAMES:
                totals[name] += arrays[name]
            longest_truth, longest_eh = keep_longest_truth_alleles(
                longest_truth, longest_eh, *sample_alleles(arrays))
        print(f"  added {os.path.basename(path)}")
    return totals, longest_truth, longest_eh


def score_longest_truth_alleles(longest_truth, longest_eh):
    """Return per-locus (n scored, shortest of them, fraction within 10%, mean |EH - truth|)."""
    present = ~np.isnan(longest_truth)
    n_scored = present.sum(axis=0)
    error = np.abs(longest_eh - longest_truth)
    within = present & (100 * error <= WITHIN_PERCENT_TOLERANCE * longest_truth)
    shortest = np.where(n_scored > 0, np.nanmin(np.where(present, longest_truth, np.inf), axis=0), np.nan)
    return (n_scored, shortest, divide(within.sum(axis=0), n_scored),
            divide(np.where(present, error, 0).sum(axis=0), n_scored))


def divide(numerator, denominator):
    """Element-wise mean, leaving nan where nothing was counted."""
    return np.where(denominator > 0, numerator / np.maximum(denominator, 1), np.nan)


def format_float(value, digits=4):
    return "" if np.isnan(value) else f"{value:.{digits}f}"


def select_cohort_files(comparison_paths, scored_labels, skipped_labels, allow_partial, parser):
    """Return the comparison files to sum: exactly one per scored sample label.

    A missing sample silently changes every denominator and a sample counted twice silently weights it
    double, so neither is allowed to pass for a complete result by accident. Files for skipped labels
    (another depth of a genome already scored) are left out.
    """
    found = [os.path.basename(path).split(".comparison.npz")[0] for path in comparison_paths]
    duplicates = sorted({label for label in found if found.count(label) > 1})
    if duplicates:
        parser.error(f"more than one comparison file for {duplicates}. Each sample must contribute "
                     f"exactly once, or its genotypes would be counted twice.")
    unexpected = sorted(set(found) - set(scored_labels) - set(skipped_labels))
    if unexpected:
        parser.error(f"comparison files for samples that are not in the sample table: {unexpected}")
    missing = sorted(set(scored_labels) - set(found))
    if missing and not allow_partial:
        parser.error(f"{len(missing)} of {len(scored_labels)} samples have no comparison file "
                     f"({', '.join(missing[:5])}{', ...' if len(missing) > 5 else ''}). Wait for them, "
                     f"or pass --allow-partial-cohort to score the cohort that is here.")
    if missing:
        print(f"WARNING: scoring {len(scored_labels) - len(missing)} of {len(scored_labels)} samples; "
              f"missing {', '.join(missing)}")
    left_out = sorted(set(found) & set(skipped_labels))
    if left_out:
        print(f"left out {', '.join(left_out)}: another sample of the same genome is scored")
    return sorted(path for path, label in zip(comparison_paths, found) if label in scored_labels)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog-bed", required=True, help="the catalog BED the comparisons are aligned to")
    parser.add_argument("--comparison-dir", default=DEFAULT_COMPARISON_DIR,
                        help="directory holding the per-sample *.comparison.npz files")
    parser.add_argument("--by-locus-path", default=DEFAULT_BY_LOCUS_PATH)
    parser.add_argument("--sample-table-path", default=DEFAULT_SAMPLE_TABLE_PATH,
                        help="the cohort this table is supposed to cover")
    parser.add_argument("--allow-partial-cohort", action="store_true",
                        help="build the table even though some samples are missing. Rates then come "
                             "from a different cohort than the sample table describes, so this has to "
                             "be asked for.")
    args = parser.parse_args()

    locus_ids, chroms, starts, ends, motifs = read_catalog(args.catalog_bed)
    print(f"{len(locus_ids):,} locus definitions in {args.catalog_bed}")

    comparison_paths = sorted(glob.glob(os.path.join(args.comparison_dir, "**", "*.comparison.npz"),
                                        recursive=True))
    if not comparison_paths:
        parser.error(f"no *.comparison.npz files under {args.comparison_dir}. Run "
                     f"run_comparison_on_selected_samples.py first, then download its output.")
    comparison_paths = select_cohort_files(comparison_paths, *read_cohort(args.sample_table_path),
                                           args.allow_partial_cohort, parser)

    print(f"summing {len(comparison_paths)} samples:")
    totals, longest_truth, longest_eh = sum_comparisons(comparison_paths, len(locus_ids))

    mean_pok = divide(totals["pok_sum"], totals["has_call"])
    mean_pok_compared = divide(totals["pok_sum_compared"], totals["has_truth"])
    fraction_exact = divide(totals["is_exact"], totals["has_truth"])
    fraction_within1 = divide(totals["is_within1"], totals["has_truth"])
    mean_abs_error = divide(totals["abs_error"], totals["has_truth"])
    fraction_within_10_percent = divide(totals["n_alleles_within_10_percent"], totals["n_alleles_compared"])
    (n_longest_scored, shortest_longest, fraction_longest_within_10_percent,
     mean_abs_error_longest) = score_longest_truth_alleles(longest_truth, longest_eh)

    os.makedirs(os.path.dirname(args.by_locus_path), exist_ok=True)
    columns = ["locus_id", "chrom", "start_0based", "end_1based", "motif", "motif_size",
               "n_samples_with_call", "mean_pOk", "n_compared", "n_exact", "n_within1",
               "fraction_exact", "fraction_within1", "mean_abs_error", "mean_pOk_on_compared",
               "n_alleles_compared", "n_alleles_within_10_percent", "fraction_alleles_within_10_percent",
               "n_longest_truth_alleles_scored", "shortest_of_the_longest_truth_alleles",
               "fraction_of_longest_truth_alleles_within_10_percent", "mean_abs_error_of_longest_truth_alleles"]
    n_rows = 0
    with gzip.open(args.by_locus_path, "wt") as f:
        f.write("\t".join(columns) + "\n")
        for row in range(len(locus_ids)):
            if totals["has_call"][row] == 0:
                continue  # never genotyped in any sample, so there is nothing to report
            f.write("\t".join(str(v) for v in [
                locus_ids[row], chroms[row], starts[row], ends[row], motifs[row], len(motifs[row]),
                totals["has_call"][row], format_float(mean_pok[row]),
                totals["has_truth"][row], totals["is_exact"][row], totals["is_within1"][row],
                format_float(fraction_exact[row]), format_float(fraction_within1[row]),
                format_float(mean_abs_error[row]), format_float(mean_pok_compared[row]),
                totals["n_alleles_compared"][row], totals["n_alleles_within_10_percent"][row],
                format_float(fraction_within_10_percent[row]),
                n_longest_scored[row], format_float(shortest_longest[row], digits=0),
                format_float(fraction_longest_within_10_percent[row]),
                format_float(mean_abs_error_longest[row]),
            ]) + "\n")
            n_rows += 1
    print(f"wrote {n_rows:,} rows to {args.by_locus_path}")

    n_compared = totals["has_truth"].sum()
    if n_compared:
        print(f"\n{n_compared:,} scored (locus, sample) observations")
        print(f"  exact match:       {100 * totals['is_exact'].sum() / n_compared:.2f}%")
        print(f"  within 1 copy:     {100 * totals['is_within1'].sum() / n_compared:.2f}%")
        print(f"  mean |EH - truth|: {totals['abs_error'].sum() / n_compared:.4f} repeats")
        print(f"  mean pOk on those: {totals['pok_sum_compared'].sum() / n_compared:.4f}")
        n_alleles = totals["n_alleles_compared"].sum()
        print(f"  alleles within 10% of truth: {100 * totals['n_alleles_within_10_percent'].sum() / n_alleles:.2f}% "
              f"of {n_alleles:,}")


if __name__ == "__main__":
    main()
