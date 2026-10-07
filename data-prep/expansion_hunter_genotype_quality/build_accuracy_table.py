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
    EH_allele_quality_largest_10, mean_abs_error_of_longest_truth_alleles
                          the per-locus accuracy score: across all scored genomes, the
                          N_LONGEST_TRUTH_ALLELES longest truth alleles at the locus and how well
                          ExpansionHunter called each of them. Long alleles are where a short-read
                          caller is most likely to fail, so this says how far its calls at this locus
                          can be trusted when an allele is long. Ties in truth length go to the sample
                          that sorts first by sample label. An allele is called correctly when
                          |ExpansionHunter - truth| is at most 1 repeat or 10% of the truth size,
                          whichever is larger, and ExpansionHunter did not call the reference size for
                          a truth allele that differs from it (a missed variant is an error even when
                          it is off by only one repeat). The score and its mean error are left blank
                          when fewer than N_LONGEST_TRUTH_ALLELES alleles were scored, so every score
                          comes from the same number of alleles; this mostly happens where the
                          assemblies give usable truth in only one or two genomes.
    EH_allele_quality_largest_10_truth_range
                          the truth sizes those alleles span, in repeats ("33-36"), blank with the score
    EH_allele_quality_all_non_ref, n_non_reference_truth_alleles
                          the same rule applied to every scored allele whose truth size differs from the
                          reference (longer or shorter): the fraction called correctly, and how many there
                          were. Blank where there were none.
    EH_allele_quality_<N>_haplotypes_distribution
                          every scored allele as distinct (truth, ExpansionHunter) repeat-count pairs with
                          their counts, most common first: "10,10x260;11,10x4". N is twice the number of
                          scored genomes; a haploid call contributes one allele.

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
# The per-locus score also accepts an error of this many repeats, so short alleles are not held to an
# exact match while long ones get 10%.
WITHIN_REPEATS_TOLERANCE = 1

ARRAY_NAMES = ("has_call", "pok_sum", "pok_sum_compared", "has_truth", "is_exact", "is_within1",
               "abs_error", "n_alleles_compared", "n_alleles_within_10_percent")
FLOAT_ARRAYS = ("pok_sum", "pok_sum_compared", "abs_error")
PER_ALLELE_ARRAYS = ("truth_short_allele", "truth_long_allele", "eh_short_allele", "eh_long_allele")
LARGEST_10_SCORE_COLUMN = "EH_allele_quality_largest_10"
LARGEST_10_TRUTH_RANGE_COLUMN = "EH_allele_quality_largest_10_truth_range"
NON_REFERENCE_SCORE_COLUMN = "EH_allele_quality_all_non_ref"
NON_REFERENCE_COUNTS = ("n_non_reference_truth_alleles", "n_non_reference_truth_alleles_called_correctly")
# Bits for each repeat count in an allele pair key (see allele_pair_keys); counts must be below 65,536.
PAIR_VALUE_BITS = 16


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
    """Return (starts, ends, motif sizes) as arrays in file order, which is the array order.

    The text columns are not kept: 5.7 million locus ids and motifs as Python strings take ~3 GB, so
    read_catalog_text_columns reads them again while the table is written.
    """
    starts, ends, motif_sizes = [], [], []
    with gzip.open(bed_path, "rt") as f:
        for line in f:
            _, start, end, motif = line.split("\t")[:4]
            starts.append(int(start))
            ends.append(int(end))
            motif_sizes.append(len(motif))
    return np.array(starts), np.array(ends), np.array(motif_sizes)


def read_catalog_text_columns(bed_path):
    """Yield (locus_id, chrom, start, end, motif) for each catalog row, in file order."""
    with gzip.open(bed_path, "rt") as f:
        for line in f:
            chrom, start, end, motif = line.split("\t")[:4]
            yield f"{chrom}-{start}-{end}-{motif}", chrom, start, end, motif


def keep_longest_truth_alleles(longest_truth, longest_eh, truth_alleles, eh_alleles, loci_per_chunk=500_000):
    """Merge one sample's alleles into the running per-locus N_LONGEST_TRUTH_ALLELES longest.

    All arrays are (rows, n_loci) with nan for an empty slot. Returns the updated (truth, eh) pair. The
    sort is stable and the running set comes first, so a tie in truth length keeps the allele seen
    earlier. Loci are sorted loci_per_chunk at a time: sorting all 5.7 million at once needs ~3.5 GB of
    temporary arrays per sample.
    """
    truth = np.vstack([longest_truth, truth_alleles])
    eh = np.vstack([longest_eh, eh_alleles])
    n_kept = min(len(truth), N_LONGEST_TRUTH_ALLELES)
    kept_truth = np.empty((n_kept, truth.shape[1]), dtype=truth.dtype)
    kept_eh = np.empty((n_kept, eh.shape[1]), dtype=eh.dtype)
    for start in range(0, truth.shape[1], loci_per_chunk):
        chunk = slice(start, start + loci_per_chunk)
        order = np.argsort(np.where(np.isnan(truth[:, chunk]), np.inf, -truth[:, chunk]), axis=0,
                           kind="stable")[:n_kept]
        kept_truth[:, chunk] = np.take_along_axis(truth[:, chunk], order, axis=0)
        kept_eh[:, chunk] = np.take_along_axis(eh[:, chunk], order, axis=0)
    return kept_truth, kept_eh


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


def is_called_correctly(truth, eh, reference_copies):
    """Return which alleles ExpansionHunter called correctly: |EH - truth| is at most
    WITHIN_REPEATS_TOLERANCE repeats or WITHIN_PERCENT_TOLERANCE % of truth, whichever is larger, and EH
    did not call the reference size for a truth allele that differs from it. False where truth is nan."""
    error = np.abs(eh - truth)
    within_tolerance = (error <= WITHIN_REPEATS_TOLERANCE) | (100 * error <= WITHIN_PERCENT_TOLERANCE * truth)
    missed_variant = (eh == reference_copies) & (truth != reference_copies)
    return within_tolerance & ~missed_variant


def allele_pair_keys(truth, eh):
    """Encode each scored (truth, EH) allele pair of one sample as one int64: locus row, truth, EH.

    truth and eh are (rows, n_loci) with nan where nothing was scored.
    """
    present = ~np.isnan(truth)
    locus_rows = np.nonzero(present)[1].astype(np.int64)
    truth_values = truth[present].astype(np.int64)
    eh_values = eh[present].astype(np.int64)
    for name, values in (("truth", truth_values), ("ExpansionHunter", eh_values)):
        if len(values) and (values.min() < 0 or values.max() >= 1 << PAIR_VALUE_BITS):
            raise ValueError(f"a {name} repeat count is outside 0 to {(1 << PAIR_VALUE_BITS) - 1}")
    return (locus_rows << (2 * PAIR_VALUE_BITS)) | (truth_values << PAIR_VALUE_BITS) | eh_values


def merge_allele_pair_counts(keys, counts, new_keys):
    """Add new_keys (one per allele) to the running sorted distinct keys and their counts."""
    unique_keys, inverse = np.unique(np.concatenate([keys, new_keys]), return_inverse=True)
    counts = np.bincount(inverse, weights=np.concatenate([counts, np.ones(len(new_keys))]))
    return unique_keys, counts.astype(np.int64)


def allele_pair_distribution_by_locus(keys, counts, n_loci):
    """Return a function giving one locus row's distribution string, e.g. "10,10x260;11,10x4": each
    distinct (truth, ExpansionHunter) repeat-count pair and how many alleles had it, most common first,
    ties by truth then ExpansionHunter repeat count."""
    locus_rows = keys >> (2 * PAIR_VALUE_BITS)
    order = np.lexsort((keys, -counts, locus_rows))
    truth = (keys[order] >> PAIR_VALUE_BITS) & ((1 << PAIR_VALUE_BITS) - 1)
    eh = keys[order] & ((1 << PAIR_VALUE_BITS) - 1)
    counts = counts[order]
    boundaries = np.searchsorted(locus_rows[order], np.arange(n_loci + 1))

    def distribution(row):
        start, end = boundaries[row], boundaries[row + 1]
        return ";".join(f"{t},{e}x{n}" for t, e, n in zip(truth[start:end], eh[start:end], counts[start:end]))
    return distribution


def sum_comparisons(comparison_paths, reference_copies):
    """Add up the per-sample arrays. Returns ({array name: per-locus total}, longest truth alleles,
    the ExpansionHunter calls paired with them, the distinct allele pair keys and their counts).

    The totals also hold n_non_reference_truth_alleles and n_non_reference_truth_alleles_called_correctly.
    """
    n_loci = len(reference_copies)
    totals = {name: np.zeros(n_loci, dtype=np.float64 if name in FLOAT_ARRAYS else np.int64)
              for name in ARRAY_NAMES + NON_REFERENCE_COUNTS}
    longest_truth = np.full((0, n_loci), np.nan, dtype=np.float32)
    longest_eh = np.full((0, n_loci), np.nan, dtype=np.float32)
    pair_keys, pair_counts = np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.int64)
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
            truth, eh = sample_alleles(arrays)
        longest_truth, longest_eh = keep_longest_truth_alleles(longest_truth, longest_eh, truth, eh)
        non_reference = ~np.isnan(truth) & (truth != reference_copies)
        totals["n_non_reference_truth_alleles"] += non_reference.sum(axis=0)
        totals["n_non_reference_truth_alleles_called_correctly"] += (
            non_reference & is_called_correctly(truth, eh, reference_copies)).sum(axis=0)
        pair_keys, pair_counts = merge_allele_pair_counts(pair_keys, pair_counts, allele_pair_keys(truth, eh))
        print(f"  added {os.path.basename(path)}")
    return totals, longest_truth, longest_eh, pair_keys, pair_counts


def score_longest_truth_alleles(longest_truth, longest_eh, reference_copies):
    """Return per-locus (n scored, shortest of them, longest of them, fraction called correctly,
    mean |EH - truth|).

    reference_copies holds each locus's reference repeat count. The last two are nan unless all
    N_LONGEST_TRUTH_ALLELES were scored.
    """
    present = ~np.isnan(longest_truth)
    n_scored = present.sum(axis=0)
    error = np.abs(longest_eh - longest_truth)
    correct = present & is_called_correctly(longest_truth, longest_eh, reference_copies)
    shortest = np.where(n_scored > 0, np.nanmin(np.where(present, longest_truth, np.inf), axis=0), np.nan)
    longest = np.where(n_scored > 0, np.nanmax(np.where(present, longest_truth, -np.inf), axis=0), np.nan)
    n_scored_if_complete = np.where(n_scored == N_LONGEST_TRUTH_ALLELES, n_scored, 0)
    return (n_scored, shortest, longest, divide(correct.sum(axis=0), n_scored_if_complete),
            divide(np.where(present, error, 0).sum(axis=0), n_scored_if_complete))


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

    starts, ends, motif_sizes = read_catalog(args.catalog_bed)
    n_loci = len(starts)
    print(f"{n_loci:,} locus definitions in {args.catalog_bed}")
    # Rounded down, as the truth genotyper's NumRepeatsInReference is (identical at every HG002 locus).
    reference_copies = ((ends - starts) // motif_sizes).astype(np.float32)

    comparison_paths = sorted(glob.glob(os.path.join(args.comparison_dir, "**", "*.comparison.npz"),
                                        recursive=True))
    if not comparison_paths:
        parser.error(f"no *.comparison.npz files under {args.comparison_dir}. Run "
                     f"run_comparison_on_selected_samples.py first, then download its output.")
    comparison_paths = select_cohort_files(comparison_paths, *read_cohort(args.sample_table_path),
                                           args.allow_partial_cohort, parser)

    print(f"summing {len(comparison_paths)} samples:")
    totals, longest_truth, longest_eh, pair_keys, pair_counts = sum_comparisons(comparison_paths, reference_copies)

    mean_pok = divide(totals["pok_sum"], totals["has_call"])
    mean_pok_compared = divide(totals["pok_sum_compared"], totals["has_truth"])
    fraction_exact = divide(totals["is_exact"], totals["has_truth"])
    fraction_within1 = divide(totals["is_within1"], totals["has_truth"])
    mean_abs_error = divide(totals["abs_error"], totals["has_truth"])
    fraction_within_10_percent = divide(totals["n_alleles_within_10_percent"], totals["n_alleles_compared"])
    (n_longest_scored, shortest_longest, longest_longest, fraction_longest_called_correctly,
     mean_abs_error_longest) = score_longest_truth_alleles(longest_truth, longest_eh, reference_copies)
    fraction_non_reference_called_correctly = divide(totals["n_non_reference_truth_alleles_called_correctly"],
                                                     totals["n_non_reference_truth_alleles"])
    distribution = allele_pair_distribution_by_locus(pair_keys, pair_counts, n_loci)

    os.makedirs(os.path.dirname(args.by_locus_path), exist_ok=True)
    columns = ["locus_id", "chrom", "start_0based", "end_1based", "motif", "motif_size",
               "n_samples_with_call", "mean_pOk", "n_compared", "n_exact", "n_within1",
               "fraction_exact", "fraction_within1", "mean_abs_error", "mean_pOk_on_compared",
               "n_alleles_compared", "n_alleles_within_10_percent", "fraction_alleles_within_10_percent",
               "n_longest_truth_alleles_scored", "shortest_of_the_longest_truth_alleles",
               LARGEST_10_SCORE_COLUMN, LARGEST_10_TRUTH_RANGE_COLUMN, "mean_abs_error_of_longest_truth_alleles",
               "n_non_reference_truth_alleles", NON_REFERENCE_SCORE_COLUMN,
               f"EH_allele_quality_{2 * len(comparison_paths)}_haplotypes_distribution"]
    n_rows = 0
    with gzip.open(args.by_locus_path, "wt") as f:
        f.write("\t".join(columns) + "\n")
        for row, (locus_id, chrom, start, end, motif) in enumerate(read_catalog_text_columns(args.catalog_bed)):
            if totals["has_call"][row] == 0:
                continue  # never genotyped in any sample, so there is nothing to report
            has_largest_10_score = not np.isnan(fraction_longest_called_correctly[row])
            f.write("\t".join(str(v) for v in [
                locus_id, chrom, start, end, motif, len(motif),
                totals["has_call"][row], format_float(mean_pok[row]),
                totals["has_truth"][row], totals["is_exact"][row], totals["is_within1"][row],
                format_float(fraction_exact[row]), format_float(fraction_within1[row]),
                format_float(mean_abs_error[row]), format_float(mean_pok_compared[row]),
                totals["n_alleles_compared"][row], totals["n_alleles_within_10_percent"][row],
                format_float(fraction_within_10_percent[row]),
                n_longest_scored[row], format_float(shortest_longest[row], digits=0),
                format_float(fraction_longest_called_correctly[row]),
                f"{shortest_longest[row]:.0f}-{longest_longest[row]:.0f}" if has_largest_10_score else "",
                format_float(mean_abs_error_longest[row]),
                totals["n_non_reference_truth_alleles"][row],
                format_float(fraction_non_reference_called_correctly[row]),
                distribution(row),
            ]) + "\n")
            n_rows += 1
    if row != n_loci - 1:
        raise ValueError(f"{args.catalog_bed} changed while the table was being written")
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
        has_score = ~np.isnan(fraction_longest_called_correctly)
        if has_score.any():
            print(f"  {N_LONGEST_TRUTH_ALLELES}-longest-truth-allele score: mean "
                  f"{fraction_longest_called_correctly[has_score].mean():.4f} over {has_score.sum():,} loci, "
                  f"{(fraction_longest_called_correctly[has_score] == 1).mean():.1%} of them 1.0")
        has_non_reference = ~np.isnan(fraction_non_reference_called_correctly)
        if has_non_reference.any():
            print(f"  non-reference truth alleles called correctly: "
                  f"{totals['n_non_reference_truth_alleles_called_correctly'].sum() / totals['n_non_reference_truth_alleles'].sum():.2%} "
                  f"of {totals['n_non_reference_truth_alleles'].sum():,}, at {has_non_reference.sum():,} loci")


if __name__ == "__main__":
    main()
