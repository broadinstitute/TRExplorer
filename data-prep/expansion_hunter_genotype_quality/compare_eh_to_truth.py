"""Compare one sample's ExpansionHunter genotypes to its assembly-derived truth genotypes, and write
the per-locus result as arrays aligned to the catalog's own row order.

Runs on one sample. run_comparison_on_selected_samples.py runs it once per sample on Hail Batch, and
build_accuracy_table.py sums the resulting arrays across samples.

Why arrays rather than a TSV: the catalog holds 5,657,854 locus definitions, so a per-sample TSV would
be ~5.7 hundred million rows across 100 samples, and re-parsing that locally costs hours. Positions in
these arrays are catalog row numbers, so summing 100 samples is 100 numpy adds.

What a locus contributes:
  has_call         ExpansionHunter produced a genotype
  pok_sum          that genotype's per-allele pOk, averaged over alleles, where has_call
  ptooshort_sum    the same for pTooShort, and ptoolong_sum for pTooLong. pOk pools two opposite
                   failure modes, so the direction probabilities are kept separately: pTooShort is
                   the score for a missed expansion, pTooLong for a missed contraction.
  has_truth        the locus is COMPARABLE in this sample: ExpansionHunter called it, a truth
                   genotype exists, the locus lies wholly inside this sample's DipCall
                   high-confidence regions, and the two sides agree about ploidy. A locus with truth
                   but no call is not comparable and contributes nothing, the same convention the
                   earlier boundary-comparison work used.
  pok_sum_compared the same pOk, but only where has_truth, so a rate and a confidence can be read off
                   the same set of observations. ptooshort_sum_compared and ptoolong_sum_compared
                   are the matching restrictions of the other two.
  is_quick         ExpansionHunter genotyped this locus on its fast path ("QuickGenotype": true)
                   rather than through the full genotyping module, which is ~80x slower per locus
  is_exact         every allele matches truth exactly            (only where has_truth)
  is_within1       every allele is within 1 repeat of truth      (only where has_truth)
  n_alleles_compared           alleles scored: 2, or 1 for a haploid call (only where has_truth)
  n_alleles_within_10_percent  of those, alleles within 10% of the truth allele size, i.e.
                               |ExpansionHunter - truth| <= 0.1 x truth (only where has_truth). A locus
                               has every allele within 10% when the two counts are equal. Below 10
                               truth copies this only allows an exact match. Measuring in bp instead
                               gives the same result, since the motif size cancels.
  abs_error        mean |ExpansionHunter - truth| over alleles   (only where has_truth)
  truth_short_allele, truth_long_allele   the truth repeat counts, and
  eh_short_allele, eh_long_allele         the ExpansionHunter repeat counts they were paired with
                       (shorter with shorter), so a per-allele score, such as how well the longest truth
                       alleles were called, can be computed downstream. A haploid call has one allele,
                       stored in both. nan unless has_truth.
  eh_diff_from_ref     how far the CALL sits from a homozygous-reference genotype, in repeats:
                       max over alleles of |allele - reference copies|. Separates "it called an
                       expansion" from "it called the reference", which the rates alone conflate,
                       since most loci are reference length in most samples.
  truth_short_diff_from_ref  the TRUTH short allele minus the reference copies, and
  truth_long_diff_from_ref   the same for the long allele. Both signed, so a missed expansion (call
                       at the reference, truth above it) can be told apart both from a false one and
                       from a missed contraction, which is what pTooShort and pTooLong score
                       separately. max of the two absolute values is the distance eh_diff_from_ref
                       reports for the call side. Two arrays rather than one signed distance because
                       ~3% of the loci more than 2 repeats from the reference straddle it, with the
                       short allele contracted and the long one expanded; a single signed number
                       would have to name one direction and hide the other.
  coverage             locus coverage ExpansionHunter reported
  n_spanning_reads     reads spanning the whole repeat. An expansion longer than a read removes
                       these, so a reference-length call with few of them is the signature of an
                       expansion the caller did not see.
  n_flanking_reads     reads overlapping one boundary only
  n_irr_reads          in-repeat reads

The high-confidence check is not optional: the truth genotyper reports the reference length both for a
locus that genuinely matches the reference and for one DipCall could never call, and those two are
indistinguishable in its output. Scoring the second kind would credit ExpansionHunter for agreeing with
an invented genotype.

A haploid ExpansionHunter call on a sex chromosome is compared against truth only when truth is
homozygous; otherwise the two sides disagree about ploidy and the comparison is not meaningful.

Usage:
    python3 compare_eh_to_truth.py --catalog-bed catalog.bed.gz --eh-json sample.json.gz \\
        --truth-tsv sample.tandem_repeat_genotypes.tsv.gz --high-confidence-bed sample.dip.bed.gz \\
        --output-npz sample.comparison.npz
"""

import argparse
import gzip
import re

import numpy as np

# A genotype counts as near-concordant when every allele is within this many repeat copies of truth.
NEAR_CONCORDANCE_TOLERANCE = 1
# An allele counts as within 10% when its error is at most this percentage of the truth allele size.
# Kept as an integer percentage so the check stays exact integer arithmetic on repeat counts.
WITHIN_PERCENT_TOLERANCE = 10

LOCUS_RESULTS_PATTERN = re.compile(r'"LocusResults":\s*\{')
LOCUS_ID_PATTERN = re.compile(r'"LocusId":\s*"([^"]+)"')
QUICK_GENOTYPE_PATTERN = re.compile(r'"QuickGenotype":\s*(true|false)')
GENOTYPE_PATTERN = re.compile(r'"Genotype":\s*"([^"]+)"')
# Per-allele genotype-quality probabilities, averaged over a locus's alleles. Each name yields two
# output arrays: <name>_sum over every call, and <name>_sum_compared over the calls scored against
# truth. Each is also the record key prefix, so the parser and the arrays stay in step.
PROBABILITY_PATTERNS = {
    "pok": re.compile(r'"pOk":\s*([0-9.eE+-]+)'),
    "ptooshort": re.compile(r'"pTooShort":\s*([0-9.eE+-]+)'),
    "ptoolong": re.compile(r'"pTooLong":\s*([0-9.eE+-]+)'),
}
COVERAGE_PATTERN = re.compile(r'"Coverage":\s*([0-9.eE+-]+)')
READ_COUNT_PATTERNS = {
    "n_spanning_reads": re.compile(r'"CountsOfSpanningReads":\s*"([^"]*)"'),
    "n_flanking_reads": re.compile(r'"CountsOfFlankingReads":\s*"([^"]*)"'),
    "n_irr_reads": re.compile(r'"CountsOfInrepeatReads":\s*"([^"]*)"'),
}
READ_COUNT_PAIR_PATTERN = re.compile(r"\((\d+),\s*(\d+)\)")


def total_reads(counts):
    """Total reads in an ExpansionHunter count string like "(9, 20), (10, 1)" -> 21."""
    return sum(int(reads) for _, reads in READ_COUNT_PAIR_PATTERN.findall(counts or ""))


def open_maybe_gzipped(path, mode="rt"):
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode)


def locus_id_without_chr(locus_id):
    """Return the LocusId with any leading "chr" removed: "chr1-100-130-CAG" -> "1-100-130-CAG".

    The three sides do not agree on the prefix: the catalog BED and filter_vcf_to_tandem_repeats write
    "chr1-...", while the TRExplorer v2.1 ExpansionHunter catalog uses "1-...". Every LocusId is passed
    through this before matching, so either spelling joins.
    """
    return locus_id[3:] if locus_id.startswith("chr") else locus_id


def catalog_row_numbers(catalog_bed_path):
    """Return {LocusId without "chr": row number} for the catalog, in file order.

    The LocusId is rebuilt the same way convert_bed_to_expansion_hunter_catalog builds it, which is
    also what filter_vcf_to_tandem_repeats writes, so all three sides agree without a coordinate join
    once the "chr" prefix is dropped (see locus_id_without_chr).
    """
    row_number_by_locus_id = {}
    with open_maybe_gzipped(catalog_bed_path) as f:
        for row_number, line in enumerate(f):
            chrom, start, end, motif = line.split("\t")[:4]
            row_number_by_locus_id[locus_id_without_chr(f"{chrom}-{start}-{end}-{motif}")] = row_number
    return row_number_by_locus_id


def parse_eh_json(eh_json_path):
    """Yield (LocusId, record) per locus, where record holds everything this script reads.

    Parsed line by line rather than with json.load because this file holds one record per catalog
    locus and is hundreds of megabytes; only three fields per locus are needed.

    The repeat counts come from the "Genotype" string ("9/9"), not from AlleleQualityMetrics, because
    the latter is deduplicated by allele size and so lists a homozygote's allele only once. pOk is
    averaged over the AlleleQualityMetrics entries, which is what the earlier reports used.

    A locus record starts where its object opens inside "LocusResults", found by counting braces, not
    at its "LocusId" line: ExpansionHunter writes a locus's keys in alphabetical order, so "Coverage"
    comes before "LocusId", and starting the record at the LocusId would credit each locus's coverage
    to the locus before it. This relies on one locus ending and the next starting on different lines,
    which is how ExpansionHunter writes it (pretty-printed). No string value in the file contains a
    brace.
    """
    locus_results_depth = None  # brace depth just inside the LocusResults object, once it opens
    depth = 0
    locus_id, record = None, None
    with open_maybe_gzipped(eh_json_path) as f:
        for line in f:
            depth_before = depth
            if locus_results_depth is None and LOCUS_RESULTS_PATTERN.search(line):
                locus_results_depth = depth_before + 1
            depth += line.count("{") - line.count("}")
            if locus_results_depth is None:
                continue
            if depth_before <= locus_results_depth < depth:
                # This line opens an object directly inside LocusResults: the next locus starts.
                if locus_id is not None:
                    yield locus_id, record
                locus_id, record = None, new_eh_record()
            elif depth < locus_results_depth:
                # LocusResults closed; nothing after it (RunInfo) belongs to a locus.
                if locus_id is not None:
                    yield locus_id, record
                locus_id, record = None, None
            if record is None:
                continue
            match = LOCUS_ID_PATTERN.search(line)
            if match and locus_id is None:
                locus_id = match.group(1)
                # No `continue`: a writer that puts more than one field on a line would otherwise
                # have every field after the LocusId silently dropped.
            match = GENOTYPE_PATTERN.search(line)
            if match:
                record["genotype"] = match.group(1)
                continue
            match = QUICK_GENOTYPE_PATTERN.search(line)
            if match:
                record["is_quick"] = match.group(1) == "true"
                continue
            match = COVERAGE_PATTERN.search(line)
            if match:
                record["coverage"] = float(match.group(1))
                continue
            for name, pattern in READ_COUNT_PATTERNS.items():
                match = pattern.search(line)
                if match:
                    record[name] = total_reads(match.group(1))
                    break
            else:
                for name, pattern in PROBABILITY_PATTERNS.items():
                    for match in pattern.finditer(line):
                        record[f"{name}_values"].append(float(match.group(1)))
    if locus_id is not None:
        yield locus_id, record


def new_eh_record():
    record = {"genotype": None, "is_quick": False, "coverage": 0.0,
              "n_spanning_reads": 0, "n_flanking_reads": 0, "n_irr_reads": 0}
    record.update({f"{name}_values": [] for name in PROBABILITY_PATTERNS})
    return record


def parse_genotype(genotype):
    """Repeat count of each allele in an ExpansionHunter "Genotype" string, sorted.

    "9/9" -> (9, 9); "6/171" -> (6, 171); "4" -> (4,) for a haploid call on a sex chromosome.
    Returns None when the genotype is absent or null.
    """
    if not genotype or genotype == "null":
        return None
    return tuple(sorted(int(allele) for allele in genotype.split("/")))


def mean_or_none(values):
    return float(np.mean(values)) if values else None


def load_high_confidence_regions(bed_path):
    """Return {chrom: (starts, ends)} as sorted numpy arrays of the regions DipCall could call."""
    starts_and_ends = {}
    with open_maybe_gzipped(bed_path) as f:
        for line in f:
            chrom, start, end = line.split("\t")[:3]
            starts_and_ends.setdefault(chrom, ([], []))
            starts_and_ends[chrom][0].append(int(start))
            starts_and_ends[chrom][1].append(int(end))
    return {chrom: (np.array(starts), np.array(ends)) for chrom, (starts, ends) in starts_and_ends.items()}


def is_fully_callable(regions_by_chrom, chrom, start, end):
    """Is one interval wholly inside a single high-confidence region?

    DipCall's regions are disjoint, so it is enough to find the last region starting at or before the
    locus and check that it also covers the locus's end.
    """
    if chrom not in regions_by_chrom:
        return False
    starts, ends = regions_by_chrom[chrom]
    index = int(np.searchsorted(starts, start, side="right")) - 1
    return index >= 0 and end <= ends[index]


def compare_alleles(eh_counts, truth_short, truth_long):
    """Return (is_exact, is_within1, mean_absolute_error, n_alleles_within_10_percent), or None if the
    ploidies disagree.

    An allele is within 10% when |ExpansionHunter - truth| <= 10% of the truth allele size. Repeat counts
    are whole numbers, so a truth allele under 10 copies (including 0) counts only an exact match.
    """
    truth_counts = (truth_short, truth_long)
    if len(eh_counts) == 1:
        if truth_short != truth_long:
            return None
        truth_counts = (truth_short,)
    errors = [abs(eh - truth) for eh, truth in zip(eh_counts, truth_counts)]
    n_within_10_percent = sum(100 * error <= WITHIN_PERCENT_TOLERANCE * truth
                              for error, truth in zip(errors, truth_counts))
    return (max(errors) == 0, max(errors) <= NEAR_CONCORDANCE_TOLERANCE, float(np.mean(errors)),
            n_within_10_percent)


def load_truth(truth_tsv_path, regions_by_chrom):
    """Return {LocusId without "chr": (short allele, long allele, reference copies)} inside the
    high-confidence regions.

    Rows whose repeat-count columns are empty are dropped: those are loci where the assembly's
    overlapping variants could not be resolved into a genotype, which is missing truth rather than a
    reference-length call.
    """
    truth_by_locus_id = {}
    with open_maybe_gzipped(truth_tsv_path) as f:
        header = next(f).rstrip("\n").split("\t")
        column = {name: i for i, name in enumerate(header)}
        for line in f:
            fields = line.rstrip("\n").split("\t")
            short_allele = fields[column["NumRepeatsShortAllele"]]
            long_allele = fields[column["NumRepeatsLongAllele"]]
            if not short_allele or not long_allele:
                continue
            if not is_fully_callable(regions_by_chrom, fields[column["Chrom"]],
                                     int(fields[column["Start0Based"]]), int(fields[column["End"]])):
                continue
            reference_copies = fields[column["NumRepeatsInReference"]]
            truth_by_locus_id[locus_id_without_chr(fields[column["LocusId"]])] = (
                int(float(short_allele)), int(float(long_allele)),
                float(reference_copies) if reference_copies else float("nan"))
    return truth_by_locus_id


def compare_sample(catalog_bed_path, eh_json_path, truth_tsv_path, high_confidence_bed_path):
    """Return the per-locus comparison arrays for one sample, aligned to catalog row order."""
    row_number_by_locus_id = catalog_row_numbers(catalog_bed_path)
    n_loci = len(row_number_by_locus_id)
    print(f"{n_loci:,} locus definitions in {catalog_bed_path}")

    truth_by_locus_id = load_truth(truth_tsv_path, load_high_confidence_regions(high_confidence_bed_path))
    print(f"{len(truth_by_locus_id):,} loci have a truth genotype inside the high-confidence regions")

    arrays = {
        "has_call": np.zeros(n_loci, dtype=np.uint8),
        "is_quick": np.zeros(n_loci, dtype=np.uint8),
        "pok_sum": np.zeros(n_loci, dtype=np.float32),
        "pok_sum_compared": np.zeros(n_loci, dtype=np.float32),
        "ptooshort_sum": np.zeros(n_loci, dtype=np.float32),
        "ptooshort_sum_compared": np.zeros(n_loci, dtype=np.float32),
        "ptoolong_sum": np.zeros(n_loci, dtype=np.float32),
        "ptoolong_sum_compared": np.zeros(n_loci, dtype=np.float32),
        "has_truth": np.zeros(n_loci, dtype=np.uint8),
        "is_exact": np.zeros(n_loci, dtype=np.uint8),
        "is_within1": np.zeros(n_loci, dtype=np.uint8),
        "n_alleles_compared": np.zeros(n_loci, dtype=np.uint8),
        "n_alleles_within_10_percent": np.zeros(n_loci, dtype=np.uint8),
        "abs_error": np.zeros(n_loci, dtype=np.float32),
        "truth_short_allele": np.full(n_loci, np.nan, dtype=np.float32),
        "truth_long_allele": np.full(n_loci, np.nan, dtype=np.float32),
        "eh_short_allele": np.full(n_loci, np.nan, dtype=np.float32),
        "eh_long_allele": np.full(n_loci, np.nan, dtype=np.float32),
        "eh_diff_from_ref": np.full(n_loci, np.nan, dtype=np.float32),
        "truth_short_diff_from_ref": np.full(n_loci, np.nan, dtype=np.float32),
        "truth_long_diff_from_ref": np.full(n_loci, np.nan, dtype=np.float32),
        "coverage": np.zeros(n_loci, dtype=np.float32),
        "n_spanning_reads": np.zeros(n_loci, dtype=np.uint16),
        "n_flanking_reads": np.zeros(n_loci, dtype=np.uint16),
        "n_irr_reads": np.zeros(n_loci, dtype=np.uint16),
    }
    n_unknown_locus_ids = 0
    for locus_id, record in parse_eh_json(eh_json_path):
        locus_id = locus_id_without_chr(locus_id)
        row_number = row_number_by_locus_id.get(locus_id)
        if row_number is None:
            n_unknown_locus_ids += 1
            continue
        eh_counts = parse_genotype(record["genotype"])
        if eh_counts is None:
            continue
        mean_probabilities = {name: mean_or_none(record[f"{name}_values"])
                              for name in PROBABILITY_PATTERNS}
        arrays["has_call"][row_number] = 1
        arrays["is_quick"][row_number] = int(record["is_quick"])
        arrays["coverage"][row_number] = record["coverage"]
        for name in ("n_spanning_reads", "n_flanking_reads", "n_irr_reads"):
            arrays[name][row_number] = min(record[name], np.iinfo(np.uint16).max)
        for name, mean_value in mean_probabilities.items():
            if mean_value is not None:
                arrays[f"{name}_sum"][row_number] = mean_value
        truth = truth_by_locus_id.get(locus_id)
        if truth is None:
            continue
        truth_short, truth_long, reference_copies = truth
        if np.isfinite(reference_copies):
            arrays["eh_diff_from_ref"][row_number] = max(abs(allele - reference_copies)
                                                         for allele in eh_counts)
            arrays["truth_short_diff_from_ref"][row_number] = truth_short - reference_copies
            arrays["truth_long_diff_from_ref"][row_number] = truth_long - reference_copies
        comparison = compare_alleles(eh_counts, truth_short, truth_long)
        if comparison is None:
            continue
        arrays["has_truth"][row_number] = 1
        for name, mean_value in mean_probabilities.items():
            if mean_value is not None:
                arrays[f"{name}_sum_compared"][row_number] = mean_value
        (arrays["is_exact"][row_number], arrays["is_within1"][row_number], arrays["abs_error"][row_number],
         arrays["n_alleles_within_10_percent"][row_number]) = comparison
        arrays["n_alleles_compared"][row_number] = len(eh_counts)
        arrays["truth_short_allele"][row_number] = truth_short
        arrays["truth_long_allele"][row_number] = truth_long
        arrays["eh_short_allele"][row_number] = eh_counts[0]
        arrays["eh_long_allele"][row_number] = eh_counts[-1]

    if n_unknown_locus_ids:
        raise ValueError(f"{n_unknown_locus_ids:,} LocusIds in {eh_json_path} are not in "
                         f"{catalog_bed_path}, so the two were built from different catalogs")
    print(f"{arrays['has_call'].sum():,} loci genotyped, {arrays['has_truth'].sum():,} of them scored "
          f"against truth, {arrays['is_exact'].sum():,} exact")
    return arrays


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog-bed", required=True, help="BED the ExpansionHunter catalog was built from")
    parser.add_argument("--eh-json", required=True, help="ExpansionHunter output for this sample")
    parser.add_argument("--truth-tsv", required=True, help="filter_vcf_to_tandem_repeats output for this sample")
    parser.add_argument("--high-confidence-bed", required=True, help="this sample's DipCall .dip.bed.gz")
    parser.add_argument("--output-npz", required=True, help="where to write the per-locus arrays")
    args = parser.parse_args()

    np.savez_compressed(args.output_npz, **compare_sample(
        args.catalog_bed, args.eh_json, args.truth_tsv, args.high_confidence_bed))
    print(f"wrote {args.output_npz}")


if __name__ == "__main__":
    main()
