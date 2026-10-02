"""Summarize the full run's SpliceAI scores into one row per locus, for the TRExplorer BigQuery table.

Reads the full run's input chunks (full_run_inputs/<gene_set>/, which list every simulated allele of
every locus and the target sizes it stands for) and its downloaded outputs (the run folder from
`modal volume get spliceai-tr-full-run /outputs/<gene_set>/<run_id> <existing dir>`, which writes
<existing dir>/<run_id>/); only the chunks in the
input chunks' manifest.json with a summary for the same input SHA-256 are used, and every one of them must
have one), and writes a TSV with:

    LocusId                                    the TRExplorer locus ID, e.g. 1-12345-12400-CAG
    SpliceAI_MaxDeltaScore                     largest SpliceAI delta score (acceptor or donor, gain or
                                               loss, on the transcript the server selects) of any
                                               simulated allele at the locus; 0 if none reached 0.01
    SpliceAI_MaxDeltaScoreAlleleSize           which simulated size gave it: 2.5pct, 97.5pct, 99.5pct
                                               (HPRC256 total allele length percentiles), or +1x, +2x,
                                               +3x (the 99.5th percentile plus that multiple of the
                                               motif's range); empty if no allele reached 0.01
    SpliceAI_MinAlleleSizeThatAffectsSplicing  the smallest simulated size whose delta score is at least
                                               AFFECTS_SPLICING_MIN_DELTA_SCORE (0.2), with the same
                                               values; empty if none
    SpliceAI_DeltaScoreByRepeatCount           every simulated allele, sorted by repeat count, as
                                               repeat_count:max_delta_score, followed by the type of the
                                               largest change (AG acceptor gain, AL acceptor loss, DG
                                               donor gain, DL donor loss) when it is at least 0.01, e.g.
                                               "18:0.000,25:0.030AG,40:0.310AG,95:0.880AL" (3 decimals,
                                               as SpliceAI-lookup reports the scores)

The sizes are in increasing length order (2.5pct < 97.5pct < 99.5pct < +1x < +2x < +3x). One allele
can stand for several sizes when they round to the same repeat count; it is then named by the smallest.
Repeat counts are the hg38 tract length divided by the motif length (whole units, as in the allele
design) plus the allele's change in units. An allele with no stored record scored below 0.01 on the
selected transcript (the pipeline stores only alleles at or above it), so its delta score is 0. An
allele whose scoring returned an error (listed in its chunk's summary.json) has no score and is left out
of every column; a locus with no allele scored is left out of the output.

Loci the run did not score (not polymorphic in HPRC256, or outside every GENCODE v50 basic transcript)
are not in the output, so their columns stay empty (NULL) in BigQuery.

Usage:
    python3 summarize_full_run_spliceai_scores_per_locus.py --gene-set basic --outputs-dir full_run_outputs/basic/<run_id>
"""
import argparse
import collections
import gzip
import json
import os

from make_tr_expansion_benchmark_alleles import HIGHEST_PERCENTILE, TARGET_LABELS

HERE = os.path.dirname(os.path.abspath(__file__))
AFFECTS_SPLICING_MIN_DELTA_SCORE = 0.2
MIN_DELTA_SCORE_TO_NAME_THE_CHANGE = 0.01
DELTA_SCORES = ("DS_AG", "DS_AL", "DS_DG", "DS_DL")
OUTPUT_COLUMNS = ("LocusId", "SpliceAI_MaxDeltaScore", "SpliceAI_MaxDeltaScoreAlleleSize",
                  "SpliceAI_MinAlleleSizeThatAffectsSplicing", "SpliceAI_DeltaScoreByRepeatCount")


def allele_size_name(target_label):
    """Returns the column value for a target label, e.g. "97.5th percentile" -> "97.5pct", "99.5th + 2x motif range" -> "+2x"."""
    if target_label.startswith(f"{HIGHEST_PERCENTILE}th + "):
        return "+" + target_label.split(" + ")[1].split()[0]
    return target_label.split("th percentile")[0] + "pct"


def summarize_locus(locus_id, alleles, target_labels, records_by_variant, unscored_variants=frozenset()):
    """Returns the output row (a dict of OUTPUT_COLUMNS) for one locus of the full run, or None if none of its
    alleles was scored.

    Args:
        locus_id (str): the full run's locus ID, e.g. chr1-12345-12400-CAG
        alleles (list): the locus's simulated variants, chrom-pos-ref-alt
        target_labels (list): for each allele, the TARGET_LABELS it stands for, smallest first
        records_by_variant (dict): the stored records of this locus's alleles, by variant
        unscored_variants (set): variants whose scoring returned an error; they are left out
    """
    chrom, start, end, motif = locus_id.rsplit("-", 3)
    reference_count = (int(end) - int(start)) // len(motif)
    size_order = {label: i for i, label in enumerate(TARGET_LABELS)}
    scored = []
    for variant, labels in zip(alleles, target_labels):
        if variant in unscored_variants:
            continue
        _, _, ref, alt = variant.split("-")
        repeat_count = reference_count + (len(alt) - len(ref)) // len(motif)
        record = records_by_variant.get(variant)
        delta, change = 0.0, ""
        if record:
            change = max(DELTA_SCORES, key=lambda k: float(record[k]))
            delta = float(record[change])
            change = change[3:] if delta >= MIN_DELTA_SCORE_TO_NAME_THE_CHANGE else ""
        smallest_label = min(labels, key=size_order.__getitem__)
        scored.append((size_order[smallest_label], repeat_count, delta, change, allele_size_name(smallest_label)))
    if not scored:
        return None
    scored.sort()
    max_delta = max(s[2] for s in scored)
    # Ties go to the smallest size
    max_size = next(s[4] for s in scored if s[2] == max_delta) if max_delta >= MIN_DELTA_SCORE_TO_NAME_THE_CHANGE else ""
    min_affecting_size = next((s[4] for s in scored if s[2] >= AFFECTS_SPLICING_MIN_DELTA_SCORE), "")
    by_repeat_count = ",".join(f"{count}:{delta:.3f}{change}" for _, count, delta, change, _ in sorted(scored, key=lambda s: s[1]))
    return {
        "LocusId": locus_id[3:] if locus_id.startswith("chr") else locus_id,
        "SpliceAI_MaxDeltaScore": f"{max_delta:.3f}",
        "SpliceAI_MaxDeltaScoreAlleleSize": max_size,
        "SpliceAI_MinAlleleSizeThatAffectsSplicing": min_affecting_size,
        "SpliceAI_DeltaScoreByRepeatCount": by_repeat_count,
    }


def read_current_chunk_names(inputs_dir, outputs_dir):
    """Returns the manifest's chunk names, after checking the outputs hold a summary for each chunk's exact input."""
    manifest = json.load(open(os.path.join(inputs_dir, "manifest.json")))
    missing = []
    for chunk in manifest["chunks"]:
        summary_path = os.path.join(outputs_dir, f"{chunk['chunk']}.summary.json")
        if not os.path.exists(summary_path) or json.load(open(summary_path)).get("input_sha256") != chunk["sha256"]:
            missing.append(chunk["chunk"])
    if missing:
        raise SystemExit(f"{len(missing):,} of {len(manifest['chunks']):,} chunks have no output for their current input "
                         f"in {outputs_dir} (e.g. {missing[0]}); finish or download the run first.")
    return [chunk["chunk"] for chunk in manifest["chunks"]]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--gene-set", choices=["basic", "comprehensive"], required=True)
    parser.add_argument("--outputs-dir", required=True, help="The downloaded run folder (chunk .jsonl.gz and .summary.json files)")
    parser.add_argument("--inputs-dir", help="The run's input chunks and manifest. Default: full_run_inputs/<gene_set> in this folder")
    parser.add_argument("--output-tsv", help="Default: TRExplorer_SpliceAI_summary.<gene_set>.tsv.gz in this folder")
    args = parser.parse_args()

    inputs_dir = args.inputs_dir or os.path.join(HERE, "full_run_inputs", args.gene_set)
    output_tsv = args.output_tsv or os.path.join(HERE, f"TRExplorer_SpliceAI_summary.{args.gene_set}.tsv.gz")
    counts = collections.Counter()
    with gzip.open(output_tsv, "wt") as out:
        out.write("\t".join(OUTPUT_COLUMNS) + "\n")
        for chunk in read_current_chunk_names(inputs_dir, args.outputs_dir):
            with gzip.open(os.path.join(inputs_dir, f"{chunk}.json.gz"), "rt") as f:
                loci = json.load(f)["loci"]
            records_by_locus = collections.defaultdict(dict)
            with gzip.open(os.path.join(args.outputs_dir, f"{chunk}.jsonl.gz"), "rt") as f:
                for line in f:
                    record = json.loads(line)
                    records_by_locus[record["locus"]][record["variant"]] = record
            unscored_variants = {variant for variant, _ in json.load(open(os.path.join(args.outputs_dir, f"{chunk}.summary.json")))["errors"]}
            counts["unscored alleles"] += len(unscored_variants)
            for locus in loci:
                row = summarize_locus(locus["locus"], locus["alleles"], locus["target_labels"], records_by_locus[locus["locus"]],
                                      unscored_variants)
                if row is None:
                    counts["unscored loci"] += 1
                    continue
                out.write("\t".join(row[c] for c in OUTPUT_COLUMNS) + "\n")
                counts["loci"] += 1
                counts["affected"] += bool(row["SpliceAI_MinAlleleSizeThatAffectsSplicing"])
    print(f"Wrote {counts['loci']:,} loci to {output_tsv}; {counts['affected']:,} with a delta score >= "
          f"{AFFECTS_SPLICING_MIN_DELTA_SCORE}; left out {counts['unscored alleles']:,} alleles whose scoring returned an "
          f"error, and {counts['unscored loci']:,} loci with no allele scored")


if __name__ == "__main__":
    main()
