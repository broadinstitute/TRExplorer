"""Parse TRGT VCF to compute joint distributions of allele size vs. purity (AP) and
allele size vs. methylation (AM) for each single-motif locus.

Two gzipped TSVs are emitted (one per metric). One row is emitted per VCF record (=
per TRGT catalog interval). Each row carries three identifier columns: ``locus_id``
(the INFO/TRID), ``interval`` = ``{chrom}:{vcf_start_0based}-{vcf_end_1based}``
(always set), and ``vc`` = the inner span from ``INFO/STRUC`` if the row was
genotyped as part of a variation cluster (``<VC:...>``) or the empty string for
an isolated TR (``<TR:...>``). These three columns together uniquely identify
the TRGT interval; a ``locus_id`` that appears in multiple intervals (once as a
standalone TR and once or more inside a VC) emits multiple rows that differ in
``interval`` and ``vc``.

The first distribution column is always the unsuffixed "All samples" value.
Stratification options mirror those of
str_analysis/convert_multisample_LPS_table_to_allele_frequency_histograms.py: when both
--stratify-by-population and --stratify-by-sex are set, each sample additionally
contributes to its Pop_Sex cell plus the Pop row-marginal and the Sex column-marginal,
adding 17 stratum-suffixed columns (5 populations + 2 sexes + 10 Pop_Sex cells).

Output column format:
    locus_id, interval, vc, AlleleSizeAndPurityDistribution, AlleleSizeAndPurityDistribution__<stratum>, ...
    locus_id, interval, vc, AlleleSizeAndMethylationDistribution, AlleleSizeAndMethylationDistribution__<stratum>, ...

Output filenames (suffixes match the LPS script):
    <stem>.allele_size_purity.stratified[.only_<pop>][.only_<sex>][.by_population][.by_sex].{N}_samples.tsv.gz
    <stem>.methylation.stratified[.only_<pop>][.only_<sex>][.by_population][.by_sex].{N}_samples.tsv.gz
"""

import argparse
import collections
import gzip
from pathlib import Path

import pandas as pd


def format_distribution(counter):
    """Format a Counter of (repeat_count, binned_value) pairs as 'rc/val:count,...' string,
    sorted by repeat_count then value."""
    if not counter:
        return ""
    return ",".join(
        f"{rc}/{val}:{count}"
        for (rc, val), count in sorted(counter.items())
    )


def bin_value(v):
    """Bin a float value to nearest 0.02 increment."""
    return round(round(v / 0.02) * 0.02, 2)


def parse_struc_vc_span(struc):
    """Returns the inner VC span from ``INFO/STRUC`` or ``""`` for non-VC rows.

    For an isolated TR row, ``STRUC`` is ``<TR:locus-id>`` and this returns
    ``""``. For a variation-cluster row, ``STRUC`` is ``<VC:chrom:start-end>``
    and this returns ``chrom:start-end``. The VC span follows the no-``chr``
    convention used by TRIDs in this catalog (e.g. ``13:102161564-102161724``).
    """
    if not struc.startswith("<VC:"):
        return ""
    end_idx = struc.find(">", 4)
    if end_idx == -1:
        return ""
    return struc[4:end_idx]


def _strip_chr(chrom):
    """Strips a leading ``chr`` from a chromosome name (matching TRID/STRUC convention)."""
    return chrom[3:] if chrom.startswith("chr") else chrom


def get_vcf_sample_ids(input_vcf):
    """Read the VCF header and return the list of sample IDs in column order."""
    opener = gzip.open if str(input_vcf).endswith(".gz") else open
    with opener(input_vcf, "rt") as f:
        for line in f:
            if line.startswith("#CHROM"):
                return line.rstrip("\n").split("\t")[9:]
            if not line.startswith("#"):
                break
    raise ValueError(f"No #CHROM header line found in {input_vcf}")


def build_sample_id_to_strata(df_metadata, sample_ids_to_include, stratify_by_population, stratify_by_sex):
    """Build a {sample_id: [stratum_label, ...]} mapping. Every sample always contributes
    to the unsuffixed "All" stratum (label ""). With both stratify flags, each sample
    additionally contributes to its Pop_Sex cell plus the Pop and Sex marginals."""
    df_for_strata = df_metadata[df_metadata.SampleId.isin(set(sample_ids_to_include))]
    if stratify_by_population and stratify_by_sex:
        return {
            sid: ["", f"{pop}_{sex}", pop, sex]
            for sid, pop, sex in zip(df_for_strata.SampleId, df_for_strata.Population, df_for_strata.Sex)
        }
    if stratify_by_population:
        return {sid: ["", pop] for sid, pop in zip(df_for_strata.SampleId, df_for_strata.Population)}
    if stratify_by_sex:
        return {sid: ["", sex] for sid, sex in zip(df_for_strata.SampleId, df_for_strata.Sex)}
    return {sid: [""] for sid in sample_ids_to_include}


def parse_vcf_and_compute_distributions(input_vcf, sample_index_to_strata, num_loci=None):
    """Parse VCF and yield (locus_ids, interval, vc, purity_counters, methylation_counters) per locus.

    ``locus_ids`` is the list of LocusIds (each ``chrom-start-end-motif``)
    from ``INFO/TRID`` whose motif suffix matches the record's single
    ``INFO/MOTIFS`` value — exactly one LocusId for a standalone TR row, but
    multiple LocusIds for a compound TRID that lists several loci with the
    same motif. The caller writes one output row per LocusId; all LocusIds
    of one VCF record share the same purity/methylation distributions plus
    interval/vc, differing only in the LocusId column.

    purity_counters and methylation_counters are dicts mapping stratum label -> Counter
    of (repeat_count, binned_value) -> count. ``interval`` is
    ``{chrom}:{vcf_start_0based}-{vcf_end_1based}`` (no ``chr`` prefix); ``vc`` is
    the inner ``<VC:...>`` span (e.g. ``13:102161564-102161724``) or ``""`` for an
    isolated TR row.
    """
    opener = gzip.open if str(input_vcf).endswith(".gz") else open
    with opener(input_vcf, "rt") as f:
        rows_processed = 0
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                continue

            if num_loci is not None and rows_processed >= num_loci:
                break

            fields = line.rstrip("\n").split("\t")
            chrom = fields[0]
            pos = fields[1]
            info = fields[7]

            motifs = None
            trid = None
            end = None
            struc = ""
            for kv in info.split(";"):
                if kv.startswith("MOTIFS="):
                    motifs = kv[7:]
                elif kv.startswith("TRID="):
                    trid = kv[5:]
                elif kv.startswith("END="):
                    end = kv[4:]
                elif kv.startswith("STRUC="):
                    struc = kv[6:]

            if motifs is None or trid is None:
                continue
            if "," in motifs:
                continue  # skip multi-motif loci

            try:
                vcf_start_0based = int(pos)
                vcf_end_1based = int(end) if end is not None else None
            except ValueError:
                continue
            if vcf_end_1based is None:
                continue

            interval = f"{_strip_chr(chrom)}:{vcf_start_0based}-{vcf_end_1based}"
            vc = parse_struc_vc_span(struc)

            motif_length = len(motifs)

            fmt_fields = fields[8].split(":")
            try:
                al_idx = fmt_fields.index("AL")
                ap_idx = fmt_fields.index("AP")
            except ValueError:
                continue

            am_idx = fmt_fields.index("AM") if "AM" in fmt_fields else None

            purity_counters = collections.defaultdict(collections.Counter)
            methylation_counters = collections.defaultdict(collections.Counter)

            for sample_idx, sample_field in enumerate(fields[9:]):
                strata = sample_index_to_strata.get(sample_idx)
                if not strata:
                    continue

                parts = sample_field.split(":")
                gt = parts[0]
                if gt == "." or gt == "./.":
                    continue

                try:
                    al_values = parts[al_idx].split(",")
                    ap_values = parts[ap_idx].split(",")
                except IndexError:
                    continue

                am_values = None
                if am_idx is not None:
                    try:
                        am_values = parts[am_idx].split(",")
                    except IndexError:
                        pass

                for i in range(len(al_values)):
                    try:
                        al = int(al_values[i])
                    except (ValueError, IndexError):
                        continue

                    repeat_count = al // motif_length

                    try:
                        ap = float(ap_values[i])
                    except (ValueError, IndexError):
                        pass
                    else:
                        binned_ap = bin_value(ap)
                        for stratum in strata:
                            purity_counters[stratum][(repeat_count, binned_ap)] += 1

                    if am_values is not None and i < len(am_values) and am_values[i] != ".":
                        try:
                            am = float(am_values[i])
                        except ValueError:
                            pass
                        else:
                            binned_am = bin_value(am)
                            for stratum in strata:
                                methylation_counters[stratum][(repeat_count, binned_am)] += 1

            rows_processed += 1
            if rows_processed % 100000 == 0:
                print(f"Processed {rows_processed:,d} rows...", flush=True)

            # Expand TRID into the LocusIds whose suffix matches the row's motif.
            # For an isolated TR this is the TRID itself; for a single-motif
            # compound TRID this is every listed LocusId ending in -{motif}.
            suffix = f"-{motifs}"
            matching_locus_ids = [lid for lid in trid.split(",") if lid.endswith(suffix)]
            if not matching_locus_ids:
                continue
            yield matching_locus_ids, interval, vc, purity_counters, methylation_counters

    print(f"Done. Processed {rows_processed:,d} rows total.", flush=True)


def main():
    parser = argparse.ArgumentParser(
        description="Compute allele size vs. purity and methylation distributions from TRGT VCF.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--input-vcf", default="trgt-hprc.vcf.gz",
                        help="Path to TRGT VCF file")
    parser.add_argument("--sample-metadata-tsv",
                        default="https://storage.googleapis.com/tandem-repeat-catalog/1kGP_metadata.tsv",
                        help="Sample ancestry metadata TSV file (local path or URL)")
    parser.add_argument("--population", choices=["AFR", "AMR", "EAS", "EUR", "SAS"],
                        help="If specified, only include samples from this population")
    parser.add_argument("--sex", choices=["male", "female"],
                        help="If specified, only include samples from this sex")
    parser.add_argument("--stratify-by-population", action="store_true",
                        help="Add per-population distribution columns")
    parser.add_argument("--stratify-by-sex", action="store_true",
                        help="Add per-sex distribution columns")
    parser.add_argument("-n", "--num-samples", type=int, default=None,
                        help="Limit number of samples to process")
    parser.add_argument("-l", "--num-loci", type=int, default=None,
                        help="Limit number of VCF data rows (loci) to process")
    args = parser.parse_args()

    if args.stratify_by_population and args.population:
        parser.error("--stratify-by-population and --population are mutually exclusive")
    if args.stratify_by_sex and args.sex:
        parser.error("--stratify-by-sex and --sex are mutually exclusive")

    print(f"Reading sample IDs from {args.input_vcf}")
    vcf_sample_ids = get_vcf_sample_ids(args.input_vcf)
    print(f"Found {len(vcf_sample_ids):,d} samples in VCF")

    df_metadata = pd.read_table(args.sample_metadata_tsv, dtype={"SampleId": str})
    metadata_ids = set(df_metadata.SampleId)
    unexpected = set(vcf_sample_ids) - metadata_ids
    if unexpected:
        parser.error(f"{len(unexpected):,d} VCF sample id(s) not found in {args.sample_metadata_tsv}: "
                     f"{', '.join(sorted(unexpected))}")

    df_included = df_metadata[df_metadata.SampleId.isin(set(vcf_sample_ids))]
    if args.population:
        df_included = df_included[df_included.Population == args.population]
        print(f"Kept {len(df_included):,d} samples from population {args.population}")
    if args.sex:
        df_included = df_included[df_included.Sex == args.sex]
        print(f"Kept {len(df_included):,d} samples from sex {args.sex}")

    valid_ids = set(df_included.SampleId)
    sample_ids_to_include = [s for s in vcf_sample_ids if s in valid_ids]
    if args.num_samples is not None and len(sample_ids_to_include) > args.num_samples:
        sample_ids_to_include = sample_ids_to_include[:args.num_samples]

    if not sample_ids_to_include:
        parser.error("No samples remain after applying --population/--sex/--num-samples filters; nothing to do.")

    sample_id_to_strata = build_sample_id_to_strata(
        df_included, sample_ids_to_include, args.stratify_by_population, args.stratify_by_sex
    )
    strata_labels = sorted({label for labels in sample_id_to_strata.values() for label in labels})

    sample_index_to_strata = {}
    for sample_idx, sample_id in enumerate(vcf_sample_ids):
        if sample_id in sample_id_to_strata:
            sample_index_to_strata[sample_idx] = sample_id_to_strata[sample_id]

    base_stem = str(Path(args.input_vcf)).replace(".vcf.gz", "").replace(".vcf", "")
    suffix = ".stratified"
    if args.population:
        suffix += f".only_{args.population}"
    if args.sex:
        suffix += f".only_{args.sex}"
    if args.stratify_by_population:
        suffix += ".by_population"
    if args.stratify_by_sex:
        suffix += ".by_sex"
    suffix += f".{len(sample_ids_to_include)}_samples.tsv.gz"
    purity_path = Path(base_stem + ".allele_size_purity" + suffix)
    methylation_path = Path(base_stem + ".methylation" + suffix)

    purity_header = ["locus_id", "interval", "vc"] + [
        f"AlleleSizeAndPurityDistribution__{label}" if label else "AlleleSizeAndPurityDistribution"
        for label in strata_labels
    ]
    methylation_header = ["locus_id", "interval", "vc"] + [
        f"AlleleSizeAndMethylationDistribution__{label}" if label else "AlleleSizeAndMethylationDistribution"
        for label in strata_labels
    ]

    print(f"Writing data from {len(sample_ids_to_include):,d} samples to:")
    print(f"  {purity_path}")
    print(f"  {methylation_path}")

    purity_rows_written = 0
    methylation_rows_written = 0
    with gzip.open(purity_path, "wt") as purity_out, gzip.open(methylation_path, "wt") as methylation_out:
        purity_out.write("\t".join(purity_header) + "\n")
        methylation_out.write("\t".join(methylation_header) + "\n")

        for locus_ids, interval, vc, purity_counters, methylation_counters in parse_vcf_and_compute_distributions(
            args.input_vcf, sample_index_to_strata, args.num_loci
        ):
            purity_values = [format_distribution(purity_counters.get(label, collections.Counter())) for label in strata_labels]
            has_purity = any(purity_values)
            methylation_values = [format_distribution(methylation_counters.get(label, collections.Counter())) for label in strata_labels]
            has_methylation = any(methylation_values)
            if not has_purity and not has_methylation:
                continue
            for locus_id in locus_ids:
                row_prefix = f"{locus_id}\t{interval}\t{vc}"
                if has_purity:
                    purity_out.write(f"{row_prefix}\t" + "\t".join(purity_values) + "\n")
                    purity_rows_written += 1
                if has_methylation:
                    methylation_out.write(f"{row_prefix}\t" + "\t".join(methylation_values) + "\n")
                    methylation_rows_written += 1

    print(f"Wrote {purity_rows_written:,d} rows to {purity_path}")
    print(f"Wrote {methylation_rows_written:,d} rows to {methylation_path}")


if __name__ == "__main__":
    main()
