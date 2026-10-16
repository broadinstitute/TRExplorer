"""Convert a TRGT multisample VCF to a wide table of total allele lengths, in repeat units, per sample.

The output has the same layout as the trgt-lps table that
str_analysis/convert_multisample_LPS_table_to_allele_frequency_histograms.py reads, so the same converter
turns it into allele size histograms, genotype (biallelic) histograms and summary statistics:

    trid    motif    <sample1>    <sample2>    ...
    4-3074876-3074933-CAG    CAG    17,21    12,18    ...

Each sample value is the allele lengths from FORMAT/AL (the record's length in bp, without the padding base)
divided by the motif length and rounded down, in genotype order: "17,21" for a diploid call, "4" for a
hemizygous one, "." for a no-call, and "." in place of an allele TRGT did not measure (e.g. "17,.").

Which record a locus's values come from:
- A locus that is part of a variation cluster gets the cluster's record (INFO/STRUC "<VC..."), so its values
  are the total length of the whole cluster, divided by the locus's own motif length. Its single-locus record
  is skipped: TRGT's single-locus allele lengths for a repeat inside a cluster often absorb the variation of
  neighbouring repeats (e.g. at 4-3074939-3074966-CCG they track the HTT CAG expansion next to it).
- Any other locus gets its single-locus record (INFO/STRUC "<TR:...>").
A record gets one output row per motif in INFO/MOTIFS, which the converter then assigns to the record's loci
with that motif.

The script also writes the subset of the VCF's (trid, locus_id, motif, interval, vc) metadata TSV that covers
exactly the converted records. The converter requires every record in its --vcf-trid-metadata-tsv to have
a row in the input table, so the full metadata TSV (which also lists the skipped single-locus records) can't
be used.

Usage:
    python3 convert_TRGT_vcf_to_multisample_total_allele_length_table.py \\
        --input-vcf ../hprc-lps_2026-05-19/trgt-hprc.unique_trids.vcf.gz \\
        --input-interval-metadata-tsv ../hprc-lps_2026-05-19/trgt-hprc.unique_trids.interval_metadata.tsv.gz \\
        --output-table ../hprc-lps_2026-05-19/hprc-total-allele-length.unique_trids.txt.gz
"""

import argparse
import gzip
import multiprocessing
import shutil
import subprocess
from pathlib import Path

CHROMS = [f"chr{c}" for c in list(range(1, 23)) + ["X", "Y", "M"]]


def convert_sample_field_to_total_lengths_in_repeats(sample_field, motif_length):
    """Returns one sample's value for the output table from its "GT:AL" field.

    Args:
        sample_field (str): e.g. "0/1:51,63", "0:9" (hemizygous), "./.:." (no-call) or "./1:.,30"
        motif_length (int): the motif length in bp

    Returns:
        str: the allele lengths in whole repeats (AL // motif_length) in genotype order, joined by ",",
            with "." for an allele that was not measured, or "." if no allele was called at all
    """
    genotype, lengths = sample_field.split(":")
    alleles = genotype.replace("|", "/").split("/")
    if all(allele == "." for allele in alleles):
        return "."
    lengths = lengths.split(",")
    return ",".join(
        "." if allele == "." or length == "." else str(int(length) // motif_length)
        for allele, length in zip(alleles, lengths))


def convert_record_to_rows(trid, motifs, structure, sample_fields, locus_ids_in_variation_clusters):
    """Returns the output rows (lists of strings) for one VCF record: one per motif, or none if it is skipped.

    Args:
        trid (str): INFO/TRID, e.g. "4-3074876-3074933-CAG" or "VC:4:3074876-3074946"
        motifs (str): INFO/MOTIFS, e.g. "CAG" or "CAG,CAA"
        structure (str): INFO/STRUC, e.g. "<TR:4-3074876-3074933-CAG>" or "<VC:...>"
        sample_fields (list): each sample's "GT:AL" field
        locus_ids_in_variation_clusters (set): LocusIds that are part of a variation cluster record

    Returns:
        list: [trid, motif, <one value per sample>] rows. A single-locus record of a locus that is part of a
            variation cluster gets none, since the locus's values come from the cluster's record.
    """
    if structure.startswith("<TR:") and trid in locus_ids_in_variation_clusters:
        return []
    return [[trid, motif] + [convert_sample_field_to_total_lengths_in_repeats(f, len(motif)) for f in sample_fields]
            for motif in motifs.split(",")]


def convert_chromosome(vcf_path, chrom, output_path, locus_ids_in_variation_clusters):
    """Writes one chromosome's output rows to a gzipped table chunk (without a header).

    Returns:
        set: the TRIDs of the records that produced output rows
    """
    command = ["bcftools", "query", "-r", chrom, "-f", "%INFO/TRID\t%INFO/MOTIFS\t%INFO/STRUC[\t%GT:%AL]\n",
               str(vcf_path)]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, text=True)
    converted_trids = set()
    with gzip.open(output_path, "wt") as out:
        for line in process.stdout:
            trid, motifs, structure, *sample_fields = line.rstrip("\n").split("\t")
            rows = convert_record_to_rows(trid, motifs, structure, sample_fields, locus_ids_in_variation_clusters)
            for row in rows:
                out.write("\t".join(row) + "\n")
            if rows:
                converted_trids.add(trid)
    if process.wait() != 0:
        raise RuntimeError(f"bcftools query failed for {chrom}")
    print(f"{chrom}: converted {len(converted_trids):,} records", flush=True)
    return converted_trids


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input-vcf", required=True, type=Path)
    parser.add_argument("--input-interval-metadata-tsv", required=True, type=Path,
                        help="The VCF's (trid, locus_id, motif, interval, vc) TSV from "
                             "str_analysis.extract_trid_metadata_from_TRGT_vcf")
    parser.add_argument("--output-table", required=True, type=Path, help="Output .txt.gz path")
    parser.add_argument("--output-interval-metadata-tsv", type=Path,
                        help="Default: the output table path with .txt.gz replaced by .interval_metadata.tsv.gz")
    parser.add_argument("--n-processes", type=int, default=8)
    args = parser.parse_args()

    output_metadata_tsv = args.output_interval_metadata_tsv or Path(
        str(args.output_table).removesuffix(".txt.gz") + ".interval_metadata.tsv.gz")
    with gzip.open(args.input_interval_metadata_tsv, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        if header != ["trid", "locus_id", "motif", "interval", "vc"]:
            raise ValueError(f"Unexpected columns in {args.input_interval_metadata_tsv}: {header}")
        locus_ids_in_variation_clusters = {fields[1] for fields in (line.rstrip("\n").split("\t") for line in f)
                                           if fields[4]}
    print(f"{len(locus_ids_in_variation_clusters):,} loci are part of a variation cluster record")

    sample_ids = subprocess.run(["bcftools", "query", "-l", str(args.input_vcf)], check=True, capture_output=True,
                                text=True).stdout.split()
    print(f"Found {len(sample_ids):,} samples in {args.input_vcf}")

    chunk_paths = [Path(f"{args.output_table}.{chrom}.chunk.gz") for chrom in CHROMS]
    with multiprocessing.Pool(args.n_processes) as pool:
        converted_trids_by_chrom = pool.starmap(
            convert_chromosome,
            [(args.input_vcf, chrom, path, locus_ids_in_variation_clusters) for chrom, path in zip(CHROMS, chunk_paths)],
            chunksize=1)

    # Concatenated gzip members are themselves a valid gzip file, so the chunks are appended as they are.
    with open(args.output_table, "wb") as out:
        out.write(gzip.compress(("\t".join(["trid", "motif"] + sample_ids) + "\n").encode()))
        for path in chunk_paths:
            with open(path, "rb") as chunk:
                shutil.copyfileobj(chunk, out)
            path.unlink()
    converted_trids = set().union(*converted_trids_by_chrom)
    print(f"Wrote rows for {len(converted_trids):,} records to {args.output_table}")

    metadata_rows_written = 0
    with gzip.open(args.input_interval_metadata_tsv, "rt") as f, gzip.open(output_metadata_tsv, "wt") as out:
        out.write(f.readline())
        for line in f:
            if line.split("\t", 1)[0] in converted_trids:
                out.write(line)
                metadata_rows_written += 1
    print(f"Wrote {metadata_rows_written:,} metadata rows to {output_metadata_tsv}")


if __name__ == "__main__":
    main()
