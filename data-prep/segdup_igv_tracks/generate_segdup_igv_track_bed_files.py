"""Write the hg38 segmental duplications (segdups) of UCSC and of Vollger et al. 2022 as BED files for TRExplorer's
UCSC SegDups and Vollger 2022 SegDups IGV tracks.

The inputs (see run.sh for where they come from) list each duplicated pair twice, once from each copy's side:
- genomicSuperDups.txt.gz: UCSC's hg38 genomicSuperDups table, made with the WGAC method (Bailey et al. 2001, 2002)
- hg38.chr_only.SDs.bed.gz: the hg38 segdups of Vollger et al. 2022, made with SEDEF

Each output row is one copy. Its color is the one that the UCSC Genome Browser's segdup track uses for its sequence
identity (the fracMatch column: the fraction of aligned bases that match, not counting indels): gray (#AAAAAA) below
0.98, yellow (#D1C800) from 0.98 to 0.99, and orange (#E78200) from 0.99 up. Its strand is the other copy's orientation
relative to this one. Its name column holds GFF-style tags, which the #gffTags header line tells IGV to list in the
popup of a clicked segdup: the other copy's location (1-based) and the sequence identity.

The rows kept are the segdups that TRExplorer's UCSC_SegDupMaxSequenceIdentity and Vollger2022_SegDupMaxSequenceIdentity
columns come from: rows whose copy is on a main chromosome, as are all TRExplorer loci, leaving out rows whose other
copy overlaps a chrY pseudoautosomal region (PAR). The hg38 reference masks the chrY PARs with Ns, so reads can't align
to a copy there. Rows whose own copy lies entirely inside a chrY PAR are left out too, since no locus can be there, but
a chrY copy that extends past the end of a PAR is kept, as it is for the columns. Only the Vollger et al. 2022 file has
rows with a copy in a chrY PAR.
"""

import argparse
import collections
import gzip
import os
from pathlib import Path

MAIN_CHROMOSOMES = {f"chr{chrom}" for chrom in list(range(1, 23)) + ["X", "Y", "M"]}

# The hg38 chrY PARs as (chromosome, 0-based start, end), from the Genome Reference Consortium's GRCh38 assembly regions
# file: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_assembly_regions.txt
HG38_CHRY_PSEUDOAUTOSOMAL_REGIONS = [("chrY", 10000, 2781479), ("chrY", 56887902, 57217415)]

# For each source of segdups: the BED file it's written to, the 0-based column index of each value read from its file,
# and the number of columns the file's rows have. "strand" is the column with the other copy's orientation relative to
# this copy (the Vollger et al. 2022 files' strand1 column is always +, and their strand2 column has the orientation).
SEGDUP_SOURCES = {
    "ucsc": {"output_bed_filename": "UCSC_genomicSuperDups.hg38.bed", "chrom": 1, "start_0based": 2, "end_1based": 3,
             "strand": 6, "other_chrom": 7, "other_start_0based": 8, "other_end_1based": 9, "fracMatch": 26,
             "num_columns": 30},
    "vollger2022": {"output_bed_filename": "Vollger2022_SDs.hg38.bed", "chrom": 0, "start_0based": 1, "end_1based": 2,
                    "strand": 13, "other_chrom": 9, "other_start_0based": 10, "other_end_1based": 11, "fracMatch": 23,
                    "num_columns": 44},
}

# The colors of the UCSC Genome Browser's segdup track as BED itemRgb values, each with the lowest sequence identity
# that gets it, from the highest cutoff down
SEQUENCE_IDENTITY_COLORS = [(0.99, "231,130,0"), (0.98, "209,200,0"), (0, "170,170,170")]


def get_color_for_sequence_identity(sequence_identity):
    """Returns the BED itemRgb color for a segdup's sequence identity.

    Args:
        sequence_identity (float): the segdup's fracMatch

    Returns:
        str: the itemRgb value, eg. "231,130,0"
    """
    return next(color for min_sequence_identity, color in SEQUENCE_IDENTITY_COLORS
                if sequence_identity >= min_sequence_identity)


def overlaps_a_chry_pseudoautosomal_region(chrom, start_0based, end_1based):
    """Returns whether an interval overlaps one of the hg38 chrY PARs by at least 1 bp.

    Args:
        chrom (str): chromosome of the interval, eg. chrY
        start_0based (int): start of the interval
        end_1based (int): end of the interval

    Returns:
        bool: True if the interval overlaps a chrY PAR
    """
    return any(chrom == par_chrom and start_0based < par_end and end_1based > par_start
               for par_chrom, par_start, par_end in HG38_CHRY_PSEUDOAUTOSOMAL_REGIONS)


def is_inside_a_chry_pseudoautosomal_region(chrom, start_0based, end_1based):
    """Returns whether an interval lies entirely inside one of the hg38 chrY PARs.

    Args:
        chrom (str): chromosome of the interval, eg. chrY
        start_0based (int): start of the interval
        end_1based (int): end of the interval

    Returns:
        bool: True if the interval is inside a chrY PAR
    """
    return any(chrom == par_chrom and start_0based >= par_start and end_1based <= par_end
               for par_chrom, par_start, par_end in HG38_CHRY_PSEUDOAUTOSOMAL_REGIONS)


def read_segdup_bed_rows(segdups_path, segdups_source):
    """Reads a UCSC genomicSuperDups table or a Vollger et al. 2022 SD BED file and returns the output BED rows of the
    segdups that are kept.

    Args:
        segdups_path (str): path of the gzipped segdup file
        segdups_source (str): a key of SEGDUP_SOURCES

    Returns:
        tuple: (list of BED rows sorted by chromosome and position, each a list of 9 strings; Counter with the number of
            rows left out for each reason)
    """
    columns = SEGDUP_SOURCES[segdups_source]
    bed_rows = []
    counts_of_rows_left_out = collections.Counter()
    with gzip.open(segdups_path, "rt") as f:
        for line_number, line in enumerate(f, start=1):
            fields = line.rstrip("\n").split("\t")
            if line.startswith("#"):
                # the Vollger et al. 2022 files' header line names their columns, so check that they're where expected
                if fields[columns["fracMatch"]] != "fracMatch" or fields[columns["strand"]] != "strand2":
                    raise ValueError(f"{segdups_path} header doesn't have fracMatch in column "
                                     f"{columns['fracMatch'] + 1} and strand2 in column {columns['strand'] + 1}")
                continue

            if len(fields) != columns["num_columns"]:
                raise ValueError(f"{segdups_path} line {line_number} has {len(fields)} columns rather than the "
                                 f"{columns['num_columns']} of the {segdups_source} segdup files")

            # rounded to the 6 decimals of the files' values, since the Vollger et al. 2022 files write some of them with
            # floating point noise (eg. 0.9857549999999999 for 0.985755)
            sequence_identity = round(float(fields[columns["fracMatch"]]), 6)
            if not 0 < sequence_identity <= 1:
                raise ValueError(f"{segdups_path} line {line_number} has fracMatch {sequence_identity}, which is not a "
                                 f"fraction")
            strand = fields[columns["strand"]]
            if strand not in ("+", "-"):
                raise ValueError(f"{segdups_path} line {line_number} has strand {strand} rather than + or -")

            chrom = fields[columns["chrom"]]
            start_0based = int(fields[columns["start_0based"]])
            end_1based = int(fields[columns["end_1based"]])
            other_chrom = fields[columns["other_chrom"]]
            other_start_0based = int(fields[columns["other_start_0based"]])
            other_end_1based = int(fields[columns["other_end_1based"]])
            if chrom not in MAIN_CHROMOSOMES:
                counts_of_rows_left_out["not on a main chromosome"] += 1
                continue
            if overlaps_a_chry_pseudoautosomal_region(other_chrom, other_start_0based, other_end_1based):
                counts_of_rows_left_out["whose other copy overlaps a chrY PAR"] += 1
                continue
            if is_inside_a_chry_pseudoautosomal_region(chrom, start_0based, end_1based):
                counts_of_rows_left_out["inside a chrY PAR"] += 1
                continue

            bed_rows.append([
                chrom, str(start_0based), str(end_1based),
                f"OtherCopy={other_chrom}:{other_start_0based + 1}-{other_end_1based};"
                f"SequenceIdentity={sequence_identity:g}",
                "0", strand, str(start_0based), str(end_1based), get_color_for_sequence_identity(sequence_identity),
            ])

    bed_rows.sort(key=lambda bed_row: (bed_row[0], int(bed_row[1]), int(bed_row[2])))
    return bed_rows, counts_of_rows_left_out


def write_bgzipped_and_indexed_bed(bed_rows, output_bed_path):
    """Writes BED rows after a #gffTags header line, then bgzips the file and indexes it with tabix.

    Args:
        bed_rows (list): BED rows sorted by chromosome and position, each a list of strings
        output_bed_path (Path): path of the BED file, which becomes <output_bed_path>.gz and <output_bed_path>.gz.tbi
    """
    with open(output_bed_path, "wt") as f:
        f.write("#gffTags\n")
        for bed_row in bed_rows:
            f.write("\t".join(bed_row) + "\n")

    for command in f"bgzip -f {output_bed_path}", f"tabix -f -p bed {output_bed_path}.gz":
        if os.system(command) != 0:
            raise RuntimeError(f"Command failed: {command}")


def main():
    parser = argparse.ArgumentParser(description="Write the UCSC and Vollger et al. 2022 hg38 segdups as BED files for "
                                                 "TRExplorer's IGV tracks")
    parser.add_argument("--ucsc-genomic-super-dups-path", default="genomicSuperDups.txt.gz",
                        help="UCSC's hg38 genomicSuperDups table")
    parser.add_argument("--vollger2022-segdups-path", default="hg38.chr_only.SDs.bed.gz",
                        help="The hg38 segdups of Vollger et al. 2022")
    parser.add_argument("--output-dir", default=".", help="Directory to write the bgzipped BED files and their indexes to")
    args = parser.parse_args()

    for segdups_source, segdups_path in [("ucsc", args.ucsc_genomic_super_dups_path),
                                         ("vollger2022", args.vollger2022_segdups_path)]:
        bed_rows, counts_of_rows_left_out = read_segdup_bed_rows(segdups_path, segdups_source)
        output_bed_path = Path(args.output_dir) / SEGDUP_SOURCES[segdups_source]["output_bed_filename"]
        write_bgzipped_and_indexed_bed(bed_rows, output_bed_path)

        print(f"Wrote {len(bed_rows):,d} segdups from {segdups_path} to {output_bed_path}.gz")
        for reason, count in sorted(counts_of_rows_left_out.items()):
            print(f"    left out {count:,d} {reason}")
        counts_by_color = collections.Counter(bed_row[8] for bed_row in bed_rows)
        next_higher_cutoffs = [None] + [min_sequence_identity for min_sequence_identity, _ in SEQUENCE_IDENTITY_COLORS]
        for (min_sequence_identity, color), next_higher_cutoff in zip(SEQUENCE_IDENTITY_COLORS, next_higher_cutoffs):
            sequence_identity_range = (f"{min_sequence_identity} and up" if next_higher_cutoff is None
                                       else f"{min_sequence_identity} to {next_higher_cutoff}")
            print(f"    {counts_by_color[color]:10,d}  sequence identity {sequence_identity_range}  (color {color})")


if __name__ == "__main__":
    main()
