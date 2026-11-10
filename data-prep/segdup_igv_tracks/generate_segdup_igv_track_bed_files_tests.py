import gzip
import tempfile
import unittest
from pathlib import Path

from generate_segdup_igv_track_bed_files import get_color_for_sequence_identity, read_segdup_bed_rows

GRAY = "170,170,170"
YELLOW = "209,200,0"
ORANGE = "231,130,0"


def make_ucsc_segdup_line(chrom, start_0based, end_1based, strand, other_chrom, other_start_0based, other_end_1based,
                          frac_match):
    """Returns a line of UCSC's genomicSuperDups table, with filler in the columns that aren't read."""
    fields = ["585", chrom, str(start_0based), str(end_1based), f"{other_chrom}:{other_start_0based}", "0", strand,
              other_chrom, str(other_start_0based), str(other_end_1based)] + ["0"] * 16 + [frac_match] + ["0"] * 3
    return "\t".join(fields) + "\n"


def make_vollger2022_segdup_lines(rows):
    """Returns the header line and rows of a Vollger et al. 2022 SD BED file, with filler in the columns that aren't
    read.

    Args:
        rows (list): (chrom, start_0based, end_1based, other_chrom, other_start_0based, other_end_1based, strand2,
            frac_match) tuples
    """
    header = ["#chr1", "start1", "end1", "name", "fakeScore", "strand1", "start1", "end1", "color", "chr2", "start2",
              "end2", "score", "strand2"] + [f"column{i}" for i in range(14, 23)] + ["fracMatch"] + [
                 f"column{i}" for i in range(24, 44)]
    lines = ["\t".join(header) + "\n"]
    for chrom, start_0based, end_1based, other_chrom, other_start_0based, other_end_1based, strand2, frac_match in rows:
        fields = [chrom, str(start_0based), str(end_1based), "name", "0", "+", str(start_0based), str(end_1based),
                  "0,0,0", other_chrom, str(other_start_0based), str(other_end_1based), "0", strand2] + ["0"] * 9 + [
                     frac_match] + ["0"] * 20
        lines.append("\t".join(fields) + "\n")
    return lines


class GetColorForSequenceIdentityTests(unittest.TestCase):

    def test_colors_at_and_around_the_cutoffs(self):
        for sequence_identity, expected_color in [(0.9, GRAY), (0.979999, GRAY), (0.98, YELLOW), (0.989999, YELLOW),
                                                  (0.99, ORANGE), (1.0, ORANGE)]:
            self.assertEqual(get_color_for_sequence_identity(sequence_identity), expected_color, sequence_identity)


class ReadSegdupBedRowsTests(unittest.TestCase):

    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp_dir.cleanup)

    def write_gzipped_lines(self, lines):
        path = str(Path(self.temp_dir.name) / "segdups.gz")
        with gzip.open(path, "wt") as f:
            f.writelines(lines)
        return path

    def test_ucsc_rows_are_converted_sorted_and_filtered(self):
        path = self.write_gzipped_lines([
            make_ucsc_segdup_line("chr2", 500, 2000, "-", "chr1", 9, 1500, "0.985"),
            # the other copy is on an unplaced contig, which is part of the reference, so the row is kept
            make_ucsc_segdup_line("chr1", 3000, 4000, "+", "chrUn_KI270746v1", 100, 1100, "0.94308"),
            make_ucsc_segdup_line("chr1", 9, 1500, "-", "chr2", 500, 2000, "0.985"),
            # this copy is on an unplaced contig
            make_ucsc_segdup_line("chrUn_KI270746v1", 100, 1100, "+", "chr1", 3000, 4000, "0.94308"),
        ])
        bed_rows, counts_of_rows_left_out = read_segdup_bed_rows(path, "ucsc")
        self.assertEqual(bed_rows, [
            ["chr1", "9", "1500", "OtherCopy=chr2:501-2000;SequenceIdentity=0.985", "0", "-", "9", "1500", YELLOW],
            ["chr1", "3000", "4000", "OtherCopy=chrUn_KI270746v1:101-1100;SequenceIdentity=0.94308", "0", "+", "3000",
             "4000", GRAY],
            ["chr2", "500", "2000", "OtherCopy=chr1:10-1500;SequenceIdentity=0.985", "0", "-", "500", "2000", YELLOW],
        ])
        self.assertEqual(counts_of_rows_left_out, {"not on a main chromosome": 1})

    def test_vollger2022_rows_at_chry_pars_are_kept_or_left_out_as_for_the_columns(self):
        path = self.write_gzipped_lines(make_vollger2022_segdup_lines([
            # a pair whose chrY copy is inside chrY PAR1: both rows are left out
            ("chrX", 100000, 120000, "chrY", 100000, 120000, "+", "0.999"),
            ("chrY", 100000, 120000, "chrX", 100000, 120000, "+", "0.999"),
            # a pair whose chrY copy extends past the end of chrY PAR1 (2781479): the chrY row is kept, since the
            # columns use it for the chrY loci past the PAR, while the chrX row is left out
            ("chrX", 2780000, 2781700, "chrY", 2780000, 2782000, "+", "0.99995"),
            ("chrY", 2780000, 2782000, "chrX", 2780000, 2781700, "+", "0.99995"),
            # a chrY segdup that ends right where chrY PAR1 starts doesn't overlap it
            ("chrY", 5000, 10000, "chr7", 2000, 7000, "-", "0.995"),
            # an identity written with floating point noise is shown, and colored, as the value it stands for
            ("chr7", 2000, 7000, "chrY", 5000, 10000, "-", "0.9899999999999999"),
        ]))
        bed_rows, counts_of_rows_left_out = read_segdup_bed_rows(path, "vollger2022")
        self.assertEqual(bed_rows, [
            ["chr7", "2000", "7000", "OtherCopy=chrY:5001-10000;SequenceIdentity=0.99", "0", "-", "2000", "7000",
             ORANGE],
            ["chrY", "5000", "10000", "OtherCopy=chr7:2001-7000;SequenceIdentity=0.995", "0", "-", "5000", "10000",
             ORANGE],
            ["chrY", "2780000", "2782000", "OtherCopy=chrX:2780001-2781700;SequenceIdentity=0.99995", "0", "+",
             "2780000", "2782000", ORANGE],
        ])
        self.assertEqual(counts_of_rows_left_out,
                         {"whose other copy overlaps a chrY PAR": 2, "inside a chrY PAR": 1})

    def test_vollger2022_header_with_columns_in_unexpected_places_is_an_error(self):
        lines = make_vollger2022_segdup_lines([])
        path = self.write_gzipped_lines([lines[0].replace("fracMatch", "fracMatchIndel")])
        with self.assertRaisesRegex(ValueError, "header"):
            read_segdup_bed_rows(path, "vollger2022")

    def test_row_with_the_wrong_number_of_columns_is_an_error(self):
        path = self.write_gzipped_lines([make_ucsc_segdup_line("chr1", 9, 1500, "-", "chr2", 500, 2000, "0.985")])
        with self.assertRaisesRegex(ValueError, "columns"):
            read_segdup_bed_rows(path, "vollger2022")


if __name__ == "__main__":
    unittest.main()
