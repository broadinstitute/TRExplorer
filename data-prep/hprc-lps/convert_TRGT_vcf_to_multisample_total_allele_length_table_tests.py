import unittest

from convert_TRGT_vcf_to_multisample_total_allele_length_table import (
    convert_record_to_rows, convert_sample_field_to_total_lengths_in_repeats)


class ConvertSampleFieldToTotalLengthsInRepeatsTests(unittest.TestCase):

    def test_diploid_calls_keep_genotype_order(self):
        self.assertEqual(convert_sample_field_to_total_lengths_in_repeats("0/1:51,63", 3), "17,21")
        self.assertEqual(convert_sample_field_to_total_lengths_in_repeats("23|12:63,51", 3), "21,17")

    def test_lengths_are_rounded_down_to_whole_repeats(self):
        # CAGCAGCA (8 bp) is 2 whole CAG repeats
        self.assertEqual(convert_sample_field_to_total_lengths_in_repeats("0/0:8,8", 3), "2,2")
        self.assertEqual(convert_sample_field_to_total_lengths_in_repeats("1/1:2,0", 3), "0,0")

    def test_hemizygous_call_has_one_value(self):
        self.assertEqual(convert_sample_field_to_total_lengths_in_repeats("0:9", 2), "4")

    def test_no_call_is_a_single_dot(self):
        self.assertEqual(convert_sample_field_to_total_lengths_in_repeats("./.:.", 3), ".")
        self.assertEqual(convert_sample_field_to_total_lengths_in_repeats(".:.", 3), ".")

    def test_unmeasured_allele_of_a_partial_call_is_a_dot(self):
        self.assertEqual(convert_sample_field_to_total_lengths_in_repeats("./1:.,30", 3), ".,10")
        self.assertEqual(convert_sample_field_to_total_lengths_in_repeats("0/1:51,.", 3), "17,.")


class ConvertRecordToRowsTests(unittest.TestCase):

    LOCUS_IDS_IN_VARIATION_CLUSTERS = {"4-3074876-3074933-CAG", "4-3074927-3074936-CAA"}

    def test_single_locus_record_of_an_isolated_locus(self):
        rows = convert_record_to_rows("1-100-112-AT", "AT", "<TR:1-100-112-AT>", ["0/1:12,14", "./.:."],
                                      self.LOCUS_IDS_IN_VARIATION_CLUSTERS)
        self.assertEqual(rows, [["1-100-112-AT", "AT", "6,7", "."]])

    def test_single_locus_record_of_a_locus_in_a_variation_cluster_is_skipped(self):
        rows = convert_record_to_rows("4-3074876-3074933-CAG", "CAG", "<TR:4-3074876-3074933-CAG>", ["0/1:51,63"],
                                      self.LOCUS_IDS_IN_VARIATION_CLUSTERS)
        self.assertEqual(rows, [])

    def test_variation_cluster_record_gets_one_row_per_motif(self):
        rows = convert_record_to_rows("VC:4:3074876-3074946", "CAG,CAAT",
                                      "<VC:4-3074876-3074933-CAG,4-3074927-3074936-CAAT>", ["36|19:64,76"],
                                      self.LOCUS_IDS_IN_VARIATION_CLUSTERS)
        self.assertEqual(rows, [["VC:4:3074876-3074946", "CAG", "21,25"], ["VC:4:3074876-3074946", "CAAT", "16,19"]])


if __name__ == "__main__":
    unittest.main()
