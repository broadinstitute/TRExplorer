import json
import unittest

from spliceai_l4_cost_and_disk_benchmark import (
    add_full_run_record_metadata, find_sai10k_problem, full_run_record, make_storage_records, thousandths)


def transcript(t_id, delta_scores, rows, priority="N"):
    ds_ag, ds_al, ds_dg, ds_dl = delta_scores
    return {"t_id": t_id, "t_priority": priority, "DS_AG": ds_ag, "DS_AL": ds_al, "DS_DG": ds_dg, "DS_DL": ds_dl,
            "DP_AG": 1, "DP_AL": 2, "DP_DG": 3, "DP_DL": 4, "ALL_NON_ZERO_SCORES": rows}


def row(pos, ra, aa, rd, ad):
    return {"pos": pos, "ref": "A", "alt": "A", "RA": ra, "AA": aa, "RD": rd, "AD": ad}


class MakeStorageRecordsTests(unittest.TestCase):

    def test_keeps_a_position_exactly_at_the_cutoff(self):
        # Observed in a production response: "0.026" - "0.016" as floats is 0.00999...
        result = {"scores": [transcript("T1", ("0.000", "0.000", "0.021", "0.000"), [
            row(100, "0.016", "0.026", "0.000", "0.000"),
            row(101, "0.016", "0.025", "0.000", "0.000"),
        ])], "allNonZeroScoresTranscriptId": "T1", "sai10kPredictions": None}
        records = make_storage_records("chr1-100-A-AT", result)
        self.assertEqual(json.loads(records["deltas_only"])["positions"], [[100, "0.010", "0.000"]])
        self.assertEqual(json.loads(records["deltas_with_ref_alt"])["positions"], [[100, "0.016", "0.026", "0.000", "0.000"]])

    def test_stores_the_transcript_the_server_selected_for_sai10k(self):
        result = {"scores": [
            transcript("mane", ("0.020", "0.000", "0.000", "0.000"), [], priority="MS"),
            transcript("higher_scoring_other", ("0.000", "0.500", "0.000", "0.000"), []),
        ], "allNonZeroScoresTranscriptId": "mane", "sai10kPredictions": {"transcript_id": "mane"}}
        record = json.loads(make_storage_records("chr1-100-A-AT", result)["deltas_only"])
        self.assertEqual(record["transcript"], "mane")
        self.assertEqual(record["DS_AG"], "0.020")
        self.assertEqual(record["sai10k"], {"transcript_id": "mane"})

    def test_nothing_stored_when_the_selected_transcript_is_below_the_cutoff(self):
        # Another transcript clears the cutoff, but it is not the one SAI-10k describes
        result = {"scores": [
            transcript("mane", ("0.009", "0.000", "0.000", "0.000"), [], priority="MS"),
            transcript("other", ("0.500", "0.000", "0.000", "0.000"), []),
        ], "allNonZeroScoresTranscriptId": "mane"}
        self.assertEqual(make_storage_records("v", result), {"deltas_only": None, "deltas_with_ref_alt": None})

    def test_flags_a_selected_transcript_missing_from_the_structure_table(self):
        # Transcripts found in the table carry EXON_STARTS from the server's structure lookup
        with_structure = {**transcript("T1", ("0.200", "0.000", "0.000", "0.000"), []), "EXON_STARTS": [1]}
        without_structure = transcript("T2", ("0.200", "0.000", "0.000", "0.000"), [])
        self.assertIsNone(find_sai10k_problem({"scores": [with_structure], "allNonZeroScoresTranscriptId": "T1"}))
        result = {"scores": [with_structure, without_structure], "allNonZeroScoresTranscriptId": "T2"}
        self.assertEqual(find_sai10k_problem(result), "T2 is not in the transcript-structure table")
        record = json.loads(make_storage_records("v", result)["deltas_with_ref_alt"])
        self.assertEqual(record["sai10k_problem"], "T2 is not in the transcript-structure table")

    def test_flags_a_sai10k_exception(self):
        result = {"scores": [{**transcript("T1", ("0.200", "0.000", "0.000", "0.000"), []), "EXON_STARTS": [1]}],
                  "allNonZeroScoresTranscriptId": "T1", "sai10kPredictions": None,
                  "sai10kPredictionsError": "Internal error computing SAI-10k predictions."}
        self.assertEqual(find_sai10k_problem(result), "Internal error computing SAI-10k predictions.")

    def test_no_sai10k_problem_field_when_complete(self):
        result = {"scores": [{**transcript("T1", ("0.200", "0.000", "0.000", "0.000"), []), "EXON_STARTS": [1]}],
                  "allNonZeroScoresTranscriptId": "T1", "sai10kPredictions": {"x": 1}}
        self.assertNotIn("sai10k_problem", json.loads(make_storage_records("v", result)["deltas_only"]))

    def test_nothing_stored_without_scores(self):
        self.assertEqual(make_storage_records("v", {"error": "no scores"}), {"deltas_only": None, "deltas_with_ref_alt": None})

    def test_thousandths_is_exact(self):
        self.assertEqual(thousandths("0.026") - thousandths("0.016"), 10)


class FullRunRecordMetadataTests(unittest.TestCase):

    def test_benchmark_records_get_the_full_runs_locus_and_target_labels(self):
        record = json.dumps({"variant": "chr1-10-C-CCA", "positions": []}, separators=(",", ":"))
        shards = [{"results_by_locus": [[
            {"variant": "chr1-10-C-CCA", "storage_records": {"deltas_only": record, "deltas_with_ref_alt": record}},
            {"variant": "chr1-10-CCA-C", "storage_records": {"deltas_only": None, "deltas_with_ref_alt": None}}]]}]
        sample = {"locus_ids": ["chr1-10-20-CA"], "alleles_by_locus": [["chr1-10-C-CCA", "chr1-10-CCA-C"]],
                  "target_labels_by_locus": [[["99.5th + 1x motif range"], ["2.5th percentile"]]]}
        add_full_run_record_metadata(shards, sample)
        stored = shards[0]["results_by_locus"][0][0]["storage_records"]["deltas_with_ref_alt"]
        self.assertEqual(stored, full_run_record("chr1-10-20-CA", ["99.5th + 1x motif range"], record))
        self.assertEqual(list(json.loads(stored))[:2], ["locus", "target_labels"])
        self.assertIsNone(shards[0]["results_by_locus"][0][1]["storage_records"]["deltas_with_ref_alt"])


if __name__ == "__main__":
    unittest.main()
