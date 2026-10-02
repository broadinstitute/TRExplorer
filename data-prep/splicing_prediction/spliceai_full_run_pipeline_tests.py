import gzip
import json
import os
import tempfile
import unittest

from make_full_run_input_chunks import allele_design_files, chunk_name, starts_a_new_chunk, write_chunk
from spliceai_full_run_pipeline import (
    RESERVED_CORES, RESERVED_MEMORY_MIB, SCORING_FILES, compute_run_id, compute_scoring_fingerprint, dollars_for_chunk,
    find_other_running_apps_of_this_pipeline, find_other_runs_with_these_inputs, find_pending_chunks, make_output_lines,
    write_atomically)
from spliceai_l4_cost_and_disk_benchmark import CORE_PRICE_PER_SECOND, GIB_PRICE_PER_SECOND, L4_PRICE_PER_SECOND


def result(variant, record):
    return {"variant": variant, "storage_records": {"deltas_only": record, "deltas_with_ref_alt": record}}


class MakeOutputLinesTests(unittest.TestCase):

    def test_keeps_chunk_order_skips_unstored_and_adds_the_locus_and_target_labels(self):
        loci = [{"locus": "chr1-10-20-CA", "alleles": ["a1", "a2"], "target_labels": [["2.5th percentile"], ["+1x", "+2x"]]},
                {"locus": "chr2-5-9-T", "alleles": ["b1"], "target_labels": [["+3x"]]}]
        # Results arrive in completion order, keyed by the locus's index in the chunk
        results_by_locus_index = {
            1: [result("b1", json.dumps({"variant": "b1", "positions": []}))],
            0: [result("a1", None), result("a2", json.dumps({"variant": "a2", "positions": [[1, "0.1", "0.2", "0.0", "0.0"]]}))],
        }
        lines = [json.loads(line) for line in make_output_lines(loci, results_by_locus_index)]
        self.assertEqual([line["variant"] for line in lines], ["a2", "b1"])
        self.assertEqual(lines[0], {"locus": "chr1-10-20-CA", "target_labels": ["+1x", "+2x"], "variant": "a2",
                                    "positions": [[1, "0.1", "0.2", "0.0", "0.0"]]})
        self.assertEqual(lines[1]["target_labels"], ["+3x"])
        self.assertEqual(list(lines[0])[0], "locus")


class FindOtherRunsWithTheseInputsTests(unittest.TestCase):

    def test_counts_only_other_runs_whose_summaries_match_these_input_chunks(self):
        manifest_chunks = [{"chunk": "chunk_00000", "sha256": "a"}, {"chunk": "chunk_00001", "sha256": "b"}]
        summaries_by_run_id = {
            "this": {"chunk_00000": {"input_sha256": "a"}},
            "earlier_same_inputs": {"chunk_00000": {"input_sha256": "a"}, "chunk_00001": {"input_sha256": "b"}},
            "earlier_other_inputs": {"chunk_00000": {"input_sha256": "old"}},
        }
        self.assertEqual(find_other_runs_with_these_inputs(summaries_by_run_id, "this", manifest_chunks),
                         {"earlier_same_inputs": 2})


class AlleleDesignFilesTests(unittest.TestCase):

    def test_every_allele_design_file_exists(self):
        for path in allele_design_files("basic"):
            self.assertTrue(os.path.exists(path), path)


class FindPendingChunksTests(unittest.TestCase):

    def test_done_only_when_the_summary_has_the_same_input_sha256_and_scoring_fingerprint(self):
        manifest_chunks = [{"chunk": f"chunk_0000{i}", "sha256": sha256} for i, sha256 in enumerate(["aaa", "bbb", "ccc", "ddd", "eee"])]
        summaries = {
            "chunk_00000": {"input_sha256": "aaa", "scoring_fingerprint": "v50"},
            "chunk_00001": {"input_sha256": "old", "scoring_fingerprint": "v50"},
            "chunk_00003": {"input_sha256": "ddd", "scoring_fingerprint": "v49"},
            "chunk_00004": {"input_sha256": "eee"},  # written before summaries had a fingerprint
        }
        self.assertEqual([c["chunk"] for c in find_pending_chunks(manifest_chunks, summaries, "v50")],
                         ["chunk_00001", "chunk_00002", "chunk_00003", "chunk_00004"])


class ComputeScoringFingerprintTests(unittest.TestCase):

    def test_changes_with_file_contents_and_settings_only(self):
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "table.tsv")
            write_atomically(path, b"ENST1\t+\n")
            first = compute_scoring_fingerprint([path], {"distance": 10000})
            self.assertEqual(compute_scoring_fingerprint([path], {"distance": 10000}), first)
            self.assertNotEqual(compute_scoring_fingerprint([path], {"distance": 5000}), first)
            write_atomically(path, b"ENST1\t-\n")
            self.assertNotEqual(compute_scoring_fingerprint([path], {"distance": 10000}), first)

    def test_run_id_depends_only_on_the_scoring_fingerprint(self):
        fingerprint = compute_scoring_fingerprint([], {"distance": 10000})
        run_id = compute_run_id(fingerprint)
        self.assertEqual(len(run_id), 16)
        self.assertNotEqual(compute_run_id(compute_scoring_fingerprint([], {"distance": 500})), run_id)

    def test_regenerated_inputs_rescore_only_changed_chunks(self):
        # Same folder (same fingerprint); the summary's input SHA-256 decides which chunks are done
        summaries = {"chunk_00000": {"input_sha256": "aaa", "scoring_fingerprint": "f"},
                     "chunk_00001": {"input_sha256": "bbb", "scoring_fingerprint": "f"}}
        regenerated = [{"chunk": "chunk_00000", "sha256": "aaa"}, {"chunk": "chunk_00001", "sha256": "ccc"}]
        self.assertEqual([c["chunk"] for c in find_pending_chunks(regenerated, summaries, "f")], ["chunk_00001"])

    def test_every_scoring_file_exists(self):
        self.assertEqual([str(p) for p in SCORING_FILES if not p.exists()], [])


class FindOtherRunningAppsTests(unittest.TestCase):

    def test_ignores_this_app_stopped_apps_and_other_pipelines(self):
        app_list = [
            {"app_id": "ap-this", "description": "spliceai-tr-full-run", "state": "ephemeral"},
            {"app_id": "ap-old-running", "description": "spliceai-tr-full-run", "state": "ephemeral (detached)"},
            {"app_id": "ap-old-stopped", "description": "spliceai-tr-full-run", "state": "stopped"},
            {"app_id": "ap-other", "description": "spliceai-l4-cost-and-disk-benchmark", "state": "ephemeral"},
        ]
        self.assertEqual([a["app_id"] for a in find_other_running_apps_of_this_pipeline(app_list, "ap-this")],
                         ["ap-old-running"])


class DollarsForChunkTests(unittest.TestCase):

    def test_bills_the_higher_of_reserved_and_used_cores_and_memory(self):
        under = {"container_seconds": 100, "cores_used_while_busy": 2.0, "memory_gib": 21.0}
        self.assertAlmostEqual(dollars_for_chunk(under), 100 * (
            L4_PRICE_PER_SECOND + RESERVED_CORES * CORE_PRICE_PER_SECOND + RESERVED_MEMORY_MIB / 1024 * GIB_PRICE_PER_SECOND))
        over = {"container_seconds": 100, "cores_used_while_busy": 6.0, "memory_gib": 40.0}
        self.assertAlmostEqual(dollars_for_chunk(over), 100 * (
            L4_PRICE_PER_SECOND + 6.0 * CORE_PRICE_PER_SECOND + 40.0 * GIB_PRICE_PER_SECOND))


class FileWritingTests(unittest.TestCase):

    def test_write_atomically_leaves_no_temporary_file(self):
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "out.json")
            write_atomically(path, b"{}")
            self.assertEqual(os.listdir(d), ["out.json"])

    def test_identical_chunks_give_identical_checksums(self):
        loci = [{"locus": "chr1-10-20-CA", "alleles": ["chr1-10-C-CCA"]}]
        with tempfile.TemporaryDirectory() as d:
            first = write_chunk(d, "basic", loci)
            second = write_chunk(d, "basic", loci)
            with gzip.open(os.path.join(d, f"{chunk_name('chr1-10-20-CA')}.json.gz"), "rt") as f:
                content = json.load(f)
        self.assertEqual(first, second)
        self.assertEqual(first["n_alleles"], 1)
        self.assertTrue(first["chunk"].startswith("chunk_chr1_000000010_"))
        self.assertEqual(content, {"gene_set": "basic", "chunk": chunk_name("chr1-10-20-CA"), "loci": loci})


class ChunkBoundaryTests(unittest.TestCase):

    def chunk(self, locus_ids, mean_loci_per_chunk):
        chunks = []
        for locus_id in locus_ids:
            if not chunks or starts_a_new_chunk(locus_id, len(chunks[-1]), mean_loci_per_chunk):
                chunks.append([])
            chunks[-1].append(locus_id)
        return chunks

    def test_removing_a_locus_changes_only_its_own_chunk(self):
        locus_ids = [f"chr1-{i * 100}-{i * 100 + 20}-CA" for i in range(1, 3000)]
        before = self.chunk(locus_ids, 50)
        after = self.chunk(locus_ids[:10] + locus_ids[11:], 50)
        changed = [c for c in after if c not in before]
        self.assertEqual(len(changed), 1)
        self.assertGreater(len(before), 20)

    def test_removing_a_locus_inside_a_capped_chunk_changes_only_its_own_chunk(self):
        locus_ids = [f"chr1-{i * 100}-{i * 100 + 20}-CA" for i in range(1, 40000)]
        before = self.chunk(locus_ids, 200)
        capped = [i for i, c in enumerate(before) if len(c) > 600]
        self.assertTrue(capped)
        removed = locus_ids.index(before[capped[0]][300])
        after = self.chunk(locus_ids[:removed] + locus_ids[removed + 1:], 200)
        self.assertEqual(len([c for c in after if c not in before]), 1)

    def test_capped_chunks_end_soon_after_three_times_the_mean(self):
        sizes = [len(c) for c in self.chunk([f"chr1-{i}-{i + 5}-A" for i in range(100000)], 200)]
        self.assertLess(max(sizes), 600 + 200)


if __name__ == "__main__":
    unittest.main()
