# SpliceAI scores for simulated tandem repeat alleles

This folder scores simulated contractions and expansions at TRExplorer v2.1 loci with SpliceAI and
SAI-10k, using the production SpliceAI-lookup code on Modal L4 GPUs. It then summarizes the scores into
four `SpliceAI_*` columns of the TRExplorer BigQuery table. The goal is to show, for each locus, whether
biologically plausible repeat-length changes are predicted to affect splicing.

## Which loci and alleles

- **Loci:** every TRExplorer v2.1 locus that lies inside a GENCODE v50 basic transcript, which is where
  SpliceAI-lookup gives scores. The locus must also be polymorphic in HPRC256, meaning it has at least 2
  distinct total allele lengths among the 256 HPRC samples' haplotypes. These lengths come from the TRGT
  VCF's single-locus records; variation-cluster records are skipped. With the current inputs this is
  2,070,141 loci.
- **Alleles per locus:** up to 6 target total lengths:
  - `2.5pct`, `97.5pct` and `99.5pct`: those percentiles of the HPRC256 allele lengths.
  - `+1x`, `+2x` and `+3x`: the 99.5th percentile plus 1, 2 or 3 times the motif's range. The motif's
    range is the 75th percentile, over HPRC256 loci with the same canonical motif and at least 5 distinct
    lengths, of each locus's 97.5th minus 2.5th percentile range. A motif with fewer than 100 such loci
    uses the value for all motifs of its length.
- **Size changes:** each target becomes a whole number of repeat units added to or removed from the hg38
  tract, with exact halves rounded away from zero.
  - Targets that round to the same change are scored once, and a change of 0 is dropped.
  - Insertions and deletions are capped at 5,000 bp, and a contraction always keeps at least one unit.
  - Loci with a base other than A, C, G or T in the anchor base or the tract are skipped.
- **Scoring settings:** distance 10,000, no masking, GENCODE v50 basic.

With the current inputs this gives 8,999,426 alleles in 2,151 chunks. The estimated cost is about $455,
with a 90% range of $445 to $465. That is about 366 L4 GPU-hours, which takes about 37 hours at
Modal's limit of 10 GPUs at a time.

## Scripts

**Pipeline steps, in the order they run:**

| Step | Script | Writes |
|---|---|---|
| 1 | `count_catalog_loci_near_gencode_splice_sites.py` | `TRExplorer_v2.1_loci_vs_GENCODE_v50_splice_sites.distances.tsv.gz` (the catalog's loci, with their distances to splice sites) |
| 2 | `compute_hprc256_total_allele_length_stats.py` | `HPRC256_total_allele_length_stats.tsv.gz` (per-locus HPRC256 allele length histograms) |
| 3 | `compute_hprc256_range_percentiles_by_motif.py` | `HPRC256_range_percentile_by_motif.tsv` (each motif's range) |
| 4 | `make_transcript_structures_table.py` | `transcripts_hg38.gencode_v50_comprehensive.tsv` (the transcript table SAI-10k reads) |
| 5 | `make_full_run_input_chunks.py` | `full_run_inputs/<gene_set>/` (allele chunks and `manifest.json`) |
| 6 | `spliceai_full_run_pipeline.py` | scores on the Modal Volume `spliceai-tr-full-run` |
| 7 | `summarize_full_run_spliceai_scores_per_locus.py` | `TRExplorer_SpliceAI_summary.<gene_set>.tsv.gz` (one row per locus) |

**Modules that are both a library and a benchmark.** Two scripts have "benchmark" in their name, but the
pipeline imports code from them:

- `make_tr_expansion_benchmark_alleles.py` defines the allele design. Steps 3, 5 and 7 import its
  constants, file paths and allele-building functions. Its own `main` only writes a random sample of
  loci for the benchmark, `benchmark/tr_benchmark_alleles.<gene_set>.json`
  (`python3 make_tr_expansion_benchmark_alleles.py --gene-set basic --n-loci 1000`).
- `spliceai_l4_cost_and_disk_benchmark.py` defines the scoring setup used by steps 5 and 6: the Modal
  image, the GPU settings, the per-locus scoring code and the output record format. Its own `main` is
  the benchmark (`modal run spliceai_l4_cost_and_disk_benchmark.py --gene-sets basic`). That run scores
  the sample and writes `benchmark/l4_benchmark_results.<gene_set>.json` and
  `benchmark/l4_benchmark_estimate.json`, the cost and disk estimates above. It refuses a sample made
  under an older allele design unless you pass `--allow-stale-sample`.

**Benchmark and design analysis only**, in `benchmark/`. These are not needed to run the pipeline:

- `analyze_population_based_design.py` reads the benchmark sample and results and writes
  `population_based_design_report.html`, comparing candidate sets of target sizes.
- `compare_simulated_sizes_to_pathogenic_thresholds.py` checks how often the simulated sizes reach
  STRchive pathogenic thresholds. `+3x` reaches the threshold at 53 of 63 comparable loci. It reads
  `~/code/STRchive/data/STRchive-loci.json`.
- `compare_simulated_sizes_to_pathogenic_thresholds_using_lps.py` does the same check using longest
  pure segment (LPS) lengths instead of total lengths. It does worse, reaching 29 of the 62 loci it can compare.

They import the pipeline modules from this folder and can be run from any directory, e.g.
`python3 benchmark/analyze_population_based_design.py`.

Each `*_tests.py` file runs with `python3 -m unittest <module>_tests` from this folder.

## Setup

Everything except the Modal steps runs locally. Data files are gitignored, so a fresh checkout has
only the code.

1. **Python packages:** `pip3 install modal numpy pyfaidx`. Install `markdown` and `matplotlib` too if you want to run
   `analyze_population_based_design.py`. `bcftools` must be on your `PATH` for step 2.
2. **Modal:** create an account and run `modal token new`. The pipeline runs on L4 GPUs, at most 10 at
   a time (`MAX_CONCURRENT_CONTAINERS` in `spliceai_full_run_pipeline.py`).
3. **Reference genome:** `~/hg38.fa` with its `.fai` index.
4. **SpliceAI-lookup checkout** at `~/code/SpliceAI-lookup`. The pipeline reads these files from it:
   - `google_cloud_run_services/docker/ref/GRCh38/hg38.fa.gz` and its index files, which are uploaded
     to the Modal image (`REFERENCE_DIR` in `spliceai_l4_cost_and_disk_benchmark.py`).
   - `google_cloud_run_services/gencode.v50.GRCh38.comprehensive.sorted.txt.gz`, the input to step 4.
5. **Production code snapshot.** Copy the production scoring code and annotation files into
   `spliceai_lookup_files_gencode_v50/`:
   ```bash
   S=~/code/SpliceAI-lookup/google_cloud_run_services
   mkdir -p spliceai_lookup_files_gencode_v50
   cp $S/server.py $S/sai10k_predictions.py \
      $S/docker/spliceai/annotations/GRCh38/gencode.v50.basic.annotation.txt.gz \
      $S/docker/spliceai/annotations/GRCh38/gencode.v50.annotation.txt.gz \
      $S/docker/ref/GRCh38/gencode.v50.basic.annotation.transcript_annotations.json.gz \
      $S/docker/ref/GRCh38/gencode.v50.annotation.transcript_annotations.json.gz \
      spliceai_lookup_files_gencode_v50/
   ```
   The run's scoring fingerprint covers these files, so the run folder name changes whenever they
   change (see "Reruns" below).
6. **HPRC256 TRGT VCF:** `../hprc-lps_2026-05-19/trgt-hprc.unique_trids.vcf.gz`, the 256 HPRC samples
   genotyped with the TRExplorer catalog (`VCF` in `compute_hprc256_total_allele_length_stats.py`).
7. **Catalog and GENCODE downloads, into this folder:**
   ```bash
   curl -LO https://github.com/broadinstitute/tandem-repeat-catalogs/releases/download/v2.1/TRExplorer.repeat_catalog_v2.1.hg38.1_to_1000bp_motifs.bed.gz
   curl -LO https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.basic.annotation.gtf.gz
   curl -LO https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.annotation.gtf.gz
   ```

## Running the pipeline

Run every command from this folder.

```bash
# 1. Distance from every catalog locus to the nearest splice site (about 1 min, 1.7 GB of memory)
python3 count_catalog_loci_near_gencode_splice_sites.py \
    --catalog-bed TRExplorer.repeat_catalog_v2.1.hg38.1_to_1000bp_motifs.bed.gz \
    --gtf gencode.v50.basic.annotation.gtf.gz --gtf gencode.v50.annotation.gtf.gz \
    --output-prefix TRExplorer_v2.1_loci_vs_GENCODE_v50_splice_sites

# 2-4. HPRC256 allele lengths, motif ranges, and the SAI-10k transcript table
python3 compute_hprc256_total_allele_length_stats.py --n-processes 8
python3 compute_hprc256_range_percentiles_by_motif.py
python3 make_transcript_structures_table.py

# 5. Input chunks (about 1,000 loci each; 2.7 GB of memory)
python3 make_full_run_input_chunks.py --gene-set basic

# 6. Score on Modal
modal run spliceai_full_run_pipeline.py::main --gene-set basic                                # dry run: what would run
modal run spliceai_full_run_pipeline.py::main --gene-set basic --max-chunks 2 --launch        # small test
modal run --detach spliceai_full_run_pipeline.py::main --gene-set basic --launch              # full run
modal run spliceai_full_run_pipeline.py::report --gene-set basic                              # progress, cost, run folder

# 7. Download the run folder that main and report print, then summarize it per locus
modal volume get spliceai-tr-full-run /outputs/basic/<run_id> ./full_run_outputs/basic
python3 summarize_full_run_spliceai_scores_per_locus.py --gene-set basic --outputs-dir full_run_outputs/basic
```

`--detach` keeps the run going if this terminal disconnects. All chunks are submitted up front, so Modal
holds the queue.

### Reruns

- **Interrupted runs.** Launching again with the same inputs and code skips the chunks that are already
  done and scores the rest. `main` refuses to launch while an earlier app of this pipeline is still
  running. Stop that app with `modal app stop <app_id>`, or pass `--force`.
- **Allele design changes.** `manifest.json` records a fingerprint of the design files, and the pipeline
  refuses inputs built under an older design. Chunk boundaries are set by the loci themselves, so after
  a design change `make_full_run_input_chunks.py` rewrites only the chunks whose loci changed, and only
  those are scored again.
- **Scoring changes.** The run folder name (`run_id`) is the first 16 characters of a fingerprint of
  the scoring code and data: the production snapshot, the transcript table, the reference files and the
  two pipeline `.py` files. Any edit to these, even a comment, starts a new run folder. If another
  folder already holds scores for the same inputs, `main` refuses to launch. Pass
  `--rescore-already-scored-chunks` only if the scores really should be redone, since that pays for
  every chunk again.
- **Which outputs are current.** The local `full_run_inputs/<gene_set>/manifest.json` decides: only the
  chunks listed there, with a matching input SHA-256 in their summary, are current. `report` lists the
  rest, and the summary script reads only the current chunks, refusing to run if any is missing. So keep
  the inputs a run was launched with until it has been summarized. Each launch also writes
  `current_chunks.json` to the run folder, a record of the chunks it launched; no script reads it.

## Outputs

**Modal Volume** `spliceai-tr-full-run`, folder `/outputs/<gene_set>/<run_id>/`:

- `<chunk>.jsonl.gz`: one JSON line per allele whose selected transcript has a delta score of at least
  0.01. The selected transcript is the one SpliceAI-lookup shows: MANE Select, then MANE Plus Clinical,
  then canonical, then the largest sum of delta scores. Each line has:
  - the locus ID and the allele's target labels;
  - the transcript ID, and the four delta scores with their positions;
  - the SAI-10k predictions, plus `sai10k_problem` when they are missing or degraded;
  - every position where the acceptor or donor delta is at least 0.01, stored as
    `[pos, REF acceptor, ALT acceptor, REF donor, ALT donor]`.

  Alleles below 0.01 are not stored, and the summary step counts them as 0. Alleles whose scoring
  returned an error (listed in the chunk's summary) are left out of the summary instead. The estimate is 0.35 GB
  gzipped for the full run.
- `<chunk>.summary.json`: counts, errors, alleles that failed twice, SAI-10k problems, timing, the input
  chunk's SHA-256 and the scoring fingerprint.
- `current_chunks.json`: the chunks of the most recent launch, with their input SHA-256s.

**Per-locus summary** `TRExplorer_SpliceAI_summary.<gene_set>.tsv.gz`, one row per scored locus:

| Column | Meaning |
|---|---|
| `LocusId` | TRExplorer locus ID without the `chr` prefix, e.g. `1-12345-12400-CAG` |
| `SpliceAI_MaxDeltaScore` | largest delta score (acceptor or donor, gain or loss, on the selected transcript) of any simulated allele; 0 if none reached 0.01 |
| `SpliceAI_MaxDeltaScoreAlleleSize` | the size that gave it: `2.5pct`, `97.5pct`, `99.5pct`, `+1x`, `+2x` or `+3x`; empty if none reached 0.01 |
| `SpliceAI_MinAlleleSizeThatAffectsSplicing` | the smallest size with a delta score of at least 0.2, with the same values; empty if none |
| `SpliceAI_DeltaScoreByRepeatCount` | every simulated allele, sorted by repeat count, as `repeat_count:max_delta_score` plus the type of the largest change when it is at least 0.01 (`AG` acceptor gain, `AL` acceptor loss, `DG` donor gain, `DL` donor loss). Example: `18:0.000,25:0.030AG,40:0.310AG,95:0.880AL` |

The sizes, in order of increasing length, are 2.5pct, 97.5pct, 99.5pct, +1x, +2x, +3x. When several sizes
round to the same allele, the allele is named by the smallest. Loci that were not scored are left out of
the summary, so their columns are NULL in BigQuery. These are loci that are not polymorphic in HPRC256,
or that lie outside every GENCODE v50 basic transcript.

## Loading into BigQuery

The column definitions are in `../../bigquery-proxy/global_constants.py`, as `SPLICEAI_COLUMN_NAMES`
in the "SpliceAI (simulated alleles)" group. There are two ways to load them, and both are run from
`../../bigquery-proxy/`.

**A. Fill the live catalog table in place**, with no rebuild:

```bash
cd ../../bigquery-proxy
TSV=../data-prep/splicing_prediction/TRExplorer_SpliceAI_summary.basic.tsv.gz
python3 add_spliceai_columns_to_catalog_table.py $TSV --check-only     # only reports how many LocusIds match
python3 add_spliceai_columns_to_catalog_table.py $TSV                  # fills the most recent catalog_* table
python3 add_spliceai_columns_to_catalog_table.py $TSV --table-id catalog_YYYYMMDD_HHMMSS   # or a specific one
```

This runs these steps:

1. Adds any missing `SpliceAI_*` columns to the table's schema.
2. Loads the TSV into a temporary staging table.
3. Sets every row's `SpliceAI_*` columns to NULL, then copies in the values by `LocusId`.
4. Deletes the staging table.

Loci missing from the TSV end up NULL, even if an earlier fill set them.

**B. Include the columns when rebuilding the table**, by adding this option to the usual
`load_bigquery_main_table.py` command (see the `load` target in `../../bigquery-proxy/Makefile`):

```bash
--spliceai-summary-tsv ../data-prep/splicing_prediction/TRExplorer_SpliceAI_summary.basic.tsv.gz
```

Either way, regenerate the website afterwards so its filters and export dialog list the new columns:

```bash
cd ../website && python3 generate_website.py
```
