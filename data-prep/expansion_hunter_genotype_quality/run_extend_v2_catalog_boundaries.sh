#!/usr/bin/env bash
# Extend the locus boundaries in the TRExplorer v2 catalog outward into their flanking sequence using
# the gap-purity rule, producing the catalog that ExpansionHunter is run against.
#
# The rule lives in str_analysis/extend_str_catalog_boundaries.py (github.com/broadinstitute/str-analysis,
# commit b73cbe3), which is installed in editable mode, so `python3 -m` resolves to that working copy
# from any directory.
#
# -k keeps each extended locus's ORIGINAL definition in the output alongside its extended definition,
# making the output a superset of the input. That is what genotyping both definitions in a single
# ExpansionHunter run and comparing them requires; without it each extended locus's original
# definition is replaced and there is nothing to compare against.
#
# Thresholds are the script's own defaults, spelled out here so a later change to those defaults does
# not silently change what this produces:
#   --max-repeats-in-gap 2           how many consecutive interrupted motif copies the scan may cross
#   --min-purity-of-new-sequence 0.9 fraction of the ADDED bases that must match a perfect motif tiling
#   --min-reference-purity 0.9       a locus is only considered if its existing definition is this pure
#
# Run 2026-08-24: ~6 minutes for all 5,599,658 loci. Output is bgzipped and tabix-indexed.
set -u

REFERENCE_FASTA=~/hg38.fa
CATALOG=~/code/tandem-repeat-catalogs/results__2026-02-01/release_draft_2026-02-01/TRExplorer.repeat_catalog_v2.hg38.1_to_1000bp_motifs.bed.gz
OUTPUT_DIR=~/code/tandem-repeat-explorer/data-prep/flank_metric_from_andrea/boundary_optimization/data
OUTPUT=${OUTPUT_DIR}/TRExplorer.repeat_catalog_v2.hg38.1_to_1000bp_motifs.with_extended_definitions.bed.gz

python3 -m str_analysis.extend_str_catalog_boundaries \
  --reference-fasta "$REFERENCE_FASTA" \
  --keep-original-definitions-of-extended-loci \
  --max-repeats-in-gap 2 \
  --min-purity-of-new-sequence 0.9 \
  --min-reference-purity 0.9 \
  --show-progress-bar \
  --output-path "$OUTPUT" \
  "$CATALOG" 2>&1 | tee "${OUTPUT_DIR}/extend_v2_catalog.log"
