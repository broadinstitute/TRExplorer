#!/usr/bin/env bash
# Convert the extended-definitions BED written by run_extend_v2_catalog_boundaries.sh into an
# ExpansionHunter catalog, and upload it to the bucket run_expansion_hunter_on_selected_samples.py
# reads from.
#
# LocusIds are each definition's own "{chrom}-{start_0based}-{end_1based}-{motif}", which is what
# convert_bed_to_expansion_hunter_catalog.py assigns. An extended definition therefore carries no id
# linking it to the locus it grew from; the two are matched by overlap when the results are analyzed.
#
# --trim is deliberately NOT passed. It would shorten every locus whose span is not a whole number of
# motif copies, which is about 62% of this catalog, and that is a boundary change of exactly the kind
# this experiment is measuring. Leaving it off keeps the original definitions byte-identical to the
# published v2 catalog's.
#
# Run 2026-08-24: 5,657,854 definitions, no duplicate locus ids.
set -u

CATALOG_DIR=~/code/tandem-repeat-explorer/data-prep/flank_metric_from_andrea/boundary_optimization/data
INPUT_BED=${CATALOG_DIR}/TRExplorer.repeat_catalog_v2.hg38.1_to_1000bp_motifs.with_extended_definitions.bed.gz
OUTPUT_JSON=${CATALOG_DIR}/TR_catalog.TRExplorer-v2.with_extended_definitions.5657854_loci.ExpansionHunter.json.gz
CLOUD_DIR=gs://tandem-repeat-catalog/v2.0

python3 -m str_analysis.convert_bed_to_expansion_hunter_catalog \
  --show-progress-bar \
  --output-path "$OUTPUT_JSON" \
  "$INPUT_BED" 2>&1 | tee "${CATALOG_DIR}/convert_extended_catalog_to_eh.log"

gsutil cp "$OUTPUT_JSON" "${CLOUD_DIR}/"
