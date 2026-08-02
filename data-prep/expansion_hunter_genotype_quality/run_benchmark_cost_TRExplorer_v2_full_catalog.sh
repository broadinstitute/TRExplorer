#!/usr/bin/env bash
# Cost benchmark: EHv5-bw2-optimized genotyping the REAL full TRExplorer-v2 catalog (5,594,988 loci,
# gs://tandem-repeat-catalog/v2.0/, downloaded 2026-08-02) vs. a 2% pilot subset (111,900 loci, every
# 50th locus genome-wide) -- unlike the prior 2026-06-25 benchmark (393k-811k loci, motif<=9bp/span<=120bp
# only), this is the ACTUAL unrestricted catalog composition (motifs up to ~1000bp, spans up to several kb;
# many such loci get cheaply skipped by EH's read-length heuristics rather than fully aligned -- confirmed
# by a local dry run -- so per-locus cost here may not match the small-motif-only benchmark's rate).
# EH genotyping step ONLY (--skip-*-step). UNSHARDED (EHV5_NUM_SHARDS=1). highmem. Samples: same 3 WGS
# replicates as the prior benchmark (HG00438, NA12878, HG00733), illumina ~31x, all confirmed PCR-free
# (NYGC 1000G 30x). Configs requested: cpu1/t2, cpu2/t4, cpu4/t8. Idempotent (skips a combo whose log
# already has "Submitted batch"). Can be run from anywhere -- cd's into run_tools/ itself.
set -u

RUN_TOOLS_DIR=~/code/str-truth-set-v2/run_tools
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

CATALOG_BASE=gs://tandem-repeat-catalog/v2.0
OUTPUT_DIR=gs://str-truth-set-v2/benchmark_cost__2026_08_01/tool_results
LOGDIR="$SCRIPT_DIR/launch_logs"
mkdir -p "$LOGDIR"

# "catalog_filename|label"
CATALOGS=(
  "TR_catalog.TRExplorer-v2.pilot_every_50th.ExpansionHunter.json.gz|pilot_111900"
  "TR_catalog.TRExplorer-v2.5594988_loci.ExpansionHunter.json.gz|full_5594988"
)
# "cpu threads"
CONFIGS=("1 2" "2 4" "4 8")
SAMPLES="HG00438 NA12878 HG00733"

launch () {
  cat_file="$1"; label="$2"; cpu="$3"; threads="$4"
  subdir="${label}/cpu${cpu}_t${threads}"
  log="${LOGDIR}/${label}.cpu${cpu}_t${threads}.log"
  if [ -f "$log" ] && grep -q "Submitted batch" "$log"; then
    echo ">>> SKIP ${subdir} (already submitted: $(grep -oE 'https://batch[^ ]+' "$log" | head -1))"
    return
  fi
  sample_args=""
  for s in $SAMPLES; do sample_args="$sample_args --sample-id $s"; done
  echo ">>> LAUNCH ${subdir} [illumina: ${SAMPLES}] -> $log"
  (cd "$RUN_TOOLS_DIR" && \
  EHV5_STREAMING_CPU=$cpu EHV5_STREAMING_THREADS=$threads EHV5_STREAMING_MEMORY=highmem EHV5_NUM_SHARDS=1 \
  python3 run_genotyping_tools.py \
    --no-wait \
    --tool EHv5-bw2-optimized \
    --data-type illumina \
    --custom-catalog-path "${CATALOG_BASE}/${cat_file}" \
    --output-dir "$OUTPUT_DIR" \
    --output-subdir "$subdir" \
    --skip-combine-expansion-hunter-step \
    --skip-add-columns-step \
    --skip-plot-accuracy-step \
    $sample_args \
  ) > "$log" 2>&1
  echo "    exit=$? batch: $(grep -oE 'https://batch[^ ]+' "$log" | head -1)"
}

for entry in "${CATALOGS[@]}"; do
  cat_file="${entry%%|*}"; label="${entry##*|}"
  for cfg in "${CONFIGS[@]}"; do
    set -- $cfg
    launch "$cat_file" "$label" "$1" "$2"
  done
done
echo "ALL DONE"
