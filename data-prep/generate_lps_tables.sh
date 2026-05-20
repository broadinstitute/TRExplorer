#!/usr/bin/env bash
# Regenerate every per-(LocusId, Interval, VC) table needed by the HPRC256
# BigQuery loaders, starting from the merged TRGT multisample VCF + the
# pre-computed trgt-lps wide-format TSV.
#
# Pipeline:
#   1. extract_vcf_interval_metadata.py
#        VCF -> small (trid, locus_id, motif, interval, vc) TSV used by step 4.
#   2. compute_allele_size_purity_and_methylation_distributions_from_vcf.py
#        VCF -> stratified allele-size-vs-purity + allele-size-vs-methylation
#        TSVs (with --stratify-by-population --stratify-by-sex). After writing,
#        we symlink the outputs into the legacy hprc-lps/ directory so the
#        purity/methylation BQ loader's defaults still work without --flag overrides.
#   3. run_decompose_hprc_alleles.py
#        Hail Batch dispatcher: VCF -> trviz decomposed-alleles parquet. The
#        orchestrator passes --force so worker-script edits always reach Batch.
#   4. convert_multisample_LPS_table_to_allele_frequency_histograms.py
#        SINGLE invocation with --stratify-by-population --stratify-by-sex.
#        Produces the all-strata file that load_bigquery_HPRC256_stratified_LPS_frequency_data.py
#        requires; also symlinked into hprc-lps/ for the loader defaults.

set -euo pipefail

# Anchor every path to this script's directory so it can be invoked from
# anywhere (e.g. `bash data-prep/generate_lps_tables.sh` from the repo root).
cd "$(dirname "${BASH_SOURCE[0]}")"

SCRIPT=hprc-lps/convert_multisample_LPS_table_to_allele_frequency_histograms.py
EXTRACT_SCRIPT=hprc-lps/extract_vcf_interval_metadata.py
PURITY_METH_SCRIPT=hprc-lps/compute_allele_size_purity_and_methylation_distributions_from_vcf.py
DECOMPOSE_SCRIPT=run_decompose_hprc_alleles.py
BATCH_DIR=hprc-lps_2026-05-19
LEGACY_DIR=hprc-lps
META=$BATCH_DIR/1kGP_metadata.tsv
LPS_TABLE=$BATCH_DIR/hprc-lps.txt.gz
VCF=$BATCH_DIR/trgt-hprc.vcf.gz
INTERVAL_TSV=$BATCH_DIR/trgt-hprc.interval_metadata.tsv.gz

for required in "$META" "$LPS_TABLE" "$VCF" "$VCF.tbi" \
                "$SCRIPT" "$EXTRACT_SCRIPT" "$PURITY_METH_SCRIPT" "$DECOMPOSE_SCRIPT"; do
    if [[ ! -e "$required" ]]; then
        echo "ERROR: required input file not found: $required" >&2
        exit 1
    fi
done

echo "============================================================"
echo "Step 1/4: extract_vcf_interval_metadata.py"
echo "============================================================"
python3 "$EXTRACT_SCRIPT" --input-vcf "$VCF" --output-tsv "$INTERVAL_TSV"
echo "[done] step 1: $INTERVAL_TSV"

echo
echo "============================================================"
echo "Step 2/4: compute_allele_size_purity_and_methylation_distributions_from_vcf.py"
echo "============================================================"
python3 "$PURITY_METH_SCRIPT" \
    --input-vcf "$VCF" \
    --sample-metadata-tsv "$META" \
    --stratify-by-population \
    --stratify-by-sex

# The compute script writes outputs next to --input-vcf (in $BATCH_DIR). The
# BQ purity/methylation loader's defaults look in $LEGACY_DIR. Symlink the
# produced files into $LEGACY_DIR so the loader works without --flag overrides.
N_SAMPLES=256
PURITY_OUT="$BATCH_DIR/trgt-hprc.allele_size_purity.stratified.by_population.by_sex.${N_SAMPLES}_samples.tsv.gz"
METH_OUT="$BATCH_DIR/trgt-hprc.methylation.stratified.by_population.by_sex.${N_SAMPLES}_samples.tsv.gz"
mkdir -p "$LEGACY_DIR"
ln -sf "../$PURITY_OUT" "$LEGACY_DIR/trgt-hprc.allele_size_purity.stratified.by_population.by_sex.${N_SAMPLES}_samples.tsv.gz"
ln -sf "../$METH_OUT"   "$LEGACY_DIR/trgt-hprc.methylation.stratified.by_population.by_sex.${N_SAMPLES}_samples.tsv.gz"
echo "[done] step 2: purity + methylation stratified TSVs (symlinked into $LEGACY_DIR for BQ loader defaults)"

echo
echo "============================================================"
echo "Step 3/4: run_decompose_hprc_alleles.py  (Hail Batch dispatch, with --force)"
echo "============================================================"
echo "Dispatching the decompose pipeline to Hail Batch. The orchestrator passes"
echo "--force so any local edits to data-prep/decompose_hprc_alleles.py reach"
echo "the Batch workers (otherwise a stale GCS-cached copy may be reused)."
python3 "$DECOMPOSE_SCRIPT" --force
echo "[done] step 3: decompose pipeline submitted"

echo
echo "============================================================"
echo "Step 4/4: convert_multisample_LPS_table_to_allele_frequency_histograms.py  (single stratified run)"
echo "============================================================"
# Single invocation with --stratify-by-population --stratify-by-sex. This is
# the exact file that load_bigquery_HPRC256_stratified_LPS_frequency_data.py
# requires (named .per_locus_and_motif.by_population.by_sex.<N>_samples.tsv.gz).
# Replaces the prior 18-way fan-out, which loaded the full ~1.5 GB VCF
# metadata map per process (peak ~18-30 GB on a workstation) and did NOT
# produce the file the loader expects anyway.
python3 "$SCRIPT" \
    --sample-metadata-tsv "$META" \
    --input-table "$LPS_TABLE" \
    --vcf-interval-tsv "$INTERVAL_TSV" \
    --stratify-by-population \
    --stratify-by-sex

# Symlink the produced file into $LEGACY_DIR so the LPS BQ loader's default
# --input-tsv path resolves.
CONVERT_OUT_BASENAME="hprc-lps.per_locus_and_motif.by_population.by_sex.${N_SAMPLES}_samples.tsv.gz"
CONVERT_OUT="$BATCH_DIR/$CONVERT_OUT_BASENAME"
if [[ ! -e "$CONVERT_OUT" ]]; then
    echo "WARNING: expected convert output not found at $CONVERT_OUT" >&2
else
    ln -sf "../$CONVERT_OUT" "$LEGACY_DIR/$CONVERT_OUT_BASENAME"
fi
echo "[done] step 4: stratified LPS TSV ready (symlinked into $LEGACY_DIR for BQ loader defaults)"

echo
echo "============================================================"
echo "Local steps complete: extract + purity/methylation + convert finished."
echo
echo "IMPORTANT: step 3 (decompose) was dispatched to Hail Batch and is still"
echo "running asynchronously. Do NOT run load_bigquery_trviz_allele_decompositions.py"
echo "until the Batch pipeline reports completion (monitor in the Batch web UI)."
echo "============================================================"
