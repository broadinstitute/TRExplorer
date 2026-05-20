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
#        Two invocations sharing the same VCF interval map + LPS table:
#          4a) --stratify-by-population --stratify-by-sex -> the all-strata file
#              required by load_bigquery_HPRC256_stratified_LPS_frequency_data.py
#          4b) no stratification -> the per_locus_and_motif file required by
#              load_bigquery_main_table.py via --hprc256-tsv.
#        Both outputs are symlinked into hprc-lps/ for the loader defaults.

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
# rm -f first so we never leave a stale real file from a prior run shadowing
# the canonical symlink path. ln -sf alone replaces an existing entry, but
# being explicit makes the "hprc-lps/ holds only current-snapshot symlinks"
# convention obvious.
PURITY_LINK="$LEGACY_DIR/trgt-hprc.allele_size_purity.stratified.by_population.by_sex.${N_SAMPLES}_samples.tsv.gz"
METH_LINK="$LEGACY_DIR/trgt-hprc.methylation.stratified.by_population.by_sex.${N_SAMPLES}_samples.tsv.gz"
rm -f "$PURITY_LINK" "$METH_LINK"
ln -s "../$PURITY_OUT" "$PURITY_LINK"
ln -s "../$METH_OUT"   "$METH_LINK"
echo "[done] step 2: purity + methylation stratified TSVs (symlinked into $LEGACY_DIR for BQ loader defaults)"

echo
echo "============================================================"
echo "Step 3/4: run_decompose_hprc_alleles.py  (Hail Batch dispatch)"
echo "============================================================"
echo "Dispatching the decompose pipeline to Hail Batch. Each invocation stages"
echo "shard outputs + worker script under a fresh timestamped subdir of"
echo "gs://...decomposed_alleles/, so concurrent runs can't clobber each other."
echo "Worker script is re-uploaded unconditionally so local edits always reach"
echo "the cluster. To resume a prior interrupted run, pass --run-id <timestamp>."
python3 "$DECOMPOSE_SCRIPT"
echo "[done] step 3: decompose pipeline submitted"

echo
echo "============================================================"
echo "Step 4/4: convert_multisample_LPS_table_to_allele_frequency_histograms.py"
echo "============================================================"
# Two sequential invocations of the convert script against the same VCF
# interval map + LPS table:
#   4a) --stratify-by-population --stratify-by-sex: required by
#       load_bigquery_HPRC256_stratified_LPS_frequency_data.py
#       (.per_locus_and_motif.by_population.by_sex.<N>_samples.tsv.gz)
#   4b) no stratification: required by load_bigquery_main_table.py via
#       --hprc256-tsv (.per_locus_and_motif.<N>_samples.tsv.gz)
# Two sequential invocations is still much cheaper than the prior 18-way
# fan-out (peak ~18-30 GB on a workstation), and guarantees the stratified
# and non-stratified files are derived from the same source snapshot.
echo "  4a) stratified by population x sex"
python3 "$SCRIPT" \
    --sample-metadata-tsv "$META" \
    --input-table "$LPS_TABLE" \
    --vcf-interval-tsv "$INTERVAL_TSV" \
    --stratify-by-population \
    --stratify-by-sex

echo "  4b) non-stratified"
python3 "$SCRIPT" \
    --sample-metadata-tsv "$META" \
    --input-table "$LPS_TABLE" \
    --vcf-interval-tsv "$INTERVAL_TSV"

# Symlink both produced files into $LEGACY_DIR so the BQ loaders' default
# --input paths resolve. rm -f first so any stale real file from a prior run
# (legacy hprc-lps/ contents predate this script's symlink convention) does
# not shadow the canonical path.
STRATIFIED_BASENAME="hprc-lps.per_locus_and_motif.by_population.by_sex.${N_SAMPLES}_samples.tsv.gz"
STRATIFIED_OUT="$BATCH_DIR/$STRATIFIED_BASENAME"
UNSTRATIFIED_BASENAME="hprc-lps.per_locus_and_motif.${N_SAMPLES}_samples.tsv.gz"
UNSTRATIFIED_OUT="$BATCH_DIR/$UNSTRATIFIED_BASENAME"
for f in "$STRATIFIED_OUT" "$UNSTRATIFIED_OUT"; do
    if [[ ! -e "$f" ]]; then
        echo "WARNING: expected convert output not found at $f" >&2
    fi
done
if [[ -e "$STRATIFIED_OUT" ]]; then
    rm -f "$LEGACY_DIR/$STRATIFIED_BASENAME"
    ln -s "../$STRATIFIED_OUT" "$LEGACY_DIR/$STRATIFIED_BASENAME"
fi
if [[ -e "$UNSTRATIFIED_OUT" ]]; then
    rm -f "$LEGACY_DIR/$UNSTRATIFIED_BASENAME"
    ln -s "../$UNSTRATIFIED_OUT" "$LEGACY_DIR/$UNSTRATIFIED_BASENAME"
fi
echo "[done] step 4: stratified + non-stratified LPS TSVs ready (symlinked into $LEGACY_DIR for BQ loader defaults)"

echo
echo "============================================================"
echo "Local steps complete: extract + purity/methylation + convert finished."
echo
echo "IMPORTANT: step 3 (decompose) was dispatched to Hail Batch and is still"
echo "running asynchronously. Do NOT run load_bigquery_trviz_allele_decompositions.py"
echo "until the Batch pipeline reports completion (monitor in the Batch web UI)."
echo "============================================================"
