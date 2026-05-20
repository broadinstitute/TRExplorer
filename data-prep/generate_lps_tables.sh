#!/usr/bin/env bash
# Regenerate all the per-(LocusId, Interval, VC) tables used by the HPRC256
# BigQuery loaders, starting from the merged TRGT multisample VCF + the
# pre-computed trgt-lps wide-format TSV.
#
# Pipeline:
#   1. extract_vcf_interval_metadata.py
#        VCF -> small (trid, locus_id, motif, interval, vc) TSV used by step 4.
#   2. compute_allele_size_purity_and_methylation_distributions_from_vcf.py
#        VCF -> stratified allele-size-vs-purity + allele-size-vs-methylation
#        TSVs consumed by load_bigquery_HPRC256_stratified_allele_purity_and_methylation_data.py.
#   3. run_decompose_hprc_alleles.py
#        Hail Batch dispatcher: VCF -> trviz decomposed-alleles parquet
#        consumed by load_bigquery_trviz_allele_decompositions.py.
#   4. convert_multisample_LPS_table_to_allele_frequency_histograms.py
#        18 parallel invocations: trgt-lps TSV + interval-metadata TSV ->
#        stratified per-(LocusId, Interval, VC) allele frequency histograms
#        consumed by load_bigquery_HPRC256_stratified_LPS_frequency_data.py.

set -euo pipefail

SCRIPT=hprc-lps/convert_multisample_LPS_table_to_allele_frequency_histograms.py
EXTRACT_SCRIPT=hprc-lps/extract_vcf_interval_metadata.py
PURITY_METH_SCRIPT=hprc-lps/compute_allele_size_purity_and_methylation_distributions_from_vcf.py
DECOMPOSE_SCRIPT=run_decompose_hprc_alleles.py
META=hprc-lps_2026-05-19/1kGP_metadata.tsv
LPS_TABLE=hprc-lps_2026-05-19/hprc-lps.txt.gz
VCF=hprc-lps_2026-05-19/trgt-hprc.vcf.gz
INTERVAL_TSV=hprc-lps_2026-05-19/trgt-hprc.interval_metadata.tsv.gz

for required in "$META" "$LPS_TABLE" "$VCF"; do
    if [[ ! -f "$required" ]]; then
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
echo "[done] step 2: purity + methylation stratified TSVs"

echo
echo "============================================================"
echo "Step 3/4: run_decompose_hprc_alleles.py  (Hail Batch dispatch)"
echo "============================================================"
echo "Dispatching the decompose pipeline to Hail Batch. Monitor progress in"
echo "the Hail Batch web UI; this step returns once the orchestrator has"
echo "submitted the worker jobs."
python3 "$DECOMPOSE_SCRIPT"
echo "[done] step 3: decompose pipeline submitted"

echo
echo "============================================================"
echo "Step 4/4: convert_multisample_LPS_table_to_allele_frequency_histograms.py (18 parallel)"
echo "============================================================"
pids=()
labels=()

start_convert() {
    local label="$1"; shift
    python3 "$SCRIPT" \
        --sample-metadata-tsv "$META" \
        --input-table "$LPS_TABLE" \
        --vcf-interval-tsv "$INTERVAL_TSV" \
        "$@" &
    pids+=($!)
    labels+=("$label")
}

start_convert "all"
for s in male female; do
    start_convert "sex=$s" --sex "$s"
done
for p in AFR AMR EAS EUR SAS; do
    start_convert "pop=$p" --population "$p"
    for s in male female; do
        start_convert "pop=$p,sex=$s" --population "$p" --sex "$s"
    done
done

failed=()
for i in "${!pids[@]}"; do
    if ! wait "${pids[$i]}"; then
        failed+=("${labels[$i]} (pid=${pids[$i]})")
    fi
done

if [[ ${#failed[@]} -gt 0 ]]; then
    echo "ERROR: ${#failed[@]} convert invocation(s) failed:" >&2
    for f in "${failed[@]}"; do
        echo "    - $f" >&2
    done
    exit 1
fi
echo "[done] step 4: all ${#pids[@]} convert invocations completed"

echo
echo "============================================================"
echo "All steps completed successfully."
echo "============================================================"
