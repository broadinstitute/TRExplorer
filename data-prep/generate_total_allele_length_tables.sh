#!/usr/bin/env bash
# Build the HPRC256 total allele length tables (the HPRC256TotalLength_* columns and their population x sex
# table) from the same unique-TRID TRGT VCF that generate_lps_tables.sh uses.
#
# Pipeline:
#   1. hprc-lps/convert_TRGT_vcf_to_multisample_total_allele_length_table.py
#        VCF -> wide (trid, motif, <sample>...) table of total allele lengths in repeat units (FORMAT/AL //
#        motif length): a variation cluster's record for loci in a cluster, the single-locus record for the
#        rest. Also writes the matching subset of the interval metadata TSV.
#   2. convert_multisample_LPS_table_to_allele_frequency_histograms.py, the same converter as for the LPS
#      table, run twice in parallel on that table:
#        2a) --stratify-by-population --stratify-by-sex -> the file loaded into the population x sex table
#        2b) no stratification -> the file load_bigquery_main_table.py reads via --hprc256-total-length-tsv
#      Both outputs are symlinked into hprc-lps/ under stable names for the loaders' defaults.
#
# Step 1 reads the interval metadata TSV written by step 1 of generate_lps_tables.sh, so run that first if the
# VCF changed.

set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")"

BATCH_DIR=hprc-lps_2026-05-19
LEGACY_DIR=hprc-lps
META=$BATCH_DIR/1kGP_metadata.tsv
VCF=$BATCH_DIR/trgt-hprc.unique_trids.vcf.gz
INTERVAL_TSV=$BATCH_DIR/trgt-hprc.unique_trids.interval_metadata.tsv.gz
TOTAL_LENGTH_TABLE=$BATCH_DIR/hprc-total-allele-length.unique_trids.txt.gz
TOTAL_LENGTH_INTERVAL_TSV=$BATCH_DIR/hprc-total-allele-length.unique_trids.interval_metadata.tsv.gz
N_SAMPLES=256

for required in "$META" "$VCF" "$VCF.tbi" "$INTERVAL_TSV"; do
    if [[ ! -e "$required" ]]; then
        echo "ERROR: required input file not found: $required" >&2
        exit 1
    fi
done

echo "Step 1/2: convert_TRGT_vcf_to_multisample_total_allele_length_table.py"
python3 hprc-lps/convert_TRGT_vcf_to_multisample_total_allele_length_table.py \
    --input-vcf "$VCF" \
    --input-interval-metadata-tsv "$INTERVAL_TSV" \
    --output-table "$TOTAL_LENGTH_TABLE" \
    --output-interval-metadata-tsv "$TOTAL_LENGTH_INTERVAL_TSV"

echo "Step 2/2: convert_multisample_LPS_table_to_allele_frequency_histograms.py (stratified and not, in parallel)"
# The converter writes each output next to its input table, i.e. into $BATCH_DIR.
python3 -m str_analysis.convert_multisample_LPS_table_to_allele_frequency_histograms \
    --sample-metadata-tsv "$META" \
    --input-table "$TOTAL_LENGTH_TABLE" \
    --vcf-trid-metadata-tsv "$TOTAL_LENGTH_INTERVAL_TSV" \
    --stratify-by-population \
    --stratify-by-sex &
STRATIFIED_PID=$!
python3 -m str_analysis.convert_multisample_LPS_table_to_allele_frequency_histograms \
    --sample-metadata-tsv "$META" \
    --input-table "$TOTAL_LENGTH_TABLE" \
    --vcf-trid-metadata-tsv "$TOTAL_LENGTH_INTERVAL_TSV" &
UNSTRATIFIED_PID=$!
wait "$STRATIFIED_PID"
wait "$UNSTRATIFIED_PID"

TOTAL_LENGTH_STEM=$(basename "$TOTAL_LENGTH_TABLE" .txt.gz)
for kind in "per_locus_and_motif.by_population.by_sex" "per_locus_and_motif"; do
    produced="$BATCH_DIR/${TOTAL_LENGTH_STEM}.${kind}.${N_SAMPLES}_samples.tsv.gz"
    link="$LEGACY_DIR/hprc-total-allele-length.${kind}.${N_SAMPLES}_samples.tsv.gz"
    if [[ ! -e "$produced" ]]; then
        echo "ERROR: expected converter output not found at $produced" >&2
        exit 1
    fi
    rm -f "$link"
    ln -s "../$produced" "$link"
    echo "Linked $link -> ../$produced"
done
echo "Done."
