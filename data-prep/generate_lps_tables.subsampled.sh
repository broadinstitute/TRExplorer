#!/usr/bin/env bash
# Subsampling experiment: rebuild the per-(LocusId, Interval, VC) table for a single
# population at a range of sample counts, to see how the allele frequency histograms
# settle as N grows. Not part of the published pipeline; generate_lps_tables.sh is.
set -euo pipefail

# Anchor every path to this script's directory so it can be invoked from anywhere,
# matching generate_lps_tables.sh.
cd "$(dirname "${BASH_SOURCE[0]}")"

SCRIPT=hprc-lps/convert_multisample_LPS_table_to_allele_frequency_histograms.py
EXTRACT_SCRIPT=hprc-lps/extract_vcf_interval_metadata.py
META=hprc-lps_2026-05-19/1kGP_metadata.tsv
# The unique-TRID pair, matching generate_lps_tables.sh. The originals give a variation
# cluster the TRID of a repeat it contains, which the convert script now rejects rather
# than pairing by stream order, so pointing these back at them aborts every run below.
# See https://github.com/PacificBiosciences/trgt-lps/issues/5.
LPS_TABLE=hprc-lps_2026-05-19/hprc-lps.unique_trids.txt.gz
VCF=hprc-lps_2026-05-19/trgt-hprc.unique_trids.vcf.gz
INTERVAL_TSV=hprc-lps_2026-05-19/trgt-hprc.unique_trids.interval_metadata.tsv.gz

python3 "$EXTRACT_SCRIPT" --input-vcf "$VCF" --output-tsv "$INTERVAL_TSV"

# Each convert loads the whole interval map (~2.5 GB for the HPRC256 input) before it
# starts streaming, regardless of --num-samples, so the fan-out is bounded by RAM rather
# than by cores. MAX_JOBS=3 keeps the peak near 8 GB; running all 9 at once needs ~23 GB.
MAX_JOBS=3
population=AFR
pids=()
reaped=0
failed=0

# `wait -n` would be the natural throttle but it needs bash 4.3, and the only bash on macOS
# is 3.2. Wait on the oldest outstanding job instead: slightly cruder (it can block on a
# slow job while a later one has already finished) but it keeps at most MAX_JOBS running
# and, unlike a bare `wait`, it reads every child's exit status exactly once.
reap_one() {
    wait "${pids[$reaped]}" || failed=1
    reaped=$((reaped + 1))
}

for n in 10 20 30 40 50 60 70 80 90; do
    while (( ${#pids[@]} - reaped >= MAX_JOBS )); do reap_one; done
    python3 "$SCRIPT" --sample-metadata-tsv "$META" --input-table "$LPS_TABLE" \
        --vcf-trid-metadata-tsv "$INTERVAL_TSV" --population "$population" --num-samples "$n" &
    pids[${#pids[@]}]=$!
done

while (( reaped < ${#pids[@]} )); do reap_one; done
if (( failed )); then
    echo "ERROR: at least one convert invocation failed" >&2
    exit 1
fi
echo "[done] subsampled tables for $population at 10..90 samples"
