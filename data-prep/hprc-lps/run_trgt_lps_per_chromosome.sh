#!/usr/bin/env bash
# Run trgt-lps one chromosome at a time and concatenate the results.
#
# trgt-lps 0.11.0 aborts intermittently on a whole-genome run with a double-free inside
# its rayon worker threads (macOS reports
# ___BUG_IN_CLIENT_OF_LIBMALLOC_POINTER_BEING_FREED_WAS_NOT_ALLOCATED while dropping a
# trgt_lps::vcf::Locus). It is a race, not a data problem: the same input runs clean on
# smaller slices. Going chromosome by chromosome keeps each job short enough to retry
# cheaply, and only one slice is on disk at a time.

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."

BATCH_DIR=hprc-lps_2026-05-19
VCF=$BATCH_DIR/trgt-hprc.unique_trids.vcf.gz
OUT=$BATCH_DIR/hprc-lps.unique_trids.txt.gz
WORK=$BATCH_DIR/lps_per_chrom
# 12 of the machine's 14 cores. Dropping to 4 was tried and did not stop the run being killed from
# outside the script, so the lighter footprint bought nothing while making each chromosome roughly
# three times slower. The retry loop below covers the rayon race that more threads make likelier.
THREADS=12
MAX_ATTEMPTS=4

CHROMS=$(for i in $(seq 1 22); do echo chr$i; done; echo chrX; echo chrY; echo chrM)

for required in "$VCF" "$VCF.tbi"; do
    if [[ ! -f "$required" ]]; then
        echo "ERROR: required input not found: $required" >&2
        exit 1
    fi
done

mkdir -p "$WORK"

# Take an exclusive lock on the work directory. Every invocation writes the same fingerprint, the
# same chromosome slices and the same temporary tables, so two overlapping runs (the ordinary way
# this script gets resumed after a failure) could truncate a slice another trgt-lps is reading, or
# concatenate a table the other is still writing. mkdir is atomic on every filesystem we use.
LOCK=$WORK/.lock
if ! mkdir "$LOCK" 2>/dev/null; then
    echo "ERROR: $LOCK exists, so another run of this script is already working in $WORK." >&2
    echo "       Wait for it to finish, or remove $LOCK if no such process is running." >&2
    exit 1
fi
trap 'rmdir "$LOCK" 2>/dev/null' EXIT

# Cheap identity for what produced the cached tables: the input VCF's byte size, a checksum of its
# tabix index (which changes whenever any record position does), and the trgt-lps version. Hashing
# the 16 GB VCF itself would cost minutes. Without this, a re-run after the VCF is regenerated, or
# after trgt-lps is upgraded, would silently concatenate tables computed from two different inputs
# or two different versions of the tool into one published file.
FINGERPRINT_FILE=$WORK/input.fingerprint
fingerprint="$(wc -c < "$VCF" | tr -d ' ') $(shasum -a 256 "$VCF.tbi" | cut -d' ' -f1) $(trgt-lps --version)"
# A missing fingerprint file counts as a mismatch: tables left by a run from before this check
# existed cannot be vouched for either, and adopting them is the same silent mixing.
if [[ "$(cat "$FINGERPRINT_FILE" 2>/dev/null)" != "$fingerprint" ]]; then
    if compgen -G "$WORK/*.txt" > /dev/null; then
        echo "=== $(date '+%H:%M:%S')  cached chromosome tables do not match $VCF; discarding them"
        echo "    so the output cannot mix results computed from two different VCFs"
        rm -f "$WORK"/*.txt
    fi
    # $OUT was concatenated from the tables just discarded, so it describes the previous
    # input. Move it aside now rather than at the end: a chromosome that exhausts its
    # retries, or an interrupted run, exits long before the concatenation below, and would
    # otherwise leave a table from the old VCF sitting at the canonical path next to the
    # new one, looking current.
    if [[ -e "$OUT" ]]; then
        echo "=== $(date '+%H:%M:%S')  $OUT was built from the previous input; moving it to $OUT.stale"
        mv "$OUT" "$OUT.stale"
    fi
fi
printf '%s\n' "$fingerprint" > "$FINGERPRINT_FILE"

for chrom in $CHROMS; do
    out_txt=$WORK/$chrom.txt
    if [[ -s "$out_txt" ]]; then
        echo "=== $(date '+%H:%M:%S')  $chrom: already done ($(wc -l < "$out_txt" | tr -d ' ') rows), skipping"
        continue
    fi

    slice=$WORK/$chrom.vcf.gz
    echo "=== $(date '+%H:%M:%S')  $chrom: slicing"
    { tabix -H "$VCF"; tabix "$VCF" "$chrom"; } | bgzip -@ $THREADS > "$slice"
    tabix -f -p vcf "$slice"

    ok=0
    for attempt in $(seq 1 $MAX_ATTEMPTS); do
        echo "=== $(date '+%H:%M:%S')  $chrom: trgt-lps attempt $attempt"
        if trgt-lps --vcf "$slice" --threads $THREADS > "$out_txt.tmp"; then
            mv "$out_txt.tmp" "$out_txt"
            ok=1
            echo "    $chrom: $(wc -l < "$out_txt" | tr -d ' ') rows"
            break
        fi
        echo "    $chrom: attempt $attempt failed"
        rm -f "$out_txt.tmp"
    done
    rm -f "$slice" "$slice.tbi"
    if [[ $ok -ne 1 ]]; then
        echo "ERROR: $chrom failed $MAX_ATTEMPTS times" >&2
        exit 1
    fi
done

echo "=== $(date '+%H:%M:%S')  concatenating"
# Keep the header from the first chromosome only; every run emits the same one.
{
    first=1
    for chrom in $CHROMS; do
        if [[ $first -eq 1 ]]; then
            cat "$WORK/$chrom.txt"
            first=0
        else
            tail -n +2 "$WORK/$chrom.txt"
        fi
    done
} | bgzip -@ $THREADS > "$OUT.tmp"
mv "$OUT.tmp" "$OUT"
ls -l "$OUT"
echo "=== $(date '+%H:%M:%S')  done ==="
