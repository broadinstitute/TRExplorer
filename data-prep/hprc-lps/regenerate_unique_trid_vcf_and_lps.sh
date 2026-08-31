#!/usr/bin/env bash
# Rebuild the published HPRC256 TRGT VCF and LPS table so that variation clusters carry
# their own TRIDs, fixing https://github.com/PacificBiosciences/trgt-lps/issues/5.
#
#   1. rewrite_vc_trids_in_vcf.py  swaps TRID and STRUC on variation cluster records
#      (TRID=VC:{chrom}:{start}-{end}, STRUC=<VC:{old trid}>), matching the v2.1 catalog.
#   2. tabix indexes the rewritten VCF.
#   3. run_trgt_lps_per_chromosome.sh re-derives the LPS table from it, so its TRIDs are
#      unique by construction and no row needs to be paired back to a VCF record after the
#      fact. That driver exists because a single whole-genome trgt-lps run aborts
#      intermittently; see its header for details.
#
# Nothing is re-genotyped: only the TRID and STRUC INFO fields change.

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."

BATCH_DIR=hprc-lps_2026-05-19
IN_VCF=$BATCH_DIR/trgt-hprc.vcf.gz
OUT_VCF=$BATCH_DIR/trgt-hprc.unique_trids.vcf.gz
THREADS=12

# Leave no partial artifact behind if a step dies mid-write.
trap 'rm -f "$OUT_VCF.tmp"' EXIT

echo "=== $(date '+%H:%M:%S')  step 1/3: rewriting variation cluster TRIDs ==="
gzcat "$IN_VCF" \
    | python3 hprc-lps/rewrite_vc_trids_in_vcf.py \
    | bgzip -@ $THREADS > "$OUT_VCF.tmp"
mv "$OUT_VCF.tmp" "$OUT_VCF"
ls -l "$OUT_VCF"

echo "=== $(date '+%H:%M:%S')  step 2/3: indexing ==="
tabix -f -p vcf "$OUT_VCF"

echo "=== $(date '+%H:%M:%S')  verifying the rewrite before the expensive step ==="
# Confirm on a few chromosomes that nothing but TRID and STRUC changed. set -e aborts here rather
# than spending many hours re-deriving the LPS table from a VCF that is already wrong.
python3 hprc-lps/verify_trid_rewrite.py \
    --original-vcf "$IN_VCF" --rewritten-vcf "$OUT_VCF" chr1 chr19 chrX chrM

echo "=== $(date '+%H:%M:%S')  step 3/3: re-running trgt-lps, one chromosome at a time ==="
bash hprc-lps/run_trgt_lps_per_chromosome.sh

echo "=== $(date '+%H:%M:%S')  done ==="
