"""Boundary-extension rule (v2): tolerate a short run of imperfect motif units in the flank rather
than requiring an immediate exact match, as long as the whole "newly added sequence" -- the
imperfect run PLUS the pure repeat run immediately following it -- is itself highly pure.

Per side, independently, repeat:
  1. Scan outward from the current boundary in whole-motif-unit steps, looking for the next EXACT
     copy of the motif, allowing at most `max_repeats_in_gap` consecutive non-exact units before it.
     If no exact unit is found within that budget, stop -- this side is done.
  2. From that exact unit, consume the entire contiguous run of further exact units that follows.
  3. "Newly added sequence" = everything from the current boundary through the end of that pure run
     (i.e. the imperfect gap PLUS the pure run after it). Compute its purity (fraction of bases
     matching an in-frame perfect tiled repeat of the motif).
  4. If purity >= `min_purity_of_new_sequence`, accept: move the boundary to the end of the pure run
     and go back to step 1 (try to extend further). Otherwise, reject and stop -- the boundary stays
     at its pre-step-1 position.

Example: boundary sits at the end of a perfect AT repeat. Flank = "CTCT" (2 imperfect AT-frame units)
followed by "ATATATATAT" (5 pure AT units). max_repeats_in_gap=2 lets the scan reach the first exact
unit inside "ATATATATAT". Newly added sequence = "CTCTATATATATAT" (14bp, 2 mismatched bases from the
CTCT gap) -- purity = 12/14 = 0.857. With the default min_purity_of_new_sequence=0.9 this specific
example would be REJECTED (below threshold); a less degenerate gap or a longer pure run would pass.

Reuses boundary_metrics.cumulative_mismatched_bases for the per-position phase/match arrays rather
than re-deriving the reversed-coordinate math for the left side by hand.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from boundary_metrics import cumulative_mismatched_bases

MAX_REPEATS_IN_GAP = 2
MIN_PURITY_OF_NEW_SEQUENCE = 0.9
MAX_OFFSET = 300
MAX_RESCAN_ATTEMPTS = 5  # each attempt triples the fetch window, so the last is 81x the first. The
# first window is motif-dependent (see scan_side) and at least MAX_OFFSET, which for a motif of 100bp
# or less gives 300bp, 900bp, 2700bp, 8100bp, 24300bp. Five is one attempt beyond the offset where a
# real locus's extension was observed to stabilize, so convergence (two consecutive equal values) has
# a chance to be confirmed rather than merely reached on the last attempt.


def scan_side(fasta, chrom, side, start_0based, end_1based, motif, core_length):
    """Extend one side, tripling the fetched flank window on each retry until the returned extension
    stops changing. A truncated fetch window can cut short
    either the anchor search or the pure-repeat run that follows an anchor inside
    extend_by_gap_purity, in either case producing a value that looks like a genuine,
    rule-determined stop but isn't -- the only way to tell the two apart from outside that function
    is to check whether a bigger window changes the answer. Returns (extension,
    still_capped_after_max_attempts).

    Lives here rather than in either scan script because both of them need it and they must agree:
    when the two produced their own copies, one was fixed to test for convergence and the other was
    left testing whether the extension merely landed short of the fetch window, which accepts a
    truncated answer as final.

    The first window is at least the anchor search's full reach, (max_repeats_in_gap + 1) motif
    copies. A fixed 300bp start is smaller than that for motifs over 100bp, and two windows too small
    to finish the search can agree on a value neither of them was in a position to determine, which
    convergence alone cannot tell apart from a genuine stop.
    """
    offset = max(MAX_OFFSET, (MAX_REPEATS_IN_GAP + 1) * len(motif))
    prev_ext = None
    for _ in range(MAX_RESCAN_ATTEMPTS):
        if side == "left":
            flank = fasta.fetch(chrom, max(0, start_0based - offset), start_0based)[::-1]
        else:
            flank = fasta.fetch(chrom, end_1based, end_1based + offset)
        ext = extend_by_gap_purity(flank, motif, core_length, side)
        if ext == prev_ext:
            return ext, False
        prev_ext = ext
        offset *= 3
    return prev_ext, True


def extend_by_gap_purity(
    flank_seq,
    motif,
    core_length,
    side,
    max_repeats_in_gap=MAX_REPEATS_IN_GAP,
    min_purity_of_new_sequence=MIN_PURITY_OF_NEW_SEQUENCE,
    verbose=False,
):
    """flank_seq: forward genomic order for 'right', reversed (boundary-adjacent base first) for
    'left', which is the same convention pytrf_style_extension.py uses. Returns the accepted
    extension length.
    """
    motif = motif.upper()
    flank_seq = flank_seq.upper()
    motif_length = len(motif)
    max_offset = len(flank_seq)

    mismatched = cumulative_mismatched_bases(flank_seq, motif, core_length, side)

    def mismatches_in(start, end):
        if start >= end:
            return 0
        return int(mismatched[end - 1]) - (int(mismatched[start - 1]) if start > 0 else 0)

    def is_exact_unit(unit_start):
        unit_end = unit_start + motif_length
        return unit_end <= max_offset and mismatches_in(unit_start, unit_end) == 0

    accepted = 0
    while True:
        pos = accepted
        gap_units = 0
        anchor_start = None
        while gap_units <= max_repeats_in_gap:
            if pos + motif_length > max_offset:
                break
            if is_exact_unit(pos):
                anchor_start = pos
                break
            pos += motif_length
            gap_units += 1

        if anchor_start is None:
            if verbose:
                print(f"    no exact anchor within {max_repeats_in_gap} imperfect unit(s) -- stop at {accepted}")
            return accepted

        run_end = anchor_start
        while is_exact_unit(run_end):
            run_end += motif_length

        n_bases = run_end - accepted
        n_mismatches = mismatches_in(accepted, run_end)
        purity = 1 - n_mismatches / n_bases

        if purity >= min_purity_of_new_sequence:
            if verbose:
                n_pure_units = (run_end - anchor_start) // motif_length
                print(
                    f"    accept: newly-added={flank_seq[accepted:run_end]!r} ({n_bases}bp = "
                    f"{gap_units} imperfect unit(s) + {n_pure_units} exact unit(s)), "
                    f"purity={purity:.3f} >= {min_purity_of_new_sequence} -- extend to {run_end}"
                )
            accepted = run_end
            continue

        if verbose:
            print(
                f"    reject: newly-added={flank_seq[accepted:run_end]!r} ({n_bases}bp), "
                f"purity={purity:.3f} < {min_purity_of_new_sequence} -- stop at {accepted}"
            )
        return accepted


def compute_candidate_definitions(
    chrom, start, end, motif, fasta,
    max_repeats_in_gap=MAX_REPEATS_IN_GAP,
    min_purity_of_new_sequence=MIN_PURITY_OF_NEW_SEQUENCE,
):
    """Given a locus's current [start, end) (0-based half-open) and motif, apply the gap-purity
    rule independently on each side and return the 4 candidate definitions (dedup'd): baseline,
    left-only-extended, right-only-extended, both-extended.
    """
    core_seq = fasta.fetch(chrom, start, end)
    core_length = len(core_seq)

    left_flank = fasta.fetch(chrom, max(0, start - MAX_OFFSET), start)[::-1]
    right_flank = fasta.fetch(chrom, end, end + MAX_OFFSET)

    left_extension = extend_by_gap_purity(left_flank, motif, core_length, "left", max_repeats_in_gap, min_purity_of_new_sequence)
    right_extension = extend_by_gap_purity(right_flank, motif, core_length, "right", max_repeats_in_gap, min_purity_of_new_sequence)

    definitions = {
        "baseline": (start, end),
        "left_only": (start - left_extension, end),
        "right_only": (start, end + right_extension),
        "both": (start - left_extension, end + right_extension),
    }
    seen = {}
    for label, region in definitions.items():
        seen.setdefault(region, []).append(label)

    return {
        "left_extension": left_extension,
        "right_extension": right_extension,
        "definitions": definitions,
        "unique_regions": seen,
    }


LOCI = [
    ("EP400", "chr12", 132062548, 132062611, "CAG", "starting from narrower def, per gnomAD blog"),
    ("C9ORF72", "chr9", 27573528, 27573546, "GGCCCC", "KnownDiseaseAssociatedLoci current def"),
    ("RFC1", "chr4", 39348424, 39348479, "AARRG", "KnownDiseaseAssociatedLoci current def"),
    ("BEAN1", "chr16", 66490398, 66490453, "TAAAA", "KnownDiseaseAssociatedLoci current def"),
    ("RUNX2", "chr6", 45422750, 45422792, "GCN", "KnownDiseaseAssociatedLoci current def"),
    ("AT_locus_chr1_210718331", "chr1", 210718331, 210718351, "AT", "flank-repetitiveness round4 example"),
]


def main():
    import pysam

    fasta = pysam.FastaFile("/Users/weisburd/hg38.fa")

    for locus_id, chrom, start, end, motif, note in LOCI:
        print(f"=== {locus_id}  {chrom}:{start}-{end}  motif={motif}  ({note}) ===")
        result = compute_candidate_definitions(chrom, start, end, motif, fasta)
        print(f"  left_extension={result['left_extension']}bp  right_extension={result['right_extension']}bp")
        for region, labels in result["unique_regions"].items():
            new_start, new_end = region
            n_units = (new_end - new_start) // len(motif)
            print(f"  {chrom}:{new_start}-{new_end}  ({n_units}x{motif})  <- {', '.join(labels)}")
        print()


if __name__ == "__main__":
    main()
