"""Alternative boundary-extension rule: scan outward in whole-motif-unit steps looking for the next
EXACT copy of the motif. Accept a jump out to that copy (encompassing the gap in between) if the gap's
hamming distance to a perfect in-frame repeat is <= 1/3 of the gap length (gap length is always a
multiple of the motif length). Keep scanning as long as the cumulative gap fraction stays <= 1/3;
stop (no further extension) the moment it exceeds 1/3 without having found an anchoring exact copy.

This differs from boundary_metrics.py's purity-threshold curves: here the trigger for stopping is
finding a genuinely bad, unanchored gap, not the position purity crosses a fixed value -- so it can
tolerate an interruption (like EP400's CAA/CAA) as long as clean motif copies resume soon enough.

Reuses boundary_metrics.cumulative_purity_and_hamming for the actual per-position phase/match logic
(already validated against hand-computed EP400 numbers) rather than re-deriving it, since getting the
left-side reversed-coordinate phase math right by hand is easy to get subtly wrong (first draft of
this script did).

CAVEAT: compute_candidate_definitions() below only ever fetches MAX_OFFSET (300) bases of flank, so an
extension that reaches exactly 300bp can't be distinguished from a genuine rule-determined stop -- it
may just mean the fetched window ran out before the rule did. Callers that treat a returned extension
as final (e.g. threshold_sweep_trexplorer_v1.py) should treat values at/near 300bp as a lower bound,
not a confirmed stopping point. This module implements the exploratory motif-reanchoring rule that was
superseded by gap_purity_extension.py for the currently-deployed catalog.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from boundary_metrics import cumulative_purity_and_hamming

HAMMING_FRACTION_THRESHOLD = 1 / 3
MAX_OFFSET = 300


def extend_by_motif_reanchoring(
    flank_seq, motif, core_length, side, threshold=HAMMING_FRACTION_THRESHOLD, min_anchor_units=1, verbose=False
):
    """flank_seq: forward genomic order for 'right', reversed (boundary-adjacent base first) for 'left'
    -- same convention as boundary_metrics.py. Returns the accepted extension length in bp.

    threshold: max allowed mismatch fraction of the gap since the last accepted point.
    min_anchor_units: how many consecutive exact-match units are required at the end of a gap to
        "anchor" it and accept the jump (default 1, matching the original rule; >1 is a stricter
        variant tested in the rule-adjustment analysis).
    """
    motif = motif.upper()
    flank_seq = flank_seq.upper()
    motif_length = len(motif)
    max_offset = len(flank_seq)

    cumulative = cumulative_purity_and_hamming(flank_seq, motif, core_length, side)
    mismatched = cumulative["mismatched_bases"]  # mismatched[g-1] = mismatches in flank_seq[:g]

    def cumulative_mismatches_up_to(g):
        return int(mismatched[g - 1]) if g > 0 else 0

    def anchor_mismatches(g):
        """Mismatches within the last `min_anchor_units` units ending at offset g."""
        anchor_len = min_anchor_units * motif_length
        return cumulative_mismatches_up_to(g) - cumulative_mismatches_up_to(g - anchor_len)

    accepted = 0
    while True:
        found_new_accept = False
        n_units = 1
        while True:
            g = accepted + n_units * motif_length
            if g > max_offset:
                if verbose:
                    print(f"    ran out of flank sequence at g={g}")
                return accepted
            gap_len = g - accepted
            gap_mismatches = cumulative_mismatches_up_to(g) - cumulative_mismatches_up_to(accepted)
            frac = gap_mismatches / gap_len
            if frac > threshold:
                if verbose:
                    print(f"    gap of {gap_len}bp ({n_units} units) exceeds threshold: {frac:.3f} > {threshold:.3f} -- stop")
                return accepted
            if n_units >= min_anchor_units and anchor_mismatches(g) == 0:
                if verbose:
                    print(
                        f"    accept: gap={gap_len}bp ({n_units} units), frac={frac:.3f} <= {threshold:.3f}, "
                        f"last {min_anchor_units} unit(s) exact motif copy"
                    )
                accepted = g
                found_new_accept = True
                break
            n_units += 1
        if not found_new_accept:
            return accepted


def compute_candidate_definitions(chrom, start, end, motif, fasta, threshold=HAMMING_FRACTION_THRESHOLD, min_anchor_units=1):
    """Given a locus's current [start, end) (0-based half-open) and motif, apply the motif-reanchoring
    rule independently on each side and return the 4 candidate definitions (dedup'd): baseline,
    left-only-extended, right-only-extended, both-extended.
    """
    core_seq = fasta.fetch(chrom, start, end)
    core_length = len(core_seq)

    left_flank = fasta.fetch(chrom, max(0, start - MAX_OFFSET), start)[::-1]
    right_flank = fasta.fetch(chrom, end, end + MAX_OFFSET)

    left_extension = extend_by_motif_reanchoring(left_flank, motif, core_length, "left", threshold, min_anchor_units)
    right_extension = extend_by_motif_reanchoring(right_flank, motif, core_length, "right", threshold, min_anchor_units)

    definitions = {
        "baseline": (start, end),
        "left_only": (start - left_extension, end),
        "right_only": (start, end + right_extension),
        "both": (start - left_extension, end + right_extension),
    }
    # dedup identical regions (e.g. a side that doesn't extend at all)
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
    # locus_id, chrom, start_0based, end, motif, baseline_note
    ("EP400", "chr12", 132062548, 132062611, "CAG", "starting from narrower def, per gnomAD blog"),
    ("C9ORF72", "chr9", 27573528, 27573546, "GGCCCC", "KnownDiseaseAssociatedLoci current def"),
    ("RFC1", "chr4", 39348424, 39348479, "AARRG", "KnownDiseaseAssociatedLoci current def"),
    ("BEAN1", "chr16", 66490398, 66490453, "TAAAA", "KnownDiseaseAssociatedLoci current def"),
    ("RUNX2", "chr6", 45422750, 45422792, "GCN", "KnownDiseaseAssociatedLoci current def"),
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
