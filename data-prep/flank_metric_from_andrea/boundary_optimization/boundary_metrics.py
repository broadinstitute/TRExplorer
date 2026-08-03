"""Metrics for evaluating how far a TR locus boundary should be extended into its flanking sequence.

All metrics are computed as a function of `offset` = distance (bp) from the current annotated
boundary, independently for the left and right flank. Convention used throughout: the annotated
repeat region is assumed to be phase-aligned such that core_seq[0] == motif[0] (this matches how
ExpansionHunter's own computeRepeatSequencePurity anchors reference-repeat purity in
RepeatPurity.cpp / io/LocusSpecDecoding.cpp -- no rotation search is performed).

Given that convention, the motif phase at any genomic position relative to the core's start
(core start = position 0) is simply `relative_position % motif_length` (Python's modulo already
returns a non-negative result for negative relative_position, which is what left-flank positions
need).
"""

import numpy as np

# IUPAC ambiguity codes: the set of literal bases each code represents. A motif position using one
# of these (e.g. RFC1's "AARRG", R = A-or-G; RUNX2's "GCN", N = any base) should match any base in
# its set, not just its own literal letter -- otherwise a real, biologically-correct repeat unit gets
# counted as mismatching at every degenerate position. This mirrors what ExpansionHunter's own
# LocusStructure motif syntax means, though EH's own computeRepeatSequencePurity (RepeatPurity.cpp)
# does NOT special-case this (a literal-char comparison, so 'N' in the motif never matches a real
# base there either) -- so this is a deliberate improvement over EH's own convention, not a
# reproduction of it. Values are frozensets of ASCII codes for fast membership use in numpy.
IUPAC_CODE_TO_BASES = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "GC", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}


def _build_iupac_match_table():
    """256x256 boolean lookup: table[seq_base_ord][motif_base_ord] = True if seq_base is a literal
    base consistent with motif_base's IUPAC code (seq_base is always a literal A/C/G/T in practice).
    """
    table = np.zeros((256, 256), dtype=bool)
    for motif_code, represented_bases in IUPAC_CODE_TO_BASES.items():
        for seq_base in represented_bases:
            table[ord(seq_base), ord(motif_code)] = True
    return table


_IUPAC_MATCH_TABLE = _build_iupac_match_table()


def _flank_phases(core_length, motif_length, side, max_offset):
    """Motif-phase (index into motif) for each of the first `max_offset` bases of the given flank.

    side='right': flank[k] is at relative position core_length + k.
    side='left': flank[k] (k=0 is the base immediately left of the core) is at relative position -(k + 1).
    """
    offsets = np.arange(max_offset)
    if side == "right":
        relative_positions = core_length + offsets
    elif side == "left":
        relative_positions = -(offsets + 1)
    else:
        raise ValueError(f"side must be 'left' or 'right', got {side!r}")
    return relative_positions % motif_length


def cumulative_purity_and_hamming(flank_seq, motif, core_length, side):
    """Cumulative purity/hamming-distance of the region [boundary, boundary + offset) for offset in 1..len(flank_seq).

    Returns a dict of numpy arrays (index i = offset i+1): matched_bases, mismatched_bases (= Hamming
    distance to a perfect tiled repeat), total_bases, purity (matched / total).
    This is the same matched-base convention ExpansionHunter uses for ReferenceRepeatPurity.
    """
    max_offset = len(flank_seq)
    phases = _flank_phases(core_length, len(motif), side, max_offset)
    motif_array = np.frombuffer(motif.upper().encode(), dtype=np.uint8)
    flank_array = np.frombuffer(flank_seq.upper().encode(), dtype=np.uint8)
    motif_at_phase = motif_array[phases]
    # IUPAC ambiguity codes in the MOTIF (not the sequence) mark degenerate positions -- e.g. the
    # polyalanine loci GCN/NGC (N = any base) or RFC1's AARRG (R = A-or-G) -- match any base in
    # their represented set, not just their own literal letter.
    is_match = _IUPAC_MATCH_TABLE[flank_array, motif_at_phase]
    matched_cumulative = np.cumsum(is_match)
    total = np.arange(1, max_offset + 1)
    return {
        "matched_bases": matched_cumulative,
        "mismatched_bases": total - matched_cumulative,
        "total_bases": total,
        "purity": matched_cumulative / total,
    }
