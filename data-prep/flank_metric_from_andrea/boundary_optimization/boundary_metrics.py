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

import itertools

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


def iupac_match(seq_base, motif_base):
    """True if literal seq_base is consistent with motif_base's IUPAC ambiguity code."""
    return bool(_IUPAC_MATCH_TABLE[ord(seq_base.upper()), ord(motif_base.upper())])


def _expand_iupac_kmer(kmer):
    """All literal ACGT k-mers consistent with `kmer`, which may itself contain IUPAC codes."""
    options = [IUPAC_CODE_TO_BASES.get(base.upper(), base.upper()) for base in kmer]
    return {"".join(combo) for combo in itertools.product(*options)}


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


def cumulative_edit_distance(flank_seq, motif, core_length, side):
    """Cumulative DP edit (Levenshtein) distance between flank_seq[:offset] and a perfect tiled
    repeat of the same length, for every offset in 1..len(flank_seq), extracted from a single DP table.

    Standard edit-distance DP has the property that dp[i][j] = edit_distance(s1[:i], s2[:j]) for all
    i, j -- so dp[i][i] for i=1..max_offset gives the whole cumulative curve from ONE O(n^2) table,
    instead of recomputing a fresh DP per offset (which would cost O(n^3)).
    """
    max_offset = len(flank_seq)
    # tile must follow the same per-position phase sequence as cumulative_purity_and_hamming/
    # _flank_phases -- for side='left' the phase DECREASES with offset (flank_seq is stored
    # boundary-adjacent-base-first, i.e. reversed genomic order), so a simple forward-tiled motif
    # starting at phase0 is wrong for the left side; build the tile position-by-position instead.
    phases = _flank_phases(core_length, len(motif), side, max_offset)
    motif_upper = motif.upper()
    tile = "".join(motif_upper[p] for p in phases)

    s1 = flank_seq.upper()
    s2 = tile.upper()
    n = max_offset

    previous_row = list(range(n + 1))
    diagonal = [0] * (n + 1)
    for i in range(1, n + 1):
        current_row = [i] + [0] * n
        s1_base = s1[i - 1]
        for j in range(1, n + 1):
            substitution_cost = 0 if iupac_match(s1_base, s2[j - 1]) else 1  # IUPAC code in the motif tile = wildcard
            current_row[j] = min(
                previous_row[j] + 1,  # deletion
                current_row[j - 1] + 1,  # insertion
                previous_row[j - 1] + substitution_cost,  # substitution / match
            )
        diagonal[i] = current_row[i]
        previous_row = current_row

    cumulative = np.array(diagonal[1:])
    total = np.arange(1, max_offset + 1)
    return {
        "edit_distance": cumulative,
        "normalized_edit_distance": cumulative / total,
    }


def local_linguistic_complexity_curve(flank_seq, window_size, k=3):
    """Andrea metric #1 (linguistic complexity), as a sliding-window curve.

    value at position i = (# distinct k-mers in flank_seq[i : i+window_size]) / min(window_size - k + 1, 4**k)
    0 = perfectly repetitive window, 1 = maximally complex (all k-mers unique, within the cap of 4**k).
    Uses an O(1)-amortized rolling window (add/drop one k-mer per step) rather than recomputing per offset.
    """
    n_positions = len(flank_seq) - window_size + 1
    if n_positions <= 0:
        return np.array([])
    max_possible = min(window_size - k + 1, 4**k)

    kmers = [flank_seq[i : i + k] for i in range(len(flank_seq) - k + 1)]
    counts = {}
    values = np.empty(n_positions)

    def add(kmer):
        counts[kmer] = counts.get(kmer, 0) + 1

    def remove(kmer):
        counts[kmer] -= 1
        if counts[kmer] == 0:
            del counts[kmer]

    n_kmers_in_window = window_size - k + 1
    for kmer in kmers[:n_kmers_in_window]:
        add(kmer)
    values[0] = len(counts) / max_possible

    for i in range(1, n_positions):
        remove(kmers[i - 1])
        add(kmers[i + n_kmers_in_window - 1])
        values[i] = len(counts) / max_possible

    return values


def local_kmer_jaccard_curve(flank_seq, motif, window_size, side, k=3):
    """Andrea metric #2 (k-mer Jaccard similarity to the motif), as a sliding-window curve.

    value at position i = Jaccard(kmer_set(flank_seq[i:i+window_size]), kmer_set(tiled motif)).
    The tiled-motif k-mer set is periodic (at most len(motif) distinct k-mers for a pure repeat), so
    it's computed once, not per offset.
    """
    n_positions = len(flank_seq) - window_size + 1
    if n_positions <= 0:
        return np.array([])

    # flank_seq for side='left' is stored boundary-adjacent-base-first (reversed genomic order, same
    # convention as cumulative_purity_and_hamming/cumulative_edit_distance), so a genuinely perfect
    # left-extension reads as the motif REVERSED and tiled, not forward-tiled.
    tiling_motif = motif[::-1] if side == "left" else motif
    # a pure repeat has at most len(motif) distinct k-mers (one per phase); tile a couple of periods
    # so wraparound k-mers are included, then take exactly one period's worth of window positions.
    motif_tile = tiling_motif * (2 + k // max(1, len(tiling_motif)))
    # expand any IUPAC ambiguity codes in the motif (e.g. 'N' in GCN/NGC) into every literal ACGT
    # possibility -- flank k-mers are always literal genomic sequence, so a literal-vs-ambiguity-code
    # set intersection would otherwise always be empty for any degenerate motif.
    motif_kmers = set()
    for i in range(len(tiling_motif)):
        motif_kmers |= _expand_iupac_kmer(motif_tile[i : i + k])

    kmers = [flank_seq[i : i + k] for i in range(len(flank_seq) - k + 1)]
    counts = {}
    values = np.empty(n_positions)

    def add(kmer):
        counts[kmer] = counts.get(kmer, 0) + 1

    def remove(kmer):
        counts[kmer] -= 1
        if counts[kmer] == 0:
            del counts[kmer]

    n_kmers_in_window = window_size - k + 1
    for kmer in kmers[:n_kmers_in_window]:
        add(kmer)

    def jaccard():
        window_kmer_set = set(counts.keys())
        union = window_kmer_set | motif_kmers
        if not union:
            return 0.0
        return len(window_kmer_set & motif_kmers) / len(union)

    values[0] = jaccard()
    for i in range(1, n_positions):
        remove(kmers[i - 1])
        add(kmers[i + n_kmers_in_window - 1])
        values[i] = jaccard()

    return values


def compute_boundary_extension_profile(
    left_flank_seq, core_seq, right_flank_seq, motif, max_offset=300, local_window_size=24, compute_edit_distance=True
):
    """Compute all boundary-extension metrics for both sides of one locus.

    left_flank_seq / right_flank_seq: raw genomic sequence immediately outside the annotated core,
        in forward genomic orientation, each at least `max_offset` bases long.
    core_seq: the annotated repeat region sequence (only its length is used).
    motif: the annotated repeat unit, phase-aligned so core_seq[0] == motif[0].

    Returns {'left': {...}, 'right': {...}}, each side a dict of the metric curves described in the
    individual metric functions above (indexed by offset - 1, i.e. curve[0] is offset=1bp from boundary).
    """
    core_length = len(core_seq)
    core_seq = core_seq.upper()
    # core purity: same phase formula as the right-side flank with core_length=0, since core[0] == motif[0]
    core_matched_total = int(cumulative_purity_and_hamming(core_seq, motif, 0, "right")["matched_bases"][-1])

    profile = {}
    for side, raw_flank in (("left", left_flank_seq), ("right", right_flank_seq)):
        flank = raw_flank[::-1] if side == "left" else raw_flank
        flank = flank[:max_offset].upper()

        purity_and_hamming = cumulative_purity_and_hamming(flank, motif, core_length, side)
        side_profile = dict(purity_and_hamming)

        # full-region curves: same convention ExpansionHunter uses for ReferenceRepeatPurity, i.e. purity
        # of the WHOLE region [core + extension], not just the newly-added flank segment. Since matching
        # is per-position, this is just the core's fixed matched-base count plus the flank's running total.
        total_bases = core_length + side_profile["total_bases"]
        total_matched = core_matched_total + side_profile["matched_bases"]
        side_profile["full_region_purity"] = total_matched / total_bases
        side_profile["full_region_mismatched_bases"] = total_bases - total_matched
        side_profile["full_region_total_bases"] = total_bases

        if compute_edit_distance:
            side_profile.update(cumulative_edit_distance(flank, motif, core_length, side))

        side_profile["linguistic_complexity"] = local_linguistic_complexity_curve(flank, local_window_size)
        side_profile["kmer_jaccard_similarity"] = local_kmer_jaccard_curve(flank, motif, local_window_size, side)

        profile[side] = side_profile

    return profile
