"""Decompose TRGT multisample VCF alleles into motif compositions.

Reads VCF lines from stdin (header `#` lines are ignored) and writes one
parquet row per ``INFO/TRID`` sub-entry to ``--output``. A simple VCF record
emits one row; a compound record (TRID with comma-separated sub-TRIDs that
TRGT uses to cluster adjacent catalog entries into one VCF row) emits one
row per sub-TRID, each decomposed under the motif parsed from that sub-TRID.
Coordinates and motif come from the sub-TRID itself, not from VCF POS/END/
MOTIFS, since TRGT sometimes extends VCF POS to anchor at an adjacent locus
while keeping each sub-TRID faithful to its originating catalog entry.

Decomposition uses trviz when feasible, with bypass paths for short motifs
and skip reasons for inputs that are too long or contain non-ACGT bases.
"""

import argparse
import collections
import re
import sys

import pyarrow as pa
import pyarrow.parquet as pq


# Hard cap on per-allele sequence length passed to trviz / mono-di bypass.
MAX_SEQ_LEN = 50_000

# trviz dynamic-programming memory cap: len(seq) * num_motifs * max_motif_len.
MAX_DP_PRODUCT = 50_000_000

# Parquet write buffer size (rows). Kept small because deeply-compound records
# can emit per-row allele_data in the megabytes; with 10k rows the in-memory
# batch would exceed multi-GB on chr1-like regions and OOM the standard worker.
PARQUET_CHUNK_ROWS = 1_000

# Regex used to split GT field on `/` or `|`.
_GT_SPLIT_RE = re.compile(r"[/|]")

# Note: there is no "interruption sentinel" token — when the mono/di bypass
# encounters a non-motif run it emits the literal substring (e.g. "GA"). The
# frontend distinguishes motif tokens from non-motif tokens by membership in
# the locus's INFO/MOTIFS list.

# Defer the trviz import — let it raise at first use so a missing dep does
# not prevent --help from working in environments without trviz installed.
_decomposer = None


def _get_decomposer():
    """Returns a cached trviz Decomposer instance (lazy-imported)."""
    global _decomposer
    if _decomposer is None:
        from trviz.decomposer import Decomposer
        _decomposer = Decomposer()
    return _decomposer


def chrom_to_index(chrom):
    """Maps a chromosome name to an integer index.

    Args:
        chrom: Chromosome name with or without the ``chr`` prefix.

    Returns:
        An integer where 1..22 map to themselves, X=23, Y=24, MT/M=25.
        Returns 0 for any other contig.
    """
    name = chrom[3:] if chrom.startswith("chr") else chrom
    if name == "X":
        return 23
    if name == "Y":
        return 24
    if name == "MT" or name == "M":
        return 25
    if name.isdigit():
        return int(name)
    return 0


def parse_info(info):
    """Parses a VCF INFO column into a dict.

    Args:
        info: The raw INFO field text (e.g. ``"TRID=...;END=...;MOTIFS=..."``).

    Returns:
        A dict mapping key to value. Flag-style keys (no ``=``) map to ``""``.
    """
    out = {}
    for entry in info.split(";"):
        if not entry:
            continue
        if "=" in entry:
            key, value = entry.split("=", 1)
            out[key] = value
        else:
            out[entry] = ""
    return out


def mono_di_bypass(seq, motifs):
    """Decomposes ``seq`` using only mono/di-nucleotide motifs (no trviz).

    Greedy longest-motif match at each position. When no motif matches at
    position ``i``, walks forward to the next position where some motif does
    match and emits the entire intervening substring as a single literal
    token. Concatenating all returned tokens reproduces ``seq`` exactly.

    Args:
        seq: Allele sequence (ACGT or any character — bypass tolerates
            non-ACGT bases by emitting them as literal tokens).
        motifs: Iterable of motif strings, all of length 1 or 2.

    Returns:
        List of tokens. Tokens that are members of ``motifs`` are decomposed
        repeats; all other tokens are literal non-motif substrings (which the
        frontend renders with ``OTHER_MOTIFS_COLOR``).

    Example:
        >>> mono_di_bypass("AAAAA", ["A"])
        ['A', 'A', 'A', 'A', 'A']
        >>> mono_di_bypass("GTGTGTGTGTGA", ["GT"])
        ['GT', 'GT', 'GT', 'GT', 'GT', 'GA']
    """
    # Sort motifs longest-first so a 2-bp motif beats a 1-bp motif when both match.
    sorted_motifs = sorted(set(motifs), key=lambda m: -len(m))
    tokens = []
    i = 0
    n = len(seq)
    while i < n:
        matched = None
        for m in sorted_motifs:
            if seq.startswith(m, i):
                matched = m
                break
        if matched is None:
            j = i + 1
            while j < n and not any(seq.startswith(m, j) for m in sorted_motifs):
                j += 1
            tokens.append(seq[i:j])
            i = j
        else:
            tokens.append(matched)
            i += len(matched)
    return tokens


def decompose_or_skip(seq, motifs):
    """Decomposes ``seq`` into a list of motif tokens or returns a skip reason.

    Checks are evaluated in order; the FIRST matching branch wins:

      1. ``len(seq) > MAX_SEQ_LEN`` -> skip (``too_long``).
      2. ``seq`` contains non-ACGT chars AND motifs are not all mono/di
         -> skip (``contains_invalid_bases``).
      3. All motifs are length 1 or 2 -> mono/di bypass.
      4. DP-table size estimate exceeds ``MAX_DP_PRODUCT`` -> skip (``too_long``).
      5. Otherwise -> trviz ``Decomposer().decompose(seq, list(motifs))``;
         any exception is treated as ``contains_invalid_bases``.

    Args:
        seq: Allele sequence (post-anchor-strip, no leading anchor base).
        motifs: Iterable of motif strings from the locus' INFO/MOTIFS.

    Returns:
        Tuple ``(composition_csv, skip_reason)`` where ``composition_csv`` is
        a comma-joined motif sequence (empty string if skipped) and
        ``skip_reason`` is one of ``"", "too_long", "contains_invalid_bases"``.
    """
    if len(seq) > MAX_SEQ_LEN:
        return "", "too_long"

    motifs_list = list(motifs)
    all_mono_di = all(len(m) <= 2 for m in motifs_list) if motifs_list else False

    # Only screen for non-ACGT when we'd be feeding the trviz path; the
    # mono/di bypass itself tolerates non-motif bases via interruption tokens.
    if not all_mono_di:
        if any(c not in "ACGT" for c in seq):
            return "", "contains_invalid_bases"

    if all_mono_di:
        return ",".join(mono_di_bypass(seq, motifs_list)), ""

    max_motif_len = max(len(m) for m in motifs_list) if motifs_list else 0
    if len(seq) * len(motifs_list) * max_motif_len > MAX_DP_PRODUCT:
        return "", "too_long"

    try:
        tokens = _get_decomposer().decompose(seq, motifs_list)
    except Exception:
        return "", "contains_invalid_bases"
    return ",".join(tokens), ""


def build_seq_counts(ref_repeat, alts, sample_fields, gt_idx):
    """Tallies haplotype counts per unique allele sequence.

    Args:
        ref_repeat: REF sequence with TRGT anchor base stripped.
        alts: List of ALT sequences (anchor-stripped except for ``*``).
        sample_fields: Per-sample FORMAT fields from the VCF record.
        gt_idx: Index of ``GT`` within the FORMAT colon-separated keys.

    Returns:
        Tuple ``(seq_counts, total_called)`` where ``seq_counts`` is an
        ``OrderedDict`` mapping sequence to integer haplotype count (REF
        inserted first even when its count is 0) and ``total_called`` is the
        sum of all called haplotype counts across samples.
    """
    counts = collections.Counter()
    total_called = 0
    for s in sample_fields:
        # Pull out just the GT subfield without splitting the whole record.
        gt = s.split(":", gt_idx + 1)[gt_idx]
        for h in _GT_SPLIT_RE.split(gt):
            if h == "." or h == "":
                continue
            counts[int(h)] += 1
            total_called += 1

    seq_counts = collections.OrderedDict()
    seq_counts[ref_repeat] = counts.get(0, 0)
    for i, alt in enumerate(alts, start=1):
        seq_counts[alt] = seq_counts.get(alt, 0) + counts.get(i, 0)
    return seq_counts, total_called


def format_allele_data(seq_counts, ref_repeat, motif):
    """Builds the delimited ``allele_data`` STRING for one (sub-)TRID.

    Args:
        seq_counts: Dict mapping sequence to integer haplotype count.
        ref_repeat: REF sequence (anchor-stripped) — used to set ``is_ref``.
        motif: The single motif string parsed from this row's TRID. Each
            sub-TRID of a compound record gets its own motif here, so trviz
            sees just that motif when decomposing the full compound allele.

    Returns:
        A string of ``;``-separated allele records sorted by descending count.
        Each record is ``sequence|composition|count|is_ref`` for the common
        decomposition-succeeded case (4 fields), or
        ``sequence|composition|count|is_ref|skip_reason`` when the allele was
        skipped (5 fields, e.g. ``too_long``, ``contains_invalid_bases``,
        ``spanning_del``). Motif units are separated by ``,`` inside
        ``composition``. The trailing ``|skip_reason`` is omitted for
        non-skipped alleles to save ~25 MB across the BQ table.
    """
    motifs = [motif]
    records = []
    for seq, cnt in sorted(seq_counts.items(), key=lambda kv: -kv[1]):
        if seq == "*":
            comp_csv, skip_reason = "", "spanning_del"
        else:
            comp_csv, skip_reason = decompose_or_skip(seq, motifs)
        is_ref = "1" if seq == ref_repeat else "0"
        record = f"{seq}|{comp_csv}|{cnt}|{is_ref}"
        if skip_reason:
            record += f"|{skip_reason}"
        records.append(record)
    return ";".join(records)


def parse_trid(trid):
    """Splits a TRID like ``"4-3074876-3074933-CAG"`` into its components.

    Args:
        trid: TRID string in ``chrom-start_0based-end_1based-motif`` form.

    Returns:
        Tuple ``(chrom, start_0based, end_1based, motif)`` or ``None`` if the
        TRID doesn't have exactly 4 dash-separated parts or the start/end
        components aren't integers. Motifs are nucleotide strings that never
        contain a dash, so a 4-way split is safe.
    """
    parts = trid.split("-", 3)
    if len(parts) != 4:
        return None
    chrom, start_s, end_s, motif = parts
    try:
        return (chrom, int(start_s), int(end_s), motif)
    except ValueError:
        return None


def process_record(line):
    """Parses one VCF data line and returns a list of dicts for parquet emission.

    Each VCF record may carry a simple TRID (one sub-TRID) or a compound TRID
    (multiple comma-separated sub-TRIDs, when TRGT clusters adjacent catalog
    entries). We emit one parquet row per sub-TRID. Coordinates and motif are
    parsed from the sub-TRID itself rather than from VCF POS/END/MOTIFS, since
    TRGT can extend VCF POS to anchor at an adjacent locus while leaving the
    TRID faithful to the originating catalog entry.

    Args:
        line: One non-header VCF line (no trailing newline required).

    Returns:
        A list of dicts matching the disk schema (one per sub-TRID), or an
        empty list if the line is empty / malformed / has no TRID / has a TRID
        that can't be parsed as ``chrom-start-end-motif``.
    """
    line = line.rstrip("\n")
    if not line:
        return []
    fields = line.split("\t")
    if len(fields) < 10:
        return []
    chrom, _pos, _id, ref, alt_field, _qual, _filter, info, fmt = fields[:9]
    sample_fields = fields[9:]

    info_dict = parse_info(info)
    trid_field = info_dict.get("TRID", "")
    if not trid_field:
        return []

    # Strip the TRGT anchor base; spanning deletion `*` is preserved as-is.
    # The same anchor-stripped allele sequences are shared across every
    # sub-TRID for a compound record (each sub-TRID re-decomposes the full
    # compound allele under its own motif).
    ref_repeat = ref[1:]
    raw_alts = [] if alt_field == "." else alt_field.split(",")
    alts = [a[1:] if a != "*" else a for a in raw_alts]

    gt_idx = fmt.split(":").index("GT")
    seq_counts, total_called = build_seq_counts(ref_repeat, alts, sample_fields, gt_idx)
    total_unique = len(seq_counts)

    chrom_idx_for_vcf = chrom_to_index(chrom)

    rows = []
    for sub_trid in trid_field.split(","):
        parsed = parse_trid(sub_trid)
        if parsed is None:
            continue
        sub_chrom, sub_start, sub_end, sub_motif = parsed
        if not sub_motif:
            continue
        rows.append({
            "locus_id": sub_trid,
            "chrom": sub_chrom,
            "chrom_index": chrom_idx_for_vcf,
            "start_0based": sub_start,
            "end_1based": sub_end,
            "motifs": sub_motif,
            "total_allele_number": total_called,
            "total_unique_alleles": total_unique,
            "allele_data": format_allele_data(seq_counts, ref_repeat, sub_motif),
        })
    return rows


def get_schema():
    """Returns the pyarrow schema for the output parquet file."""
    return pa.schema([
        ("locus_id", pa.string()),
        ("chrom", pa.string()),
        ("chrom_index", pa.int64()),
        ("start_0based", pa.int64()),
        ("end_1based", pa.int64()),
        ("motifs", pa.string()),
        ("total_allele_number", pa.int64()),
        ("total_unique_alleles", pa.int64()),
        ("allele_data", pa.string()),
    ])


def _flush_batch(writer, batch):
    """Writes ``batch`` (a dict of column -> list) to ``writer`` and clears it."""
    if not batch["locus_id"]:
        return
    table = pa.Table.from_pydict(batch, schema=writer.schema)
    writer.write_table(table)
    for column in batch.values():
        column.clear()


def stream_vcf_to_parquet(input_stream, output_path, max_rows=None):
    """Streams VCF text into a parquet file, one row per sub-TRID.

    Args:
        input_stream: Iterable of lines (e.g. ``sys.stdin``).
        output_path: Destination parquet path.
        max_rows: If not ``None``, stop after emitting this many output rows
            (sub-TRIDs, not VCF lines).

    Returns:
        Tuple ``(rows_written, vcf_lines_skipped)``. ``vcf_lines_skipped``
        counts non-header VCF lines that ``process_record`` returned an empty
        list for (empty/malformed lines, missing TRID, unparseable TRID).
    """
    schema = get_schema()
    batch = {name: [] for name in schema.names}
    rows_written = 0
    vcf_lines_skipped = 0

    with pq.ParquetWriter(output_path, schema) as writer:
        for line in input_stream:
            if line.startswith("#"):
                continue
            rows = process_record(line)
            if not rows:
                vcf_lines_skipped += 1
                continue
            for row in rows:
                for name in schema.names:
                    batch[name].append(row[name])
                rows_written += 1
                if len(batch["locus_id"]) >= PARQUET_CHUNK_ROWS:
                    _flush_batch(writer, batch)
                if max_rows is not None and rows_written >= max_rows:
                    break
            if max_rows is not None and rows_written >= max_rows:
                break
        _flush_batch(writer, batch)

    return rows_written, vcf_lines_skipped


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, help="Output parquet path.")
    parser.add_argument(
        "--max-rows",
        type=int,
        default=None,
        help="Optional cap on number of VCF records processed (for local testing).",
    )
    args = parser.parse_args()

    rows_written, vcf_lines_skipped = stream_vcf_to_parquet(sys.stdin, args.output, max_rows=args.max_rows)
    print(f"Wrote {rows_written} rows to {args.output} (skipped {vcf_lines_skipped} unparseable VCF lines)")


if __name__ == "__main__":
    main()
