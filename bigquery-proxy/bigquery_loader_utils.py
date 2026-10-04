"""Helpers shared by the scripts that load per-locus HPRC256 tables into BigQuery."""

import gzip

import tqdm


def scan_for_smallest_trid(input_tsv, col_idx, locus_id_column="locus_id", interval_column="interval",
                           trid_column="trid"):
    """Find the row to keep for each (locus id, interval) in a per-locus TSV.

    Two VCF records can describe the same repeat over the same interval under different TRIDs, so
    the upstream TSVs now emit a row for each. These BigQuery schemas have no TRID column to tell
    those rows apart, and the website fetches them with an unordered LIMIT 1, so a loader keeps only
    one row per (locus id, interval): the one with the smallest TRID, which is the same row
    load_bigquery_main_table.py picks when it builds the main catalog table.

    Only the key and the TRID are held in memory, since the histogram and distribution columns these
    TSVs carry are large.

    Args:
        input_tsv: path of the gzipped TSV.
        col_idx: column name -> index, from the TSV header.
        locus_id_column: name of the locus id column ("locus_id" or "LocusId").
        interval_column: name of the interval column ("interval" or "Interval").
        trid_column: name of the TRID column ("trid" or "TRID").

    Returns:
        A dictionary of "<locus id>\\t<interval>" -> the smallest TRID seen for that key.
    """
    locus_i, interval_i = col_idx[locus_id_column], col_idx[interval_column]
    trid_i = col_idx[trid_column]
    smallest_trid_by_key = {}
    with gzip.open(input_tsv, "rt") as f:
        f.readline()
        for line in tqdm.tqdm(f, unit=" rows", unit_scale=True, desc="Finding the smallest TRID per locus"):
            fields = line.rstrip("\n").split("\t")
            key = f"{fields[locus_i]}\t{fields[interval_i]}"
            trid = fields[trid_i]
            if key not in smallest_trid_by_key or trid < smallest_trid_by_key[key]:
                smallest_trid_by_key[key] = trid
    return smallest_trid_by_key


def is_smallest_trid_row(fields, col_idx, smallest_trid_by_key, locus_id_column="locus_id",
                         interval_column="interval", trid_column="trid"):
    """Return whether this row is the one scan_for_smallest_trid chose for its (locus id, interval).

    Always True when `smallest_trid_by_key` is empty, which is how a TSV written before the TRID
    tie-break was added (no TRID column, already one row per key) loads unchanged.
    """
    if not smallest_trid_by_key:
        return True
    key = f"{fields[col_idx[locus_id_column]]}\t{fields[col_idx[interval_column]]}"
    return fields[col_idx[trid_column]] == smallest_trid_by_key[key]
