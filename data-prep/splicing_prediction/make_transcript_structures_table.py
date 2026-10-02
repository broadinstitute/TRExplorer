"""Build the transcript-structure table SAI-10k reads, for every GENCODE v50 comprehensive transcript.

Production server.py reads this table (transcripts_hg38) from Cloud SQL. Since SpliceAI-lookup commit
78030ec, build_and_deploy.py's update_transcript_tables loads it from the comprehensive genePred file
gencode.v50.GRCh38.comprehensive.sorted.txt.gz. This script reads that same file and applies the same
conversions (version suffix stripped, the unsuffixed copy of PAR / duplicate names wins,
cdsStart == cdsEnd means non-coding), so the benchmark and the full run see what production sees,
without a database connection. The output has the format of a transcripts_hg38 export:

    transcript_id  strand  cds_start  cds_end  exon_starts  exon_ends     (0-based starts, \\N for no CDS)

(The GENCODE v49 production table held only 280,000 transcripts, because the basic genePred file
overwrote the comprehensive one; the v49 version of this script rebuilt the missing ones from the GTF.)

Usage:
    python3 make_transcript_structures_table.py
"""
import gzip
import os

HERE = os.path.dirname(os.path.abspath(__file__))
GENE_PRED = "/Users/weisburd/code/SpliceAI-lookup/google_cloud_run_services/gencode.v50.GRCh38.comprehensive.sorted.txt.gz"
OUTPUT = os.path.join(HERE, "transcripts_hg38.gencode_v50_comprehensive.tsv")


def gene_pred_rows_to_table(gene_pred_lines):
    """Converts sorted genePredExt lines to {transcript_id: tab-joined table row}, as update_transcript_tables loads them.

    Each line starts with the index column the sorted file adds, then the genePredExt columns.
    """
    rows = [line.rstrip("\n").split("\t")[1:] for line in gene_pred_lines if line.strip()]
    # Stripping the version maps PAR copies (ENST..._PAR_Y) and a few others onto the same ID; the
    # production loader keeps the last row inserted and moves the unsuffixed names to the end
    # (a stable sort, so the file order decides among the rest).
    rows.sort(key=lambda fields: "_" not in fields[0])
    table = {}
    for name, chrom, strand, tx_start, tx_end, cds_start, cds_end, exon_count, exon_starts, exon_ends, *_ in rows:
        if cds_start == cds_end:
            cds_start = cds_end = "\\N"
        transcript_id = name.split(".")[0]
        table[transcript_id] = "\t".join([transcript_id, strand, cds_start, cds_end, exon_starts, exon_ends])
    return table


def main():
    with gzip.open(GENE_PRED, "rt") as f:
        table = gene_pred_rows_to_table(f)
    with open(OUTPUT, "w") as f:
        for transcript_id in sorted(table):
            f.write(table[transcript_id] + "\n")
    print(f"Wrote {len(table):,} transcripts to {OUTPUT}")


if __name__ == "__main__":
    main()
