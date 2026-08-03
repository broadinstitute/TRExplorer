"""Apply the motif-reanchoring rule (thresholds 0.15 and 0.20) across the existing stratified
TRExplorer v1 sample (19,994 loci, data/tier3_trexplorer_v1_sample.tsv), to see how the rule behaves
on a broad, non-curated locus set (all literal motifs, no IUPAC degeneracy, unlike the known-disease
loci tested previously) and to find "disagreement" loci where 0.20 extends but 0.15 doesn't (or
extends less) -- these are the interesting cases for a follow-up real-read validation.
"""

import sys
from pathlib import Path

import pandas as pd
import pysam

sys.path.insert(0, str(Path(__file__).parent))
from motif_reanchor_extension import compute_candidate_definitions

FASTA_PATH = "/Users/weisburd/hg38.fa"
INPUT_TSV = Path(__file__).parent / "data" / "tier3_trexplorer_v1_sample.tsv"
OUTPUT_TSV = Path(__file__).parent / "data" / "threshold_sweep_trexplorer_v1.tsv"


def main():
    loci_df = pd.read_csv(INPUT_TSV, sep="\t")
    loci_df = loci_df[~loci_df["too_close_to_chrom_end"]]
    print(f"processing {len(loci_df)} loci")

    fasta = pysam.FastaFile(FASTA_PATH)
    rows = []
    for i, row in enumerate(loci_df.itertuples()):
        if i % 2000 == 0:
            print(f"  {i}/{len(loci_df)}")
        chrom, start, end, motif = row.chrom, int(row.start_0based), int(row.end), row.motif
        core_seq = fasta.fetch(chrom, start, end)
        if not core_seq:
            continue
        try:
            r15 = compute_candidate_definitions(chrom, start, end, motif, fasta, threshold=0.15)
            r20 = compute_candidate_definitions(chrom, start, end, motif, fasta, threshold=0.20)
        except Exception:
            continue
        ext15 = r15["left_extension"] + r15["right_extension"]
        ext20 = r20["left_extension"] + r20["right_extension"]
        if ext15 == 0 and ext20 == 0:
            continue  # only keep loci where at least one threshold extends -- that's the interesting set
        rows.append(
            {
                "locus_id": row.locus_id,
                "chrom": chrom,
                "start_0based": start,
                "end": end,
                "motif": motif,
                "motif_length": len(motif),
                "left_ext_015": r15["left_extension"],
                "right_ext_015": r15["right_extension"],
                "left_ext_020": r20["left_extension"],
                "right_ext_020": r20["right_extension"],
                "total_015": ext15,
                "total_020": ext20,
                "agree": ext15 == ext20,
            }
        )

    result_df = pd.DataFrame(rows)
    result_df.to_csv(OUTPUT_TSV, sep="\t", index=False)
    print(f"\nwrote {len(result_df)} rows (loci with nonzero extension under 0.15 or 0.20) to {OUTPUT_TSV}")

    n_015 = (result_df["total_015"] > 0).sum()
    n_020 = (result_df["total_020"] > 0).sum()
    n_disagree = (~result_df["agree"]).sum()
    print(f"extended under 0.15: {n_015}/{len(loci_df)}")
    print(f"extended under 0.20: {n_020}/{len(loci_df)}")
    print(f"disagreement (different extension amount): {n_disagree}")
    print("\nmotif length distribution among extended loci (0.20):")
    print(result_df[result_df["total_020"] > 0]["motif_length"].value_counts().sort_index())


if __name__ == "__main__":
    main()
