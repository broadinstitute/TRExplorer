import pandas as pd
import numpy as np
from str_analysis.utils.canonical_repeat_unit import compute_canonical_motif

def get_individual_repeat_id_that_matches_longest_pure_segment_motif(row):
    trid = row["TRID"]
    if "," not in trid:  
        return trid
    
    trid_parts = trid.split(",")
    trid2 = None
    for trid_part in trid_parts:
        if trid_part.count("-") != 3:
            return None
        trid_part_canonical_motif = compute_canonical_motif(trid_part.split("-")[3])
        if trid_part_canonical_motif == row["longestPureSegmentMotif"]:
            if trid2 is None:
                trid2 = trid_part
            else:
                return None

    return trid2

df = pd.read_table("AoULR_phase1_TRGT_Weisburd_v1_lpsStats.txt.gz")
df = df[df.Stdev.notna()]
#df = df.sample(10**6)

df["longestPureSegmentMotif"] = df["longestPureSegmentMotif"].apply(compute_canonical_motif)
df["TRID2"] = df.apply(get_individual_repeat_id_that_matches_longest_pure_segment_motif, axis=1)

before = len(df)
df = df[df.TRID2.notna()]
print(f"Kept {len(df):,} out of {before:,} ({100 * len(df) / before:.2f}%) rows where TRID2 is not None")

df["motif"] = df["TRID2"].str.split("-").str[3]
df["canonical_motif"] = df["motif"].apply(compute_canonical_motif)

before = len(df)
df = df[df.canonical_motif == df.longestPureSegmentMotif]
print(f"Kept {len(df):,} out of {before:,} ({100 * len(df) / before:.2f}%) rows where canonical motif == longestPureSegmentMotif")

df["chrom"] = df["TRID2"].str.split("-").str[0]
df["motif_size"] = df["motif"].str.len()

df["Mode"] = df.apply(lambda m: m["Mode"] / len(m["longestPureSegmentMotif"]), axis=1)
df["Stdev"] = df.apply(lambda m: m["Stdev"] / len(m["longestPureSegmentMotif"]), axis=1)
df["isVariationCluster"] = np.where(df["TRID"].str.contains(","), 1, 0)

df = df[["TRID2", "Mode", "Stdev", "isVariationCluster", "chrom", "motif_size"]]

df.rename(columns={"Mode": "mode_AoU", "Stdev": "stdev_AoU"}, inplace=True)
df = df.set_index("TRID2")


df2 = pd.read_table("hprc_lps.2025_05.grouped_by_locus_and_motif.with_biallelic_histogram.tsv.gz")
df2 = df2[["locus_id", "mode_allele", "stdev"]]
df2.rename(columns={"mode_allele": "mode_HPRC", "stdev": "stdev_HPRC"}, inplace=True)
df2 = df2.set_index("locus_id")

df = df.join(df2, how="inner").reset_index()
df["diff_AoU_HPRC_mode"] = df["mode_AoU"] - df["mode_HPRC"]
df["diff_AoU_HPRC_stdev"] = df["stdev_AoU"] - df["stdev_HPRC"]
df["ratio_AoU_stdev_overHPRC_stdev"] = np.where(df["stdev_HPRC"] == 0, 0, df["stdev_AoU"] / df["stdev_HPRC"])  

df = df[sorted(df.columns)]
output_filename = "AoU_vs_HPRC256_comparison_table.tsv"
df.to_csv(output_filename, sep="\t", index=False, header=True)
print(f"Wrote {len(df):,} rows to {output_filename}")
