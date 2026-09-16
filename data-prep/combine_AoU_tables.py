"""Combines multiple AoU (All of Us) data tables into a single output table.

This script merges three input tables:
- AoULR_phase1_TRGT_Weisburd_v1_lpsStats.txt.gz
- AoULR_phase1_TRGT_Weisburd_v1_distinctAndTotalAlleles.txt.gz
- AoULR_phase1_TRGT_Weisburd_v1.0.1_TRConstraint.txt.gz

It filters to keep only loci where the canonical motif matches the longest pure segment motif,
computes additional statistics, and outputs a combined table with OE_len_percentile and
StdevRankByMotif columns.
"""

import pandas as pd
from str_analysis.utils.canonical_repeat_unit import compute_canonical_motif

filenames = [
    "AoULR_phase1_TRGT_Weisburd_v1_lpsStats.txt.gz",
    "AoULR_phase1_TRGT_Weisburd_v1_distinctAndTotalAlleles.txt.gz",
    "AoULR_phase1_TRGT_Weisburd_v1.0.1_TRConstraint.txt.gz"
]

df_list = []

print("-" * 80)
for filename in filenames:
    print(f"Reading {filename}")
    df = pd.read_table(filename)
    print(f"Read {len(df):,} rows with {len(df.TRID.unique()):,} unique TRIDs from {filename} and {len(df.columns):,} columns:\n\t\t{', '.join(df.columns)}")
    df_list.append(df)
    print("-" * 80)

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

df = df_list[0]
df = df[df.Stdev.notna()].copy()

# Keep the motif exactly as it arrived, before the canonicalization below merges strand variants.
# It is what tells apart a row pair that collapsed together here from one that was already
# duplicated in the input table.
df["rawLongestPureSegmentMotif"] = df["longestPureSegmentMotif"]
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

# compute_canonical_motif folds strand and rotation together, so motifs that differ only by strand
# ("CT" and "AG", "A" and "T") or by strand plus a rotation ("ATC" and "ATG") become the same value
# above. Those arrive as separate rows, one per motif, so a locus whose alleles were not all called
# on the same strand is left with more than one row for the same TRID2.
# Collapse those rows only when they genuinely divide the locus's alleles between them, which is what
# their N_motif adding up to numCalledAlleles shows. N_motif counts the alleles whose longest pure
# segment had that motif, so in that case the minority strand is a small slice of the cohort and the
# row the most alleles support is the one to keep. When the rows do not add up, each is its own
# summary of the whole cohort and they simply disagree: their N_motif is often identical, which would
# leave the choice to input row order. Those are left for the conflicting-rows check below, which is
# the same treatment the rows that already shared a motif before canonicalization get.
alleles_per_locus = df_list[1].set_index("TRID")["numCalledAlleles"]

before = len(df)
df = df.sort_values("N_motif", ascending=False, kind="mergesort")
dominant_raw_motif = df.groupby("TRID2")["rawLongestPureSegmentMotif"].transform("first")
distinct_motif_rows = df.drop_duplicates(subset=["TRID2", "rawLongestPureSegmentMotif", "N_motif"])
divides_the_locus = (df.TRID2.map(distinct_motif_rows.groupby("TRID2")["N_motif"].sum())
                     == df.TRID.map(alleles_per_locus))
df = df[~divides_the_locus | (df.rawLongestPureSegmentMotif == dominant_raw_motif)].sort_index()
print(f"Kept {len(df):,} out of {before:,} ({100 * len(df) / before:.2f}%) rows after dropping "
      f"minority-strand rows from loci whose strand-split rows account for every called allele")

df = df.drop(columns=["rawLongestPureSegmentMotif"])

# A row that repeats another row in every column is the same summary of the same locus written
# twice, so keeping one copy loses nothing.
before = len(df)
df = df.drop_duplicates()
print(f"Kept {len(df):,} out of {before:,} ({100 * len(df) / before:.2f}%) rows after dropping rows "
      f"that repeat another row exactly")

# Whatever still shares a TRID2 is one locus and one motif summarized more than once with conflicting
# numbers, and nothing here can tell which summary is right. These rows do not divide the locus's
# alleles between them either: summing their N_motif overshoots numCalledAlleles by a whole cohort,
# so the row with the larger N_motif is not the more complete one, just a differently wrong one.
# Drop the locus rather than pick a row, so that it ends up with no All of Us statistics at all
# instead of arbitrary ones. Fixing this for real needs a corrected input table.
ambiguous = df.TRID2.duplicated(keep=False)
if ambiguous.any():
    examples = sorted(df.TRID2[ambiguous].unique())[:3]
    print(f"WARNING: dropping {df.TRID2[ambiguous].nunique():,} loci ({ambiguous.sum():,} rows) that "
          f"are summarized more than once with conflicting values, e.g. {examples}")
    df = df[~ambiguous]

# Every consumer of this table keys on TRID2, so it has to identify a row on its own.
if df.TRID2.duplicated().any():
    raise ValueError(f"TRID2 is not unique: {df.TRID2.duplicated().sum():,} rows share one")

# join on the full TRID column since the other 2 tables only have 1 row per full TRID
df_list[0] = df

for df in df_list:
    df.set_index("TRID", inplace=True)

print(f"Joining dataframes: {filenames[0]} and {filenames[1]}")
df = df_list[0].join(df_list[1], how="left")
df.reset_index(inplace=True)
print(f"Resulting dataframe has {len(df):,} rows with {len(df.TRID.unique()):,} unique TRIDs and {len(df.columns):,} columns")

print(f"Joining dataframes: df and {filenames[2]}")
df.set_index("TRID", inplace=True)
df = df.join(df_list[2], how="left")
df.reset_index(inplace=True)
print(f"Resulting dataframe has {len(df):,} rows with {len(df.TRID.unique()):,} unique TRIDs and {len(df.columns):,} columns")

df["OE_len_percentile"] = df["OE_len"].rank(pct=True)

df_grouped_by_motif = df.groupby("canonical_motif")
df['StdevRankByMotif'] = df_grouped_by_motif["Stdev"].rank(ascending=False, method='min').astype(int)
df['StdevRankTotalNumberByMotif'] = df_grouped_by_motif["TRID"].transform("count").astype(int)

output_filename = "AoULR_phase1_TRGT_Weisburd_v1_combined.txt.gz"
df.to_csv(output_filename, sep="\t", index=False, header=True)
print(f"Wrote {len(df):,} rows to {output_filename}")
