"""Combines multiple AoU (All of Us) data tables into a single output table.

This script merges three input tables:
- AoULR_phase1_TRGT_Weisburd_v1_lpsStats.txt.gz
- AoULR_phase1_TRGT_Weisburd_v1_distinctAndTotalAlleles.txt.gz
- AoULR_phase1_TRGT_Weisburd_v1.0.1_TRConstraint.txt.gz

It filters to keep only rows whose longest pure segment motif matches the locus motif on the same
strand, and only rows that cover more than half of their locus's called alleles, computes additional
statistics, and outputs a combined table with OE_len_percentile and StdevRankByMotif columns.
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

def compute_same_strand_motif(motif):
    """Returns the motif's alphabetically first rotation, without considering its reverse complement."""
    return compute_canonical_motif(motif, include_reverse_complement=False)

def get_individual_repeat_id_that_matches_longest_pure_segment_motif(row):
    trid = row["TRID"]
    if "," not in trid:
        return trid

    trid_parts = trid.split(",")
    trid2 = None
    for trid_part in trid_parts:
        if trid_part.count("-") != 3:
            return None
        if compute_same_strand_motif(trid_part.split("-")[3]) == row["sameStrandLongestPureSegmentMotif"]:
            if trid2 is None:
                trid2 = trid_part
            else:
                return None

    return trid2

df = df_list[0]
df = df[df.Stdev.notna()].copy()

# The input's longestPureSegmentMotif is read from the plus strand, so a row whose motif is only the
# reverse complement of a locus's motif measures a different, neighbouring run: at 1-6695162-6695173-T,
# the "A" row describes the 19 bp A run right before the T run. Rows are therefore matched to loci on
# the same strand (rotations only), while the output keeps the strand-merged canonical motif.
df["sameStrandLongestPureSegmentMotif"] = df["longestPureSegmentMotif"].apply(compute_same_strand_motif)
df["longestPureSegmentMotif"] = df["longestPureSegmentMotif"].apply(compute_canonical_motif)
df["TRID2"] = df.apply(get_individual_repeat_id_that_matches_longest_pure_segment_motif, axis=1)


before = len(df)
df = df[df.TRID2.notna()]
print(f"Kept {len(df):,} out of {before:,} ({100 * len(df) / before:.2f}%) rows where TRID2 is not None")

df["motif"] = df["TRID2"].str.split("-").str[3]
df["canonical_motif"] = df["motif"].apply(compute_canonical_motif)

before = len(df)
df = df[df.motif.apply(compute_same_strand_motif) == df.sameStrandLongestPureSegmentMotif]
print(f"Kept {len(df):,} out of {before:,} ({100 * len(df) / before:.2f}%) rows whose longest pure segment motif "
      f"matches the locus motif on the same strand")
df = df.drop(columns=["sameStrandLongestPureSegmentMotif"])

# Each input row summarizes only the N_motif alleles whose longest pure segment had its motif, so a
# locus can have several rows, e.g. one for a single allele with an odd pure run hundreds of bp long.
# Keep a row only if it covers more than half of the locus's called alleles, since only then does it
# describe the cohort. Rows that split a locus's alleles between them can't both pass. The input also
# sometimes summarizes the same locus more than once for the whole cohort (its rows' N_motif add up to
# more than numCalledAlleles). When only one of those summaries has a row covering more than half of
# the alleles, that row is kept: in the cases checked, it is the summary of the locus itself, while the
# other summary is dominated by a longer neighbouring run. When several do, the conflicting-rows check
# below drops the locus. A locus with no such row, e.g. a homopolymer inside a variation cluster whose
# longest pure segment is another locus's motif in most alleles, gets no All of Us statistics rather
# than a minority's.
alleles_per_locus = df_list[1].set_index("TRID")["numCalledAlleles"]
if not alleles_per_locus.index.is_unique:
    raise ValueError(f"{filenames[1]} has more than one row for some TRIDs")
if not df.TRID.isin(alleles_per_locus.index).all():
    raise ValueError(f"{(~df.TRID.isin(alleles_per_locus.index)).sum():,} rows have a TRID that is not in {filenames[1]}")

n_motif_sum_per_trid = df_list[0].drop_duplicates().groupby("TRID")["N_motif"].sum()
trids_summarized_more_than_once = n_motif_sum_per_trid.index[
    n_motif_sum_per_trid > alleles_per_locus.reindex(n_motif_sum_per_trid.index)]

before = len(df)
df = df[df.N_motif * 2 > df.TRID.map(alleles_per_locus)]
print(f"Kept {len(df):,} out of {before:,} ({100 * len(df) / before:.2f}%) rows that cover more than half "
      f"of their locus's called alleles")

# A row that repeats another row in every column is the same summary of the same locus written
# twice, so keeping one copy loses nothing.
before = len(df)
df = df.drop_duplicates()
print(f"Kept {len(df):,} out of {before:,} ({100 * len(df) / before:.2f}%) rows after dropping rows "
      f"that repeat another row exactly")

# Whatever still shares a TRID2 is one locus summarized more than once for the whole cohort, with a row
# covering more than half of the alleles in each summary and conflicting numbers (e.g. the same N_motif
# but a different Stdev), and nothing here can tell which summary is right. Drop the locus rather than
# pick a row, so that it ends up with no All of Us statistics at all instead of arbitrary ones. Fixing
# this for real needs a corrected input table.
ambiguous = df.TRID2.duplicated(keep=False)
if ambiguous.any():
    examples = sorted(df.TRID2[ambiguous].unique())[:3]
    print(f"WARNING: dropping {df.TRID2[ambiguous].nunique():,} loci ({ambiguous.sum():,} rows) that "
          f"are summarized more than once with conflicting values, e.g. {examples}")
    df = df[~ambiguous]

print(f"Kept the only summary with a row covering more than half of the called alleles at "
      f"{df.TRID.isin(trids_summarized_more_than_once).sum():,} loci that the input summarizes more than once")

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
