import pandas as pd
from str_analysis.utils.canonical_repeat_unit import compute_canonical_motif

filenames = [
    "AoULR_phase1_TRGT_Weisburd_v1.0.1_lpsStats.txt.gz",
    "AoULR_phase1_TRGT_Weisburd_v1.0.1_distinctAndTotalAlleles.txt.gz",
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


output_filename = "AoULR_phase1_TRGT_Weisburd_v1.0.1_combined.txt.gz"
df.to_csv(output_filename, sep="\t", index=False, header=True)
print(f"Wrote {len(df):,} rows to {output_filename}")



