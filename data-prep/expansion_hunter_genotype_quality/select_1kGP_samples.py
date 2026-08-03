"""Select up to 3 1kGP samples per (Population, Gender) combination, prioritizing samples that have an
OpenHGL (HPRC) de novo assembly (https://github.com/lh3/OpenHGL) -- this sort applies to every group,
but only visibly matters (changes which samples are kept) for groups with more available samples than
SAMPLES_PER_GROUP.

Inputs:
  ~/code/str-truth-set-v2/20130606_sample_info_1kGP.tsv                      -- 1kGP phase-3 sample metadata
  ~/code/str-truth-set-v2/run_tools/broad_short_read_cram_paths_for_all_1kGP.txt  -- the 698-sample short-read CRAM list
  https://raw.githubusercontent.com/lh3/OpenHGL/master/human579.meta.tsv     -- OpenHGL assembly metadata

Output:
  selected_1kGP_samples.tsv (in this script's own directory)
"""
import os
import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SAMPLE_INFO_PATH = os.path.expanduser("~/code/str-truth-set-v2/20130606_sample_info_1kGP.tsv")
CRAM_LIST_PATH = os.path.expanduser("~/code/str-truth-set-v2/run_tools/broad_short_read_cram_paths_for_all_1kGP.txt")
OPENHGL_META_URL = "https://raw.githubusercontent.com/lh3/OpenHGL/master/human579.meta.tsv"
OUTPUT_PATH = os.path.join(SCRIPT_DIR, "selected_1kGP_samples.tsv")

SAMPLES_PER_GROUP = 3

openhgl = pd.read_table(OPENHGL_META_URL, header=None,
    names=["assembly", "sex_chrom", "sample_name", "sex", "sgdp_region", "kg_pop", "country"])
openhgl_1kg_samples = set(openhgl.loc[openhgl["kg_pop"] != ".", "sample_name"].unique())

meta = pd.read_table(SAMPLE_INFO_PATH)
crams = pd.read_table(CRAM_LIST_PATH)
crams = crams.dropna(subset=["sample_id"])

merged = crams.merge(meta, left_on="sample_id", right_on="Sample", how="left")
merged = merged.dropna(subset=["Population", "Gender"])  # drops samples absent from the 1kGP metadata sheet

merged["has_t2t_assembly"] = merged["sample_id"].isin(openhgl_1kg_samples)

# has_t2t_assembly=True sorts first within each (Population, Gender) group; sample_id breaks ties deterministically
merged = merged.sort_values(["Population", "Gender", "has_t2t_assembly", "sample_id"], ascending=[True, True, False, True])
selected = merged.groupby(["Population", "Gender"], group_keys=False).head(SAMPLES_PER_GROUP)

out_cols = ["sample_id", "Population", "Population Description", "Gender", "has_t2t_assembly", "cram_path", "crai_path"]
selected[out_cols].to_csv(OUTPUT_PATH, sep="\t", index=False)

print(f"selected {len(selected)} of {len(merged)} samples ({selected['has_t2t_assembly'].sum()} have OpenHGL assemblies)")
print(f"wrote {OUTPUT_PATH}")
