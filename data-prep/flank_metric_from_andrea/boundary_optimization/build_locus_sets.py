"""Build tier2 (loci sharing a canonical motif with a known disease locus) and tier3 (stratified
sample of the full TRExplorer v1 catalog) locus sets, to compare against tier1 (the 73 known
disease-associated loci, already built directly from KnownDiseaseAssociatedLoci_July2024.json).

Usage: python3 build_locus_sets.py
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, "/Users/weisburd/code/str-analysis")
from str_analysis.utils.canonical_repeat_unit import compute_canonical_motif
from str_analysis.utils.fasta_utils import get_chromosome_sizes

TREXPLORER_V1_BED = "/Users/weisburd/code/tandem-repeat-explorer/catalogs/data/TR_Explorer.repeat_catalog_v1.hg38.1_to_1000bp_motifs.bed.gz"
TIER1_TSV = Path(__file__).parent / "data" / "tier1_known_disease_loci.tsv"
DATA_DIR = Path(__file__).parent / "data"
FASTA_PATH = "/Users/weisburd/hg38.fa"

MAX_PER_MOTIF = 150
TIER3_SAMPLE_SIZE = 20_000
MOTIF_LENGTH_BUCKETS = [(1, 2), (3, 3), (4, 4), (5, 5), (6, 6), (7, 10), (11, 20), (21, 50), (51, 1000)]
FLANK_MARGIN = 500  # matches compute_boundary_profiles.py's max --max-offset headroom
RNG_SEED = 42


def canonical_motif_cache():
    cache = {}

    def get(motif):
        if motif not in cache:
            cache[motif] = compute_canonical_motif(motif, include_reverse_complement=True)
        return cache[motif]

    return get


def load_trexplorer():
    print(f"Loading {TREXPLORER_V1_BED} ...")
    df = pd.read_csv(
        TREXPLORER_V1_BED, sep="\t", header=None, names=["chrom", "start", "end", "motif", "_col5"], compression="gzip"
    )
    df = df.drop(columns=["_col5"])
    print(f"  {len(df):,} rows loaded")
    return df


def add_canonical_motif(df, get_canonical):
    distinct_motifs = df["motif"].unique()
    print(f"  {len(distinct_motifs):,} distinct motif strings; computing canonical form for each (cached)...")
    motif_to_canonical = {m: get_canonical(m) for m in distinct_motifs}
    df["canonical_motif"] = df["motif"].map(motif_to_canonical)
    return df


def add_too_close_to_chrom_end(df, chrom_sizes):
    df["too_close_to_chrom_end"] = df.apply(
        lambda r: r["start_0based"] < FLANK_MARGIN
        or r["end"] > chrom_sizes.get(r["chrom"], float("inf")) - FLANK_MARGIN,
        axis=1,
    )
    return df


def build_tier2(trexplorer, tier1, get_canonical):
    tier1 = tier1.copy()
    tier1["canonical_motif"] = tier1["motif"].apply(get_canonical)

    tier1_locus_key = set(zip(tier1["chrom"], tier1["start_0based"], tier1["end"]))
    canonical_to_tier1_ids = tier1.groupby("canonical_motif")["locus_id"].apply(list).to_dict()

    matched = trexplorer[trexplorer["canonical_motif"].isin(canonical_to_tier1_ids)].copy()
    matched = matched[~matched.apply(lambda r: (r["chrom"], r["start"], r["end"]) in tier1_locus_key, axis=1)]
    print(f"  {len(matched):,} TRExplorer loci share a canonical motif with a tier1 locus (pre-cap)")

    rng = np.random.default_rng(RNG_SEED)
    parts = []
    for canonical_motif, group in matched.groupby("canonical_motif"):
        if len(group) > MAX_PER_MOTIF:
            group = group.iloc[rng.choice(len(group), size=MAX_PER_MOTIF, replace=False)]
        parts.append(group)
    tier2 = pd.concat(parts, ignore_index=True) if parts else matched.iloc[0:0].copy()

    def source_note(canonical_motif):
        ids = canonical_to_tier1_ids[canonical_motif]
        note = ",".join(ids[:5])
        if len(ids) > 5:
            note += f"...({len(ids)} total)"
        return f"matches {canonical_motif} via {note}"

    tier2["source_note"] = tier2["canonical_motif"].apply(source_note)
    tier2["locus_id"] = tier2.apply(lambda r: f"{r['chrom']}-{r['start']}-{r['end']}-{r['motif']}", axis=1)
    tier2["motif_length"] = tier2["motif"].str.len()
    tier2["tier"] = "same_canonical_motif"
    tier2 = tier2.rename(columns={"start": "start_0based"})
    print(f"  {len(tier2):,} tier2 loci after capping at {MAX_PER_MOTIF}/motif")
    return tier2[
        ["locus_id", "chrom", "start_0based", "end", "motif", "canonical_motif", "motif_length", "tier", "source_note"]
    ]


def build_tier3(trexplorer, tier1, tier2):
    exclude_keys = set(zip(tier1["chrom"], tier1["start_0based"], tier1["end"])) | set(
        zip(tier2["chrom"], tier2["start_0based"], tier2["end"])
    )
    df = trexplorer.copy()
    df["motif_length"] = df["motif"].str.len()

    def bucket_label(length):
        for lo, hi in MOTIF_LENGTH_BUCKETS:
            if lo <= length <= hi:
                return f"{lo}-{hi}"
        return "51-1000"

    df["length_bucket"] = df["motif_length"].apply(bucket_label)

    print("  motif-length-bucket proportions in the full catalog:")
    full_props = df["length_bucket"].value_counts(normalize=True).sort_index()
    print(full_props.to_string())

    rng = np.random.default_rng(RNG_SEED)
    parts = []
    for bucket, group in df.groupby("length_bucket"):
        target_n = round(TIER3_SAMPLE_SIZE * len(group) / len(df))
        target_n = min(target_n, len(group))
        sampled = group.iloc[rng.choice(len(group), size=target_n, replace=False)]
        parts.append(sampled)
    tier3 = pd.concat(parts, ignore_index=True)

    n_before_dedup = len(tier3)
    tier3 = tier3[~tier3.apply(lambda r: (r["chrom"], r["start"], r["end"]) in exclude_keys, axis=1)]
    print(f"  sampled {n_before_dedup:,} loci, {len(tier3):,} remain after excluding tier1/tier2 overlaps")

    print("  motif-length-bucket proportions in the sample:")
    print(tier3["length_bucket"].value_counts(normalize=True).sort_index().to_string())

    tier3 = tier3.rename(columns={"start": "start_0based"})
    tier3["locus_id"] = tier3.apply(lambda r: f"{r['chrom']}-{r['start_0based']}-{r['end']}-{r['motif']}", axis=1)
    tier3["tier"] = "trexplorer_v1_sample"
    tier3["source_note"] = "stratified random sample, seed=" + str(RNG_SEED)
    return tier3[
        ["locus_id", "chrom", "start_0based", "end", "motif", "canonical_motif", "motif_length", "tier", "source_note"]
    ]


def main():
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    get_canonical = canonical_motif_cache()

    tier1 = pd.read_csv(TIER1_TSV, sep="\t")

    trexplorer = load_trexplorer()
    trexplorer = add_canonical_motif(trexplorer, get_canonical)

    print("\nBuilding tier2 (same canonical motif as a known disease locus)...")
    tier2 = build_tier2(trexplorer, tier1, get_canonical)

    print("\nBuilding tier3 (stratified sample of all TRExplorer v1 loci)...")
    tier3 = build_tier3(trexplorer, tier1, tier2)

    print("\nFlagging loci too close to a chromosome end...")
    chrom_sizes = get_chromosome_sizes(FASTA_PATH)
    tier2 = add_too_close_to_chrom_end(tier2, chrom_sizes)
    tier3 = add_too_close_to_chrom_end(tier3, chrom_sizes)

    tier2.to_csv(DATA_DIR / "tier2_same_canonical_motif.tsv", sep="\t", index=False)
    tier3.to_csv(DATA_DIR / "tier3_trexplorer_v1_sample.tsv", sep="\t", index=False)

    print(f"\n=== Summary ===")
    print(f"tier1 (known disease loci): {len(tier1)} rows")
    print(f"tier2 (same canonical motif): {len(tier2)} rows, {int(tier2['too_close_to_chrom_end'].sum())} flagged near chrom end")
    print(f"tier3 (TRExplorer v1 sample): {len(tier3)} rows, {int(tier3['too_close_to_chrom_end'].sum())} flagged near chrom end")
    print(f"distinct tier1 canonical motifs: {tier1['motif'].apply(get_canonical).nunique()}")
    print("Wrote data/tier2_same_canonical_motif.tsv and data/tier3_trexplorer_v1_sample.tsv")


if __name__ == "__main__":
    main()
