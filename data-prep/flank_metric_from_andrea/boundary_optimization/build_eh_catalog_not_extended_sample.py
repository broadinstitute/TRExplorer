"""Build an ExpansionHunter catalog of a random sample of 15,000 loci (among Andrea's catalog loci
whose original definition exactly matches TRExplorer v2) that the gap-purity rule does NOT extend
(left_ext==0 and right_ext==0) -- a comparison population for the "extended vs. not extended"
absolute-pOk violin plots, sampled from data/gap_purity_scan_v2_exact_matches.tsv with a fixed seed
for reproducibility.

Usage:
    python3 build_eh_catalog_not_extended_sample.py
"""

import json
from pathlib import Path

import pandas as pd

DATA_DIR = Path(__file__).parent / "data"
INPUT_TSV = DATA_DIR / "gap_purity_scan_v2_exact_matches.tsv"
SAMPLE_TSV = DATA_DIR / "not_extended_sample_15000.tsv"
OUTPUT_PATH = DATA_DIR / "not_extended_sample_15000.EH_catalog.json"
N_SAMPLE = 15_000
RANDOM_SEED = 42


def main():
    df = pd.read_csv(INPUT_TSV, sep="\t")
    not_extended = df[(df["left_ext"] == 0) & (df["right_ext"] == 0)]
    print(f"{len(not_extended):,} not-extended loci in the v2-exact-match pool")

    sample = not_extended.sample(n=N_SAMPLE, random_state=RANDOM_SEED)
    sample.to_csv(SAMPLE_TSV, sep="\t", index=False)

    catalog = []
    for row in sample.itertuples():
        chrom_no_prefix = row.chrom.removeprefix("chr")
        locus_id = f"{chrom_no_prefix}-{row.start_0based}-{row.end_1based}-{row.motif}__not_extended"
        catalog.append({
            "LocusId": locus_id,
            "ReferenceRegion": f"{row.chrom}:{row.start_0based}-{row.end_1based}",
            "LocusStructure": f"({row.motif})*",
            "VariantType": "Repeat",
        })

    with open(OUTPUT_PATH, "w") as f:
        json.dump(catalog, f, indent=2)
    print(f"wrote {len(catalog):,} catalog entries to {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
