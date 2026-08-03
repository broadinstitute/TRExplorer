"""Build a validation catalog (baseline / reanchored_015 / reanchored_020, deduped) for the loci among
the top 20,000 (ranked by flank repetitiveness, data/top20k_flank_repetitive_extensions.tsv) that
actually get extended by the motif-reanchoring rule at threshold 0.15 or 0.20. Same catalog-entry
convention as build_reanchor_validation_catalog_3.py (round 3).

Usage:
    python3 build_flank_repetitive_validation_catalog.py
"""

import json
from pathlib import Path

import pandas as pd

INPUT_TSV = Path(__file__).parent / "data" / "top20k_flank_repetitive_extensions.tsv"
OUTPUT_PATH = Path(__file__).parent / "data" / "flank_repetitive_validation_catalog.json"


def make_entry(locus_id, chrom, start_0based, end, motif):
    return {
        "LocusId": locus_id,
        "ReferenceRegion": f"{chrom}:{start_0based}-{end}",
        "LocusStructure": f"({motif})*",
        "VariantType": "Repeat",
    }


def main():
    df = pd.read_csv(INPUT_TSV, sep="\t")
    extended = df[(df["left_ext_015"] > 0) | (df["right_ext_015"] > 0) | (df["left_ext_020"] > 0) | (df["right_ext_020"] > 0)]
    print(f"{len(extended)} / {len(df)} loci extended at 0.15 or 0.20")

    catalog = []
    seen_ids = set()
    summary_rows = []
    for row in extended.itertuples():
        chrom, start, end, motif = row.chrom, row.start_0based, row.end_1based, row.motif
        base_id = row.LocusId

        variants = {
            "baseline": (start, end),
            "reanchored_015": (start - row.left_ext_015, end + row.right_ext_015),
            "reanchored_020": (start - row.left_ext_020, end + row.right_ext_020),
        }
        region_to_labels = {}
        for label, region in variants.items():
            region_to_labels.setdefault(region, []).append(label)

        for region, labels in region_to_labels.items():
            entry_id = f"{base_id}__{'+'.join(labels)}"
            if entry_id in seen_ids:
                continue
            seen_ids.add(entry_id)
            vstart, vend = region
            catalog.append(make_entry(entry_id, chrom, vstart, vend, motif))
            summary_rows.append((row.LocusId, labels, chrom, vstart, vend, motif))

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT_PATH, "w") as f:
        json.dump(catalog, f, indent=2)

    print(f"wrote {len(catalog)} catalog entries ({len(extended)} loci x up to 3 variants, deduped) to {OUTPUT_PATH}")
    for locus_id, labels, chrom, vstart, vend, motif in summary_rows[:15]:
        print(f"  {locus_id:35s} {'+'.join(labels):25s} {chrom}:{vstart}-{vend}")
    print("  ...")


if __name__ == "__main__":
    main()
