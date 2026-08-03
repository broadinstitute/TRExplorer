"""Build a validation catalog comparing threshold=0.15 vs threshold=0.20 for the motif-reanchoring
rule, across the 98 TRExplorer v1 sample loci that get nonzero extension under either threshold
(from threshold_sweep_trexplorer_v1.py). Each locus gets up to 3 entries (deduped): baseline,
reanchored_015, reanchored_020.

NOTE: this writes reanchor_validation_catalog_3.json, NOT the reanchor_validation_catalog_3_male.json
that run_hail_validation_3.py and analyze_round3_threshold_sweep.py actually consume -- that file is a
hand-augmented copy (adds 4 ARX_1/ARX_2 baseline/reanchored *_MALE_FIX entries for the corrected --sex
male round-3 rerun) that no script in this repo produces. Re-running this builder does NOT update the
file the round-3 pipeline actually reads.
"""

import json
from pathlib import Path

import pandas as pd

INPUT_TSV = Path(__file__).parent / "data" / "threshold_sweep_trexplorer_v1.tsv"
OUTPUT_PATH = Path(__file__).parent / "data" / "reanchor_validation_catalog_3.json"


def make_entry(locus_id, chrom, start_0based, end, motif):
    return {
        "LocusId": locus_id,
        "ReferenceRegion": f"{chrom}:{start_0based}-{end}",
        "LocusStructure": f"({motif})*",
        "VariantType": "Repeat",
    }


def main():
    df = pd.read_csv(INPUT_TSV, sep="\t")
    print(f"loaded {len(df)} candidate loci")

    catalog = []
    seen_ids = set()
    summary_rows = []
    for row in df.itertuples():
        chrom, start, end, motif = row.chrom, row.start_0based, row.end, row.motif
        base_id = row.locus_id.replace(":", "_")  # sanitize for LocusId use

        variants = {
            "baseline": (start, end),
            "reanchored_015": (start - row.left_ext_015, end + row.right_ext_015),
            "reanchored_020": (start - row.left_ext_020, end + row.right_ext_020),
        }
        # dedup identical regions per locus (e.g. 0.15 and 0.20 agree, or 0.15 == baseline)
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
            summary_rows.append((row.locus_id, labels, chrom, vstart, vend, motif))

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT_PATH, "w") as f:
        json.dump(catalog, f, indent=2)

    print(f"wrote {len(catalog)} catalog entries ({len(df)} loci x up to 3 variants, deduped) to {OUTPUT_PATH}")
    for locus_id, labels, chrom, vstart, vend, motif in summary_rows[:20]:
        print(f"  {locus_id:50s} {'+'.join(labels):25s} {chrom}:{vstart}-{vend}")
    print("  ...")


if __name__ == "__main__":
    main()
