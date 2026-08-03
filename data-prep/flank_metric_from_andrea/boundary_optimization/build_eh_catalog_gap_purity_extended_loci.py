"""Build an ExpansionHunter variant catalog covering both the original and gap-purity-rule-extended
definition of every locus (among Andrea's catalog loci whose original definition exactly matches
TRExplorer v2) that the rule actually extends -- data/gap_purity_scan_v2_exact_matches.tsv, filtered
to left_ext>0 or right_ext>0 (16,264 loci -> 32,528 entries, one pair per locus).

LocusId convention: "{source locus's own stable LocusId}__{original|extended}", e.g.
"6-45422750-45422792-GCN__original" / "6-45422750-45422792-GCN__extended" for RUNX2 -- keyed to the
ORIGINAL locus's own coordinates (not the extended region's own coordinates) so the two variants of
one locus always share a stable base id, and two different loci can never collide by both extending
to the same resulting region.

Usage:
    python3 build_eh_catalog_gap_purity_extended_loci.py
"""

import json
from pathlib import Path

import pandas as pd

DATA_DIR = Path(__file__).parent / "data"
INPUT_TSV = DATA_DIR / "gap_purity_scan_v2_exact_matches.tsv"
OUTPUT_PATH = DATA_DIR / "boundary_optimization_gap_purity_extended_loci.EH_catalog.json"


def make_entry(base_locus_id, chrom, start, end, motif, label):
    return {
        "LocusId": f"{base_locus_id}__{label}",
        "ReferenceRegion": f"{chrom}:{start}-{end}",
        "LocusStructure": f"({motif})*",
        "VariantType": "Repeat",
    }


def main():
    df = pd.read_csv(INPUT_TSV, sep="\t")
    extended = df[(df["left_ext"] > 0) | (df["right_ext"] > 0)]
    print(f"{len(extended):,} extended loci")

    catalog = []
    for row in extended.itertuples():
        catalog.append(make_entry(row.LocusId, row.chrom, row.start_0based, row.end_1based, row.motif, "original"))
        new_start = row.start_0based - row.left_ext
        new_end = row.end_1based + row.right_ext
        catalog.append(make_entry(row.LocusId, row.chrom, new_start, new_end, row.motif, "extended"))

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT_PATH, "w") as f:
        json.dump(catalog, f, indent=2)

    print(f"wrote {len(catalog):,} catalog entries to {OUTPUT_PATH}")
    for entry in catalog[:4]:
        print(f"  {entry['LocusId']:40s} {entry['ReferenceRegion']}")
    print("  ...")


if __name__ == "__main__":
    main()
