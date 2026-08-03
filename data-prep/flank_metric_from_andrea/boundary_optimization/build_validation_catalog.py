"""Build a small ExpansionHunter catalog with several boundary-definition variants of 5 hand-picked
loci, for a real-read Hail Batch validation run (see run_hail_validation.py). Each variant is a
separate LocusId so one EH invocation genotypes all of them in a single pass.

Variant choice per locus:
  - "current": the real, currently-annotated KnownDiseaseAssociatedLoci boundary (same LocusStructure
    string as the production catalog, including any IUPAC-degenerate repeat units).
  - "purity90": both boundaries extended outward by the largest amount that keeps full-region
    reference-repeat purity >= 0.90 (ExpansionHunter's own IRR-read purity cutoff; see
    ehunter/locus/IrrPairFinder.hh), per compute_boundary_profiles.py's threshold-crossing table,
    ROUNDED DOWN to a whole number of motif-length units so the extended LocusStructure's implied
    repeat count stays well-defined (an unaligned offset, e.g. +7bp of a 3bp motif, leaves the
    boundary mid-unit and produces a confusing genotype).
  - EP400 additionally gets "wider" and "narrower": the two real alternative definitions gnomAD
    reported choosing between (see 2026-07 gnomAD STR data update, EP400 section).
  - RUNX2 additionally gets BOTH "purity90_naive" and "purity90_corrected": RUNX2's motif is the
    IUPAC-degenerate "GCN" (matches KnownDiseaseAssociatedLoci's real LocusStructure). The naive
    purity calculation treats the wildcard 'N' position as an automatic match against ANY base --
    which inflates purity for any random downstream sequence by ~1/motif_length regardless of
    whether it's really repeat-like, so the naive threshold crossing overshoots (right side: 21bp /
    7 units). "purity90_corrected" excludes wildcard motif positions from both the numerator and
    denominator of the purity calculation (right side: 18bp / 6 units) -- see the report for the
    full worked comparison; this is a real limitation of the plain purity metric for degenerate
    motifs (~14% of the 73 known disease loci use one), not fixed in boundary_metrics.py itself
    (that would require reprocessing all 3 tiers), just demonstrated here as a documented caveat.

EP400 isn't a KnownDiseaseAssociatedLoci entry -- it's included as the motivating worked example
from the gnomAD blog post, computed by hand from hg38 (verified against boundary_metrics.py).
"""

import json
from pathlib import Path

OUTPUT_PATH = Path(__file__).parent / "data" / "validation_catalog.json"

# (locus_id, chrom, current_start_0based, current_end, repeat_unit_for_locus_structure, extra_variants)
# extra_variants: {variant_name: (start_0based, end)}
LOCI = [
    {
        "locus_id": "EP400",
        "chrom": "chr12",
        "motif": "CAG",
        "current": (132062524, 132062611),  # gnomAD's "wider" (=de-facto current) definition, 29xCAG
        "variants": {
            "narrower": (132062548, 132062611),  # gnomAD's alternative narrower definition, 21xCAG
            "purity90": (132062524 - 6, 132062611 + 6),  # full_region_purity >= 0.90 rule, unit-aligned (2 units/side)
        },
    },
    {
        "locus_id": "C9ORF72",
        "chrom": "chr9",
        "motif": "GGCCCC",
        "current": (27573528, 27573546),
        "variants": {"purity90": (27573528 - 6, 27573546 + 0)},  # unit-aligned: 1 unit left, 0 right
    },
    {
        "locus_id": "RFC1",
        "chrom": "chr4",
        "motif": "AARRG",  # matches the real catalog's IUPAC-degenerate LocusStructure, not the plain RepeatUnit field
        "current": (39348424, 39348479),
        "variants": {"purity90": (39348424 - 5, 39348479 + 10)},  # unit-aligned: 1 unit left, 2 units right
    },
    {
        "locus_id": "BEAN1",
        "chrom": "chr16",
        "motif": "TAAAA",
        "current": (66490398, 66490453),
        "variants": {"purity90": (66490398 - 10, 66490453 + 15)},  # unit-aligned: 2 units left, 3 units right
    },
    {
        "locus_id": "RUNX2",
        "chrom": "chr6",
        "motif": "GCN",
        "current": (45422750, 45422792),
        "variants": {
            "purity90_naive": (45422750 - 9, 45422792 + 21),  # wildcard 'N' auto-matches -- inflated, see docstring
            "purity90_corrected": (45422750 - 6, 45422792 + 18),  # wildcard positions excluded from purity calc
        },
    },
]


def make_entry(locus_id, chrom, start_0based, end, motif):
    return {
        "LocusId": locus_id,
        "ReferenceRegion": f"{chrom}:{start_0based}-{end}",
        "LocusStructure": f"({motif})*",
        "VariantType": "Repeat",
    }


def main():
    catalog = []
    for locus in LOCI:
        start, end = locus["current"]
        catalog.append(make_entry(f"{locus['locus_id']}_current", locus["chrom"], start, end, locus["motif"]))
        for variant_name, (vstart, vend) in locus["variants"].items():
            catalog.append(
                make_entry(f"{locus['locus_id']}_{variant_name}", locus["chrom"], vstart, vend, locus["motif"])
            )

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT_PATH, "w") as f:
        json.dump(catalog, f, indent=2)

    print(f"Wrote {len(catalog)} catalog entries ({len(LOCI)} loci x variants) to {OUTPUT_PATH}")
    for entry in catalog:
        print(f"  {entry['LocusId']:20s} {entry['ReferenceRegion']:30s} {entry['LocusStructure']}")


if __name__ == "__main__":
    main()
