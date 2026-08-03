"""Build a validation catalog for the motif-reanchoring extension rule (motif_reanchor_extension.py),
applied to all known-disease loci that (a) get nonzero extension from the rule and (b) are NOT
compound/multi-region loci (ReferenceRegion is a list in the source catalog) -- compound loci are
excluded because a "flank" extension there risks just re-absorbing the locus's OWN other sub-region
(as documented for PRNP in the original report), not genuinely new repeat-like sequence, which would
make the before/after comparison meaningless. Each locus gets 2 entries: baseline (current annotated
definition) and "reanchored" (both sides extended per the rule, using its exact ReferenceRegion /
RepeatUnit as in KnownDiseaseAssociatedLoci_July2024.json).

RUNX2, C9ORF72, RFC1, BEAN1 are excluded here even though tested before, since C9ORF72/RFC1/BEAN1
get ZERO extension from this rule (nothing new to test) and RUNX2 already has this exact "reanchored"
region tested in validation_catalog.json under a different variant name. EP400 is also excluded --
already fully covered (both its baseline and reanchored/"wider" definitions were tested previously).
"""

import json
import sys
from pathlib import Path

import pysam

sys.path.insert(0, str(Path(__file__).parent))
from motif_reanchor_extension import compute_candidate_definitions

CATALOG_PATH = "/Users/weisburd/code/tandem-repeat-explorer/catalogs/data/KnownDiseaseAssociatedLoci_July2024.json"
FASTA_PATH = "/Users/weisburd/hg38.fa"
OUTPUT_PATH = Path(__file__).parent / "data" / "reanchor_validation_catalog_2.json"

# Loci with nonzero extension from the rule, excluding compound-region loci (ATXN8OS, FRA10AC1, HTT,
# TMEM185A) and loci already covered by the first validation round (RUNX2, EP400, C9ORF72, RFC1, BEAN1).
SELECTED_LOCUS_IDS = [
    "HOXA13_3", "NIPA1", "PABPN1", "ARX_1", "ATXN3",
    "GIPC1", "NOTCH2NLA", "ZIC2", "ARX_2", "HOXA13_1", "LRP12",
    "RUNX2",  # exact reanchored region (right+9bp only) wasn't among the variants tested in round 1
]


def make_entry(locus_id, chrom, start_0based, end, motif):
    return {
        "LocusId": locus_id,
        "ReferenceRegion": f"{chrom}:{start_0based}-{end}",
        "LocusStructure": f"({motif})*",
        "VariantType": "Repeat",
    }


def main():
    fasta = pysam.FastaFile(FASTA_PATH)
    known_loci = {locus["LocusId"]: locus for locus in json.load(open(CATALOG_PATH))}

    catalog = []
    summary_rows = []
    for locus_id in SELECTED_LOCUS_IDS:
        locus = known_loci[locus_id]
        chrom, coords = locus["ReferenceRegion"].split(":")
        start, end = map(int, coords.split("-"))
        motif = locus["RepeatUnit"]

        result = compute_candidate_definitions(chrom, start, end, motif, fasta)
        new_start = start - result["left_extension"]
        new_end = end + result["right_extension"]

        catalog.append(make_entry(f"{locus_id}_baseline", chrom, start, end, motif))
        catalog.append(make_entry(f"{locus_id}_reanchored", chrom, new_start, new_end, motif))

        summary_rows.append(
            (locus_id, chrom, start, end, motif, result["left_extension"], result["right_extension"], new_start, new_end)
        )

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT_PATH, "w") as f:
        json.dump(catalog, f, indent=2)

    print(f"Wrote {len(catalog)} catalog entries ({len(SELECTED_LOCUS_IDS)} loci x 2) to {OUTPUT_PATH}\n")
    print(f"{'LocusId':12s} {'motif':6s} {'baseline':30s} {'left_ext':9s} {'right_ext':9s} {'reanchored':30s}")
    for locus_id, chrom, start, end, motif, left_ext, right_ext, new_start, new_end in summary_rows:
        print(
            f"{locus_id:12s} {motif:6s} {chrom}:{start}-{end:<20d} {left_ext:9d} {right_ext:9d} "
            f"{chrom}:{new_start}-{new_end}"
        )


if __name__ == "__main__":
    main()
