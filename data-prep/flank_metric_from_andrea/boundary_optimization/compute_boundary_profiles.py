"""Compute boundary-extension metric profiles for every locus in a locus-set TSV produced by
build_locus_sets.py, using boundary_metrics.compute_boundary_extension_profile.

Usage: python3 compute_boundary_profiles.py <tier_tsv> <output_prefix> [--max-offset N]
    [--no-edit-distance] [--local-window-size N]

Writes:
  <output_prefix>.profiles.pkl   -- {locus_id: {"left": {...curves...}, "right": {...curves...}, ...}}
  <output_prefix>.purity_crossings.tsv -- one row per locus x side x threshold, the first offset at
      which purity would drop below that threshold if the boundary were extended that far. This is
      the compact table the threshold-finding / plotting step actually uses.
"""

import argparse
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

sys.path.insert(0, str(Path(__file__).parent))
from boundary_metrics import compute_boundary_extension_profile

FASTA_PATH = "/Users/weisburd/hg38.fa"

PURITY_THRESHOLDS = [1.00, 0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.60, 0.50]


def best_phase_motif(core_seq, motif):
    """Rotate `motif` so that motif[0] is in the reading frame that best matches core_seq[0].

    TRExplorer's bed file does not guarantee the motif string starts in phase with the reference
    region's first base (unlike EH's own hand-curated LocusStructure catalogs, where it does by
    construction) -- so the phase is determined empirically here rather than assumed.
    """
    core_seq = core_seq.upper()
    motif = motif.upper()
    motif_length = len(motif)
    best_phase, best_matches = 0, -1
    for phase in range(motif_length):
        rotated = motif[phase:] + motif[:phase]
        matches = sum(
            1 for i, base in enumerate(core_seq) if base == rotated[i % motif_length] or rotated[i % motif_length] == "N"
        )
        if matches > best_matches:
            best_matches, best_phase = matches, phase
    rotated_motif = motif[best_phase:] + motif[:best_phase]
    core_purity = best_matches / max(1, len(core_seq))
    return rotated_motif, core_purity


def first_offset_below_threshold(purity_curve, threshold):
    """Smallest offset (1-indexed) at which purity_curve drops below `threshold`, else None."""
    below = np.where(purity_curve < threshold)[0]
    return int(below[0] + 1) if len(below) else None


def compute_profiles_for_tier(fasta, loci_df, max_offset, compute_edit_distance, local_window_size):
    profiles = {}
    crossing_rows = []
    low_core_purity_count = 0
    skipped_near_chrom_end = 0
    skipped_no_sequence = 0

    for row in loci_df.itertuples():
        if getattr(row, "too_close_to_chrom_end", False):
            skipped_near_chrom_end += 1
            continue

        chrom, start, end, motif = row.chrom, int(row.start_0based), int(row.end), row.motif
        core_seq = fasta.fetch(chrom, start, end)
        if not core_seq:
            skipped_no_sequence += 1
            continue

        rotated_motif, core_purity = best_phase_motif(core_seq, motif)
        if core_purity < 0.7:
            low_core_purity_count += 1

        flank_start = max(0, start - max_offset)
        left_flank = fasta.fetch(chrom, flank_start, start)
        right_flank = fasta.fetch(chrom, end, end + max_offset)
        this_max_offset = min(max_offset, len(left_flank), len(right_flank))
        if this_max_offset < 10:
            skipped_no_sequence += 1
            continue

        profile = compute_boundary_extension_profile(
            left_flank,
            core_seq,
            right_flank,
            rotated_motif,
            max_offset=this_max_offset,
            local_window_size=local_window_size,
            compute_edit_distance=compute_edit_distance,
        )
        profiles[row.locus_id] = {
            "chrom": chrom,
            "start_0based": start,
            "end": end,
            "motif": motif,
            "rotated_motif": rotated_motif,
            "core_purity": core_purity,
            "core_length": end - start,
            "tier": row.tier,
            "max_offset_used": this_max_offset,
            "profile": profile,
        }

        for side in ("left", "right"):
            full_purity = profile[side]["full_region_purity"]
            local_purity = profile[side]["purity"]
            for threshold in PURITY_THRESHOLDS:
                crossing_rows.append(
                    {
                        "locus_id": row.locus_id,
                        "tier": row.tier,
                        "motif_length": len(motif),
                        "side": side,
                        "core_length": end - start,
                        "core_purity": core_purity,
                        "threshold": threshold,
                        "full_region_first_offset_below": first_offset_below_threshold(full_purity, threshold),
                        "local_flank_first_offset_below": first_offset_below_threshold(local_purity, threshold),
                        "max_offset_used": this_max_offset,
                    }
                )

    print(
        f"  processed {len(profiles)} loci; skipped {skipped_near_chrom_end} (near chrom end), "
        f"{skipped_no_sequence} (no/short sequence); {low_core_purity_count} had best-phase core "
        f"purity < 0.7 (flagged, not dropped)"
    )
    return profiles, pd.DataFrame(crossing_rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("tier_tsv")
    parser.add_argument("output_prefix")
    parser.add_argument("--max-offset", type=int, default=300)
    parser.add_argument("--no-edit-distance", action="store_true")
    parser.add_argument("--local-window-size", type=int, default=24)
    args = parser.parse_args()

    loci_df = pd.read_csv(args.tier_tsv, sep="\t")
    print(f"loaded {len(loci_df)} loci from {args.tier_tsv}")

    fasta = pysam.FastaFile(FASTA_PATH)
    profiles, crossings_df = compute_profiles_for_tier(
        fasta,
        loci_df,
        max_offset=args.max_offset,
        compute_edit_distance=not args.no_edit_distance,
        local_window_size=args.local_window_size,
    )

    output_prefix = Path(args.output_prefix)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)

    with open(f"{output_prefix}.profiles.pkl", "wb") as pickle_file:
        pickle.dump(profiles, pickle_file)
    crossings_df.to_csv(f"{output_prefix}.purity_crossings.tsv", sep="\t", index=False)

    print(f"wrote {output_prefix}.profiles.pkl ({len(profiles)} loci)")
    print(f"wrote {output_prefix}.purity_crossings.tsv ({len(crossings_df)} rows)")


if __name__ == "__main__":
    main()
