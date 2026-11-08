set -ex

# Inputs: the supplementary tables of Tanudisastro et al. 2026 (https://doi.org/10.1038/s41586-026-10438-2).
# The zip holds them in a folder named 2025-03-08020C-s2, so -j drops that folder and puts the files where the script
# looks for them.
curl -L -o 41586_2026_10438_MOESM3_ESM.zip "https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-026-10438-2/MediaObjects/41586_2026_10438_MOESM3_ESM.zip"
unzip -o -j 41586_2026_10438_MOESM3_ESM.zip -d 41586_2026_10438_MOESM3_ESM

python3 generate_Tanudisastro_2024_lookup_json.py --catalog-path ~/code/tandem-repeat-catalogs/results__2026-10-29/release_draft_2026-10-29/TRExplorer.repeat_catalog_v2.1.hg38.1_to_1000bp_motifs.EH.json.gz
