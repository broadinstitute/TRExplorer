set -ex

# Inputs. PMC answers command-line downloads with a web page instead of the file, so download these in a browser and
# put them in this directory:
#   Garg et al. 2019 (https://pmc.ncbi.nlm.nih.gov/articles/PMC6801017/), Supplementary Tables 1 and 2:
#     https://pmc.ncbi.nlm.nih.gov/articles/instance/6801017/bin/NIHMS1054049-supplement-Supp_TableS1.xlsx
#     https://pmc.ncbi.nlm.nih.gov/articles/instance/6801017/bin/NIHMS1054049-supplement-Supp_TableS2.xlsx
#   Garg et al. 2020 (https://pmc.ncbi.nlm.nih.gov/articles/PMC7536611/), Table S15:
#     https://pmc.ncbi.nlm.nih.gov/articles/instance/7536611/bin/mmc16.xlsx

python3 generate_Garg_2019_and_2020_DMR_lookup_json.py --trexplorer-catalog ~/code/tandem-repeat-catalogs/results__2026-10-29/release_draft_2026-10-29/TRExplorer.repeat_catalog_v2.1.hg38.1_to_1000bp_motifs.EH.json.gz
