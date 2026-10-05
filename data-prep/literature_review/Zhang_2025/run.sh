set -ex

# Inputs. NCBI's FTP server returned empty responses on 2026-10-04; if that happens, retry later.
# The design files are pinned to str_mpra_design commit 5157fe4bcafe98c7d6b5385891a94d843d2e76a0.
curl -L -o GSE306816_hSTR1_HEK293_linear_regression.csv.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE306nnn/GSE306816/suppl/GSE306816_hSTR1_HEK293_linear_regression.csv.gz
curl -L -o tss_str_pairs.tsv https://raw.githubusercontent.com/gymreklab/str_mpra_design/5157fe4bcafe98c7d6b5385891a94d843d2e76a0/design/tss_str_pairs.tsv
curl -L -o array_probes.tsv https://raw.githubusercontent.com/gymreklab/str_mpra_design/5157fe4bcafe98c7d6b5385891a94d843d2e76a0/design/array_probes.tsv

python3 generate_Zhang_2025_lookup_json.py --trexplorer-catalog ~/code/tandem-repeat-catalogs/results__2026-08-30/release_draft_2026-08-30/TRExplorer.repeat_catalog_v2.1.hg38.1_to_1000bp_motifs.EH.json.gz
