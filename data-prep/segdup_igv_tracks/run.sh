set -ex

# Inputs:
#   UCSC's hg38 genomicSuperDups table (unchanged by UCSC since 2014)
curl -L -o genomicSuperDups.txt.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz

#   The hg38 segdups of Vollger et al. 2022 (https://doi.org/10.1126/science.abj6965), from annotation_bed_files.zip in
#   the paper's Zenodo record (https://zenodo.org/records/5502036). The zip is 492 MB and Zenodo can be slow and break
#   off downloads, so curl retries and resumes them.
if [ ! -f hg38.chr_only.SDs.bed.gz ]; then
    curl -L --retry 20 --retry-all-errors -C - -o annotation_bed_files.zip "https://zenodo.org/api/records/5502036/files/annotation_bed_files.zip/content"
    unzip -o -j annotation_bed_files.zip annotation_bed_files/hg38.chr_only/hg38.chr_only.SDs.bed.gz
    rm annotation_bed_files.zip
fi

python3 generate_segdup_igv_track_bed_files.py

gsutil cp UCSC_genomicSuperDups.hg38.bed.gz UCSC_genomicSuperDups.hg38.bed.gz.tbi Vollger2022_SDs.hg38.bed.gz Vollger2022_SDs.hg38.bed.gz.tbi gs://tgg-viewer/ref/GRCh38/segdups/
