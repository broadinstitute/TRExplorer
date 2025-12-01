"""
Filter the big file of Talkowski lab non-coding annotations (ncAnnot.v0.14.jul2024.bed.gz) described in 
https://docs.google.com/spreadsheets/d/1y8KxzNVqRlZSTcImkNLxRMuO8xDy33Xdz-KHHT67_VQ/edit?gid=1468174182#gid=1468174182


Categories to keep:

                category      size    size_Mb
                  zooCHAR     19689   0.019689
        gencode.v39.tRNAs     47738   0.047738
                   zooHAR     49173   0.049173
  gencode.v39.gene.snoRNA    108990   0.108990
                      uce    126007   0.126007
   gencode.v39.gene.miRNA    147935   0.147935
   gencode.v39.gene.snRNA    208089   0.208089
    fetal_cortex_sil_gw18    751680   0.751680
                      har    852699   0.852699
              uce_primate   1009099   1.009099
                     ucne   1415142   1.415142
                    vista   1495209   1.495209
                  en3_pLS   2162675   2.162675
      en3_dnase_k4m3_ctcf   2255591   2.255591
       cFootprint_primate   3607629   3.607629
            TF_footprints   3974490   3.974490
           en3_dnase_k4m3   4521752   4.521752
    fetal_cortex_enh_gw18   6538860   6.538860
             en3_pLS_ctcf   7797549   7.797549
                   PROcap   8294693   8.294693
           abc_loeuf_0.15   8712151   8.712151
         abc_pTriplo_0.84  75492825  75.492825
           abc_pHaplo_0.9  82798783  82.798783


Keep most annotations that span <=12MB + abc_loeuf_0.15 + abc_pTriplo_0.84 +  abc_pHaplo_0.9 (based on discussions with Hope)
"""

import argparse
import collections
import gzip
import os
import tqdm


parser = argparse.ArgumentParser()
parser.add_argument("--non-coding-annotations-bed", default="ncAnnot.v0.14.jul2024.bed.gz")
args = parser.parse_args()

output_path = args.non_coding_annotations_bed.replace(".bed.gz", "") + ".filtered.bed.gz"
print(f"Filtering {args.non_coding_annotations_bed} to {output_path}")

categories_to_keep = {
    "gencode.v39.tRNAs",
    "gencode.v39.gene.snoRNA",
    "gencode.v39.gene.miRNA",
    "gencode.v39.gene.snRNA",

    "zooCHAR",
    "zooHAR",
    "har",
    "uce",
    "uce_primate",
    "ucne",

    "cFootprint_primate",
    "TF_footprints",
    "zoonomia_tfbs",
    "UCSC_TF_ChIP",

    "vista",
    "fetal_cortex_sil_gw18",
    "fetal_cortex_enh_gw18",

    "en3_dnase_k4m3",
    "en3_pLS",
    "en3_dnase_k4m3_ctcf",
    "en3_pLS_ctcf",
    "PROcap",

    "abc_loeuf_0.15",
    "abc_pTriplo_0.84",
    "abc_pHaplo_0.9",
}


with gzip.open(args.non_coding_annotations_bed, "rt") as f:
    output_records = []
    counter = 0
    category_item_counters = collections.Counter()
    category_span_counters = collections.Counter()
    for line in tqdm.tqdm(f, unit=" lines", unit_scale=True, total=85_823_915):
        counter += 1
        fields = line.rstrip().split("\t")
        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        category = fields[3]
        category_item_counters[category] += 1
        category_span_counters[category] += end - start
        if category in categories_to_keep:
            output_records.append((chrom, start, end, category))
    
    print(f"Kept {len(output_records):,d} out of {counter:,d} ({100 * len(output_records) / counter:.2f}%) annotations")

print("Category counters:")
for category, count in sorted(category_item_counters.items(), key=lambda x: -x[1]):
    print(f" {count:10,d}   {category}")

print("Category size counters:")
for category, size in sorted(category_span_counters.items(), key=lambda x: -x[1]):
    print(f" {size:10,d} MB   {category}")

with open(output_path, "wt") as f:
    for chrom, start, end, category in sorted(output_records):
        f.write(f"{chrom}\t{start}\t{end}\t{category}\n")

os.system(f"bgzip -f {output_path}")
os.system(f"tabix -f {output_path}")

spanning_size = sum(category_span_counters[category] for category in categories_to_keep)
print(f"Wrote {len(output_records):,d} annotations to {output_path} spanning {spanning_size  / 1024 / 1024:.2f} MB:")
for category_to_keep in sorted(categories_to_keep, key=lambda x: -category_span_counters[x]):
    print(f" {category_item_counters[category_to_keep]:10,d} intervals  spanning {category_span_counters[category_to_keep] / 1024 / 1024:5.2f} MB  ({category_span_counters[category_to_keep] / spanning_size:5.1%}) {category_to_keep} ")



"""
Filtered results, sorted by span in genome (in MB):
     81,533 intervals  spanning 78.96 MB  (39.0%) abc_pHaplo_0.9
     73,698 intervals  spanning 72.00 MB  (35.5%) abc_pTriplo_0.84
      7,594 intervals  spanning  8.31 MB  ( 4.1%) abc_loeuf_0.15
     11,191 intervals  spanning  7.91 MB  ( 3.9%) PROcap
     27,134 intervals  spanning  7.44 MB  ( 3.7%) en3_pLS_ctcf
     24,218 intervals  spanning  6.24 MB  ( 3.1%) fetal_cortex_enh_gw18
     16,733 intervals  spanning  4.31 MB  ( 2.1%) en3_dnase_k4m3
    209,348 intervals  spanning  3.79 MB  ( 1.9%) TF_footprints
    208,280 intervals  spanning  3.44 MB  ( 1.7%) cFootprint_primate
      8,793 intervals  spanning  2.15 MB  ( 1.1%) en3_dnase_k4m3_ctcf
      7,578 intervals  spanning  2.06 MB  ( 1.0%) en3_pLS
      4,242 intervals  spanning  1.43 MB  ( 0.7%) vista
      4,351 intervals  spanning  1.35 MB  ( 0.7%) ucne
     33,368 intervals  spanning  0.96 MB  ( 0.5%) uce_primate
      3,167 intervals  spanning  0.81 MB  ( 0.4%) har
      2,784 intervals  spanning  0.72 MB  ( 0.4%) fetal_cortex_sil_gw18
      1,901 intervals  spanning  0.20 MB  ( 0.1%) gencode.v39.gene.snRNA
      1,822 intervals  spanning  0.14 MB  ( 0.1%) gencode.v39.gene.miRNA
        481 intervals  spanning  0.12 MB  ( 0.1%) uce
        926 intervals  spanning  0.10 MB  ( 0.1%) gencode.v39.gene.snoRNA
        312 intervals  spanning  0.05 MB  ( 0.0%) zooHAR
        648 intervals  spanning  0.05 MB  ( 0.0%) gencode.v39.tRNAs
        141 intervals  spanning  0.02 MB  ( 0.0%) zooCHAR
"""


"""
All categories:

zooCHAR                                 0.02
gencode.v39.tRNAs                       0.05
zooHAR                                  0.05
gencode.v39.gene.snoRNA                 0.11
gencode.v39.exon.snoRNA                 0.11
fantom5_ubiq_cells                      0.11
uce                                     0.13
fantom5_ubiq_organs                     0.13
gencode.v39.gene.miRNA                  0.15
gencode.v39.exon.miRNA                  0.15
gencode.v39.gene.snRNA                  0.21
gencode.v39.exon.snRNA                  0.21
gencode.v39.gene.miscRNA                0.45
gencode.v39.exon.miscRNA                0.45
fantom5_clust                           0.66
fetal_cortex_sil_gw18                   0.75
har                                     0.85
uce_primate                             1.01
fantom5_all_organs                      1.24
fantom5_spec_organs                     1.24
ucne                                    1.42
vista                                   1.50
en3_pLS                                 2.16
en3_dnase_k4m3_ctcf                     2.26
cFootprint_primate                      3.61
TF_footprints                           3.97
en3_dnase_k4m3                          4.52
chromHMM_insul                          6.19
fetal_cortex_enh_gw18                   6.54
fantom5_spec_cells                      6.55
fantom5_all_cells                       6.55
en3_pLS_ctcf                            7.80
PROcap                                  8.29
abc_loeuf_0.15                          8.71
fantom5_robust                         11.45
fantom5_perm                           12.38
atac-seq_pREs                          12.51
gencode.v39.exon.pseudogene            13.37
en3_ctcf                               14.32
cFootprint_mamm                        14.90
en3_pELS                               17.43
en4_bival                              17.62
en4_CTCF                               17.70
fantom5_ncz                            18.64
en3_pELS_ctcf                          19.17
en4_zinc                               23.80
zoonomia_group3                        25.53
zoonomia_tfbs                          27.07
en4_prom                               29.01
zoonomia_group2                        30.47
z_4.0                                  31.41
gencode.v39.exon.lncRNA                33.91
UCSC_TF_ChIP                           36.78
en4_prom_flank                         37.39
chromHMM_prom                          41.52
cDHS_primate                           53.98
en3_dELS_ctcf                          58.33
gencode.v39.UTR                        58.76
OCRs                                   60.51
gencode.v39.gene.pseudogene            64.04
TAD                                    73.05
abc_pTriplo_0.84                       75.49
phyloP_2                               81.91
abc_pHaplo_0.9                         82.80
en4_enh_low                            96.10
devcisreg                              98.60
en4_enh                               104.35
gencode.v39.exon.protein_coding       106.28
hg38_DNase                            116.95
zoonomia_group1                       120.55
en3_dELS                              127.32
gencode.v39.exon                      153.27
fetal_brain_atac_6_13gw               164.76
gencode.v39.intron.lncRNA             166.40
cDHS_mamm                             189.12
z_2.18                                200.99
abc_0.015                             231.91
chromHMM_polycomb                     283.86
phastCons_0.2                         341.78
chromHMM_enh                          353.03
en4_tx                                379.82
gencode.v39.gene.lncRNA               517.60
gencode.v39.intron.protein_coding     911.87
genedesert                           1039.81
gencode.v39.gene.protein_coding      1300.53
gencode.v39.integrated               1767.04
gencode.v39.annotation.genes         1767.86
"""