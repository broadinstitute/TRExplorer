# ENCODE4 and IGVF columns for TRExplorer loci

This folder computes a column for every TRExplorer v2.1 catalog locus from functional assays in the final
phase of ENCODE (ENCODE4, described in "The Encyclopedia of DNA Elements", bioRxiv 2026,
doi:10.64898/2026.07.06.731365) and in its successor consortium IGVF (https://data.igvf.org/). The output
is one TSV keyed by `LocusId`, in the same form as the SpliceAI summary TSV in `../splicing_prediction/`.
`../../bigquery-proxy/load_bigquery_main_table.py` loads it by default (`--encode4-and-igvf-validated-elements-tsv`).
A few loci in the TSV (17 in v2.1) are not in the BigQuery catalog JSON, so they are not loaded.

## Column

| Column | Meaning |
|---|---|
| `ENCODE4_and_IGVF_CRISPR_or_MPRA_validated_elements` | `data source:assay:biosample type` triples in which an element called functional overlaps the locus, e.g. `ENCODE:MPRA:HepG2,IGVF:CRISPR:HCT116` |

"Overlap" means at least 1 bp in common. All intervals are 0-based half-open. A biosample type is the
portal's biosample term name, such as `K562`. Several files from the same biosample type count once.

**An empty value does not mean the locus was tested.** The column only lists positive calls. Most loci
were never tested by any of these assays.

### What counts as "validated"

Each call follows a rule from the source portal or the producing lab:

- **ENCODE CRISPR:** rows with `Significant=TRUE` in the GRCh38 CRISPR "element quantifications" TSVs.
  Most of these rows are in the harmonized K562 dataset used to train ENCODE-rE2G
  ([ENCSR998YDI](https://www.encodeproject.org/annotations/ENCSR998YDI/)). The other files are from WTC11.
- **ENCODE MPRA:** elements with log2 fold change > 0 and FDR < 0.05 in the GRCh38 MPRA "element
  enrichments" BEDs. Files are skipped with a warning when they don't have that format's 11 columns, or
  when they have no FDR. 9 of the 53 files have no FDR: 8 have the placeholder -1 in every q-value, and
  ENCFF230JYM has 0. These 9 include all three WTC11 MPRA files, so there are no `ENCODE:MPRA:WTC11` calls.
- **IGVF CRISPR:** the 46 released "differential element quantifications" files come in about 17 layouts.
  Each layout is called with the rule its lab published:

  | Layout | Lab | Rule |
  |---|---|---|
  | `Significant` / `significant` column (Flow-FISH, TAP-seq) | Engreitz | Rows marked TRUE |
  | `growth_significant` / `mig_significant` (proliferation, migration screens) | Gersbach | Rows marked TRUE |
  | `FRACTEL_pval_fdr_corr` (Perturb-seq, FACS screens) | Gersbach | FDR < 0.05, as in the FRACTEL paper (bioRxiv 2026) |
  | `adj.pval.EnhancerEffect.noAux` without a `Significant` column (cohesin Flow-FISH) | Engreitz | FDR < 0.05 with cohesin present, as in Guckelberger et al. (bioRxiv 2024) |
  | `minus_auxin_padj` with `element_hg38` (cohesin TAP-seq) | Engreitz | FDR < 0.05 with cohesin present, as above |

  8 files are skipped with a warning:
  - 3 Perturb-seq files with only `p_val_adj`, and 1 FACS file with `FDR`/`Padj_Mean`: no published rule
    was found for these.
  - 2 Huangfu lab FACS files: the lab's rule (IDR ≤ 0.001 and a Z-score of at least 3 in each
    replicate) needs per-replicate Z-scores that the files don't have.
  - 1 file with no significance column.
  - 1 TAP-seq file with only hg19 coordinates.

  Rows without coordinates (negative controls such as `nc_MYC`, and gene-level targets such as
  `NA:NA-NA_ABL1`) are also dropped.

**Known gap:** most ENCODE CRISPR screens (Flow-FISH, FACS, proliferation; about 190 files) publish only
per-guide counts, with no per-element calls, so they are not used. Using them would mean re-running a
screen analysis.

**Coverage (2026-10-04):**

| Data | Loci with a call |
|---|---|
| Any | 118,894 (2.1%) |
| ENCODE MPRA | 116,790 |
| ENCODE CRISPR | 794 |
| IGVF CRISPR | 1,481 |

## ENCODE4 data that was tried and dropped

Each of these annotated too many loci to tell them apart, or didn't show that the TR itself is functional:

- **cCRE classes:** they overlapped 19.4% of loci, while cCREs cover 20.2% of the genome. So TRs are not
  enriched in them.
- **Number of biosample types with an intact Hi-C loop anchor at the locus:** this was non-zero for 82%
  of loci, because the anchors are 1 to 10 kb wide.
- **ENCODE-rE2G predicted target genes:** these covered 30.8% of loci (29.7% without links from a gene's
  own promoter to that gene).
- **STARR-seq peaks (32 files):** these overlapped 4.0% of loci, more than MPRA and CRISPR combined, but
  agreed poorly with MPRA (7,353 loci in both). A STARR-seq peak is hundreds of bp wide, and its fragments
  are tested on a plasmid with whatever repeat length the cell line carries. So a peak doesn't show that
  the TR itself is active.

## Caveats

- **Absence is not evidence of inactivity:** MPRA and CRISPR tested only chosen elements, mostly in a few
  cell lines. Most loci were never tested.
- **The TR may not be the active part:** an element is often a few hundred bp, so the activity can come
  from the sequence around the repeat.

## Running

Run every command from this folder.

```bash
# 1. Download the ENCODE and IGVF files (~100 files, 0.17 GB) and write encode4_source_files_manifest.tsv.
#    Files already present with the right md5 are skipped, so this can be rerun after an interruption.
python3 download_encode4_source_files.py --dry-run     # print file counts and sizes only
python3 download_encode4_source_files.py

# 2. Compute the column for every catalog locus. Prints its peak memory use after each step.
python3 compute_encode4_columns_for_catalog_loci.py \
    --catalog-bed ../splicing_prediction/TRExplorer.repeat_catalog_v2.1.hg38.1_to_1000bp_motifs.bed.gz

# Tests
python3 -m unittest compute_encode4_columns_for_catalog_loci_tests
```

The catalog BED comes from
https://github.com/broadinstitute/tandem-repeat-catalogs/releases/download/v2.1/TRExplorer.repeat_catalog_v2.1.hg38.1_to_1000bp_motifs.bed.gz

**Output:** `TRExplorer_v2.1_ENCODE4_columns.tsv.gz`, with one row per locus that has a value.

**Memory:** the catalog is held as numpy arrays rather than interval trees, and the ENCODE files are
read one at a time.
