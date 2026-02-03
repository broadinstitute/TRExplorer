"""Global constants shared between BigQuery data loading and website generation.

This module contains the BigQuery schema definition with column descriptions that serve
as the single source of truth for column documentation across the project.
"""

# Gene region priority description used in multiple column descriptions
GENE_REGION_PRIORITY = (
    "with significance defined as coding > 5' UTR > 3' UTR > exon in non-coding gene > "
    "intron > promoter region spanning 1,000bp upstream of the 5' UTR"
)

# Export group names for organizing columns in the export dialog
EXPORT_GROUP_CORE = "Core"
EXPORT_GROUP_ADDITIONAL_LAYERS = "Additional Layers"
EXPORT_GROUP_GENE_ANNOTATIONS = "Gene Annotations"
EXPORT_GROUP_POLYMORPHISM_HPRC256 = "Polymorphism (HPRC256)"
EXPORT_GROUP_POLYMORPHISM_AOU1027 = "Polymorphism (AoU1027)"
EXPORT_GROUP_POLYMORPHISM_TENK10K = "Polymorphism (TenK10K)"
EXPORT_GROUP_VAMOS = "Vamos"
EXPORT_GROUP_SC_ETRS = "sc-eTRs (Tanudisastro 2024)"
EXPORT_GROUP_PHEWAS = "PheWAS (Manigbas 2024)"


BIGQUERY_COLUMNS = [
    # Core locus identifiers
    {
        "type": "INTEGER",
        "name": "id",
        "mode": "REQUIRED",
        "description": "Auto-incrementing unique identifier for each row in the database.",
    },
    {
        "type": "INTEGER",
        "name": "chrom_index",
        "mode": "REQUIRED",
        "description": "Numeric chromosome index (1-22, 23=X, 24=Y, 25=MT) used for partitioning and sorting.",
    },
    {
        "type": "STRING",
        "name": "chrom",
        "mode": "REQUIRED",
        "description": "Chromosome name without 'chr' prefix (e.g., '1', '2', 'X', 'Y', 'MT').",
        "displayName": "Chromosome",
        "group": EXPORT_GROUP_CORE,
        "exportColumn": "CONCAT('chr', chrom) AS chrom",
    },
    {
        "type": "INTEGER",
        "name": "start_0based",
        "mode": "REQUIRED",
        "description": "0-based start coordinate of the TR locus in hg38.",
        "displayName": "Start (0-based)",
        "group": EXPORT_GROUP_CORE,
    },
    {
        "type": "INTEGER",
        "name": "end_1based",
        "mode": "REQUIRED",
        "description": "1-based end coordinate of the TR locus in hg38 (half-open interval).",
        "displayName": "End (1-based)",
        "group": EXPORT_GROUP_CORE,
    },
    {
        "type": "STRING",
        "name": "ReferenceRegion",
        "mode": "REQUIRED",
        "description": "Genomic coordinates of the TR locus in the hg38 reference genome. Example: chr1:94418421-94418442",
        "displayName": "Reference Region",
        "group": EXPORT_GROUP_CORE,
    },
    {
        "type": "STRING",
        "name": "LocusId",
        "description": "A unique identifier for this TR locus in hg38, formatted as chrom-start-end-motif.",
        "displayName": "Locus ID",
        "group": EXPORT_GROUP_CORE,
    },
    {
        "type": "INTEGER",
        "name": "ReferenceRegionSize",
        "description": "Size of the TR locus in base pairs (end_1based - start_0based).",
        "displayName": "Locus Size (bp)",
        "group": EXPORT_GROUP_CORE,
    },
    {
        "type": "INTEGER",
        "name": "MotifSize",
        "description": "Length of the repeat motif in base pairs.",
        "displayName": "Motif Size",
        "group": EXPORT_GROUP_CORE,
    },
    {
        "type": "STRING",
        "name": "ReferenceMotif",
        "description": "The repeat motif sequence as it appears in hg38.",
        "displayName": "Motif",
        "group": EXPORT_GROUP_CORE,
    },
    {
        "type": "STRING",
        "name": "CanonicalMotif",
        "description": (
            "The repeat motif in hg38, normalized by computing all cyclic rotations of the motif "
            "and its reverse complement, and then selecting the one that is alphabetically first. "
            "For example, the canonical version of CTG is AGC."
        ),
        "displayName": "Canonical Motif",
        "group": EXPORT_GROUP_CORE,
    },
    {
        "type": "INTEGER",
        "name": "NumRepeatsInReference",
        "description": (
            "Number of repeats calculated by taking the width of the TR locus reference interval "
            "and dividing it by the motif size, then rounding down to the nearest integer."
        ),
        "displayName": "# Repeats in Reference",
        "group": EXPORT_GROUP_CORE,
    },
    {
        "type": "FLOAT",
        "name": "ReferenceRepeatPurity",
        "description": (
            "The fraction of bases within the TR locus reference sequence that match perfect repeats "
            "of the specified motif. For example, CAG.CAG.CAG has purity = 1.0, while CAG.CAA.CAG "
            "has purity = 0.89. Range: 0 to 1."
        ),
        "displayName": "Repeat Purity",
        "group": EXPORT_GROUP_CORE,
    },

    # Highest purity motif fields
    {
        "type": "STRING",
        "name": "HighestPurityMotif",
        "description": "The motif that achieves the highest purity score at this locus (may differ from ReferenceMotif).",
    },
    {
        "type": "FLOAT",
        "name": "HighestPurityMotifPurity",
        "description": "The purity score achieved by the HighestPurityMotif. Range: 0 to 1.",
    },
    {
        "type": "FLOAT",
        "name": "HighestPurityMotifQuality",
        "description": "Quality score for the HighestPurityMotif assignment.",
    },

    # Locus quality and context
    {
        "type": "INTEGER",
        "name": "NsInFlanks",
        "description": "Number of N bases within the flanking regions of this locus. Used to filter loci for tools like ExpansionHunter that don't support Ns in flanks.",
        "displayName": "Ns in Flanks",
        "group": EXPORT_GROUP_CORE,
    },
    {
        "type": "INTEGER",
        "name": "TRsInRegion",
        "description": "Number of TR loci in the vicinity of this locus that are separated from each other by no more than 6bp of spacer sequence.",
        "displayName": "Nearby TRs",
        "group": EXPORT_GROUP_CORE,
    },
    {
        "type": "STRING",
        "name": "Source",
        "description": "TRExplorer catalog version and the source catalog that contributed this TR locus to the TRExplorer catalog.",
        "displayName": "Source",
        "group": EXPORT_GROUP_CORE,
    },

    # Non-overlapping locus flags
    {
        "type": "INTEGER",
        "name": "NonOverlappingLongestLocus",
        "description": "Whether this is the longest locus among overlapping loci at this position (1 = yes, 0 = no).",
        "displayName": "Non-Overlapping Longest",
        "group": EXPORT_GROUP_ADDITIONAL_LAYERS,
    },
    {
        "type": "INTEGER",
        "name": "NonOverlappingPurestLocus",
        "description": "Whether this is the purest locus among overlapping loci at this position (1 = yes, 0 = no).",
        "displayName": "Non-Overlapping Purest",
        "group": EXPORT_GROUP_ADDITIONAL_LAYERS,
    },

    # Source catalog membership flags
    {
        "type": "STRING",
        "name": "FoundInKnownDiseaseAssociatedLoci",
        "description": "Whether this locus was found in the known disease-associated loci catalog.",
    },
    {
        "type": "STRING",
        "name": "FoundInIllumina174kPolymorphicTRs",
        "description": "Whether this locus was found in the Illumina 174k polymorphic TR catalog.",
    },
    {
        "type": "STRING",
        "name": "FoundInPerfectRepeatsInReference",
        "description": "Whether this locus was found in the perfect repeats in reference catalog.",
    },
    {
        "type": "STRING",
        "name": "FoundInPolymorphicTRsInT2TAssemblies",
        "description": "Whether this locus was found in the polymorphic TRs in T2T assemblies catalog.",
    },
    {
        "type": "INTEGER",
        "name": "FoundInGangSTRCatalog",
        "description": "Whether this locus was found in the GangSTR v17 catalog (1 = yes, NULL = no).",
    },

    # Mappability
    {
        "type": "FLOAT",
        "name": "LeftFlankMappability",
        "description": "UCSC 36-mer mappability score for the left flanking region. Range: 0 to 1.",
    },
    {
        "type": "FLOAT",
        "name": "RightFlankMappability",
        "description": "UCSC 36-mer mappability score for the right flanking region. Range: 0 to 1.",
    },
    {
        "type": "FLOAT",
        "name": "FlanksAndLocusMappability",
        "description": "UCSC 36-mer mappability scores averaged across the TR locus and +/- 150bp of its flanking sequences. Range: 0 to 1.",
        "displayName": "Mappability",
        "group": EXPORT_GROUP_CORE,
    },

    # Variation clusters
    {
        "type": "STRING",
        "name": "VariationCluster",
        "description": (
            "Genomic interval of the variation cluster containing this locus, if any. "
            "Variation Clusters (VCs) are genomic intervals that contain one or more TRs "
            "embedded within a wider region that harbors multiple common polymorphisms."
        ),
        "displayName": "Variation Cluster",
        "group": EXPORT_GROUP_ADDITIONAL_LAYERS,
    },
    {
        "type": "INTEGER",
        "name": "VariationClusterSizeDiff",
        "description": "Size difference in bp between the variation cluster interval and the reference repeat interval.",
        "displayName": "VC Size Diff",
        "group": EXPORT_GROUP_ADDITIONAL_LAYERS,
    },
    {
        "type": "STRING",
        "name": "VariationClusterFilterReason",
        "description": "Reason this locus was filtered from the variation cluster analysis, if applicable.",
    },

    # Disease association
    {
        "type": "STRING",
        "name": "KnownDiseaseAssociatedLocus",
        "description": "Name of the known disease-associated locus if this TR overlaps one.",
        "displayName": "Known Disease Locus",
        "group": EXPORT_GROUP_ADDITIONAL_LAYERS,
    },
    {
        "type": "STRING",
        "name": "KnownDiseaseAssociatedMotif",
        "description": "The pathogenic repeat motif at the known disease-associated locus.",
        "displayName": "Known Disease Motif",
        "group": EXPORT_GROUP_ADDITIONAL_LAYERS,
    },
    {
        "type": "STRING",
        "name": "DiseaseInfo",
        "description": "JSON object containing detailed disease information including LocusId, STRchiveId, STRipyId, and disease names.",
    },

    # Gene annotations - Gencode
    {
        "type": "STRING",
        "name": "GencodeGeneRegion",
        "description": f"The most significant gene region that overlaps the TR locus in Gencode v48, {GENE_REGION_PRIORITY}.",
        "displayName": "Gencode Gene Region",
        "group": EXPORT_GROUP_GENE_ANNOTATIONS,
    },
    {
        "type": "STRING",
        "name": "GencodeGeneName",
        "description": "Gencode v48 gene name",
        "displayName": "Gencode Gene Name",
        "group": EXPORT_GROUP_GENE_ANNOTATIONS,
    },
    {
        "type": "STRING",
        "name": "GencodeGeneId",
        "description": "Gencode v48 gene ID (ENSG)",
        "displayName": "Gencode Gene ID",
        "group": EXPORT_GROUP_GENE_ANNOTATIONS,
    },
    {
        "type": "STRING",
        "name": "GencodeTranscriptId",
        "description": "Gencode v48 transcript ID (ENST)",
        "displayName": "Gencode Transcript ID",
        "group": EXPORT_GROUP_GENE_ANNOTATIONS,
    },

    # Gene annotations - RefSeq
    {
        "type": "STRING",
        "name": "RefseqGeneRegion",
        "description": f"The most significant RefSeq gene region that overlaps the TR locus, {GENE_REGION_PRIORITY}.",
        "displayName": "RefSeq Gene Region",
        "group": EXPORT_GROUP_GENE_ANNOTATIONS,
    },
    {
        "type": "STRING",
        "name": "RefseqGeneName",
        "description": "RefSeq gene name",
        "displayName": "RefSeq Gene Name",
        "group": EXPORT_GROUP_GENE_ANNOTATIONS,
    },
    {
        "type": "STRING",
        "name": "RefseqGeneId",
        "description": "RefSeq gene ID",
        "displayName": "RefSeq Gene ID",
        "group": EXPORT_GROUP_GENE_ANNOTATIONS,
    },
    {
        "type": "STRING",
        "name": "RefseqTranscriptId",
        "description": "RefSeq transcript ID",
        "displayName": "RefSeq Transcript ID",
        "group": EXPORT_GROUP_GENE_ANNOTATIONS,
    },

    # Gene annotations - MANE
    {
        "type": "STRING",
        "name": "ManeGeneRegion",
        "description": f"The most significant gene region that overlaps the TR locus in MANE v1.4, {GENE_REGION_PRIORITY}.",
        "displayName": "MANE Gene Region",
        "group": EXPORT_GROUP_GENE_ANNOTATIONS,
    },
    {
        "type": "STRING",
        "name": "ManeGeneName",
        "description": "MANE v1.4 gene name",
        "displayName": "MANE Gene Name",
        "group": EXPORT_GROUP_GENE_ANNOTATIONS,
    },
    {
        "type": "STRING",
        "name": "ManeGeneId",
        "description": "MANE v1.4 gene ID",
        "displayName": "MANE Gene ID",
        "group": EXPORT_GROUP_GENE_ANNOTATIONS,
    },
    {
        "type": "STRING",
        "name": "ManeTranscriptId",
        "description": "MANE v1.4 transcript ID",
        "displayName": "MANE Transcript ID",
        "group": EXPORT_GROUP_GENE_ANNOTATIONS,
    },

    # Illumina174k / 1kGP data
    {
        "type": "STRING",
        "name": "AlleleFrequenciesFromIllumina174k",
        "description": "Allele frequencies based on ExpansionHunter calls in 2,504 short-read genomes from 1kGP.",
        "displayName": "1kGP Allele Frequencies",
        "group": EXPORT_GROUP_ADDITIONAL_LAYERS,
    },
    {
        "type": "FLOAT",
        "name": "StdevFromIllumina174k",
        "description": "Standard deviation of the allele frequency distribution based on ExpansionHunter calls in 2,504 short-read genomes from 1kGP.",
        "displayName": "1kGP Stdev",
        "group": EXPORT_GROUP_ADDITIONAL_LAYERS,
    },

    # T2T assemblies data
    {
        "type": "STRING",
        "name": "AlleleFrequenciesFromT2TAssemblies",
        "description": "Allele frequencies based on assembly-to-hg38 alignments of 78 diploid T2T assemblies.",
        "displayName": "T2T Allele Frequencies",
        "group": EXPORT_GROUP_ADDITIONAL_LAYERS,
    },
    {
        "type": "FLOAT",
        "name": "StdevFromT2TAssemblies",
        "description": "Standard deviation of the allele frequency distribution based on assembly-to-hg38 alignments of 78 diploid T2T assemblies.",
        "displayName": "T2T Stdev",
        "group": EXPORT_GROUP_ADDITIONAL_LAYERS,
    },

    # TenK10K data
    {
        "type": "STRING",
        "name": "TenK10K_AlleleHistogram",
        "description": "Allele frequencies based on ExpansionHunter calls in 1,925 short-read genomes from the TenK10K Phase 1 dataset.",
        "displayName": "TenK10K Allele Histogram",
        "group": EXPORT_GROUP_POLYMORPHISM_TENK10K,
    },
    {
        "type": "STRING",
        "name": "TenK10K_BiallelicHistogram",
        "description": "Biallelic genotype histogram from TenK10K Phase 1 short-read genomes.",
        "displayName": "TenK10K Biallelic Histogram",
        "group": EXPORT_GROUP_POLYMORPHISM_TENK10K,
    },
    {
        "type": "INTEGER",
        "name": "TenK10K_MinAllele",
        "description": "Minimum allele size observed in TenK10K Phase 1.",
        "displayName": "TenK10K Min Allele",
        "group": EXPORT_GROUP_POLYMORPHISM_TENK10K,
    },
    {
        "type": "INTEGER",
        "name": "TenK10K_ModeAllele",
        "description": "Mode (most common) allele size in TenK10K Phase 1.",
        "displayName": "TenK10K Mode Allele",
        "group": EXPORT_GROUP_POLYMORPHISM_TENK10K,
    },
    {
        "type": "FLOAT",
        "name": "TenK10K_Stdev",
        "description": "Standard deviation of allele sizes in TenK10K Phase 1.",
        "displayName": "TenK10K Stdev",
        "group": EXPORT_GROUP_POLYMORPHISM_TENK10K,
    },
    {
        "type": "FLOAT",
        "name": "TenK10K_Median",
        "description": "Median allele size in TenK10K Phase 1.",
        "displayName": "TenK10K Median",
        "group": EXPORT_GROUP_POLYMORPHISM_TENK10K,
    },
    {
        "type": "FLOAT",
        "name": "TenK10K_99thPercentile",
        "description": "99th percentile allele size in TenK10K Phase 1.",
        "displayName": "TenK10K 99th Percentile",
        "group": EXPORT_GROUP_POLYMORPHISM_TENK10K,
    },
    {
        "type": "INTEGER",
        "name": "TenK10K_MaxAllele",
        "description": "Maximum allele size observed in TenK10K Phase 1.",
        "displayName": "TenK10K Max Allele",
        "group": EXPORT_GROUP_POLYMORPHISM_TENK10K,
    },
    {
        "type": "INTEGER",
        "name": "TenK10K_UniqueAlleleLengths",
        "description": "Number of unique allele lengths observed in TenK10K Phase 1.",
        "displayName": "TenK10K Unique Allele Lengths",
        "group": EXPORT_GROUP_POLYMORPHISM_TENK10K,
    },
    {
        "type": "INTEGER",
        "name": "TenK10K_NumCalledAlleles",
        "description": "Total number of alleles called in TenK10K Phase 1.",
        "displayName": "TenK10K Num Called",
        "group": EXPORT_GROUP_POLYMORPHISM_TENK10K,
    },
    {
        "type": "INTEGER",
        "name": "TenK10K_StdevRankByMotif",
        "description": "Rank of this locus by standard deviation among all loci with the same motif in TenK10K Phase 1.",
        "displayName": "TenK10K Stdev Rank",
        "group": EXPORT_GROUP_POLYMORPHISM_TENK10K,
    },
    {
        "type": "INTEGER",
        "name": "TenK10K_StdevRankTotalNumberByMotif",
        "description": "Total number of loci with the same motif used for ranking in TenK10K Phase 1.",
        "displayName": "TenK10K Stdev Rank Total",
        "group": EXPORT_GROUP_POLYMORPHISM_TENK10K,
    },

    # HPRC256 data
    {
        "type": "STRING",
        "name": "HPRC256_AlleleHistogram",
        "description": "Allele frequencies based on TRGT calls in 256 PacBio HiFi samples from the Human Pan-genome Reference Consortium (HPRC). Allele sizes computed using the longest pure segment (LPS).",
        "displayName": "HPRC256 Allele Histogram",
        "group": EXPORT_GROUP_POLYMORPHISM_HPRC256,
    },
    {
        "type": "STRING",
        "name": "HPRC256_BiallelicHistogram",
        "description": "Biallelic genotype histogram from HPRC 256 PacBio samples.",
        "displayName": "HPRC256 Biallelic Histogram",
        "group": EXPORT_GROUP_POLYMORPHISM_HPRC256,
    },
    {
        "type": "INTEGER",
        "name": "HPRC256_MinAllele",
        "description": "Minimum allele size (LPS) observed in HPRC256.",
        "displayName": "HPRC256 Min Allele",
        "group": EXPORT_GROUP_POLYMORPHISM_HPRC256,
    },
    {
        "type": "INTEGER",
        "name": "HPRC256_ModeAllele",
        "description": "Mode (most common) allele size (LPS) in HPRC256.",
        "displayName": "HPRC256 Mode Allele",
        "group": EXPORT_GROUP_POLYMORPHISM_HPRC256,
    },
    {
        "type": "FLOAT",
        "name": "HPRC256_Stdev",
        "description": "Standard deviation of allele sizes (LPS) in HPRC256.",
        "displayName": "HPRC256 Stdev",
        "group": EXPORT_GROUP_POLYMORPHISM_HPRC256,
    },
    {
        "type": "FLOAT",
        "name": "HPRC256_Median",
        "description": "Median allele size (LPS) in HPRC256.",
        "displayName": "HPRC256 Median",
        "group": EXPORT_GROUP_POLYMORPHISM_HPRC256,
    },
    {
        "type": "FLOAT",
        "name": "HPRC256_99thPercentile",
        "description": "99th percentile allele size (LPS) in HPRC256.",
        "displayName": "HPRC256 99th Percentile",
        "group": EXPORT_GROUP_POLYMORPHISM_HPRC256,
    },
    {
        "type": "INTEGER",
        "name": "HPRC256_MaxAllele",
        "description": "Maximum allele size (LPS) observed in HPRC256.",
        "displayName": "HPRC256 Max Allele",
        "group": EXPORT_GROUP_POLYMORPHISM_HPRC256,
    },
    {
        "type": "INTEGER",
        "name": "HPRC256_UniqueAlleleLengths",
        "description": "Number of unique LPS allele lengths observed in HPRC256.",
        "displayName": "HPRC256 Unique Allele Lengths",
        "group": EXPORT_GROUP_POLYMORPHISM_HPRC256,
    },
    {
        "type": "INTEGER",
        "name": "HPRC256_NumCalledAlleles",
        "description": "Total number of alleles called in HPRC256.",
        "displayName": "HPRC256 Num Called",
        "group": EXPORT_GROUP_POLYMORPHISM_HPRC256,
    },
    {
        "type": "INTEGER",
        "name": "HPRC256_StdevRankByMotif",
        "description": "Rank of this locus by standard deviation among all loci with the same motif in HPRC256.",
        "displayName": "HPRC256 Stdev Rank",
        "group": EXPORT_GROUP_POLYMORPHISM_HPRC256,
    },
    {
        "type": "INTEGER",
        "name": "HPRC256_StdevRankTotalNumberByMotif",
        "description": "Total number of loci with the same motif used for ranking in HPRC256.",
        "displayName": "HPRC256 Stdev Rank Total",
        "group": EXPORT_GROUP_POLYMORPHISM_HPRC256,
    },

    # AoU1027 data
    {
        "type": "INTEGER",
        "name": "AoU1027_MinAllele",
        "description": "Minimum allele size (LPS) observed in AoU1027 (1,027 African American participants from All of Us).",
        "displayName": "AoU1027 Min Allele",
        "group": EXPORT_GROUP_POLYMORPHISM_AOU1027,
    },
    {
        "type": "INTEGER",
        "name": "AoU1027_ModeAllele",
        "description": "Mode (most common) allele size (LPS) in AoU1027.",
        "displayName": "AoU1027 Mode Allele",
        "group": EXPORT_GROUP_POLYMORPHISM_AOU1027,
    },
    {
        "type": "FLOAT",
        "name": "AoU1027_Stdev",
        "description": "Standard deviation of allele sizes (LPS) in AoU1027.",
        "displayName": "AoU1027 Stdev",
        "group": EXPORT_GROUP_POLYMORPHISM_AOU1027,
    },
    {
        "type": "FLOAT",
        "name": "AoU1027_Median",
        "description": "Median allele size (LPS) in AoU1027.",
        "displayName": "AoU1027 Median",
        "group": EXPORT_GROUP_POLYMORPHISM_AOU1027,
    },
    {
        "type": "FLOAT",
        "name": "AoU1027_99thPercentile",
        "description": "99th percentile allele size (LPS) in AoU1027.",
        "displayName": "AoU1027 99th Percentile",
        "group": EXPORT_GROUP_POLYMORPHISM_AOU1027,
    },
    {
        "type": "INTEGER",
        "name": "AoU1027_MaxAllele",
        "description": "Maximum allele size (LPS) observed in AoU1027.",
        "displayName": "AoU1027 Max Allele",
        "group": EXPORT_GROUP_POLYMORPHISM_AOU1027,
    },
    {
        "type": "INTEGER",
        "name": "AoU1027_NumCalledAlleles",
        "description": "Total number of alleles called in AoU1027.",
        "displayName": "AoU1027 Num Called",
        "group": EXPORT_GROUP_POLYMORPHISM_AOU1027,
    },
    {
        "type": "INTEGER",
        "name": "AoU1027_StdevRankByMotif",
        "description": "Rank of this locus by standard deviation among all loci with the same motif in AoU1027.",
        "displayName": "AoU1027 Stdev Rank",
        "group": EXPORT_GROUP_POLYMORPHISM_AOU1027,
    },
    {
        "type": "INTEGER",
        "name": "AoU1027_StdevRankTotalNumberByMotif",
        "description": "Total number of loci with the same motif used for ranking in AoU1027.",
        "displayName": "AoU1027 Stdev Rank Total",
        "group": EXPORT_GROUP_POLYMORPHISM_AOU1027,
    },
    {
        "type": "FLOAT",
        "name": "AoU1027_OE_Length",
        "description": "Observed / Expected TR length metric from the TR constraint model (Danzi et al. 2025).",
        "displayName": "Constraint: O/E Length",
        "group": EXPORT_GROUP_ADDITIONAL_LAYERS,
    },
    {
        "type": "FLOAT",
        "name": "AoU1027_OE_LengthPercentile",
        "description": "Observed / Expected TR length metric displayed as a percentile rank relative to all other TR loci in the catalog.",
        "displayName": "Constraint: O/E Length Percentile",
        "group": EXPORT_GROUP_ADDITIONAL_LAYERS,
    },

    # RepeatMasker and non-coding annotations
    {
        "type": "STRING",
        "name": "RepeatMaskerIntervals",
        "description": "UCSC Repeat Masker intervals that overlap the TR locus",
        "displayName": "RepeatMasker Intervals",
        "group": EXPORT_GROUP_ADDITIONAL_LAYERS,
    },
    {
        "type": "STRING",
        "name": "NonCodingAnnotations",
        "description": "Non-coding annotations from various sources collected by the Talkowski Lab",
        "displayName": "Non-coding Annotations",
        "group": EXPORT_GROUP_ADDITIONAL_LAYERS,
    },

    # Vamos data
    {
        "type": "STRING",
        "name": "VamosUniqueMotifs",
        "description": "Unique motifs detected at this locus by the Vamos v2.1 pipeline when applied to HPRC assemblies.",
        "displayName": "Vamos Unique Motifs",
        "group": EXPORT_GROUP_VAMOS,
    },
    {
        "type": "STRING",
        "name": "VamosEfficientMotifs",
        "description": "Efficient motif representation for Vamos genotyping.",
        "displayName": "Vamos Efficient Motifs",
        "group": EXPORT_GROUP_VAMOS,
    },
    {
        "type": "STRING",
        "name": "VamosMotifFrequencies",
        "description": "Different motifs and their frequencies detected at this locus by Vamos v2.1 in 94 HPRC assemblies.",
        "displayName": "Vamos Motif Frequencies",
        "group": EXPORT_GROUP_VAMOS,
    },
    {
        "type": "INTEGER",
        "name": "VamosNumUniqueMotifs",
        "description": "Number of unique motifs detected by Vamos at this locus.",
        "displayName": "Vamos Num Unique Motifs",
        "group": EXPORT_GROUP_VAMOS,
    },
    {
        "type": "INTEGER",
        "name": "IncludeInVamosCatalog",
        "description": "Whether this locus should be included in Vamos catalog exports (1 = yes). Excludes overlapping loci since Vamos doesn't support them.",
        "displayName": "Include In Vamos Catalog",
        "group": EXPORT_GROUP_VAMOS,
    },

    # Tanudisastro 2024 sc-eTR data
    {
        "type": "STRING",
        "name": "Tanudisastro2024_SignificantCellTypes",
        "description": "Comma-separated list of immune cell types with significant sc-eTR associations (p < 0.05) from Tanudisastro et al. 2024.",
        "displayName": "sc-eTRs: Cell Types",
        "group": EXPORT_GROUP_SC_ETRS,
    },
    {
        "type": "FLOAT",
        "name": "Tanudisastro2024_MinPvalue",
        "description": "Minimum (most significant) p-value across all cell types for this locus from Tanudisastro et al. 2024.",
        "displayName": "sc-eTRs: Min P-value",
        "group": EXPORT_GROUP_SC_ETRS,
    },
    {
        "type": "STRING",
        "name": "Tanudisastro2024_Details",
        "description": "JSON with per-cell-type sc-eTR details (p-value, effect size, eGene, rank) from Tanudisastro et al. 2024.",
    },

    # Manigbas 2024 PheWAS data
    {
        "type": "STRING",
        "name": "Manigbas2024_AssociatedTraits",
        "description": "Comma-separated list of human traits with fine-mapped causal associations from Manigbas et al. 2024 UK Biobank PheWAS. Empty string indicates locus was tested but had no significant associations.",
        "displayName": "PheWAS: Traits",
        "group": EXPORT_GROUP_PHEWAS,
    },
    {
        "type": "FLOAT",
        "name": "Manigbas2024_MinPvalue",
        "description": "Minimum (most significant) p-value across all associated traits from Manigbas et al. 2024. NULL for loci with no associations.",
        "displayName": "PheWAS: Min P-value",
        "group": EXPORT_GROUP_PHEWAS,
    },
    {
        "type": "STRING",
        "name": "Manigbas2024_Details",
        "description": "JSON with per-trait PheWAS statistics from Manigbas et al. 2024. Keys are trait names, values contain pvalue, log10Pvalue, traitType, sampleSize, caviarRank, caviarProbability, replicatedInAoU, manigbasLocusId.",
        "displayName": "PheWAS: Details",
        "group": EXPORT_GROUP_PHEWAS,
    },

    # Reference sequence data
    {
        "type": "STRING",
        "name": "ReferenceRepeatSequence",
        "description": "The repeat sequence in hg38",
    },
    {
        "type": "STRING",
        "name": "ReferenceLeftFlank",
        "description": "The sequence immediately to the left of the TR locus in hg38",
    },
    {
        "type": "STRING",
        "name": "ReferenceRightFlank",
        "description": "The sequence immediately to the right of the TR locus in hg38",
    },
]


def get_column_descriptions():
    """Return a dictionary mapping column names to their descriptions.

    Returns:
        dict: Mapping of column name to description string.
    """
    return {col["name"]: col.get("description", "") for col in BIGQUERY_COLUMNS}


def get_exportable_columns():
    """Return a list of columns that can be exported, with their display metadata.

    Only columns that have a 'displayName' and 'group' are considered exportable.

    Returns:
        list: List of dicts with keys: name, displayName, group, column, description
              where column is the SQL expression to use (defaults to name).
    """
    exportable = []
    for col in BIGQUERY_COLUMNS:
        if "displayName" in col and "group" in col:
            exportable.append({
                "name": col["name"],
                "displayName": col["displayName"],
                "group": col["group"],
                "column": col.get("exportColumn", col["name"]),
                "description": col.get("description", ""),
            })
    return exportable
