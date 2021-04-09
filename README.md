
# modMEMOs (Modified MEMOS)
Modified version of the previously published tool MEMOS (Mazrooei et al. 2019). 

This tool is designed to identify significant enrichment of mutations within transcription factor binding sites and flanking regions. modMEMOs can be run using either ChIP-seq or ATAC-seq data.

## Running modMEMOS
**Required data**

* Mutation calls SNVs (.vcf)
* Peaks of interest (from ChIP- or ATAC-seq or both) (.bed)
* Fasta reference file (.fa)
* Transcription factor's position weight matrix (.pfm)

**modMEMOS-wrapper.sh**: the main pipeline used to execute the modMEMOs analysis.

**shuffleTest.sh**: Will randomly shuffle transcription factor sites, overlap these regions with mutation calls and count the number of SNVs. This process is repeated ${NUM} times. 

**modMEMOS_ATAC_Plus_ChIP.sh**: A bash script used to prepare bed and fasta files required to run modMEMOs and to execute modMEMOS-wrapper.sh.


# HoSRE (HotSpot in Regulatory Elements) 
An algorithm to identify hotspots of mutations within regulatory elements against a background of global and local mutation rates.

## Running HoSRE 

**MkPromotersBed_hg19()**:Create promoters reference bed file.

**FlattenPromoters()**: Reduces the promoters from the reference genome into a collapsed/flattened bed file for each gene. Therefore, it looks across all promoters for a gene rather than each promoter individually.

**Required data**
* Peaks of interest (.bed)
* Promoter reference file (.bed)
* Mutation calls SNVs (.vcf) 

**AssembleGR()**: Convert bed files into Granges while maintaining the same style throughout the analysis.

**overlapTargetsAndCatalogue()**: Determine which promoters are overlap with ATAC-seq peaks.

**removeCodingMuts()**: Function to separate coding from non-coding SNVs.

**separateVcfByBed()**: Separate SNVs into non-coding targets in order to limit the search space for the sliding window.

**applySlidingWindow()**: Slide a window of Wkb across each promoter and calculate the number of mutations within that region.

**estimateBgMutRate()**: Estimate the background mutations rate.

**mutRatePvalue()**: Using the Binomial test calculate the P value for each mutation followed by FDR adjustment.

**RunMutRate_Merge.R**: The main pipeline used to execute the analysis.

# Authors
Rene Quevedo

Samah El Ghamrasni
