# modMEMOs (Modified MEMOS)
Modified version of the previously published tool MEMOS (Mazrooei et al. 2019). 

# CRE_MutEnrich
An algorithm to identify hotspots of mutations within regulatory elements against a background of global and local mutation rates.

## Running CRE_MutEnrich 

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


