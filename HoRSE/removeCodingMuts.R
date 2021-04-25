#' removeCodingMuts
#'
#' @param vcf A VariantAnnotation object
#' @param bed A GenomicRanges object of non-coding regions
#'
#' @return A list of VariantAnnotation objects for of 
#' coding and non-coding VCF mutations
#' @export
#'
#' @examples
removeCodingMuts <- function(vcf, bed){
  require(GenomicRanges)
  require(VariantAnnotation)
  require(TxDb.Hsapiens.UCSC.hg19.knownGene)
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene 
  
  seqlevels(vcf) <- as.character(unique(seqnames(vcf)))
  seqlevelsStyle(vcf) <- seqlevelsStyle(bed) <- 'UCSC'
  ov.idx <- findOverlaps(rowRanges(vcf), bed)
  if(length(ov.idx) > 0){
    # Identify Non-Coding and Coding mutations
    nc.muts <- unique(queryHits(ov.idx))
    c.muts <- c(1:length(vcf))[-nc.muts]
    
    message(paste0(length(nc.muts), "/", length(vcf),
                   " mutations found in target region"))
  } else {
    c.muts <- 1:length(vcf)
    nc.muts <- 0
    
    message(paste0("No (0/", length(vcf),
                   ") mutations found in target region"))
  }
  if(length(nc.muts) == 0) nc.muts <- 0
  if(length(c.muts) == 0) c.muts <- 0
  
  ## Split the non-target genes into their maximum likelihood location-classifications
  loc_all <- locateVariants(vcf[c.muts,], txdb, 
                            AllVariants(promoter=PromoterVariants(
                              upstream=2000,downstream=500)))
  .getMaxLoc <- function(i){
    names(sort(table(i$LOCATION), decreasing = T))[1]
  }
  locs <- factor(sapply(split(loc_all, loc_all$QUERYID), .getMaxLoc), 
                 levels = levels(loc_all$LOCATION))
  other.vcfs <- as.list(split(vcf[c.muts, ], locs))
  
  
  return(append(list("target"=vcf[nc.muts,]),
                other.vcfs))
}