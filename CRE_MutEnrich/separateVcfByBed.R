#' separateVcfByBed
#' @description Simple function that separates the SNVs
#' into their non-coding promoters. Intention is to
#' limit the search space for the sliding window
#' to only promoters that harbor SNVs
#' 
#' @param vcf VariantAnnotation object
#' @param bed GenomicRanges object of bed file
#'
#' @return Returns a list containing a VariantAnnotation 
#' object of each interval split out. And a Bed GRanges
#' object of only relevant promoters with their names
#' matching the VariantAnnotation list
#' @export
#'
#' @examples
separateVcfByBed <- function(vcf, bed){
  require(VariantAnnotation)
  require(GenomicRanges)
  
  ov.idx <- findOverlaps(bed, rowRanges(vcf))
  ## Separate VCF based on unique promoters
  vcf.bed <- split(vcf[subjectHits(ov.idx),],
                   f=queryHits(ov.idx))
  
  ## Identify all bed intervals with SNVs in them
  bed.vcf <- bed[unique(queryHits(ov.idx)),]
  names(bed.vcf) <- names(vcf.bed)
  
  return(list("vcf"=vcf.bed,
              "bed"=bed.vcf))
}