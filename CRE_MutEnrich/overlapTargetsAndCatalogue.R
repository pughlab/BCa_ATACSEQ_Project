#' overlapTargetsAndCatalogue
#' @description Overlaps the promoters and peaks for the sample. Flags
#' the promoters that have some ATAC-seq peak coverage on them and those
#' promoters that do not
#' 
#' @param catalogue GRanges object for peaks
#' @param targets GRanges object for targetsoters
#'
#' @return GRanges object with a cov/non-covered column to split on 
#' @export
#'
#' @examples
overlapTargetsAndCatalogue <- function(catalogue, targets){
  ov.idx <- findOverlaps(targets, catalogue)
  
  # Indexes all the targets that have catalogue coverage
  all.targets <- c(1:length(targets))
  targets.covered.by.catalogue <- unique(queryHits(ov.idx))
  
  # Sets a "coverage" flag to split on
  targets$peak <- 'non-covered'
  targets[targets.covered.by.catalogue,]$peak <- 'covered'
  
  return(targets)
}