#' assembleGR
#' @description A standardized way to take in a bed data frame
#' and convert it to GRanges object while maintaining the same
#' seqlevel style across the analysis
#'
#' @param df Data frame of bed file
#' @param chrs UCSC style chrs 1-22,X,Y
#'
#' @return GRanges of BED file for chrs 1-22,X,Y
#' @export
#'
#' @examples
assembleGR <- function(df, chrs=paste0("chr", c(1:22, "X", "Y"))){
  gr <- makeGRangesFromDataFrame(df, keep.extra.columns = TRUE)
  seqlevelsStyle(gr) <- 'UCSC'
  gr[seqnames(gr) %in% chrs,]
}