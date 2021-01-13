#' estimateBgMutRate
#' @description 
#' 
#' @param prom 
#' @param prom.snv
#' @param hotspots 
#'
#' @return
#' @export
#'
#' @examples
estimateBgMutRate <- function(prom, prom.snv, hotspots, pad=TRUE){
  require(MASS)
  
  # Identify all covered/targetted regions in promoters/peaks file
  pr.cov <- split(prom, prom$peak)[['covered']]
  
  # Set all covered/targetted regions with no mutations to a mut.cnt of 0
  ov.idx <- findOverlaps(pr.cov, hotspots)
  all.pr.cov <- c(1:length(pr.cov))
  pr.nomut <- pr.cov[all.pr.cov[-unique(queryHits(ov.idx))],]
  if(length(pr.nomut > 0)) pr.nomut$cnt <- 0
  
  # Use the MASS package to fit a poisson distribution to estimate lambda
  mut.cnt <- list("pad"=pr.nomut$cnt,
                  "null"=(prom.snv$cnt / width(prom.snv))) ## Total SNVs within each range
  # mut.cnt <- c(pr.nomut$cnt, prom.snv$cnt) 
  lambda.mle.pad <- if(!is.null(mut.cnt$pad)){
    suppressWarnings(fitdistr(c(mut.cnt$pad, mut.cnt$null), densfun="poisson"))
  } else {
    NULL
  }
  lambda.mle.null <- if(length(mut.cnt$null)>0){
    suppressWarnings(fitdistr(mut.cnt$null, densfun="poisson"))
  } else {
    NULL
  }
  
  # Summarize the counts for downstream p-value estimations
  all.max.gr <- c(pr.nomut, hotspots)
  if(any(is.na(all.max.gr$peak))) all.max.gr[which(is.na(all.max.gr$peak)),]$peak <- 'covered'
  
  # Return lambda for the estimate of background mutation rate
  return(list("lambda.pad"=lambda.mle.pad,
              "lambda.null"=lambda.mle.null,
              "max.gr"=all.max.gr))
}