#' mutRatePvalue
#' @description Calculates the binomial and poisson probability for
#' the number of max mutations in each promoter/w-bin size. Also 
#' applies a multiple hypothesis adjustment test to the pvalue
#' 
#' 
#' @param hotspot A GRanges object containing the promoter and the
#' number of SNVs in each promoter under the column $cnt
#' @param w Integer size of the bin (Default: 3)
#' @param prob Background mutation probability rate (Default: 0.000075)
#' @param kb Integer size of bin scale (Default: 1000)
#' @param stat The statistical test used to calculate a p-value, poisson or binomial (Default: poisson)
#' @param method Multiple hypothesis test correction (Default: bonferroni)
#'
#' @return Returns the hotspot GRanges object but with added columns
#' in the elementMetadata for p and q binomial and poisson values
#' @export
#'
#' @examples
mutRatePvalue  <- function(hotspot, w, prob=0.000075, 
                           kb=1000, stat = 'poisson', method='fdr'){
  # mutation rate: 2.5 x 10(-8) per nucleotide
  # mutation rate: 0.000025 per 1kb
  # mutation rate: 0.000075 per 3kb (default maxW)
  
  # Based on a binomial with average size for promoter
  if(length(hotspot) > 0){
    p.l <- lapply(prob, function(p){
      if(is.null(p)) p=0.000075
      p.mat <- sapply(stat, function(s){
        switch(s,
               poisson={
                 ppois(hotspot$cnt, lambda=p * width(hotspot), lower=FALSE)
               },
               binomial={
                 pbinom(q = hotspot$cnt, size = width(hotspot), # w*kb, 
                        prob=p, lower.tail=FALSE) # prob / (w*kb)
               })
      })
      colnames(p.mat) <- paste0("p.", stat)
      as.data.frame(p.mat)
    })
    p.df <- do.call(cbind,p.l)
    colnames(p.df) <- as.character(sapply(names(prob), function(i) paste0(i, ".", paste0("p.", c(stat)))))
    
    ## Multiple hypothesis correction
    q.df <- apply(p.df, 2, p.adjust, method=method)
    colnames(q.df) <- gsub("p\\.", "q.", colnames(q.df))
    
    pq.df <- round(cbind(p.df, q.df),4)
    mcols(hotspot) <- cbind(mcols(hotspot), pq.df)
  } else {
    col.ids <- as.character(sapply(names(prob), function(i) paste0(i, ".", paste0("p.", c(stat)))))
    col.ids <- c(col.ids, gsub("\\.p\\.", ".q.", col.ids))
    col.df <- as.data.frame(matrix(ncol=length(col.ids), nrow=0))
    colnames(col.df) <- col.ids
    
    mcols(hotspot) <- cbind(mcols(hotspot), col.df)
  }
  
  hotspot
}