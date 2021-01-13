#' applySlidingWindow
#' @description Applies a sliding window approach of window size "w"
#' to count the number of SNVs within that window. For now, it grabs the
#' window with the maximum number of SNVs and returns that region, gene,
#' and which SNVs populated that region.
#' 
#' @param vcf.bed A list output from separateVcfByBed() function, contains
#' a list of VariantAnnotation objects in the first element ([[1]]), and 
#' a GRanges object in the second element with gene name ([[2]])
#' @param w window size (Default: 3)
#' @param w.scale window size scale (Default: 1000)
#'
#' @return List containing a reduced GRanges object of windows with all SNPs
#' in that promoter range, the maximum number of SNVs within a window size
#' for that promoter, as well as SimpleList object of VariantAnnotations of those 
#' SNVs
#' @export
#'
#' @examples
applySlidingWindow <- function(vcf.bed, w=3, w.scale=1000){
  # Return the simplified w-bin to actual kb size
  w <- (w * w.scale)/2
  
  prom.idx <- names(vcf.bed[['vcf']]) 
  # p <- '2686' should be 26 total, 20 max
  summ.w <- lapply(prom.idx, function(p){
    ## Extract the promoter bed and the SNVs in that region
    prom <- vcf.bed[['bed']][p,]
    snv <- vcf.bed[['vcf']][[p]]
    anno.samples <- any(grepl("^ID$", colnames(info(snv))))
    
    ## Count all SNVs in that promoter
    ov.cnt <- suppressWarnings(countOverlaps(prom, snv))
    
    ## Annotate the entire promoter count
    prom$snv <- paste(names(rowRanges(snv)), collapse=",")
    prom$cnt <- ov.cnt
    if(anno.samples) prom$samples <- paste(info(snv)$ID, collapse=",")
    
    max.run.ov <- slidingWindow(gr.range=prom, w=w, obj=snv)
    max.run <- max.run.ov[['run']]
    max.run$gene <- prom$gene
    max.run$snv <- paste(names(rowRanges(snv))[max.run.ov[['ov.idx']][[1]]], collapse=",")
    if(anno.samples) max.run$samples <- paste(info(snv)$ID[max.run.ov$ov.idx[[1]]], collapse=",")
    
    return(list("prom"=prom,
                "w"=max.run,
                "snvs"=snv[max.run.ov[['ov.idx']][[1]],]))
  })
  
  reduced.raw.gr <- unlist(as(lapply(summ.w, function(x) x[['prom']]), "GRangesList"))
  reduced.gr <- unlist(as(lapply(summ.w, function(x) x[['w']]), "GRangesList"))
  reduced.vcf <- as(lapply(summ.w, function(x) x[['snvs']]), "SimpleList")
  return(list("t.summ"=reduced.raw.gr,
              "max.summ"=reduced.gr,
              "snvs"=reduced.vcf))
}


#' Applies sliding window overlap
#' @description Uses two GRanges objects (gr.range and obj).  Uses a sliding window of 
#' size 'w' to slide over gr.range, looking for the maximum number of overlaps of the
#' obj GRAnges object
#'
#' @param gr.range GRanges object: A single interval to look across
#' @param obj GRanges object: 1 or more intervals that fall within the gr.range
#' @param w Integer specifying the size of the search window in bp
#'
#' @return list: GRanges object with the maximum number of overlaps from obj found in gr.range
#' as well as the indices of which objs were found within gr.range
#' @export
#'
#' @examples
slidingWindow <- function(gr.range, obj, w){
  ## Create an artificial sliding window based on the bin-size
  range.idx <- c(start(gr.range):(end(gr.range) - w))
  gr <- GRanges(seqnames=as.character(seqnames(gr.range)), 
                IRanges(min(range.idx):max(range.idx),
                        width=w))
  
  ## For the SlidingWindow, count the overlaps with the VCF
  ov.cnt <- countOverlaps(gr, obj)
  max.ov <- max(ov.cnt) # Configured only for MAX at the moment
  max.idx <- which(ov.cnt %in% max.ov)
  
  ## If there are multiple ranges that have the exact value,
  ## split up the ranges and treat them individually
  # Finds continuous ranges of intervals
  cont.interval <- rep(1, length(max.idx))
  changepoints <- which(diff(max.idx) != 1)
  if(length(changepoints) > 0){
    for(cp in changepoints){
      tail.end <- (cp+1):length(cont.interval)
      cont.interval[tail.end] <- cont.interval[tail.end]+1
    }
  }
  # Splits up continuous ranges and reduces
  split.intervals <- split(max.idx, cont.interval)
  max.runs <- as(sapply(split.intervals, function(idx){
    reduce(gr[idx,])
  }), 'GRangesList')
  max.run <- unlist(max.runs) # Isolate for just w-range harboring SNVs
  
  ## Identify which SNVs are in that w-range
  ov.idx <- findOverlaps(gr,obj)
  snv.idx <- lapply(split.intervals, function(max.idx){
    ov.subset <- ov.idx[queryHits(ov.idx) %in% max.idx,] # reduce SNV range
    max.obj.idx <- unique(subjectHits(ov.subset)) # get unique SNVs
  })
  
  ## Annotate the reduced w-range
  elementMetadata(max.run)$cnt <- sapply(snv.idx, length)
  
  list("run"=max.run, "ov.idx"=snv.idx)
}