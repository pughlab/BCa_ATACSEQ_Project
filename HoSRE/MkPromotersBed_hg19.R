

#### Create a hg19 promoters reference bedfile
mkPromoterBed <- function(out.file="~/Desktop/promoters.bed", ret.type="data.frame"){
  ## Make Promoter BED file
  require(org.Hs.eg.db)
  require(TxDb.Hsapiens.UCSC.hg19.knownGene)
  
  # Get all transcripts
  hg19 <- TxDb.Hsapiens.UCSC.hg19.knownGene
  transcriptCoordsByGene.GRangesList <- transcriptsBy (hg19, by = "gene") 
  
  # Map Entrez IDs to Hugo Symbols
  symbolid <- mapIds(org.Hs.eg.db,
                     keys=names(transcriptCoordsByGene.GRangesList),
                     column="SYMBOL",
                     keytype="ENTREZID",
                     multiVals="first")
  names(transcriptCoordsByGene.GRangesList) <- symbolid
  
  # Extract all promoters, add range and account for strand
  prom <- promoters(transcriptCoordsByGene.GRangesList, 
                    upstream=2500, downstream=500)
  tss <- promoters(transcriptCoordsByGene.GRangesList, 
                   upstream=0, downstream=1)
  
  
  if(ret.type=='data.frame'){
    # Rearrange to make a bed file
    prom.df <- data.frame(prom)[,c(3, 4, 5, 2, 7, 8)]
    write.table(prom.df, file=out.file, 
                sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
    
    prom.df
  } else if(ret.type=='GRanges'){
    prom.gr <- suppressWarnings(unlist(prom))
    tss.gr <- unlist(tss)
    if(all(prom.gr$tx_name == tss.gr$tx_name)){
      st <- rep(strand(prom.gr)@values, strand(prom.gr)@lengths)
      mcols(prom.gr, use.names = F)$prom.start <- start(resize(prom.gr, 1))
      mcols(prom.gr, use.names = F)$prom.end <- end(prom.gr)
      mcols(prom.gr, use.names = F)$tss <- start(tss.gr)
      mcols(prom.gr, use.names = F)[which(st == '-'),]$prom.end <- start(prom.gr)[which(st == '-')]
      mcols(prom.gr, use.names = F)$strand <- st
      prom.gr
    }
  } else {
    prom.df
  }
}
