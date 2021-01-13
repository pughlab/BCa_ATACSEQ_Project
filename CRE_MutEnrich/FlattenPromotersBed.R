library(GenomicRanges)
wdir <- ""
promoters <- read.table(file.path(wdir, "promoters.bed"),
                        sep="\t", header=FALSE, stringsAsFactors = FALSE,
                        check.names=FALSE, quote="", comment.char="")
colnames(promoters) <- c("chr", "start", "end", "gene")
chrs <- paste0("chr", c(1:22, "X", "Y"))

## Reduce without gene information
grx <- sort(makeGRangesFromDataFrame(promoters, keep.extra.columns = TRUE))
grx.gene <- split(grx, grx$gene)
red.grx.gene <- reduce(grx.gene)

## Identify promoters in multiple chromosomes
grl.count<- sapply(grx.gene, function(x)
  length(unique(sort(as.character(seqnames(x))))))
grl.count.sort<- sort(grl.count, decreasing = TRUE)

## Concatenate, annotate gene, and convert to Data.Frame
red.grx.g <- unlist(red.grx.gene)
elementMetadata(red.grx.g)$Gene <- names(red.grx.g)
names(red.grx.g) <- NULL
red.df.g <- as.data.frame(sort(red.grx.g[seqnames(red.grx.g) %in% chrs,]))

write.table(red.df.g[,c("seqnames", "start", "end", "Gene")], 
            file=file.path(git.dir, "Gencodev19_Promoters_REDUCE.bed"),
            sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)