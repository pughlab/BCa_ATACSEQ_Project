# Samah El Ghamrasni & Rene Quevedo
# Department of medical biophysics, University of Toronto
# Princess Margaret Cancer Centre - UHN, July 27, 2020

library(scales)
library(graphics)
library(gplots)
library(RColorBrewer)
library(ggplot2)
library(data.table)
library(plotrix)
require(stats)
library(reshape)
library(dplyr)
library(pheatmap)


#### modMEMOS ####
## Directory structure:
# csvs
#  |- 50flank
#  |   |- TF1-50.csv
#  |   |- TF2-50.csv
#  |   |- ....
#  |- 100flank
#  |   |- TF1-100.csv
#  |   |- TF2-100.csv
#  |   |- ....

getMutStats <- function(tf, normality_test=FALSE, stat='cohensd'){
  # Shapiro-Wilkson test for normality
  norm_bool <- apply(tf, 2, function(i){
    norm_p <- tryCatch({
      if(normality_test) {
        if(shapiro.test(i)$p.value < 0.01) FALSE else TRUE
      } else {
        TRUE
      }
    }, error=function(e){FALSE})
  })
  
  ## get p-value from column 1 [input] to column 2 [shuffle]
  # If normally distributed, do a t-test to get a pvalue
  if(all(norm_bool)){
    switch(stat,
           cnt=median(tf[,1], na.rm=TRUE),
           cohensd=effsize::cohen.d(tf[,1], tf[,2])$estimate,
           ttest=t.test(tf[,1], tf[,2], alternative = "greater")$p.value,
           wilcox=wilcox.test(tf[,1], tf[,2], alternative = "greater")$p.value)
  } else {
    NA
  }
}

path_to_csv <- ""
flanks_v <- list.files(path_to_csv, pattern="Flank")
flanks_v <- flanks_v[order(as.integer(gsub("Flank", "", flanks_v)))]

#Variables
stat <- 'cohensd' # 'cohensd', 'ttest', 'wilcox', 'cnt'
normality_test <- FALSE

# Loop through the flank folders
cntMat <- sapply(flanks_v, function(flank_i, ...){
  print(paste0(flank_i, ": ", match(flank_i, flanks_v), "/", length(flanks_v)))
  csvs_v <- list.files(file.path(path_to_csv, flank_i), pattern="csv$")
  
  # Loop through all csvs within flank X folder
  sapply(csvs_v, function(csv_i){
    tf <- read.csv(file.path(path_to_csv, flank_i, csv_i), header=FALSE, sep=",")
    getMutStats(tf, ...)
  })
}, stat='cnt', normality_test=FALSE)

pvalueMat <- sapply(flanks_v, function(flank_i, ...){
  print(paste0(flank_i, ": ", match(flank_i, flanks_v), "/", length(flanks_v)))
  csvs_v <- list.files(file.path(path_to_csv, flank_i), pattern="csv$")
  # Loop through all csvs within flank X folder
  sapply(csvs_v, function(csv_i){
    tf <- read.csv(file.path(path_to_csv, flank_i, csv_i), header=FALSE, sep=",")
    # Checks normality of TF counts in sample and shuffle
    # Then does a statistic (cohensd, ttest, wilcox, cnt) between the two
    getMutStats(tf, ...)
  })
}, stat=stat, normality_test=normality_test)

if(stat %in% c('ttest', 'wilcox', 'cohensd')){
  pvalueMat[is.na(pvalueMat)] <- -100
  if(stat %in% c('ttest', 'wilcox')){
    qvalueMat <- matrix(p.adjust(pvalueMat, method='bonferroni'), ncol=ncol(pvalueMat))
  } else {
    qvalueMat <- pvalueMat
  }
} else {
  qvalueMat <- pvalueMat
}

### Flank and TF IDs
csv_names <-  list.files(file.path(path_to_csv, flanks_v[1]), pattern="csv$")
colnames(qvalueMat) <- colnames(pvalueMat) <- gsub("Flank", "", flanks_v) ##Flanks ID
rownames(qvalueMat) <- rownames(pvalueMat) <- gsub("_peaks-bg.*\\.csv$", "", csv_names) %>% gsub(".MCF7", "", .) ##TF ID

## Create annotation file 
sigs<- median(qvalueMat_sub[,6], probs = 0.8)
mcol <- data.frame("pval50flank" =cut(qvalueMat_sub[,6],
                                       breaks = c(-Inf, sigs, Inf),
                                       labels = c("No","Yes")))
rownames(mcol)<- rownames(qvalueMat_sub)
Yes=grep("Yes", mcol$pval100flank)
No=grep("No",mcol$pval100flank)
mat_colors <- list(pval100flank=c("grey", "red"))
names(mat_colors$pval100flank) <- c("No","Yes")


#### Plotting heatmap
col_breaks <- c(-10, c(seq(-5,0,length=100),seq(0.01,10,length=100)))
my_palette <- c("#7570b3", colorRampPalette(c("#7570b3","white","#e7298a"))(n = 200))
pheatmap(qvalueMat_sub, cluster_cols = FALSE, cluster_rows = FALSE,
         breaks = col_breaks, color = my_palette, fontsize = 7,
         cellwidth= 10,
         cellheight= 10,
         annotation_row = mcol,
         annotation_colors = mat_colors)



#### Plotting line graph
qval<-melt(qvalueMat)
colnames(qval)<- c("TFs", "Flank", "sd")
ggplot(qval, aes(x=as.numeric(as.character(Flank)), y=sd, group=TFs, colour=TFs)) + 
  geom_smooth(se = TRUE)+
  geom_hline(yintercept=c(sigs, -sigs), col="red", linetype = 2)+
  theme_bw()+xlim(0,1000)+ scale_x_continuous(breaks=seq(0,1000,100))+
  scale_color_manual(values=c("#053061","#2166ac","#4393c3","#92c5de","#b2abd2","#c51b7d","#8073ac","#542788","#2d004b"))


