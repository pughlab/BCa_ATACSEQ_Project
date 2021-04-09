#### Load libraries
library(pheatmap)
library(lsa)

#load matrix
dat<-read.csv("") #path to matrix
dat<-dat[,-1]
raw<-matrix(ncol=26, nrow=nrow(dat))
for(i in 1:ncol(raw)) {
  raw[,i]<-as.numeric(as.character(dat[,i+3]))
}
colnames(raw)<-colnames(dat)[4:ncol(dat)]
rownames(raw)<-paste0(dat[,1],":",dat[,2],"-",dat[,3])
samplenames <- substr(colnames(raw),1, 12)
colnames(raw)<- samplenames

### calculate cosine similarity
Cos.raw <- cosine(raw)

### Heatmap Poting
col_breaks <- c(seq(0,0.6,length=50),seq(0.61,1,length=50))
my_palette <- colorRampPalette(c("blue","white","red"))(n = 100)
pheatmap (mat = Cos.raw,
          col = my_palette, breaks = col_breaks,
          cluster_rows = T,
          cluster_col = T,
          fontsize = 15)



