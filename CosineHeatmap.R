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




#d<- dist(raw, diag = F, method = "canberra")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("ComplexHeatmap")
#library(ComplexHeatmap)
#ht_list = Heatmap(d, name = "mat", show_column_dend = FALSE) +
#  rowAnnotation(foo = anno_barplot(colSums(raw), width = unit(4, "cm"))) +)
#draw(ht_list, padding = unit(c(2, 2, 10, 2), "mm")) # add space for titles
#decorate_annotation("foo", { 
#  grid.text("Number of peaks", y = unit(1, "npc") + unit(2, "mm"), just = "bottom") 
#})
#decorate_heatmap_body("mat", {
#  grid.text("title for the heatmap", y = unit(1, "npc") + unit(2, "mm"), just = "bottom")
#})



#library(network)
#library(sna)
#library(ggplot2)
#library(devtools)
#install_github("ggobi/ggally")
#library(GGally)
#net = network(d, mode="circle")
#ggnet2(net, label=TRUE)
#plot(net)

#qn<-normalize.quantiles(raw)
#colnames(qn)<-colnames(raw)
#rownames(qn)<-rownames(raw)
#qn.cor <- cor(qn,method="pearson")
#library(pheatmap)
#pheatmap (mat = raw,
#col= colorRampPalette(c("blue","white","red"))(30),
#cluster_rows = T,
#cluster_col = T,
#cellheight = 10,
#cellwidth = 10,
#fontsize = 12
#)


library(reshape2)
d <- apply(raw, 2, function(i){
  apply(raw, 2, function(j) sum(rowSums(data.frame(i,j)) == 2))
})
library(igraph)
library(scales)                
g <- graph.adjacency(d, weighted = TRUE, mode ='undirected')
g<- simplify(g)
V(g)$color <- ifelse(V(g)$name == 'Catalogue.EP', alpha('Red',0.4), alpha('Blue', 0.4))
test.layout <- layout_(g,with_dh(weight.edge.lengths = edge_density(g)))
dir="/Users/selghamr/OneDrive - UHN/Documents/ATAC-results-xls/Manuscript_plots/"
pdf(file=paste(dir, "MotifSitesNetwork-Plot", ".pdf", sep="")
    , width = 10, height = 8, pointsize = 12, useDingbats = FALSE)
plot(g,
     vertex.label=gsub("\\.0\\..*", "", colnames(raw)),
     vertex.size=log(colSums(raw))*1, 
     edge.width=E(g)$weight*0.001, layout = test.layout)
dev.off()               





TF2 <- substr(jacc$file2,1, 10)

library(NMF)
for(i in 2:10) {
  print(i)
  fit<-nmf(x = raw, r = i, method = "brunet",  seed=123456, nrun = 10)
  save(fit, file = paste0("Rank",i,".rData"))
}

Stats<-matrix(ncol=6, nrow=10)
colnames(Stats)<-c("Cophenetic","BasisSparseness","CoefSparseness","Residuals","Dispersion","Deviance")
for(i in 2:10){
  fname<-paste0("Rank",i,".rData")
  load(fname)
  Stats[i,1]<-cophcor(fit)
  Stats[i,2]<-sparseness(fit)[1]
  Stats[i,3]<-sparseness(fit)[2]
  Stats[i,4]<-residuals(fit)
  Stats[i,5]<-dispersion(fit)
  Stats[i,6]<-deviance(fit)
}

Stats<-Stats[-1,]
par(mfrow=c(3,2))
plot(2:10,Stats[,2], main="Basis Sparsity", pch=20, col="red", type="l", ylab="sparseness",xlab="Number of Signatures")
points(2:10,Stats[,2],pch=20, col="red")
abline(v=4, lty=2)
plot(2:10,Stats[,3], main="Signature Sparsity", pch=20, col="green", type="l", ylab="sparseness",xlab="Number of Signatures")
points(2:10,Stats[,3],pch=20, col="green")
abline(v=4, lty=2)
plot(2:10,Stats[,1], main="Cophenetic", pch=20, col="blue", type="l", ylab="Cophenetic Coefficient",xlab="Number of Signatures")
points(2:10,Stats[,1],pch=20, col="blue")
abline(v=4, lty=2)
plot(2:10,Stats[,4], main="Residual Error", pch=20, col="orange", type="l", ylab="Residual Error",xlab="Number of Signatures")
points(2:10,Stats[,4],pch=20, col="orange")
abline(v=4, lty=2)
plot(2:10,Stats[,5], main="Dispersion", pch=20, col="violet", type="l", ylab="Dispersion",xlab="Number of Signatures")
points(2:10,Stats[,5],pch=20, col="violet")
abline(v=4, lty=2)

par(mar=c(10,15,15,15))
load("Rank6.rData")
tCo<-coef(fit)
bCo<-basis(fit)

library(gplots)
Palette<-c("white", "#FFE9E8","#FFD3D1","#FFBDB9","#FFA7A2","#FF918B","#FF7B74","#FF655D","#FF4F46","#FF392E","#FF2317","#FF0D00")
heatmap.2(tCo,trace="none", col=Palette, colsep = 0:100, rowsep=0:22, sepcolor = "gray45")
pheatmap(qn.cor)


myPCA <- prcomp(raw, scale=F, center = F)
par(mfrow=c(1,1))

#Plot PCA
df <- data.frame(myPCA$rotation)
df2 <- data.frame(df[0], df$PC1, df$PC2)
library(ggfortify); library(ggplot2)
autoplot(df2)
autoplot(kmeans(df2,3), data = df2, label=TRUE)+ theme(panel.background = element_blank(), text = element_text(size=20), axis.line = element_line(colour = "black"))
a <- predict(myPCA$rotation)



df.rene <- t(dat[,-c(1:3)])
library(Rtsne)
b <- Rtsne(df.rene, initial_config = NULL, k = 8, initial_dims = 30, perplexity = 2,
          max_iter = 1100, min_cost = 0, epoch_callback = NULL, whiten = TRUE,
          epoch=100)
tsne_plot <- data.frame(x = b$Y[,1], y = b$Y[,2], col = gsub("_rep.*", "", rownames(df)))
ggplot(tsne_plot) + geom_point(aes(x=x, y=y, color=col))


ggplot(b, aes(x=b[,1], y=b[,3])) +  
  geom_point(size=3) +
  guides(colour=guide_legend(override.aes=list(size=6))) +
  ggtitle("t-SNE") +
  theme_light(base_size=20) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank()) +
  scale_colour_brewer(palette = "Set2")
plot(a)

library(lsa)


a <- cosine(raw)

raw[raw>0.1]<- 1

library(gplots)
Palette<-c("white", "#FFE9E8","#FFD3D1","#FFBDB9","#FFA7A2","#FF918B","#FF7B74","#FF655D","#FF4F46","#FF392E","#FF2317","#FF0D00")
heatmap.2(a,trace="none", colsep = 0:100, rowsep=0:22, sepcolor = "gray45")

col_breaks <- c(seq(0,0.9,length=199),seq(0.91,1,length=100))
#red_palette <- colorRampPalette(c("white", "red"))(n = 299)
my_palette <- colorRampPalette(c("#3288bd", "#fee08b",  "#d53e4f"))(n = 299)
pheatmap (mat = raw,
          #col = my_palette, breaks = col_breaks,
          col= colorRampPalette(c("blue","white","red"))(100),
          cluster_rows = T,
          cluster_col = T,
          #cellheight = 10,
          #cellwidth = 10,
          fontsize = 12
)




dat<-read.table("~/samwise/projects/KOMBAT/Targetedseq/MutRegEle/FOXA1-IGR_sig_Down_annotated-target.bed", sep="\t")
mean(dat$V3-dat$V2)
