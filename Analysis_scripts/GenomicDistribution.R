
###Load libraries####
library(ggplot2)
library(reshape)

#####ATAC Peaks distribution#####
File_dir <- "" #CEAS directory
samples<-list.files(path=File_dir,recursive=TRUE,pattern=".R")
samples<-gsub(".R","",samples)
ceas_distrib<-data.frame(Feature = c("TSS U/S 1kb","TSS U/S 1-2kb","TSS U/S 2-3kb","TSS D/S 1kb","TSS D/S 1-2kb","TSS D/S 2-3kb" ,"5\'UTR","3\'UTR","Coding exon","Intron","Distal intergenic"), 
                         Genome = c(0.8,0.7,0.7,0.8,0.7,0.6,0.2,0.8,1.1,29.4,64.2))

for (i in 1:length(samples)) {
  tmp<-scan(paste(File_dir,samples[i],".R",sep=""),what="character",skip=102,nlines=1,quiet=T,sep="%")
  tmp2<-as.data.frame(as.numeric(unlist(regmatches(tmp, gregexpr('?[0-9.]+', tmp)))[1:11]))
  ceas_distrib<-data.frame(c(ceas_distrib,tmp2))
  rm(tmp)
  rm(tmp2) 
}

###### Values processing####
ceas_distrib.t <- t(ceas_distrib)
colnames(ceas_distrib.t) <- ceas_distrib.t[1,]

prom <- data.frame(as.numeric(ceas_distrib.t[,2]), as.numeric(ceas_distrib.t[,3]), as.numeric(ceas_distrib.t[,4]), as.numeric(ceas_distrib.t[,5]),
                   as.numeric(ceas_distrib.t[,1]), as.numeric(ceas_distrib.t[,6]))
UTRs <- data.frame(as.numeric(ceas_distrib.t[,7]), as.numeric(ceas_distrib.t[,8]))

Distribution <-data.frame(Samples=ceas_distrib.t[,0], Coding.exon= as.numeric(ceas_distrib.t[,9]), UTRs=rowSums(UTRs), 
              promoter= rowSums(prom), Intron= as.numeric(ceas_distrib.t[,11]), Distal.intergenic=as.numeric(ceas_distrib.t[,10]))
rownames(Distribution)<- c("genome", samples)

####Ploting#####

M.Distribution<- melt(Distribution)
p <- ggplot(M.Distribution, aes(x=variable, y=value)) + 
  geom_boxplot(fill=c("#8dd3c7", "#ffffb3", "#bebada","#a6cee3", "#b2df8a"))
p + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5, binwidth =1.5, size = 5)+ scale_color_grey() + theme_classic()


