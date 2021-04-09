###
# ATAC Saturation Calculation 

#####
# SETUP
#####
#load libraries 
library(GenomicRanges)
library( matrixStats)
library(reshape)
library(nlstools)
#####
# FUNCTIONS 
######
# make a catalogue
makeCatalogue <- function (atac.dir){
  files <- list.files(path = atac.dir, pattern="_peaks.narrowPeak")
  setwd(atac.dir)
  peaks.l<-list()
  number_peaks<-vector()
  for (i in 1:length(files)) {
    peaks.l[[i]]<-read.delim(paste(files[i],sep=""),sep="\t",header=F)
    peaks.l[[i]]<-peaks.l[[i]][,1:3]
    number_peaks<-c(number_peaks,nrow(peaks.l[[i]]))
  }
  names(peaks.l)<-files
  return (list(peaks.l, number_peaks))
}

# randomize sample and calculate saturation of unique peaks(bootsrapping)
calculateSaturation <- function (n=1, peak.ranges){
  draw.order <- sample (1:(length(peak.ranges)), length(peak.ranges), replace = F) 
  unique.peaks <- peak.ranges[[draw.order[1]]]
  num.unique.peaks <- length(unique.peaks)
  for (i in draw.order[-1]){
    unique.peaks2 <- c(unique.peaks, peak.ranges[[i]])
    unique.peaks <- reduce (unique.peaks2)
    num.unique.peaks <- c(num.unique.peaks, length(unique.peaks))
  }
  return(num.unique.peaks)
}

#####
# MAIN 
######

atac.dir <- "~/Samwise/projects/KOMBAT/atacseq_analysis/Sorted_ATACSeq_hg19/Peaks_KOM-Immune_BC/ATAC_Peaks_All/" 

n <- 1000 # number of simulations
peak.cat <- makeCatalogue(atac.dir)
peak.ranges <- lapply(peak.cat[[1]],
                      makeGRangesFromDataFrame, 
                      seqnames.field = "V1", 
                      start.field = "V2", 
                      end.field = "V3")
saturation.calc <- do.call (rbind.data.frame, 
                            lapply(1:n, calculateSaturation, peak.ranges = peak.ranges))
#saturation.calc <- t(saturation.calc)
colnames(saturation.calc) <- paste("Sample", c(1: (length(peak.cat[[2]]))))

#reshape matrix to allow model fitting
x <- melt(saturation.calc)
data.m3 <- x
data.m3[,1] <- as.integer(gsub("Sample ", "", data.m3[,1]))
colnames(data.m3) <- c("sampleNum", "uniquePeaks")
nlsfitSS <- nls(uniquePeaks ~ SSasymp(sampleNum, Asym, R0, lrc),
                data=data.m3)

new.sams = seq(1,200,1)
pred.3 <- predict(nlsfitSS,list(sampleNum=new.sams))
Asym_coef<-summary(nlsfitSS)$coefficients[1] #100% saturation
R0_coef<-summary(nlsfitSS)$coefficients[2]
lrc_coef<-summary(nlsfitSS)$coefficients[3]
sat95<-Asym_coef*0.95 #peaks at 95% saturation
sat99<-Asym_coef*0.99 #peaks at 99% saturation
SN_sat95<--(log((sat95-Asym_coef)/(R0_coef-Asym_coef))/exp(lrc_coef)) #number of samples to reach 95% saturation
SN_sat99<--(log((sat99-Asym_coef)/(R0_coef-Asym_coef))/exp(lrc_coef)) #number of samples to reach 99% saturation
curpc.median<-round((mean(split(data.m3, f=data.m3$sampleNum)[['25']]$uniquePeaks)/Asym_coef)*100)

#potting saturation curve

boxplot(saturation.calc, ylim=c(0,150000), xlim=c(0,80))
lines(pred.3, col="dodgerblue", lwd=1,lty="dashed",xlim=c(0,50), ylim=c(0,300000),)

