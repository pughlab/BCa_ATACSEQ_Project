# Rene Quevedo & Samah El Ghamrasni 
# Department of medical biophysics, University of Toronto
# Princess Margaret Cancer Centre - UHN

### Load libraries ###
library(optparse)
library(GenomicRanges)
library(VariantAnnotation)
library(MASS)
	
#### Required parameters ####

option_list <- list(
  make_option(c("-i", "--id"), type="character", default=NULL,
              help="VCF file to calculate mutation rate for [default= %default]", metavar="character"),
  make_option(c("-l", "--label"), action="store_true", default=TRUE,
              help="Keep the input label for each interval, else replaces it with chr:start-end [default]"),
  make_option(c("-x", "--suffix"), type="character", default="FILEID-X.vcf",
              help="Suffix or the sample ID to distinguish between Tumor and Normal sample,
              for instance, 'FILEID-X.vcf' will scan for FILEID-T.vcf and FILEID-N.vcf (case sensitive)
              in their respective directories (somatic or germline)", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL,
              help="Directory to save new files to [default= %default]", metavar="character"),
  make_option(c("-w", "--w"), type="integer", default=NULL,
              help="Bin size that is a multipler of --wscale, if NULL, will determine bin size", metavar="integer"),
  make_option(c("-z", "--wscale"), type="integer", default=1000,
              help="Multiplier of bin-size", metavar="integer"),
  make_option(c("-m", "--somaticdir"), type="character", default="",
              help="Path to directory containing somatic VCF file [default= %default]", metavar="character"),
  make_option(c("-n", "--germlinedir"), type="character", default=NULL,
              help="Path to directory containing germline VCF file [default= %default]", metavar="character"),
  make_option(c("-g", "--gitdir"), type="character", default="~/git/BCa_ATACSEQ_Project/CRE_MutEnrich",
              help="Path to Git directory [default= %default]", metavar="character"),
  make_option(c("-r", "--refdir"), type="character", default="",
              help="Path to directory containing target/catalogue files [default= %default]", metavar="character"),
  make_option(c("-t", "--target"), type="character", default="Gencodev19_Promoters_REDUCE.bed",
              help="Targets file name [default= %default]", metavar="character"),
  make_option(c("-c", "--catalogue"), type="character", default="",
              help="Catalogue file name [default= %default]", metavar="character"),
  make_option(c("-f", "--flank"), type="integer", default=0,
              help="Flanking regions (bp) around the catalogue if you don't have a targets bed
              file and want to capture mutations in the flanking regions of your catalogue", metavar="character"),
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print extra output [default]")
  )

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)


### Setup Parameters ###
opt$id=''
opt$refdir=''
opt$target=''  ## regions of interest
opt$catalogue='' ## accessible chromatin catalogue
opt$suffix='FILEID-X.vcf'
opt$label=FALSE
opt$outdir=''
opt$w=6
opt$somaticdir=''
opt$germlinedir=''
opt$verbose = TRUE
opt$wscale = 100 
opt$flank = 200 
opt$suffix <- gsub("^FILEID", "", opt$suffix)

vcf.file <- c("somatic"=paste0(opt$id, gsub("X", "T", opt$suffix)),
              "germline"=paste0(opt$id, gsub("X", "N", opt$suffix))) 
out.dir <- opt$outdir # "./"

germline.dir <- opt$germlinedir # NULL
somatic.dir <- opt$somaticdir 
ref.dir <- opt$refdir 
git.dir <-  opt$gitdir 
target.file <- opt$target 
catalogue.file <- opt$catalogue 
flank <- opt$flank

##############
#### MAIN ####
## Temporary source functions
for(i in list.files(file.path(git.dir, "R"))){
  source(file.path(git.dir, "R", i))
}

#### Load in data ####
# Catalogue
catalogue <- read.table(file.path(ref.dir, catalogue.file),
                        header=FALSE, check.names=FALSE, quote="",
                        stringsAsFactors = FALSE)
colnames(catalogue) <- c("chr", "start", "end", "gene")
pk.gr <- suppressWarnings(assembleGR(catalogue))

# Targets
targets <- read.table(file.path(ref.dir, target.file),
                      header=FALSE, check.names=FALSE, quote="",
                      stringsAsFactors = FALSE)
colnames(targets) <- c("chr", "start", "end", "gene")
if(flank > 0){
  targets$start <- targets$start - flank
  targets$end <- targets$end + flank
}

pr.gr <- suppressWarnings(assembleGR(targets))
if(!opt$label){
  message("Renaming target interval annotations to chr:start-end...")
  replace.ids <- c('^DistalRE$')
  rpl.idx <- as.integer(unlist(sapply(replace.ids, grep, x=pr.gr$gene)))
  locs <- paste0(as.character(seqnames(pr.gr)), ":",
                 start(pr.gr), "-", end(pr.gr))
  pr.gr[rpl.idx,]$gene <- locs[rpl.idx]
}

# VCF
vcf.path <- ''
vcf <- sort(suppressWarnings(filterVCF(readVcf(vcf.path))))


#### Pre-processing ####
### Flags all overlaps between the two search spaces
pr.gr <- overlapTargetsAndCatalogue(pk.gr, pr.gr) # targets covered by catalogue
pk.gr <- overlapTargetsAndCatalogue(pr.gr, pk.gr) # catalogues covered by targets
all.ref.gr <- list('target'=pr.gr, 'catalogue'=pk.gr)
null <- if(opt$verbose) compareSearchSpace(pk.gr, pr.gr, 'target')

## Split non-coding from coding mutations
vcf <- removeCodingMuts(vcf, pr.gr)

## Set the bin-size for the sliding window
#if(is.null(opt$w)) opt$w <- determineW(vcf.spl[['non-coding']], q=5)
w <- opt$w
w.scale <- opt$wscale

## Break down variants into based on their search ranges
# Somatic: Catalogue and Target, then each SNV location
s.gr <- separateVcfByBed(vcf$target, bed=all.ref.gr$target)

#### Hotspot Calling ####
## Apply sliding window
# t.summ = 'all SNPs in the given range range'
# max.summ = 'maximum number of SNVs within a window for given range'
# snvs = 'VariantAnnotations of those SNVs'
s.hotspot <- applySlidingWindow(s.gr, w=w, w.scale=w.scale)

## Estimates the background mutation rate for the sample using a poisson dist.
space.gr <- all.ref.gr$target
bg <- estimateBgMutRate(prom=space.gr,
                        prom.snv=s.hotspot$t.summ,
                        hotspots=s.hotspot$max.summ)

## Calculates the p-value for each number of mutations
nc.stat <- mutRatePvalue(hotspot = bg$max.gr,
                         prob = list('pad'=bg$lambda.pad$estimate,
                                     'null'=bg$lambda.null$estimate),
                         w=w, kb=w.scale, method='fdr', stat='poisson')
						 
nc.stat.df<- as.data.frame(nc.stat, row.names=NULL)						 

write.table(nc.stat.df, "nc.stat.df.txt", append = FALSE, sep = "\t", dec = ".",
              row.names = TRUE, col.names = TRUE)

