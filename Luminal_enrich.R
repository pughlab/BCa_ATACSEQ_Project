##Load libraries
library(DescTools)

####################
#### ENRICHMENT ####
####################

#### Functions #####
genRandData <- function(mus=c(0.95, 0.05, 0.05), sds=c(0.05, 0.05, 0.05)){
  labels <- c(rep("lum", 5), rep("her2", 3), rep("tnbc", 10))
  
  dat <- c(rnorm(n = 5, mean = mus[1], sd = sds[1]), 
           rnorm(n=3, mean=mus[2], sd=sds[2]),
           rnorm(n=10, mean=mus[3], sd=sds[3]))
  dat[dat>1] <- 1; dat[dat<0] <- 0; 
  dat <- setNames(dat, labels)
  return(dat)
}

getEnrichment <- function(dat, id, nperm=1000){
  dat_x <- dat
  
  ## Get indices of cancer-type and set non-cancer-type to f(-x)
  idx <- (names(dat_x) != id)
  dat_x[idx] <- -1*dat_x[idx]
  
  ## Use cancer-type and non-cancer type as controls to get range
  lum_auc <- DescTools::AUC(1:length(dat), as.integer(!idx))
  nonlum_auc <- DescTools::AUC(1:length(dat), -1*as.integer(idx))
  
  ## Get AUC for cancer-type and standardize it
  auc <- DescTools::AUC(1:length(dat_x), dat_x)
  auc_standardized = (auc-nonlum_auc)/(lum_auc-nonlum_auc)
  
  ## Do a permutation test to get significance
  rand_aucs <- sapply(1:nperm, function(i){
    dat_y <- dat
    idx_y <- sample(1:length(dat), size = sum(idx), replace = F)
    
    dat_y[idx_y] <- -1*dat_y[idx_y]
    DescTools::AUC(1:length(dat_y), dat_y)
  })
  pval <- (sum(rand_aucs > auc) / length(rand_aucs))
  
  return(data.frame("auc"=auc, "s_auc"=auc_standardized, "p"=pval))
}

#### Main #####

### read essentiality screen matrix
dat <- read.csv("",header = TRUE, sep = ",") 
###give cell lines lum, her2 or tnbc based on the subtype
colnames(dat) <- c(rep("lum",5), rep("her2",3),rep("tnbc", 10))


p_vals<-apply(dat_class, 1, function(dat){
  enrich_vals <- t(sapply(c('lum', 'her2', 'tnbc'), function(id){
    getEnrichment(dat, id, 1000)
  }))
  return(enrich_vals)
})

lapply(c("lum", "her2", "tnbc"), function(id){
  id_pval <- t(sapply(p_vals, function(i) as.numeric(i[id,])))
  id_pval <- id_pval[order(id_pval[,3]),]
  id_pval
})


