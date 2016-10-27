### run ml-2dqsar
### szu
### 2016-10-20

### accept args

args <- commandArgs(TRUE)
### set parameter.
jobdir <- "/data/home/szu/Lab/ml-2dqsar/code/result/"
datadir <- "/data/home/szu/Lab/ml-2dqsar/"
ml2dqsar_dir <- "/data/home/szu/Lab/ml-2dqsar/code/"

ml2dqsar_filenm <- "FuncPreComPiHierBayes.R"
wiperfunc_filenm <- "wiperdata.R"

## gpcr_names <- c("CPIs_GPCR_monoamine_receptor",
##                 "CPIs_GPCR_nucleotide-like_receptor",
##                 "CPIs_GPCR  lipid-like ligand receptor",
##                 "CPIs_Peptide GPCR")
gpcr_names <- c("lipidlike",
                "monoamine",
                "nucleotidelike",
                "peptide",
                "others")
kinase_names <- c("kinase_agc", "kinase_camk","kinase_ck1","kinase_cmgc",
                  "kinase_ste","kinase_tk","kinase_tkl")

fealist <- c("macc","phychem","fingerprint","morgan")
isbinaryf <- c(TRUE, FALSE, TRUE, TRUE)

nfold <- 5
lower_limit <- 50
upper_limit <- 1000
sigma2dot <- 100
iternum <- 1000

isGPCR <- TRUE
subd_index= c(1) # subdirs choose index
fc_index= c(1) # feature choose index

cutoff = 0.05 # to filter the features.
###------keep unchanged below----------###
### packages
library(methods)
library(glmnet)
library(data.table)
###library(ggplot2)
library(optimx)
###library(doParallel)
###library(foreach)
###library(caret)
library(matrixStats)
library(Metrics)

### set dir, and load func.
setwd(jobdir)
source(paste(ml2dqsar_dir,ml2dqsar_filenm,sep=""))
source(paste(ml2dqsar_dir,wiperfunc_filenm,sep=""))

### load features.
features <- fealist[fc_index]
isbinaryfc <- isbinaryf[fc_index]
if(isGPCR){
  subdirs <- gpcr_names[subd_index]
} else{
  subdirs <- kinase_names[subd_index]
}

### run
for(f in 1:length(features)){# f for feature
  for(k in 1:length(subdirs)){ # k for subdirs.
    ## load data.
    ldat = loadfiles(datadir,gpcr_nm, kinase_nm,isGPCR,f,k)
    save(cpi, file=paste('cpi_',subdirs[k],'_',features[f],'.data',sep=""))
    save(testResult, file=paste('testResult_',subdirs[k],'_',features[f],'.data',sep=""))
  }
}
