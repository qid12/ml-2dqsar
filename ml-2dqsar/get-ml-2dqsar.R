### run ml-2dqsar
### szu
### 2016-10-20

### accept args

## args <- commandArgs(TRUE)
### set parameter.
jobdir <- "/data/home/szu/Lab/ml-2dqsar/getmcode/result/"
datadir <- "/data/home/szu/Lab/ml-2dqsar/"
ml2dqsar_dir <- "/data/home/szu/Lab/ml-2dqsar/getmcode/"

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
fi= 1 # feature choose index

cutoff = 0.05 # to filter the features.

isaddstd <- TRUE
isusemin <- FALSE
useminv <- 0.05
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

isbinaryfc <- isbinaryf[fi]
if(isGPCR){
  subdirs <- gpcr_names[subd_index]
} else{
  subdirs <- kinase_names[subd_index]
}

### run
for(k in 1:length(subdirs)){ # k for subdirs.
  ## load data.
  ldat = loadfiles(datadir,gpcr_nm, kinase_nm,isGPCR,fi,k)
  precpi <- getpredata(ldat)
  cpi <- getcpid(precpi, fi, cutoff,lower_limit, upper_limit)

  sresult <- get_singlpar_sigma(cpi,nfold,isaddstd,isusemin,useminv)

  mresult <- get_ml2dqsarpar(cpi,sresult$sigma,sigma2dot,iternum)

  save(sresult, file=paste('s_',subdirs[k],'_',fealist[fi],'.data',sep=""))
  save(mresult, file=paste('m_',subdirs[k],'_',fealist[fi],'.data',sep=""))
}
