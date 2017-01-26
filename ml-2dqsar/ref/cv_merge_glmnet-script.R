### run for cv of all data for glmnet
### szu
### 2016-11-07
### accept args

args <- commandArgs(TRUE)
### set parameter.
jobdir <- "/data/home/szu/Lab/ml-2dqsar/cv_alld_glmnet/result/"
datadir <- "/data/home/szu/Lab/ml-2dqsar/"
ml2dqsar_dir <- "/data/home/szu/Lab/ml-2dqsar/cv_alld_glmnet/"

ml2dqsar_filenm <- "FuncPreComPiHierBayes.R"
wiperfunc_filenm <- "wiperdata.R"
mergdglm_filenm <- "cv_merge-glmnet.R"
gpcr_names <- c("lipidlike",
                "monoamine",
                "nucleotidelike",
                "peptide",
                "others")
kinase_names <- c("kinase_agc",
                  "kinase_camk",
                  "kinase_ck1",
                  "kinase_cmgc",
                  "kinase_ste",
                  "kinase_tk",
                  "kinase_all_tkl",
                  "kinase_tkl",
                  "kinase_stt_tkl")

fealist <- c("macc","phychem","fingerprint","morgan")
isbinaryf <- c(TRUE, FALSE, TRUE, TRUE)

nfold <- 5
lower_limit <- 30
upper_limit <- 1000
sigma2dot <- 100
iternum <- 1000

if(as.numeric(args[1]) > 0){
  isGPCR <- TRUE
} else{
  isGPCR <- FALSE
}
subd_index= c(as.numeric(args[2])) # subdirs choose index
fc_index= c(as.numeric(args[3])) # feature choose index

cutoff = 0.05 # to filter the features.
###------keep unchanged below----------###
### packages
library(methods)
library(glmnet)
library(data.table)
###library(ggplot2)
##library(optimx)
###library(doParallel)
###library(foreach)
###library(caret)
library(matrixStats)
library(Metrics)

### set dir, and load func.
setwd(jobdir)
##source(paste(ml2dqsar_dir,ml2dqsar_filenm,sep=""))
source(paste(ml2dqsar_dir,wiperfunc_filenm,sep=""))
source(paste(ml2dqsar_dir,mergdglm_filenm, sep=""))
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
    files = ldat$allfiles
    dirs = ldat$protnms
    oneGroup = list()
    len = list()
    oneGroup_xy = list()
    for (i in (1:length(files))){
      oneGroup[[i]] <- fread(as.character(files[i][[1]]), sep=',',header=TRUE)
      len[i] <- length(oneGroup[[i]]$aveaffinty)
      ## load data, take -log10 for y.
      oneGroup[[i]]$aveaffinty <- -log10(oneGroup[[i]]$aveaffinty)
      oneGroup_xy$X[[i]] <- subset(oneGroup[[i]], select = c(-V1,-aveaffinty))
      oneGroup_xy$y[[i]] <- subset(oneGroup[[i]], select = aveaffinty)
    }
    keepft <- bifeakeep(oneGroup_xy,cutoff)
    ##for(i in length(oneGroup_xy$y)){
    ## with = FALSE is needed since it's data.table.  f*ck!!!
    ## But it cannot change itself. f*ck!!!
    ##  oneGroup_xy$X[[i]] = oneGroup_xy$X[[i]][ ,keepft, with=FALSE]
    ##}
    ## scale data
    for (i in (1:length(files))){
      oneGroup_xy$X[[i]] <- scale(oneGroup_xy$X[[i]], center=TRUE, scale = FALSE)
      oneGroup_xy$y[[i]] <- scale(oneGroup_xy$y[[i]], center=TRUE, scale = FALSE)
    }
    ## filter based on compounds' number.
    cpi <- list()
    cpi$proteins <- list()
    cpi$rootMSE  <- c()
    j = 0
    for(i in 1:length(oneGroup_xy$y)){
      if(len[[i]]>=lower_limit & len[[i]] <= upper_limit){
        j = j + 1
        ## as.character is needed since dirs is factor. f*ck!!!
        cpi$proteins[j] = as.character(dirs[[i]])
        cpi$y[[j]] <- oneGroup_xy$y[[i]] # log of affinity
        ## after scale, not data.table. f*ck!!!
        cpi$X[[j]] <- oneGroup_xy$X[[i]][ ,keepft] # Matrix of the predictors.
      }
    }
    result_alldata_glm <- cv_merge_glm(cpi,nfold)
    save(result_alldata_glm,
         file=paste('alldata_glm_',subdirs[k],'_',features[f],'.data',
                    sep = ""))
  } ## for  subdirs
} ## for features
