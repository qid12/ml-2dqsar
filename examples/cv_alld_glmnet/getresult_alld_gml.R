## get result of all data comparison
## szu
## 20161113
## load packages
library(ggplot2)
library(reshape)
## basic information.
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
## functions.
load_obj <- function(filenm){
  ## load rdata to new variable.
  env <- new.env()
  nm <- load(filenm,env)[1]
  env[[nm]]
}
getrmse_ms<- function(mqsard,cpid){
  protnum <- length(cpid$y)
  fold <- length(mqsard)
  rmsem <- rep(0,protnum)
  for(i in 1:fold){
    rmsem = rmsem + sqrt(mseMult[[i]]$mse)
  }
  result <- list()
  result$rmsem <- rmsem/fold
  result$rmses <- cpid$rootMSE
  result$protnm <- unlist(cpid$proteins)
  return(result)
}
getrmse_alldglm <- function(alldglm){
  fold <- length(alldglm)
  protnum <- length(alldglm$protnm)
  rmsealld <- rep(0,protnum)
  for(i in 1:fold){
    rmsealld <- rmsealld + sqrt(alldglm[[i]]$mse)
  }
  result <- list()
  result$rmsealld <- rmsealld/fold
  result$protnm <- alldgml$protnum
  return(result)
}
## directory of result.
job_dir_pre <- "/Users/wangchao/home/songpeng/git-recipes/ml-2dqsar/examples/"
ref_ml_dir_rela <- "result_gpcrkinase_eqsj_autosig_fpmo/"
alld_glmnet_dir_rela <- "cv_all_glmnet/result/"
alld_result_file_pre <- "alldata_glm"
singled_result_file_pre <- "cpi_"
mqsar_result_file_pre <- "testResult_"
