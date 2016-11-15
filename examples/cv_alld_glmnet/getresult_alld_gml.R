## get result of all data comparison
## szu
## 20161113

## directory of result.
job_dir_pre <- "/Users/wangchao/home/songpeng/git-recipes/ml-2dqsar/examples/"
ref_ml_dir_rela <- "result_gpcrkinase_eqsj_autosig_fpmo/"
alld_glmnet_dir_rela <- "cv_all_glmnet/result/"
alld_result_file_pre <- "alldata_glm"
singled_result_file_pre <- "cpi_"
mqsar_result_file_pre <- "testResult_"
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
  result <- data.frame(rmsem = rmsem/fold,
                       rmses = cpid$rootMSE,
                       protnm = unlist(cpid$proteins))
  return(result)
}
getrmse_alldglm <- function(alldglm){
  fold <- length(alldglm)
  protnum <- length(alldglm$protnm)
  rmsealld <- rep(0,protnum)
  for(i in 1:fold){
    rmsealld <- rmsealld + sqrt(alldglm[[i]]$mse)
  }
  result <- data.frame(rmsea = rmesealld/fold,
                       protnm = alldgml$protnm)
  return(result)
}
getfnm <- function(isgpcr,proind,ftind){
  if(isgpcr > 0){
    pronm <- gpcr_names[proind]
  } else{
    pronm <- kinase_names[proind]
  }
  ftnm <- fealist[ftind]
  fnmall <- list()
  fnmall$singledfn <- paste(job_dir_pre,
                            ref_ml_dir_rela,
                            singled_result_file_pre,'_',
                            pronm,'_',ftnm,'.data',sep = "")
  fnmall$mqsarfn <- paste(job_dir_pre,
                          ref_ml_dir_rela,
                          mqsar_result_file_pre,'_',
                          pronm,'_',ftnm,'.data',sep = "")
  fnmall$alldfn <- paste(job_dir_pre,
                         alld_glmnet_dir_rela,
                         alld_result_file_pre,'_',
                         pronm,'_',ftnm,'.data',sep = "")
  return(fnmall)
}
load_alld <- function(ftnmall){
  singled <- load_obj(ftnmall$singledfn)
  mqsard <- load_obj(ftnmall$mqsarfn)
  alld <- load_obj(ftnmall$alldfn)

  rmsems <- getrmse_ms(mqsard,singled)
  rmsea <- getrmse_alldglm(alld)

  allrmse <- merge(rmses, rmsea, by = "pronm")
  return(allrmse)
}
## run
