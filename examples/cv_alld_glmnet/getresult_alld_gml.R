## get result of all data comparison
## szu
## 20161113

## load packages
library(ggplot2)
library(reshape)
## directory of result.
job_dir_pre <- "/Users/wangchao/home/songpeng/git-recipes/ml-2dqsar/examples/"
ref_ml_dir_rela <- "result_gpcrkinase_eqsj_autosig_fpmo/"
alld_glmnet_dir_rela <- "cv_alld_glmnet/result/"
alld_result_file_pre <- "alldata_glm_"
singled_result_file_pre <- "cpi_"
mqsar_result_file_pre <- "testResult_"
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
    rmsem <- rmsem + sqrt(mqsard[[i]]$mse)
  }
  result <- data.frame(rmsem = rmsem/fold,
                       rmses = cpid$rootMSE,
                       protnm = unlist(cpid$proteins))
  return(result)
}
getrmse_alldglm <- function(alldglm){
  fold <- length(alldglm) - 1 # leave one for protnm
  protnum <- length(alldglm$protnm)
  rmsealld <- rep(0,protnum)
  for(i in 1:fold){
    rmsealld <- rmsealld + sqrt(alldglm[[i]]$mse)
  }
  result <- data.frame(rmsea = rmsealld/fold,
                       protnm = unlist(alldglm$protnm))
  return(result)
}
getfnm <- function(isgpcr,proind,ftind){
  if(isgpcr > 0){
    pronm <- gpcr_names[proind]
  } else{
    pronm <- kinase_names[proind]
  }
  ftnm <- fealist[ftind]
  afnm <- list()
  afnm$singledfn <- paste(job_dir_pre,
                            ref_ml_dir_rela,
                            singled_result_file_pre,
                            pronm,'_',ftnm,'.data',sep = "")
  afnm$mqsarfn <- paste(job_dir_pre,
                          ref_ml_dir_rela,
                          mqsar_result_file_pre,
                          pronm,'_',ftnm,'.data',sep = "")
  afnm$alldfn <- paste(job_dir_pre,
                         alld_glmnet_dir_rela,
                         alld_result_file_pre,
                         pronm,'_',ftnm,'.data',sep = "")
  return(afnm)
}
load_alld <- function(afnm){
  singled <- load_obj(afnm$singledfn)
  mqsard <- load_obj(afnm$mqsarfn)
  alld <- load_obj(afnm$alldfn)

  rmsems <- getrmse_ms(mqsard,singled)
  rmsea <- getrmse_alldglm(alld)

  allrmse <- merge(rmsems, rmsea, by = "protnm")
  return(allrmse)
}
## run
##-----
isgpcr <- 1
proind <- 4
ftind <- 3

afnm_g_p <- getfnm(isgpcr, proind, ftind)
armse_g_p <- load_alld(afnm_g_p)

##-----
isgpcr <- 0
proind <- 1
ftind <- 3

afnm_k_a <- getfnm(isgpcr, proind, ftind)
armse_k_a <- load_alld(afnm_k_a)
