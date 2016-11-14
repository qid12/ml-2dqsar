## get result of all data comparison
## szu
## 20161113

library(ggplot2)
library(reshape)

load_obj <- function(filenm){
  ## load rdata to new variable.
  env <- new.env()
  nm <- load(filenm,env)[1]
  env[[nm]]
}

job_dir_pre <- "/Users/wangchao/home/songpeng/git-recipes/ml-2dqsar/examples/"
ref_ml_dir_rela <- "result_gpcrkinase_eqsj_autosig_fpmo"
alld_glmnet_dir_rela <- "cv_all_glmnet/result/"
