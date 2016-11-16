### get result of sigma2j on different relative values.
### szu
### 20161116

## set dir
job_dir_pre <- "/Users/wangchao/home/songpeng/git-recipes/ml-2dqsar/examples/cv_sigma2j_same_rset/"
result_subdir <- "result/"
mqsar_result_file_pre <- "mresult_"
singled_result_file_pre <- "cpi_"
wiperfunc_filenm <- "wiperdata.R"
## add libs.
source(paste(job_dir_pre,wiperfunc_filenm,sep = ""))
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
sigma2jvec_ref <- c(0.1,0.5,1,2,4,8)
## functions
load_obj <- function(filenm){
  ## load rdata to new variable.
  env <- new.env()
  nm <- load(filenm,env)[1]
  env[[nm]]
}
getrmse_ms<- function(mqsard,cpid){
  ## get rmse for single and mqsar results.
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
load_alld <- function(isgpcr,proind,ftind){
  if(isgpcr > 0){
    pronm <- gpcr_names[proind]
  } else{
    pronm <- kinase_names[proind]
  }
  ftnm <- fealist[ftind]
  singledfn <- paste(job_dir_pre,
                     result_subdir,
                     singled_result_file_pre,
                     pronm,'_',ftnm,'.data',sep = "")
  singled <- load_obj(singledfn)
  result <- data.frame(protnm = unlist(singled$proteins),
                       rmses = singled$rootMSE)
  mfiles <- list.files(path=paste(job_dir_pre,result_subdir,sep=""),
                       full.names = FALSE,no..=TRUE,
                       pattern=paste(mqsar_result_file_pre,
                                     pronm,".*",ftnm,".*",sep=""))
  floatmatch <- gregexpr("0\\.[0-9]+",mfiles)
  sigma2jvec_real <- as.numeric(unique(unlist(regmatches(mfiles,floatmatch))))
  for(i in 1:length(mfiles)){
    mqsard <- load_obj(mfiles[i])
    temp <- getrmse_ms(mqsard,singled)
    colnm <- paste("m",sigma2jvec_real[i],"_",
                   sigma2jvec_ref[i],sep="")
    result[colnm] <- temp$rmsem
  }
  return(result)
}
## run
isgpcr <- 1
proind <- 4
ftind <- 3
armse_g_p <- load_alld(isgpcr,proind,ftind)
