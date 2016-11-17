### get result of sigma2dot on different relative values.
### szu
### 20161116

## set dir
job_dir_pre <- "/Users/wangchao/home/songpeng/git-recipes/ml-2dqsar/examples/cv_sigma2dot_sigma2jsame/"
result_subdir <- "result/"
mqsar_result_file_pre <- "mresult_"
singled_result_file_pre <- "cpi_"
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
sigma2dot_ref <- c(0.1,0.2,0.4,0.8,1,2,4,8,16,32,100,500,1000)
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
  floatmatch <- gregexpr("[0-9]+\\.{0,}[0-9]{0,}",mfiles)
  ## pairs of sigma2j and sigma2dot
  sigma2dot_tmp <- as.numeric(unlist(regmatches(mfiles,floatmatch)))
  ## get sigma2dot
  sigma2dot_real_ind <- seq(2,length(sigma2dot_tmp),2)
  sigma2dot_real <- round(sigma2dot_tmp[sigma2dot_real_ind],3)
  ## sort refsize
  j <- 1
  refsize <- rep(0.0,length(mfiles))
  for(i in 1:length(mfiles)){
    refsize[i] <- round(sigma2dot_tmp[j+1]/sigma2dot_tmp[j],1)
    j <- j+2
  }
  refsizeind <- rank(refsize)
  for(i in refsizeind){
    mfnm <- paste(job_dir_pre,
                  result_subdir,
                  mfiles[1],sep="")
    mqsard <- load_obj(mfnm)
    temp <- getrmse_ms(mqsard,singled)
    colnm <- paste(sigma2dot_real[i],"_",refsize[i],sep="")
    result[colnm] <- temp$rmsem
  }
  return(result)
}
## run
isgpcr <- 1
proind <- 4
ftind <- 3
armse_g_p <- load_alld(isgpcr,proind,ftind)
