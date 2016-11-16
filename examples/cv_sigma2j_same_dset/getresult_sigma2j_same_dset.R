### get result of sigma2j on directly different values.
### szu
### 20161116

## set dir
job_dir_pre <- "/Users/wangchao/home/songpeng/git-recipes/ml-2dqsar/examples/cv_sigma2j_same_dset/"
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
sigma2jvec <- c(0.05,0.1,0.2,0.4,0.8,1,2,4,16)
## functions
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
  for(i in 1:length(sigma2jvec)){
    mqsarfn <- paste(job_dir_pre,
                     result_subdir,
                     mqsar_result_file_pre,
                     pronm,'_',ftnm,'_',
                     sigma2jvec[i],'.data',sep = "")
    mqsard <- load_obj(mqsarfn)
    temp <- getrmse_ms(mqsard,singled)
    ## assign(paste('m',sigma2jvec[i],sep=""),temp$rmsem)
    ## eval(as.name(paste(something)))
    result[paste("m",sigma2jvec[i],sep="")] = temp$rmsem
  }
  return(result)
}
## run
isgpcr <- 1
proind <- 4
ftind <- 3
armse_g_p <- load_alld(isgpcr,proind,ftind)
