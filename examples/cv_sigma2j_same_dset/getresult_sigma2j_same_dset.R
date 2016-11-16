### get result of sigma2j on directly different values.
### szu
### 20161116

## set dir
job_dir <- "/Users/wangchao/home/songpeng/git-recipes/ml-2dqsar/examples/cv_sigma2j_same_dset/"
result_subdir <- "result/"
mqsar_result_file_pre <- "mresult_"
singled_result_file_pre <- "cpi_"
wiperfunc_filenm <- "wiperdata.R"
## add libs.
source(paste(job_dir,wiperfunc_filenm,sep = ""))
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
