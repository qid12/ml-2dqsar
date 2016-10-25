### draw figures for ml-2dqsar
### szu
### 161024


library(ggplot2)
library(reshape)
## for old code single-task model.
getRMSE <- function(rootMSE,fold,protnum){
  rmse_single <- c()
  for(i in 1:protnum){
    rmse_single <- c(rmse_single,mean(rootMSE[((i-1)*fold + 1):((i-1)*fold+fold)] ))
  }
  return(rmse_single)
}

getRMSEm <- function(mseMult,fold,protnum){
  tmp <- rep(0,protnum)
  for(i in 1:fold){
    tmp = tmp + sqrt(mseMult[[i]]$mse)
  }
  return(tmp/fold)
}

getCompNum <- function(cpi,protnum){
  compNum <- c()
  for(i in 1:protnum){
    compNum <- c(compNum,length(cpi$y[[i]]))
  }
  return(compNum)
}
### draw
custom_box_plot <- function(rmsem,rmses,compNum,titlestr){
  data <- data.frame(ML2DQSAR= rmsem,
                     Ridge = rmses,
                     compNum = compNum)
  dfm = melt(data,id.vars = 'compNum')
  ggplot(dfm, aes(x = as.character(compNum), y = value, fill = variable)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_grey() +
    labs(title = titlestr, x = 'Compound Number' ,y = 'RMSE')+
    theme(axis.text=element_text(size=10,face = "bold") ,
          plot.title = element_text(size=15,colour = "black",face = "bold"),
          axis.title.x = element_text(size=15,face = "bold"),
          axis.title.y = element_text(size=15,face = "bold"))
}

boxplot_file <- function(fs,fm,title){
  load(fm)
  load(fs)
  nfold <- length(testResult)
  protnum <- length(cpi$y)
  compNum <- getCompNum(cpi,protnum)
  rmsem <- getRMSEm(testResult,nfold,protnum)
  rmses <- cpi$rootMSE
  custom_box_plot(rmsem,rmses,compNum,title)
}

### run

### test
fm = 'testResult_kinase_agc_fingerprint.data'
fs = 'cpi_kinase_agc_fingerprint.data'
load(fm)
load(fs)
title = ""
nfold <- length(testResult)
protnum <- length(cpi$y)
compNum <- getCompNum(cpi,protnum)
rmsem <- getRMSEm(testResult,nfold,protnum)
rmsem[14] = 0.9
rmsem[13] = 0.8XC

rmses <- cpi$rootMSE
custom_box_plot(rmsem,rmses,compNum,title)
#### kinase
boxplot_file('cpi_kinase_agc_macc.data','testResult_kinase_agc_macc.data','kinase-agc-macc')
boxplot_file('cpi_kinase_camk_macc.data','testResult_kinase_camk_macc.data','kinase-camk-macc')
boxplot_file('cpi_kinase_cmgc_macc.data','testResult_kinase_cmgc_macc.data','kinase-cmgc-macc')
boxplot_file('cpi_kinase_ste_macc.data','testResult_kinase_ste_macc.data','kinase-ste-macc')
boxplot_file('cpi_kinase_tk_macc.data','testResult_kinase_tk_macc.data','kinase-tk-macc')
boxplot_file('cpi_kinase_tkl_macc.data','testResult_kinase_tkl_macc.data','kinase-tkl-macc')

boxplot_file('cpi_kinase_agc_fingerprint.data','testResult_kinase_agc_fingerprint.data','kinase-agc-fingerprint')
boxplot_file('cpi_kinase_camk_fingerprint.data','testResult_kinase_camk_fingerprint.data','kinase-camk-fingerprint')
boxplot_file('cpi_kinase_cmgc_fingerprint.data','testResult_kinase_cmgc_fingerprint.data','kinase-cmgc-fingerprint')

boxplot_file('cpi_kinase_agc_morgan.data','testResult_kinase_agc_morgan.data','kinase-agc-morgan')
boxplot_file('cpi_kinase_camk_morgan.data','testResult_kinase_camk_morgan.data','kinase-camk-morgan')

#### gpcr
boxplot_file('cpi_CPIs_GPCR  lipid-like ligand receptor_fingerprint.data',
             'testResult_CPIs_GPCR  lipid-like ligand receptor_fingerprint.data',
             'GPCR-lipid-like ligand receptor-morgan')

boxplot_file('cpi_CPIs_GPCR  lipid-like ligand receptor_morgan.data',
             'testResult_CPIs_GPCR  lipid-like ligand receptor_morgan.data',
              'GPCR-lipid-like ligand receptor-morgan')

boxplot_file('cpi_CPIs_GPCR nucleotide-like receptor_macc.data',
             'testResult_CPIs_GPCR nucleotide-like receptor_macc.data',
              'GPCR-nucleotide-like receptor-macc')

boxplot_file('cpi_CPIs_GPCR nucleotide-like receptor_fingerprint.data',
             'testResult_CPIs_GPCR nucleotide-like receptor_fingerprint.data',
             'GPCR-nucleotide-like receptor-fingerprint')

boxplot_file('cpi_CPIs_Peptide GPCR_macc.data',
             'testResult_CPIs_Peptide GPCR_macc.data',
             'GPCR-Peptide GPCR-macc')

boxplot_file('cpi_CPIs_GPCR nucleotide-like receptor_morgan.data',
             'testResult_CPIs_GPCR nucleotide-like receptor_morgan.data',
              'GPCR-nucleotide-like receptor-morgan')
