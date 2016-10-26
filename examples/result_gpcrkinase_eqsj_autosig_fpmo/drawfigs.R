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
fm = 'testResult_kinase_tk_fingerprint.data'
fs = 'cpi_kinase_tk_fingerprint.data'
load(fm)
load(fs)
title = ""
nfold <- length(testResult)
protnum <- length(cpi$y)
compNum <- getCompNum(cpi,protnum)
rmsem <- getRMSEm(testResult,nfold,protnum)
rmses <- cpi$rootMSE
custom_box_plot(rmsem,rmses,compNum,title)

### draw figures.
boxplot_file('cpi_kinase_agc_morgan.data',
             'testResult_kinase_agc_morgan.data',
             'kinase agc with morgan')
boxplot_file('cpi_kinase_agc_fingerprint.data',
             'testResult_kinase_agc_fingerprint.data',
             'kinase ag with pubchem')

boxplot_file('cpi_kinase_camk_morgan.data',
             'testResult_kinase_camk_morgan.data',
             'kinase camk morgan')
boxplot_file('cpi_kinase_camk_fingerprint.data',
             'testResult_kinase_camk_fingerprint.data',
             'kinase camk fingerprint')

## sigma = 0.0076 kept fp: 343
## performance is ok, should compared with learn in group
boxplot_file('cpi_kinase_tk_morgan.data',
             'testResult_kinase_tk_morgan.data',
             'kinase tk morgan')
boxplot_file('cpi_kinase_tk_fingerprint.data',
             'testResult_kinase_tk_fingerprint.data',
             'kinase tk fingerprint')


boxplot_file('cpi_kinase_cmgc_morgan.data',
             'testResult_kinase_cmgc_morgan.data',
             'kinase cmgc morgan')

## sigma = 0.0089 kept fp: 339
## performance is bad, should compared with larger variance
boxplot_file('cpi_kinase_cmgc_fingerprint.data',
             'testResult_kinase_cmgc_fingerprint.data',
             'kinase cmgc fingerprint')

boxplot_file('cpi_lipidlike_morgan.data',
             'testResult_lipidlike_morgan.data',
             'gpcr lipidlike morgan')
boxplot_file('cpi_lipidlike_fingerprint.data',
             'testResult_lipidlike_fingerprint.data',
             'lipidlike fingerprint')

boxplot_file('cpi_monoamine_morgan.data',
             'testResult_monoamine_morgan.data',
             'gpcr monoamine morgan')
boxplot_file('cpi_monoamine_fingerprint.data',
             'testResult_monoamine_fingerprint.data',
             'monoamine fingerprint')

boxplot_file('cpi_peptide_morgan.data',
             'testResult_peptide_morgan.data',
             'gpcr peptide morgan')
boxplot_file('cpi_peptide_fingerprint.data',
             'testResult_peptide_fingerprint.data',
             'peptide fingerprint')
