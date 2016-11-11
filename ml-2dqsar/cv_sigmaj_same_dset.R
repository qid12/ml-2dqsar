### run ml-2dqsar
### szu
### 2016-10-20

### accept args
args <- commandArgs(trailingOnly = TRUE)
### set parameter.
jobdir <- "/data/home/szu/Lab/ml-2dqsar/cvk-tkl/result/"
datadir <- "/data/home/szu/Lab/ml-2dqsar/"
ml2dqsar_dir <- "/data/home/szu/Lab/ml-2dqsar/cvk-tkl/"

ml2dqsar_filenm <- "FuncPreComPiHierBayes.R"
wiperfunc_filenm <- "wiperdata.R"
nfold <- 5
lower_limit <- 30
upper_limit <- 500
sigma2dot <- 100
iternum <- 1000

if(as.numeric(args[1]) > 0){
  isGPCR <- TRUE
} else{
  isGPCR <- FALSE
}
subd_index <- c(as.numeric(args[2])) # subdirs choose index
fc_index <- c(as.numeric(args[3])) # feature choose index
sigma2j <- as.numeric(args[4]) # sigmaj size
cutoff = 0.05 # to filter the features.

###------keep unchanged below----------###
### packages
library(methods)
library(glmnet)
library(data.table)
###library(ggplot2)
library(optimx)
###library(doParallel)
###library(foreach)
###library(caret)
library(matrixStats)
library(Metrics)

### set dir, and load func.
setwd(jobdir)
source(paste(ml2dqsar_dir,ml2dqsar_filenm,sep=""))
source(paste(ml2dqsar_dir,wiperfunc_filenm,sep=""))
## subdirs, fealist, and isbinaryf.
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
                  "kinase_tkl",
                  "kinase_stt_tkl",
                  "kinase_all_tkl")

fealist <- c("macc","phychem","fingerprint","morgan")
isbinaryf <- c(TRUE, FALSE, TRUE, TRUE)
### load features.
features <- fealist[fc_index]
isbinaryfc <- isbinaryf[fc_index]
if(isGPCR){
  subdirs <- gpcr_names[subd_index]
} else{
  subdirs <- kinase_names[subd_index]
}

### run
for(f in 1:length(features)){# f for feature
  for(k in 1:length(subdirs)){ # k for subdirs.
    ## load data.
    ldat = loadfiles(datadir,gpcr_nm, kinase_nm,isGPCR,f,k)
    files = ldat$allfiles
    dirs = ldat$protnms
    oneGroup = list()
    len = list()
    oneGroup_xy = list()
    for (i in (1:length(files))){
      oneGroup[[i]] <- fread(as.character(files[i][[1]]), sep=',',header=TRUE)
      len[i] <- length(oneGroup[[i]]$aveaffinty)
      ## load data, take -log10 for y.
      oneGroup[[i]]$aveaffinty <- -log10(oneGroup[[i]]$aveaffinty)
      oneGroup_xy$X[[i]] <- subset(oneGroup[[i]], select = c(-V1,-aveaffinty))
      oneGroup_xy$y[[i]] <- subset(oneGroup[[i]], select = aveaffinty)
    }
    keepft <- bifeakeep(oneGroup_xy,cutoff)
    ##for(i in length(oneGroup_xy$y)){
    ## with = FALSE is needed since it's data.table.  f*ck!!!
    ## But it cannot change itself. f*ck!!!
    ##  oneGroup_xy$X[[i]] = oneGroup_xy$X[[i]][ ,keepft, with=FALSE]
    ##}
    ## scale data
    for (i in (1:length(files))){
      oneGroup_xy$X[[i]] <- scale(oneGroup_xy$X[[i]], center=TRUE, scale = FALSE)
      oneGroup_xy$y[[i]] <- scale(oneGroup_xy$y[[i]], center=TRUE, scale = FALSE)
    }
    ## filter based on compounds' number.
    cpi <- list()
    cpi$proteins <- list()
    cpi$rootMSE  <- c()
    j = 0
    for(i in 1:length(oneGroup_xy$y)){
      if(len[[i]]>=lower_limit & len[[i]] <= upper_limit){
        j = j + 1
        ## as.character is needed since dirs is factor. f*ck!!!
        cpi$proteins[j] = as.character(dirs[[i]])
        cpi$y[[j]] <- oneGroup_xy$y[[i]] # log of affinity
        ## after scale, not data.table. f*ck!!!
        cpi$X[[j]] <- oneGroup_xy$X[[i]][ ,keepft] # Matrix of the predictors.
      }
    }
    ## check j == 0 ; continue
    ## if(j==0){next;}
    ## filter features, and update cpi
    cpi$omega_matrix = matrix(nrow=length(cpi$y),ncol=ncol(cpi$X[[1]]))
    for(i in 1:length(cpi$y)){
      ## label data for n-fold cross validation.
      ## flds <- createFolds(cpi$y[[i]],
      ##                     k = nfold,
      ##                     list=TRUE,
      ##                     returnTrain = FALSE)
      flds <- cut(seq(1,length(cpi$y[[i]])),breaks=nfold,labels=FALSE)
      tmp_mse = 0
      for(j in 1:nfold){
        testindex <- which(flds==j,arr.ind=TRUE)
        fit <- glmnet(x = cpi$X[[i]][-testindex,],
                      y = cpi$y[[i]][-testindex],
                      lambda=cv.glmnet(cpi$X[[i]][-testindex,],
                                       y = cpi$y[[i]][-testindex],
                                       intercept = FALSE, nfold = nfold)$lambda.min,
                      intercept = FALSE)
        test_result <- predict(fit,
                               cpi$X[[i]][testindex,])
        tmp_mse = tmp_mse + sqrt(mse(cpi$y[[i]][testindex],test_result))
      } # end of nfold for j.
      cpi$rootMSE[[i]] <- tmp_mse /nfold
      fit <- glmnet(x = cpi$X[[i]],
                    y = cpi$y[[i]],
                    lambda=cv.glmnet(cpi$X[[i]],y = cpi$y[[i]],intercept = FALSE, nfold = nfold)$lambda.min,
                    intercept = FALSE)
      omega <- coef(fit)
      cpi$omega_matrix[i,] = t(omega[2:nrow(omega)])
    } # end of protein num for i.
    ## estimate sigma based on the omega metrix.
    ## if each row corresponds to one protein, then
    ## only use features not zero
    ##sigma2j <- mean(colVars(cpi$omega_matrix))
    cpi$sigma2j <- sigma2j
    for(i in 1:length(cpi$y)){
      cpi$X[[i]] = t(cpi$X[[i]])
    }
    ## run ml-2dqsar with personal code and sigma
    mresult <- crossPreComPiHier(cpi,
                                    nfold=nfold,
                                    sigma2j = sigma2j,
                                    sigma2dot = sigma2dot,
                                    maxit = TRUE,
                                    iteration = iternum,
                                    parallel = FALSE)
    save(cpi, file=paste('cpi_',subdirs[k],'_',features[f],'.data',sep=""))
    save(mresult, file=paste('mresult_',subdirs[k],'_',features[f],'_',sigma2j,'.data',sep=""))
  }
}
