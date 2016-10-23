### run ml-2dqsar
### szu
### 2016-10-20

### set parameter.
jobdir <- "D:/lab/sx/"
datadir <- "D:/lab/sx/"
ml2dqsar_dir <- "D:/lab/sx/"

ml2dqsar_filenm <- "FuncPreComPiHierBayes.R"
wiperfunc_filenm <- "wiperdata.R"

gpcr_names <- c("CPIs_GPCR_monoamine_receptor",
                "CPIs_GPCR_nucleotide-like_receptor",
                "CPIs_GPCR  lipid-like ligand receptor",
                "CPIs_Peptide GPCR")
kinase_names <- c("kinase_agc", "kinase_camk","kinase_ck1","kinase_cmgc",
                  "kinase_ste","kinase_tk","kinase_tkl")
fealist <- c("macc","phychem","fingerprint","morgan")
isbinaryf <- c(TRUE, FALSE, TRUE, TRUE)
nfold <- 5
lower_limit <- 50
upper_limit <- 1000
sigma2dot <- 100
iternum <- 1000

isGPCR <- TRUE
subd_index= c(3) # subdirs choose index
fc_index= c(4) # feature choose index

cutoff = 0.05 # to filter the features.
###------keep unchanged below----------###
### packages
library(glmnet)
library(data.table)
library(ggplot2)
library(optimx)
library(doParallel)
library(foreach)
library(caret)
library(matrixStats)
library(Metrics)

### set dir, and load func.
setwd(jobdir)
source(paste(ml2dqsar_dir,ml2dqsar_filenm,sep=""))
source(paste(ml2dqsar_dir,wiperfunc_filenm,sep=""))

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
    for(i in length(oneGroup_xy$y)){
      oneGroup_xy$X[[i]] = oneGroup_xy$X[[i]][ ,keepft]
    }
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
        cpi$proteins[j] = dirs[[i]]
        cpi$y[[j]] <- oneGroup_xy$y[[i]] # log of affinity
        cpi$X[[j]] <- oneGroup_xy$X[[i]] # Matrix of the covarietes
      }
    }
    ## check j == 0 ; continue
    if(j==0){next;}
    ## filter features, and update cpi
    cpi$omega_matrix = matrix(nrow=length(cpi$y),ncol=ncol(cpi$X[[1]]))
    for(i in 1:length(cpi$y)){
      ## label data for n-fold cross validation.
      flds <- createFolds(cpi$y[[i]],
                          k = nfold,
                          list=TRUE,
                          returnTrain = FALSE)
      tmp_mse = 0
      for(j in 1:nfold){
        fit <- glmnet(x = cpi$X[[i]][-flds[[j]],],
                      y = cpi$y[[i]][-flds[[j]]],
                      lambda=cv.glmnet(cpi$X[[i]][-flds[[j]],],
                                       y = cpi$y[[i]][-flds[[j]]],intercept = FALSE, nfold = nfold)$lambda.min,
                      intercept = FALSE)
        test_result <- predict(fit,
                               cpi$X[[i]][flds[[j]],])
        tmp_mse = tmp_mse + sqrt(mse(cpi$y[[i]][flds[[j]]],test_result))
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
    sigma <- mean(colVars(cpi$omega_matrix))
    cpi$sigma <- sigma
    for(i in 1:length(cpi$y)){
      cpi$X[[i]] = t(cpi$X[[i]])
    }
    ## run ml-2dqsar with personal code and sigma
    testResult <- crossPreComPiHier(cpi,
                                    nfold=nfold,
                                    sigma2j = sigma,
                                    sigma2dot = sigma2dot,
                                    maxit = TRUE,
                                    iteration = iternum,
                                    parallel = FALSE)
    save(cpi, file=paste('cpi_',subdirs[k],'_',features[f],'.data',sep=""))
    save(testResult, file=paste('testResult_',subdirs[k],'_',features[f],'.data',sep=""))
  }
}
