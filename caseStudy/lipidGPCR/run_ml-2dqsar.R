### run ml-2dqsar
### demo
### szu
### 2016-10-20

### packages
library(glmnet) # for ridge regression
library(data.table) # for load data quickly
library(ggplot2) # for draw figures
library(optimx) # for optimization
library(doParallel) # for parallel calculations
library(foreach)
library(caret) # for create cross validation samples
library(matrixStats) # for calculate matrix columnwise variance

setwd("c:/users/zhou/desktop/sx")
source("FuncPreComPiHierBayes.R") # for ml2dqsar


### use proteins with compound between (50 to 1000)

nfold = 5
lower_limit = 50
upper_limit = 1000
cpuNumber = 3

### read file

files <- list.files(path="c:/users/zhou/desktop/sx/GPCR/CPIs_GPCR  lipid-like ligand receptor",
                      full.names=TRUE,recursive=TRUE,
                      pattern=".*macc.csv")

dirs <- dir(path="c:/users/zhou/desktop/sx/GPCR/CPIs_GPCR  lipid-like ligand receptor",
                      full.names=FALSE,no..=FALSE)

oneGroup = list()
len = list()
oneGroup_scaled = list()
for (i in (1:length(files)))
{
  oneGroup[[i]] <- fread(files[i][[1]], sep=',',header=TRUE)

  len[i] <- length(oneGroup[[i]]$aveaffinty) #spelling error, col_name in the data is aveaffinty
  ### load data, take -log10 for y.
  
  oneGroup[[i]]$aveaffinty <- -log10(oneGroup[[i]]$aveaffinty)
  ### center data both (x and y) (scale data for some of phychem)
  ### y column
  
  oneGroup[[i]] <- scale(oneGroup[[i]], center=TRUE, scale = FALSE)
  
  oneGroup_scaled$X[[i]] <- subset(oneGroup[[i]], select = c(-V1,-aveaffinty))
  oneGroup_scaled$y[[i]] <- subset(oneGroup[[i]], select = aveaffinty)
}

### data structure
cpi <- list()
#proteins <- list()
j = 0
for(i in 1:length(oneGroup_scaled$y)){
  if(len[[i]]>=lower_limit & len[[1]] <= upper_limit){
    j = j + 1
    #proteins[j] = dirs[[i]]
    cpi$y[[j]] <- oneGroup_scaled$y[[i]] # log of affinity
    cpi$X[[j]] <- oneGroup_scaled$X[[i]] # Matrix of the covarietes
  }
}


cpi$omega_matrix = matrix(nrow=ncol(cpi$X[[1]]),ncol=nfold*length(cpi$y))
#cpi$fit = list()

for(i in 1:length(cpi$y)){

  ### label data for n-fold cross validation.
  flds <- createFolds(cpi$y[[i]],
                      k = nfold,
                      list=TRUE,
                      returnTrain = FALSE)
  
  for(j in 1:nfold){
    
    ### run every task for ridge regression with glmnet package.
    #### for every protein,run
    #### lambda chosen as nested cv.
    fit <- glmnet(x = cpi$X[[i]][-flds[[j]],],
                  y = cpi$y[[i]][-flds[[j]]],
                  lambda=cv.glmnet(cpi$X[[i]][-flds[[j]],],
                                   y = cpi$y[[i]][-flds[[j]]],intercept = FALSE)$lambda.min,
                  intercept = FALSE)
    omega <- coef(fit) # for each protein
    #### record predicted tests and omega.
    test_result <- predict(fit,
                           cpi$X[[i]][flds[[j]],])
    
    cpi$rootMSE[(i-1)*nfold+j] <- sqrt(mse(cpi$y[[i]][flds[[j]]],test_result))
    
    #cpi$fit[(i-1)*nfold+j] <- fit

    cpi$omega_matrix[,(i-1)*nfold+j] = omega[2:168]
  }
}

### estimate sigma based on the omega metrix.
#### if each row corresponds to one protein, then
sigma <- mean(colVars(cpi$omega_matrix))

for(i in 1:length(cpi$y))
{
  cpi$X[[i]] = t(cpi$X[[i]])
}


### run ml-2dqsar with personal code and sigma
#### set
#cl <- makeCluster(cpuNumber)
#registerDoparallel(cl)
#### run ml-2dqsar
testResult <- crossPreComPiHier(cpi,
                                nfold=nfold,
                                sigma2j = sigma,
                                sigma2dot = 100,
                                maxit = TRUE,
                                iteration = 1000,
                                parallel = FALSE)
### see the result.
#### load data
load('cpi.data')
load('testResult.data')

getRMSE <- function(rootMSE,fold,protnum){
  rmse_single <- c()
  for(i in 1:protnum){
    rmse_single <- c(rmse_single,mean(rootMSE[((i-1)*fold + 1):((i-1)*fold+fold)] ))
  }
  return(rmse_single)
}
rmse_single <- getRMSE(cpi$rootMSE,5,16)


getRMSEm <- function(mseMult,fold,protnum){
  tmp <- rep(0,protnum)
  for(i in 1:fold){
    tmp = tmp + mseMult[[i]]$mse
  }
  return(sqrt(tmp/fold))
}
rmse_multi <- getRMSEm(testResult,5,16)

getCompNum <- function(cpi,protnum){
  compNum <- c()
  for(i in 1:protnum){
    compNum <- c(compNum,length(cpi$y[[i]]))
  }
  return(compNum)
}
compNum <- getCompNum(cpi,16)

statSM <- data.frame(single = rmse_single,
                     multi = rmse_multi,
                     compNum = compNum)
