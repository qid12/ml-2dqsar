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
library(rowVars) # for calculate matrix columnwise variance
source("FuncPreComPiHierBayes.R") # for ml2dqsar

### use proteins with compound between (50 to 1000)

### read file
setwd("your data directory.")
oneGroup <- fread("everyDataFile",sep=',',header=TRUE)

### load data, take -log10 for y.
data$affinity <- -log10(data$affinity)
### center data both (x and y) (scale data for some of phychem)
### y column
data <- scale(data,center=TRUE, scale = FALSE)


### data structure
cpi$y[[i]]  # ith protein's affinities, vector
cpi$X[[i]]  # ith protien's features, matrix (n * p), big X for matrix.
cpi$flds[[i]] # ith group label for cv
cpi$coef[[i]] # ith protein result for omege
cpi$coefm[[i]] # ith protein result for omega modified by ml2qsar.
cpi$result[[i]] # ith protein result for ridge prediction.
cpi$resultm[[i]] # ith protein result for ml2qsar.
cpi$sigma # sigma parameter for ml2qsar.
cpi$sigmadot <-  100


### label data for n-fold cross validation.
flds <- createFolds(data$affnity,k = nfold,list=TRUE,returnTrain = FALSE)
data$matrix[-flds[[i]],] # get the train sample.
data$matrix[flds[[1]],] # get the test sample.
#### another way for n-fold sample
flds <- cut(seq(1,ncol(data)),breaks = nfold, labels = FALSE)
testIndexes <- which(flds==i,arr.ind=TRUE)
testData <- data[testIndex,]
trainData <- data[-testIndex,]

### run every task for ridge regression with glmnet package.
#### for every protein,run
fit <- glmnet(x = data$matrix[train, ], y = data$affinity[train],
   lambda=cv.glmnet(data$matrix[train,],data$afffinity[train])$lambda.min)
omega <- coef(fit) # for each protein
#### record predicted tests and omega.
test_result <- predict(fit,data$matrix[test,])
rootMSE <- sqrt(mse(y[test],test_result))


### estimate sigma based on the omega metrix.
#### if each row corresponds to one protein, then
sigma <- mean(colVars(omega_matrix))


### run ml-2dqsar with personal code and sigma
#### set
cl <- makeCluster(cpuNumber)
registerDoparallel(cl)
#### run ml-2dqsar
testResult <- crossPreComPiHier(cpi,
                                nfold=foldnumber,
                                sigma2j = sigma,
                                sigma2dot = 100,
                                maxit = TRUE,
                                iteration = 1000,
                                parallel = TRUE)
