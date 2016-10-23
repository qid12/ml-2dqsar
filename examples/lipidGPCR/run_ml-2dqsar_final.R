### run ml-2dqsar
### szu
### 2016-10-20

### set parameter.
jobdir <- "D:/lab/TransferLearning/TransferModel/Code/data_161023/"
datadir <- "D:/lab/TransferLearning/TransferModel/Code/data_161023/"
ml2dqsar_dir <- "D:/lab/TransferLearning/TransferModel/Code/data_161023/"

ml2dqsar_filenm <- "FuncPreComPiHierBayes.R"
gpcr_names <- c("CPIs_GPCR_monoamine_receptor",
                "CPIs_GPCR_nucleotide-like_receptor",
                "CPIs_GPCR  lipid-like ligand receptor",
                "CPIs_Peptide GPCR")
kinase_names <- c("kinase_agc", "kinase_camk","kinase_ck1","kinase_cmgc",
                  "kinase_ste","kinase_tk","kinase_tkl")
fealist <- c("macc","phychem","fingerprint","morgan")

nfold <- 5
lower_limit = 50
upper_limit = 1000
### packages
library(glmnet) # for ridge regression
library(data.table) # for load data quickly
library(ggplot2) # for draw figures
library(optimx) # for optimization
library(doParallel) # for parallel calculations
library(foreach)
library(caret) # for create cross validation samples
library(matrixStats) # for calculate matrix columnwise variance
library(Metrics) 
### set dir
setwd(jobdir)
source(paste(ml2dqsar_dir,ml2dqsar_filenm,sep="")) # under jobdir.
### use proteins with compound between (50 to 1000)
proclass <- "GPCR" # if kinase set it proclass <- ""
subdirs <- gpcr_names[1,2]
features <- fealist[1,2,3]
### run
for(f in 1:length(features)){# f for feature
  for(k in 1:length(subdirs)) # k for subdirs
  {
    ## read file
    files <- list.files(path=paste(jobdir,subdirs[k],sep=""),
                        full.names=TRUE,recursive=TRUE,
                        pattern=paste(".*",features[f],".csv",sep=""))
    ## record protei names
    if(proclass){
        dirs <- dir(path=paste("D:/lab/TransferLearning/TransferModel/Code/data_161023/GPCR/",subdirs[k],sep=""),
                    full.names=FALSE,no..=FALSE)
    }
    dirs <- dir(path=paste("D:/lab/TransferLearning/TransferModel/Code/data_161023/GPCR/",subdirs[k],sep=""),
                full.names=FALSE,no..=FALSE)

    oneGroup = list()
    len = list()
    oneGroup_scaled = list()
    for (i in (1:length(files)))
    {
      oneGroup[[i]] <- fread(files[i][[1]], sep=',',header=TRUE)

      len[i] <- length(oneGroup[[i]]$aveaffinty) #spelling error, col_name in the data is aveaffinty

      ## load data, take -log10 for y.
      oneGroup[[i]]$aveaffinty <- -log10(oneGroup[[i]]$aveaffinty)

      ## center data both (x and y) (scale data for some of phychem)
      ## y column
      oneGroup[[i]] <- scale(oneGroup[[i]], center=TRUE, scale = FALSE)

      oneGroup_scaled$X[[i]] <- subset(oneGroup[[i]], select = c(-V1,-aveaffinty))
      oneGroup_scaled$y[[i]] <- subset(oneGroup[[i]], select = aveaffinty)
    }

    ## data structure
    cpi <- list()
    cpi$proteins <- list()
    j = 0
    for(i in 1:length(oneGroup_scaled$y)){
      if(len[[i]]>=lower_limit & len[[i]] <= upper_limit){
        j = j + 1
        cpi$proteins[j] = dirs[[i]]
        cpi$y[[j]] <- oneGroup_scaled$y[[i]] # log of affinity
        cpi$X[[j]] <- oneGroup_scaled$X[[i]] # Matrix of the covarietes
      }
    }
    ### check j == 0 ; continue
    cpi$omega_matrix = matrix(nrow=ncol(cpi$X[[1]]),ncol=nfold*length(cpi$y))
    for(i in 1:length(cpi$y)){
      ## label data for n-fold cross validation.
      flds <- createFolds(cpi$y[[i]],
                          k = nfold,
                          list=TRUE,
                          returnTrain = FALSE)

      for(j in 1:nfold){
        ## run every task for ridge regression with glmnet package.
        ## for every protein,run
        ## lambda chosen as nested cv.
        fit <- glmnet(x = cpi$X[[i]][-flds[[j]],],
                      y = cpi$y[[i]][-flds[[j]]],
                      lambda=cv.glmnet(cpi$X[[i]][-flds[[j]],],
                                       y = cpi$y[[i]][-flds[[j]]],intercept = FALSE)$lambda.min,
                      intercept = FALSE)
        omega <- coef(fit) # for each protein

        ## record predicted tests and omega.
        test_result <- predict(fit,
                               cpi$X[[i]][flds[[j]],])

        cpi$rootMSE[(i-1)*nfold+j] <- sqrt(mse(cpi$y[[i]][flds[[j]]],test_result))
        cpi$omega_matrix[,(i-1)*nfold+j] = omega[2:nrow(omega)]
      }
    }
    ## estimate sigma based on the omega metrix.
    ## if each row corresponds to one protein, then
    ## only use features not zero
    sigma <- mean(colVars(cpi$omega_matrix))
    cpi$sigma <- sigma
    for(i in 1:length(cpi$y)){
      cpi$X[[i]] = t(cpi$X[[i]])
    }

    ## run ml-2dqsar with personal code and sigma
    ##cl <- makeCluster(cpuNumber)
    ##registerDoparallel(cl)
    ## run ml-2dqsar
    testResult <- crossPreComPiHier(cpi,
                                    nfold=nfold,
                                    sigma2j = sigma,
                                    sigma2dot = 100,
                                    maxit = TRUE,
                                    iteration = 1000,
                                    parallel = FALSE)
    save(cpi, file=paste('cpi_',subdirs[k],'_',features[f],'.data',sep=""))
    save(testResult, file=paste('testResult_',subdirs[k],'_',features[f],'.data',sep=""))
  }
}
