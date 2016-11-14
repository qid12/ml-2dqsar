### contrary to single-model, merge all data
### szu
### 2016-11-07
### need glmnet package.
mergedatalist <- function(datalist){
  totalrow <- 0
  totalcol <- ncol(datalist$X[[1]])
  tasknm <- length(datalist$X)
  for(i in 1:tasknm){
    totalrow = totalrow + nrow(datalist$X[[i]])
  }
  mergedm <- matrix(,nrow=totalrow,ncol=totalcol)
  mergedv <- rep(0.0,totalrow)
  trow <- 1
  for(i in 1:tasknm){
    tmpl <- trow + nrow(datalist$X[[i]]) - 1
    mergedm[trow:tmpl,] <- datalist$X[[i]]
    mergedv[trow:tmpl] <- datalist$y[[i]]
    trow = tmpl + 1
  }
  finald <- list(X = mergedm, y = mergedv)
  return(finald)
}

tp_merge_glm <- function(traind,testd,nfold=5){
  td <- mergedatalist(traind)
  tasknm <- length(traind$X)
  msevec <- rep(0.0, tasknm)
  fit <- glmnet(x = td$X,
                y = td$y,
                lambda = cv.glmnet(x=td$X,y=td$y,intercept =FALSE,nfold=nfold)$lambda.min,
                intercept  = FALSE)
  for(i in 1:tasknm){
    tpre <- predict(fit, testd$X[[i]])
    msevec[i] <- mse(testd$y[[i]],tpre)
  }
  omega <- coef(fit)
  finalr <- list(omega = omega[2:ncol(td$X)], mse = msevec)
  return(finalr)
}

cv_merge_glm <- function(datalist, nfold){
  tasknm <- length(datalist$X)
  traincv <- list()
  testcv <- list()
  for(j in 1:nfold){
    traincv[[j]] <- list()
    testcv[[j]] <- list()
    traincv[[j]]$X <- list()
    traincv[[j]]$y <- list()
    testcv[[j]]$X <- list()
    testcv[[j]]$y <- list()
    for(i in 1:tasknm){
      foldid <- sample(rep(seq(nfold),length=nrow(datalist$X[[i]])))
      which <- foldid==j
      traincv[[j]]$y[[i]] = datalist$y[[i]][!which]
      traincv[[j]]$X[[i]] = datalist$X[[i]][!which, ]
      testcv[[j]]$X[[i]] = datalist$X[[i]][which, ]
      testcv[[j]]$y[[i]] = datalist$y[[i]][which]
    } # end of i in group
  } # end of j in nfold
  ## init finalr.
  finalr <- list()
  for(i in 1:nfold){
    finalr[[i]] <- list()
    finalr[[i]]$omega <- rep(0.0, ncol(datalist$X[[1]]))
    finalr[[i]]$mse <- rep(0.0, tasknm)
  }
  ## record every fold result
  for(i in 1:nfold){
    tmp <- tp_merge_glm(traincv[[i]],testcv[[i]],nfold)
    finalr[[i]]$omega <- tmp$omega
    finalr[[i]]$mse <- tmp$mse
  }
  finalr$protnm <- datalist$proteins
  return(finalr)
}
