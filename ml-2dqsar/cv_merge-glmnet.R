### contrary to single-model, merge all data
### szu
### 2016-10-20

library(glmnet)

mergedatalist <- function(datalist){
  totalrow <- 0
  totalcol <- ncol(datalist[[1]]$X)
  for(i in length(datalist)){
    totalrow = totalrow + nrow(datalist[[i]]$X)
  }
  mergedm <- matrix(,nrow=totalrow,ncol=totalcol)
  mergedv <- rep(0.0,totalrow)
  trow <- 0
  for(i in length(datalist)){
    tmpl <- trow + nrow(datalist$X[[i]])
    mergedm[trow:tmpl,] <- datalist[[i]]$X
    mergedv[trow:tmpl] <- datalist[[i]]$y
    trow = tmpl + 1
  }
  finald <- data.frame(X = mergedm, y = mergedv)
  return(finald)
}

tp_merge_glm <- function(traind,testd,nfold=5){
  td <- mergedatalist(traind)
  tasknm <- rep(0.0, length(testd))
  msevec <- rep(0.0, tasknm)
  fit <- glmnet(x = td$X,
                y = td$y,
                lambda = cv.glmnet(x=td$X,y=td$y,interpret=FALSE,nfold=nfold)$lambda.min,
                intercept  = FALSE)
  for(i in tasknm){
    tpre <- predict(fit, testd[[i]]$X)
    msevec[i] <- mse(testd[[i]]$y,tpre)
  }
  omega <- coef(fit)[2:ncol(td)]
  finalr <- data.frame(omega = omega,
                       mse = msevec)
  return(finalr)
}

cv_merge_glm <- function(datalist, nfold){
  traincv <- list()
  testcv <- list()
  for(j in 1:nfold){
    traincv[[j]] <- list()
    testcv[[j]] <- list()
    traincv[[j]]$X <- list()
    traincv[[j]]$y <- list()
    testcv[[j]]$X <- list()
    testcv[[j]]$y <- list()
    for(i in 1:group){
      foldid <- sample(rep(seq(nfold),length=ncol(datalist$X[[i]])))
      which <- foldid==j
      traincv[[j]]$y[[i]] = datalist$y[[i]][!which]
      traincv[[j]]$X[[i]] = datalist$X[[i]][,!which]
      testcv[[j]]$X[[i]] = datalist$X[[i]][,which]
      testcv[[j]]$y[[i]] = datalist$y[[i]][which]
    } # end of i in group
  } # end of j in nfold
  finalr <- tp_merge_glm(traincv,testcv,nfold)
  return(finalr)
}
