# TransferLearning Function Verision
# Songpeng Zu
# 2015-01-22

trainPreComPiHier <- function(datalist,sigma2j=2,sigma2dot=1,maxit =TRUE,iteration = 30000){
  # Inition 
  m <- length(datalist$y)
  p <- nrow(datalist$X[[1]])
  
  totalnum <- 0
  for(i in 1:m){
    totalnum <- totalnum + length(unlist(datalist$y[[i]]))
  }
  
  innerProduct <- list(Xy <- list(),XX <- list())
  for(i in 1:m){
    innerProduct$Xy[[i]] <- as.vector(datalist$X[[i]]%*%datalist$y[[i]])
    innerProduct$XX[[i]] <- datalist$X[[i]]%*%t(datalist$X[[i]])
  }
  
  theta0 <- c(rep(0.01,times=(m+1)*p),0.1)
  LB0 <- c(rep(-100,times=(m+1)*p),0.00001)
  
  # optimization target
  objectfunc <- function(theta){
    sum1 = 0
    sum2 = 0
    for(i in 1:m){
      tmp1 <- datalist$y[[i]]-as.vector(t(datalist$X[[i]])%*%theta[((i-1)*p+1):(i*p)])
      sum1 <- sum1 + sum(tmp1^2)
      
      tmp2 <- theta[((i-1)*p+1):(i*p)] - theta[(m*p+1):((m+1)*p)]
      sum2 <- sum2 + sum(tmp2^2)
    }
    
    sum <- -sum1*0.5/tail(theta,n=1) - sum2*0.5/sigma2j - sum(theta[((m*p)+1):((m+1)*p)]^2)*0.5/sigma2dot -
      log(2*pi*tail(theta,n=1))*0.5*totalnum
    return(sum) 
  }
  
  # gradfunc
  gradientfunc <- function(theta){
    grad <- c()
    tmp <- rep(0,times=p)
    tmpsum <- 0
    for(i in 1:m){
      grad <- c(grad,innerProduct$Xy[[i]]/tail(theta,n=1)+theta[(m*p+1):((m+1)*p)]/sigma2j
                - as.vector(innerProduct$XX[[i]]%*%theta[((i-1)*p+1):(i*p)])/tail(theta,n=1)
                - theta[((i-1)*p+1):(i*p)]/sigma2j)
      
      tmp <- tmp + theta[((i-1)*p+1):(i*p)]
      
      tmp1 <- datalist$y[[i]]-as.vector(t(datalist$X[[i]])%*%theta[((i-1)*p+1):(i*p)])
      tmpsum <- tmpsum + sum(tmp1^2)
    }
    grad <- c(grad,-m*theta[(m*p+1):((m+1)*p)]/sigma2j+tmp/sigma2j-theta[(m*p+1):((m+1)*p)]/sigma2dot)
    grad <- c(grad,tmpsum*0.5/tail(theta,n=1)^2-totalnum*0.5/tail(theta,n=1))
    
    return(grad)
  }
  
  # This function needs the package optimx.
  tmp <- optimx(par = theta0,fn= objectfunc,gr= gradientfunc,lower=LB0,method="L-BFGS-B",
         control = list(maximize=maxit,maxit=iteration))
  return(tmp)
}

predictPreComPiHier <- function(predictopt,datalist,group){
  p <- nrow(datalist$X[[1]])
  resultpre <- list()
  m <- length(datalist$X)
  if(missing(group)){
    for(i in 1:m){
      resultpre$pre[[i]] <- as.vector(t(datalist$X[[i]])%*%as.numeric(unlist(predictopt[((i-1)*p+1):(i*p)])))
      if(!is.null(datalist$y)){
        tmparray <- resultpre$pre[[i]]-datalist$y[[i]]
        resultpre$mse[[i]] <- mean(tmparray^2)
        resultpre$var[[i]] <- var(tmparray)
        resultpre$rsquare[[i]] <- sum((resultpre$pre[[i]]-mean(datalist$y[[i]]))^2)/sum((datalist$y[[i]]-mean(datalist$y[[i]]))^2) 
        resultpre$y[[i]] <- datalist$y[[i]]
      }
    }  
  }
  else{
    i <- group
    resultpre$pre <- as.vector(t(datalist$X[[i]])%*%as.numeric(unlist(predictopt[((i-1)*p+1):(i*p)])))
    if(!is.null(datalist$y)){
      tmparray <- resultpre$pre-datalist$y[[i]]
      
      resultpre$mse <- mean(tmparray^2)
      resultpre$var <- var(tmparray)
      resultpre$rsquare <- sum((resultpre$pre[[i]]-mean(datalist$y[[i]]))^2)/sum((datalist$y[[i]]-mean(datalist$y[[i]]))^2)
      resultpre$y <- datalist$y[[i]]
    }
  }
  resultpre$converge <- predictopt$convcode
  resultpre$sigmay <- unlist(predictopt[(m+1)*p+1])
  resultpre$PreComPiHier <- predictopt
  
  return(resultpre)
}



crossPreComPiHier <- function(datalist,nfold = 10,sigma2j=2,sigma2dot=1,maxit =TRUE,iteration = 30000,parallel=FALSE){
  resultcross <- list()
  foldid <- list()
  group <- length(datalist$y)
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
    } 
  }
  
  if(parallel && require(foreach)){
    resultcross <- foreach(i = seq(nfold),.export=c("trainPreComPiHier","predictPreComPiHier"),.packages = "optimx") %dopar%{
      tmodel <- trainPreComPiHier(traincv[[i]],sigma2j,sigma2dot,maxit,iteration)     
      tpreresult <- predictPreComPiHier(tmodel,testcv[[i]])
    }
  }else{
    for(i in seq(nfold)){
      tmodel <- trainPreComPiHier(traincv[[i]],sigma2j,sigma2dot,maxit,iteration)     
      resultcross[[i]] <- predictPreComPiHier(tmodel,testcv[[i]])
    }
  }
  return(resultcross)
}





