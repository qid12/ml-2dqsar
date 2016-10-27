### wiper for run ml-2dqsar.
### szu
### 161023

loadfiles <- function(datadir,gpcr_nm,kinase_nm,isGPCR,f,k){
  if(isGPCR){
    subdirs <- gpcr_names[subd_index]
  } else{
    subdirs <- kinase_names[subd_index]
  }
  if(isGPCR){
    files <- list.files(path=paste(datadir,"GPCR/",subdirs[k],sep=""),
                        full.names=TRUE,recursive=TRUE,
                        pattern=paste(".*",features[f],".csv",sep=""))
    dirs <- dir(path=paste(datadir,"GPCR/",subdirs[k],sep=""),
                full.names=FALSE,no..=FALSE)
  } else {
    files <- list.files(path=paste(datadir,"kinase/",subdirs[k],sep=""),
                        full.names=TRUE,recursive=TRUE,
                        pattern=paste(".*",features[f],".csv",sep=""))
    dirs <- dir(path=paste(datadir,'kinase/',subdirs[k],sep=""),
                full.names = FALSE, no.. = FALSE)
  }
  loadfiles = data.frame(allfiles = files,
                         protnms = dirs)
  return(loadfiles)
}

getpredata <- function(filesl,){
  files <- filesl$allfiles
  dirs <- filesl$protnums
  oneGroup = list()
  len = list()
  oneGroup_xy = list()
  for (i in seq(files)){
    oneGroup[[i]] <- fread(as.character(files[i][[1]]), sep=',',header=TRUE)
    len[i] <- length(oneGroup[[i]]$aveaffinty)
    ## load data, take -log10 for y.
    oneGroup[[i]]$aveaffinty <- -log10(oneGroup[[i]]$aveaffinty)

    oneGroup_xy$X[[i]] <- subset(oneGroup[[i]], select = c(-V1,-aveaffinty))
    ## oneGroup_xy$X[[i]] <- scale(oneGroup_xy$X[[i]], center = TRUE, scale = FALSE)
    ## oneGroup_xy$xcolm[[i]] <- attributes(oneGroup_xy$X[[i]])$`scaled:center`

    oneGroup_xy$y[[i]] <- subset(oneGroup[[i]], select = aveaffinty)
    ## oneGroup_xy$y[[i]] <- scale(oneGroup_xy$y[[i]],center = TRUE, scale = FALSE)
    ## oneGroup_xy$ym[[i]] <- attributes(oneGroup_xy$y[[i]])x$`scaled:center`
  }
  return(oneGroup_xy)
}

bifeakeep <- function(datalist,cutoff){
  ## for binary feature only in datalist$X
  ##   each column of X is feature.
  ## return as vector for index to keep.
  feanum <- ncol(datalist$X[[1]])
  samplenum <- 0
  bifeaneed <- 1:feanum
  colsmeans <- rep(0,feanum)
  for(prot in seq(datalist$X)){
    colsmeans = colsmeans + colSums(datalist$X[[prot]])
    samplenum <- samplenum + nrow(datalist$X[[prot]])
  }
  colsmeans <- colsmeans/samplenum
  return(bifeaneed[which( (colsmeans > cutoff) & (colsmeans < (1-cutoff))) ])
}

getcpid <- function(oneGroup_xy, isbinaryf, cutoff, lower_limit, upper_limit){
  if(isbinaryf){
    keepft <- bifeakeep(oneGroup_xy,cutoff)

  } else{
    keepft <- seq(ncol(oneGroup_xy$X[[1]]))
  }

  for(i in seq(filse)){
    oneGroup_xy$X[[i]] <- scale(oneGroup_xy$X[[i]], center = TRUE, scale = FALSE)
    oneGroup_xy$xcolm[[i]] <- attributes(oneGroup_xy$X[[i]])$`scaled:center`

    oneGroup_xy$y[[i]] <- scale(oneGroup_xy$y[[i]],center = TRUE, scale = FALSE)
    oneGroup_xy$ym[[i]] <- attributes(oneGroup_xy$y[[i]])$`scaled:center`
  }

  cpi <- list()
  cpi$proteins <- list()
  cpi$rootMSE  <- c()
  j = 0
  for(i in seq(oneGroup_xy$y)){
    if(len[[i]]>=lower_limit & len[[i]] <= upper_limit){
      j = j + 1
      ## as.character is needed since dirs is factor. f*ck!!!
      cpi$proteins[j] = as.character(dirs[[i]])
      cpi$y[[j]] <- oneGroup_xy$y[[i]] # log of affinity
      cpi$ym[[j]] <- oneGroup_xy$ym[[i]]
      ## after scale, not data.table. f*ck!!!
      cpi$X[[j]] <- t(oneGroup_xy$X[[i]][ ,keepft]) # Matrix of the predictors.
      cpi$xcolm[[j]] <- oneGroup_xy$xcolm[[i]]
    }
  }
  try(if(j==0) stop('No data between "lower_limit" and "upper_limit".'))
  return(cpi)
}

get_singlpar_sigma <- function(cpi, nfold, isaddstd,isusemin,useminv){
  omega_matrix  = matrix(nrow=length(cpi$y),
                         ncol = nrow(cpi$X[[1]]))
  for(i in seq(cpi$y)){
    fit <- glmnet(x=cpi$X,
                  y=cpi$y,
                  lambda=cv.glmnet(cpi$X, cpi$y, intercept=FALSE,nfold=nfold)$lambda.min,
                  intercept = FALSE)
    omega <- coef(fit)
    omega_matrix[i, ] = t(omega[2:nrow(omega)])
  }
  sigma <- mean(colVars(omega_matrix))
  if(isaddstd){
    sigma <- sigma+2*sqrt(var(colVars(omega_matrix)))
  }
  if(isusemin){
    sigma <- max(c(sigma,useminv))
  }
  result <- list()
  result$singlepar <- omega_matrix
  result$sigma <- sigmae
  return(result)
}

get_ml2dqsarpar <- function(cpi,sigma2j,sigma2dot){
  p <- nrow(cpi$X[[1]])
  omega_matrix = matrix(nrow=length(cpi$y),
                        ncol=p)
  tmodel <- trainPreComPiHier(cpi$X,
                              sigma2j,
                              sigma2dot,
                              maxit,
                              iteration)
  numm <- length(cpi$y)
  for(i in seq(numm)){
    froml = (i-1)*p + 1
    tor = i*p
    omega_matri[i,] <- as.numeric(unlist(tmodel[froml:tor]))
  }
  froml = numm*p+1
  tor  = froml + p
  omega_dot = as.numeric(unlist(tmodel[froml:tor]))
  sigma_y = as.numeric(unlist(tmodel[tor+1]))

  result <- list()
  result$omega_matrix <- omega_matrix
  result$omega_dot <- omega_dot
  result$sigma_y <- sigma_y
  return(result)
}
