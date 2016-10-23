### wiper for run ml-2dqsar.
### szu
### 161023

bifeakeep <- function(datalist,cutoff){
  ## for binary feature only in datalist$X
  ##   each column of X is feature.
  ## return as vector for index to keep.

  feanum <- ncol(datalist$X[1])
  samplenum <- 0
  bifeakeep <- 1:feanum
  colsmeans <- rep(0,feanum)
  for(prot in length(datalist$X)){
    colsmeans = colsmeans + colSums(datalist$X[[i]])
    samplenum <- samplenum + nrow(datalist$X[[i]])
  }
  colsmeans <- colsmeans/samplenum
  return bifeakeep[which((colsmeans > cutoff) & (colsmeans < 1-cutoff))]
}

loadfiles <- function(datadir,gpcr_nm,kinase_nm,isGPCR,f,k){
  if(isGPCR){
    subdirs <- gpcr_names[subd_index]
  } else{
    subdirs <- kinase_names[subd_index]
  }
  if(isGPCR){
    files <- list.files(path=paste(datadir,"GPCR",subdirs[k],sep=""),
                        full.names=TRUE,recursive=TRUE,
                        pattern=paste(".*",features[f],".csv",sep=""))
    dirs <- dir(path=paste(datadir,"GPCR",subdirs[k],sep=""),
                full.names=FALSE,no..=FALSE)
  } else {
    files <- list.files(path=paste(datadir,subdirs[k],sep=""),
                        full.names=TRUE,recursive=TRUE,
                        pattern=paste(".*",features[f],".csv",sep=""))
    dirs <- dir(path=paste(datadir,subdirs[k],sep=""),
                full.names = FALSE, no.. = FALSE)
  }
  loadfiles = data.frame(allfiles = files,
                         protnms = dirs)
  return(loadfiles)
}
