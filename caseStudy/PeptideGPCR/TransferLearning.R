# Transfer Learning.
# Songpeng Zu
# 2015-01-19

#-- Write the object function.
objectfunc <- function(theta){
  # totalnum: the total number of the samples among all the groups
  # theta, a numeric vector, includes w, w*, and sigma_y
  # need CPI (list structure), including result, a list of IC50s, predictor, a list of samples (p*n)
  # need sigma_j and sigma_*
  # m: group number
  # p: dimension of the covariets. 
  sum1 = 0
  sum2 = 0
  for(i in 1:m){
    tmp1 <- CPI$y[[i]]-as.vector(t(CPI$X[[i]])%*%theta[((i-1)*p+1):(i*p)])
    sum1 <- sum1 + sum(tmp1^2)
    
    tmp2 <- theta[((i-1)*p+1):(i*p)] - theta[(m*p+1):((m+1)*p)]
    sum2 <- sum2 + sum(tmp2^2)
  }
  
  sum <- -sum1*0.5/tail(theta,n=1) - sum2*0.5/sigma2j - sum(theta[((m*p)+1):((m+1)*p)]^2)*0.5/sigma2dot -
    log(2*pi*tail(theta,n=1))*0.5*totalnum
  return(sum) 
}

#-- write the deriviative function.

gradientfunc <- function(theta){
  # Here we need a new list called innerProduct, which contains Xy and XX for each groups.
  grad <- c()
  tmp <- rep(0,times=p)
  tmpsum <- 0
  for(i in 1:m){
    grad <- c(grad,innerProduct$Xy[[i]]/tail(theta,n=1)+theta[(m*p+1):((m+1)*p)]/sigma2j
              - as.vector(innerProduct$XX[[i]]%*%theta[((i-1)*p+1):(i*p)])/tail(theta,n=1)
              - theta[((i-1)*p+1):(i*p)]/sigma2j)
 
    tmp <- tmp + theta[((i-1)*p+1):(i*p)]
    
    tmp1 <- CPI$y[[i]]-as.vector(t(CPI$X[[i]])%*%theta[((i-1)*p+1):(i*p)])
    tmpsum <- tmpsum + sum(tmp1^2)
  }
  grad <- c(grad,-m*theta[(m*p+1):((m+1)*p)]/sigma2j+tmp/sigma2j-theta[(m*p+1):((m+1)*p)]/sigma2dot)
  grad <- c(grad,tmpsum*0.5/tail(theta,n=1)^2-totalnum*0.5/tail(theta,n=1))
  
  return(grad)
}

#-- Simulation
#--optimization 
# Bounds constrain: sigma_y should be larger than 0.

#library(lbfgs)
#out.lbfgs <- lbfgs(objectfunc,gradientfunc,theta0)

#library(Rsolnp)
#solnp(theta0,fun = objectfunc,LB = c(-5,-5,-5,-5,-5,-5,0.001))

#-- Real Data
setwd("D:/lab/TransferLearning/TransferModel/Code")
library(data.table)
# load data 
CPIs.sub <- fread("CPIs_substructure_PeptideGPCR",sep="\t",header=TRUE)
setkeyv(CPIs.sub,c("compound","protein"))
CPIs.prop <- fread("CPIs_property_PeptideGPCR",sep="\t",header=TRUE)
setkeyv(CPIs.prop,c("compound","protein"))
CPIs.all <- merge(CPIs.prop,CPIs.sub,suffixes = c("",".y"))


sub.dyn.raw <- read.table("FeatureSelect_Substructure_PepetideGPCR",sep="\t",header=FALSE,
                          skip=1, quote="",colClasses=c("character","numeric"))
prop.dyn.raw <- read.table("FeatureSelect_Propter_PepetideGPCR",sep="\t",header=FALSE,
                           skip=1,quote="",colClasses=c("character","numeric"))

score <- c(prop.dyn.raw[,2],sub.dyn.raw[,2])
dyn.result <- data.frame(score)

row.names(dyn.result)<-c(prop.dyn.raw[,1],paste("Sub",1:881,sep=""))

feafilter <- function(dyn.result,cutoff){
  feafilter <- row.names(dyn.result)[which(dyn.result>=cutoff)]
}
subs.del <- c(paste("Sub",1:115,sep=""),paste("Sub",264:327,sep=""))
subs.select <- setdiff(feafilter(data.frame(c(sub.dyn.raw[,2]),row.names=paste("Sub",1:881,sep="")),1),subs.del)
detect.dyn <- colSums(CPIs.all[,subs.select,with=FALSE])/nrow(CPIs.all)
subs.del.freq <- c(subs.select[which(detect.dyn<=0.05)],subs.select[which(detect.dyn>=0.95)])
subs.select.feature <- setdiff(subs.select,subs.del.freq)
select.feature <- c(feafilter(data.frame(c(prop.dyn.raw[,2]),row.names=prop.dyn.raw[,1]),1),subs.select.feature)

CPIs.select.colname <- c("protein","compound","affinity","logP","MR","MW","TPSA",
                         select.feature)

CPIs.select <- CPIs.all[,CPIs.select.colname,with=FALSE]
# Add one constant column
CPIs.select <- cbind(CPIs.select,rep(1,times=nrow(CPIs.select)))
colnames(CPIs.select) <- c(CPIs.select.colname,"constant")
  
setkey(CPIs.select,"protein")
proteins <- unique(as.factor(CPIs.select$protein))

# Here we choose several proteins first.
prot.stat <- as.data.frame(table(as.factor(CPIs.select$protein)))
prot.stat.sort <- prot.stat[order(prot.stat$Freq,decreasing=TRUE),]

# we choose several proteins according to the numbers of compounds targeting them.
proteins <- c("CHEMBL4015","CHEMBL3974","CHEMBL3473","CHEMBL2049","CHEMBL252","CHEMBL4843")
proteins <- c("CHEMBL2049","CHEMBL252")
proteins <- c("CHEMBL252")
#Construct the CPIS data.
CPI <- list()
for(i in 1:length(proteins)){
  tmprows = which(CPIs.select[,1,with=FALSE]==as.character(proteins[i]))
  CPI$y[[i]] <- as.vector(unlist(log(CPIs.select[tmprows,3,with=FALSE]))) # log of affinity
  CPI$X[[i]] <- t(CPIs.select[tmprows,4:ncol(CPIs.select),with=FALSE]) # Matrix of the covarietes
}


m <- length(proteins)
p<- length(4:ncol(CPIs.select))

#totalnum <- ncol(CPIs.select) 
totalnum <- 0
for(i in 1:m){
  totalnum <- totalnum + length(unlist(CPI$y[[i]]))
}

sigma2j <- 2
sigma2dot <- 1

innerProduct <- list(Xy <- list(),XX <- list())
for(i in 1:m){
  innerProduct$Xy[[i]] <- as.vector(CPI$X[[i]]%*%CPI$y[[i]])
  innerProduct$XX[[i]] <- CPI$X[[i]]%*%t(CPI$X[[i]])
}

theta0 <- c(rep(0.01,times=(m+1)*p),0.1)
LB0 <- c(rep(-100,times=(m+1)*p),0.00001)
resultsolnp <- solnp(theta0,fun = objectfunc,LB = LB0)

library(optimx)
#resultoptimx <- optimx(par = theta0,fn= objectfunc,gr = gradientfunc,lower=LB0,method="L-BFGS-B")
resultoptimx <- optimx(par = theta0,fn= objectfunc,gr=gradientfunc,lower=LB0,method="L-BFGS-B",
                       control = list(maximize=TRUE,maxit=30000))
plotmatrix <- matrix(as.numeric(resultoptimx[1:1988]),nrow=p,byrow=FALSE)
nba_heatmap <- heatmap(plotmatrix, Rowv=NA, Colv=NA, col = cm.colors(256), scale="column", margins=c(3,20),cexRow=0.2)



#library(ggplot2)
#gp <- ggplot(as.data.frame(plotmatrix)) + geom_tile(colour = "white") 
#gp <- gp  + scale_fill_gradient(low = "white", high = "steelblue")

# Note: we need center the data or add const 1 as one predictor.

#-- Cross Validation Design
# First Strategy: partition each group and label as 1 to n-fold and compare the RMSE (average mean square error)
source("FuncPreComPiHierBayes.R")
library(doParallel)
cl <- makeCluster(2)
registerDoParallel(cl)

library(foreach)
library(optimx)

restmp <- trainPreComPiHier(CPI,sigma2j=2,sigma2dot=1,maxit =TRUE,iteration = 30000)
tmpget <- predictPreComPiHier(restmp,CPI)
cvtmp <- crossPreComPiHier(CPI,nfold=5,sigma2j=2,sigma2dot=1,maxit =TRUE,iteration = 30000,parallel=FALSE)

#saveRDS(tmpget,"tstsave")
#city <- readRDS("city.rds")

CVResultHier1 <- crossPreComPiHier(CPI,nfold=10,sigma2j=2,sigma2dot=1,maxit =TRUE,iteration = 50000,parallel=FALSE)
saveRDS(CVResultHier1,"CVResultHier1.rds")
#CVResultHier1 <- readRDS("CVResultHier1.rds")

CVResultHier2 <- crossPreComPiHier(CPI,nfold=10,sigma2j=0.1,sigma2dot=1,maxit =TRUE,iteration = 50000,parallel=FALSE)
saveRDS(CVResultHier2,"CVResultHier2.rds")


CVResultHier3 <- crossPreComPiHier(CPI,nfold=10,sigma2j=1,sigma2dot=0.1,maxit =TRUE,iteration = 50000,parallel=FALSE)

saveRDS(CVResultHier3,"CVResultHier3.rds")
save(CVResultHier3,file="saveCVResultHier3")

CVResultHier4 <- crossPreComPiHier(CPI,nfold=10,sigma2j=4,sigma2dot=1,maxit =TRUE,iteration = 50000,parallel=FALSE)
saveRDS(CVResultHier4,"CVResultHier4.rds")

CVResultHier5 <- crossPreComPiHier(CPI,nfold=10,sigma2j=5,sigma2dot=2,maxit =TRUE,iteration = 50000,parallel=FALSE)
saveRDS(CVResultHier5,"CVResultHier5.rds")

CVResultHier6 <- crossPreComPiHier(CPI,nfold=10,sigma2j=2,sigma2dot=5,maxit =TRUE,iteration = 50000,parallel=FALSE)
saveRDS(CVResultHier6,"CVResultHier6.rds")

#-- Compare with single ridge regression model.
# From the hierachical model, given one task or group, we sum the \omega_dot, 
# then \omega_j follows zero mean and \mu2_j+\mu2_dot variance.

library(glmnet)
plot(cv.glmnet(x = t(CPI$X[[6]]), y = CPI$y[[6]], family = "gaussian", alpha=0, nfold=10))
plot(cv.glmnet(x = t(CPI$X[[5]]), y = CPI$y[[5]], family = "gaussian", alpha=0, nfold=10))
plot(cv.glmnet(x = t(CPI$X[[4]]), y = CPI$y[[4]], family = "gaussian", alpha=0, nfold=10))
plot(cv.glmnet(x = t(CPI$X[[1]]), y = CPI$y[[1]], family = "gaussian", alpha=0, nfold=10))
plot(cv.glmnet(x = t(CPI$X[[2]]), y = CPI$y[[2]], family = "gaussian", alpha=0, nfold=10))
plot(cv.glmnet(x = t(CPI$X[[3]]), y = CPI$y[[3]], family = "gaussian", alpha=0, nfold=10))

cv.glmnet(x = t(CPI$X[[4]]), y = CPI$y[[4]], family = "gaussian", alpha=0, lambda = c(0.53,1),nfold=10)
cv.glmnet(x = t(CPI$X[[5]]), y = CPI$y[[5]], family = "gaussian", alpha=0, lambda = c(0.53,1),nfold=10)
# sigma2j = 2, sigma2dot = 1
sigma2j <- 2
sigma2dot <- 1
sigmay <- mean(sapply(c(1:10),function(x) CVResultHier1[[x]]$sigmay))
lamda <- sigmay/sum(sigma2j+sigma2dot)
p1 <- mean(sapply(c(1:10),function(x) CVResultHier1[[x]]$mse[1]))
p2 <- mean(sapply(c(1:10),function(x) CVResultHier1[[x]]$mse[2]))

mseHier <- function(CVresult){
  result <- matrix(,nrow=10,ncol=6)
  for (i in 1:length(CVresult)){
    result[i,] <- CVresult[[i]]$mse
  }
  return(result)
}
mse1 <- mseHier(CVResultHier1)
mse1revise <- mse1[c(1,3,5,7,8,9),]
colMeans(mse1revise)

performance <- c(1/32,1/12,4/25,6.4/11.6,1.7/3.4)

library(ggplot2)
df <- data.frame(Protein = factor(c("CHEMBL4015","CHEMBL3473","CHEMBL3974","CHEMBL4843")),
                                  CompoundNumber=c("1000","600","300","50"),
                 MSEreduce = c(1/32,4/25,1/12,1.7/3.4))
 
ggplot(data=df, aes(x=Protein, y=MSEreduce, fill=CompoundNumber,width=0.5)) + geom_bar(colour="black", stat="identity") +
  ggtitle("Hierarchial Bayesian Model v.s. Ridge Regression")+scale_fill_brewer() + 
  scale_x_discrete(name="Proteins from Peptide GRCRs")  + 
  theme(axis.title.x = element_text(size=20), axis.text.x  = element_text(vjust=0.5, size=16),
        axis.title.y = element_text(size=20), axis.text.y = element_text(vjust=0.5,size=16),
        plot.title = element_text(lineheight=.8, face="bold",size=20))
#        axis.title.y = element_text(face="bold",size=20), axis.text.y = element_text(vjust=0.5,size=12)

# Pattern
library(gplots)
#plotmatrix <- matrix(as.numeric(resultoptimx[1:1988]),ncol=p,byrow=TRUE)
resultoptimx <- trainPreComPiHier(CPI,sigma2j=2,sigma2dot=1,maxit =TRUE,iteration = 30000)
plotmatrix <- matrix(as.numeric(resultoptimx[1:((length(proteins)+1)*nrow(CPI$X[[1]]))]),nrow=p,byrow=FALSE)
distance <- dist(plotmatrix,method="euclidean")
fit <- hclust(distance, method="ward.D2") 
#redgreen(75)
#cm.colors(256)
#scale="column"
nba_heatmap <- heatmap.2(plotmatrix,scale="column",tracecol="darkblue",,Rowv=NA, Colv=NA, col = redblue(100),cexRow=0.3,dendrogram=c("row"))
nba_heatmap <- heatmap.2(plotmatrix,scale="column",tracecol="white",,Rowv=NA, Colv=NA, col = redgreen(75),cexRow=0.3)
