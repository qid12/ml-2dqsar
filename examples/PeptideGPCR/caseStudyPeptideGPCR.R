###Transfer Learning
###In this file we would use both independent linear regression and
###multiple task-baed linear regression.
###Songpeng Zu
###20141205

###set directory
setwd("D:/lab/TransferLearning/TransferModel/Code")
###load package
library(data.table)
library(DAAG)
library(glmnet)

###load data
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

###Delete some PubChem substructure's mannualy due to their chemical meaning.
###Unwanted:
###Sub1 to Sub115, they are all like  ">=H"
###Sub264 to Sub327, they are all like "Si-Cl"
subs.del <- c(paste("Sub",1:115,sep=""),paste("Sub",264:327,sep=""))
subs.select <-
setdiff(feafilter(data.frame(c(sub.dyn.raw[,2]),row.names=paste("Sub",1:881,sep="")),1),subs.del)

###Delete too low or too high frequency of substructures.
detect.dyn <- colSums(CPIs.all[,subs.select,with=FALSE])/nrow(CPIs.all)
hist(detect.dyn,breaks=40,col="red",
     main="Histogram of the frequency of detected substructures \nby dynamic slicing with equal partition method",
     xlab="Frequency of substructures in CPIs from Peptide GPCR family",warn.unused=TRUE,
     ylab="")
subs.del.freq <- c(subs.select[which(detect.dyn<=0.05)],subs.select[which(detect.dyn>=0.95)])
subs.select.feature <- setdiff(subs.select,subs.del.freq)

###filter some properties.

select.feature <- c(feafilter(data.frame(c(prop.dyn.raw[,2]),row.names=prop.dyn.raw[,1]),1),subs.select.feature)

CPIs.select.colname <- c("protein","compound","affinity","logP","MR","MW","TPSA",
                         select.feature)
CPIs.select <- CPIs.all[,CPIs.select.colname,with=FALSE]

###??? Unique the compounds with high chemical similarity.
###prot.stat <- as.data.frame(table(as.factor(CPIs.select$protein)))
###prot.stat.sort <- prot.stat[order(prot.stat$Freq,decreasing=TRUE),]

###temp case for unique the compounds.
###specprotein <- "CHEMBL344"
###setkey(CPIs.select,"protein")
###speccompounds <- as.data.frame(CPIs.select[specprotein,select.feature,with=FALSE])
###row.names(speccompounds) <- unlist(CPIs.select[specprotein,"compound",with=FALSE])
###comsim <- dist(speccompounds,method="binary")
###fit <- hclust(comsim,method="ward.D")

###pdf("temp.pdf",width=80,height=15)
###plot(fit,labels=FALSE)
###dev.off()
###??? scale the covariates.

###　Independent linear regressio
proteins <- unique(as.factor(CPIs.select$protein))
prot.stat <- as.data.frame(table(as.factor(CPIs.select$protein)))
prot.stat.sort <- prot.stat[order(prot.stat$Freq,decreasing=TRUE),]
setkey(CPIs.select,"protein")

specpro1 <- as.character(prot.stat.sort$Var1[1])
speccom1 <- as.data.frame(CPIs.select[specpro1,])
data1.reg <- data.frame(y=log(speccom1$affinity),speccom1[,colnames(speccom1)[4:ncol(speccom1)]])
###fmla <- as.formula(paste("log(affinity) ~ ",paste(colnames(speccom1)[4:ncol(speccom1)],collapse="+")))
###mod1 <- lm(formula = fmla,data = speccom1)
mod1 <- lm(y~.,data=data1.reg)
cv.lm(df=data1.reg, mod1, m=5) 

###Regression Diagnostics.
###library(car)
###outlierTest(mod1)
###qqPlot(mod1, main="QQ Plot")
###collinear <- vif(mod1)

###Try ridge-regression

x <- as.matrix(data1.reg[,2:ncol(data1.reg)])
y <- as.matrix(data1.reg[,1])
set.seed(1)
train <- sample(1:nrow(x), nrow(x)*2/3)
test <- (-train)
r1 <- glmnet(x = x[train, ], y = y[train], family = "gaussian", alpha = 0.02)
plot(r1, xvar = "lambda")
r1.cv <- cv.glmnet(x = x, y = y, family = "gaussian", alpha = 0.02, nfold = 10)
###r1.cv <- cv.glmnet(x = x[train,], y = y[train], family = "gaussian", alpha = 0, nfold = 10)
plot(r1.cv)

mte <- predict(r1, x[test, ])
mte <- apply((mte - y[test])^2, 2, mean)
points(log(r1$lambda), mte, col = "blue", pch = 19)
legend("topleft", legend = c("10 - fold CV", "Test"), col = c("red", "blue"))

###Try Elastic-net regression model.
alpha.cv.glm <- function(alpha,x,y,train.p=2/3,nfold=10){
  
  plot(cv.glmnet(x = x, y = y, family = "gaussian", alpha=alpha, nfold=nfold),main=paste("Alpha is ",alpha))
  
  train <- sample(1:nrow(x), nrow(x)*train.p)
  test <- (-train)
  r1 <- glmnet(x = x[train, ], y = y[train], family = "gaussian", alpha=alpha)
  mte <- predict(r1, x[test, ])
  mte <- apply((mte - y[test])^2, 2, mean)
  points(log(r1$lambda), mte, col = "blue", pch = 19)
  legend("topleft", legend = c("10 fold CV red", "Test blue"), col = c("red", "blue"))
}

alpha.array <- seq(from=1,to=0.51,by=-0.01)
pdf("elasticnet_alpha_1.pdf",width=20,height=50) par(mfrow=c(10,5))
lapply(alpha.array,function(t)
alpha.cv.glm(t,as.matrix(data1.reg[,2:ncol(data1.reg)]),as.matrix(data1.reg[,1])))
dev.off()

alpha.array <- seq(from=0.5,to=0.01,by=-0.01)
pdf("elasticnet_alpha_2.pdf",width=20,height=50) par(mfrow=c(10,5))
lapply(alpha.array,function(t)
alpha.cv.glm(t,as.matrix(data1.reg[,2:ncol(data1.reg)]),as.matrix(data1.reg[,1])))
dev.off()
