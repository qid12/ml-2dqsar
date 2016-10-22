# QSAR on ginkgolide
# Songpeng Zu
# 2014-10-23

setwd("D:/lab/ginkgolide")
library(ChemmineR)
library(fmcsR)
library(ChemmineOB)

# Read the profile of CHEMBL250 (P25105).
PAFRComs <- read.table("Chembl2Cid2Affinity_250.txt",sep="\t",quote="",as.is= FALSE)
PAFRComs.matrix <- data.frame(PAFRComs[,2:3])
rownames(PAFRComs.matrix) <- PAFRComs[,1]
colnames(PAFRComs.matrix) <- c("CID","Affinity")
PAFRComs.df <- data.frame(PAFRComs.matrix)


PAFRComs <- getIds(PAFRComs.matrix[,1])
GIN.compounds <- getIds(c(9909368,6324617,16211418,441296,46937025))

# Get the compounds' properties From ChemmineR.
prop.PAFRComs <- propOB(PAFRComs)
propma.PAFRComs <- data.frame(MF=MF(PAFRComs, addH=TRUE), MW=MW(PAFRComs, addH=TRUE), 
                     Ncharges=sapply(bonds(PAFRComs, type="charge"), length), 
                     atomcountMA(PAFRComs, addH=TRUE), 
                     groups(PAFRComs, type="countMA"), 
                     rings(PAFRComs, upper=6, type="count", arom=TRUE))
prop.all.PAFRComs <- data.frame(propma.PAFRComs,prop.PAFRComs[,c("HBA1","HBA2","HBD","logP","MR","nF","TPSA")])
rownames(prop.all.PAFRComs) <- rownames(PAFRComs.matrix)

prop.GIN <- propOB(GIN.compounds)
propma.GIN <- data.frame(MF=MF(GIN.compounds, addH=TRUE), MW=MW(GIN.compounds, addH=TRUE), 
                              Ncharges=sapply(bonds(GIN.compounds, type="charge"), length), 
                              atomcountMA(GIN.compounds, addH=TRUE), 
                              groups(GIN.compounds, type="countMA"), 
                              rings(GIN.compounds, upper=6, type="count", arom=TRUE))
prop.all.GIN <- data.frame(propma.GIN,prop.GIN[,c("HBA1","HBA2","HBD","logP","MR","nF","TPSA")])
rownames(prop.all.GIN) <- c("CHEMBL465161", "CHEMBL514432","CHEMBL221821","CHEMBL374004","CHEMBL2113267")


# Select the covariates according to the differences properties of GINs.
select.cov <- c("MW","O","ROH","HBA1","HBA2","HBD","logP","MR","TPSA")

prop.chemmine.GIN <- prop.all.GIN[,select.cov]
prop.chemmine.PAFRComs <- prop.all.PAFRComs[,select.cov]
prop.chemmine.PAFRComs <- prop.chemmine.PAFRComs[-which(rownames(prop.chemmine.PAFRComs)=="CHEMBL514432"),]
prop.chemmine.all <- rbind(prop.chemmine.GIN,prop.chemmine.PAFRComs)

simscale <- function(x){
  (x-min(x))/(max(x)-min(x))
}

prop.chemmine.all.scale <- apply(prop.chemmine.all,2,simscale)

#CHEMBL514432 should be deleted first.

# Regression model.

reg.comps.tmp <- data.frame(prop.chemmine.PAFRComs,affinity=PAFRComs.matrix[,2])
reg.comps <- reg.comps[-which(rownames(reg.comps.tmp)=="CHEMBL514432"),]

reg.comps.scale <- data.frame(prop.chemmine.all.scale[rownames(prop.chemmine.PAFRComs),],
                              affinity=reg.comps[,"affinity"])
fit.scale <- lm(-log10(affinity) ~ MW + O + ROH+HBA1+HBA2+HBD+logP+MR+TPSA, data=reg.comps.scale)
summary(fit.scale)
predict(fit.scale, data.frame(prop.chemmine.all.scale[rownames(prop.all.GIN),]))

fit <- lm(-log10(affinity) ~ MW + O + ROH+HBA1+HBA2+HBD+logP+MR+TPSA, data=reg.comps)
summary(fit)
predict(fit,prop.chemmine.GIN)




layout(matrix(c(1,2,3,4),2,2)) # optional 4 graphs/page 
plot(fit.scale)



library(lattice)
library(DAAG)
cv.lm(df=reg.comps, fit.scale, m=5) # 5 fold cross-validation

library(MASS)
step <- stepAIC(fit, direction="both")
step$anova # display results

fit1 <- lm(-log(affinity) ~ MW + O + ROH+HBA1+HBA2+HBD+logP+MR+TPSA, data=reg.comps)
fit2 <- lm(-log(affinity) ~ MW+ O + ROH+HBA1+HBD, data=reg.comps)
anova(fit1, fit2)

library(relaimpo)
calc.relimp(fit,type=c("lmg","last","first","pratt"),
            rela=TRUE)

boot <- boot.relimp(fit, b = 1000, type = c("lmg", 
                                            "last", "first", "pratt"), rank = TRUE, 
                    diff = TRUE, rela = TRUE)
booteval.relimp(boot) # print result
plot(booteval.relimp(boot,sort=TRUE)) # plot result


vif(fit)
sqrt(vif(fit))

library(car)
crPlots(fit)
ceresPlots(fit)
