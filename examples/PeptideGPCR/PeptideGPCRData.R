# Generate one group data for transferlearning.
# PeptideGPCR 
# Songpeng Zu 
# 2014-12-01

setwd("D:/lab/TransferLearning/TransferModel/Code")
#install.packages("dynslicing_1.0.zip",repos = NULL,type="source")
library(ggplot2)
library(data.table)
library(dynslicing)
#----------------------------------------------------------------------------------#
# Information of Protein Classification from Homo sapiens
Pro2Class <- read.table("Protein2Class_Chembl.txt",sep="\t",header=FALSE,quote="",
                        col.names=c("protein","class"))

ProClassName <- read.table("ProClassName.txt",sep="\t",header=FALSE,quote="",colClasses = "character",
                           col.names=c("class","detail"))

OrganPro <- read.table("organism.txt",sep="\t",header=FALSE,quote="",colClasses=c("character","factor"),
                       col.names=c("protein","organism"))

total <- merge(Pro2Class,OrganPro,by="protein")
total <- unique(total)
# We noticed that four proteins have multiple classes from the source file. We modified them and change 
# the data of "Protein2Class_chembl.txt".

# The four proteins: CHEMBL2683 CHEMBL3521 CHEMBL3646 CHEMBL5213
# The modified results:
# CHEMBL2683 Transporter enzyme Both are right. (Do not need to change.)
# CHEMBL3521 phosphatase
# CHEMBL3646 phosphatase
# CHEMBL5213 protein kinase

pro.human <- total[which(total$organism=="Homo sapiens"),]
#pie(table(pro.human$class))
#barplot(table(pro.human$class),cex.names = 0.5,las=2)
#barplot(table(pro.human$class),cex.names = 0.5,horiz=TRUE,las=2,col="red")
proclass.stat <- as.data.frame(table(pro.human$class))
colnames(proclass.stat) <- c("ProteinClass","Frequency")
bp <- ggplot(proclass.stat, aes(ProteinClass,Frequency)) + geom_bar(fill="red",stat="identity") 
bp +theme(axis.text.y = element_text(colour = "black",size=16),
          axis.title.y = element_blank(),axis.title.x=element_text(colour="black",
                                                                   size=16)) +coord_flip()
#--------------------------------------------------------------------------------------#
# Information of Known Protein-Compound Interactions from Chembl
CPIs <- read.table("compounds_targets_interaction_7.txt",sep="\t",header=FALSE,quote="",
                   colClasses=c("character","character","numeric","NULL"),
                   col.names=c("compound","protein","affinity","NA"))
CPIs.class <- merge(CPIs,pro.human,by="protein")
write.table(CPIs.class,file="CPIsByR.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)

rm(list=ls())
#--------------------------------------------------------------------------------------#
# We firstly see GPCRs.
CPIs.class <- read.table("CPIsByR.txt",quote="",sep="\t",header=TRUE)
CPIs.PeptideGPCR <- CPIs.class[which(CPIs.class$class=="Peptide GPCR"),1:3]
CPIs.otherGPCR <- CPIs.class[which(CPIs.class$class=="GPCR others"),1:3]
CPIs.monoGPCR <- CPIs.class[which(CPIs.class$class=="GPCR monoamine receptor"),1:3]
CPIs.lipidGPCR <- CPIs.class[which(CPIs.class$class=="GPCR  lipid-like ligand receptor"),1:3]
CPIs.nucGPCR <- CPIs.class[which(CPIs.class$class=="GPCR nucleotide-like receptor"),1:3]

# We next mainly focus on Peptide GPCRs.
# Need to only consider the unique affinities.
rm(list=c("CPIs.otherGPCR","CPIs.monoGPCR","CPIs.lipidGPCR","CPIs.nucGPCR","CPIs.class"))
z <- as.data.table(CPIs.PeptideGPCR)
keys <- colnames(z)[!grepl('affinity',colnames(z))]
CPIs.PeptideGPCR.ave <- z[,list(aveaffinty=mean(affinity)),keys]
rm(list=c("z","keys","CPIs.PeptideGPCR"))

CPIs.PeptideGPCR.stat <- as.data.frame(table(CPIs.PeptideGPCR.ave$protein))
colnames(CPIs.PeptideGPCR.stat) <- c("PeptideGPCR","Frequency")
bp <- ggplot(CPIs.PeptideGPCR.stat, aes(PeptideGPCR,Frequency)) + geom_bar(fill="blue",stat="identity") 
bp +theme(axis.text.y = element_text(size=8),
          axis.title.y = element_text(colour="black",face="bold",size=16),
          axis.title.x=element_text(colour="black",size=16)) +coord_flip()
rm("CPIs.PeptideGPCR.stat")
#-------------------------------------------------------------------------------------#
# Load Compounds Substructures (from PubChem) and other properties.
# Substructure.
substruc <- fread("SubStructure.txt",sep="\t",header=FALSE)
cid.subs <- read.table("cid.txt",header=FALSE,colClasses="character",col.names="cid")
#rownames(substruc) <- cid.subs$cid
substruc <- data.table(cbind(cid.subs,substruc))
setkey(substruc,cid)

chem2cid <- read.table("chembl_id-cid.txt",header=FALSE,colClasses = "character",
                       col.names=c("compound","cid"),sep="\t",quote="")
# Multi chembl IDs versus one cid
chem2cid <- chem2cid[!duplicated(chem2cid$cid),]
CPIs.sub.pep.cid <- merge(CPIs.PeptideGPCR.ave,chem2cid,by="compound")
substruc.pep <- substruc[CPIs.sub.pep.cid$cid,]
CPIs.sub.pep.cid <- cbind(CPIs.sub.pep.cid,substruc.pep[,2:ncol(substruc.pep),with=FALSE])
colnames(CPIs.sub.pep.cid) <- c("compound","protein","affinity","cid",paste("Sub",1:881,sep=""))
rm(list=c("substruc","cid.subs","chem2cid"))
# Try dynamic slicing method to filter some substructures.
subs.feature.total <- substruc.pep[order(CPIs.sub.pep.cid$affinity),2:ncol(substruc.pep),with=FALSE]
dim <- 2
lambda <- 1.0
alpha <- 1.0

bf.subs <- apply(subs.feature.total,2,function(x) ds_eqp_k(x,dim,lambda))
rm(list=c("substruc.pep","CPIs.PeptideGPCR"))
#Bayes factor is too slow in our sample set, so we use the equal partition. 
#bf.subs.bf <- apply(subs.feature.total,2,function(x) ds_bf_u(x,dim,lambda,alpha))


#-----------------------------------------------------------------------------------#
# Load Compounds physico-chemical features generated by ChemmineR.
propers <- fread("compoundcharacter/prop.csv",sep="\t",verbose=TRUE,colClasses = "character")
colnames(propers)[1] <- "cid"
setkey(propers,cid)
chem2cid.prop <- read.table("compoundcharacter/chembl_id-cid.txt",header=FALSE,
                            colClasses="character",sep="\t",quote="",col.names=c("compound","cid"))
chem2cid.prop <- chem2cid.prop[!duplicated(chem2cid.prop$cid),]

propers.chem <- cbind(chem2cid.prop$compound, propers[chem2cid.prop$cid,])
colnames(propers.chem)[1] <- "compound"

select.prop <- c("compound","cid","cansmi","HBA1","HBA2","HBD","logP","MR","MW","nF","TPSA","RNH2","R2NH",
                 "R3N","ROPO3","ROH","RCHO","RCOR","RCOOH","RCOOR","ROR","RCCH","RCN","RINGS","AROMATIC")
propers.chem.select <- as.data.frame(propers.chem)[,select.prop]

rm(list=c("propers","chem2cid.prop","propers.chem","select.prop"))

CPIs.prop.pep.cid <- merge(CPIs.PeptideGPCR.ave,propers.chem.select,by="compound")
rm("propers.chem.select")

colnames(CPIs.prop.pep.cid)[3] <- "affinity"

# Try bayes factor.
col.prop.feature <- c("HBA1","HBA2","HBD","nF","RNH2","R2NH",
                      "R3N","ROPO3","ROH","RCHO","RCOR","RCOOH","RCOOR","ROR","RCCH","RCN","RINGS","AROMATIC")
prop.feature.total <- CPIs.prop.pep.cid[order(CPIs.prop.pep.cid$affinity),
                                      c(6,7,8,12,14:ncol(CPIs.prop.pep.cid)),with=FALSE]
lambda <- 1.0
bf.prop <- apply(as.data.frame(prop.feature.total),2,function(x) 
  ds_eqp_k(as.integer(x),max(as.integer(x)+1),lambda))


# Export the CPIs data and Results of Dynamic slicing.
write.table(CPIs.prop.pep.cid,file="CPIs_property_PeptideGPCR",quote=FALSE,sep="\t",
            col.names=TRUE,row.names=FALSE)
write.table(CPIs.sub.pep.cid,file="CPIs_substructure_PeptideGPCR",quote=FALSE,sep="\t",
            col.names=TRUE,row.names=FALSE)
write.table(bf.subs,file="FeatureSelect_Substructure_PepetideGPCR",quote=FALSE,sep="\t",
            col.names=TRUE,row.names=TRUE)
write.table(bf.prop,file="FeatureSelect_Propter_PepetideGPCR",quote=FALSE,sep="\t",
            col.names=TRUE,row.names=TRUE)

