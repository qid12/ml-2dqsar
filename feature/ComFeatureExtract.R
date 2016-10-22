# Extract the compounds' chemical features .
# Songpeng Zu /zusongpeng@gmail.com/

#-- load package
library(data.table)
library(ChemmineR)
#library(Rcpi)

#-- load files.
phychemprop <- fread("cidProperty.txt",sep="\t",verbose=TRUE)
colnames(phychemprop)[1] <- "cid"

cids <- read.table("cid.txt",col.names = c("cid"))
# Note: Since read.table default : comment.char = "#", and "#" has a specific meaning in the SMILES format, we must use
# comment.char="" when reading SMILES from the disc using the read.table function.
# While in fread, you don't need to worry about this.
SMIcidstripsalts <- read.table("StripSalts.smiles",col.names=c("smilesStripSalts"),colClasses="character",comment.char = "")
ncharSmiles <- sapply(SMIcidstripsalts$smilesStripSalts,function(x) nchar(x))
cid2smi <- cbind(cids,SMIcidstripsalts,ncharSmiles)

# Delete duplicates the smiles.
cid2smitmp <- cid2smi[!duplicated(cid2smi[,c("smilesStripSalts")]),]
# Delete "small" molecules (small: the number of atoms is no more than three without H.)
cid2smiuniq <- cid2smitmp[which(cid2smitmp$ncharSmiles>=4),]

#-- Physicochemical properties
chemblphychemp <- merge(phychemprop,cid2smiuniq,by=c("cid"))
# 
# MWStripSalts <- MW(smiles2sdf(chemblphychemp$smilesStripSalts))
# 
# MWStripSalts <- rep(0,times=nrow(chemblphychemp))
# for (i in 1:nrow(chemblphychemp)){
#   MWStripSalts[i] <- MW(smiles2sdf(chemblphychemp$cansmi[i]))
# }
# 
# 
# sdfgetstripsalts<- smiles2sdf(chemblphychemp$smilesStripSalts[1:10])
# Bugs
chemblbugs <- chemblphychemp[which(chemblphychemp[,c("")])]

# I face some problems in this part. And try to use the python packages of RDKit and openbabel.
chem2cid <- read.table("chembl_id-cid.txt",header=FALSE,col.names=c("compound","cid"),sep="\t",quote="")
chem2cid2smilestmp <- merge(chem2cid,chemblphychemp[,c("cid","smilesStripSalts","cansmi","cansmiNS"),with=FALSE],
                            by=c("cid"))
chem2cid2smiles <- chem2cid2smilestmp[!duplicated(chem2cid2smilestmp[,c("cid")]),]
write.table(chem2cid2smiles,file="chemble2cid2smiles.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
