# Extract the compounds' chemical features .
# Songpeng Zu /zusongpeng@gmail.com/

#-- load package
library(data.table)
library(ChemmineR)
library(Rcpi)

#-- load files.
phychemprop <- fread("cidProperty.txt",sep="\t",verbose=TRUE)
colnames(phychemprop)[1] <- "cid"

cids <- read.table("cid.txt",col.names = c("cid"))
SMIcidstripsalts <- read.table("StripSalts.smiles",col.names=c("smilesStripSalts"),colClasses="character")
ncharSmiles <- sapply(SMIcidstripsalts$smilesStripSalts,function(x) nchar(x))
cid2smi <- cbind(cids,SMIcidstripsalts,ncharSmiles)

# Delete duplicates the smiles.
cid2smitmp <- cid2smi[!duplicated(cid2smi[,c("smilesStripSalts")]),]
# Delete "small" molecules (small: the number of atoms is no more than three without H.)
cid2smiuniq <- cid2smitmp[which(cid2smitmp$ncharSmiles>=3),]

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

# I face some problems in this part. And try to use the python packages of RDKit and openbabel.
chem2cid <- read.table("chembl_id-cid.txt",header=FALSE,col.names=c("compound","cid"),sep="\t",quote="")
chem2cid2smiles <- merge(chem2cid,chemblphychemp[,c("cid","smilesStripSalts","cansmi","cansmiNS"),with=FALSE],
                            by=c("cid"))
write.table(chem2cid2smiles,file="chemble2cid2smiles.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
