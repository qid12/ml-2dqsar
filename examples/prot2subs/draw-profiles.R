### draw subprofiles in a given protein family.
### szu
### 2016-10-27

## load package
library(corrplot)
library(ape)

## set 
setwd("~/home/songpeng/git-recipes/QSARMulT/examples/prot2subs")
col3 <- colorRampPalette(c("blue", "white", "red")) 
samplenum <- 50
### load data
load('m_monoamine_fingerprint.data')

### draw profles
drawcorplot <- function(fn,samplenum){
  load(fn)
  choose_index = sample(ncol(mresult$omega_matrix), samplenum)
  corrplot(mresult$omega_matrix[ ,choose_index], is.corr=FALSE, diag = FALSE,
                        method = "circle", tl.pos="n",
                        col=col3(20))
  return(choose_index)
}

getpubchemind <- function(mresult, filterind){
  return(mresult$kf[filterind]-1)
}

### run
k_mono = "m_monoamine_fingerprint.data"
g_lipid = "m_lipidlike_fingerprint.data"
g_s_lipid = "s_lipidlike_fingerprint.data"
samind = drawcorplot(g_lipid, samplenum)
gpubchemind = getpubchemind(mresult, samind)
load(g_lipid)
load(g_s_lipid)
View(mresult$omega_matrix[,samind])
corrplot(sresult$omega_matrix[ ,samind], is.corr=FALSE, diag = FALSE,
         method = "circle", tl.pos="n",
         col=col3(20))

### ginkgolide structure similarity
simgink = as.matrix(read.table("ginkgolidesim"))
row.names(simgink) = c("GA","GB","GC","GJ","GM","GK")
colnames(simgink) =  c("GA","GB","GC","GJ","GM","GK")
corrplot(simgink, is.corr = FALSE,
         method = "circle",  type="lower", tl.pos="d",tl.cex=2.5,cl.cex = 2)
### phylo tree.
mypal = c("#556270", "#4ECDC4", "#1B676B", "#FF6B6B", "#C44D58")
hc = hclust(dist(mresult$omega_matrix))
cutnum <- 5 
clus5 <- cutree(hc,cutnum)
plot(as.phylo(hc),type="unrooted",col="gray80",
     use.edge.length = TRUE,
     tip.color = mypal[clus5],
     label.offset = 1)

