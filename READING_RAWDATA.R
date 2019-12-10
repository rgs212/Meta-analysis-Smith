
rm(list=ls())
library(minfi)
library(parallel)
library(wateRmelon)
######################################################################
setwd("/mnt/data1/Ehsan/meta_CrossCortical")
######################################################################
pheno<-read.csv("/mnt/data1/Bex/Meta/SVAregressed/Allphenosub.csv")
pheno.AZ<-pheno[which(pheno$Cohort=="Az" & pheno$region=="STG"),]
rm(pheno)
######################################################################
idatPath <- c("/mnt/data1/Ehsan/meta_CrossCortical/idats")
barcodes<-as.character(pheno.AZ$sample)
AZ<- readEPIC(idatPath=idatPath, barcodes=barcodes, parallel = TRUE)
identical(colnames(AZ), as.character(pheno.AZ$sample))
######################################################################
load("tissue.MTG.rda")
rownames(tissue.MTG)<-as.character(tissue.MTG$CompleteBarcode)
pheno.MTG<-tissue.MTG[barcodes,]
identical(colnames(AZ), rownames(pheno.MTG))
pheno<-rbind(pheno.AZ,pheno.MTG)
######################################################################
save(AZ, pheno, file="AZ_methylumiSet_pheno.Rdata")
######################################################################
######################################################################
setwd("/mnt/data1/Ehsan/meta_CrossCortical")

targets <- data.frame(barcodes)
idatPath <- c("/mnt/data1/Ehsan/meta_CrossCortical/idats")
colnames(targets) <- "sentrixFull"
targets$Basename <- with(targets, file.path(idatPath, sentrixFull))

# Construct RGset from IDAT files
rgSet <- read.metharray.exp(targets = targets, verbose = TRUE, extended = FALSE)
save(rgSet, pheno, file="rgSet404.AZ_pheno_cellProp.Rdata")
######################################################################
######################################################################

setwd("/mnt/data1/Ehsan/MTG/")
pheno.az2<-read.csv("pheno192.csv")
pheno.az2<-pheno.az2[which(pheno.az2$Sample_Group=="BS"),]

idatPath <- c("/mnt/data1/Ehsan/MTG/IDATS")
barcodes<-as.character(pheno.az2$Sentrix)
AZ2<- readEPIC(idatPath=idatPath, barcodes=barcodes, parallel = TRUE)
save(AZ2, pheno.az2, file="AZ2_methylumiSet_pheno.Rdata")
######################################################################
setwd("/mnt/data1/Ehsan/meta_CrossCortical")

targets <- data.frame(barcodes)
colnames(targets) <- "sentrixFull"
targets$Basename <- with(targets, file.path(idatPath, sentrixFull))

# Construct RGset from IDAT files
rgSet <- read.metharray.exp(targets = targets, verbose = TRUE, extended = FALSE)
save(rgSet, pheno.az2, file="rgSet96.AZ2_pheno_cellProp.Rdata")
######################################################################
######################################################################

betas<-betas(AZ)
betas.rs.az1<-betas[grep("rs", rownames(betas)),]

betas<-betas(AZ2)
betas.rs.az2<-betas[grep("rs", rownames(betas)),]

betas.rs.az12<-cbind(betas.rs.az1,betas.rs.az2)

snpCor<-cor(betas.rs.az12, use="complete.obs")

for(i in 1:ncol(betas.rs.az12)){
  snpCor[i,i]<-NA
}


#finds the maximum correlation values between samples and plots a histogram
#note: correlations of ~0.5 do not mean siblings as there are too few SNPs for this
corMax<-apply(snpCor, 1, max, na.rm = TRUE)
pdf("SNPCorrelations.LBB1.LC.pdf")
hist(corMax, xlab = "Max. correlation with all other samples", main = "")
dev.off()
genetic.cor<-data.frame(corMax)
overlap<-rownames(genetic.cor)[which(genetic.cor[,1] > 0.9)]
az1<-intersect(colnames(AZ),overlap)
az2<-intersect(colnames(AZ2),overlap)

AZ1<-AZ[,-which(colnames(AZ)%in%az1)]
dim(AZ1)
rownames(pheno)<-as.character(pheno$X)
pheno.AZ1<-pheno[colnames(AZ1),]
identical(colnames(AZ1),as.character(pheno.AZ1$X))
#[1] TRUE

AZ2<-AZ2[,which(colnames(AZ2)%in%az2)]
rownames(pheno.az2)<-as.character(pheno.az2$Sentrix)
pheno.AZ2<-pheno.az2[colnames(AZ2),]
identical(colnames(AZ2),as.character(pheno.AZ2$Sentrix))
#[1] TRUE

load("rgSet404.AZ_pheno_cellProp.Rdata")
rgsetAZ1<-rgSet[, colnames(AZ1)]
identical(colnames(rgsetAZ1),as.character(pheno.AZ1$X))
#[1] TRUE

load("rgSet96.AZ2_pheno_cellProp.Rdata")
rgsetAZ2<-rgSet[, colnames(AZ2)]
identical(colnames(rgsetAZ2),as.character(pheno.AZ2$Sentrix))
#[1] TRUE

setwd("/mnt/data1/Ehsan/meta_CrossCortical/RAW_DATASETS")
save(AZ1,rgsetAZ1,pheno.AZ1,file="AZ1_RgSet_Methylumi_Pheno_309.Rdata")
save(AZ2,rgsetAZ2,pheno.AZ2,file="AZ2_RgSet_Methylumi_Pheno_95.Rdata")
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
rm(list=ls())
###reading the demographics and phenotipic aoutcomes
setwd("/mnt/data1/Ehsan/meta_CrossCortical/LBB1/PFC")
pheno<-read.csv("/mnt/data1/Bex/Meta/SVAregressed/Allphenosub.csv")
###################################
pheno.LBB1.PFC<-pheno[which(pheno$Cohort=="LBB1" & pheno$region=="PFC"),]
idatPath <- c("/mnt/data1/Bex/AD/iDATs/LBBA/")
barcodes<-as.character(pheno.LBB1.PFC$X)
LBB1.PFC<- readEPIC(idatPath=idatPath, barcodes=barcodes, parallel = TRUE)
dim(LBB1.PFC)#485577      114
identical(colnames(LBB1.PFC), as.character(pheno.LBB1.PFC$X))
#[1] TRUE
###################################
pheno.LBB1.EC<-pheno[which(pheno$Cohort=="LBB1" & pheno$region=="EC"),]
idatPath <- c("/mnt/data1/Bex/AD/iDATs/LBBE/")
barcodes<-as.character(pheno.LBB1.EC$X)
LBB1.EC<- readEPIC(idatPath=idatPath, barcodes=barcodes, parallel = TRUE)
dim(LBB1.EC)#485577      107
identical(colnames(LBB1.EC), as.character(pheno.LBB1.EC$X))
#[1] TRUE
###################################
pheno.LBB1.STG<-pheno[which(pheno$Cohort=="LBB1" & pheno$region=="STG"),]
idatPath <- c("/mnt/data1/Bex/AD/iDATs/LBBF/")
barcodes<-as.character(pheno.LBB1.STG$X)
LBB1.STG<- readEPIC(idatPath=idatPath, barcodes=barcodes, parallel = TRUE)
dim(LBB1.STG)#485577      117
identical(colnames(LBB1.STG), as.character(pheno.LBB1.STG$X))
#[1] TRUE
############################################################################
barcodes<-as.character(pheno.LBB1.PFC$X)
targets <- data.frame(barcodes)
idatPath <- c("/mnt/data1/Bex/AD/iDATs/LBBA/")
colnames(targets) <- "sentrixFull"
targets$Basename <- with(targets, file.path(idatPath, sentrixFull))

# Construct RGset from IDAT files
rgSet.PFC <- read.metharray.exp(targets = targets, verbose = TRUE, extended = FALSE)
############################################################################
############################################################################
barcodes<-as.character(pheno.LBB1.EC$X)
targets <- data.frame(barcodes)
idatPath <- c("/mnt/data1/Bex/AD/iDATs/LBBE/")
colnames(targets) <- "sentrixFull"
targets$Basename <- with(targets, file.path(idatPath, sentrixFull))

# Construct RGset from IDAT files
rgSet.EC <- read.metharray.exp(targets = targets, verbose = TRUE, extended = FALSE)
############################################################################
barcodes<-as.character(pheno.LBB1.STG$X)
targets <- data.frame(barcodes)
idatPath <- c("/mnt/data1/Bex/AD/iDATs/LBBF/")
colnames(targets) <- "sentrixFull"
targets$Basename <- with(targets, file.path(idatPath, sentrixFull))

# Construct RGset from IDAT files
rgSet.STG <- read.metharray.exp(targets = targets, verbose = TRUE, extended = FALSE)
############################################################################
setwd("/mnt/data1/Ehsan/meta_CrossCortical/RAW_DATASETS")
save(LBB1.PFC, LBB1.EC, LBB1.STG, 
     rgSet.PFC, rgSet.EC, rgSet.STG, 
     pheno.LBB1.PFC, pheno.LBB1.EC, pheno.LBB1.STG, 
     file="LBB1__RgSet_Methylumi_Pheno_EC107_PFC114_STG117.Rdata")
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
pheno<-read.csv("/mnt/data1/Bex/Meta/SVAregressed/Allphenosub.csv")
pheno.LBB2.EC<-pheno[which(pheno$Cohort=="LBB2" & pheno$region=="EC"),]
idatPath <- c("/mnt/data1/Ehsan/meta_CrossCortical/LBB2/IDATS")
barcodes<-as.character(pheno.LBB2.EC$X)
LBB2.EC<- readEPIC(idatPath=idatPath, barcodes=barcodes, parallel = TRUE)
dim(LBB2.EC)#485577      92
identical(colnames(LBB2.EC), as.character(pheno.LBB2.EC$X))
#[1] TRUE
###################################
barcodes<-as.character(pheno.LBB2.EC$X)
targets <- data.frame(barcodes)
idatPath <- c("/mnt/data1/Ehsan/meta_CrossCortical/LBB2/IDATS")
colnames(targets) <- "sentrixFull"
targets$Basename <- with(targets, file.path(idatPath, sentrixFull))

# Construct RGset from IDAT files
rgSet.EC.2 <- read.metharray.exp(targets = targets, verbose = TRUE, extended = FALSE)
setwd("/mnt/data1/Ehsan/meta_CrossCortical/RAW_DATASETS")
save(LBB2.EC,rgSet.EC.2, pheno.LBB2.EC, file="LBB2__RgSet_Methylumi_Pheno_EC92.Rdata")
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
###reading the demographics and phenotipic aoutcomes
setwd("/mnt/data1/Ehsan/meta_CrossCortical/MS/PFC")
pheno<-read.csv("/mnt/data1/Bex/Meta/SVAregressed/Allphenosub.csv")
###################################
pheno.MS<-pheno[which(pheno$Cohort=="MS"),]
idatPath <- c("/mnt/data1/Bex/Mount Sinai/iDats")
barcodes<-as.character(pheno.MS$X)
MS<- readEPIC(idatPath=idatPath, barcodes=barcodes, parallel = TRUE)
dim(MS)#
identical(colnames(MS), as.character(pheno.MS$X))
#[1] TRUE
targets <- data.frame(barcodes)
colnames(targets) <- "sentrixFull"
targets$Basename <- with(targets, file.path(idatPath, sentrixFull))
# Construct RGset from IDAT files
rgSet.MS <- read.metharray.exp(targets = targets, verbose = TRUE, extended = FALSE)
save(MS,rgSet.MS,pheno.MS, file="MS.methylumi.rgSet.pheno.Rdata")
############################################################################
rgSet.MS.PFC<-rgSet.MS[, which(pheno.MS$region=="PFC")]
rgSet.MS.STG<-rgSet.MS[, which(pheno.MS$region=="STG")]


MS.PFC<-MS[, which(pheno.MS$region=="PFC")]
MS.STG<-MS[, which(pheno.MS$region=="STG")]

pheno.MS.PFC<-pheno.MS[which(pheno.MS$region=="PFC"),]
pheno.MS.STG<-pheno.MS[which(pheno.MS$region=="STG"),]
setwd("/mnt/data1/Ehsan/meta_CrossCortical/RAW_DATASETS")
save(MS.PFC,MS.STG,
     rgSet.MS.PFC,rgSet.MS.STG,
     pheno.MS.PFC,pheno.MS.STG,
     file="MS__RgSet_Methylumi_Pheno_PFC142_STG144.Rdata")
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
pheno<-read.csv("/mnt/data1/Bex/Meta/SVAregressed/Allphenosub.csv")
pheno.RM<-pheno[which(pheno$Cohort=="ROSMAP" & pheno$region=="PFC"),]
rm(pheno)
###
idatPath <- c("/mnt/data1/Ehsan/meta_CrossCortical/ROSMAP/")
barcodes<-as.character(pheno.RM$X)
RM<- readEPIC(idatPath=idatPath, barcodes=barcodes, parallel = TRUE)
dim(RM)#485577      720
identical(colnames(RM), as.character(pheno.RM$X))

tissue.PFC<-read.csv("samplesROSMAP.csv")
rownames(tissue.PFC)<-as.character(tissue.PFC$Full_Sentrix)
pheno<-tissue.PFC[barcodes,]
pheno<-cbind(pheno,pheno.RM )
identical(colnames(RM), rownames(pheno))
#[1] TRUE

targets <- data.frame(barcodes)
idatPath <- c("/mnt/data1/Ehsan/meta_CrossCortical/ROSMAP/IDAT")
colnames(targets) <- "sentrixFull"
targets$Basename <- with(targets, file.path(idatPath, sentrixFull))

# Construct RGset from IDAT files
rgSet <- read.metharray.exp(targets = targets, verbose = TRUE, extended = FALSE)

setwd("/mnt/data1/Ehsan/meta_CrossCortical/RAW_DATASETS")
save(RM,rgSet,pheno,
     file="ROSMAP__RgSet_Methylumi_Pheno_PFC720.Rdata")



######################################################################
######################################################################
pheno<-read.csv("/mnt/data1/Bex/Meta/SVAregressed/Allphenosub.csv")
pheno.RM<-pheno[which(pheno$Cohort=="ROSMAP" & pheno$region=="PFC"),]
rm(pheno)
###
idatPath <- c("/mnt/data1/Ehsan/meta_CrossCortical/ROSMAP/")
barcodes<-as.character(pheno.RM$X)
RM<- readEPIC(idatPath=idatPath, barcodes=barcodes, parallel = TRUE)
dim(RM)#485577      720
identical(colnames(RM), as.character(pheno.RM$X))

tissue.PFC<-read.csv("samplesROSMAP.csv")
rownames(tissue.PFC)<-as.character(tissue.PFC$Full_Sentrix)
pheno<-tissue.PFC[barcodes,]
pheno<-cbind(pheno,pheno.RM )
identical(colnames(RM), rownames(pheno))
#[1] TRUE

targets <- data.frame(barcodes)
idatPath <- c("/mnt/data1/Ehsan/meta_CrossCortical/ROSMAP/IDAT")
colnames(targets) <- "sentrixFull"
targets$Basename <- with(targets, file.path(idatPath, sentrixFull))

# Construct RGset from IDAT files
rgSet <- read.metharray.exp(targets = targets, verbose = TRUE, extended = FALSE)

setwd("/mnt/data1/Ehsan/meta_CrossCortical/RAW_DATASETS")
save(RM,rgSet,pheno,
     file="ROSMAP__RgSet_Methylumi_Pheno_PFC720.Rdata")
######################################################################
pheno.LBB2.CER<-read.csv("/mnt/data1/Bex/Meta/CB/LBB2Hpheno.csv")
idatPath <- c("/mnt/data1/Ehsan/meta_CrossCortical/LBB1/CER/idats/")

barcodes<-as.character(pheno.LBB2.CER$X)
LBB2.CER<- readEPIC(idatPath=idatPath, barcodes=barcodes, parallel = F)
dim(LBB2.CER)#485577      720
identical(colnames(LBB2.CER), as.character(pheno.LBB2.CER$X))

targets <- data.frame(barcodes)
colnames(targets) <- "sentrixFull"
targets$Basename <- with(targets, file.path(idatPath, sentrixFull))

# Construct RGset from IDAT files
rgSet <- read.metharray.exp(targets = targets, verbose = TRUE, extended = FALSE)

setwd("/mnt/data1/Ehsan/meta_CrossCortical/RAW_DATASETS")
save(LBB2.CER,rgSet,pheno.LBB2.CER,
     file="LBB2__RgSet_Methylumi_Pheno_CER95.Rdata")

######################################################################
rm(list=ls())
pheno.LBB1.CER<-read.csv("/mnt/data1/Bex/Meta/CB/LBB1Hpheno.csv")
idatPath <- c("/mnt/data1/Bex/AD/iDATs/LBBH")
barcodes<-as.character(pheno.LBB1.CER$X.1)
LBB1.CER<- readEPIC(idatPath=idatPath, barcodes=barcodes, parallel = F)
dim(LBB1.CER)#485577      112
identical(colnames(LBB1.CER), as.character(pheno.LBB1.CER$X.1))

targets <- data.frame(barcodes)
colnames(targets) <- "sentrixFull"
targets$Basename <- with(targets, file.path(idatPath, sentrixFull))

# Construct RGset from IDAT files
rgSet <- read.metharray.exp(targets = targets, verbose = TRUE, extended = FALSE)

setwd("/mnt/data1/Ehsan/meta_CrossCortical/RAW_DATASETS")
save(LBB1.CER,rgSet,pheno.LBB1.CER,
     file="LBB1__RgSet_Methylumi_Pheno_CER112.Rdata")
######################################################################


######################################################################
rm(list=ls())
pheno<-read.csv("/mnt/data1/Bex/AD/iDATs/Arizona/idat_files/ArizonaPhenoData.csv")
pheno.AZ.CER<-pheno[which(pheno$Region=="CBL"),]

idatPath <- c("/mnt/data1/Ehsan/meta_CrossCortical/AZ1/IDATS/")
barcodes<-as.character(pheno.AZ.CER$Complete.Barcode)
AZ.CER<- readEPIC(idatPath=idatPath, barcodes=barcodes, parallel = F)
dim(AZ.CER)#485577      404
identical(colnames(AZ.CER), as.character(pheno.AZ.CER$Complete.Barcode))

targets <- data.frame(barcodes)
colnames(targets) <- "sentrixFull"
targets$Basename <- with(targets, file.path(idatPath, sentrixFull))

# Construct RGset from IDAT files
rgSet <- read.metharray.exp(targets = targets, verbose = TRUE, extended = FALSE)

setwd("/mnt/data1/Ehsan/meta_CrossCortical/RAW_DATASETS")
save(AZ.CER,rgSet,pheno.AZ.CER,
     file="AZ__RgSet_Methylumi_Pheno_CER404.Rdata")
######################################################################








rm(list=ls())
library(minfi)
library(parallel)
library(wateRmelon)
######################################################################
setwd("/mnt/data1/Ehsan/meta_CrossCortical")
######################################################################
pheno<-read.table("/mnt/data1/Ehsan/meta_CrossCortical/SAARLAND/SAARLAND/pheno/SampleSheet_BulkTissue.txt", sep="", header=TRUE)
pheno.FC<-pheno[which(pheno$BrainRegion=="FrontalCortex"),]
rm(pheno)
######################################################################
idatPath <- c("/mnt/data1/Ehsan/meta_CrossCortical/SAARLAND/SAARLAND/rawdata_BulkTissue/")
barcodes<-as.character(pheno.FC$SentrixAdress)
SAAR.FC<- readEPIC(idatPath=idatPath, barcodes=barcodes, parallel = TRUE)
identical(colnames(SAAR.FC), as.character(pheno.FC$SentrixAdress))
######################################################################
######################################################################
targets <- data.frame(barcodes)
colnames(targets) <- "sentrixFull"
targets$Basename <- with(targets, file.path(idatPath, sentrixFull))

# Construct RGset from IDAT files
rgSet <- read.metharray.exp(targets = targets, verbose = TRUE, extended = FALSE)
setwd("/mnt/data1/Ehsan/meta_CrossCortical/RAW_DATASETS")
save(SAAR.FC,rgSet,pheno.FC,
     file="SAAR__RgSet_Methylumi_Pheno_FC63.Rdata")
######################################################################
######################################################################

rm(list=ls())
library(minfi)
library(parallel)
library(wateRmelon)
######################################################################
######################################################################
pheno<-read.table("/mnt/data1/Ehsan/meta_CrossCortical/SAARLAND/SAARLAND/pheno/SampleSheet_BulkTissue.txt", sep="", header=TRUE)
pheno.TC<-pheno[which(pheno$BrainRegion=="TemporalCortex"),]
rm(pheno)
######################################################################
idatPath <- c("/mnt/data1/Ehsan/meta_CrossCortical/SAARLAND/SAARLAND/rawdata_BulkTissue/")
barcodes<-as.character(pheno.TC$SentrixAdress)
SAAR.TC<- readEPIC(idatPath=idatPath, barcodes=barcodes, parallel = TRUE)
identical(colnames(SAAR.TC), as.character(pheno.TC$SentrixAdress))
######################################################################
######################################################################
targets <- data.frame(barcodes)
colnames(targets) <- "sentrixFull"
targets$Basename <- with(targets, file.path(idatPath, sentrixFull))

# Construct RGset from IDAT files
rgSet <- read.metharray.exp(targets = targets, verbose = TRUE, extended = FALSE)
setwd("/mnt/data1/Ehsan/meta_CrossCortical/RAW_DATASETS")
save(SAAR.TC,rgSet,pheno.TC,
     file="SAAR__RgSet_Methylumi_Pheno_TC65.Rdata")
######################################################################
######################################################################




rm(list=ls())
library(minfi)
library(parallel)
library(wateRmelon)
######################################################################
######################################################################
pheno<-read.table("/mnt/data1/Ehsan/meta_CrossCortical/SAARLAND/SAARLAND/pheno/SampleSheet_NeuNsorts.txt", sep="", header=TRUE)
######################################################################
idatPath <- c("/mnt/data1/Ehsan/meta_CrossCortical/SAARLAND/SAARLAND/rawdata_NeuNTissue/")
barcodes<-as.character(pheno$SentrixAdress)
SAAR.N<- readEPIC(idatPath=idatPath, barcodes=barcodes, parallel = TRUE)
identical(colnames(SAAR.N), as.character(pheno$SentrixAdress))
######################################################################
######################################################################
targets <- data.frame(barcodes)
colnames(targets) <- "sentrixFull"
targets$Basename <- with(targets, file.path(idatPath, sentrixFull))

# Construct RGset from IDAT files
rgSet <- read.metharray.exp(targets = targets, verbose = TRUE, extended = FALSE)
save(rgSet, pheno, file="rgSet62.SAAR.N_pheno_.Rdata")
setwd("/mnt/data1/Ehsan/meta_CrossCortical/RAW_DATASETS")
save(SAAR.N,rgSet,pheno,
     file="SAAR__RgSet_Methylumi_Pheno_N62.Rdata")
######################################################################
######################################################################

