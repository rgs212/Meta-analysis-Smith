rm(list=ls())
############################################################################
############################################################################
###loading required libraries
library(minfi)
library(parallel)
library(wateRmelon)
library(corrplot)
library(sva)
###reading the demographics and phenotipic aoutcomes
setwd("/mnt/data1/Ehsan/meta_CrossCortical/RAW_DATASETS")
load("SAAR__RgSet_Methylumi_Pheno_FC63.Rdata")
###################################
#######
setwd("/mnt/data1/Ehsan/meta_CrossCortical/SAARLAND/FC/QC")
###### QC.1#######
rgSet<-rgSet[,which(pheno.FC$Age > 60)]
Methylumiset<-SAAR.FC[,which(pheno.FC$Age > 60)]
pheno<-pheno.FC[which(pheno.FC$Age > 60),]
identical(colnames(rgSet), as.character(pheno$SentrixAdress))
rm(pheno.FC, SAAR.FC)


con.red<-getRed(rgSet)[getProbeInfo(rgSet, type = "Control")$Address,]
rownames(con.red)<-getProbeInfo(rgSet, type = "Control")$ExtendedType
con.green<-getGreen(rgSet)[getProbeInfo(rgSet, type = "Control")$Address,]
rownames(con.green)<-getProbeInfo(rgSet, type = "Control")$ExtendedType

neg.red.bs<-con.red[which(getProbeInfo(rgSet, type = "Control")$Type == "NEGATIVE"),]
neg.green.bs<-con.green[which(getProbeInfo(rgSet, type = "Control")$Type == "NEGATIVE"),]

neg.red.mean.bs.lc<-colMeans(neg.red.bs)
neg.green.mean.bs.lc<-colMeans(neg.green.bs)
green.mean.bs.lc<-colMeans(getGreen(rgSet))
red.mean.bs.lc<-colMeans(getRed(rgSet))

pdf("NegativeCtrlbckg.pdf", height=8, width=8)
par(mfrow=c(1,2))
plot(neg.red.mean.bs.lc, xlab = "Sample", ylab = "Intensity-red", col="grey", pch = 16, ylim = c(0, 12000))
par(new=TRUE)
plot(red.mean.bs.lc, xlab = "Sample", ylab = "", col="red", pch = 16, ylim = c(0, 12000))

plot(neg.green.mean.bs.lc, xlab = "Sample", ylab = "Intensity-green", col="grey", pch = 16, ylim = c(0, 12000))
par(new=TRUE)
plot(green.mean.bs.lc, xlab = "Sample", ylab = "", col="Green", pch = 16, ylim = c(0, 12000))
dev.off()

BAD.BckgIntensityDiff<-c(rownames(data.frame(which(neg.red.mean.bs.lc >1000))),rownames(data.frame(which(red.mean.bs.lc < 2000))))

avgPval<- colMeans(pvals(Methylumiset))

pdf("BckgIntensityDiffPval.pdf")
barplot(avgPval,xlab=NULL)
dev.off()

BAD.BckgIntensityDiffPval<-c(as.character(colnames(Methylumiset[,which(avgPval> 0.005)])),BAD.BckgIntensityDiff)
identical(colnames(Methylumiset), as.character(pheno$SentrixAdress))
#[1] TRUE

###### QC xxx.  MULTIDIMENSIONAL SCALING OF SEX CHROMOSOMES #######
meth      <- methylated(Methylumiset)
unmeth   <- unmethylated(Methylumiset)

M.mean<-apply(meth, 2, mean)
U.mean<-apply(unmeth, 2, mean)

Mean.M<-mean(M.mean)
Mean.U<-mean(U.mean)
sd.M<-sd(M.mean)
sd.U<-sd(U.mean)

BAD.samples.Intensity<-as.character(pheno[which(Mean.M +3*sd.M < M.mean |
                                                  M.mean < Mean.M-3*sd.M | 
                                                  Mean.U +3*sd.U < U.mean |
                                                  U.mean< Mean.U-3*sd.U),]$SentrixAdress)

pdf("IntensityByFactor.pdf", height=6, width=12)
par(mfrow=c(1,3))

#plots of median sample intensites and threshold lines to see if signals were sufficiently high
Sex<-as.factor(pheno$Gender)
plot(M.mean, U.mean, pch = 16, xlab = "Mean M intensity", ylab = "Mean U intensity", col = rainbow(nlevels(Sex))[factor(Sex)], main="BS Signal Intensities Coloured by Sex")
legend("topleft", legend=levels(factor(Sex)), col = rainbow(nlevels(Sex)), pch = 16, cex=0.6)

Position<-as.factor(substr(pheno$SentrixAdress,12,18))
plot(M.mean, U.mean, pch = 16, xlab = "Mean M intensity", ylab = "Mean U intensity", col = rainbow(nlevels(Position))[factor(Position)], main="BS Signal Intensities Coloured by Position")
legend("topleft", legend=levels(factor(Position)), col = rainbow(nlevels(Position)), pch = 16, cex=0.6)

Chip<-as.factor(substr(pheno$SentrixAdress,1,11))
plot(M.mean, U.mean, pch = 16, xlab = "Mean M intensity", ylab = "Mean U intensity", col = rainbow(nlevels(Chip))[factor(Chip)], main="BS Signal Intensities Coloured by Chip")
legend("topleft", legend=levels(factor(Chip)), col = rainbow(nlevels(Chip)), pch = 16, cex=0.6)
dev.off()

############################################################################
############################################################################
###### QC xxx.  MULTIDIMENSIONAL SCALING OF SEX CHROMOSOMES #######
#calculate bisulphite conversion statistic using bscon (from wateRmelon) and plots histogram
bs<-bscon(Methylumiset)
pdf("HistogramBisulphiteConversionStatistics.pdf")
hist(bs, xlab = "Mean % BS conversion", main = "")
dev.off()
BAD.samples.Bisulifite<-as.character(colnames(Methylumiset[,which(bs< 80)]))
############################################################################
############################################################################
###### QC xxx.  MULTIDIMENSIONAL SCALING OF SEX CHROMOSOMES #######
#this is a check that gender is correct (ie for sample mix up!)
###identify which PCAs correlate with sex
findGenderPC<-function(betas, sex, npcs = 20){
  
  betas.com<-betas[complete.cases(betas),]
  pca<-prcomp(betas.com)
  
  pca.cor<-rep(NA, npcs)
  for(i in 1:npcs){
    pca.cor[i]<-cor(pca$rotation[,i], as.numeric(as.factor(sex)))
  }
  top<-order(abs(pca.cor), decreasing = TRUE)[1]
  second<-order(abs(pca.cor), decreasing = TRUE)[2]
  print(paste("Top correlated principal components with sex:", top, ",", second))
  plot(pca$rotation[,top], pca$rotation[,second], pch = 16, col = c("magenta", "blue")[as.factor(sex)], xlab = paste("PC", top), ylab = paste("PC", second))
  legend("topright", levels(as.factor(sex)), pch = 16, col = c("magenta", "blue"))
  predSex<-rep(NA, length(sex))
  options.sex<-levels(as.factor(sex))
  if(abs(pca.cor[top]) > 0.9){
    print("Top PC has r > 0.9 with sex so good enough to confirm reported sexes")
    if(sign(pca.cor[top]) == 1){
      
      predSex[which(pca$rotation[,top] < 0)]<-options.sex[1]
      predSex[which(pca$rotation[,top] > 0)]<-options.sex[2]
    } else {
      
      predSex[which(pca$rotation[,top] < 0)]<-options.sex[2]
      predSex[which(pca$rotation[,top] > 0)]<-options.sex[1]
    }
  }
  return(predSex)
}


clusterGender<-function(betas, sex, npcs = 20, thres = 0.5){
  if(npcs > ncol(betas)){
    npcs<-ncol(betas)
    print(paste("As only", ncol(betas), "samples, can look look at that number of PCs", sep = " "))
  }
  betas.com<-betas[complete.cases(betas),]
  pca<-prcomp(betas.com)
  
  pca.cor<-rep(NA, npcs)
  for(i in 1:npcs){
    pca.cor[i]<-cor(pca$rotation[,i], as.numeric(as.factor(sex)))
  }
  predSex<-kmeans(pca$rotation[,which(abs(pca.cor) > thres)],2)$cluster
  return(predSex)
}

#runs the two gender predicting functions
pdf("PredictGender.pdf")
predSex<-findGenderPC(betas(Methylumiset), pheno$Gender)
dev.off()
BAD.samples.sexMissmatch<-as.character(pheno[which(predSex!=pheno$Gender),]$SentrixAdress)
############################################################################
############################################################################
############################################################################
#Uses SNP probes on the array to find genetic correlations between samples
betas<-betas(Methylumiset)
betas.rs<-betas[grep("rs", rownames(betas)),]
snpCor<-cor(betas.rs, use="complete.obs")

#ignores correlations between a sample and itself
for(i in 1:ncol(betas.rs)){
  snpCor[i,i]<-NA
}
#finds the maximum correlation values between samples and plots a histogram
#note: correlations of ~0.5 do not mean siblings as there are too few SNPs for this
corMax<-apply(snpCor, 1, max, na.rm = TRUE)
pdf("SNPCorrelation.MS.LC.pdf")
hist(corMax, xlab = "Max. correlation with all other samples", main = "")
dev.off()
genetic.cor<-data.frame(corMax)
rownames(genetic.cor)<-colnames(betas.rs)
#Visscher, P.M., et al. PLoS Genet, 2006. 2(3): p. e41.
############################################################################
pheno[which(genetic.cor[,1] > 0.65),]
genetic.cor[which(genetic.cor[,1] > 0.65),]
BAD.genetic.relationship<-rownames(genetic.cor)[which(genetic.cor[,1] > 0.65)]

############################################################################
############################################################################
############################################################################
bad<-c(BAD.samples.Bisulifite,BAD.BckgIntensityDiffPval,BAD.BckgIntensityDiff,BAD.samples.Intensity,BAD.genetic.relationship, BAD.samples.sexMissmatch)
bad<-unique(bad)
save(bad,BAD.samples.Bisulifite,BAD.BckgIntensityDiffPval,BAD.BckgIntensityDiff,BAD.samples.Intensity,BAD.genetic.relationship, BAD.samples.sexMissmatch,
     file="bad.samples.Rdata")
############################################################################
dat.pf<-pfilter(Methylumiset[,-which(colnames(Methylumiset)%in%bad)], perc=5)
#0 samples having 5 % of sites with a detection p-value greater than 0.05 were removed 
#Samples removed:  
#   176 sites were removed as beadcount <3 in 5 % of samples 
#   2894 sites having 1 % of samples with a detection p-value greater than 0.05 were removed 
dat.dasen<-dasen(dat.pf)
betas.dasen<-betas(dat.dasen)
pheno<-pheno[-which(colnames(Methylumiset)%in%bad),]
identical(colnames(betas.dasen), as.character(pheno$SentrixAdress))
#[1] TRUE

tiff("density.betas.tiff")
plot(density(betas.dasen[,1],na.rm=T),main="Distribution of Beta Values", xlim=c(0,1), col=c("black"), lty=1)
for (n in 2:ncol(betas.dasen)){
  lines(density(betas.dasen[,n],na.rm=T))}
dev.off()

cross<-read.csv("/mnt/data1/450K_reference/CrossHybridisingProbesPriceORWeksberg.csv")
crosslist<-cross[,2]
length(crosslist)
#[1]43233  

d<-betas.dasen[ ! rownames(betas.dasen) %in% crosslist, ]
dim(d) # [1] 440480      4   

#Made a file from the Weksburg_CommonSNP_EuAf_within10bpSBE column and Weksburg_CommonSNP_EuAf_inProbeSeq in the SNPs inprobeANno.csv. Called this snpprobes.csv. It contains 
Af_SNP<-read.csv("/mnt/data1/450K_reference/snp_chyb_Af_out.csv")
dim(Af_SNP)
# [1] 53794     1
SNPlist<-as.character(Af_SNP[,2]) 
d1<- d[ ! rownames(d) %in% SNPlist, ]
dim(d1) # [1] 425721    141

EU_SNP<-read.csv("/mnt/data1/450K_reference/snp_chyb_EU_out.csv")
dim(EU_SNP)
# [1] 53794     1
SNPlist<-as.character(EU_SNP[,2])
d2<- d1[ ! rownames(d1) %in% SNPlist, ]
dim(d2)
#[1] 424958    141
betas<-data.frame(d2)
identical(substr(colnames(betas),2,18),as.character(pheno$SentrixAdress))
colnames(betas)<-as.character(pheno$SentrixAdress)
setwd("/mnt/data1/Ehsan/meta_CrossCortical/NORMALISED_BETA")
save(betas,pheno, file="SAARLAND.FC.48.beta.pheno.Rdata")
#TRUE

resid <- function( row, age, sex, cet ){
  
  fit <- try (
    lm( row ~ age + sex + cet),
    silent=TRUE
  )
  if(inherits(fit,'try-error')) return(rep(NA, length(sex)))
  fit$residuals
} 

cl<- makeCluster(32)

pheno<-pheno[!is.na(pheno$Braak),]
betas<-betas[,as.character(pheno$SentrixAdress)]


colData(rgSet)@listData$sentrixFull <- as.character(colData(rgSet)@listData$sentrixFull)
library(minfi)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("FlowSorted.DLPFC.450k")

library('FlowSorted.DLPFC.450k')
cet.saar.fc <- estimateCellCounts(rgSet[,colnames(betas)], 
                              compositeCellType = "DLPFC",
                              processMethod = "auto", 
                              probeSelect = "auto",
                              cellTypes = c("NeuN_neg", "NeuN_pos"),
                              referencePlatform = NULL,
                              returnAll = FALSE, 
                              meanPlot = TRUE, 
                              verbose = TRUE)
cellprop<-(cet.saar.fc[,2]/(cet.saar.fc[,2]+cet.saar.fc[,1]))*100
pheno<-cbind(pheno, cellprop)
pheno$Braak<-pheno$Braak


dat.reg<- {              
  sex    	<- 	as.factor(pheno$Gender)
  age    	<- 	pheno$Age 
  cet		<- as.numeric(pheno$cellprop)
  t(parApply(cl, betas, 1, resid, age, sex, cet))
}
head(dat.reg[1:5,1:5])
identical(colnames(dat.reg), as.character(pheno$SentrixAdress))
stopCluster(cl)
setwd("/mnt/data1/Ehsan/meta_CrossCortical/REGRESSED_BETA")
save(dat.reg,pheno, file="SAAR.FC.45.pheno.REGRESSED.Rdata")

mod = model.matrix( ~ as.numeric(Braak),data=pheno)
mod0 = model.matrix(~1,data=pheno)

mydf<-na.omit(dat.reg)

sva<- sva::sva(dat = as.matrix(mydf), mod = mod, mod0 = mod0)

setwd("/mnt/data1/Ehsan/meta_CrossCortical/SVA/")
save(sva, file="sva.regressed.SAAR.fc.Rdata")
load("sva.regressed.SAAR.fc.Rdata")

svs<-sva$sv
pheno.SAAR.fc<-cbind(pheno, svs)

setwd("/mnt/data1/Ehsan/meta_CrossCortical/pheno/")
write.csv(pheno.SAAR.fc,file="pheno.SAAR.fc.45.csv")








rm(list=ls())
############################################################################
############################################################################
###loading required libraries
library(minfi)
library(parallel)
library(wateRmelon)
library(corrplot)
library(sva)
###reading the demographics and phenotipic aoutcomes
setwd("/mnt/data1/Ehsan/meta_CrossCortical/RAW_DATASETS")
load("SAAR__RgSet_Methylumi_Pheno_N62.Rdata")
###################################
#######
setwd("/mnt/data1/Ehsan/meta_CrossCortical/SAARLAND/NeuN/QC")
###### QC.1#######
rgSet<-rgSet[,which(pheno$Age > 60)]
Methylumiset<-SAAR.N[,which(pheno$Age > 60)]
pheno<-pheno[which(pheno$Age > 60),]
identical(colnames(rgSet), as.character(pheno$SentrixAdress))
rm(SAAR.N)

con.red<-getRed(rgSet)[getProbeInfo(rgSet, type = "Control")$Address,]
rownames(con.red)<-getProbeInfo(rgSet, type = "Control")$ExtendedType
con.green<-getGreen(rgSet)[getProbeInfo(rgSet, type = "Control")$Address,]
rownames(con.green)<-getProbeInfo(rgSet, type = "Control")$ExtendedType

neg.red.bs<-con.red[which(getProbeInfo(rgSet, type = "Control")$Type == "NEGATIVE"),]
neg.green.bs<-con.green[which(getProbeInfo(rgSet, type = "Control")$Type == "NEGATIVE"),]

neg.red.mean.bs.lc<-colMeans(neg.red.bs)
neg.green.mean.bs.lc<-colMeans(neg.green.bs)
green.mean.bs.lc<-colMeans(getGreen(rgSet))
red.mean.bs.lc<-colMeans(getRed(rgSet))

pdf("NegativeCtrlbckg.pdf", height=8, width=8)
par(mfrow=c(1,2))
plot(neg.red.mean.bs.lc, xlab = "Sample", ylab = "Intensity-red", col="grey", pch = 16, ylim = c(0, 12000))
par(new=TRUE)
plot(red.mean.bs.lc, xlab = "Sample", ylab = "", col="red", pch = 16, ylim = c(0, 12000))

plot(neg.green.mean.bs.lc, xlab = "Sample", ylab = "Intensity-green", col="grey", pch = 16, ylim = c(0, 12000))
par(new=TRUE)
plot(green.mean.bs.lc, xlab = "Sample", ylab = "", col="Green", pch = 16, ylim = c(0, 12000))
dev.off()

BAD.BckgIntensityDiff<-c(rownames(data.frame(which(neg.red.mean.bs.lc >1000))),rownames(data.frame(which(red.mean.bs.lc < 2000))))

avgPval<- colMeans(pvals(Methylumiset))

pdf("BckgIntensityDiffPval.pdf")
barplot(avgPval,xlab=NULL)
dev.off()

BAD.BckgIntensityDiffPval<-c(as.character(colnames(Methylumiset[,which(avgPval> 0.005)])),BAD.BckgIntensityDiff)
identical(colnames(Methylumiset), as.character(pheno$SentrixAdress))
#[1] TRUE

###### QC xxx.  MULTIDIMENSIONAL SCALING OF SEX CHROMOSOMES #######
meth      <- methylated(Methylumiset)
unmeth   <- unmethylated(Methylumiset)

M.mean<-apply(meth, 2, mean)
U.mean<-apply(unmeth, 2, mean)

Mean.M<-mean(M.mean)
Mean.U<-mean(U.mean)
sd.M<-sd(M.mean)
sd.U<-sd(U.mean)

BAD.samples.Intensity<-as.character(pheno[which(Mean.M +3*sd.M < M.mean |
                                                  M.mean < Mean.M-3*sd.M | 
                                                  Mean.U +3*sd.U < U.mean |
                                                  U.mean< Mean.U-3*sd.U),]$SentrixAdress)

pdf("IntensityByFactor.pdf", height=6, width=12)
par(mfrow=c(1,3))

#plots of median sample intensites and threshold lines to see if signals were sufficiently high
Sex<-as.factor(pheno$Gender)
plot(M.mean, U.mean, pch = 16, xlab = "Mean M intensity", ylab = "Mean U intensity", col = rainbow(nlevels(Sex))[factor(Sex)], main="BS Signal Intensities Coloured by Sex")
legend("topleft", legend=levels(factor(Sex)), col = rainbow(nlevels(Sex)), pch = 16, cex=0.6)

Position<-as.factor(substr(pheno$SentrixAdress,12,18))
plot(M.mean, U.mean, pch = 16, xlab = "Mean M intensity", ylab = "Mean U intensity", col = rainbow(nlevels(Position))[factor(Position)], main="BS Signal Intensities Coloured by Position")
legend("topleft", legend=levels(factor(Position)), col = rainbow(nlevels(Position)), pch = 16, cex=0.6)

Chip<-as.factor(substr(pheno$SentrixAdress,1,11))
plot(M.mean, U.mean, pch = 16, xlab = "Mean M intensity", ylab = "Mean U intensity", col = rainbow(nlevels(Chip))[factor(Chip)], main="BS Signal Intensities Coloured by Chip")
legend("topleft", legend=levels(factor(Chip)), col = rainbow(nlevels(Chip)), pch = 16, cex=0.6)
dev.off()

############################################################################
############################################################################
###### QC xxx.  MULTIDIMENSIONAL SCALING OF SEX CHROMOSOMES #######
#calculate bisulphite conversion statistic using bscon (from wateRmelon) and plots histogram
bs<-bscon(Methylumiset)
pdf("HistogramBisulphiteConversionStatistics.pdf")
hist(bs, xlab = "Mean % BS conversion", main = "")
dev.off()
BAD.samples.Bisulifite<-as.character(colnames(Methylumiset[,which(bs< 80)]))
############################################################################
############################################################################
###### QC xxx.  MULTIDIMENSIONAL SCALING OF SEX CHROMOSOMES #######
#this is a check that gender is correct (ie for sample mix up!)
###identify which PCAs correlate with sex
findGenderPC<-function(betas, sex, npcs = 20){
  
  betas.com<-betas[complete.cases(betas),]
  pca<-prcomp(betas.com)
  
  pca.cor<-rep(NA, npcs)
  for(i in 1:npcs){
    pca.cor[i]<-cor(pca$rotation[,i], as.numeric(as.factor(sex)))
  }
  top<-order(abs(pca.cor), decreasing = TRUE)[1]
  second<-order(abs(pca.cor), decreasing = TRUE)[2]
  print(paste("Top correlated principal components with sex:", top, ",", second))
  plot(pca$rotation[,top], pca$rotation[,second], pch = 16, col = c("magenta", "blue")[as.factor(sex)], xlab = paste("PC", top), ylab = paste("PC", second))
  legend("topright", levels(as.factor(sex)), pch = 16, col = c("magenta", "blue"))
  predSex<-rep(NA, length(sex))
  options.sex<-levels(as.factor(sex))
  if(abs(pca.cor[top]) > 0.9){
    print("Top PC has r > 0.9 with sex so good enough to confirm reported sexes")
    if(sign(pca.cor[top]) == 1){
      
      predSex[which(pca$rotation[,top] < 0)]<-options.sex[1]
      predSex[which(pca$rotation[,top] > 0)]<-options.sex[2]
    } else {
      
      predSex[which(pca$rotation[,top] < 0)]<-options.sex[2]
      predSex[which(pca$rotation[,top] > 0)]<-options.sex[1]
    }
  }
  return(predSex)
}


clusterGender<-function(betas, sex, npcs = 20, thres = 0.5){
  if(npcs > ncol(betas)){
    npcs<-ncol(betas)
    print(paste("As only", ncol(betas), "samples, can look look at that number of PCs", sep = " "))
  }
  betas.com<-betas[complete.cases(betas),]
  pca<-prcomp(betas.com)
  
  pca.cor<-rep(NA, npcs)
  for(i in 1:npcs){
    pca.cor[i]<-cor(pca$rotation[,i], as.numeric(as.factor(sex)))
  }
  predSex<-kmeans(pca$rotation[,which(abs(pca.cor) > thres)],2)$cluster
  return(predSex)
}

#runs the two gender predicting functions
pdf("PredictGender.pdf")
predSex<-findGenderPC(betas(Methylumiset), pheno$Gender)
dev.off()
BAD.samples.sexMissmatch<-NULL
############################################################################
############################################################################
############################################################################
#Uses SNP probes on the array to find genetic correlations between samples
betas<-betas(Methylumiset)
betas.rs<-betas[grep("rs", rownames(betas)),]
snpCor<-cor(betas.rs, use="complete.obs")

#ignores correlations between a sample and itself
for(i in 1:ncol(betas.rs)){
  snpCor[i,i]<-NA
}
#finds the maximum correlation values between samples and plots a histogram
#note: correlations of ~0.5 do not mean siblings as there are too few SNPs for this
corMax<-apply(snpCor, 1, max, na.rm = TRUE)
pdf("SNPCorrelation.MS.LC.pdf")
hist(corMax, xlab = "Max. correlation with all other samples", main = "")
dev.off()
genetic.cor<-data.frame(corMax)
rownames(genetic.cor)<-colnames(betas.rs)
#Visscher, P.M., et al. PLoS Genet, 2006. 2(3): p. e41.
############################################################################
pheno[which(genetic.cor[,1] > 0.65),]
genetic.cor[which(genetic.cor[,1] > 0.65),]
BAD.genetic.relationship<-NULL

############################################################################
############################################################################
############################################################################
bad<-c(BAD.samples.Bisulifite,BAD.BckgIntensityDiffPval,BAD.BckgIntensityDiff,BAD.samples.Intensity,BAD.genetic.relationship, BAD.samples.sexMissmatch)
bad<-unique(bad)
save(bad,BAD.samples.Bisulifite,BAD.BckgIntensityDiffPval,BAD.BckgIntensityDiff,BAD.samples.Intensity,BAD.genetic.relationship, BAD.samples.sexMissmatch,
     file="bad.samples.Rdata")
############################################################################
dat.pf<-pfilter(Methylumiset[,-which(colnames(Methylumiset)%in%bad)], perc=5)
#0 samples having 5 % of sites with a detection p-value greater than 0.05 were removed 
#Samples removed:  
#   7396 sites were removed as beadcount <3 in 5 % of samples 
#   2162 sites having 1 % of samples with a detection p-value greater than 0.05 were removed 
dat.dasen<-dasen(dat.pf)
betas.dasen<-betas(dat.dasen)
pheno<-pheno[-which(colnames(Methylumiset)%in%bad),]
identical(colnames(betas.dasen), as.character(pheno$SentrixAdress))
#[1] TRUE

tiff("density.betas.tiff")
plot(density(betas.dasen[,1],na.rm=T),main="Distribution of Beta Values", xlim=c(0,1), col=c("black"), lty=1)
for (n in 2:ncol(betas.dasen)){
  lines(density(betas.dasen[,n],na.rm=T))}
dev.off()

cross<-read.csv("/mnt/data1/450K_reference/CrossHybridisingProbesPriceORWeksberg.csv")
crosslist<-cross[,2]
length(crosslist)
#[1]43233  

d<-betas.dasen[ ! rownames(betas.dasen) %in% crosslist, ]
dim(d) # [1] 440480      4   

#Made a file from the Weksburg_CommonSNP_EuAf_within10bpSBE column and Weksburg_CommonSNP_EuAf_inProbeSeq in the SNPs inprobeANno.csv. Called this snpprobes.csv. It contains 
Af_SNP<-read.csv("/mnt/data1/450K_reference/snp_chyb_Af_out.csv")
dim(Af_SNP)
# [1] 53794     1
SNPlist<-as.character(Af_SNP[,2]) 
d1<- d[ ! rownames(d) %in% SNPlist, ]
dim(d1) # [1] 425721    141

EU_SNP<-read.csv("/mnt/data1/450K_reference/snp_chyb_EU_out.csv")
dim(EU_SNP)
# [1] 53794     1
SNPlist<-as.character(EU_SNP[,2])
d2<- d1[ ! rownames(d1) %in% SNPlist, ]
dim(d2)
#[1] 424958    141
betas<-data.frame(d2)
identical(substr(colnames(betas),2,18),as.character(pheno$SentrixAdress))
colnames(betas)<-as.character(pheno$SentrixAdress)
setwd("/mnt/data1/Ehsan/meta_CrossCortical/NORMALISED_BETA")
save(betas,pheno, file="SAARLAND.N.51.beta.pheno.Rdata")
#TRUE

resid <- function( row, age, sex){
  
  fit <- try (
    lm( row ~ age + sex),
    silent=TRUE
  )
  if(inherits(fit,'try-error')) return(rep(NA, length(sex)))
  fit$residuals
} 

cl<- makeCluster(32)

phenoNP<-pheno[which(pheno$NeuNsort=="NeuNnegative"),]
betasNP<-betas[,as.character(phenoNP$SentrixAdress)]

phenoNN<-pheno[which(pheno$NeuNsort=="NeuNpositive"),]
betasNN<-betas[,as.character(phenoNN$SentrixAdress)]


dat.reg.NP<- {              
  sex    	<- 	as.factor(phenoNP$Gender)
  age    	<- 	phenoNP$Age 
  t(parApply(cl, betasNP, 1, resid, age, sex))
}
head(dat.reg.NP[1:5,1:5])
identical(colnames(dat.reg.NP), as.character(phenoNP$SentrixAdress))

dat.reg.NN<- {              
  sex    	<- 	as.factor(phenoNN$Gender)
  age    	<- 	phenoNN$Age 
  t(parApply(cl, betasNN, 1, resid, age, sex))
}
head(dat.reg.NN[1:5,1:5])
identical(colnames(dat.reg.NN), as.character(phenoNN$SentrixAdress))



setwd("/mnt/data1/Ehsan/meta_CrossCortical/REGRESSED_BETA")
save(dat.reg.NN,dat.reg.NP ,phenoNP, phenoNN, file="SAAR.NN26.NP25.pheno.REGRESSED.Rdata")


mod = model.matrix( ~ as.numeric(Braak),data=phenoNN)
mod0 = model.matrix(~1,data=phenoNN)
mydf<-na.omit(dat.reg.NN)
sva.NN<- sva::sva(dat = as.matrix(mydf), mod = mod, mod0 = mod0)

mod = model.matrix( ~ as.numeric(Braak),data=phenoNP)
mod0 = model.matrix(~1,data=phenoNP)
mydf<-na.omit(dat.reg.NP)
sva.NP<- sva::sva(dat = as.matrix(mydf), mod = mod, mod0 = mod0)

setwd("/mnt/data1/Ehsan/meta_CrossCortical/SVA/")
save(sva.NN, sva.NP, file="sva.regressed.SAAR.NN.NP.Rdata")

svs<-sva.NN$sv
pheno.SAAR.NN<-cbind(phenoNN, svs)

svs<-sva.NP$sv
pheno.SAAR.NP<-cbind(phenoNP, svs)


setwd("/mnt/data1/Ehsan/meta_CrossCortical/pheno/")
write.csv(pheno.SAAR.NN,file="pheno.SAAR.NN.26.csv")
write.csv(pheno.SAAR.NP,file="pheno.SAAR.NP.25.csv")






rm(list=ls())
############################################################################
############################################################################
###loading required libraries
library(minfi)
library(parallel)
library(wateRmelon)
library(corrplot)
library(sva)
###reading the demographics and phenotipic aoutcomes
setwd("/mnt/data1/Ehsan/meta_CrossCortical/RAW_DATASETS")
load("SAAR__RgSet_Methylumi_Pheno_TC65.Rdata")
###################################
#######
setwd("/mnt/data1/Ehsan/meta_CrossCortical/SAARLAND/FC/QC")
###### QC.1#######
rgSet<-rgSet[,which(pheno.TC$Age > 60)]
Methylumiset<-SAAR.TC[,which(pheno.TC$Age > 60)]
pheno<-pheno.TC[which(pheno.TC$Age > 60),]
identical(colnames(rgSet), as.character(pheno$SentrixAdress))
rm(pheno.TC, SAAR.TC)


con.red<-getRed(rgSet)[getProbeInfo(rgSet, type = "Control")$Address,]
rownames(con.red)<-getProbeInfo(rgSet, type = "Control")$ExtendedType
con.green<-getGreen(rgSet)[getProbeInfo(rgSet, type = "Control")$Address,]
rownames(con.green)<-getProbeInfo(rgSet, type = "Control")$ExtendedType

neg.red.bs<-con.red[which(getProbeInfo(rgSet, type = "Control")$Type == "NEGATIVE"),]
neg.green.bs<-con.green[which(getProbeInfo(rgSet, type = "Control")$Type == "NEGATIVE"),]

neg.red.mean.bs.lc<-colMeans(neg.red.bs)
neg.green.mean.bs.lc<-colMeans(neg.green.bs)
green.mean.bs.lc<-colMeans(getGreen(rgSet))
red.mean.bs.lc<-colMeans(getRed(rgSet))

pdf("NegativeCtrlbckg.pdf", height=8, width=8)
par(mfrow=c(1,2))
plot(neg.red.mean.bs.lc, xlab = "Sample", ylab = "Intensity-red", col="grey", pch = 16, ylim = c(0, 12000))
par(new=TRUE)
plot(red.mean.bs.lc, xlab = "Sample", ylab = "", col="red", pch = 16, ylim = c(0, 12000))

plot(neg.green.mean.bs.lc, xlab = "Sample", ylab = "Intensity-green", col="grey", pch = 16, ylim = c(0, 12000))
par(new=TRUE)
plot(green.mean.bs.lc, xlab = "Sample", ylab = "", col="Green", pch = 16, ylim = c(0, 12000))
dev.off()

BAD.BckgIntensityDiff<-c(rownames(data.frame(which(neg.red.mean.bs.lc >1000))),rownames(data.frame(which(red.mean.bs.lc < 2000))))

avgPval<- colMeans(pvals(Methylumiset))

pdf("BckgIntensityDiffPval.pdf")
barplot(avgPval,xlab=NULL)
dev.off()

BAD.BckgIntensityDiffPval<-c(as.character(colnames(Methylumiset[,which(avgPval> 0.005)])),BAD.BckgIntensityDiff)
identical(colnames(Methylumiset), as.character(pheno$SentrixAdress))
#[1] TRUE

###### QC xxx.  MULTIDIMENSIONAL SCALING OF SEX CHROMOSOMES #######
meth      <- methylated(Methylumiset)
unmeth   <- unmethylated(Methylumiset)

M.mean<-apply(meth, 2, mean)
U.mean<-apply(unmeth, 2, mean)

Mean.M<-mean(M.mean)
Mean.U<-mean(U.mean)
sd.M<-sd(M.mean)
sd.U<-sd(U.mean)

BAD.samples.Intensity<-as.character(pheno[which(Mean.M +3*sd.M < M.mean |
                                                  M.mean < Mean.M-3*sd.M | 
                                                  Mean.U +3*sd.U < U.mean |
                                                  U.mean< Mean.U-3*sd.U),]$SentrixAdress)

pdf("IntensityByFactor.pdf", height=6, width=12)
par(mfrow=c(1,3))

#plots of median sample intensites and threshold lines to see if signals were sufficiently high
Sex<-as.factor(pheno$Gender)
plot(M.mean, U.mean, pch = 16, xlab = "Mean M intensity", ylab = "Mean U intensity", col = rainbow(nlevels(Sex))[factor(Sex)], main="BS Signal Intensities Coloured by Sex")
legend("topleft", legend=levels(factor(Sex)), col = rainbow(nlevels(Sex)), pch = 16, cex=0.6)

Position<-as.factor(substr(pheno$SentrixAdress,12,18))
plot(M.mean, U.mean, pch = 16, xlab = "Mean M intensity", ylab = "Mean U intensity", col = rainbow(nlevels(Position))[factor(Position)], main="BS Signal Intensities Coloured by Position")
legend("topleft", legend=levels(factor(Position)), col = rainbow(nlevels(Position)), pch = 16, cex=0.6)

Chip<-as.factor(substr(pheno$SentrixAdress,1,11))
plot(M.mean, U.mean, pch = 16, xlab = "Mean M intensity", ylab = "Mean U intensity", col = rainbow(nlevels(Chip))[factor(Chip)], main="BS Signal Intensities Coloured by Chip")
legend("topleft", legend=levels(factor(Chip)), col = rainbow(nlevels(Chip)), pch = 16, cex=0.6)
dev.off()

############################################################################
############################################################################
###### QC xxx.  MULTIDIMENSIONAL SCALING OF SEX CHROMOSOMES #######
#calculate bisulphite conversion statistic using bscon (from wateRmelon) and plots histogram
bs<-bscon(Methylumiset)
pdf("HistogramBisulphiteConversionStatistics.pdf")
hist(bs, xlab = "Mean % BS conversion", main = "")
dev.off()
BAD.samples.Bisulifite<-as.character(colnames(Methylumiset[,which(bs< 80)]))
############################################################################
############################################################################
###### QC xxx.  MULTIDIMENSIONAL SCALING OF SEX CHROMOSOMES #######
#this is a check that gender is correct (ie for sample mix up!)
###identify which PCAs correlate with sex
findGenderPC<-function(betas, sex, npcs = 20){
  
  betas.com<-betas[complete.cases(betas),]
  pca<-prcomp(betas.com)
  
  pca.cor<-rep(NA, npcs)
  for(i in 1:npcs){
    pca.cor[i]<-cor(pca$rotation[,i], as.numeric(as.factor(sex)))
  }
  top<-order(abs(pca.cor), decreasing = TRUE)[1]
  second<-order(abs(pca.cor), decreasing = TRUE)[2]
  print(paste("Top correlated principal components with sex:", top, ",", second))
  plot(pca$rotation[,top], pca$rotation[,second], pch = 16, col = c("magenta", "blue")[as.factor(sex)], xlab = paste("PC", top), ylab = paste("PC", second))
  legend("topright", levels(as.factor(sex)), pch = 16, col = c("magenta", "blue"))
  predSex<-rep(NA, length(sex))
  options.sex<-levels(as.factor(sex))
  if(abs(pca.cor[top]) > 0.9){
    print("Top PC has r > 0.9 with sex so good enough to confirm reported sexes")
    if(sign(pca.cor[top]) == 1){
      
      predSex[which(pca$rotation[,top] < 0)]<-options.sex[1]
      predSex[which(pca$rotation[,top] > 0)]<-options.sex[2]
    } else {
      
      predSex[which(pca$rotation[,top] < 0)]<-options.sex[2]
      predSex[which(pca$rotation[,top] > 0)]<-options.sex[1]
    }
  }
  return(predSex)
}


clusterGender<-function(betas, sex, npcs = 20, thres = 0.5){
  if(npcs > ncol(betas)){
    npcs<-ncol(betas)
    print(paste("As only", ncol(betas), "samples, can look look at that number of PCs", sep = " "))
  }
  betas.com<-betas[complete.cases(betas),]
  pca<-prcomp(betas.com)
  
  pca.cor<-rep(NA, npcs)
  for(i in 1:npcs){
    pca.cor[i]<-cor(pca$rotation[,i], as.numeric(as.factor(sex)))
  }
  predSex<-kmeans(pca$rotation[,which(abs(pca.cor) > thres)],2)$cluster
  return(predSex)
}

#runs the two gender predicting functions
pdf("PredictGender.pdf")
predSex<-findGenderPC(betas(Methylumiset), pheno$Gender)
dev.off()
BAD.samples.sexMissmatch<-as.character(pheno[which(predSex!=pheno$Gender),]$SentrixAdress)
############################################################################
############################################################################
############################################################################
#Uses SNP probes on the array to find genetic correlations between samples
betas<-betas(Methylumiset)
betas.rs<-betas[grep("rs", rownames(betas)),]
snpCor<-cor(betas.rs, use="complete.obs")

#ignores correlations between a sample and itself
for(i in 1:ncol(betas.rs)){
  snpCor[i,i]<-NA
}
#finds the maximum correlation values between samples and plots a histogram
#note: correlations of ~0.5 do not mean siblings as there are too few SNPs for this
corMax<-apply(snpCor, 1, max, na.rm = TRUE)
pdf("SNPCorrelation.MS.LC.pdf")
hist(corMax, xlab = "Max. correlation with all other samples", main = "")
dev.off()
genetic.cor<-data.frame(corMax)
rownames(genetic.cor)<-colnames(betas.rs)
#Visscher, P.M., et al. PLoS Genet, 2006. 2(3): p. e41.
############################################################################
pheno[which(genetic.cor[,1] > 0.65),]
genetic.cor[which(genetic.cor[,1] > 0.65),]
BAD.genetic.relationship<-rownames(genetic.cor)[which(genetic.cor[,1] > 0.65)]

############################################################################
############################################################################
############################################################################
bad<-c(BAD.samples.Bisulifite,BAD.BckgIntensityDiffPval,BAD.BckgIntensityDiff,BAD.samples.Intensity,BAD.genetic.relationship, BAD.samples.sexMissmatch)
bad<-unique(bad)
save(bad,BAD.samples.Bisulifite,BAD.BckgIntensityDiffPval,BAD.BckgIntensityDiff,BAD.samples.Intensity,BAD.genetic.relationship, BAD.samples.sexMissmatch,
     file="bad.samples.Rdata")
############################################################################
dat.pf<-pfilter(Methylumiset[,-which(colnames(Methylumiset)%in%bad)], perc=5)
#0 samples having 5 % of sites with a detection p-value greater than 0.05 were removed 
#Samples removed:  
#   125 sites were removed as beadcount <3 in 5 % of samples 
#   2897 sites having 1 % of samples with a detection p-value greater than 0.05 were removed 
dat.dasen<-dasen(dat.pf)
betas.dasen<-betas(dat.dasen)
pheno<-pheno[-which(colnames(Methylumiset)%in%bad),]
identical(colnames(betas.dasen), as.character(pheno$SentrixAdress))
#[1] TRUE

tiff("density.betas.tiff")
plot(density(betas.dasen[,1],na.rm=T),main="Distribution of Beta Values", xlim=c(0,1), col=c("black"), lty=1)
for (n in 2:ncol(betas.dasen)){
  lines(density(betas.dasen[,n],na.rm=T))}
dev.off()

cross<-read.csv("/mnt/data1/450K_reference/CrossHybridisingProbesPriceORWeksberg.csv")
crosslist<-cross[,2]
length(crosslist)
#[1]43233  

d<-betas.dasen[ ! rownames(betas.dasen) %in% crosslist, ]
dim(d) # [1] 440480      4   

#Made a file from the Weksburg_CommonSNP_EuAf_within10bpSBE column and Weksburg_CommonSNP_EuAf_inProbeSeq in the SNPs inprobeANno.csv. Called this snpprobes.csv. It contains 
Af_SNP<-read.csv("/mnt/data1/450K_reference/snp_chyb_Af_out.csv")
dim(Af_SNP)
# [1] 53794     1
SNPlist<-as.character(Af_SNP[,2]) 
d1<- d[ ! rownames(d) %in% SNPlist, ]
dim(d1) # [1] 425721    141

EU_SNP<-read.csv("/mnt/data1/450K_reference/snp_chyb_EU_out.csv")
dim(EU_SNP)
# [1] 53794     1
SNPlist<-as.character(EU_SNP[,2])
d2<- d1[ ! rownames(d1) %in% SNPlist, ]
dim(d2)
#[1] 424958    141
betas<-data.frame(d2)
identical(substr(colnames(betas),2,18),as.character(pheno$SentrixAdress))
colnames(betas)<-as.character(pheno$SentrixAdress)
setwd("/mnt/data1/Ehsan/meta_CrossCortical/NORMALISED_BETA")
save(betas,pheno, file="SAARLAND.TC.51.beta.pheno.Rdata")
#TRUE

resid <- function( row, age, sex, cet ){
  
  fit <- try (
    lm( row ~ age + sex + cet),
    silent=TRUE
  )
  if(inherits(fit,'try-error')) return(rep(NA, length(sex)))
  fit$residuals
} 

cl<- makeCluster(32)

pheno<-pheno[!is.na(pheno$Braak),]
betas<-betas[,as.character(pheno$SentrixAdress)]


colData(rgSet)@listData$sentrixFull <- as.character(colData(rgSet)@listData$sentrixFull)
library(minfi)


library('FlowSorted.DLPFC.450k')
cet.saar.tc <- estimateCellCounts(rgSet[,colnames(betas)], 
                                  compositeCellType = "DLPFC",
                                  processMethod = "auto", 
                                  probeSelect = "auto",
                                  cellTypes = c("NeuN_neg", "NeuN_pos"),
                                  referencePlatform = NULL,
                                  returnAll = FALSE, 
                                  meanPlot = TRUE, 
                                  verbose = TRUE)
cellprop<-(cet.saar.tc[,2]/(cet.saar.tc[,2]+cet.saar.tc[,1]))*100
pheno<-cbind(pheno, cellprop)
pheno$Braak<-pheno$Braak


dat.reg<- {              
  sex    	<- 	as.factor(pheno$Gender)
  age    	<- 	pheno$Age 
  cet		<- as.numeric(pheno$cellprop)
  t(parApply(cl, betas, 1, resid, age, sex, cet))
}
head(dat.reg[1:5,1:5])
identical(colnames(dat.reg), as.character(pheno$SentrixAdress))
stopCluster(cl)

setwd("/mnt/data1/Ehsan/meta_CrossCortical/REGRESSED_BETA")
save(dat.reg,pheno, file="SAAR.TC.51.pheno.REGRESSED.Rdata")

mod = model.matrix( ~ as.numeric(Braak),data=pheno)
mod0 = model.matrix(~1,data=pheno)

mydf<-na.omit(dat.reg)

sva<- sva::sva(dat = as.matrix(mydf), mod = mod, mod0 = mod0)

setwd("/mnt/data1/Ehsan/meta_CrossCortical/SVA/")
save(sva, file="sva.regressed.SAAR.tc.Rdata")
load("sva.regressed.SAAR.tc.Rdata")

svs<-sva$sv
pheno.SAAR.tc<-cbind(pheno, svs)

setwd("/mnt/data1/Ehsan/meta_CrossCortical/pheno/")
write.csv(pheno.SAAR.tc,file="pheno.SAAR.tc.51.csv")










































