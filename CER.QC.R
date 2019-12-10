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
load("LBB2__RgSet_Methylumi_Pheno_CER95.Rdata")
###################################
#######LBB2.CER
setwd("/mnt/data1/Ehsan/meta_CrossCortical/CER/LBB2/QC")
###### QC.1#######
rgSet<-rgSet
Methylumiset<-LBB2.CER
pheno<-pheno.LBB2.CER

### Identify control probes

con.red<-getRed(rgSet)[getProbeInfo(rgSet, type = "Control")$Address,]
rownames(con.red)<-getProbeInfo(rgSet, type = "Control")$ExtendedType
con.green<-getGreen(rgSet)[getProbeInfo(rgSet, type = "Control")$Address,]
rownames(con.green)<-getProbeInfo(rgSet, type = "Control")$ExtendedType

### Identify just the negative control probes
neg.red.bs<-con.red[which(getProbeInfo(rgSet, type = "Control")$Type == "NEGATIVE"),]
neg.green.bs<-con.green[which(getProbeInfo(rgSet, type = "Control")$Type == "NEGATIVE"),]

neg.red.mean.bs<-colMeans(neg.red.bs)
neg.green.mean.bs<-colMeans(neg.green.bs)
green.mean.bs<-colMeans(getGreen(rgSet))
red.mean.bs<-colMeans(getRed(rgSet))

pdf("NegativeCtrlbckg.pdf", height=8, width=8)
par(mfrow=c(1,2))
plot(neg.red.mean.bs, xlab = "Sample", ylab = "Intensity-red", col="grey", pch = 16, ylim = c(0, 12000))
par(new=TRUE)
plot(red.mean.bs, xlab = "Sample", ylab = "", col="red", pch = 16, ylim = c(0, 12000))

plot(neg.green.mean.bs, xlab = "Sample", ylab = "Intensity-green", col="grey", pch = 16, ylim = c(0, 12000))
par(new=TRUE)
plot(green.mean.bs, xlab = "Sample", ylab = "", col="Green", pch = 16, ylim = c(0, 12000))
dev.off()

###################################
###################################
###################################
BAD.BckgIntensityDiff<-c(rownames(data.frame(which(neg.red.mean.bs >1000))),rownames(data.frame(which(red.mean.bs < 2000))))
###################################
###################################
###################################


avgPval<- colMeans(pvals(Methylumiset))

pdf("BckgIntensityDiffPval.pdf")
barplot(avgPval,xlab=NULL)
dev.off()

BAD.BckgIntensityDiffPval<-as.character(colnames(Methylumiset[,which(avgPval> 0.005)]))
identical(colnames(Methylumiset), as.character(pheno$X))
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
                                                  U.mean< Mean.U-3*sd.U),]$X)

pdf("IntensityByFactor.pdf", height=6, width=12)
par(mfrow=c(1,3))

#plots of median sample intensites and threshold lines to see if signals were sufficiently high
Sex<-as.factor(pheno$Gender)
plot(M.mean, U.mean, pch = 16, xlab = "Mean M intensity", ylab = "Mean U intensity", col = rainbow(nlevels(Sex))[factor(Sex)], main="BS Signal Intensities Coloured by Sex")
legend("topleft", legend=levels(factor(Sex)), col = rainbow(nlevels(Sex)), pch = 16, cex=0.6)

Position<-as.factor(substr(pheno$X,12,18))
plot(M.mean, U.mean, pch = 16, xlab = "Mean M intensity", ylab = "Mean U intensity", col = rainbow(nlevels(Position))[factor(Position)], main="BS Signal Intensities Coloured by Position")
legend("topleft", legend=levels(factor(Position)), col = rainbow(nlevels(Position)), pch = 16, cex=0.6)

Chip<-as.factor(substr(pheno$X,1,11))
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
BAD.samples.sexMissmatch<-as.character(pheno[which(predSex!=pheno$Gender),]$X)
#if you look at the plot you see that this is a false detection
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
data.frame(pheno[which(genetic.cor[,1] > 0.65),], genetic.cor[which(genetic.cor[,1] > 0.65),])
#                    Age        Gender Braak                                   %
#5815381006_R06C02 85.49760      1     3                                     0.9953073
#5815381017_R04C02 87.84668      2     5                                     0.6673011
#5822020002_R02C02 84.30390      1     3                                     0.6673011
#5822062011_R02C02 79.10472      1     1                                     0.9953073
#seems there are a twin pair or a replicate that we keep one awith shorter PMI
BAD.genetic.relationship<-c("3998567094_R03C01")
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
#   106 sites were removed as beadcount <3 in 5 % of samples 
#   4195 sites having 1 % of samples with a detection p-value greater than 0.05 were removed 
dat.dasen<-dasen(dat.pf)
betas.dasen<-betas(dat.dasen)
pheno<-pheno[-which(colnames(Methylumiset)%in%bad),]
identical(colnames(betas.dasen), as.character(pheno$X))
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
identical(substr(colnames(betas),2,18),as.character(pheno$X))
colnames(betas)<-as.character(pheno$X)
setwd("/mnt/data1/Ehsan/meta_CrossCortical/NORMALISED_BETA")
save(betas,pheno, file="LBB2.CER.93.beta.pheno.Rdata")
#TRUE
load("LBB2.CER.93.beta.pheno.Rdata")


resid <- function( row, age, sex){
  
  fit <- try (
    lm( row ~ age + sex),
    silent=TRUE
  )
  if(inherits(fit,'try-error')) return(rep(NA, length(sex)))
  fit$residuals
} 

cl<- makeCluster(32)


pheno<-pheno[!is.na(pheno$Braak_Stage),]
betas<-data.frame(betas[,as.character(pheno$X)])

dat.reg<- {              
  sex    	<- 	as.factor(pheno$Gender)
  age    	<- 	pheno$Age 
  t(parApply(cl, betas, 1, resid, age, sex))
}


head(dat.reg[1:5,1:5])

identical(colnames(dat.reg), as.character(pheno$X))

setwd("/mnt/data1/Ehsan/meta_CrossCortical/REGRESSED_BETA")
save(dat.reg,pheno, file="LBB2.CER.92.pheno.REGRESSED.Rdata")



mod = model.matrix( ~ as.numeric(Braak_Stage),data=pheno)
mod0 = model.matrix(~1,data=pheno)

mydf<-na.omit(dat.reg)

sva<- sva::sva(dat = as.matrix(mydf), mod = mod, mod0 = mod0)
setwd("/mnt/data1/Ehsan/meta_CrossCortical/SVA/")
save(sva, file="sva.regressed.lbb2.cer.Rdata")

svs<-sva$sv
pheno.LBB2.CER<-cbind(pheno, svs)

setwd("/mnt/data1/Ehsan/meta_CrossCortical/pheno/")
write.csv(pheno.LBB2.CER,file="pheno.LBB2.CER.92.csv")






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
load("LBB1__RgSet_Methylumi_Pheno_CER112.Rdata")
###################################
#######LBB2.CER
setwd("/mnt/data1/Ehsan/meta_CrossCortical/CER/LBB1/QC")
###### QC.1#######
rgSet<-rgSet
Methylumiset<-LBB1.CER
pheno<-pheno.LBB1.CER
pheno$X<- as.character(pheno$X.1)
### Identify control probes

con.red<-getRed(rgSet)[getProbeInfo(rgSet, type = "Control")$Address,]
rownames(con.red)<-getProbeInfo(rgSet, type = "Control")$ExtendedType
con.green<-getGreen(rgSet)[getProbeInfo(rgSet, type = "Control")$Address,]
rownames(con.green)<-getProbeInfo(rgSet, type = "Control")$ExtendedType

### Identify just the negative control probes
neg.red.bs<-con.red[which(getProbeInfo(rgSet, type = "Control")$Type == "NEGATIVE"),]
neg.green.bs<-con.green[which(getProbeInfo(rgSet, type = "Control")$Type == "NEGATIVE"),]

neg.red.mean.bs<-colMeans(neg.red.bs)
neg.green.mean.bs<-colMeans(neg.green.bs)
green.mean.bs<-colMeans(getGreen(rgSet))
red.mean.bs<-colMeans(getRed(rgSet))

pdf("NegativeCtrlbckg.pdf", height=8, width=8)
par(mfrow=c(1,2))
plot(neg.red.mean.bs, xlab = "Sample", ylab = "Intensity-red", col="grey", pch = 16, ylim = c(0, 12000))
par(new=TRUE)
plot(red.mean.bs, xlab = "Sample", ylab = "", col="red", pch = 16, ylim = c(0, 12000))

plot(neg.green.mean.bs, xlab = "Sample", ylab = "Intensity-green", col="grey", pch = 16, ylim = c(0, 12000))
par(new=TRUE)
plot(green.mean.bs, xlab = "Sample", ylab = "", col="Green", pch = 16, ylim = c(0, 12000))
dev.off()

###################################
###################################
###################################
BAD.BckgIntensityDiff<-c(rownames(data.frame(which(neg.red.mean.bs >1000))),rownames(data.frame(which(red.mean.bs < 2000))))
###################################
###################################
###################################


avgPval<- colMeans(pvals(Methylumiset))

pdf("BckgIntensityDiffPval.pdf")
barplot(avgPval,xlab=NULL)
dev.off()

BAD.BckgIntensityDiffPval<-as.character(colnames(Methylumiset[,which(avgPval> 0.005)]))
identical(colnames(Methylumiset), as.character(pheno$X.1))
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
                                                  U.mean< Mean.U-3*sd.U),]$X.1)

pdf("IntensityByFactor.pdf", height=6, width=12)
par(mfrow=c(1,3))

#plots of median sample intensites and threshold lines to see if signals were sufficiently high
Sex<-as.factor(pheno$Gender)
plot(M.mean, U.mean, pch = 16, xlab = "Mean M intensity", ylab = "Mean U intensity", col = rainbow(nlevels(Sex))[factor(Sex)], main="BS Signal Intensities Coloured by Sex")
legend("topleft", legend=levels(factor(Sex)), col = rainbow(nlevels(Sex)), pch = 16, cex=0.6)

Position<-as.factor(substr(pheno$X.1,12,18))
plot(M.mean, U.mean, pch = 16, xlab = "Mean M intensity", ylab = "Mean U intensity", col = rainbow(nlevels(Position))[factor(Position)], main="BS Signal Intensities Coloured by Position")
legend("topleft", legend=levels(factor(Position)), col = rainbow(nlevels(Position)), pch = 16, cex=0.6)

Chip<-as.factor(substr(pheno$X.1,1,11))
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
BAD.samples.sexMissmatch<-as.character(pheno[which(predSex!=pheno$Gender),]$X)
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
data.frame(pheno[which(genetic.cor[,1] > 0.65),], genetic.cor[which(genetic.cor[,1] > 0.65),])
#                    Age        Gender Braak                                   %
#5815381006_R06C02 85.49760      1     3                                     0.9953073
#5815381017_R04C02 87.84668      2     5                                     0.6673011
#5822020002_R02C02 84.30390      1     3                                     0.6673011
#5822062011_R02C02 79.10472      1     1                                     0.9953073
#seems there are a twin pair or a replicate that we keep one awith shorter PMI
BAD.genetic.relationship<-c("6042316127_R06C02")
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
#   106 sites were removed as beadcount <3 in 5 % of samples 
#   4195 sites having 1 % of samples with a detection p-value greater than 0.05 were removed 
dat.dasen<-dasen(dat.pf)
betas.dasen<-betas(dat.dasen)
pheno<-pheno[-which(colnames(Methylumiset)%in%bad),]
identical(colnames(betas.dasen), as.character(pheno$X))
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
identical(substr(colnames(betas),2,18),as.character(pheno$X))
colnames(betas)<-as.character(pheno$X)
setwd("/mnt/data1/Ehsan/meta_CrossCortical/NORMALISED_BETA")
save(betas,pheno, file="LBB1.CER.111.beta.pheno.Rdata")
#TRUE

resid <- function( row, age, sex){
  
  fit <- try (
    lm( row ~ age + sex),
    silent=TRUE
  )
  if(inherits(fit,'try-error')) return(rep(NA, length(sex)))
  fit$residuals
} 




library(parallel)
load("LBB1.CER.111.beta.pheno.Rdata")
pheno.pfc<-read.csv("/mnt/data1/Bex/Meta/SVAregressed/Allphenosub.csv")
pheno.LBB1<-pheno.pfc[which(pheno.pfc$Cohort=="LBB1" & pheno.pfc$region=="PFC"),]
pheno$sample<-as.character(pheno$Brain.Bank.ID)
pheno.pfc$sample<-as.character(pheno.pfc$sample)

pheno<-merge(pheno.LBB1,pheno,by="sample")

pheno<-pheno[!is.na(pheno$Braak),]
betas<-betas[,as.character(pheno$X.1)]
identical(colnames(betas),as.character(pheno$X.1))


cl<- makeCluster(32)

dat.reg<- {              
  sex    	<- 	as.factor(pheno$Gender2)
  age    	<- 	pheno$Age 
  t(parApply(cl, betas, 1, resid, age, sex))
}


head(dat.reg[1:5,1:5])
identical(colnames(dat.reg), as.character(pheno$X.1))

setwd("/mnt/data1/Ehsan/meta_CrossCortical/REGRESSED_BETA")
save(dat.reg,pheno, file="LBB1.CER.104.pheno.REGRESSED.Rdata")



mod = model.matrix( ~ as.numeric(Braak),data=pheno)
mod0 = model.matrix(~1,data=pheno)

mydf<-na.omit(dat.reg)

sva<- sva::sva(dat = as.matrix(mydf), mod = mod, mod0 = mod0)
setwd("/mnt/data1/Ehsan/meta_CrossCortical/SVA/")
save(sva, file="sva.regressed.lbb1.cer.Rdata")
load("sva.regressed.lbb1.cer.Rdata")
svs<-sva$sv
pheno.LBB1.CER<-cbind(pheno, svs)

setwd("/mnt/data1/Ehsan/meta_CrossCortical/pheno/")
write.csv(pheno.LBB1.CER,file="pheno.LBB1.CER.104.csv")






























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
load("AZ__RgSet_Methylumi_Pheno_CER404.Rdata")
###################################
#######LBB2.CER
setwd("/mnt/data1/Ehsan/meta_CrossCortical/CER/AZ/QC")
###### QC.1#######
rgSet<-rgSet
Methylumiset<-AZ.CER
pheno<-pheno.AZ.CER
pheno$X<- as.character(pheno$Complete.Barcode)
### Identify control probes

con.red<-getRed(rgSet)[getProbeInfo(rgSet, type = "Control")$Address,]
rownames(con.red)<-getProbeInfo(rgSet, type = "Control")$ExtendedType
con.green<-getGreen(rgSet)[getProbeInfo(rgSet, type = "Control")$Address,]
rownames(con.green)<-getProbeInfo(rgSet, type = "Control")$ExtendedType

### Identify just the negative control probes
neg.red.bs<-con.red[which(getProbeInfo(rgSet, type = "Control")$Type == "NEGATIVE"),]
neg.green.bs<-con.green[which(getProbeInfo(rgSet, type = "Control")$Type == "NEGATIVE"),]

neg.red.mean.bs<-colMeans(neg.red.bs)
neg.green.mean.bs<-colMeans(neg.green.bs)
green.mean.bs<-colMeans(getGreen(rgSet))
red.mean.bs<-colMeans(getRed(rgSet))

pdf("NegativeCtrlbckg.pdf", height=8, width=8)
par(mfrow=c(1,2))
plot(neg.red.mean.bs, xlab = "Sample", ylab = "Intensity-red", col="grey", pch = 16, ylim = c(0, 12000))
par(new=TRUE)
plot(red.mean.bs, xlab = "Sample", ylab = "", col="red", pch = 16, ylim = c(0, 12000))

plot(neg.green.mean.bs, xlab = "Sample", ylab = "Intensity-green", col="grey", pch = 16, ylim = c(0, 12000))
par(new=TRUE)
plot(green.mean.bs, xlab = "Sample", ylab = "", col="Green", pch = 16, ylim = c(0, 12000))
dev.off()

###################################
###################################
###################################
BAD.BckgIntensityDiff<-c(rownames(data.frame(which(neg.red.mean.bs >1000))),rownames(data.frame(which(red.mean.bs < 2000))))
###################################
###################################
###################################


avgPval<- colMeans(pvals(Methylumiset))

pdf("BckgIntensityDiffPval.pdf")
barplot(avgPval,xlab=NULL)
dev.off()

BAD.BckgIntensityDiffPval<-as.character(colnames(Methylumiset[,which(avgPval> 0.005)]))
identical(colnames(Methylumiset), as.character(pheno$X))
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
                                                  U.mean< Mean.U-3*sd.U),]$X)

pdf("IntensityByFactor.pdf", height=6, width=12)
par(mfrow=c(1,3))

#plots of median sample intensites and threshold lines to see if signals were sufficiently high
Sex<-as.factor(pheno$Gender2)
plot(M.mean, U.mean, pch = 16, xlab = "Mean M intensity", ylab = "Mean U intensity", col = rainbow(nlevels(Sex))[factor(Sex)], main="BS Signal Intensities Coloured by Sex")
legend("topleft", legend=levels(factor(Sex)), col = rainbow(nlevels(Sex)), pch = 16, cex=0.6)

Position<-as.factor(substr(pheno$X,12,18))
plot(M.mean, U.mean, pch = 16, xlab = "Mean M intensity", ylab = "Mean U intensity", col = rainbow(nlevels(Position))[factor(Position)], main="BS Signal Intensities Coloured by Position")
legend("topleft", legend=levels(factor(Position)), col = rainbow(nlevels(Position)), pch = 16, cex=0.6)

Chip<-as.factor(substr(pheno$X,1,11))
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
BAD.samples.sexMissmatch<-as.character(pheno[which(predSex!=pheno$Gender),]$X)
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
data.frame(pheno[which(genetic.cor[,1] > 0.65),], genetic.cor[which(genetic.cor[,1] > 0.65),])
#                    Age        Gender Braak                                   %
#5815381006_R06C02 85.49760      1     3                                     0.9953073
#5815381017_R04C02 87.84668      2     5                                     0.6673011
#5822020002_R02C02 84.30390      1     3                                     0.6673011
#5822062011_R02C02 79.10472      1     1                                     0.9953073
#seems there are a twin pair or a replicate that we keep one awith shorter PMI
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
#   792 sites were removed as beadcount <3 in 5 % of samples 
#   8741 sites having 1 % of samples with a detection p-value greater than 0.05 were removed 
dat.dasen<-dasen(dat.pf)
betas.dasen<-betas(dat.dasen)
pheno<-pheno[-which(colnames(Methylumiset)%in%bad),]
identical(colnames(betas.dasen), as.character(pheno$X))
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
identical(substr(colnames(betas),2,18),as.character(pheno$X))
colnames(betas)<-as.character(pheno$X)
setwd("/mnt/data1/Ehsan/meta_CrossCortical/NORMALISED_BETA")
save(betas,pheno, file="AZ.CER.340.beta.pheno.Rdata")
#TRUE

resid <- function( row, age, sex){
  
  fit <- try (
    lm( row ~ age + sex),
    silent=TRUE
  )
  if(inherits(fit,'try-error')) return(rep(NA, length(sex)))
  fit$residuals
} 




library(parallel)

pheno<-pheno[!is.na(pheno$Braak.score.1),]
betas<-betas[,as.character(pheno$X)]
identical(colnames(betas),as.character(pheno$X))


cl<- makeCluster(32)

dat.reg<- {              
  sex    	<- 	as.factor(pheno$Gender2)
  age    	<- 	pheno$Age 
  t(parApply(cl, betas, 1, resid, age, sex))
}


head(dat.reg[1:5,1:5])
identical(colnames(dat.reg), as.character(pheno$X))

setwd("/mnt/data1/Ehsan/meta_CrossCortical/REGRESSED_BETA")
save(dat.reg,pheno, file="AZ.CER.337.pheno.REGRESSED.Rdata")



mod = model.matrix( ~ as.numeric(Braak.score.1),data=pheno)
mod0 = model.matrix(~1,data=pheno)

mydf<-na.omit(dat.reg)

sva<- sva::sva(dat = as.matrix(mydf), mod = mod, mod0 = mod0)

setwd("/mnt/data1/Ehsan/meta_CrossCortical/SVA/")
save(sva, file="sva.regressed.az.cer.Rdata")

svs<-sva$sv
pheno.AZ.CER<-cbind(pheno, svs)

setwd("/mnt/data1/Ehsan/meta_CrossCortical/pheno/")
write.csv(pheno.AZ.CER,file="pheno.AZ.CER.337.csv")


