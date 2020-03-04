rm(list=ls())

###------------------------------
###loading required libraries
###------------------------------
library(minfi)
library(parallel)
library(wateRmelon)
library(corrplot)
library(sva)

###------------------------------
###reading the demographics and phenotipic aoutcomes
###------------------------------
setwd("")
pheno<-read.csv(".csv")

###------------------------------
###Loading the raw idat files to creat a MethylumiSet object 
###------------------------------
idatPath <- c("")
barcodes<-as.character(pheno$Sentix_IDs)
methylumiset<- readEPIC(idatPath=idatPath, barcodes=barcodes, parallel = TRUE)# Construct methylumiSet from IDAT files
identical(colnames(methylumiset), as.character(pheno$Sentix_IDs))

###------------------------------
###Loading the raw idat files to creat a RGChannelSet object 
###------------------------------
targets <- data.frame(barcodes)
colnames(targets) <- "sentrixFull"
targets$Basename <- with(targets, file.path(idatPath, sentrixFull))
rgSet<- read.metharray.exp(targets = targets, verbose = TRUE, extended = FALSE)# Construct RGset from IDAT files

###-------------------------------------------------------------------------------------------------------------------
### QC.a : Detecting samples with extreme intesities using the negative control probes as references
###-------------------------------------------------------------------------------------------------------------------
con.red<-getRed(rgSet)[getProbeInfo(rgSet, type = "Control")$Address,]
rownames(con.red)<-getProbeInfo(rgSet, type = "Control")$ExtendedType
con.green<-getGreen(rgSet)[getProbeInfo(rgSet, type = "Control")$Address,]
rownames(con.green)<-getProbeInfo(rgSet, type = "Control")$ExtendedType

### Identify just the negative control probes
neg.red<-con.red[which(getProbeInfo(rgSet, type = "Control")$Type == "NEGATIVE"),]
neg.green<-con.green[which(getProbeInfo(rgSet, type = "Control")$Type == "NEGATIVE"),]

neg.red.mean<-colMeans(neg.red)
neg.green.mean<-colMeans(neg.green)
green.mean<-colMeans(getGreen(rgSet))
red.mean<-colMeans(getRed(rgSet))

par(mfrow=c(1,2))
plot(neg.red.mean, xlab = "Sample", ylab = "Intensity-red", col="grey", pch = 16, ylim = c(0, 12000))
par(new=TRUE)
plot(red.mean, xlab = "Sample", ylab = "", col="red", pch = 16, ylim = c(0, 12000))
plot(neg.green.mean, xlab = "Sample", ylab = "Intensity-green", col="grey", pch = 16, ylim = c(0, 12000))
par(new=TRUE)
plot(green.mean, xlab = "Sample", ylab = "", col="Green", pch = 16, ylim = c(0, 12000))

###Detecting the bad samples 
BAD.BckgIntensityDiff<-c(rownames(data.frame(which(neg.red.mean >1000))),rownames(data.frame(which(red.mean < 2000))))

###-------------------------------------------------------------------------------------------------------------------
### QC.b : Detecting samples with low background to signal ratio
###-------------------------------------------------------------------------------------------------------------------
avgPval<- colMeans(pvals(methylumiset))
barplot(avgPval,xlab=NULL)

###Detecting the bad samples 
BAD.BckgIntensityDiffPval<-as.character(colnames(methylumiset[,which(avgPval> 0.005)]))

###-------------------------------------------------------------------------------------------------------------------
### QC.c : detecting the samples with the mean intensity of methylated or unmethylated signals were three standard deviations above or below the mean
###-------------------------------------------------------------------------------------------------------------------
meth      <- methylated(methylumiset)
unmeth   <- unmethylated(methylumiset)

M.mean<-apply(meth, 2, mean)
U.mean<-apply(unmeth, 2, mean)

Mean.M<-mean(M.mean)
Mean.U<-mean(U.mean)
sd.M<-sd(M.mean)
sd.U<-sd(U.mean)
###Detecting the bad samples 
BAD.samples.Intensity<-as.character(pheno[which(Mean.M +3*sd.M < M.mean |
                                                  M.mean < Mean.M-3*sd.M | 
                                                  Mean.U +3*sd.U < U.mean |
                                                  U.mean< Mean.U-3*sd.U),]$sentrixFull)
												  
		par(mfrow=c(1,3))

#plots of average sample intensites and threshold lines to see if signals were sufficiently high
Sex<-as.factor(pheno$Gender)
plot(M.mean, U.mean, pch = 16, xlab = "Mean M intensity", ylab = "Mean U intensity", col = rainbow(nlevels(Sex))[factor(Sex)], main="BS Signal Intensities Coloured by Sex")
legend("topleft", legend=levels(factor(Sex)), col = rainbow(nlevels(Sex)), pch = 16, cex=0.6)

Position<-as.factor(substr(pheno$X,12,18))
plot(M.mean, U.mean, pch = 16, xlab = "Mean M intensity", ylab = "Mean U intensity", col = rainbow(nlevels(Position))[factor(Position)], main="BS Signal Intensities Coloured by Position")
legend("topleft", legend=levels(factor(Position)), col = rainbow(nlevels(Position)), pch = 16, cex=0.6)

Chip<-as.factor(substr(pheno$X,1,11))
plot(M.mean, U.mean, pch = 16, xlab = "Mean M intensity", ylab = "Mean U intensity", col = rainbow(nlevels(Chip))[factor(Chip)], main="BS Signal Intensities Coloured by Chip")
legend("topleft", legend=levels(factor(Chip)), col = rainbow(nlevels(Chip)), pch = 16, cex=0.6)										  

###-------------------------------------------------------------------------------------------------------------------
### QC.d : Detecting the samples with bisulfite conversion efficiency < 80%
###-------------------------------------------------------------------------------------------------------------------
#calculate bisulphite conversion statistic using bscon (from wateRmelon) and plots histogram
bs<-bscon(methylumiset)
pdf("HistogramBisulphiteConversionStatistics.pdf")
hist(bs, xlab = "Mean % BS conversion", main = "")
dev.off()
###Detecting the bad samples 
BAD.samples.Bisulifite<-as.character(colnames(Methylumiset[,which(bs< 80)]))										  

###-------------------------------------------------------------------------------------------------------------------
### QC.e : Detecting a mismatch between reported and predicted sex
###-------------------------------------------------------------------------------------------------------------------
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
predSex<-findGenderPC(betas(Methylumiset), pheno$Gender)

###Detecting the bad samples 
BAD.samples.sexMissmatch<-as.character(pheno[which(predSex!=pheno$Gender),]$sentrixFull)

###-------------------------------------------------------------------------------------------------------------------
### QC.f : Uses SNP probes on the array to find genetic correlations between samples
###-------------------------------------------------------------------------------------------------------------------
betas<-betas(methylumiset)
betas.rs<-betas[grep("rs", rownames(betas)),]
snpCor<-cor(betas.rs, use="complete.obs")

#ignores correlations between a sample and itself
for(i in 1:ncol(betas.rs)){
  snpCor[i,i]<-NA
}

#finds the maximum correlation values between samples and plots a histogram
#note: correlations of ~0.5 do not mean siblings as there are too few SNPs for this
corMax<-apply(snpCor, 1, max, na.rm = TRUE)
hist(corMax, xlab = "Max. correlation with all other samples", main = "")
dev.off()
genetic.cor<-data.frame(corMax)
rownames(genetic.cor)<-colnames(betas.rs)
#Visscher, P.M., et al. PLoS Genet, 2006. 2(3): p. e41.

###Detecting the bad samples 
BAD.genetic.relationship<-data.frame(pheno[which(genetic.cor[,1] > 0.65),c("Age","Gender","Braak")], genetic.cor[which(genetic.cor[,1] > 0.65),])$sentrixFull

###-------------------------------------------------------------------------------------------------------------------
###Detecting the unique bad samples identified at the different QC steps
###-------------------------------------------------------------------------------------------------------------------
bad<-unique((c(BAD.samples.Bisulifite,BAD.BckgIntensityDiffPval,BAD.BckgIntensityDiff,BAD.samples.Intensity,BAD.genetic.relationship, BAD.samples.sexMissmatch))

###-------------------------------------------------------------------------------------------------------------------
###pfilter function within the wateRmelon package29, with the following criteria used for exclusion: samples with a detection  P > 0.05 in more than 5% of probes, 
###probes with > three beadcount in 5% of samples and  probes having 1% of samples with a detection P value > 0.05).
###-------------------------------------------------------------------------------------------------------------------
dat.pf<-pfilter(methylumiset[,-which(colnames(Methylumiset)%in%bad)], perc=5)

###-------------------------------------------------------------------------------------------------------------------
###Quantile normalization was applied using the dasen function in the wateRmelon package and generaeting beta values
###-------------------------------------------------------------------------------------------------------------------
dat.dasen<-dasen(dat.pf)
betas.dasen<-betas(dat.dasen)
pheno<-pheno[-which(colnames(Methylumiset)%in%bad),]
identical(colnames(betas.dasen), as.character(pheno$sentrixFull))


###-------------------------------------------------------------------------------------------------------------------
###Removing the probes with evidence for cross-hybridizition
###-------------------------------------------------------------------------------------------------------------------
cross<-read.csv("")
crosslist<-cross[,2]
df<-betas.dasen[ ! rownames(betas.dasen) %in% crosslist, ]
SNPlist<-as.character(Af_SNP[,2]) 
df1<- d[ ! rownames(df) %in% SNPlist, ]

EU_SNP<-read.csv("/mnt/data1/450K_reference/snp_chyb_EU_out.csv")
SNPlist<-as.character(EU_SNP[,2])
df2<- df1[ ! rownames(df1) %in% SNPlist, ]

betas<-data.frame(df2)
identical((colnames(betas),as.character(pheno$sentrixFull))

###-------------------------------------------------------------------------------------------------------------------
### CET
###-------------------------------------------------------------------------------------------------------------------

load("CETS_Image.RData")

refProfile <- getReference(brain, modelIdx)
head(refProfile)
prop <- estProportion(dilution, profile = refProfile)

prop<- as.data.frame(estProportion(betas.dasen, profile = refProfile))
colnames(prop)<-"cellprop"
pheno<- merge(pheno, prop, by="row.names")
rownames(pheno)<-[,1]

###-------------------------------------------------------------------------------------------------------------------
### harmonisation: regressing out the effects of age and sex in all samples in each cohort and tissue separately, with neuron/glia proportions included as an additional covariate in cortical regions
###-------------------------------------------------------------------------------------------------------------------
resid <- function( row, age, sex, cet ){
  
  fit <- try (
    lm( row ~ age + sex + cet),
    silent=TRUE
  )
  if(inherits(fit,'try-error')) return(rep(NA, length(sex)))
  fit$residuals
} 

cl<- makeCluster()
pheno<-pheno[!is.na(pheno$Braak),]
betas<-betas[,!is.na(pheno$Braak)]
dat.reg<- {              
  sex    	<- 	as.factor(pheno$Gender)
  age    	<- 	pheno$Age 
  cet		<- as.numeric(pheno$cellprop)
  t(parApply(cl, betas, 1, resid, age, sex, cet))
}
identical(colnames(dat.reg), as.character(pheno$sentrixFull))
###-------------------------------------------------------------------------------------------------------------------
### Surrogate variables (SVs) were calculated using the sva function in the SVA package
###-------------------------------------------------------------------------------------------------------------------
mod = model.matrix( ~ as.numeric(Braak),data=pheno)
mod0 = model.matrix(~1,data=pheno)
mydf<-na.omit(dat.reg)

sva<- sva::sva(dat = as.matrix(mydf), mod = mod, mod0 = mod0)
############################################################################

