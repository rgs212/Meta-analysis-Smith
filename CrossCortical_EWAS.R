
rm(list=ls())
library(metafor)
library(parallel)

setwd("/mnt/data1/Ehsan/meta_CrossCortical/REGRESSED_BETA/")
load("LBB1.PFC.109.pheno.REGRESSED.Rdata")
lbb1.pfc<-dat.reg
rm(dat.reg, pheno)

load("LBB1.EC.101.pheno.REGRESSED.Rdata")
lbb1.ec<-dat.reg
rm(dat.reg, pheno)

load("LBB1.STG.112.pheno.REGRESSED.Rdata")
lbb1.stg<-dat.reg
rm(dat.reg, pheno)

load("LBB2.EC.87.pheno.REGRESSED.Rdata")
lbb2.ec<-dat.reg
rm(dat.reg, pheno)

load("MS.PFC.141.pheno.REGRESSED.Rdata")
ms.pfc<-dat.reg
rm(dat.reg, pheno)

load("MS.STG.143.pheno.REGRESSED.Rdata")
ms.stg<-dat.reg
rm(dat.reg, pheno)

load("AZ1.MTG.265.pheno.REGRESSED.Rdata")
az.mtg.1<-dat.reg
rm(dat.reg, pheno)

load("AZ2.MTG.88.pheno.REGRESSED.Rdata")
az.mtg.2<-dat.reg
rm(dat.reg, pheno)

load("ROSMAP.711.pheno.REGRESSED.Rdata")
rosmap.pfc<-dat.reg
rm(dat.reg, pheno)
dim(rosmap.pfc)

r1<-rownames(lbb1.pfc)
r2<-rownames(lbb1.ec)
r3<-rownames(lbb1.stg)
r4<-rownames(lbb2.ec)
r5<-rownames(ms.pfc)
r6<-rownames(ms.stg)
r7<-rownames(az.mtg.1)
r8<-rownames(az.mtg.2)
r9<-rownames(rosmap.pfc)


x<-intersect(intersect(intersect(intersect(intersect(intersect(intersect(intersect(r1,r2),r3),r4),r5),r6),r7),r8),r9)
length(x)
##[1] 412999
rm(r1,r2,r3,r4,r5,r6,r7,r8,r9)

lbb1.pfc  <-	lbb1.pfc	[x,]
lbb1.ec	  <-	lbb1.ec	[x,]
lbb1.stg  <-	lbb1.stg	[x,]
lbb2.ec	  <-	lbb2.ec	[x,]
ms.pfc	  <-	ms.pfc	[x,]
ms.stg	  <-	ms.stg	[x,]
az.mtg.1  <-	az.mtg.1	[x,]
az.mtg.2  <-	az.mtg.2	[x,]
rosmap.pfc<-	rosmap.pfc	[x,]

rm(x)

setwd("/mnt/data1/Ehsan/meta_CrossCortical/pheno/")
pheno.lbb1.pfc<-read.csv("pheno.LBB1.PFC.109.csv")
pheno.lbb1.ec<-read.csv("pheno.LBB1.EC.101.csv")
pheno.lbb1.stg<-read.csv("pheno.LBB1.STG.112.csv")
pheno.lbb2.ec<-read.csv("pheno.LBB2.EC.88.csv")
pheno.ms.pfc<-read.csv("pheno.MS.PFC.141.csv")
pheno.ms.stg<-read.csv("pheno.MS.STG.143.csv")
pheno.az.mtg.1<-read.csv("pheno.AZ1.265.csv")
pheno.az.mtg.2<-read.csv("pheno.AZ1.88.csv")
pheno.rosmap.pfc<-read.csv("pheno.ROSMAP.711.csv")


EWAS <- function(x,Braak){
  res<-lm(x ~ Braak)
  return(coef(summary(res))[2,])
} 

cl<- makeCluster(32)

Braak<-pheno.lbb1.pfc$Braak
res.lbb1.pfc<-t(parApply(cl, lbb1.pfc, 1, EWAS,Braak))
head(res.lbb1.pfc)
chisq <- qchisq(1-(res.lbb1.pfc[,4]),1)
median(chisq)/qchisq(0.5,1)
#[1] 1.21789


Braak<-pheno.lbb1.ec$Braak
res.lbb1.ec<-t(parApply(cl, lbb1.ec, 1, EWAS,Braak))
head(res.lbb1.ec)
chisq <- qchisq(1-(res.lbb1.ec[,4]),1)
median(chisq)/qchisq(0.5,1)
#[1] 1.045444

Braak<-pheno.lbb1.stg$Braak
res.lbb1.stg<-t(parApply(cl, lbb1.stg, 1, EWAS,Braak))
head(res.lbb1.stg)
chisq <- qchisq(1-(res.lbb1.stg[,4]),1)
median(chisq)/qchisq(0.5,1)
#[1] 0.9434773

Braak<-pheno.lbb2.ec$Braak
res.lbb2.ec<-t(parApply(cl, lbb2.ec, 1, EWAS,Braak))
head(res.lbb2.ec)
chisq <- qchisq(1-(res.lbb2.ec[,4]),1)
median(chisq)/qchisq(0.5,1)
# [1] 0.9482733


Braak<-pheno.ms.pfc$Braak
res.ms.pfc<-t(parApply(cl, ms.pfc, 1, EWAS,Braak))
head(res.ms.pfc)
chisq <- qchisq(1-(res.ms.pfc[,4]),1)
median(chisq)/qchisq(0.5,1)
# [1] 1.02681


Braak<-pheno.ms.stg$Braak
res.ms.stg<-t(parApply(cl, ms.stg, 1, EWAS,Braak))
head(res.ms.stg)
chisq <- qchisq(1-(res.ms.stg[,4]),1)
median(chisq)/qchisq(0.5,1)
# 1.08183


Braak<-pheno.az.mtg.1$Braak
res.az.mtg.1<-t(parApply(cl, az.mtg.1, 1, EWAS,Braak))
head(res.az.mtg.1)
chisq <- qchisq(1-(res.az.mtg.1[,4]),1)
median(chisq)/qchisq(0.5,1)
# [1] 1.457281


Braak<-pheno.az.mtg.2$Braak
res.az.mtg.2<-t(parApply(cl, az.mtg.2, 1, EWAS,Braak))
head(res.az.mtg.2)
chisq <- qchisq(1-(res.az.mtg.2[,4]),1)
median(chisq)/qchisq(0.5,1)
# [1] 1.161352


Braak<-pheno.rosmap.pfc$Braak
res.rosmap.pfc<-t(parApply(cl, rosmap.pfc, 1, EWAS,Braak))
head(res.rosmap.pfc)
chisq <- qchisq(1-(res.rosmap.pfc[,4]),1)
median(chisq)/qchisq(0.5,1)
#[1] 1.349907



res.ewas<-data.frame(res.lbb1.pfc,
                res.lbb1.ec,
                res.lbb1.stg,
                res.lbb2.ec,
                res.ms.pfc,
                res.ms.stg,
                res.az.mtg.1,
                res.az.mtg.2,
                res.rosmap.pfc)
setwd("/mnt/data1/Ehsan/meta_CrossCortical/")
save(res.ewas, file="res.ewas.noSV.rda")
#===========================================================================
EWAS <- function(x,Braak, sv1,sv2,sv3,sv4,sv5){
  res<-lm(x ~ Braak + sv1 + sv2 + sv3 + sv4 + sv5)
  return(coef(summary(res))[2,])
} 

cl<- makeCluster(32)
Braak<-pheno.lbb1.pfc$Braak
sv1<-pheno.lbb1.pfc$X1
sv2<-pheno.lbb1.pfc$X2
sv3<-pheno.lbb1.pfc$X3
sv4<-pheno.lbb1.pfc$X4
sv5<-pheno.lbb1.pfc$X5

res.lbb1.pfc<-t(parApply(cl, lbb1.pfc, 1, EWAS,Braak, sv1,sv2,sv3,sv4,sv5))
head(res.lbb1.pfc)
chisq <- qchisq(1-(res.lbb1.pfc[,4]),1)
median(chisq)/qchisq(0.5,1)
#[1] 1.083204
#===========================================================================
EWAS <- function(x){
  res<-lm(x ~ Braak)
  return(coef(summary(res))[2,])
} 

Braak<-pheno.lbb1.ec$Braak
res.lbb1.ec<-t(parApply(cl, lbb1.ec, 1, EWAS,Braak))
head(res.lbb1.ec)
chisq <- qchisq(1-(res.lbb1.ec[,4]),1)
median(chisq)/qchisq(0.5,1)
#[1] 1.045444
#===========================================================================
EWAS <- function(x,Braak, sv1,sv2,sv3,sv4,sv5){
  res<-lm(x ~ Braak + sv1 + sv2 + sv3 + sv4 + sv5)
  return(coef(summary(res))[2,])
} 
Braak<-pheno.lbb1.stg$Braak
sv1<-pheno.lbb1.stg$X1
sv2<-pheno.lbb1.stg$X2
sv3<-pheno.lbb1.stg$X3
sv4<-pheno.lbb1.stg$X4
sv5<-pheno.lbb1.stg$X5

res.lbb1.stg<-t(parApply(cl, lbb1.stg, 1, EWAS,Braak, sv1,sv2,sv3,sv4,sv5))
head(res.lbb1.stg)
chisq <- qchisq(1-(res.lbb1.stg[,4]),1)
median(chisq)/qchisq(0.5,1)
#[1] 1.007998
#===========================================================================
EWAS <- function(x){
  res<-lm(x ~ Braak)
  return(coef(summary(res))[2,])
} 
Braak<-pheno.ms.pfc$Braak
res.ms.pfc<-t(parApply(cl, ms.pfc, 1, EWAS,Braak))
head(res.ms.pfc)
chisq <- qchisq(1-(res.ms.pfc[,4]),1)
median(chisq)/qchisq(0.5,1)
# [1] 1.02681
#===========================================================================
EWAS <- function(x){
  res<-lm(x ~ Braak)
  return(coef(summary(res))[2,])
} 
Braak<-pheno.ms.stg$Braak
res.ms.stg<-t(parApply(cl, ms.stg, 1, EWAS,Braak))
head(res.ms.stg)
chisq <- qchisq(1-(res.ms.stg[,4]),1)
median(chisq)/qchisq(0.5,1)
# 1.08183
#===========================================================================
EWAS <- function(x,Braak, sv1,sv2){
  res<-lm(x ~ Braak + sv1 + sv2)
  return(coef(summary(res))[2,])
} 
Braak<-pheno.lbb2.ec$Braak
sv1<-pheno.lbb2.ec$X1
sv2<-pheno.lbb2.ec$X2

res.lbb2<-t(parApply(cl, lbb2.ec, 1, EWAS,Braak, sv1,sv2))
head(res.lbb2)
chisq <- qchisq(1-(res.lbb2[,4]),1)
median(chisq)/qchisq(0.5,1)
#[1] 1.077476
save(res.lbb2, file="/mnt/data1/Ehsan/meta_CrossCortical/LBB2/res.lbb2.pfc.2sv.Rdata")

#=============================================================================
#=============================================================================
#=============================================================================

EWAS <- function(x,Braak, sv1,sv2,sv3,sv4,sv5,sv6,sv7,sv8,sv9,sv10){
  res<-lm(x ~ Braak + sv1 + sv2 + sv3 + sv4 + sv5 + sv6 + sv7 +sv8 +sv9 + sv10)
  return(coef(summary(res))[2,])
} 

sv1<-pheno.az.mtg.1$X1
sv2<-pheno.az.mtg.1$X2
sv3<-pheno.az.mtg.1$X3
sv4<-pheno.az.mtg.1$X4
sv5<-pheno.az.mtg.1$X5
sv6<-pheno.az.mtg.1$X6
sv7<-pheno.az.mtg.1$X7
sv8<-pheno.az.mtg.1$X8
sv9<-pheno.az.mtg.1$X9
sv10<-pheno.az.mtg.1$X10

Braak<-pheno.az.mtg.1$Braak
res.az.mtg.1<-t(parApply(cl, az.mtg.1, 1, EWAS,Braak, sv1,sv2,sv3,sv4,sv5,sv6,sv7,sv8,sv9,sv10))
head(res.az.mtg.1)
chisq <- qchisq(1-(res.az.mtg.1[,4]),1)
median(chisq)/qchisq(0.5,1)
# [1] 0.994314
save(res.az.mtg.1, file="/mnt/data1/Ehsan/meta_CrossCortical/AZ1/res.az1.mtg.Rdata")

#===========================================================================
#===========================================================================
EWAS <- function(x,Braak, sv1,sv2,sv3,sv4,sv5,sv6,sv7,sv8){
  res<-lm(x ~ Braak + sv1 + sv2 + sv3 + sv4 + sv5 + sv6 + sv7 +sv8 )
  return(coef(summary(res))[2,])
} 


sv1<-pheno.az.mtg.2$X1
sv2<-pheno.az.mtg.2$X2
sv3<-pheno.az.mtg.2$X3
sv4<-pheno.az.mtg.2$X4
sv5<-pheno.az.mtg.2$X5
sv6<-pheno.az.mtg.2$X6
sv7<-pheno.az.mtg.2$X7
sv8<-pheno.az.mtg.2$X8


Braak<-pheno.az.mtg.2$Braak
res.az.mtg.2<-t(parApply(cl, az.mtg.2, 1, EWAS,Braak, sv1,sv2,sv3,sv4,sv5,sv6,sv7,sv8))
head(res.az.mtg.2)
chisq <- qchisq(1-(res.az.mtg.2[,4]),1)
median(chisq)/qchisq(0.5,1)
#[1] 1.097746

save(res.az.mtg.2, file="/mnt/data1/Ehsan/meta_CrossCortical/AZ2/res.az2.mtg.Rdata")


#=============================================================================
#=============================================================================
#=============================================================================
#=============================================================================
EWAS <- function(x,Braak, sv1,sv2,sv3,sv4,sv5){
  res<-lm(x ~ Braak + sv1 + sv2 + sv3 + sv4 + sv5)
  return(coef(summary(res))[2,])
} 

sv1<-pheno.rosmap.pfc$X1
sv2<-pheno.rosmap.pfc$X2
sv3<-pheno.rosmap.pfc$X3
sv4<-pheno.rosmap.pfc$X4
sv5<-pheno.rosmap.pfc$X5

Braak<-pheno.rosmap.pfc$Braak
res.rosmap.pfc<-t(parApply(cl, rosmap.pfc, 1, EWAS,Braak, sv1,sv2,sv3,sv4,sv5))
head(res.rosmap.pfc)
chisq <- qchisq(1-(res.rosmap.pfc[,4]),1)
median(chisq)/qchisq(0.5,1)
#[1] 1.062303


save(res.rosmap.pfc,file="/mnt/data1/Ehsan/meta_CrossCortical/ROSMAP/res.rosmap.pfc.5sv.Rdata")
#=============================================================================
#=============================================================================
#=============================================================================
#=============================================================================
rm(list=ls())
library(metafor)
library(parallel)

setwd("/mnt/data1/Ehsan/meta_CrossCortical/REGRESSED_BETA/")
load("LBB1.PFC.109.pheno.REGRESSED.Rdata")
lbb1.pfc<-dat.reg
rm(dat.reg, pheno)

load("LBB1.EC.101.pheno.REGRESSED.Rdata")
lbb1.ec<-dat.reg
rm(dat.reg, pheno)

load("LBB1.STG.112.pheno.REGRESSED.Rdata")
lbb1.stg<-dat.reg
rm(dat.reg, pheno)

load("LBB2.EC.87.pheno.REGRESSED.Rdata")
lbb2.ec<-dat.reg
rm(dat.reg, pheno)

load("MS.PFC.141.pheno.REGRESSED.Rdata")
ms.pfc<-dat.reg
rm(dat.reg, pheno)

load("MS.STG.143.pheno.REGRESSED.Rdata")
ms.stg<-dat.reg
rm(dat.reg, pheno)

load("AZ1.MTG.265.pheno.REGRESSED.Rdata")
az.mtg.1<-dat.reg
rm(dat.reg, pheno)

load("AZ2.MTG.88.pheno.REGRESSED.Rdata")
az.mtg.2<-dat.reg
rm(dat.reg, pheno)

load("ROSMAP.709.pheno.REGRESSED.Rdata")
rosmap.pfc<-dat.reg
rm(dat.reg, pheno)

r1<-rownames(lbb1.pfc)
r2<-rownames(lbb1.ec)
r3<-rownames(lbb1.stg)
r4<-rownames(lbb2.ec)
r5<-rownames(ms.pfc)
r6<-rownames(ms.stg)
r7<-rownames(az.mtg.1)
r8<-rownames(az.mtg.2)
r9<-rownames(rosmap.pfc)


x<-intersect(intersect(intersect(intersect(intersect(intersect(intersect(intersect(r1,r2),r3),r4),r5),r6),r7),r8),r9)
length(x)
##[1] 413003
rm(r1,r2,r3,r4,r5,r6,r7,r8,r9)

lbb1.pfc  <-	lbb1.pfc	[x,]
lbb1.ec	  <-	lbb1.ec	[x,]
lbb1.stg  <-	lbb1.stg	[x,]
lbb2.ec	  <-	lbb2.ec	[x,]
ms.pfc	  <-	ms.pfc	[x,]
ms.stg	  <-	ms.stg	[x,]
az.mtg.1  <-	az.mtg.1	[x,]
az.mtg.2  <-	az.mtg.2	[x,]
rosmap.pfc<-	rosmap.pfc	[x,]



rm(x)

setwd("/mnt/data1/Ehsan/meta_CrossCortical/pheno/")
pheno.lbb1.pfc<-read.csv("pheno.LBB1.PFC.109.csv")
pheno.lbb1.ec<-read.csv("pheno.LBB1.EC.101.csv")
pheno.lbb1.stg<-read.csv("pheno.LBB1.STG.112.csv")
pheno.lbb2.ec<-read.csv("pheno.LBB2.EC.88.csv")
pheno.ms.pfc<-read.csv("pheno.MS.PFC.141.csv")
pheno.ms.stg<-read.csv("pheno.MS.STG.143.csv")
pheno.az.mtg.1<-read.csv("pheno.AZ1.265.csv")
pheno.az.mtg.2<-read.csv("pheno.AZ1.88.csv")
pheno.rosmap.pfc<-read.csv("pheno.ROSMAP.709.csv")




rownames(pheno.lbb1.pfc)<-as.character(pheno.lbb1.pfc$sample)
rownames(pheno.lbb1.ec)<-as.character(pheno.lbb1.ec$sample)
rownames(pheno.lbb1.stg)<-as.character(pheno.lbb1.stg$sample)

colnames(lbb1.pfc)<-rownames(pheno.lbb1.pfc)
colnames(lbb1.ec)<-rownames(pheno.lbb1.ec)
colnames(lbb1.stg)<-rownames(pheno.lbb1.stg)

y<-intersect(rownames(lbb1.pfc),intersect(rownames(lbb1.ec),rownames(lbb1.stg)))


length(y)
#[1] 413003
dat.pfc  <-	lbb1.pfc	[y,]
dat.ec	  <-	lbb1.ec	[y,]
dat.stg  <-	lbb1.stg	[y,]

pheno.pfc<-pheno.lbb1.pfc[,1:20]
pheno.ec<-pheno.lbb1.ec[,1:20]
pheno.stg<-pheno.lbb1.stg[,1:20]



identical(colnames(dat.pfc),rownames(pheno.pfc))
#[1] TRUE
identical(colnames(dat.ec),rownames(pheno.ec))
#[1] TRUE
identical(colnames(dat.stg),rownames(pheno.stg))
#[1] TRUE


dat.LBB1<-cbind(dat.pfc,dat.ec,dat.stg)
pheno.LBB1<-rbind(pheno.pfc, pheno.ec, pheno.stg )

identical(as.character(colnames(dat.LBB1)), as.character(pheno.LBB1$sample))



pheno<-pheno.LBB1[!is.na(pheno.LBB1$Braak),]
betas<-dat.LBB1[,as.character(pheno$sample)]

identical(as.character(colnames(betas)), as.character(pheno$sample))

Braak<-pheno$Braak
pheno$id<-as.character(pheno$sample)
id<-pheno$id
region<-as.character(pheno$region)
methyl<-as.numeric(betas["cg11823178",])
svs<-data.frame(pheno[, 11:15])


data<-data.frame(id,region,Braak,methyl,svs)

res<-nlme::lme(methyl ~ Braak, random=~1|id, data=data)
summary(res)

res<-nlme::lme(methyl ~ Braak + X1  + X2 + X3 + X4 + X5 , random=~1|id, data=data)
summary(res)

cl<- makeCluster(32)

meta <- function(x, Braak,id, sv1, sv2, sv3, sv4, sv5){
  library(nlme)
  res<-try(nlme::lme(as.numeric(x) ~ Braak + sv1  + sv2 + sv3 + sv4 + sv5 , random=~1|id), silent=T)
  if(inherits(res,'try-error')) return(rep(NA,5))
  return(coef(summary(res))[2,])
} 


sv1<-pheno$X1
sv2<-pheno$X2
sv3<-pheno$X3
sv4<-pheno$X4
sv5<-pheno$X5


res<-t(parApply(cl, betas, 1, meta, Braak,id, sv1, sv2, sv3, sv4, sv5))

chisq <- qchisq(1-(res[,5]),1)
median(chisq)/qchisq(0.5,1)

setwd("/mnt/data1/Ehsan/meta_CrossCortical/LBB1/CrossRegion")
save(res, file="lme.regressed.lbb1.322.rda")

#=============================================================================
#=============================================================================
#=============================================================================
#=============================================================================
########################################################################################################################
########################################################################################################################
########################################################################################################################


rm(list=ls())
library(metafor)
library(parallel)

setwd("/mnt/data1/Ehsan/meta_CrossCortical/REGRESSED_BETA/")
load("LBB1.PFC.109.pheno.REGRESSED.Rdata")
lbb1.pfc<-dat.reg
rm(dat.reg, pheno)

load("LBB1.EC.101.pheno.REGRESSED.Rdata")
lbb1.ec<-dat.reg
rm(dat.reg, pheno)

load("LBB1.STG.112.pheno.REGRESSED.Rdata")
lbb1.stg<-dat.reg
rm(dat.reg, pheno)

load("LBB2.EC.87.pheno.REGRESSED.Rdata")
lbb2.ec<-dat.reg
rm(dat.reg, pheno)

load("MS.PFC.141.pheno.REGRESSED.Rdata")
ms.pfc<-dat.reg
rm(dat.reg, pheno)

load("MS.STG.143.pheno.REGRESSED.Rdata")
ms.stg<-dat.reg
rm(dat.reg, pheno)

load("AZ1.MTG.265.pheno.REGRESSED.Rdata")
az.mtg.1<-dat.reg
rm(dat.reg, pheno)

load("AZ2.MTG.88.pheno.REGRESSED.Rdata")
az.mtg.2<-dat.reg
rm(dat.reg, pheno)

load("ROSMAP.709.pheno.REGRESSED.Rdata")
rosmap.pfc<-dat.reg
rm(dat.reg, pheno)

r1<-rownames(lbb1.pfc)
r2<-rownames(lbb1.ec)
r3<-rownames(lbb1.stg)
r4<-rownames(lbb2.ec)
r5<-rownames(ms.pfc)
r6<-rownames(ms.stg)
r7<-rownames(az.mtg.1)
r8<-rownames(az.mtg.2)
r9<-rownames(rosmap.pfc)


x<-intersect(intersect(intersect(intersect(intersect(intersect(intersect(intersect(r1,r2),r3),r4),r5),r6),r7),r8),r9)
length(x)
##[1] 413003
rm(r1,r2,r3,r4,r5,r6,r7,r8,r9)

lbb1.pfc  <-	lbb1.pfc	[x,]
lbb1.ec	  <-	lbb1.ec	[x,]
lbb1.stg  <-	lbb1.stg	[x,]
lbb2.ec	  <-	lbb2.ec	[x,]
ms.pfc	  <-	ms.pfc	[x,]
ms.stg	  <-	ms.stg	[x,]
az.mtg.1  <-	az.mtg.1	[x,]
az.mtg.2  <-	az.mtg.2	[x,]
rosmap.pfc<-	rosmap.pfc	[x,]



rm(x)

setwd("/mnt/data1/Ehsan/meta_CrossCortical/pheno/")
pheno.lbb1.pfc<-read.csv("pheno.LBB1.PFC.109.csv")
pheno.lbb1.ec<-read.csv("pheno.LBB1.EC.101.csv")
pheno.lbb1.stg<-read.csv("pheno.LBB1.STG.112.csv")
pheno.lbb2.ec<-read.csv("pheno.LBB2.EC.88.csv")
pheno.ms.pfc<-read.csv("pheno.MS.PFC.141.csv")
pheno.ms.stg<-read.csv("pheno.MS.STG.143.csv")
pheno.az.mtg.1<-read.csv("pheno.AZ1.265.csv")
pheno.az.mtg.2<-read.csv("pheno.AZ1.88.csv")
pheno.rosmap.pfc<-read.csv("pheno.ROSMAP.709.csv")




rownames(pheno.ms.pfc)<-as.character(pheno.ms.pfc$sample)
rownames(pheno.ms.stg)<-as.character(pheno.ms.stg$sample)

colnames(ms.pfc)<-rownames(pheno.ms.pfc)
colnames(ms.stg)<-rownames(pheno.ms.stg)

y<-intersect(rownames(ms.pfc),rownames(ms.stg))


length(y)
#[1] 413003
dat.pfc  <-	ms.pfc	[y,]
dat.stg  <-	ms.stg	[y,]



pheno.pfc<-pheno.ms.pfc[,1:20]
pheno.stg<-pheno.ms.stg[,1:20]


identical(colnames(dat.pfc),rownames(pheno.pfc))
#[1] TRUE
identical(colnames(dat.stg),rownames(pheno.stg))
#[1] TRUE


dat.MS<-cbind(dat.pfc,dat.stg)
pheno.MS<-rbind(pheno.pfc, pheno.stg )

identical(as.character(colnames(dat.MS)), as.character(pheno.MS$sample))



pheno<-pheno.MS[!is.na(pheno.MS$Braak),]
betas<-dat.MS[,as.character(pheno$sample)]

identical(as.character(colnames(betas)), as.character(pheno$sample))


Braak<-pheno$Braak
pheno$id<-as.character(pheno$sample)
id<-pheno$id
region<-as.character(pheno$region)
methyl<-as.numeric(betas["cg11823178",])
svs<-data.frame(pheno[, 11:15])


data<-data.frame(id,region,Braak,methyl,svs)

res<-nlme::lme(methyl ~ Braak, random=~1|id, data=data)
summary(res)

res<-nlme::lme(methyl ~ Braak + X1  + X2 + X3 + X4 + X5 , random=~1|id, data=data)
summary(res)

cl<- makeCluster(32)

meta <- function(x, Braak,id, sv1, sv2, sv3, sv4, sv5){
  library(nlme)
  res<-try(nlme::lme(as.numeric(x) ~ Braak + sv1  + sv2 + sv3 + sv4 + sv5 , random=~1|id), silent=T)
  if(inherits(res,'try-error')) return(rep(NA,5))
  return(coef(summary(res))[2,])
} 


sv1<-pheno$X1
sv2<-pheno$X2
sv3<-pheno$X3
sv4<-pheno$X4
sv5<-pheno$X5


res<-t(parApply(cl, betas, 1, meta, Braak,id, sv1, sv2, sv3, sv4, sv5))
head(res)

res.ms<-na.omit(res)
chisq <- qchisq(1-(res.ms[,5]),1)
median(chisq)/qchisq(0.5,1)

setwd("/mnt/data1/Ehsan/meta_CrossCortical/MS/CrossRegion")
save(res, file="lme.regressed.ms.284.rda")

##############################################################################
##############################################################################
setwd("/mnt/data1/Ehsan/meta_CrossCortical/REGRESSED_BETA/")

load("LBB1.CER.104.pheno.REGRESSED.Rdata")
lbb1.cer<-dat.reg
rm(dat.reg, pheno)

load("LBB2.CER.92.pheno.REGRESSED.Rdata")
lbb2.cer<-dat.reg
rm(dat.reg, pheno)

load("AZ.CER.337.pheno.REGRESSED.Rdata")
az.cer<-dat.reg
rm(dat.reg, pheno)

setwd("/mnt/data1/Ehsan/meta_CrossCortical/pheno/")
pheno.lbb1.cer<-read.csv("pheno.LBB1.CER.104.csv")
pheno.lbb2.cer<-read.csv("pheno.LBB2.CER.92.csv")
pheno.az.cer<-read.csv("pheno.AZ.CER.337.csv")


r1<-rownames(lbb1.cer)
r2<-rownames(lbb2.cer)
r3<-rownames(az.cer)

x<-intersect(intersect(r1,r2),r3)
#419,220

lbb1.cer  <-	lbb1.cer[x,]
lbb2.cer	<-	lbb2.cer[x,]
az.cer    <-  az.cer[x,]

EWAS <- function(x,Braak){
  res<-lm(x ~ Braak)
  return(coef(summary(res))[2,])
} 

cl<- makeCluster(32)

Braak<-pheno.lbb1.cer$Braak
res.lbb1.cer<-t(parApply(cl, lbb1.cer, 1, EWAS,Braak))
head(res.lbb1.cer)
chisq <- qchisq(1-(res.lbb1.cer[,4]),1)
median(chisq)/qchisq(0.5,1)
#[1] 1.357345


EWAS <- function(x,Braak, sv1,sv2,sv3,sv4,sv5){
  res<-lm(x ~ Braak + sv1 + sv2 + sv3 + sv4 + sv5)
  return(coef(summary(res))[2,])
} 

sv1<-pheno.lbb1.cer$X1
sv2<-pheno.lbb1.cer$X2
sv3<-pheno.lbb1.cer$X3
sv4<-pheno.lbb1.cer$X4
sv5<-pheno.lbb1.cer$X5

res.lbb1.cer<-t(parApply(cl, lbb1.cer, 1, EWAS,Braak, sv1,sv2,sv3,sv4,sv5))
head(res.lbb1.cer)
chisq <- qchisq(1-(res.lbb1.cer[,4]),1)
median(chisq)/qchisq(0.5,1)

#0.94221


Braak<-pheno.lbb2.cer$Braak

sv1<-pheno.lbb2.cer$X1
sv2<-pheno.lbb2.cer$X2
sv3<-pheno.lbb2.cer$X3
sv4<-pheno.lbb2.cer$X4
sv5<-pheno.lbb2.cer$X5

res.lbb2.cer<-t(parApply(cl, lbb2.cer, 1, EWAS,Braak, sv1,sv2,sv3,sv4,sv5))
head(res.lbb2.cer)
chisq <- qchisq(1-(res.lbb2.cer[,4]),1)
median(chisq)/qchisq(0.5,1)
#[1]1.055837

Braak<-pheno.az.cer$Braak.score.1

sv1<-pheno.az.cer$X1
sv2<-pheno.az.cer$X2
sv3<-pheno.az.cer$X3
sv4<-pheno.az.cer$X4
sv5<-pheno.az.cer$X5

res.az.cer<-t(parApply(cl, az.cer, 1, EWAS,Braak, sv1,sv2,sv3,sv4,sv5))
head(res.az.cer)
chisq <- qchisq(1-(res.az.cer[,4]),1)
median(chisq)/qchisq(0.5,1)

#[1] 0.9757663

setwd("/mnt/data1/Ehsan/meta_CrossCortical/CER/")
save(res.lbb1.cer, res.lbb2.cer, res.az.cer, file="res.cer.lbb1_2_5svs.rda")




rm(list=ls())
library(metafor)
library(parallel)


setwd("/mnt/data1/Ehsan/meta_CrossCortical/pheno/")
pheno.saar.pfc<-read.csv("pheno.SAAR.fc.45.csv")
pheno.saar.tc<-read.csv("pheno.SAAR.tc.51.csv")

pheno.saar.np<-read.csv("pheno.SAAR.NN.26.csv")
pheno.saar.nn<-read.csv("pheno.SAAR.NP.25.csv")


setwd("/mnt/data1/Ehsan/meta_CrossCortical/REGRESSED_BETA/")
load("SAAR.FC.45.pheno.REGRESSED.Rdata")
dat.reg.fc<-dat.reg
load("SAAR.TC.51.pheno.REGRESSED.Rdata")
dat.reg.tc<-dat.reg
rm(dat.reg, pheno)
identical(colnames(dat.reg.fc), as.character(pheno.saar.pfc$SentrixAdress))
identical(colnames(dat.reg.tc), as.character(pheno.saar.tc$SentrixAdress))


load("SAAR.NN26.NP25.pheno.REGRESSED.Rdata")

setwd("/mnt/data1/Ehsan/meta_CrossCortical/")
load("RESULTS_CROSS.CORTICAL.META.Rda")
res.meta.bonf<-results[which(results$bonf.FE<0.05),]
saar.tc.220<-data.frame(dat.reg.tc[rownames(res.meta.bonf),])
saar.pfc.220<-data.frame(dat.reg.fc[rownames(res.meta.bonf),])

saar.np.220<-data.frame(dat.reg.NN[rownames(res.meta.bonf),])
saar.nn.220<-data.frame(dat.reg.NP[rownames(res.meta.bonf),])

cl<- makeCluster(4)
identical(substr(colnames(saar.pfc.220),2,20), as.character(pheno.saar.pfc$SentrixAdress))


EWAS <- function(x,Braak){
  res<-lm(x ~ Braak)
  return(coef(summary(res))[2,])
} 
res.saar.pfc<-t(parApply(cl, saar.pfc.220, 1, EWAS,pheno.saar.pfc$Braak))
head(res.saar.pfc)
res.saar.tc<-t(parApply(cl, saar.tc.220, 1, EWAS,pheno.saar.tc$Braak))
head(res.saar.tc)
res.saar.nn<-t(parApply(cl, saar.nn.220, 1, EWAS,pheno.saar.nn$Braak))
head(res.saar.nn)
res.saar.np<-t(parApply(cl, saar.np.220, 1, EWAS,pheno.saar.np$Braak))
head(res.saar.np)

meta <- function(x, Braak,id){
  library(nlme)
  res<-try(nlme::lme(as.numeric(x) ~ Braak, random=~1|id), silent=T)
  if(inherits(res,'try-error')) return(rep(NA,5))
  return(coef(summary(res))[2,])
} 

pheno.tc<-pheno.saar.tc[, c(2,6)]
pheno.fc<-pheno.saar.pfc[, c(2,6)]
pheno<-rbind(pheno.fc,pheno.tc)
dat<-cbind(saar.pfc.220, saar.tc.220)

res.mixed.tc.fc<-t(parApply(cl, dat, 1, meta, pheno$BraakStaging,pheno$DonorID))
head(res.mixed.tc.fc)


res.220<-cbind(res.saar.pfc,res.saar.tc,res.saar.nn,res.saar.np,res.mixed.tc.fc)
colnames(res.220)<-c("MUNICH.FC.Estimate","MUNICH.FC.SE","MUNICH.FC.t","MUNICH.FC.P",
                     "MUNICH.TC.Estimate","MUNICH.TC.SE","MUNICH.TC.t","MUNICH.TC.P",
                     "MUNICH.Nneg.Estimate","MUNICH.Nneg.SE","MUNICH.Nneg.t","MUNICH.Nneg.P",
                     "MUNICH.Npos.Estimate","MUNICH.Npos.SE","MUNICH.Npos.t","MUNICH.Npos.P",
                     "MUNICH.cc.Estimate","MUNICH.cc.SE","DF","MUNICH.cc.t","MUNICH.cc.P")
                     

setwd("/mnt/data1/Ehsan/meta_CrossCortical/")
load("RESULTS_CROSS.CORTICAL.META.Rda")
res.cc<-results
res.all<-read.csv("/mnt/data1/Bex/Meta/New/alldataresultswithCB.csv")
rownames(res.cc)<-rownames(results)
head(res.cc)
rownames(res.all)<-res.all[,1]
x<-intersect(rownames(res.all), rownames(res.cc))
res.all.f<-res.all[x,]
res.cc.f<-res.cc[x,]
identical(rownames(res.all.f),rownames(res.cc.f))
res.all.f<-cbind(res.all.f,res.cc.f)
head(res.all.f)
resultsall220<-cbind(res.all.f[rownames(res.220),1:124],res.220)
write.table(resultsall220, file="results.220.txt", sep="\t")



