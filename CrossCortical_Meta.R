rm(list=ls())
library(metafor)
library(parallel)

load("/mnt/data1/Ehsan/meta_CrossCortical/LBB1/CrossRegion/lme.regressed.lbb1.322.rda")
res.lbb1<-res
rm(res)
load("/mnt/data1/Ehsan/meta_CrossCortical/LBB2/res.lbb2.pfc.2sv.Rdata")
load("/mnt/data1/Ehsan/meta_CrossCortical/MS/CrossRegion/lme.regressed.ms.284.rda")
res.ms<-res
rm(res)
load("/mnt/data1/Ehsan/meta_CrossCortical/AZ1/res.az1.mtg.Rdata")
res.az1<-res.az.mtg.1
rm(res.az.mtg.1)
load("/mnt/data1/Ehsan/meta_CrossCortical/AZ2/res.az2.mtg.Rdata")
res.az2<-res.az.mtg.2
rm(res.az.mtg.2)
load("/mnt/data1/Ehsan/meta_CrossCortical/ROSMAP/res.rosmap.pfc.5sv.Rdata")
res.rosmap<-res.rosmap.pfc
rm(res.rosmap.pfc)
                                
x<-intersect(intersect(intersect(intersect(intersect(rownames(res.az1), rownames(res.az2)),rownames(res.ms)),rownames(res.lbb2)),rownames(res.lbb1)),rownames(res.rosmap))


ANNOT450<-read.csv("/mnt/data1/450K_reference/AdditionalAnnotation450K_Price_Great_SNPsinProbeSequence.csv")
rownames(ANNOT450)<-ANNOT450[,1]
XY<-rownames(ANNOT[which(ANNOT$CHR=="X"|ANNOT$CHR=="Y"),])
length(intersect(x,XY))

x<-x[-which(x%in%XY)]

res.lbb1<-res.lbb1[x,]
res.lbb2<-res.lbb2[x,]
res.az1<-res.az1[x,]
res.az2<-res.az2[x,]
res.lbb2<-res.lbb2[x,]
res.ms<-res.ms[x,]
res.rosmap<-res.rosmap[x,]


es<-cbind(res.lbb1[,1],res.ms[,1],res.lbb2[,1],res.az1[,1],res.az2[,1],res.rosmap[,1])
se<-cbind(res.lbb1[,2],res.ms[,2],res.lbb2[,2],res.az1[,2],res.az2[,2],res.rosmap[,2])
es1<-na.omit(es)
se1<-na.omit(se)
identical(rownames(es), rownames(se))

library(BiocParallel)
register(MulticoreParam(32, log=TRUE))
colnames(es1) <- colnames(se1) <- c("LBB1","MS","LBB2","AZ1","AZ2","ROSMAP")
rownames(es1) <- rownames(se1) <-rownames(es1)
head(colnames(es1))

library(bacon)
bc <- bacon(NULL, es1, se1)


knitr::kable(estimates(bc))

inflation(bc)
#LBB1       MS     LBB2      AZ1      AZ2   ROSMAP 
#1.242364 1.085899 1.035447 0.996604 1.044030 1.027027 
bias(bc)

knitr::kable(es(bc)[1:5,])
knitr::kable(se(bc)[1:5,])


bcm <- bacon::meta(na.omit(bc))
head(pval(bcm))

setwd( "/mnt/data1/Ehsan/meta_CrossCortical/")
bmp("meta.cc.bmp", width=17,height=10, unit="cm", res=600)
plot(bcm, type="qq")
dev.off()

p<-data.frame(pval(bcm))

chisq <- qchisq(1-(as.numeric(p$meta)),1)
lambda = median(chisq)/qchisq(0.5,1)
#[1] 1.184592


es.b<-es(bc)
se.b<-se(bc)

identical(rownames(es.b), rownames(se.b))
data<-cbind(es.b,se.b)

save(data, file="CrosscorticalMetaanalysis.Rda")
load("CrosscorticalMetaanalysis.Rda")

out<-metagen(data[1,1:6],data[1,7:12])


for(i in 1:nrow(data)){
  out<-metagen(data[i,1:6],data[i,7:12])
  res.meta[i,]<-c(out$TE.fixed, out$seTE.fixed, out$pval.fixed, out$TE.random, out$seTE.random, out$pval.random, out$I2, 1-pchisq(out$Q, out$df.Q))
}

library(parallel)
cl<-parallel::makeCluster(32)
res.meta<-t(parApply(cl,data,1,meta))
colnames(res.meta)<-c("TE.fixed", "seTE.fixed", "pval.fixed", "TE.random", "seTE.random", "pval.random", "I2", "1-pchisq(Q, df.Q)")
head(res.meta)

results<-res.meta
chisq <- qchisq(1-(results[,3]),1)
median(chisq)/qchisq(0.5,1)
#[1] 1.184279
chisq <- qchisq(1-(results[,6]),1)
median(chisq)/qchisq(0.5,1)
#[1] 0.9858139

bonf.FE<-p.adjust(results[,3], method="bonferroni")
length(which(bonf.FE<0.05))
#[1] 220
fdr.FE<-p.adjust(results[,3], method="fdr")
length(which(fdr.FE<0.05))
#[1] 2408

bonf.RE<-p.adjust(results[,6], method="bonferroni")
length(which(bonf.RE<0.05))
#[1] 107
fdr.RE<-p.adjust(results[,6], method="fdr")
length(which(fdr.RE<0.05))
#[1] 1189

resultsfdr<-cbind(results, bonf.FE, bonf.RE, fdr.FE, fdr.RE)

annot<-ANNOT
head(annot)
rownames(annot)<-annot[,1]
x<-intersect(rownames(annot),rownames(resultsfdr))
annot.sub<-annot[x,]
resultsfdr<-resultsfdr[x,]

res.meta.bacon<-cbind(resultsfdr,annot.sub)
head(res.meta.bacon)

res.meta.bacon.ordered<-res.meta.bacon[order(res.meta.bacon$pval.fixed),]

res.meta.bonf<-res.meta.bacon.ordered[which(res.meta.bacon.ordered$bonf.FE < 0.05),]
write.table(res.meta.bonf, file="res.meta.bonf.txt", sep="\t")

results<-res.meta.bacon.ordered

save(results, file="RESULTS_CROSS.CORTICAL.META.Rda")



dmr <- data.frame("chrom" = 	paste0("chr", results$CHR), 
                  "start" = 	results[rownames(results), "MAPINFO.1"],
                  "end" 	= 	results[rownames(results), "MAPINFO.1.1"],
                  "pvalue" = 	results$pval.fixed)

dmr<-dmr[order(dmr[,1],dmr[,2]),]
head(dmr)
setwd("/mnt/data1/Ehsan/meta_CrossCortical/")
write.table(dmr, file="meta.cc.DMR.txt",sep = "\t", row.names = FALSE,col.names = TRUE, quote = FALSE)

#comb-p pipeline -c 4 --dist 500 --seed 1.0e-4 --anno hg19 -p  meta.cc.DMR meta.cc.DMR.txt
res.meta.bonf.f<-res.meta.bonf[-which(res.meta.bonf$UCSC_RefGene_Name==""),]
length(unique(res.meta.bonf.f$Closest_TSS_gene_name))
smithList<-as.character(unique(res.meta.bonf.f$Closest_TSS_gene_name))
save(smithList, file="CCBONFERRONIUNIQUEGENENAME.rda")
res.meta.f<-results[-which(results$UCSC_RefGene_Name==""),]
length(unique(res.meta.f$Closest_TSS_gene_name))
GeneList<-as.character(unique(res.meta.bonf.f$Closest_TSS_gene_name))






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

results<-res.all.f
mitogene<-read.table("mito.txt", sep="\t")
mitogene<-as.character(mitogene[,1])

results.mito<-results[which(results$gene %in% mitogene),]


results.mito.fc<-results.mito[which(results.mito$PFC.P_Fixed<0.05),]
write.table(results.mito.fc,file="results.mito.fc.txt", sep="\t")
results.mito.cc<-results.mito[which(results.mito$pval.fixed<0.05),]
write.table(results.mito.cc,file="results.mito.cc.txt", sep="\t")

length(intersect(rownames(results.mito.fc), rownames(results.mito.cc)))





res220all<-results[rownames(res.meta.bonf),]
Gene<-results[which(results$gene=="ADAMTS4"),]
ADAMTS4<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
           "1" ,"NA", Gene$pval.cc)

Gene<-results[which(results$gene=="CR1"),]
CR1<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
       ((sumlog(Gene$pval.cc))$df)/2,(sumlog(Gene$pval.cc))$chisq, 
       (sumlog(Gene$pval.cc))$p)

Gene<-results[which(results$gene=="BIN1"),]
BIN1<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
        ((sumlog(Gene$pval.cc))$df)/2,(sumlog(Gene$pval.cc))$chisq, 
        (sumlog(Gene$pval.cc))$p)


Gene<-results[which(results$gene=="INPPD5"),]
Gene<-results[which(results$gene=="HESX1"),]
HESX1<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
         ((sumlog(Gene$pval.cc))$df)/2,(sumlog(Gene$pval.cc))$chisq, 
         (sumlog(Gene$pval.cc))$p)



Gene<-results[which(results$gene=="CLNK"),]
CLNK<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
        ((sumlog(Gene$pval.cc))$df)/2,(sumlog(Gene$pval.cc))$chisq, 
        (sumlog(Gene$pval.cc))$p)

Gene<-results[which(results$gene=="HS3ST1"),]

HS3ST1<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
          ((sumlog(Gene$pval.cc))$df)/2,(sumlog(Gene$pval.cc))$chisq, 
          (sumlog(Gene$pval.cc))$p)

Gene<-results[which(results$gene=="HLA-DRB1"),]


Gene<-results[which(results$gene=="CD2AP"),]
CD2AP<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
         ((sumlog(Gene$pval.cc))$df)/2,(sumlog(Gene$pval.cc))$chisq, 
         (sumlog(Gene$pval.cc))$p)

Gene<-results[which(results$gene=="ZCWPW1"),]
ZCWPW1<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
          ((sumlog(Gene$pval.cc))$df)/2,(sumlog(Gene$pval.cc))$chisq, 
          (sumlog(Gene$pval.cc))$p)

Gene<-results[which(results$gene=="EPHA1"),]
EPHA1<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
         ((sumlog(Gene$pval.cc))$df)/2,(sumlog(Gene$pval.cc))$chisq, 
         (sumlog(Gene$pval.cc))$p)

Gene<-results[which(results$gene=="CNTNAP2"),]
CNTNAP2<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
           ((sumlog(Gene$pval.cc))$df)/2,(sumlog(Gene$pval.cc))$chisq, 
           (sumlog(Gene$pval.cc))$p)

Gene<-results[which(results$gene=="CLU"),]
CLU<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
       ((sumlog(Gene$pval.cc))$df)/2,(sumlog(Gene$pval.cc))$chisq, 
       (sumlog(Gene$pval.cc))$p)

Gene<-results[which(results$gene=="ECHDC3"),]
ECHDC3<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
          ((sumlog(Gene$pval.cc))$df)/2,(sumlog(Gene$pval.cc))$chisq, 
          (sumlog(Gene$pval.cc))$p)

Gene<-results[which(results$gene=="MS4A6A"),]
MS4A6A<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
          ((sumlog(Gene$pval.cc))$df)/2,(sumlog(Gene$pval.cc))$chisq, 
          (sumlog(Gene$pval.cc))$p)

Gene<-results[which(results$gene=="PICALM"),]
PICALM<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
          ((sumlog(Gene$pval.cc))$df)/2,(sumlog(Gene$pval.cc))$chisq, 
          (sumlog(Gene$pval.cc))$p)

Gene<-results[which(results$gene=="SORL1"),]
SORL1<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
         ((sumlog(Gene$pval.cc))$df)/2,(sumlog(Gene$pval.cc))$chisq, 
         (sumlog(Gene$pval.cc))$p)

Gene<-results[which(results$gene=="SLC24A4"),]
SLC24A4<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
           ((sumlog(Gene$pval.cc))$df)/2,(sumlog(Gene$pval.cc))$chisq, 
           (sumlog(Gene$pval.cc))$p)

Gene<-results[which(results$gene=="ADAM10"),]
ADAM10<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
          ((sumlog(Gene$pval.cc))$df)/2,(sumlog(Gene$pval.cc))$chisq, 
          (sumlog(Gene$pval.cc))$p)

Gene<-results[which(results$gene=="APH1B"),]
APH1B<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
         ((sumlog(Gene$pval.cc))$df)/2,(sumlog(Gene$pval.cc))$chisq, 
         (sumlog(Gene$pval.cc))$p)

Gene<-results[which(results$gene=="KAT8"),]


Gene<-results[which(results$gene=="ECHDC3"),]
ECHDC3<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
          ((sumlog(Gene$pval.cc))$df)/2,(sumlog(Gene$pval.cc))$chisq, 
          (sumlog(Gene$pval.cc))$p)

Gene<-results[which(results$gene=="ABI3"),]
ABI3<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
        ((sumlog(Gene$pval.cc))$df)/2,(sumlog(Gene$pval.cc))$chisq, 
        (sumlog(Gene$pval.cc))$p)

Gene<-results[which(results$gene=="SUZ12P1"),]

Gene<-results[which(results$gene=="ALPK2"),]
ALPK2<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
         ((sumlog(Gene$pval.cc))$df)/2,(sumlog(Gene$pval.cc))$chisq, 
         (sumlog(Gene$pval.cc))$p)

Gene<-results[which(results$gene=="ABCA7"),]
ABCA7<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
         ((sumlog(Gene$pval.cc))$df)/2,(sumlog(Gene$pval.cc))$chisq, 
         (sumlog(Gene$pval.cc))$p)

Gene<-results[which(results$gene=="APOE"),]
APOE<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
        ((sumlog(Gene$pval.cc))$df)/2,(sumlog(Gene$pval.cc))$chisq, 
        (sumlog(Gene$pval.cc))$p)

Gene<-results[which(results$gene=="AC074212.3"),]

Gene<-results[which(results$gene=="CD33"),]
CD33<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
        ((sumlog(Gene$pval.cc))$df)/2,(sumlog(Gene$pval.cc))$chisq, 
        (sumlog(Gene$pval.cc))$p)
Gene<-results[which(results$gene=="CASS4"),]
CASS4<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
         ((sumlog(Gene$pval.cc))$df)/2,(sumlog(Gene$pval.cc))$chisq, 
         (sumlog(Gene$pval.cc))$p)

Gene<-results[which(results$gene=="TREM2"),]
TREM2<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
         ((sumlog(Gene$pval.cc))$df)/2,(sumlog(Gene$pval.cc))$chisq, 
         (sumlog(Gene$pval.cc))$p)


res.GWAS<-rbind(ABCA7,ABI3,ADAM10,ADAMTS4,ALPK2,APH1B,APOE,BIN1,CASS4,CD2AP,CD33,CLNK,CNTNAP2,CR1,ECHDC3,EPHA1,HESX1,PICALM,
                TREM2,SLC24A4,SORL1,ZCWPW1,CLU,HS3ST1,MS4A6A)
res.GWAS<-data.frame(res.GWAS)
str(res.GWAS)

res.GWAS[,1]<-NULL
res.GWAS[,1]<-as.numeric(as.character(res.GWAS[,1]))
res.GWAS[,2]<-as.numeric(as.character(res.GWAS[,2]))
res.GWAS[,3]<-as.numeric(as.character(res.GWAS[,3]))
res.GWAS[,4]<-as.numeric(as.character(res.GWAS[,4]))
res.GWAS[,5]<-as.numeric(as.character(res.GWAS[,5]))
res.GWAS[,6]<-as.numeric(as.character(res.GWAS[,6]))
res.GWAS$p.BH<-p.adjust(res.GWAS[,6],"BH")
res.GWAS$p.bonf<-p.adjust(res.GWAS[,6],"bonferroni")
colnames(res.GWAS)[1:6]<-c("CHR", "start","end","nProbes","X2","p.Fisher´s test")
res.GWAS.CC<-res.GWAS
res.GWAS.CC

rm(res.GWAS, ABCA7,ABI3,ADAM10,ADAMTS4,ALPK2,APH1B,APOE,BIN1,CASS4,CD2AP,CD33,CLNK,CNTNAP2,CR1,ECHDC3,EPHA1,HESX1,PICALM,
   TREM2,SLC24A4,SORL1,ZCWPW1,CLU,HS3ST1,MS4A6A)




Gene<-results[which(results$gene=="ADAMTS4"),]
ADAMTS4<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
           "1" ,"NA", Gene$PFC.P_Fixed)

Gene<-results[which(results$gene=="CR1"),]
CR1<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
       ((sumlog(Gene$PFC.P_Fixed))$df)/2,(sumlog(Gene$PFC.P_Fixed))$chisq, 
       (sumlog(Gene$PFC.P_Fixed))$p)

Gene<-results[which(results$gene=="BIN1"),]
BIN1<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
        ((sumlog(Gene$PFC.P_Fixed))$df)/2,(sumlog(Gene$PFC.P_Fixed))$chisq, 
        (sumlog(Gene$PFC.P_Fixed))$p)


Gene<-results[which(results$gene=="INPPD5"),]
Gene<-results[which(results$gene=="HESX1"),]
HESX1<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
         ((sumlog(Gene$PFC.P_Fixed))$df)/2,(sumlog(Gene$PFC.P_Fixed))$chisq, 
         (sumlog(Gene$PFC.P_Fixed))$p)



Gene<-results[which(results$gene=="CLNK"),]
CLNK<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
        ((sumlog(Gene$PFC.P_Fixed))$df)/2,(sumlog(Gene$PFC.P_Fixed))$chisq, 
        (sumlog(Gene$PFC.P_Fixed))$p)

Gene<-results[which(results$gene=="HS3ST1"),]

HS3ST1<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
          ((sumlog(Gene$PFC.P_Fixed))$df)/2,(sumlog(Gene$PFC.P_Fixed))$chisq, 
          (sumlog(Gene$PFC.P_Fixed))$p)

Gene<-results[which(results$gene=="HLA-DRB1"),]


Gene<-results[which(results$gene=="CD2AP"),]
CD2AP<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
         ((sumlog(Gene$PFC.P_Fixed))$df)/2,(sumlog(Gene$PFC.P_Fixed))$chisq, 
         (sumlog(Gene$PFC.P_Fixed))$p)

Gene<-results[which(results$gene=="ZCWPW1"),]
ZCWPW1<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
          ((sumlog(Gene$PFC.P_Fixed))$df)/2,(sumlog(Gene$PFC.P_Fixed))$chisq, 
          (sumlog(Gene$PFC.P_Fixed))$p)

Gene<-results[which(results$gene=="EPHA1"),]
EPHA1<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
         ((sumlog(Gene$PFC.P_Fixed))$df)/2,(sumlog(Gene$PFC.P_Fixed))$chisq, 
         (sumlog(Gene$PFC.P_Fixed))$p)

Gene<-results[which(results$gene=="CNTNAP2"),]
CNTNAP2<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
           ((sumlog(Gene$PFC.P_Fixed))$df)/2,(sumlog(Gene$PFC.P_Fixed))$chisq, 
           (sumlog(Gene$PFC.P_Fixed))$p)

Gene<-results[which(results$gene=="CLU"),]
CLU<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
       ((sumlog(Gene$PFC.P_Fixed))$df)/2,(sumlog(Gene$PFC.P_Fixed))$chisq, 
       (sumlog(Gene$PFC.P_Fixed))$p)

Gene<-results[which(results$gene=="ECHDC3"),]
ECHDC3<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
           ((sumlog(Gene$PFC.P_Fixed))$df)/2,(sumlog(Gene$PFC.P_Fixed))$chisq, 
           (sumlog(Gene$PFC.P_Fixed))$p)

Gene<-results[which(results$gene=="MS4A6A"),]
MS4A6A<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
          ((sumlog(Gene$PFC.P_Fixed))$df)/2,(sumlog(Gene$PFC.P_Fixed))$chisq, 
          (sumlog(Gene$PFC.P_Fixed))$p)

Gene<-results[which(results$gene=="PICALM"),]
PICALM<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
          ((sumlog(Gene$PFC.P_Fixed))$df)/2,(sumlog(Gene$PFC.P_Fixed))$chisq, 
          (sumlog(Gene$PFC.P_Fixed))$p)

Gene<-results[which(results$gene=="SORL1"),]
SORL1<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
         ((sumlog(Gene$PFC.P_Fixed))$df)/2,(sumlog(Gene$PFC.P_Fixed))$chisq, 
         (sumlog(Gene$PFC.P_Fixed))$p)

Gene<-results[which(results$gene=="SLC24A4"),]
SLC24A4<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
           ((sumlog(Gene$PFC.P_Fixed))$df)/2,(sumlog(Gene$PFC.P_Fixed))$chisq, 
           (sumlog(Gene$PFC.P_Fixed))$p)

Gene<-results[which(results$gene=="ADAM10"),]
ADAM10<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
          ((sumlog(Gene$PFC.P_Fixed))$df)/2,(sumlog(Gene$PFC.P_Fixed))$chisq, 
          (sumlog(Gene$PFC.P_Fixed))$p)

Gene<-results[which(results$gene=="APH1B"),]
APH1B<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
         ((sumlog(Gene$PFC.P_Fixed))$df)/2,(sumlog(Gene$PFC.P_Fixed))$chisq, 
         (sumlog(Gene$PFC.P_Fixed))$p)

Gene<-results[which(results$gene=="KAT8"),]


Gene<-results[which(results$gene=="ECHDC3"),]
ECHDC3<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
           ((sumlog(Gene$PFC.P_Fixed))$df)/2,(sumlog(Gene$PFC.P_Fixed))$chisq, 
           (sumlog(Gene$PFC.P_Fixed))$p)

Gene<-results[which(results$gene=="ABI3"),]
ABI3<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
        ((sumlog(Gene$PFC.P_Fixed))$df)/2,(sumlog(Gene$PFC.P_Fixed))$chisq, 
        (sumlog(Gene$PFC.P_Fixed))$p)

Gene<-results[which(results$gene=="SUZ12P1"),]

Gene<-results[which(results$gene=="ALPK2"),]
ALPK2<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
         ((sumlog(Gene$PFC.P_Fixed))$df)/2,(sumlog(Gene$PFC.P_Fixed))$chisq, 
         (sumlog(Gene$PFC.P_Fixed))$p)

Gene<-results[which(results$gene=="ABCA7"),]
ABCA7<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
         ((sumlog(Gene$PFC.P_Fixed))$df)/2,(sumlog(Gene$PFC.P_Fixed))$chisq, 
         (sumlog(Gene$PFC.P_Fixed))$p)

Gene<-results[which(results$gene=="APOE"),]
APOE<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
        ((sumlog(Gene$PFC.P_Fixed))$df)/2,(sumlog(Gene$PFC.P_Fixed))$chisq, 
        (sumlog(Gene$PFC.P_Fixed))$p)

Gene<-results[which(results$gene=="AC074212.3"),]

Gene<-results[which(results$gene=="CD33"),]
CD33<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
        ((sumlog(Gene$PFC.P_Fixed))$df)/2,(sumlog(Gene$PFC.P_Fixed))$chisq, 
        (sumlog(Gene$PFC.P_Fixed))$p)
Gene<-results[which(results$gene=="CASS4"),]
CASS4<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
         ((sumlog(Gene$PFC.P_Fixed))$df)/2,(sumlog(Gene$PFC.P_Fixed))$chisq, 
         (sumlog(Gene$PFC.P_Fixed))$p)

Gene<-results[which(results$gene=="TREM2"),]
TREM2<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
         ((sumlog(Gene$PFC.P_Fixed))$df)/2,(sumlog(Gene$PFC.P_Fixed))$chisq, 
         (sumlog(Gene$PFC.P_Fixed))$p)


res.GWAS<-rbind(ABCA7,ABI3,ADAM10,ADAMTS4,ALPK2,APH1B,APOE,BIN1,CASS4,CD2AP,CD33,CLNK,CNTNAP2,CR1,ECHDC3,EPHA1,HESX1,PICALM,
                TREM2,SLC24A4,SORL1,ZCWPW1,CLU,HS3ST1,MS4A6A)


res.GWAS<-data.frame(res.GWAS)
str(res.GWAS)

res.GWAS[,1]<-NULL
res.GWAS[,1]<-as.numeric(as.character(res.GWAS[,1]))
res.GWAS[,2]<-as.numeric(as.character(res.GWAS[,2]))
res.GWAS[,3]<-as.numeric(as.character(res.GWAS[,3]))
res.GWAS[,4]<-as.numeric(as.character(res.GWAS[,4]))
res.GWAS[,5]<-as.numeric(as.character(res.GWAS[,5]))
res.GWAS[,6]<-as.numeric(as.character(res.GWAS[,6]))
res.GWAS$p.BH<-p.adjust(res.GWAS[,6],"BH")
res.GWAS$p.bonf<-p.adjust(res.GWAS[,6],"bonferroni")
colnames(res.GWAS)[1:6]<-c("CHR", "start","end","nProbes","X2","p.Fisher´s test")
res.GWAS.PFC<-res.GWAS


rm(res.GWAS,ABCA7,ABI3,ADAM10,ADAMTS4,ALPK2,APH1B,APOE,BIN1,CASS4,CD2AP,CD33,CLNK,CNTNAP2,CR1,ECHDC3,EPHA1,HESX1,PICALM,
             TREM2,SLC24A4,SORL1,ZCWPW1,CLU,HS3ST1,MS4A6A)


Gene<-results[which(results$gene=="ADAMTS4"),]
ADAMTS4<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
            "1" ,"NA", Gene$EC.P_Fixed)

Gene<-results[which(results$gene=="CR1"),]
CR1<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
       ((sumlog(Gene$EC.P_Fixed))$df)/2,(sumlog(Gene$EC.P_Fixed))$chisq, 
       (sumlog(Gene$EC.P_Fixed))$p)

Gene<-results[which(results$gene=="BIN1"),]
BIN1<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
        ((sumlog(Gene$EC.P_Fixed))$df)/2,(sumlog(Gene$EC.P_Fixed))$chisq, 
        (sumlog(Gene$EC.P_Fixed))$p)


Gene<-results[which(results$gene=="INPPD5"),]
Gene<-results[which(results$gene=="HESX1"),]
HESX1<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
         ((sumlog(Gene$EC.P_Fixed))$df)/2,(sumlog(Gene$EC.P_Fixed))$chisq, 
         (sumlog(Gene$EC.P_Fixed))$p)



Gene<-results[which(results$gene=="CLNK"),]
CLNK<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
        ((sumlog(Gene$EC.P_Fixed))$df)/2,(sumlog(Gene$EC.P_Fixed))$chisq, 
        (sumlog(Gene$EC.P_Fixed))$p)

Gene<-results[which(results$gene=="HS3ST1"),]

HS3ST1<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
          ((sumlog(Gene$EC.P_Fixed))$df)/2,(sumlog(Gene$EC.P_Fixed))$chisq, 
          (sumlog(Gene$EC.P_Fixed))$p)

Gene<-results[which(results$gene=="HLA-DRB1"),]


Gene<-results[which(results$gene=="CD2AP"),]
CD2AP<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
         ((sumlog(Gene$EC.P_Fixed))$df)/2,(sumlog(Gene$EC.P_Fixed))$chisq, 
         (sumlog(Gene$EC.P_Fixed))$p)

Gene<-results[which(results$gene=="ZCWPW1"),]
ZCWPW1<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
          ((sumlog(Gene$EC.P_Fixed))$df)/2,(sumlog(Gene$EC.P_Fixed))$chisq, 
          (sumlog(Gene$EC.P_Fixed))$p)

Gene<-results[which(results$gene=="EPHA1"),]
EPHA1<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
         ((sumlog(Gene$EC.P_Fixed))$df)/2,(sumlog(Gene$EC.P_Fixed))$chisq, 
         (sumlog(Gene$EC.P_Fixed))$p)

Gene<-results[which(results$gene=="CNTNAP2"),]
CNTNAP2<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
           ((sumlog(Gene$EC.P_Fixed))$df)/2,(sumlog(Gene$EC.P_Fixed))$chisq, 
           (sumlog(Gene$EC.P_Fixed))$p)

Gene<-results[which(results$gene=="CLU"),]
CLU<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
       ((sumlog(Gene$EC.P_Fixed))$df)/2,(sumlog(Gene$EC.P_Fixed))$chisq, 
       (sumlog(Gene$EC.P_Fixed))$p)

Gene<-results[which(results$gene=="ECHDC3"),]
ECHDC3<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
          ((sumlog(Gene$EC.P_Fixed))$df)/2,(sumlog(Gene$EC.P_Fixed))$chisq, 
          (sumlog(Gene$EC.P_Fixed))$p)

Gene<-results[which(results$gene=="MS4A6A"),]
MS4A6A<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
          ((sumlog(Gene$EC.P_Fixed))$df)/2,(sumlog(Gene$EC.P_Fixed))$chisq, 
          (sumlog(Gene$EC.P_Fixed))$p)

Gene<-results[which(results$gene=="PICALM"),]
PICALM<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
          ((sumlog(Gene$EC.P_Fixed))$df)/2,(sumlog(Gene$EC.P_Fixed))$chisq, 
          (sumlog(Gene$EC.P_Fixed))$p)

Gene<-results[which(results$gene=="SORL1"),]
SORL1<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
         ((sumlog(Gene$EC.P_Fixed))$df)/2,(sumlog(Gene$EC.P_Fixed))$chisq, 
         (sumlog(Gene$EC.P_Fixed))$p)

Gene<-results[which(results$gene=="SLC24A4"),]
SLC24A4<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
           ((sumlog(Gene$EC.P_Fixed))$df)/2,(sumlog(Gene$EC.P_Fixed))$chisq, 
           (sumlog(Gene$EC.P_Fixed))$p)

Gene<-results[which(results$gene=="ADAM10"),]
ADAM10<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
          ((sumlog(Gene$EC.P_Fixed))$df)/2,(sumlog(Gene$EC.P_Fixed))$chisq, 
          (sumlog(Gene$EC.P_Fixed))$p)

Gene<-results[which(results$gene=="APH1B"),]
APH1B<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
         ((sumlog(Gene$EC.P_Fixed))$df)/2,(sumlog(Gene$EC.P_Fixed))$chisq, 
         (sumlog(Gene$EC.P_Fixed))$p)

Gene<-results[which(results$gene=="KAT8"),]


Gene<-results[which(results$gene=="ECHDC3"),]
ECHDC3<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
          ((sumlog(Gene$EC.P_Fixed))$df)/2,(sumlog(Gene$EC.P_Fixed))$chisq, 
          (sumlog(Gene$EC.P_Fixed))$p)

Gene<-results[which(results$gene=="ABI3"),]
ABI3<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
        ((sumlog(Gene$EC.P_Fixed))$df)/2,(sumlog(Gene$EC.P_Fixed))$chisq, 
        (sumlog(Gene$EC.P_Fixed))$p)

Gene<-results[which(results$gene=="SUZ12P1"),]

Gene<-results[which(results$gene=="ALPK2"),]
ALPK2<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
         ((sumlog(Gene$EC.P_Fixed))$df)/2,(sumlog(Gene$EC.P_Fixed))$chisq, 
         (sumlog(Gene$EC.P_Fixed))$p)

Gene<-results[which(results$gene=="ABCA7"),]
ABCA7<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
         ((sumlog(Gene$EC.P_Fixed))$df)/2,(sumlog(Gene$EC.P_Fixed))$chisq, 
         (sumlog(Gene$EC.P_Fixed))$p)

Gene<-results[which(results$gene=="APOE"),]
APOE<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
        ((sumlog(Gene$EC.P_Fixed))$df)/2,(sumlog(Gene$EC.P_Fixed))$chisq, 
        (sumlog(Gene$EC.P_Fixed))$p)

Gene<-results[which(results$gene=="AC074212.3"),]

Gene<-results[which(results$gene=="CD33"),]
CD33<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
        ((sumlog(Gene$EC.P_Fixed))$df)/2,(sumlog(Gene$EC.P_Fixed))$chisq, 
        (sumlog(Gene$EC.P_Fixed))$p)
Gene<-results[which(results$gene=="CASS4"),]
CASS4<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
         ((sumlog(Gene$EC.P_Fixed))$df)/2,(sumlog(Gene$EC.P_Fixed))$chisq, 
         (sumlog(Gene$EC.P_Fixed))$p)

Gene<-results[which(results$gene=="TREM2"),]
TREM2<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
         ((sumlog(Gene$EC.P_Fixed))$df)/2,(sumlog(Gene$EC.P_Fixed))$chisq, 
         (sumlog(Gene$EC.P_Fixed))$p)


res.GWAS<-rbind(ABCA7,ABI3,ADAM10,ADAMTS4,ALPK2,APH1B,APOE,BIN1,CASS4,CD2AP,CD33,CLNK,CNTNAP2,CR1,ECHDC3,EPHA1,HESX1,PICALM,
                TREM2,SLC24A4,SORL1,ZCWPW1,CLU,HS3ST1,MS4A6A)

res.GWAS<-data.frame(res.GWAS)
str(res.GWAS)

res.GWAS[,1]<-NULL
res.GWAS[,1]<-as.numeric(as.character(res.GWAS[,1]))
res.GWAS[,2]<-as.numeric(as.character(res.GWAS[,2]))
res.GWAS[,3]<-as.numeric(as.character(res.GWAS[,3]))
res.GWAS[,4]<-as.numeric(as.character(res.GWAS[,4]))
res.GWAS[,5]<-as.numeric(as.character(res.GWAS[,5]))
res.GWAS[,6]<-as.numeric(as.character(res.GWAS[,6]))
res.GWAS$p.BH<-p.adjust(res.GWAS[,6],"BH")
res.GWAS$p.bonf<-p.adjust(res.GWAS[,6],"bonferroni")
colnames(res.GWAS)[1:6]<-c("CHR", "start","end","nProbes","X2","p.Fisher´s test")
res.GWAS.EC<-res.GWAS

rm(res.GWAS, ABCA7,ABI3,ADAM10,ADAMTS4,ALPK2,APH1B,APOE,BIN1,CASS4,CD2AP,CD33,CLNK,CNTNAP2,CR1,ECHDC3,EPHA1,HESX1,PICALM,
   TREM2,SLC24A4,SORL1,ZCWPW1,CLU,HS3ST1,MS4A6A)




Gene<-results[which(results$gene=="ADAMTS4"),]
ADAMTS4<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
           "1" ,"NA", Gene$TG.P_Fixed)

Gene<-results[which(results$gene=="CR1"),]
CR1<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
       ((sumlog(Gene$TG.P_Fixed))$df)/2,(sumlog(Gene$TG.P_Fixed))$chisq, 
       (sumlog(Gene$TG.P_Fixed))$p)

Gene<-results[which(results$gene=="BIN1"),]
BIN1<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
        ((sumlog(Gene$TG.P_Fixed))$df)/2,(sumlog(Gene$TG.P_Fixed))$chisq, 
        (sumlog(Gene$TG.P_Fixed))$p)


Gene<-results[which(results$gene=="INPPD5"),]
Gene<-results[which(results$gene=="HESX1"),]
HESX1<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
         ((sumlog(Gene$TG.P_Fixed))$df)/2,(sumlog(Gene$TG.P_Fixed))$chisq, 
         (sumlog(Gene$TG.P_Fixed))$p)



Gene<-results[which(results$gene=="CLNK"),]
CLNK<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
        ((sumlog(Gene$TG.P_Fixed))$df)/2,(sumlog(Gene$TG.P_Fixed))$chisq, 
        (sumlog(Gene$TG.P_Fixed))$p)

Gene<-results[which(results$gene=="HS3ST1"),]

HS3ST1<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
          ((sumlog(Gene$TG.P_Fixed))$df)/2,(sumlog(Gene$TG.P_Fixed))$chisq, 
          (sumlog(Gene$TG.P_Fixed))$p)

Gene<-results[which(results$gene=="HLA-DRB1"),]


Gene<-results[which(results$gene=="CD2AP"),]
CD2AP<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
         ((sumlog(Gene$TG.P_Fixed))$df)/2,(sumlog(Gene$TG.P_Fixed))$chisq, 
         (sumlog(Gene$TG.P_Fixed))$p)

Gene<-results[which(results$gene=="ZCWPW1"),]
ZCWPW1<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
          ((sumlog(Gene$TG.P_Fixed))$df)/2,(sumlog(Gene$TG.P_Fixed))$chisq, 
          (sumlog(Gene$TG.P_Fixed))$p)

Gene<-results[which(results$gene=="EPHA1"),]
EPHA1<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
         ((sumlog(Gene$TG.P_Fixed))$df)/2,(sumlog(Gene$TG.P_Fixed))$chisq, 
         (sumlog(Gene$TG.P_Fixed))$p)

Gene<-results[which(results$gene=="CNTNAP2"),]
CNTNAP2<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
           ((sumlog(Gene$TG.P_Fixed))$df)/2,(sumlog(Gene$TG.P_Fixed))$chisq, 
           (sumlog(Gene$TG.P_Fixed))$p)

Gene<-results[which(results$gene=="CLU"),]
CLU<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
       ((sumlog(Gene$TG.P_Fixed))$df)/2,(sumlog(Gene$TG.P_Fixed))$chisq, 
       (sumlog(Gene$TG.P_Fixed))$p)

Gene<-results[which(results$gene=="ECHDC3"),]
ECHDC3<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
          ((sumlog(Gene$TG.P_Fixed))$df)/2,(sumlog(Gene$TG.P_Fixed))$chisq, 
          (sumlog(Gene$TG.P_Fixed))$p)

Gene<-results[which(results$gene=="MS4A6A"),]
MS4A6A<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
          ((sumlog(Gene$TG.P_Fixed))$df)/2,(sumlog(Gene$TG.P_Fixed))$chisq, 
          (sumlog(Gene$TG.P_Fixed))$p)

Gene<-results[which(results$gene=="PICALM"),]
PICALM<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
          ((sumlog(Gene$TG.P_Fixed))$df)/2,(sumlog(Gene$TG.P_Fixed))$chisq, 
          (sumlog(Gene$TG.P_Fixed))$p)

Gene<-results[which(results$gene=="SORL1"),]
SORL1<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
         ((sumlog(Gene$TG.P_Fixed))$df)/2,(sumlog(Gene$TG.P_Fixed))$chisq, 
         (sumlog(Gene$TG.P_Fixed))$p)

Gene<-results[which(results$gene=="SLC24A4"),]
SLC24A4<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
           ((sumlog(Gene$TG.P_Fixed))$df)/2,(sumlog(Gene$TG.P_Fixed))$chisq, 
           (sumlog(Gene$TG.P_Fixed))$p)

Gene<-results[which(results$gene=="ADAM10"),]
ADAM10<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
          ((sumlog(Gene$TG.P_Fixed))$df)/2,(sumlog(Gene$TG.P_Fixed))$chisq, 
          (sumlog(Gene$TG.P_Fixed))$p)

Gene<-results[which(results$gene=="APH1B"),]
APH1B<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
         ((sumlog(Gene$TG.P_Fixed))$df)/2,(sumlog(Gene$TG.P_Fixed))$chisq, 
         (sumlog(Gene$TG.P_Fixed))$p)

Gene<-results[which(results$gene=="KAT8"),]


Gene<-results[which(results$gene=="ECHDC3"),]
ECHDC3<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
          ((sumlog(Gene$TG.P_Fixed))$df)/2,(sumlog(Gene$TG.P_Fixed))$chisq, 
          (sumlog(Gene$TG.P_Fixed))$p)

Gene<-results[which(results$gene=="ABI3"),]
ABI3<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
        ((sumlog(Gene$TG.P_Fixed))$df)/2,(sumlog(Gene$TG.P_Fixed))$chisq, 
        (sumlog(Gene$TG.P_Fixed))$p)

Gene<-results[which(results$gene=="SUZ12P1"),]

Gene<-results[which(results$gene=="ALPK2"),]
ALPK2<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
         ((sumlog(Gene$TG.P_Fixed))$df)/2,(sumlog(Gene$TG.P_Fixed))$chisq, 
         (sumlog(Gene$TG.P_Fixed))$p)

Gene<-results[which(results$gene=="ABCA7"),]
ABCA7<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
         ((sumlog(Gene$TG.P_Fixed))$df)/2,(sumlog(Gene$TG.P_Fixed))$chisq, 
         (sumlog(Gene$TG.P_Fixed))$p)

Gene<-results[which(results$gene=="APOE"),]
APOE<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
        ((sumlog(Gene$TG.P_Fixed))$df)/2,(sumlog(Gene$TG.P_Fixed))$chisq, 
        (sumlog(Gene$TG.P_Fixed))$p)

Gene<-results[which(results$gene=="AC074212.3"),]

Gene<-results[which(results$gene=="CD33"),]
CD33<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
        ((sumlog(Gene$TG.P_Fixed))$df)/2,(sumlog(Gene$TG.P_Fixed))$chisq, 
        (sumlog(Gene$TG.P_Fixed))$p)
Gene<-results[which(results$gene=="CASS4"),]
CASS4<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
         ((sumlog(Gene$TG.P_Fixed))$df)/2,(sumlog(Gene$TG.P_Fixed))$chisq, 
         (sumlog(Gene$TG.P_Fixed))$p)

Gene<-results[which(results$gene=="TREM2"),]
TREM2<-c(unique(as.character(Gene$gene)), unique(as.character(Gene$CHR)), min(Gene$MAPINFO),max(Gene$MAPINFO),
         ((sumlog(Gene$TG.P_Fixed))$df)/2,(sumlog(Gene$TG.P_Fixed))$chisq, 
         (sumlog(Gene$TG.P_Fixed))$p)


res.GWAS<-rbind(ABCA7,ABI3,ADAM10,ADAMTS4,ALPK2,APH1B,APOE,BIN1,CASS4,CD2AP,CD33,CLNK,CNTNAP2,CR1,ECHDC3,EPHA1,HESX1,PICALM,
                TREM2,SLC24A4,SORL1,ZCWPW1,CLU,HS3ST1,MS4A6A)

res.GWAS<-data.frame(res.GWAS)
str(res.GWAS)

res.GWAS[,1]<-NULL
res.GWAS[,1]<-as.numeric(as.character(res.GWAS[,1]))
res.GWAS[,2]<-as.numeric(as.character(res.GWAS[,2]))
res.GWAS[,3]<-as.numeric(as.character(res.GWAS[,3]))
res.GWAS[,4]<-as.numeric(as.character(res.GWAS[,4]))
res.GWAS[,5]<-as.numeric(as.character(res.GWAS[,5]))
res.GWAS[,6]<-as.numeric(as.character(res.GWAS[,6]))
res.GWAS$p.BH<-p.adjust(res.GWAS[,6],"BH")
res.GWAS$p.bonf<-p.adjust(res.GWAS[,6],"bonferroni")
colnames(res.GWAS)[1:6]<-c("CHR", "start","end","nProbes","X2","p.Fisher´s test")
res.GWAS.TG<-res.GWAS

res.GWAS.PFC[order(rownames(res.GWAS.PFC)),]
res.GWAS.TG[order(rownames(res.GWAS.TG)),]
res.GWAS.EC[order(rownames(res.GWAS.EC)),]


res.gwas.enrichment<-cbind(res.GWAS.CC,res.GWAS.PFC,res.GWAS.EC,res.GWAS.TG)
write.csv(res.gwas.enrichment,file="GWAS_Janssen_Coverage.csv")


data.qq<- data.frame(rownames(res.all.f),res.all.f$CB.P_Fixed, res.all.f$TG.P_Fixed, 
                    res.all.f$PFC.P_Fixed,res.all.f$EC.P_Fixed, res.all.f$pval.cc, res.all.f$gene) 
head(data.qq)
rownames(data.qq)<-data.qq[,1]
data.qq[,1]<-NULL
head(data.qq)

str(data.qq)
dat<-na.omit(data.qq)

library(Haplin)

chisq <- qchisq(1-(as.numeric(dat[,1])),1)
lambda = median(chisq)/qchisq(0.5,1)
p<-dat[,1]
names(p)<-paste0(rownames(dat),"_",dat$res.all.f.gene)

bmp("CB.bmp" ,width=150, height=150, units="mm", res=300)
pQQ(p, nlabs = 3, conf = 0.95,  mark = 0.05,lim=c(0,10), sub= paste0("Lambda=",round(lambda , digits=3)), main="CB")
dev.off()


chisq <- qchisq(1-(as.numeric(dat[,2])),1)
lambda = median(chisq)/qchisq(0.5,1)
p<-dat[,2]
names(p)<-paste0(rownames(dat),"_",dat$res.all.f.gene)

bmp("TG.bmp" ,width=150, height=150, units="mm", res=300)
pQQ(p, nlabs = 2, conf = 0.95,  mark = 0.05,lim=c(0,17), sub= paste0("Lambda=",round(lambda , digits=3)), main="TG")
dev.off()



chisq <- qchisq(1-(as.numeric(dat[,3])),1)
lambda = median(chisq)/qchisq(0.5,1)
p<-dat[,3]
names(p)<-paste0(rownames(dat),"_",dat$res.all.f.gene)

bmp("PFC.bmp" ,width=150, height=150, units="mm", res=300)
pQQ(p, nlabs = 2, conf = 0.95,  mark = 0.05,lim=c(0,17), sub= paste0("Lambda=",round(lambda , digits=3)), main="PFC")
dev.off()


chisq <- qchisq(1-(as.numeric(dat[,4])),1)
lambda = median(chisq)/qchisq(0.5,1)
p<-dat[,4]
names(p)<-paste0(rownames(dat),"_",dat$res.all.f.gene)

bmp("EC.bmp" ,width=150, height=150, units="mm", res=300)
pQQ(p, nlabs =2, conf = 0.95,  mark = 0.05,lim=c(0,17), sub= paste0("Lambda=",round(lambda , digits=3)), main="EC")
dev.off()

chisq <- qchisq(1-(as.numeric(dat[,5])),1)
lambda = median(chisq)/qchisq(0.5,1)
p<-dat[,5]
names(p)<-paste0(rownames(dat),"_",dat$res.all.f.gene)

bmp("CC.bmp" ,width=150, height=150, units="mm", res=300)
pQQ(p, nlabs = 3, conf = 0.95,  mark = 0.05,lim=c(0,20), sub= paste0("Lambda=",round(lambda , digits=3)), main="Cross cortical")
dev.off()






setwd("/mnt/data1/Ehsan/meta_CrossCortical/")
load("RESULTS_CROSS.CORTICAL.META.Rda")
res.meta.bonf<-results[which(results$fdr.FE<0.05),]
res.meta.bacon<-results
CHR.enrich<-data.frame(table(res.meta.bonf$CHR),table(res.meta.bacon$CHR))
CHR.enrich$Var1.1<-NULL
CHR.enrich<-CHR.enrich[-c(1,24,25),]
CHR.enrich
CHR.enrich$Var1<-as.numeric(as.character(CHR.enrich$Var1))
CHR.enrich<-CHR.enrich[order(CHR.enrich$Var1),]
CHR.enrich
colnames(CHR.enrich)[1]<-"CHR"

sum<-t(CHR.enrich)
colnames(sum)<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12",
                 "chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22")
sum<-sum[-1,]
rownames(sum)<-c("GS_Probes","All_Probes")
sum<-data.frame(sum)
sum$TOTAL<-data.frame(c(length(res.meta.bonf$CHR),length(res.meta.bacon$CHR)))
colnames(sum)[23]<-c("TOTAL")


out<-matrix(data=NA, nrow=22, ncol=2)
rownames(out)<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12",
                 "chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22")
colnames(out)<-c("OR","p")
for(i in 1:22){
  out[i,]<-c((fisher.test(sum[c(1,2),c(i,23)], alternative = "greater"))$estimate, (fisher.test(sum[c(1,2),c(i,23)], alternative = "greater"))$p.value)
  }
write.csv(out, file="ChromosomeenrichmentanalysesFDR.csv")


GS_Probes <-((CHR.enrich[,2]/length(res.meta.bonf$CHR))*100) 
All_Probes <-((CHR.enrich[,3]/length(res.meta.bacon$CHR))*100) 


totalper<-rbind(All_Probes,GS_Probes)
totalper
colnames(totalper)<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12",
                      "chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22")



LOC1<-results[(which(results$CHR ==1 & results$MAPINFO > 207653395-1 & results$MAPINFO < 207850539+1)),]
LOC2<-results[(which(results$CHR ==2 & results$MAPINFO > 127839781-1 & results$MAPINFO < 127894615+1)),]
LOC3<-results[(which(results$CHR ==4 & results$MAPINFO > 11700350-1 & results$MAPINFO < 11723235+1)),]
LOC4<-results[(which(results$CHR ==6 & results$MAPINFO > 47432637-1 & results$MAPINFO < 47432637+1)),]
LOC5<-results[(which(results$CHR ==7 & results$MAPINFO > 99971313-1 & results$MAPINFO < 99971834+1)),]
LOC6<-results[(which(results$CHR ==7 & results$MAPINFO > 143099107-1 & results$MAPINFO < 143129306+1)),]
LOC7<-results[(which(results$CHR ==8 & results$MAPINFO > 27456253-1 & results$MAPINFO < 27467821+1)),]
LOC8<-results[(which(results$CHR ==10 & results$MAPINFO > 11721102-1 & results$MAPINFO < 11721477+1)),]
LOC9<-results[(which(results$CHR ==11 & results$MAPINFO > 59857581-1 & results$MAPINFO < 60103385+1)),]
LOC10<-results[(which(results$CHR ==11 & results$MAPINFO > 85655105-1 & results$MAPINFO < 85869737+1)),]
LOC11<-results[(which(results$CHR ==17 & results$MAPINFO > 56404349-1 & results$MAPINFO < 56410041+1)),]
LOC12<-results[(which(results$CHR ==19 & results$MAPINFO > 1056492-1 & results$MAPINFO < 1063443+1)),]
LOC13<-results[(which(results$CHR ==19 & results$MAPINFO > 45065601-1 & results$MAPINFO < 45733897+1)),]


Gene<-LOC1
LOC1<-c(unique(as.character(Gene$CHR)), min(Gene$MAPINFO), max(Gene$MAPINFO),
         ((sumlog(Gene$pval.fixed))$df)/2,(sumlog(Gene$pval.fixed))$chisq, 
         (sumlog(Gene$pval.fixed))$p)
Gene<-LOC2
LOC2<-c(unique(as.character(Gene$CHR)), min(Gene$MAPINFO), max(Gene$MAPINFO),
        ((sumlog(Gene$pval.fixed))$df)/2,(sumlog(Gene$pval.fixed))$chisq, 
        (sumlog(Gene$pval.fixed))$p)
Gene<-LOC6
LOC6<-c(unique(as.character(Gene$CHR)), min(Gene$MAPINFO), max(Gene$MAPINFO),
        ((sumlog(Gene$pval.fixed))$df)/2,(sumlog(Gene$pval.fixed))$chisq, 
        (sumlog(Gene$pval.fixed))$p)
Gene<-LOC7
LOC7<-c(unique(as.character(Gene$CHR)), min(Gene$MAPINFO), max(Gene$MAPINFO),
        ((sumlog(Gene$pval.fixed))$df)/2,(sumlog(Gene$pval.fixed))$chisq, 
        (sumlog(Gene$pval.fixed))$p)
Gene<-LOC9
LOC9<-c(unique(as.character(Gene$CHR)), min(Gene$MAPINFO), max(Gene$MAPINFO),
        ((sumlog(Gene$pval.fixed))$df)/2,(sumlog(Gene$pval.fixed))$chisq, 
        (sumlog(Gene$pval.fixed))$p)
Gene<-LOC10
LOC10<-c(unique(as.character(Gene$CHR)), min(Gene$MAPINFO), max(Gene$MAPINFO),
        ((sumlog(Gene$pval.fixed))$df)/2,(sumlog(Gene$pval.fixed))$chisq, 
        (sumlog(Gene$pval.fixed))$p)
Gene<-LOC11
LOC11<-c(unique(as.character(Gene$CHR)), min(Gene$MAPINFO), max(Gene$MAPINFO),
        ((sumlog(Gene$pval.fixed))$df)/2,(sumlog(Gene$pval.fixed))$chisq, 
        (sumlog(Gene$pval.fixed))$p)
Gene<-LOC13
LOC13<-c(unique(as.character(Gene$CHR)), min(Gene$MAPINFO), max(Gene$MAPINFO),
         ((sumlog(Gene$pval.fixed))$df)/2,(sumlog(Gene$pval.fixed))$chisq, 
         (sumlog(Gene$pval.fixed))$p)

CC.LOC<-data.frame(rbind(LOC1, LOC2, LOC6, LOC7, LOC9, LOC10,LOC11,LOC13))
CC.LOC$X6<-as.numeric(as.character(CC.LOC$X6))
CC.LOC$paadj<-p.adjust(as.numeric(CC.LOC[,6]),"fdr")

LOC1<-results[(which(results$CHR ==1 & results$MAPINFO > 207653395-1 & results$MAPINFO < 207850539+1)),]
LOC2<-results[(which(results$CHR ==2 & results$MAPINFO > 127839781-1 & results$MAPINFO < 127894615+1)),]
LOC3<-results[(which(results$CHR ==4 & results$MAPINFO > 11700350-1 & results$MAPINFO < 11723235+1)),]
LOC4<-results[(which(results$CHR ==6 & results$MAPINFO > 47432637-1 & results$MAPINFO < 47432637+1)),]
LOC5<-results[(which(results$CHR ==7 & results$MAPINFO > 99971313-1 & results$MAPINFO < 99971834+1)),]
LOC6<-results[(which(results$CHR ==7 & results$MAPINFO > 143099107-1 & results$MAPINFO < 143129306+1)),]
LOC7<-results[(which(results$CHR ==8 & results$MAPINFO > 27456253-1 & results$MAPINFO < 27467821+1)),]
LOC8<-results[(which(results$CHR ==10 & results$MAPINFO > 11721102-1 & results$MAPINFO < 11721477+1)),]
LOC9<-results[(which(results$CHR ==11 & results$MAPINFO > 59857581-1 & results$MAPINFO < 60103385+1)),]
LOC10<-results[(which(results$CHR ==11 & results$MAPINFO > 85655105-1 & results$MAPINFO < 85869737+1)),]
LOC11<-results[(which(results$CHR ==17 & results$MAPINFO > 56404349-1 & results$MAPINFO < 56410041+1)),]
LOC12<-results[(which(results$CHR ==19 & results$MAPINFO > 1056492-1 & results$MAPINFO < 1063443+1)),]
LOC13<-results[(which(results$CHR ==19 & results$MAPINFO > 45065601-1 & results$MAPINFO < 45733897+1)),]


Gene<-LOC1
LOC1<-c(unique(as.character(Gene$CHR)), min(Gene$MAPINFO), max(Gene$MAPINFO),
        ((sumlog(Gene$PFC.P_Fixed))$df)/2,(sumlog(Gene$PFC.P_Fixed))$chisq, 
        (sumlog(Gene$PFC.P_Fixed))$p)
Gene<-LOC2
LOC2<-c(unique(as.character(Gene$CHR)), min(Gene$MAPINFO), max(Gene$MAPINFO),
        ((sumlog(Gene$PFC.P_Fixed))$df)/2,(sumlog(Gene$PFC.P_Fixed))$chisq, 
        (sumlog(Gene$PFC.P_Fixed))$p)
Gene<-LOC6
LOC6<-c(unique(as.character(Gene$CHR)), min(Gene$MAPINFO), max(Gene$MAPINFO),
        ((sumlog(Gene$PFC.P_Fixed))$df)/2,(sumlog(Gene$PFC.P_Fixed))$chisq, 
        (sumlog(Gene$PFC.P_Fixed))$p)
Gene<-LOC7
LOC7<-c(unique(as.character(Gene$CHR)), min(Gene$MAPINFO), max(Gene$MAPINFO),
        ((sumlog(Gene$PFC.P_Fixed))$df)/2,(sumlog(Gene$PFC.P_Fixed))$chisq, 
        (sumlog(Gene$PFC.P_Fixed))$p)
Gene<-LOC9
LOC9<-c(unique(as.character(Gene$CHR)), min(Gene$MAPINFO), max(Gene$MAPINFO),
        ((sumlog(Gene$PFC.P_Fixed))$df)/2,(sumlog(Gene$PFC.P_Fixed))$chisq, 
        (sumlog(Gene$PFC.P_Fixed))$p)
Gene<-LOC10
LOC10<-c(unique(as.character(Gene$CHR)), min(Gene$MAPINFO), max(Gene$MAPINFO),
         ((sumlog(Gene$PFC.P_Fixed))$df)/2,(sumlog(Gene$PFC.P_Fixed))$chisq, 
         (sumlog(Gene$PFC.P_Fixed))$p)
Gene<-LOC11
LOC11<-c(unique(as.character(Gene$CHR)), min(Gene$MAPINFO), max(Gene$MAPINFO),
         ((sumlog(Gene$PFC.P_Fixed))$df)/2,(sumlog(Gene$PFC.P_Fixed))$chisq, 
         (sumlog(Gene$PFC.P_Fixed))$p)
Gene<-LOC13
LOC13<-c(unique(as.character(Gene$CHR)), min(Gene$MAPINFO), max(Gene$MAPINFO),
         ((sumlog(Gene$PFC.P_Fixed))$df)/2,(sumlog(Gene$PFC.P_Fixed))$chisq, 
         (sumlog(Gene$PFC.P_Fixed))$p)

PFC.LOC<-data.frame(rbind(LOC1, LOC2, LOC6, LOC7, LOC9, LOC10,LOC11,LOC13))
PFC.LOC$X6<-as.numeric(as.character(PFC.LOC$X6))
PFC.LOC$paadj<-p.adjust(as.numeric(PFC.LOC[,6]),"fdr")

PFC.LOC



LOC1<-results[(which(results$CHR ==1 & results$MAPINFO > 207653395-1 & results$MAPINFO < 207850539+1)),]
LOC2<-results[(which(results$CHR ==2 & results$MAPINFO > 127839781-1 & results$MAPINFO < 127894615+1)),]
LOC3<-results[(which(results$CHR ==4 & results$MAPINFO > 11700350-1 & results$MAPINFO < 11723235+1)),]
LOC4<-results[(which(results$CHR ==6 & results$MAPINFO > 47432637-1 & results$MAPINFO < 47432637+1)),]
LOC5<-results[(which(results$CHR ==7 & results$MAPINFO > 99971313-1 & results$MAPINFO < 99971834+1)),]
LOC6<-results[(which(results$CHR ==7 & results$MAPINFO > 143099107-1 & results$MAPINFO < 143129306+1)),]
LOC7<-results[(which(results$CHR ==8 & results$MAPINFO > 27456253-1 & results$MAPINFO < 27467821+1)),]
LOC8<-results[(which(results$CHR ==10 & results$MAPINFO > 11721102-1 & results$MAPINFO < 11721477+1)),]
LOC9<-results[(which(results$CHR ==11 & results$MAPINFO > 59857581-1 & results$MAPINFO < 60103385+1)),]
LOC10<-results[(which(results$CHR ==11 & results$MAPINFO > 85655105-1 & results$MAPINFO < 85869737+1)),]
LOC11<-results[(which(results$CHR ==17 & results$MAPINFO > 56404349-1 & results$MAPINFO < 56410041+1)),]
LOC12<-results[(which(results$CHR ==19 & results$MAPINFO > 1056492-1 & results$MAPINFO < 1063443+1)),]
LOC13<-results[(which(results$CHR ==19 & results$MAPINFO > 45065601-1 & results$MAPINFO < 45733897+1)),]


Gene<-LOC1
LOC1<-c(unique(as.character(Gene$CHR)), min(Gene$MAPINFO), max(Gene$MAPINFO),
        ((sumlog(Gene$TG.P_Fixed))$df)/2,(sumlog(Gene$TG.P_Fixed))$chisq, 
        (sumlog(Gene$TG.P_Fixed))$p)
Gene<-LOC2
LOC2<-c(unique(as.character(Gene$CHR)), min(Gene$MAPINFO), max(Gene$MAPINFO),
        ((sumlog(Gene$TG.P_Fixed))$df)/2,(sumlog(Gene$TG.P_Fixed))$chisq, 
        (sumlog(Gene$TG.P_Fixed))$p)
Gene<-LOC6
LOC6<-c(unique(as.character(Gene$CHR)), min(Gene$MAPINFO), max(Gene$MAPINFO),
        ((sumlog(Gene$TG.P_Fixed))$df)/2,(sumlog(Gene$TG.P_Fixed))$chisq, 
        (sumlog(Gene$TG.P_Fixed))$p)
Gene<-LOC7
LOC7<-c(unique(as.character(Gene$CHR)), min(Gene$MAPINFO), max(Gene$MAPINFO),
        ((sumlog(Gene$TG.P_Fixed))$df)/2,(sumlog(Gene$TG.P_Fixed))$chisq, 
        (sumlog(Gene$TG.P_Fixed))$p)
Gene<-LOC9
LOC9<-c(unique(as.character(Gene$CHR)), min(Gene$MAPINFO), max(Gene$MAPINFO),
        ((sumlog(Gene$TG.P_Fixed))$df)/2,(sumlog(Gene$TG.P_Fixed))$chisq, 
        (sumlog(Gene$TG.P_Fixed))$p)
Gene<-LOC10
LOC10<-c(unique(as.character(Gene$CHR)), min(Gene$MAPINFO), max(Gene$MAPINFO),
         ((sumlog(Gene$TG.P_Fixed))$df)/2,(sumlog(Gene$TG.P_Fixed))$chisq, 
         (sumlog(Gene$TG.P_Fixed))$p)
Gene<-LOC11
LOC11<-c(unique(as.character(Gene$CHR)), min(Gene$MAPINFO), max(Gene$MAPINFO),
         ((sumlog(Gene$TG.P_Fixed))$df)/2,(sumlog(Gene$TG.P_Fixed))$chisq, 
         (sumlog(Gene$TG.P_Fixed))$p)
Gene<-LOC13
LOC13<-c(unique(as.character(Gene$CHR)), min(Gene$MAPINFO), max(Gene$MAPINFO),
         ((sumlog(Gene$TG.P_Fixed))$df)/2,(sumlog(Gene$TG.P_Fixed))$chisq, 
         (sumlog(Gene$TG.P_Fixed))$p)

TG.LOC<-data.frame(rbind(LOC1, LOC2, LOC6, LOC7, LOC9, LOC10,LOC11,LOC13))
TG.LOC$X6<-as.numeric(as.character(TG.LOC$X6))
TG.LOC$paadj<-p.adjust(as.numeric(TG.LOC[,6]),"fdr")

TG.LOC




LOC1<-results[(which(results$CHR ==1 & results$MAPINFO > 207653395-1 & results$MAPINFO < 207850539+1)),]
LOC2<-results[(which(results$CHR ==2 & results$MAPINFO > 127839781-1 & results$MAPINFO < 127894615+1)),]
LOC3<-results[(which(results$CHR ==4 & results$MAPINFO > 11700350-1 & results$MAPINFO < 11723235+1)),]
LOC4<-results[(which(results$CHR ==6 & results$MAPINFO > 47432637-1 & results$MAPINFO < 47432637+1)),]
LOC5<-results[(which(results$CHR ==7 & results$MAPINFO > 99971313-1 & results$MAPINFO < 99971834+1)),]
LOC6<-results[(which(results$CHR ==7 & results$MAPINFO > 143099107-1 & results$MAPINFO < 143129306+1)),]
LOC7<-results[(which(results$CHR ==8 & results$MAPINFO > 27456253-1 & results$MAPINFO < 27467821+1)),]
LOC8<-results[(which(results$CHR ==10 & results$MAPINFO > 11721102-1 & results$MAPINFO < 11721477+1)),]
LOC9<-results[(which(results$CHR ==11 & results$MAPINFO > 59857581-1 & results$MAPINFO < 60103385+1)),]
LOC10<-results[(which(results$CHR ==11 & results$MAPINFO > 85655105-1 & results$MAPINFO < 85869737+1)),]
LOC11<-results[(which(results$CHR ==17 & results$MAPINFO > 56404349-1 & results$MAPINFO < 56410041+1)),]
LOC12<-results[(which(results$CHR ==19 & results$MAPINFO > 1056492-1 & results$MAPINFO < 1063443+1)),]
LOC13<-results[(which(results$CHR ==19 & results$MAPINFO > 45065601-1 & results$MAPINFO < 45733897+1)),]
LOC<-data.frame(rbind(LOC1, LOC2, LOC6, LOC7, LOC9, LOC10,LOC11,LOC13))


Gene<-LOC1
LOC1<-c(unique(as.character(Gene$CHR)), min(Gene$MAPINFO), max(Gene$MAPINFO),
        ((sumlog(Gene$EC.P_Fixed))$df)/2,(sumlog(Gene$EC.P_Fixed))$chisq, 
        (sumlog(Gene$EC.P_Fixed))$p)
Gene<-LOC2
LOC2<-c(unique(as.character(Gene$CHR)), min(Gene$MAPINFO), max(Gene$MAPINFO),
        ((sumlog(Gene$EC.P_Fixed))$df)/2,(sumlog(Gene$EC.P_Fixed))$chisq, 
        (sumlog(Gene$EC.P_Fixed))$p)
Gene<-LOC6
LOC6<-c(unique(as.character(Gene$CHR)), min(Gene$MAPINFO), max(Gene$MAPINFO),
        ((sumlog(Gene$EC.P_Fixed))$df)/2,(sumlog(Gene$EC.P_Fixed))$chisq, 
        (sumlog(Gene$EC.P_Fixed))$p)
Gene<-LOC7
LOC7<-c(unique(as.character(Gene$CHR)), min(Gene$MAPINFO), max(Gene$MAPINFO),
        ((sumlog(Gene$EC.P_Fixed))$df)/2,(sumlog(Gene$EC.P_Fixed))$chisq, 
        (sumlog(Gene$EC.P_Fixed))$p)
Gene<-LOC9
LOC9<-c(unique(as.character(Gene$CHR)), min(Gene$MAPINFO), max(Gene$MAPINFO),
        ((sumlog(Gene$EC.P_Fixed))$df)/2,(sumlog(Gene$EC.P_Fixed))$chisq, 
        (sumlog(Gene$EC.P_Fixed))$p)
Gene<-LOC10
LOC10<-c(unique(as.character(Gene$CHR)), min(Gene$MAPINFO), max(Gene$MAPINFO),
         ((sumlog(Gene$EC.P_Fixed))$df)/2,(sumlog(Gene$EC.P_Fixed))$chisq, 
         (sumlog(Gene$EC.P_Fixed))$p)
Gene<-LOC11
LOC11<-c(unique(as.character(Gene$CHR)), min(Gene$MAPINFO), max(Gene$MAPINFO),
         ((sumlog(Gene$EC.P_Fixed))$df)/2,(sumlog(Gene$EC.P_Fixed))$chisq, 
         (sumlog(Gene$EC.P_Fixed))$p)
Gene<-LOC13
LOC13<-c(unique(as.character(Gene$CHR)), min(Gene$MAPINFO), max(Gene$MAPINFO),
         ((sumlog(Gene$EC.P_Fixed))$df)/2,(sumlog(Gene$EC.P_Fixed))$chisq, 
         (sumlog(Gene$EC.P_Fixed))$p)

EC.LOC<-data.frame(rbind(LOC1, LOC2, LOC6, LOC7, LOC9, LOC10,LOC11,LOC13))
EC.LOC$X6<-as.numeric(as.character(EC.LOC$X6))
EC.LOC$paadj<-p.adjust(as.numeric(EC.LOC[,6]),"fdr")

EC.LOC


as.character(unique(LOC1$UCSC_RefGene_Name))
as.character(unique(LOC2$UCSC_RefGene_Name))
as.character(unique(LOC6$UCSC_RefGene_Name))
as.character(unique(LOC7$UCSC_RefGene_Name))
as.character(unique(LOC9$UCSC_RefGene_Name))
as.character(unique(LOC10$UCSC_RefGene_Name))
as.character(unique(LOC11$UCSC_RefGene_Name))
as.character(unique(LOC13$UCSC_RefGene_Name))




setwd("/mnt/data1/Ehsan/meta_CrossCortical/")
load("RESULTS_CROSS.CORTICAL.META.Rda")

res.meta.bonf<-results[which(results$fdr.FE<0.05),]
LOC1<-results[(which(results$CHR ==1 & results$MAPINFO > 207653395-1 & results$MAPINFO < 207850539+1)),]


loc1.enrich<-data.frame(table(res.meta.bonf$CHR),table(LOC1$CHR))
CHR.enrich$Var1.1<-NULL
CHR.enrich<-CHR.enrich[-c(1,24,25),]
CHR.enrich
CHR.enrich$Var1<-as.numeric(as.character(CHR.enrich$Var1))
CHR.enrich<-CHR.enrich[order(CHR.enrich$Var1),]
CHR.enrich
colnames(CHR.enrich)[1]<-"CHR"

sum<-t(CHR.enrich)
colnames(sum)<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12",
                 "chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22")
sum<-sum[-1,]
rownames(sum)<-c("GS_Probes","All_Probes")
sum<-data.frame(sum)
sum$TOTAL<-data.frame(c(length(res.meta.bonf$CHR),length(res.meta.bacon$CHR)))
colnames(sum)[23]<-c("TOTAL")


out<-matrix(data=NA, nrow=22, ncol=2)
rownames(out)<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12",
                 "chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22")
colnames(out)<-c("OR","p")
for(i in 1:22){
  out[i,]<-c((fisher.test(sum[c(1,2),c(i,23)], alternative = "greater"))$estimate, (fisher.test(sum[c(1,2),c(i,23)], alternative = "greater"))$p.value)
}
write.csv(out, file="ChromosomeenrichmentanalysesFDR.csv")


GS_Probes <-((CHR.enrich[,2]/length(res.meta.bonf$CHR))*100) 
All_Probes <-((CHR.enrich[,3]/length(res.meta.bacon$CHR))*100) 


totalper<-rbind(All_Probes,GS_Probes)
totalper
colnames(totalper)<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12",
                      "chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22")
length(which(results$pval.fixed<0.05))

matrix(c(62, 34318, 928, 402412),nrow = 2)
fisher.test(matrix(c(62, 34318, 928, 402412),nrow = 2), alternative = "greater")



setwd("/mnt/data1/Ehsan/meta_CrossCortical/")
load("RESULTS_CROSS.CORTICAL.META.Rda")

res.meta.bonf<-results[which(results$bonf.FE<0.05),]

sig.cpg<-rownames(res.meta.bonf)
all.cpg<-rownames(results)

go <- crystalmeth(sig.cpg = sig.cpg, 
                  all.cpg = all.cpg,
                  collection = "GO",array.type = "450K" )

write.table(go[which(go$P.DE<0.05),], file="GO.220.CC.txt", sep="\t", col.names=NA)



res.meta<-results[order(results$bonf.FE),]

sig.cpg<-rownames(res.meta[1:1000,])
all.cpg<-rownames(results)

go <- crystalmeth(sig.cpg = sig.cpg, 
                  all.cpg = all.cpg,
                  collection = "GO",array.type = "450K" )

write.table(go[which(go$P.DE<0.05),], file="GO.1000.CC.txt", sep="\t", col.names=NA)


## Perform pathway analysis

This vignette deomstarates how to perform a pathway analysis that controls from the number of probes annotated to each gene. It uses the Illumina UCSC gene annotation, which is derived from the genomic overlap of probes with RefSeq genes or up to 1500bp of the transcription start site of a gene, to annotate probes to genes. Therefore it will ignore any probe that is not annotated to any gene by this method and where probes are annotated to multiple genes all will be included. This method uses a logistic regression approach to test if genes in your specificed test list predicted pathway membership, while controlling for the number of probes annotated to each gene.

Pathways were downloaded from the [GO website](http://geneontology.org/) on the 31st August 2017 and mapped to genes including all parent ontology terms. 


```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

### Load EWAS results 


```{r,eval=FALSE}
setwd("/mnt/data1/Ehsan/meta_CrossCortical/")
load("RESULTS_CROSS.CORTICAL.META.Rda")
res<-results

uniqueAnno<-function(row){
  if(row != ""){
    return(paste(unique(unlist(strsplit(row, "\\;"))), collapse = ";"))
  } else {
    return(row)
  }
}

res$UCSC_RefGene_Name<-as.character(res$UCSC_RefGene_Name)
## reduce gene annotation to unique genes
res$UCSC_RefGene_Name<-unlist(lapply(res$UCSC_RefGene_Name, uniqueAnno))

## exclude sites with no gene anno
res<-res[-which(res$UCSC_RefGene_Name == ""),]

index<-which(res$bonf.FE < 0.05) ## this is a vector of indexes for the DMPs to include in your test list.

gene_test<-table(unlist(strsplit(res$UCSC_RefGene_Name[index], "\\;")))
gene_size<-table(unlist(strsplit(res$UCSC_RefGene_Name, "\\;")))


backgroundGenes<-names(gene_size)
extractGoTerms<-function(backgroundGenes, EPIC = TRUE){
  
  if(EPIC){
    gene_go<-gene_go.epic
  } else {
    gene_go<-gene_go.450K
  }
  backgroundGenes<-intersect(bg_genes, gene_go[,1])
  bg_gene_go<-gene_go[match(backgroundGenes, gene_go[,1]),]
  
  return(bg_gene_go)
}

bg_gene_go<-extractGoTerms(backgroundGenes, EPIC = FALSE)

minPathwaySize<-9 ## adjust if required; set to 0 if you wish to analyse all terms regardless of size.
maxPathwaySize<-2000 ## adjust if required ; set to 10000000 if you wish to analyse all terms regardless of size.

termCount<-table(unlist(strsplit(as.character(bg_gene_go[,2]), "\\|")))
terms<-names(termCount)[which(termCount > minPathwaySize & termCount <= maxPathwaySize)]

## create vector that contains the number of probes per gene that were tested
gene_size<-gene_size[bg_gene_go[,1]]
gene_size<-as.numeric(gene_size)

```
There are two ways you can define your test vector. 

1) provide the number of probes annotated to that gene that are in your test list. This is demonstrated here:
  
  ```{r,eval=FALSE}


gene_test<-gene_test[bg_gene_go[,1]]
gene_test[is.na(gene_test)]<-0

```

2) treat all genes in the test list equally regardless of how many probes were annotated to them, essentially this is a binary variable of the first option. This is demonstrated here:
  
  ```{r,eval=FALSE}


gene_test_ind<-gene_test
gene_test_ind[which(gene_test_ind > 0)]<-1


```

The analyses below use gene_test (i.e. the first option), you will need to swap gene_test_ind in place of gene_test if you wish to run the second option. 

If you don't have many DMPs they this is unlikely to make too much of a difference.


To make this more efficient set up parallel processes using the [doParallel package](https://cran.r-project.org/web/packages/doParallel/vignettes/gettingstartedParallel.pdf).


```{r,eval=FALSE}

library(doParallel)
cl<-makeCluster(30) ## adjust the number in the bracklets to use less or more resources depending on your system limitations. 
registerDoParallel(cl)


```

```{r,eval=FALSE}


r<-foreach(i=1:length(terms), .combine=rbind) %dopar%{
pathwayAnalysis(terms[i], gene_test, gene_size, bg_gene_go)
}
colnames(r)<-c("ID", "Name", "Type", "nProbesinPathway", "nGenesinPathway", "nTestListProbesinPathway",  "nTestListGenesinPathway", "P:GenesinTestList", "OR", "P:GeneSize", "Beta:GeneSize","GenesinTestListAndPathway")

write.csv(r, "") ## put appropriate filename between ""



```

NB nProbesinPathway will only be correct if gene_test is a count of the number of probes per gene in test list.

### Identify independent pathways	

Given the hierarchical structure of the [Gene Ontology](http://www.geneontology.org/) annotations, many associated terms are not independent and are associated by virtue of their overlapping membership. To reduce the results to those that represent independent associations, the code below groups terms where the significant enrichment is explained by the overlap with a more significant term.  This was achieved by taking the most significant pathway, and retesting all remaining significant pathways while controlling additionally for the best term. If the test genes no longer predicted the pathway, the term is said to be explained by the more significant pathway, and hence these pathways are grouped together. This algorithm is repeated, taking the next most significant term, until all pathways are considered as the most significant or found to be explained by a more significant term.

```{r,eval=FALSE}

## reformat r so variables are numeric

r<-as.data.frame(r, row.names = NULL, stringsAsFactors = FALSE)
for(i in c(4:11)){
mode(r[,i])<-"numeric"
}

p1.thres<-0.05 ## what threshold to define pathways as significantly enriched.
p2.thres<-0.05 ## what threshold to define pathways as no longer significant.

## sort results by order of significance
r<-r[order(r[,8]),]
r<-r[which(r[,8] < p1.thres & r[,9] > 0),]

output<-c()
while(!is.null(r)){
if(nrow(r) > 1){
### for all terms repeat analysis controlling for most significant terms
best_term<-vector(length = nrow(bg_gene_go))
best_term[grep(r[1,1], bg_gene_go[,2])]<-1
merge.id<-c()
merge.name<-c()
remove<-c()
for(j in 2:nrow(r)){
pathway<-vector(length = nrow(bg_gene_go))
pathway[grep(r[j,1], bg_gene_go[,2])]<-1
## automatically merge identical pathways
if(!identical(gene_test, best_term)){
model<-glm(pathway ~ gene_test + gene_size + best_term)
if(summary(model)$coefficients["gene_test", "Pr(>|t|)"] > p2.thres){
merge.id<-append(merge.id, r[j,1])
merge.name<-append(merge.name, r[j,2])
remove<-append(remove, j)
}
} else {
merge.id<-append(merge.id, r[j,1])
merge.name<-append(merge.name, r[j,2])
remove<-append(remove, j)
}
}
merge.id<-paste(unlist(merge.id), collapse = "|")
merge.name<-paste(unlist(merge.name), collapse = "|")
output<-rbind(output, c(r[1,], merge.id, merge.name))
r<-r[-c(1, remove),]
} else {
output<-rbind(output, c(r, "", ""))
r<-NULL
}
}

write.csv(output, "")	## put appropriate filename between ""



betas.rm.pfc[rownames(LOC1),]
)
betas.loc2<-cbind(betas.lbb1.pfc[rownames(LOC2),],
betas.ms.pfc[rownames(LOC2),],
betas.rm.pfc[rownames(LOC2),]
)
betas.loc6<-cbind(betas.lbb1.pfc[rownames(LOC6),],
betas.ms.pfc[rownames(LOC6),],
betas.rm.pfc[rownames(LOC6),]
)
betas.loc6<-na.omit(betas.loc6)
betas.loc7<-cbind(betas.lbb1.pfc[rownames(LOC7),],
betas.ms.pfc[rownames(LOC7),],
betas.rm.pfc[rownames(LOC7),]
)
betas.loc9<-cbind(betas.lbb1.pfc[rownames(LOC9),],
betas.ms.pfc[rownames(LOC9),],
betas.rm.pfc[rownames(LOC9),]
)
betas.loc10<-cbind(betas.lbb1.pfc[rownames(LOC10),],
betas.ms.pfc[rownames(LOC10),],
betas.rm.pfc[rownames(LOC10),]
)
betas.loc11<-cbind(betas.lbb1.pfc[rownames(LOC11),],
betas.ms.pfc[rownames(LOC11),],
betas.rm.pfc[rownames(LOC11),]
)
betas.loc13<-cbind(betas.lbb1.pfc[rownames(LOC13),],
betas.ms.pfc[rownames(LOC13),],
betas.rm.pfc[rownames(LOC13),]
)
betas.loc13<-na.omit(betas.loc13)
p.val<-as.numeric(LOC1$pval.fixed)
browns.Loc1.cc<-empiricalBrownsMethod(betas.loc1, p.val, extra_info=T)
p.val<-as.numeric(LOC2$pval.fixed)
browns.Loc2.cc<-empiricalBrownsMethod(betas.loc2, p.val, extra_info=T)
p.val<-as.numeric((LOC6[rownames(betas.loc6),])$pval.fixed)
browns.Loc6.cc<-empiricalBrownsMethod(betas.loc6, p.val, extra_info=T)
p.val<-as.numeric(LOC7$pval.fixed)
browns.Loc7.cc<-empiricalBrownsMethod(betas.loc7, p.val, extra_info=T)
p.val<-as.numeric(LOC9$pval.fixed)
browns.Loc9.cc<-empiricalBrownsMethod(betas.loc9, p.val, extra_info=T)
p.val<-as.numeric(LOC10$pval.fixed)
browns.Loc10.cc<-empiricalBrownsMethod(betas.loc10, p.val, extra_info=T)
p.val<-as.numeric(LOC11$pval.fixed)
browns.Loc11.cc<-empiricalBrownsMethod(betas.loc11, p.val, extra_info=T)
p.val<-as.numeric((LOC13[rownames(betas.loc13),])$pval.fixed)
browns.Loc13.cc<-empiricalBrownsMethod(betas.loc13, p.val, extra_info=T)
res.pfc<-matrix(nrow=8, ncol=3)
colnames(res.cc)<-c("P_Brown","P_Fisher","nProbes")
res.pfc[1,]<- c(browns.Loc1.cc$P_test, browns.Loc1.cc$P_Fisher, nrow(betas.loc1))
res.pfc[2,]<- c(browns.Loc2.cc$P_test, browns.Loc2.cc$P_Fisher, nrow(betas.loc2))
res.pfc[3,]<- c(browns.Loc6.cc$P_test, browns.Loc6.cc$P_Fisher, nrow(betas.loc6))
res.pfc[4,]<- c(browns.Loc7.cc$P_test, browns.Loc7.cc$P_Fisher, nrow(betas.loc7))
res.pfc[5,]<- c(browns.Loc9.cc$P_test, browns.Loc9.cc$P_Fisher, nrow(betas.loc9))
res.pfc[6,]<- c(browns.Loc10.cc$P_test, browns.Loc10.cc$P_Fisher, nrow(betas.loc10))
res.pfc[7,]<- c(browns.Loc11.cc$P_test, browns.Loc11.cc$P_Fisher, nrow(betas.loc11))
res.pfc[8,]<- c(browns.Loc13.cc$P_test, browns.Loc13.cc$P_Fisher, nrow(betas.loc13))
View(res.pfc)
colnames(res.pfc)<-c("P_Brown","P_Fisher","nProbes")
res.GWAS.enrichment<-cbind(res.cc,res.pfc, res.ec, res.tg)
write.table(res.GWAS.enrichment, file="res.GWAS.enrichment.txt", sep="\t", col.names=NA)
getwd()
min(LOC1$pval.fixed)
res.cc<-matrix(nrow=8, ncol=4)
colnames(res.cc)<-c("P_Brown","P_Fisher","nProbes")
res.cc[1,]<- c(min(LOC1$pval.fixed),browns.Loc1.cc$P_test, browns.Loc1.cc$P_Fisher, nrow(betas.loc1))
res.cc[2,]<- c(min(LOC2$pval.fixed), browns.Loc2.cc$P_Fisher, nrow(betas.loc2))
res.cc[3,]<- c(min(LOC6$pval.fixed),browns.Loc6.cc$P_test, browns.Loc6.cc$P_Fisher, nrow(betas.loc6))
res.cc[4,]<- c(min(LOC7$pval.fixed),browns.Loc7.cc$P_test, browns.Loc7.cc$P_Fisher, nrow(betas.loc7))
res.cc[5,]<- c(min(LOC9$pval.fixed),browns.Loc9.cc$P_test, browns.Loc9.cc$P_Fisher, nrow(betas.loc9))
res.cc[6,]<- c(min(LOC10$pval.fixed),browns.Loc10.cc$P_test, browns.Loc10.cc$P_Fisher, nrow(betas.loc10))
res.cc[7,]<- c(min(LOC11$pval.fixed),browns.Loc11.cc$P_test, browns.Loc11.cc$P_Fisher, nrow(betas.loc11))
res.cc[8,]<- c(min(LOC13$pval.fixed),browns.Loc13.cc$P_test, browns.Loc13.cc$P_Fisher, nrow(betas.loc13))
res.cc<-matrix(nrow=8, ncol=4)
colnames(res.cc)<-c("minP","P_Brown","P_Fisher","nProbes")
res.cc[1,]<- c(min(LOC1$pval.fixed),browns.Loc1.cc$P_test, browns.Loc1.cc$P_Fisher, nrow(betas.loc1))
res.cc[2,]<- c(min(LOC2$pval.fixed), browns.Loc2.cc$P_Fisher, nrow(betas.loc2))
res.cc[3,]<- c(min(LOC6$pval.fixed),browns.Loc6.cc$P_test, browns.Loc6.cc$P_Fisher, nrow(betas.loc6))
res.cc[4,]<- c(min(LOC7$pval.fixed),browns.Loc7.cc$P_test, browns.Loc7.cc$P_Fisher, nrow(betas.loc7))
res.cc[5,]<- c(min(LOC9$pval.fixed),browns.Loc9.cc$P_test, browns.Loc9.cc$P_Fisher, nrow(betas.loc9))
res.cc[6,]<- c(min(LOC10$pval.fixed),browns.Loc10.cc$P_test, browns.Loc10.cc$P_Fisher, nrow(betas.loc10))
res.cc[7,]<- c(min(LOC11$pval.fixed),browns.Loc11.cc$P_test, browns.Loc11.cc$P_Fisher, nrow(betas.loc11))
res.cc[8,]<- c(min(LOC13$pval.fixed),browns.Loc13.cc$P_test, browns.Loc13.cc$P_Fisher, nrow(betas.loc13))
res.cc<-matrix(nrow=8, ncol=4)
colnames(res.cc)<-c("minP","P_Brown","P_Fisher","nProbes")
res.cc[1,]<- c(min(LOC1$pval.fixed),browns.Loc1.cc$P_test, browns.Loc1.cc$P_Fisher, nrow(betas.loc1))
res.cc[2,]<- c(min(LOC2$pval.fixed),browns.Loc2.cc$P_test, browns.Loc2.cc$P_Fisher, nrow(betas.loc2))
res.cc[3,]<- c(min(LOC6$pval.fixed),browns.Loc6.cc$P_test, browns.Loc6.cc$P_Fisher, nrow(betas.loc6))
res.cc[4,]<- c(min(LOC7$pval.fixed),browns.Loc7.cc$P_test, browns.Loc7.cc$P_Fisher, nrow(betas.loc7))
res.cc[5,]<- c(min(LOC9$pval.fixed),browns.Loc9.cc$P_test, browns.Loc9.cc$P_Fisher, nrow(betas.loc9))
res.cc[6,]<- c(min(LOC10$pval.fixed),browns.Loc10.cc$P_test, browns.Loc10.cc$P_Fisher, nrow(betas.loc10))
res.cc[7,]<- c(min(LOC11$pval.fixed),browns.Loc11.cc$P_test, browns.Loc11.cc$P_Fisher, nrow(betas.loc11))
res.cc[8,]<- c(min(LOC13$pval.fixed),browns.Loc13.cc$P_test, browns.Loc13.cc$P_Fisher, nrow(betas.loc13))
LOC1$PFC.P_Fixed
betas.loc1<-cbind(betas.az1.tg[rownames(LOC1),],
betas.az2.tg[rownames(LOC1),],
betas.lbb1.tg[rownames(LOC1),],
betas.ms.tg[rownames(LOC1),],
betas.lbb1.ec[rownames(LOC1),],
betas.lbb2.ec[rownames(LOC1),],
betas.lbb1.pfc[rownames(LOC1),],
betas.ms.pfc[rownames(LOC1),],
betas.rm.pfc[rownames(LOC1),]
)
betas.loc2<-cbind(betas.az1.tg[rownames(LOC2),],
betas.az2.tg[rownames(LOC2),],
betas.lbb1.tg[rownames(LOC2),],
betas.ms.tg[rownames(LOC2),],
betas.lbb1.ec[rownames(LOC2),],
betas.lbb2.ec[rownames(LOC2),],
betas.lbb1.pfc[rownames(LOC2),],
betas.ms.pfc[rownames(LOC2),],
betas.rm.pfc[rownames(LOC2),]
)
betas.loc6<-cbind(betas.az1.tg[rownames(LOC6),],
betas.az2.tg[rownames(LOC6),],
betas.lbb1.tg[rownames(LOC6),],
betas.ms.tg[rownames(LOC6),],
betas.lbb1.ec[rownames(LOC6),],
betas.lbb2.ec[rownames(LOC6),],
betas.lbb1.pfc[rownames(LOC6),],
betas.ms.pfc[rownames(LOC6),],
betas.rm.pfc[rownames(LOC6),]
)
betas.loc6<-na.omit(betas.loc6)
betas.loc7<-cbind(betas.az1.tg[rownames(LOC7),],
betas.az2.tg[rownames(LOC7),],
betas.lbb1.tg[rownames(LOC7),],
betas.ms.tg[rownames(LOC7),],
betas.lbb1.ec[rownames(LOC7),],
betas.lbb2.ec[rownames(LOC7),],
betas.lbb1.pfc[rownames(LOC7),],
betas.ms.pfc[rownames(LOC7),],
betas.rm.pfc[rownames(LOC7),]
)
betas.loc9<-cbind(betas.az1.tg[rownames(LOC9),],
betas.az2.tg[rownames(LOC9),],
betas.lbb1.tg[rownames(LOC9),],
betas.ms.tg[rownames(LOC9),],
betas.lbb1.ec[rownames(LOC9),],
betas.lbb2.ec[rownames(LOC9),],
betas.lbb1.pfc[rownames(LOC9),],
betas.ms.pfc[rownames(LOC9),],
betas.rm.pfc[rownames(LOC9),]
)
betas.loc10<-cbind(betas.az1.tg[rownames(LOC10),],
betas.az2.tg[rownames(LOC10),],
betas.lbb1.tg[rownames(LOC10),],
betas.ms.tg[rownames(LOC10),],
betas.lbb1.ec[rownames(LOC10),],
betas.lbb2.ec[rownames(LOC10),],
betas.lbb1.pfc[rownames(LOC10),],
betas.ms.pfc[rownames(LOC10),],
betas.rm.pfc[rownames(LOC10),]
)
betas.loc11<-cbind(betas.az1.tg[rownames(LOC11),],
betas.az2.tg[rownames(LOC11),],
betas.lbb1.tg[rownames(LOC11),],
betas.ms.tg[rownames(LOC11),],
betas.lbb1.ec[rownames(LOC11),],
betas.lbb2.ec[rownames(LOC11),],
betas.lbb1.pfc[rownames(LOC11),],
betas.ms.pfc[rownames(LOC11),],
betas.rm.pfc[rownames(LOC11),]
)
betas.loc13<-cbind(betas.az1.tg[rownames(LOC13),],
betas.az2.tg[rownames(LOC13),],
betas.lbb1.tg[rownames(LOC13),],
betas.ms.tg[rownames(LOC13),],
betas.lbb1.ec[rownames(LOC13),],
betas.lbb2.ec[rownames(LOC13),],
betas.lbb1.pfc[rownames(LOC13),],
betas.ms.pfc[rownames(LOC13),],
betas.rm.pfc[rownames(LOC13),]
)
betas.loc13<-na.omit(betas.loc13)
p.val<-as.numeric(LOC1$pval.fixed)
browns.Loc1.cc<-empiricalBrownsMethod(betas.loc1, p.val, extra_info=T)
p.val<-as.numeric(LOC2$pval.fixed)
browns.Loc2.cc<-empiricalBrownsMethod(betas.loc2, p.val, extra_info=T)
p.val<-as.numeric((LOC6[rownames(betas.loc6),])$pval.fixed)
browns.Loc6.cc<-empiricalBrownsMethod(betas.loc6, p.val, extra_info=T)
p.val<-as.numeric(LOC7$pval.fixed)
browns.Loc7.cc<-empiricalBrownsMethod(betas.loc7, p.val, extra_info=T)
p.val<-as.numeric(LOC9$pval.fixed)
browns.Loc9.cc<-empiricalBrownsMethod(betas.loc9, p.val, extra_info=T)
p.val<-as.numeric(LOC10$pval.fixed)
browns.Loc10.cc<-empiricalBrownsMethod(betas.loc10, p.val, extra_info=T)
p.val<-as.numeric(LOC11$pval.fixed)
browns.Loc11.cc<-empiricalBrownsMethod(betas.loc11, p.val, extra_info=T)
p.val<-as.numeric((LOC13[rownames(betas.loc13),])$pval.fixed)
browns.Loc13.cc<-empiricalBrownsMethod(betas.loc13, p.val, extra_info=T)
res.cc<-matrix(nrow=8, ncol=4)
colnames(res.cc)<-c("minP","P_Brown","P_Fisher","nProbes")
res.cc[1,]<- c(min(LOC1$pval.fixed),browns.Loc1.cc$P_test, browns.Loc1.cc$P_Fisher, nrow(betas.loc1))
res.cc[2,]<- c(min(LOC2$pval.fixed),browns.Loc2.cc$P_test, browns.Loc2.cc$P_Fisher, nrow(betas.loc2))
res.cc[3,]<- c(min(LOC6$pval.fixed),browns.Loc6.cc$P_test, browns.Loc6.cc$P_Fisher, nrow(betas.loc6))
res.cc[4,]<- c(min(LOC7$pval.fixed),browns.Loc7.cc$P_test, browns.Loc7.cc$P_Fisher, nrow(betas.loc7))
res.cc[5,]<- c(min(LOC9$pval.fixed),browns.Loc9.cc$P_test, browns.Loc9.cc$P_Fisher, nrow(betas.loc9))
res.cc[6,]<- c(min(LOC10$pval.fixed),browns.Loc10.cc$P_test, browns.Loc10.cc$P_Fisher, nrow(betas.loc10))
res.cc[7,]<- c(min(LOC11$pval.fixed),browns.Loc11.cc$P_test, browns.Loc11.cc$P_Fisher, nrow(betas.loc11))
res.cc[8,]<- c(min(LOC13$pval.fixed),browns.Loc13.cc$P_test, browns.Loc13.cc$P_Fisher, nrow(betas.loc13))
#########################################
#########################################
#########################################
betas.loc1<-cbind(betas.az1.tg[rownames(LOC1),],
betas.az2.tg[rownames(LOC1),],
betas.lbb1.tg[rownames(LOC1),],
betas.ms.tg[rownames(LOC1),]
)
betas.loc2<-cbind(betas.az1.tg[rownames(LOC2),],
betas.az2.tg[rownames(LOC2),],
betas.lbb1.tg[rownames(LOC2),],
betas.ms.tg[rownames(LOC2),]
)
betas.loc6<-cbind(betas.az1.tg[rownames(LOC6),],
betas.az2.tg[rownames(LOC6),],
betas.lbb1.tg[rownames(LOC6),],
betas.ms.tg[rownames(LOC6),]
)
betas.loc6<-na.omit(betas.loc6)
betas.loc7<-cbind(betas.az1.tg[rownames(LOC7),],
betas.az2.tg[rownames(LOC7),],
betas.lbb1.tg[rownames(LOC7),],
betas.ms.tg[rownames(LOC7),]
)
betas.loc9<-cbind(betas.az1.tg[rownames(LOC9),],
betas.az2.tg[rownames(LOC9),],
betas.lbb1.tg[rownames(LOC9),],
betas.ms.tg[rownames(LOC9),]
)
betas.loc10<-cbind(betas.az1.tg[rownames(LOC10),],
betas.az2.tg[rownames(LOC10),],
betas.lbb1.tg[rownames(LOC10),],
betas.ms.tg[rownames(LOC10),]
)
betas.loc11<-cbind(betas.az1.tg[rownames(LOC11),],
betas.az2.tg[rownames(LOC11),],
betas.lbb1.tg[rownames(LOC11),],
betas.ms.tg[rownames(LOC11),]
)
betas.loc13<-cbind(betas.az1.tg[rownames(LOC13),],
betas.az2.tg[rownames(LOC13),],
betas.lbb1.tg[rownames(LOC13),],
betas.ms.tg[rownames(LOC13),]
)
betas.loc13<-na.omit(betas.loc13)
p.val<-as.numeric(LOC1$TG.P_Fixed)
browns.Loc1.cc<-empiricalBrownsMethod(betas.loc1, p.val, extra_info=T)
p.val<-as.numeric(LOC2$TG.P_Fixed)
browns.Loc2.cc<-empiricalBrownsMethod(betas.loc2, p.val, extra_info=T)
p.val<-as.numeric((LOC6[rownames(betas.loc6),])$TG.P_Fixed)
browns.Loc6.cc<-empiricalBrownsMethod(betas.loc6, p.val, extra_info=T)
p.val<-as.numeric(LOC7$TG.P_Fixed)
browns.Loc7.cc<-empiricalBrownsMethod(betas.loc7, p.val, extra_info=T)
p.val<-as.numeric(LOC9$TG.P_Fixed)
browns.Loc9.cc<-empiricalBrownsMethod(betas.loc9, p.val, extra_info=T)
p.val<-as.numeric(LOC10$TG.P_Fixed)
browns.Loc10.cc<-empiricalBrownsMethod(betas.loc10, p.val, extra_info=T)
p.val<-as.numeric(LOC11$TG.P_Fixed)
browns.Loc11.cc<-empiricalBrownsMethod(betas.loc11, p.val, extra_info=T)
p.val<-as.numeric((LOC13[rownames(betas.loc13),])$TG.P_Fixed)
browns.Loc13.cc<-empiricalBrownsMethod(betas.loc13, p.val, extra_info=T)
res.tg<-matrix(nrow=8, ncol=4)
colnames(res.tg)<-c("minP","P_Brown","P_Fisher","nProbes")
res.tg[1,]<- c(min(LOC1$TG.P_Fixed),browns.Loc1.cc$P_test, browns.Loc1.cc$P_Fisher, nrow(betas.loc1))
res.tg[2,]<- c(min(LOC2$TG.P_Fixed),browns.Loc2.cc$P_test, browns.Loc2.cc$P_Fisher, nrow(betas.loc2))
res.tg[3,]<- c(min(LOC6$TG.P_Fixed),browns.Loc6.cc$P_test, browns.Loc6.cc$P_Fisher, nrow(betas.loc6))
res.tg[4,]<- c(min(LOC7$TG.P_Fixed),browns.Loc7.cc$P_test, browns.Loc7.cc$P_Fisher, nrow(betas.loc7))
res.tg[5,]<- c(min(LOC9$TG.P_Fixed),browns.Loc9.cc$P_test, browns.Loc9.cc$P_Fisher, nrow(betas.loc9))
res.tg[6,]<- c(min(LOC10$TG.P_Fixed),browns.Loc10.cc$P_test, browns.Loc10.cc$P_Fisher, nrow(betas.loc10))
res.tg[7,]<- c(min(LOC11$TG.P_Fixed),browns.Loc11.cc$P_test, browns.Loc11.cc$P_Fisher, nrow(betas.loc11))
res.tg[8,]<- c(min(LOC13$TG.P_Fixed),browns.Loc13.cc$P_test, browns.Loc13.cc$P_Fisher, nrow(betas.loc13))
############################################################################
############################################################################
#########################################
betas.loc1<-cbind(betas.lbb1.ec[rownames(LOC1),],
betas.lbb2.ec[rownames(LOC1),]
)
betas.loc2<-cbind(betas.lbb1.ec[rownames(LOC2),],
betas.lbb2.ec[rownames(LOC2),]
)
betas.loc6<-cbind(betas.lbb1.ec[rownames(LOC6),],
betas.lbb2.ec[rownames(LOC6),]
)
betas.loc6<-na.omit(betas.loc6)
betas.loc7<-cbind(betas.lbb1.ec[rownames(LOC7),],
betas.lbb2.ec[rownames(LOC7),]
)
betas.loc9<-cbind(betas.lbb1.ec[rownames(LOC9),],
betas.lbb2.ec[rownames(LOC9),]
)
betas.loc10<-cbind(betas.lbb1.ec[rownames(LOC10),],
betas.lbb2.ec[rownames(LOC10),]
)
betas.loc11<-cbind(betas.lbb1.ec[rownames(LOC11),],
betas.lbb2.ec[rownames(LOC11),]
)
betas.loc13<-cbind(betas.lbb2.ec[rownames(LOC13),],
betas.lbb1.pfc[rownames(LOC13),]
)
betas.loc13<-na.omit(betas.loc13)
p.val<-as.numeric(LOC1$EC.P_Fixed)
browns.Loc1.cc<-empiricalBrownsMethod(betas.loc1, p.val, extra_info=T)
p.val<-as.numeric(LOC2$EC.P_Fixed)
browns.Loc2.cc<-empiricalBrownsMethod(betas.loc2, p.val, extra_info=T)
p.val<-as.numeric((LOC6[rownames(betas.loc6),])$EC.P_Fixed)
browns.Loc6.cc<-empiricalBrownsMethod(betas.loc6, p.val, extra_info=T)
p.val<-as.numeric(LOC7$EC.P_Fixed)
browns.Loc7.cc<-empiricalBrownsMethod(betas.loc7, p.val, extra_info=T)
p.val<-as.numeric(LOC9$EC.P_Fixed)
browns.Loc9.cc<-empiricalBrownsMethod(betas.loc9, p.val, extra_info=T)
p.val<-as.numeric(LOC10$EC.P_Fixed)
browns.Loc10.cc<-empiricalBrownsMethod(betas.loc10, p.val, extra_info=T)
p.val<-as.numeric(LOC11$EC.P_Fixed)
browns.Loc11.cc<-empiricalBrownsMethod(betas.loc11, p.val, extra_info=T)
p.val<-as.numeric((LOC13[rownames(betas.loc13),])$EC.P_Fixed)
browns.Loc13.cc<-empiricalBrownsMethod(betas.loc13, p.val, extra_info=T)
res.ec<-matrix(nrow=8, ncol=4)
colnames(res.ec)<-c("minP","P_Brown","P_Fisher","nProbes")
res.ec[1,]<- c(min(LOC1$EC.P_Fixed),browns.Loc1.cc$P_test, browns.Loc1.cc$P_Fisher, nrow(betas.loc1))
res.ec[2,]<- c(min(LOC2$EC.P_Fixed),browns.Loc2.cc$P_test, browns.Loc2.cc$P_Fisher, nrow(betas.loc2))
res.ec[3,]<- c(min(LOC6$EC.P_Fixed),browns.Loc6.cc$P_test, browns.Loc6.cc$P_Fisher, nrow(betas.loc6))
res.ec[4,]<- c(min(LOC7$EC.P_Fixed),browns.Loc7.cc$P_test, browns.Loc7.cc$P_Fisher, nrow(betas.loc7))
res.ec[5,]<- c(min(LOC9$EC.P_Fixed),browns.Loc9.cc$P_test, browns.Loc9.cc$P_Fisher, nrow(betas.loc9))
res.ec[6,]<- c(min(LOC10$EC.P_Fixed),browns.Loc10.cc$P_test, browns.Loc10.cc$P_Fisher, nrow(betas.loc10))
res.ec[7,]<- c(min(LOC11$EC.P_Fixed),browns.Loc11.cc$P_test, browns.Loc11.cc$P_Fisher, nrow(betas.loc11))
res.ec[8,]<- c(min(LOC13$EC.P_Fixed),browns.Loc13.cc$P_test, browns.Loc13.cc$P_Fisher, nrow(betas.loc13))
#########################################
#########################################
#########################################
#########################################
betas.loc1<-cbind(betas.lbb1.pfc[rownames(LOC1),],
betas.ms.pfc[rownames(LOC1),],
betas.rm.pfc[rownames(LOC1),]
)
betas.loc2<-cbind(betas.lbb1.pfc[rownames(LOC2),],
betas.ms.pfc[rownames(LOC2),],
betas.rm.pfc[rownames(LOC2),]
)
betas.loc6<-cbind(betas.lbb1.pfc[rownames(LOC6),],
betas.ms.pfc[rownames(LOC6),],
betas.rm.pfc[rownames(LOC6),]
)
betas.loc6<-na.omit(betas.loc6)
betas.loc7<-cbind(betas.lbb1.pfc[rownames(LOC7),],
betas.ms.pfc[rownames(LOC7),],
betas.rm.pfc[rownames(LOC7),]
)
betas.loc9<-cbind(betas.lbb1.pfc[rownames(LOC9),],
betas.ms.pfc[rownames(LOC9),],
betas.rm.pfc[rownames(LOC9),]
)
betas.loc10<-cbind(betas.lbb1.pfc[rownames(LOC10),],
betas.ms.pfc[rownames(LOC10),],
betas.rm.pfc[rownames(LOC10),]
)
betas.loc11<-cbind(betas.lbb1.pfc[rownames(LOC11),],
betas.ms.pfc[rownames(LOC11),],
betas.rm.pfc[rownames(LOC11),]
)
betas.loc13<-cbind(betas.lbb1.pfc[rownames(LOC13),],
betas.ms.pfc[rownames(LOC13),],
betas.rm.pfc[rownames(LOC13),]
)
betas.loc13<-na.omit(betas.loc13)
p.val<-as.numeric(LOC1$PFC.P_Fixed)
browns.Loc1.cc<-empiricalBrownsMethod(betas.loc1, p.val, extra_info=T)
p.val<-as.numeric(LOC2$PFC.P_Fixed)
browns.Loc2.cc<-empiricalBrownsMethod(betas.loc2, p.val, extra_info=T)
p.val<-as.numeric((LOC6[rownames(betas.loc6),])$PFC.P_Fixed)
browns.Loc6.cc<-empiricalBrownsMethod(betas.loc6, p.val, extra_info=T)
p.val<-as.numeric(LOC7$PFC.P_Fixed)
browns.Loc7.cc<-empiricalBrownsMethod(betas.loc7, p.val, extra_info=T)
p.val<-as.numeric(LOC9$PFC.P_Fixed)
browns.Loc9.cc<-empiricalBrownsMethod(betas.loc9, p.val, extra_info=T)
p.val<-as.numeric(LOC10$PFC.P_Fixed)
browns.Loc10.cc<-empiricalBrownsMethod(betas.loc10, p.val, extra_info=T)
p.val<-as.numeric(LOC11$PFC.P_Fixed)
browns.Loc11.cc<-empiricalBrownsMethod(betas.loc11, p.val, extra_info=T)
p.val<-as.numeric((LOC13[rownames(betas.loc13),])$PFC.P_Fixed)
browns.Loc13.cc<-empiricalBrownsMethod(betas.loc13, p.val, extra_info=T)
res.pfc<-matrix(nrow=8, ncol=4)
colnames(res.pfc)<-c("minP","P_Brown","P_Fisher","nProbes")
res.pfc[1,]<- c(min(LOC1$PFC.P_Fixed),browns.Loc1.cc$P_test, browns.Loc1.cc$P_Fisher, nrow(betas.loc1))
res.pfc[2,]<- c(min(LOC2$PFC.P_Fixed),browns.Loc2.cc$P_test, browns.Loc2.cc$P_Fisher, nrow(betas.loc2))
res.pfc[3,]<- c(min(LOC6$PFC.P_Fixed),browns.Loc6.cc$P_test, browns.Loc6.cc$P_Fisher, nrow(betas.loc6))
res.pfc[4,]<- c(min(LOC7$PFC.P_Fixed),browns.Loc7.cc$P_test, browns.Loc7.cc$P_Fisher, nrow(betas.loc7))
res.pfc[5,]<- c(min(LOC9$PFC.P_Fixed),browns.Loc9.cc$P_test, browns.Loc9.cc$P_Fisher, nrow(betas.loc9))
res.pfc[6,]<- c(min(LOC10$PFC.P_Fixed),browns.Loc10.cc$P_test, browns.Loc10.cc$P_Fisher, nrow(betas.loc10))
res.pfc[7,]<- c(min(LOC11$PFC.P_Fixed),browns.Loc11.cc$P_test, browns.Loc11.cc$P_Fisher, nrow(betas.loc11))
res.pfc[8,]<- c(min(LOC13$PFC.P_Fixed),browns.Loc13.cc$P_test, browns.Loc13.cc$P_Fisher, nrow(betas.loc13))
res.GWAS.enrichment<-cbind(res.cc,res.pfc, res.ec, res.tg)
res.GWAS.enrichment
write.table(res.GWAS.enrichment, file="res.GWAS.enrichment.txt", sep="\t", col.names=NA)

