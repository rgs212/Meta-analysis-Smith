library(bacon)
library(parallel)
library(nlme)
library(meta)

##-------------------------------------------------------
##Linear regression analyses were performed with respect to Braak stage (modelled as a continuous variable) 
##using residuals and a variable number of SVs for each study until the inflation index (lambda) fell below 1.2
##-------------------------------------------------------
EWAS <- function(x,Braak, sv1,sv2,...){
  res<-lm(x ~ Braak + sv1 + sv2 + ...)
  return(coef(summary(res))[2,])
}
  
##-------------------------------------------------------
##As multiple cortical brain regions were available for the “London 1” and “Mount Sinai” cohorts, 
##a mixed model with nested IDs was performed using the lme function within the nlme package.
##-------------------------------------------------------
EWAS <- function(x, Braak,id, sv1, sv2, ...){
  library(nlme)
  res<-try(nlme::lme(as.numeric(x) ~ Braak + sv1  + sv2 + sv3 + ... , random=~1|id), silent=T)
  if(inherits(res,'try-error')) return(rep(NA,5))
  return(coef(summary(res))[2,])
} 

##-------------------------------------------------------
##Running the EWAS using parallel processing codes
##-------------------------------------------------------
cl<- makeCluster(n)
res.cohort.i<-t(parApply(cl, cohort.i(harmonised methylation data), 1, EWAS,Braak, sv1,sv2,....))

##-------------------------------------------------------
##Calculating inflation index for each EWAS
##-------------------------------------------------------
chisq <- qchisq(1-(res.cohort.i[,""p]),1)
median(chisq)/qchisq(0.5,1)

##-------------------------------------------------------
##Intra-tissue meta-analysis 
##Estimate coefficients and SEs from each EWAS
##-------------------------------------------------------

braaklm <- function( row, braak, Age, Gender, CETs ){
     
   fit <- try (
      lm( row ~ Age + Gender + CETs + braak),
      silent=TRUE
   )
   if(inherits(fit,'try-error')) return(rep(NA,8))
   as.numeric(summary(fit)$coeff[,c(1,2,4)])
} 

results  <- {
      Gender    <- pheno[ colnames(data), 'Gender' ]
      Age    	<- pheno[ colnames(data), 'Age'])
      CETs    	<- pheno[ colnames(data), 'cellprop' ]
      braak 	<- pheno[ colnames(data), 'Braak'   ])
     
      t( apply(data, 1, braaklm, braak, Age, Gender, CETs))
   }

##-------------------------------------------------------
##A fixed-effect inverse variance meta-analysis was performed across all discovery cohorts using the metagen function
##-------------------------------------------------------

res<-matrix(data = NA, nrow = nrow(res.tissue), ncol = 6)
rownames(res)<-rownames(res.tissue)
colnames(res)<-c("TE.fixed", "seTE.fixed", "pval.fixed", "TE.random", "seTE.random", "pval.random", "I2", "1-pchisq(Q, df.Q)")

for(i in 1:nrow(res.tissue)){
out<-metagen(c(res.tissue$es.coh.1[i], res.tissue$es.coh.2[i] ...), c(res.tissue$se.coh.1[i], res.tissue$se.coh.2[i] ...))
res[i,]<-c(out$TE.fixed, out$seTE.fixed, out$pval.fixed, out$TE.random, out$seTE.random, out$pval.random, out$I2, 1-pchisq(out$Q, out$df.Q))
}

##-------------------------------------------------------
##Cross-cortex AD-associated DMPs:
##Estimate coefficients and SEs from each EWAS were extracted and were subjected to bacon to control for bias and inflation
##-------------------------------------------------------
es<-cbind(res.coh.1[,"estimate"],res.coh.2[,"estimate"], ....res.coh.6[,"estimate"])
se<-cbind(res.coh.1[,"se"],res.coh.2[,"se"], ....res.coh.6[,"se"])
identical(rownames(es), rownames(se))

library(BiocParallel)
register(MulticoreParam(32, log=TRUE))
colnames(es) <- colnames(se) <- c("coh.1","coh.2",...,"coh.6")
rownames(es) <- rownames(se) <-rownames(es)

bc <- bacon(NULL, es, se)
es.b<-es(bc)
se.b<-se(bc)

##-------------------------------------------------------
##A fixed-effect inverse variance meta-analysis was performed across all discovery cohorts using the metagen function
##-------------------------------------------------------

for(i in 1:nrow(data)){
  dat<-cbind(es.b,se.b)
  out<-metagen(data[i,1:6],data[i,7:12])
  res.meta[i,]<-c(out$TE.fixed, out$seTE.fixed, out$pval.fixed, out$TE.random, out$seTE.random, out$pval.random, out$I2, 1-pchisq(out$Q, out$df.Q))
}

cl<-makeCluster(32)
res.meta<-t(parApply(cl,data,1,meta))
colnames(res.meta)<-c("TE.fixed", "seTE.fixed", "pval.fixed", "TE.random", "seTE.random", "pval.random", "I2", "1-pchisq(Q, df.Q)")

##-------------------------------------------------------
##Creating a bed file as an input for the comb-p tool
##-------------------------------------------------------
dmr <- data.frame("chrom" = 	paste0("chr", results$CHR), 
                  "start" = 	results[rownames(results), "MAPINFO.1"],
                  "end" 	= 	results[rownames(results), "MAPINFO.1.1"],
                  "pvalue" = 	results$pval.fixed)

dmr<-dmr[order(dmr[,1],dmr[,2]),]
write.table(dmr, file="meta.DMR.bed",sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

##-------------------------------------------------------
##Identification of Differentially methylated regions using the comb-p tool
##-------------------------------------------------------
# comb-p pipeline -c 4 --dist 500 --seed 1.0e-4 --anno hg19 -p  out meta.DMR.bed

##-------------------------------------------------------
## AD GWAS enrichment analysis (Brown's method)
##-------------------------------------------------------
##Creating a matrix of LD regions according to Kunkle et al 2019 

matrix(c("chr1","207679307","207850539","CR1",
  "chr2","127882182","127894615","BIN1",
  "chr2","233976593","233981912","INPP5D",
  "chr6","32395036","32636434","HLA-DRB1",
  "chr6","40706366","41365821","TREM2",
  "chr6","47412916","47628558","CD2AP",
  "chr7","99932049","100190116","NYAP1",
  "chr7","143099107","143109208","EPHA1",
  "chr8","27195121","27238052","PTK2B",
  "chr8","27456253","27468503","CLU",
  "chr10","11703491","11723257","ECHDC3",
  "chr11","47372377","47466790","SPI1",
  "chr11","59856028","60097777","MS4A2",
  "chr11","85670385","85868640","PICALM",
  "chr11","121433926","121461593","SORL1",
  "chr14","53293307","53462216","FERMT2",
  "chr14","92926952","92957176","SLC24A4",
  "chr15","58873555","59120077","ADAM10",
  "chr16","19706199","19867021","IQCK",
  "chr16","79355847","79355847","WWOX",
  "chr17","61499732","61543566","ACE",
  "chr19","1050130","1075979","ABCA7",
  "chr20","54979828","55025377","CASS4",
  "chr21","28146668","28166355","ADAMTS1"),ncol=4,byrow=T)


LD.Probes<-data.frame(paste0("chr",results$CHR), 
                   as.numeric(results$MAPINFO.1),
                   as.numeric(results$MAPINFO.1.1),
                   results$IlmnID)
write.table(LD.Probes, file="LD.probes.bed", sep="\t", col.names= FALSE,row.names=F, quote=FALSE )


##Identification of the 450K probes falling into the defined LD regions
bedtools sort -i LD.probes.bed > LD.probes.sorted.bed
bedtools intersect -a LD.probes.sorted.bed -b LD_AD.sorted_bed -wo > LD.meta.bed

p.val<-as.numeric(Loc.i$pval.fixed)

##Conducting Brown's method
browns.Loc.i<-empiricalBrownsMethod(betas.loc.i, p.val, extra_info=T)


##-------------------------------------------------------
## GO enrichment analysis
##-------------------------------------------------------
#This vignette deomstarates how to perform a pathway analysis that controls from the number of probes annotated to each gene. It uses the Illumina UCSC gene annotation, 
#which is derived from the genomic overlap of probes with RefSeq genes or up to 1500bp of the transcription start site of a gene, to annotate probes to genes. 
#Therefore it will ignore any probe that is not annotated to any gene by this method and where probes are annotated to multiple genes all will be included. 
#This method uses a logistic regression approach to test if genes in your specificed test list predicted pathway membership, while controlling for the number of probes annotated to each gene.
#Pathways were downloaded from the [GO website](http://geneontology.org/) on the 31st August 2017 and mapped to genes including all parent ontology terms. 

### Load the results results 
setwd("")
load("RESULTS")
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

#The analyses below use gene_test (i.e. the first option), you will need to swap gene_test_ind in place of gene_test if you wish to run the second option. 

#If you don't have many DMPs they this is unlikely to make too much of a difference.


#To make this more efficient set up parallel processes using the [doParallel package](https://cran.r-project.org/web/packages/doParallel/vignettes/gettingstartedParallel.pdf).

library(doParallel)
cl<-makeCluster(30) ## adjust the number in the bracklets to use less or more resources depending on your system limitations. 
registerDoParallel(cl)

r<-foreach(i=1:length(terms), .combine=rbind) %dopar%{
  pathwayAnalysis(terms[i], gene_test, gene_size, bg_gene_go)
}
colnames(r)<-c("ID", "Name", "Type", "nProbesinPathway", "nGenesinPathway", "nTestListProbesinPathway",  "nTestListGenesinPathway", "P:GenesinTestList", "OR", "P:GeneSize", "Beta:GeneSize","GenesinTestListAndPathway")

write.csv(r, "") ## put appropriate filename between ""


#NB nProbesinPathway will only be correct if gene_test is a count of the number of probes per gene in test list.

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

##-------------------------------------------------------
## Quantifying variance explained by PES and PRS 
##-------------------------------------------------------

#Added soon

##-------------------------------------------------------
## Developing a classifier to predict Braak pathology
##-------------------------------------------------------

# GLMNET - PENALIZED REGRESSION METHODS TO PREDICT BRAAK STAGE

library(glmnet)
require(caTools)
library(wateRmelon)
library(ggplot2)
library(reshape2)

### LOAD DATA 
setwd("")
load("data.rda") ## INCLUDES data
betas<-data

## READ IN FILE WITH BRAAK STAGE AND CASE CONTROL STATUS
# Pheno file should have rownames as sample IDs/Sentrix
pheno<-read.csv("phenonew.csv")
pheno<-pheno[complete.cases(pheno$CC),]

## match to betas 
betas1<-as.matrix(betas)

## Load 220 probes only
probes<-read.table("meta_sig_probes.txt")
betas1<-betas1[which(rownames(betas1) %in% probes$V1),]
betas1<-betas1[,match(pheno1$Basename, colnames(betas1))]

# ALL SAMPLES

#### SPLIT .75 AND .25 INTO TRAINING AND TESTING DATA, RESPECTIVELY 
set.seed(101) # Set Seed so that same sample can be reproduced in future also
# # Now Selecting 75% of data as sample from total 'n' rows of the data
sample <- sample.int(n = nrow(pheno1), size = floor(.75*nrow(pheno1)), replace = F)
train <- pheno1[sample, ]
test  <- pheno1[-sample, ]

train$status<-ifelse(train$CC==0, "control","case")
train$status<-as.numeric(as.factor(train$status))

test$status<-ifelse(test$CC==0, "control","case")
test$status<-as.numeric(as.factor(test$status))

## MATCH TO BETAS
# training
betasTrain<-betas1[,match(train$Basename, colnames(betas1))]
dim(betasTrain)
# testing
betasTest<-betas1[,match(test$Basename, colnames(betas1))]
dim(betasTest) 

# RUN GLMNET
## need a model matrix when using factors in glmnet 
## to note that status is the factor, so it will be in reference to the case (as it will be alphabetical) so in results <0.5 = case, >0.5 = control
status_matrix <-model.matrix(status ~ ., train) 

# use 10 fold cross validation to estimate the lambda parameter 
# in the training data
# AS IT IS ELAST NET REGRESSION, IT USES BOTH RIDGE REGRESSION (ALPHA=0) AND LASSO REGRESSION (ALPHA=1)
alpha<-0.5
glmnet.Training.CV = cv.glmnet(t(betasTrain), train$status, nfolds=10,alpha=alpha,family="binomial") ## CV = cross-validation
# The definition of the lambda parameter:
lambda.glmnet.Training = glmnet.Training.CV$lambda.min
# Fit the elastic net predictor to the training data
glmnet.Training = glmnet(t(betasTrain),train$status, family="binomial", alpha=0.5, nlambda=100)
#Arrive at an estimate of of DNAmBraakStage
DNAmBraakStageBrainTraining=predict(glmnet.Training,t(betasTrain),type="response",s=lambda.glmnet.Training)
DNAmBraakStageBrainTesting=predict(glmnet.Training,t(betasTest),type="response",s=lambda.glmnet.Training)

#### EXTRACT THE COEFFIEICENTS FOR THE PROBES

tmp_coeffs <- coef(glmnet.Training.CV, s = "lambda.min")
myCoeff<-data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)
dim(myCoeff)
head(myCoeff) # 96

train$StatusPred<-DNAmBraakStageBrainTraining[,1]
train[,1]<-NULL
test$StatusPred<-DNAmBraakStageBrainTesting[,1]
test[,1]<-NULL
