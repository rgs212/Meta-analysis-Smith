<b>Meta-analysis of epigenome-wide association studies in Alzheimer’s disease highlights 220 differentially methylated loci across cortex</b>

<i>Rebecca G. Smith<sup>+</sup>, Ehsan Pishva<sup>+</sup>, Gemma Shireby, Adam R. Smith, Janou A.Y. Roubroeks, Eilis Hannon, Gregory Wheildon, Diego Mastroeni, Gilles Gasparoni, Matthias Riemenschneider, Armin Giese, Andrew J. Sharp, Leonard Schalkwyk, Vahram Haroutunian, Wolfgang Viechtbauer, Daniel L.A. van den Hove, Michael Weedon, Jörn Walter, Paul D. Coleman, David A. Bennett, Philip L. De Jager, Jonathan Mill, Katie Lunnon<sup>*</sup></i>

<sup>+</sup> Joint first author

<sup>*</sup> Corresponding author - k.lunnon@exeter.ac.uk

Epigenome-wide association studies of Alzheimer’s disease have highlighted neuropathology-associated DNA methylation differences, although existing studies have been limited in sample size and utilized different brain regions. Here, we combine data from six methylomic studies of Alzheimer’s disease (N=1,453 unique individuals) to identify differential methylation associated with Braak stage in different brain regions and across cortex. At an experiment-wide significance threshold (P < 1.238 x 10<sup>−7</sup>) we identified 236 CpGs in the prefrontal cortex, 95 CpGs in the temporal gyrus and ten CpGs in the entorhinal cortex, with none in the cerebellum. Our cross-cortex meta-analysis (N=1,408 donors) identified 220 CpGs associated with neuropathology, annotated to 121 genes, of which 96 genes had not been previously reported at experiment-wide significance. Polyepigenic scores derived from these 220 CpGs explain 24.7% of neuropathological variance, whilst polygenic scores accounted for 20.2% of variance in these samples. The meta-analysis summary statistics are available in our online data resource (www.epigenomicslab.com/ad-meta-analysis/).

Paper uploaded to Biorxiv (https://www.biorxiv.org/content/10.1101/2020.02.28.957894v1)/

This page contains scripts for the above paper
1. Data quality control and harmonization.r - QC and harmonization for all samples included in analysis
   - Loading the raw idat files
   - Quality control
      - Detecting samples with extreme intesities using the negative control probes as references
      - Detecting samples with low background to signal ratio
      - detecting the samples with the extreme mean intensity of methylated or unmethylated signals 
      - Detecting the samples with bisulfite conversion efficiency < 80%
      - Detecting a mismatch between reported and predicted sex
      - Uses SNP probes on the array to find genetic correlations between samples
   - WateRmelon - pfilter
   - Quantile normalization
   - Removing the probes with evidence for cross-hybridizition
   - CET
   - Harmonisation
   - Surrogate variables (SVs)armonisation
   
2. Main Analyses.r - All analyses used
   - Linear regression analyses were performed with respect to Braak stage
   - Intra-tissue meta-analysis
   - Fixed-effect inverse variance meta-analysis
   - Cross-cortex AD-associated DMPs
   - Comb-p
   - AD GWAS enrichment analysis
   - GO enrichment analysis
   - Quantifying variance explained by PES and PRS
   - Developing a classifier to predict Braak pathology

#Statistical analysis was performed in R 3.5.2 and Bioconductor 3.8, Python and PLINK 1.9. 
#R Packages used were wateRmelon 1.26.0, minfi 1.28.4, CETS, sva 3.30.1, Meta 4.10.0, 
#nlme 3.1.142, bacon 1.10.1, pROC 1.16.1, glmnet 2.0-18, PRSice 2.2.12. 
#Python package comb-p 33.1.1.

#All programs and packages installed via Anaconda on cluster computing linux system
# Operating System: CentOS Linux 7 (Core)
# CPE OS Name: cpe:/o:centos:centos:7
# Kernel: Linux 3.10.0-693.17.1.el7.x86_64
# Architecture: x86-64
