# Genetic association studies: Pipelines and useful functions

The following repository contains pipelines and functions I found useful while
working on genome-wide association studies. This work includes genotype quality
control, phasing and imputation, genome-wide association studies and analysis
thereof.

## 1. GWAS:
1. plots.r: 
    * Function for manhattan of GWAS results
    * Function for qqplot of GWAS results
1. metaanalysis.r: 
    * Function for estimating the effective number of hypothesis tests conducted, based on the eigenvalues of
    the correlation matrix of the phenotypes (based on [Galwey (2009) Genetic Epidemiology](https://onlinelibrary.wiley.com/doi/abs/10.1002/gepi.20408)).
    * Function for combining univariate association results into a meta-analysis p-values (based on [Boloorma (2014) Plos Genetics](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004198)).
1. relatedness.r:
    * Function filtering related samples from list of samples in cohort with the aim of retaining maximum number of samples in the cohort.

## 2. Genotyping quality control
 
 Coming soon...

## 3. Imputation and phasing
  
  Coming soon...
