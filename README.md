# Genetic association studies: Pipelines and useful functions

The following repository contains pipelines and functions I found useful while
working on genome-wide association studies. This work includes links to detailed
genotype quality control with **plinkQC**, reference panel download,
phasing and imputation, genome-wide association studies and analysis
thereof.

All scripts can simply be sourced to access the relevant functions. Alternatively,
the [modules package](https://github.com/klmr/modules) offers flexible and tidy
integration of R source files. Application of functions in this repository via modules in e.g.
[ukbb-fd](https://github.com/HannahVMeyer/ukbb-fd/blob/master/association/association_results.R).


## 1. Genotyping quality control
Genotype quality control for genotyping arrays with
[plinkQC](http://meyer-lab-cshl.github.io/plinkQC/)

## 2. Imputation and phasing
Phasing and imputation of plink formated genotype calls using phaseit and impute2:
[here](https://github.com/HannahVMeyer/Genetic-association-studies/tree/master/imputation).

## 3. reference
Shell scripts for downloading and processing of reference data sets used in
human genetic association studies.
1. [1000Genomes](https://github.com/HannahVMeyer/Genetic-association-studies/tree/master/reference)

## 4. GWAS
Functions for processing and depicting of genome-wide association analysis
results, find them [here](https://github.com/HannahVMeyer/Genetic-association-studies/tree/master/GWAS)

## 5. utils
Generally useful functions for plotting and data processing,
find them [here](https://github.com/HannahVMeyer/Genetic-association-studies/tree/master/utils).


