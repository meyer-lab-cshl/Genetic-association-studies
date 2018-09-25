# GWAS

The following directory contains plotting and analysis functions for depicting
and processing results from genome-wide association studies.

1. bgenieResults.r:
    * Function for reading chromosome-wide bgenie results, including optional
    posthoc minor allele frequence, imputation info score and biallelic variants
    filtering options.
1. metaanalysis.r:
    * Function for estimating the effective number of hypothesis tests conducted,
    based on the eigenvalues of the correlation matrix of the phenotypes (based
    on [Galwey (2009) Genetic Epidemiology](https://onlinelibrary.wiley.com/doi/abs/10.1002/gepi.20408)).
    * Function for combining univariate association results into a meta-analysis
    p-values (based on [Boloorma (2014) Plos Genetics](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004198)).
1. plots.r:
    * Function for manhattan of GWAS results
    * Function for qqplot of GWAS results
1. prepGarfield.r:
    * Function formating and writing association results for analyses with
    [GARFILED](https://www.ebi.ac.uk/birney-srv/GARFIELD/).
1. relatedness.r:
    * Function filtering related samples from list of samples in cohort with the aim of retaining maximum number of samples in the cohort.
