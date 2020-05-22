#' Meta-analyses of summary statistics from uni-variate GWAS
#'
#' The summary statistics of uni-variate GWAS are combined to
#' 'pseudo-multitrait' p-values as described in Boloorma et al 2013
#'
#' @param summarystats [NrSNP x P ] matrix of summary stastics from GWAS of
#' NrSNP number of SNPs for P traits.
#' @return named list with i) [NrSNPs x 1] vector of meta-analysis p-values
#' (multitrait_p) and ii) [NrSNPs x 1] vector of meta-analysis of chi square
#' statistic (chi2stats).
#' @references Bolormaa S, Pryce JE, Reverter A, Zhang Y, Barendse W, et al.
#' (2014) A Multi-Trait, Meta-analysis for Detecting Pleiotropic Polymorphisms
#' for Stature, Fatness and Reproduction in Beef Cattle. PLoS Genet 10(3):
#' e1004198. doi:10.1371/journal.pgen.1004198
pseudoMultitrait <- function(summarystats) {
    summarystats <- as.matrix(summarystats)
    df <- dim(summarystats)[2]
    cor_sstats<- cor(summarystats)
    cor_sstats_inv <- chol2inv(chol(cor_sstats))
    chi2stats <- apply(summarystats, 1, function(z) {
                         t(z) %*% cor_sstats_inv %*% z
            })
    multitrait_p <- as.vector(pchisq(chi2stats, df, lower.tail=FALSE))
    return(list(multitrait_p=multitrait_p, chi2stats=chi2stats))
}

#' Estimate effective number of tests
#'
#' Use eigenvalues of phenotype correlation matrix to estimate the effective
#' number of tests conducted. Result can be taken to adjust for Bonferoni-type
#' multiple hypothesis, with p/m_effective rather than p/m_total.
#' @param mat [N x P] matrix of phenotypes [double] from N samples and P tested
#' traits.
#' return Effective number of tests [double]
#' @references Galwey (2009) A new measure of the effective number of tests,
#' a practical tool for comparing families of non-independent significance
#' tests. Genetic Epidemiology.
Teff <- function(mat) {
    # 1. get correlation matrix
    corr_matrix <-  cor(mat,  method="spearman")
    # 2. Get eigenvalues of correlation matrix:
    eigenval <- eigen(corr_matrix, only.value=TRUE, symmetric=TRUE)$values
    # 3. Determine effective number of tests:
    t <- sum(sqrt(eigenval))^2/sum(eigenval)
    return(t)
}
