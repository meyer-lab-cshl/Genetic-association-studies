#' Meta-analyses of summary statistics from uni-variate GWAS
#'
#' The summary statistics of uni-variate GWAS are combined to
#' 'pseudo-multitrait' p-values as described in Boloorma et al 2013
#'
#' @param summarystats [NrSNP x P ] matrix of summary stastics from GWAS of
#' NrSNP number of SNPs for P traits.
#' @return [NrSNPs x 1] vector of meta-analysis p-values
#' @references Bolormaa S, Pryce JE, Reverter A, Zhang Y, Barendse W, et al.
#' (2014) A Multi-Trait, Meta-analysis for Detecting Pleiotropic Polymorphisms
#' for Stature, Fatness and Reproduction in Beef Cattle. PLoS Genet 10(3):
#' e1004198. doi:10.1371/journal.pgen.1004198
pseudoMultitrait <- function(summarystats) {
    df <- dim(summarystats)[2]
    cor_sstats<- cor(summarystats)
    cor_sstats_inv <- chol2inv(chol(cor_sstats))
    chi2stats <- apply(summarystats, 1, function(z) {
                         t(z) %*% cor_sstats_inv %*% z
            })
    multitrait_p <- as.vector(pchisq(chi2stats, df, lower.tail=FALSE))
    return(multitrait_p)
}
