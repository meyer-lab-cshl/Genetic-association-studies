#' Filter bgenie output files for biallelic variants
#'
#' @param o_bgenie [data.frame] of bgenie results, with bgenie results header as
#' column names of the data.frame
#' @return data.frame of bgenie results containing biallelic genetic variants
#' only.
biallelic <- function(o_bgenie){
    multi_rsID <- o_bgenie$rsid[duplicated(o_bgenie$rsid)]
    tmp <- o_bgenie[!(o_bgenie$rsid %in% multi_rsID),]
    multi_pos <- tmp$pos[duplicated(tmp$pos)]
    tmp <- tmp[!(tmp$rsid %in% multi_pos),]
    return(tmp)
}

#' Read and process bgenie output files
#'
#' Read the .gz bgenie output files of the specified chromosome and optionally
#' apply filters on minor allele frequency and imputation info score. In
#' addition, one can choose to select for biallelic variants only.
#'
#' @param chr [integer] chromosome number of output file.
#' @param directory [character] full path to directory where bgenie results are
#' saved.
#' @param name [character] name of bgenie file; the full filename is
#' reconstructed as {directory}/{name}{chr}.gz. Alternatively, prefix can
#' directly be provided.
#' @param prefix [character] specifies the full path the bgenie file before the
#' chromosome number, ie the full file name is {prefix}{chr}.gz.
#' @param maf [double] optional minor allele frequency filtering threshold. If
#' not desired set to NULL.
#' @param info [double] optional info score filtering threshold. If not desired
#' set to NULL.
#' @param biallelicOnly [logical] if TRUE genetic variants are filtered for
#' biallelic variants only.
#' @return data.frame of (filtered) bgenie results.

readBgenieOutput <- function(chr, directory, name, maf=0.01, info=0.4,
                             biallelicOnly=TRUE) {
    chrFile <- paste(directory, "/", name, chr,".gz",
                     sep="")
    if (!file.exists(chrFile)) {
        message("Bgenie association output for chr", chr, " cannot be found, ",
        "skip to next chromosome")
        return(NULL)
    }
    readString <- paste("zcat", chrFile)
    if (verbose) message("Reading association results from chr", chr)
    tmp <- data.table::fread(readString, sep=" ", stringsAsFactors=FALSE,
                             data.table=FALSE, header=TRUE)
    # Filter based on MAF
    if (!is.null(maf)) {
        tmp <- tmp[tmp$af > maf,]
    }
    tmp <- tmp[!tmp$af %in% c(0,1),]

    # Filter on info criterion
    if (!is.null(info)) {
        tmp <- tmp[tmp$info > info,]
    }

    # Filter for biallelic
    if (biallelicOnly) {
        tmp <- biallelic(tmp)
    }
    if (all(is.na(tmp$chr))) tmp$chr <- chr
    return(tmp)
}

#' Format bgenie output files for ld score regression
#'
#' @param o_bgenie [data.frame] of bgenie results, with bgenie results header as
#' column names of the data.frame
#' @param sumstat [character] name of column used for summary statistic. If
#' effect size estimates (beta) are used, output column will be beta, else
#' SUMSTAT.
#' @param ldshub_snps [data.frame] with at least 'SNP' column containing the
#' SNP IDs relevant for analysis on LDhub
#' \url{http://ldsc.broadinstitute.org/ldhub/}.
#' @return [data.frame] of bgenie results formated for munge_sumstats.py and
#' LDhub. See Details for munge_sumstats.py format.
#' @details  munge_sumstats.py from ldsc regression takes a file with the
#' following column names: \itemize{
#' \item SNP: Variant ID
#' \item P: p-Value,
#' \item A1: Allele 1, interpreted as ref allele for signed sumstat.
#' \item A2: Allele 2, interpreted as non-ref allele for signed sumstat
#' \item N: Sample size
#' \item BETA: [linear/logistic] regression coefficient (0 --> no effect;
#' above 0 --> A1 is trait/risk increasing)
#' \item INFO: INFO score (imputation quality; higher --> better imputation)
#' \item FRQ: Allele frequency
#' \item SIGNED_SUMSTAT: Directional summary statistic as specified by
#' --signed-sumstats.
#' }
bgenie2ldsc <- function(o_bgenie, sumstat, ldshub_snps){
        data.table::setnames(o_bgenie, new="SNP", old="rsid")
        data.table::setnames(o_bgenie, new="A1", old="a_0")
        data.table::setnames(o_bgenie, new="A2", old="a_1")
        data.table::setnames(o_bgenie, new="INFO", old="info")
        data.table::setnames(o_bgenie, new="FRQ", old="af")
        if (grepl("beta", tolower(sumstat))) {
            data.table::setnames(o_bgenie, new="BETA",
                                 old=colnames(o_bgenie)[grepl(sumstat,
                                                          colnames(o_bgenie))])
        } else {
            data.table::setnames(o_bgenie, new="SUMSTAT",
                                 old=colnames(o_bgenie)[grepl(sumstat,
                                                          colnames(o_bgenie))])
        }
        o_bgenie$P <- 10^(-o_bgenie[, grepl("-log10p", colnames(o_bgenie))])
        o_bgenie$N <- nrow(o_bgenie)
        o_bgenie <- o_bgenie[,!grepl("[_log10psechrt]{2,6}", colnames(o_bgenie))]
        if (!is.null(ldshub_snps)) {
            o_bgenie <- o_bgenie[o_bgenie$SNP %in% ldshub_snps$SNP, ]
        }
        return(o_bgenie)
}
