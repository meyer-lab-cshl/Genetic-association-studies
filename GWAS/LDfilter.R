#' Select significant loci from association analyssis and filter for LD
#'
#' Results of genetic association study are filtered by applying the specified
#' threshold to the association p-values. Variants passing the threshold are
#' pruned for LD and only the most significant variant per locus is retained.
#' Information about tag SNPs has to be provided in external files that the user
#' needs to have read access to. Tags should be in plink tag-list format
#' (\url{https://www.cog-genomics.org/plink/1.9/ld#show_tags}), and an example
#' tag set for European samples is created here: \url{}.
#'
#' @param gwas [NrSNP x 3 (+x)] [data.frame] of association results with genetic
#' variant identifier (snp_name column), chromosome identifier (chr_name column)
#' and association p-value (p_name column) for all NrSNP. Possible additional
#' columns are the allele frequency of the genetic variant (af_name column), the
#' imputation quality score (info_name column) and the chromsomal position of
#' the variant (bp_name column).
#' @param tags_dir [character] full/path/to/directory/with/genetic/tags. Tags
#' should be in plink tag-list format
#' (\url{https://www.cog-genomics.org/plink/1.9/ld#show_tags}) as generated by
#' plink --bfile bfile --show-tags all --tag-r2 r2 --tag-kb kb --out outfile.
#' Outfile should follow the format tags_dir/tags_prefix[chrID]tags_suffix,
#' where the chromosome ID will be dynamically inserted depending on the
#' chromsome ID of significant variants.
#' @param tags_prefix [character] Character string of tags file before
#' chromosome identifier, see tags_dir.
#' @param tags_suffix [character] Character string of tags file after
#' chromosome identifier, see tags_dir.
#' @param threshold [double] threshold for considering association significant.
#' @param snp_name [character] Name of genetic variant identifier column.
#' @param chr_name [character] Name of chromosome identifier column. Chromosomes
#' can be specified as simple numbers i.e. e.g. 1 or '1' or as 'chr1'.
#' @param p_name [character] Name of p-value column. Entries have to be floats.
#' @param bp_name [character, optional] Name of base pair position column.
#' Entries have to be integers.
#' @param af_name [character, optional] Name of allele frequency column. Entries
#' have to be floats.
#' @param info_name [character, optional] Name of imputation info criterion
#' column. Entries have to be floats.
#' @param is.negLog [logical] indicates if [p_name] column is converted to
#' -log10(p-value)
#' @return named [list] with all significant loci (sig) and most significant
#' genetic variant per locus (sig_no_ld).

filterSigLoci4LD <- function(gwas, tags_dir, tags_prefix, tags_suffix,
                          threshold=5*10^(-8), chr_name="CHR",
                          p_name="P", snp_name="SNP", bp_name="BP",
                          af_name="AF", info_name="INFO", is.negLog=FALSE) {

    if (!(chr_name %in% names(gwas))) {
        stop(paste("Column", chr_name, "not found!"))
    }
    if (!(p_name %in% names(gwas))) {
        stop(paste("Column", p_name, "not found!"))
    }
    if (!(snp_name %in% names(gwas))) {
        stop(paste("Column", snp_name, "not found!"))
    }

    colnames(gwas)[colnames(gwas) == snp_name] <- "SNP"
    colnames(gwas)[colnames(gwas) == chr_name] <- "CHR"
    colnames(gwas)[colnames(gwas) == p_name] <- "P"

    if (bp_name %in% names(gwas)) {
        colnames(gwas)[colnames(gwas) == bp_name] <- "BP"
    }
    if (af_name %in% names(gwas)) {
        colnames(gwas)[colnames(gwas) == af_name] <- "AF"
    }
    if (info_name %in% names(gwas)) {
        colnames(gwas)[colnames(gwas) == info_name] <- "INFO"
    }

    if (!is.negLog) {
        sig <- gwas[gwas$P < threshold,]
    } else {
        sig <- gwas[gwas$P > -log10(threshold),]
    }

    if (nrow(sig) == 0) {
        return(list(sig=NULL, sig_wo_ld=NULL))
    } else if (nrow(sig) == 1) {
        return(list(sig=sig, sig_wo_ld=sig))
    } else {
        sigChr <- unique(sig$CHR)
        ld <- lapply(sigChr, function(chr) {
            tags_file <- paste(tags_dir, "/", tags_prefix, chr, tags_suffix,
                               sep="")
            write.table(dplyr::filter(sig, CHR==chr)$SNP, sig_file, sep="\t",
                        row.names=FALSE, col.names=FALSE, quote=FALSE)
            chr <- gsub("chr", "", chr)
            sig_file <- paste(tempdir(), "/sigChr", chr, ".txt", sep="")
            out_file <- paste(tempdir(), "/sigChr", chr, "_ld.txt", sep="")
            system(paste("bash /homes/hannah/projects/GWAS/helperLDfilter.sh ",
                         sig_file, " ", tags_file, " ",
                     out_file, sep=""), wait=TRUE)
            if (file.info(out_file)$size != 0) {
                ld <- read.table(out_file, header=FALSE, stringsAsFactors=FALSE)
                colnames(ld) <- c("SNP","CHR","BP", "NTAG", "LEFT", "RIGHT",
                              "KBSPAN", "TAGS")
                system(paste("rm", sig_file, out_file))
            } else {
                ld <- NULL
                system(paste("rm", sig_file))
            }
            return(ld)
        })
        ld <- do.call(rbind, ld)

        snp2ld <- lapply(sig$SNP, findLD,  LD=ld)
        snp2ld <- do.call(rbind, snp2ld)

        sig_wo_ld <- lapply(split(snp2ld, 1:nrow(snp2ld)),
                            filterLD, allLD=snp2ld,
                            is.negLog=is.negLog, gwas=gwas)
        sig_wo_ld <- do.call(rbind, sig_wo_ld)
        sig_wo_ld <- sig_wo_ld[!duplicated(sig_wo_ld$SNP),]
        return(list(sig=sig, sig_wo_ld=sig_wo_ld))
    }
}

#' Write results of filterSigLoci4LD
#'
#' @param results_filterSigLoci4LD named [list] returned from filterSigLoci4LD
#' containing all significant genetic variants (sig) and the most significant
#' variant per locus (sig_no_ld).
#' @param threshold [double] threshold applied for significance filtering;
#' appended to output file name; for details see [name].
#' @param directory [character] full/path/to/output/directory; needs user
#' writing permissions; for output file name see [name].
#' @param name [character] name prepanded to output file; output file name
#' generated as [directory]/[name]_sig[threshold].txt and
#' [directory]/[name]_sig[threshold]_ldFiltered.txt
#' @param bed [logical] indicates if results should additonally be saved in UCSC
#' bed format (\url{https://genome.ucsc.edu/FAQ/FAQformat.html#format1}); useful
#' for displaying results in USCS genome browser.
#' return Null, simplies writes results to [directory]/[name]_sig[threshold].txt
#' and [directory]/[name]_sig[threshold]_ldFiltered.txt
writeSig <- function(results_filterSigLoci4LD, threshold, directory,
                     name="Results", bed=FALSE) {
    sig <- results_filterSigLoci4LD$sig
    sig_wo_ld <- results_filterSigLoci4LD$sig_wo_ld
    threshold = gsub("-", "", threshold)
    if (!is.null(sig)) {
        write.table(sig, paste(directory, "/", name, "_sig", threshold,".txt",
                               sep=""),
                    row.names=FALSE, col.names=TRUE, quote=FALSE)
        if (bed) {
            sig_bed  <- data.frame(paste("chr", sig$CHR, sep=""),
                                   sig$BP -1, sig$BP, sig$SNP)
            write.table(pqtlExX_bed, paste(directory, "/", name, "_sig",
                                           threshold, ".bed", sep=""),
                    row.names=FALSE, col.names=FALSE, quote=FALSE)
        }
    }
    if (!is.null(sig_wo_ld)) {
        write.table(sig_wo_ld, paste(directory, "/", name, "_sig", threshold,
                                     "_ldFiltered.txt", sep=""),
                    row.names=FALSE, col.names=TRUE, quote=FALSE)
        if (bed) {
            sig__wo_ld_bed  <- data.frame(paste("chr", sig_wo_ld$CHR, sep=""),
                                   sig_wo_ld$BP -1, sig_wo_ld$BP, sig_wo_ld$SNP)
            write.table(pqtlExX_bed, paste(directory, "/", name, "_sig",
                                           threshold, "_ldFiltered.bed", sep=""),
                    row.names=FALSE, col.names=FALSE, quote=FALSE)
        }
    }
}

#' Find LD buddies of genetic variants
#'
#' Called from filterSigLoci4LD.
#'
#' From a list containing genetic variants and their respective tag variants
#' (LD), look up provided genetic variant ID (SNP) and return vector with
#' variant ID (SNP), chromosome ID (CHR) and tag variants (TAGS). If no tags
#' exists, CHR and TAGS are set to NONE.
#' @param SNP [character] genetic variant identifier.
#' @param LD [list] with genetic variant identifier and their tag variants in
#' plink tag-list format
#' (\url{https://www.cog-genomics.org/plink/1.9/ld#show_tags}) as generated by
#' plink --bfile bfile --show-tags all --tag-r2 r2 --tag-kb kb --out outfile.
#' @return data.frame with SNP (genetic variant ID), CHR (chromosome ID of
#' variant) and TAGS (pipe separated string of genetic variants tagged by SNP).

findLD <- function(SNP, LD) {
    if (SNP %in% LD$SNP && LD$TAGS[LD$SNP == SNP] != "NONE") {
        tmp <-  dplyr::select(LD[LD$SNP == SNP,], SNP, CHR, TAGS)
    } else {
        tmp <- c(SNP, "NONE","NONE")
    }
    names(tmp) <- c("SNP", "CHR", "TAGS")
    return(tmp)
}

#' Filter GWAS results for one variant per loci
#'
#' Called from filterSigLoci4LD.
#'
#' Takes a genetic variant ID, the output of findLD (data.frame with SNP
#' (genetic variant ID), CHR (chromosome ID of variant) and TAGS (pipe separated
#' string of genetic variants tagged by SNP) and the overall association results
#' to select the genetic variant per locus with the minimum association p-value.
#' @param snp [character] genetic variant ID.
#' @param allLD output from findLD ie [NrSigSNP x 3] [data.frame] with SNP
#' (genetic variant ID), CHR (chromosome ID of variant) and TAGS (pipe separated
#' string of genetic variants tagged by SNP) for all significant SNPs
#' (NRSigSNPs).
#' @param is.negLog [logical] indicates if [p_name] column is converted to
#' -log10(p-value)
#' @return [data.frame] with one row from input gwas, containing the genetic
#' variant with the minimal p-value for that locus (and all original columns of
#' the gwas input data.frame).

filterLD <- function(snp, allLD, gwas, is.negLog=FALSE) {
    tags <- c(strsplit(allLD$TAGS[allLD$SNP == snp], split="|", fixed=TRUE)[[1]],
              snp$SNP)
    if ("NONE" %in% tags) {
        minP_ld_snp <- gwas[gwas$SNP == snp$SNP,]
    } else {
        ld_snps <- gwas[gwas$SNP %in% tags,]
        if (nrow(ld_snps) > 1) {
            if (is.negLog) {
                minP_ld_snp_tmp <- ld_snps[which.max(ld_snps$P),]
            } else {
                minP_ld_snp_tmp <- ld_snps[which.min(ld_snps$P),]
            }
            tags_minP_ld_snp <-
                c(strsplit(allLD$TAGS[allLD$SNP == minP_ld_snp_tmp$SNP],
                           split="|", fixed=TRUE)[[1]], minP_ld_snp_tmp$SNP)
            ld_snps_minP <- gwas[gwas$SNP %in% tags_minP_ld_snp,]
            if (is.negLog) {
                minP_ld_snp <- ld_snps_minP[which.max(ld_snps_minP$P),]
            } else {
                minP_ld_snp <- ld_snps_minP[which.min(ld_snps_minP$P),]
            }
        } else {
            minP_ld_snp <- ld_snps
        }
    }
    return(minP_ld_snp)
}


