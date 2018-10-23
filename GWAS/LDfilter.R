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
#' @param tags_dir [character] full/path/to/directory/with/genetic/tags. Should
#' contain genetic variants and their tags. Variants and their tags should be
#' saved in per-chromosome files, such that all variants for [chr] can be found
#' in the file tags_dir/tags_prefix[chrID]tags_suffix. Each
#' tags_dir/tags_prefix[chrID]tags_suffix file should contain at least
#' two columns, one with the genetic variant identifier [tags_columnID] and one
#' with a string of its tags [tags_columnTags], which are separated by
#' [tagsSplit]. Files can be generated in plink tag-list format via
#' (\url{https://www.cog-genomics.org/plink/1.9/ld#show_tags})
#' plink --bfile bfile --show-tags all --tag-r2 r2 --tag-kb kb --out outfile.
#' Alternatively, precomputed LD-list for European Samples for strict and
#' lenient LD pruning can be found as part of the supplementary data of the
#' GARFIELD package which can be found here:
#' https://www.ebi.ac.uk/birney-srv/GARFIELD/
#' Column names for [tags_columnID] and [tags_columnTags] will have to be added.
#' @param tags_ID [character] Column name in
#' tags_dir/tags_prefix[chr]tags_suffix containing the genetic variant
#' identifier.
#' @param tags_Tags [character] Column name in
#' tags_dir/tags_prefix[chr]tags_suffix containing the list of tags separated by
#' [tagsSplit].
#' @param tagsSplit [character] Separator of tag IDs in
#' tags_dir/tags_prefix[chr]tags_suffix.
#' @param tags_prefix [character] Character string of tags file before
#' chromosome identifier, see tags_dir.
#' @param tags_suffix [character] Character string of tags file after
#' chromosome identifier, see tags_dir.
#' @param matchBy [character]  Choice depends on format of pre-computed LD file:
#' if tags_dir/prefix[chr]suffix contains SNP TAG information based on
#' positions, set to 'BP', if SNP TAG information is based on SNP names (ie
#' rsIDs), set to 'SNP'.
#' @param threshold [double] Threshold for considering association significant.
#' @param gwas_snp [character] Name of genetic variant identifier column in
#' gwas.
#' @param gwas_chr [character] Name of chromosome identifier column in gwas.
#' Chromosomes can be specified as simple numbers i.e. e.g. 1 or '1' or as
#' 'chr1'.
#' @param gwas_p [character] Name of p-value column in gwas. Entries have to be
#' floats.
#' @param gwas_bp [character, optional] Name of base pair position column in
#' gwas. Entries have to be integers.
#' @param gwas_af [character, optional] Name of allele frequency column in gwas.
#' Entries have to be floats.
#' @param gwas_info [character, optional] Name of imputation info criterion
#' column in gwas. Entries have to be floats.
#' @param is.negLog [logical] Indicates if [p_name] column is converted to
#' -log10(p-value)
#' @param greedy [logical] Greedy pruning i.e. multiple passes of LDfilter till
#' every variant in LD is removed.
#' @return Named [list] with all significant loci (sig) and most significant
#' genetic variant per locus (sig_no_ld).

filterSigLoci4LD <- function(gwas, tags_dir, tags_prefix, tags_suffix,
                             tags_ID="ID", tags_Tags="TAGS", tagsSplit=",",
                             threshold=5*10^(-8), gwas_chr="CHR",
                             gwas_p="P", gwas_snp="SNP", gwas_bp="BP",
                             gwas_af="AF", gwas_info="INFO", is.negLog=FALSE,
                             matchBy=c("BP", "SNP"), greedy=TRUE) {

    matchBy <- match.arg(matchBy)

    if (!(gwas_chr %in% names(gwas))) {
        stop(paste("Column", gwas_chr, "not found!"))
    }
    if (!(gwas_p %in% names(gwas))) {
        stop(paste("Column", gwas_p, "not found!"))
    }
    if (!(gwas_snp %in% names(gwas))) {
        stop(paste("Column", gwas_snp, "not found!"))
    }

    colnames(gwas)[colnames(gwas) == gwas_snp] <- "SNP"
    colnames(gwas)[colnames(gwas) == gwas_chr] <- "CHR"
    colnames(gwas)[colnames(gwas) == gwas_p] <- "P"

    if (gwas_bp %in% names(gwas)) {
        colnames(gwas)[colnames(gwas) == gwas_bp] <- "BP"
    }
    if (gwas_af %in% names(gwas)) {
        colnames(gwas)[colnames(gwas) == gwas_af] <- "AF"
    }
    if (gwas_info %in% names(gwas)) {
        colnames(gwas)[colnames(gwas) == gwas_info] <- "INFO"
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
        message("Significant SNPs on chromosomes ", paste(sigChr, collapse=","))
        sig_wo_ld <- lapply(sigChr, function(chr) {
            message("Getting tags for chr", chr, "...")
            chr_str <- gsub("chr", "", chr)
            sig_on_chr <- dplyr::filter_(sig, ~CHR == chr)
            tags_file <- paste(tags_dir, "/", tags_prefix, chr_str, tags_suffix,
                               sep="")
            sig_file <- paste(tempdir(), "/sigChr", chr_str, ".txt", sep="")
            out_file <- paste(tempdir(), "/sigChr", chr_str, "_ld.txt", sep="")
            if (matchBy == "SNP") {
                write.table(sig_on_chr$SNP, sig_file, sep="\t",
                            row.names=FALSE, col.names=FALSE, quote=FALSE)
            } else if (matchBy == "BP") {
                write.table(sig_on_chr$BP, sig_file, sep="\t",
                            row.names=FALSE, col.names=FALSE, quote=FALSE)
            }
            system(paste("bash /homes/hannah/projects/GWAS/helperLDfilter.sh ",
                         sig_file, " ", tags_file, " ",  out_file,
                         sep=""), wait=TRUE)
            if (file.info(out_file)$size != 0) {
                ld <- read.table(out_file, header=TRUE, stringsAsFactors=FALSE)
                if (!(tags_Tags %in% names(ld))) {
                    stop(paste("Column", tags_Tags, "not found!"))
                }
                if (!(tags_ID %in% names(ld))) {
                    stop(paste("Column", tags_ID, "not found!"))
                }
                colnames(ld)[colnames(ld) == tags_ID] <- "ID"
                colnames(ld)[colnames(ld) == tags_Tags] <- "TAGS"
                ld$CHR <- chr_str
                system(paste("rm", sig_file, out_file))
            } else {
                ld <- NULL
                system(paste("rm", sig_file))
            }
            message("Match variants to their ld tags...")
            snp2ld <- lapply(1:nrow(sig_on_chr), function(id) {
                                 findLD(SNP=sig_on_chr$SNP[id],
                                        BP=sig_on_chr$BP[id],
                                        CHR=sig_on_chr$CHR[id],
                                        LD=ld, matchBy=matchBy,
                                        id_name='ID')
                            })
            snp2ld <- do.call(rbind, snp2ld)

            message("Filter loci for variant with minimal p-value...")
            sigSNP_wo_ld <- sapply(snp2ld$SNP, filterLD, allLD=snp2ld,
                                is.negLog=is.negLog, gwas=gwas,
                                tagsSplit=tagsSplit, matchBy=matchBy)
            sigSNP_wo_ld <- sigSNP_wo_ld[!duplicated(sigSNP_wo_ld)]

            old_count <- length(sigSNP_wo_ld)
            new_count <- length(sigSNP_wo_ld) + 1

            if (greedy) {
                while(old_count != new_count) {
                    old_count <- length(sigSNP_wo_ld)
                    newSNP <- sapply(sigSNP_wo_ld, filterLD, allLD=snp2ld,
                                  is.negLog=is.negLog, gwas=gwas,
                                  tagsSplit=tagsSplit)
                    newSNP <- newSNP[!duplicated(newSNP)]
                    new_count <- length(newSNP)
                    sigSNP_wo_ld <- newSNP
                }
            }
            sigAll_wo_ld <- gwas[gwas$SNP %in% sigSNP_wo_ld,]
            return(sigAll_wo_ld)
        })
    sig_wo_ld <- do.call(rbind, sig_wo_ld)
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
#' (LD), look up provided genetic variant ID by either variant name (SNP) or
#' unique chromosomal position (via CHR and BP) and return vector with
#' variant ID (SNP), chromosome ID (CHR), chromosomal position (BP) and tag
#' variants (TAGS). If no tags exist, TAGS are set to NONE.
#' @param SNP [character] genetic variant identifier.
#' @param CHR [character] chromosome of genetic variant; either as e.g. chr1 or
#' 1.
#' @param BP [character] chromosomal position of genetic variant.
#' @param LD [list] with genetic variant identifier and their tag variants in
#' named (see options *_name) plink tag-list format
#' (\url{https://www.cog-genomics.org/plink/1.9/ld#show_tags}) as generated by
#' plink --bfile bfile --show-tags all --tag-r2 r2 --tag-kb kb --out outfile and
#' @param id_name [character] Name of variant identifier column in LD (either
#' column with genetic variant identifier if matchBy=='SNP', or chromosomal
#' position if matchBy=='BP').
#' @param chr_name [character] Name of chromosome identifier column in LD.
#' @param tags_name [character] Name of tags column in LD.
#' @param matchBy [character]  Choice depends on format of pre-computed LD file:
#' if tags_dir/prefix[chr]suffix contains SNP TAG information based on
#' positions, set to 'BP', if SNP TAG information is based on SNP names (ie
#' rsIDs), set to 'SNP'.
#' @return data.frame with SNP (genetic variant ID), CHR (chromosome ID of
#' variant) and TAGS (pipe separated string of genetic variants tagged by SNP).

findLD <- function(SNP, CHR, BP, LD, id_name, chr_name="CHR", tags_name="TAGS",
                   matchBy=c('BP', 'SNP')) {
    matchBy <- match.arg(matchBy)
    if (!(chr_name %in% names(LD))) {
        stop(paste("Column", chr_name, "not found in LD!"))
    }
    if (!(tags_name %in% names(LD))) {
        stop(paste("Column", tags_name, "not found in LD!"))
    }
    if (!(id_name %in% names(LD))) {
        stop(paste("Column", id_name, "not found!"))
    }
    colnames(LD)[colnames(LD) == id_name] <- "ID"
    colnames(LD)[colnames(LD) == chr_name] <- "CHR"
    colnames(LD)[colnames(LD) == tags_name] <- "TAGS"
    tmp <- NULL
    if (matchBy == 'SNP') {
        if (SNP %in% LD$ID && LD$TAGS[LD$ID == SNP] != "NONE") {
            tmp <- LD[LD$ID == SNP,]
        }
    } else if (matchBy == 'BP') {
        if (BP %in% LD$ID && LD$TAGS[LD$ID == BP] != "NONE") {
            tmp <- LD[LD$ID == BP & LD$CHR == CHR,]
        }
    } else {
        stop("matchBy has to be SNP or BP, but ", matchBy, "provided.")
    }
    if (is.null(tmp)) {
        tmp <- data.frame(SNP=SNP, CHR=CHR, BP=BP, TAGS="NONE",
                          stringsAsFactors=FALSE)
    } else {
        tmp$SNP <- SNP
        tmp$BP <- BP
        tmp <- dplyr::select_(tmp, ~SNP, ~CHR, ~BP, ~TAGS)
        names(tmp) <- c("SNP", "CHR", "BP", "TAGS")
    }
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
#' @param snp_name_ld [character] Name of genetic variant identifier column in
#' LD.
#' @param tags_name_ld [character] Name of tags column in LD.
#' @param snp_name_gwas [character] Name of genetic variant identifier column in
#' gwas.
#' @param p_name_gwas [character] Name of p value column in gwas.
#' @param tagsSplit [character] separator used for tags list.
#' @return variant ID with minimal p-value of locus.

filterLD <- function(snp, allLD, gwas, is.negLog=FALSE, tags_name_ld="TAGS",
                     snp_name_ld="SNP", chr_name_ld="CHR", p_name_gwas="P",
                     snp_name_gwas="SNP", tagsSplit=",",
                     matchBy=c("BP", "SNP")) {
    matchBy <- match.arg(matchBy)

    if (!(tags_name_ld %in% names(allLD))) {
        stop(paste("Column", tags_name_ld, "not found in allLD!"))
    }
    if (!(snp_name_ld %in% names(allLD))) {
        stop(paste("Column", snp_name_ld, "not found in allLD!"))
    }
    if (!(chr_name_ld %in% names(allLD))) {
        stop(paste("Column", chr_name_ld, "not found in allLD!"))
    }
    if (!(p_name_gwas %in% names(gwas))) {
        stop(paste("Column", p_name_gwas, "not found in gwas!"))
    }
    if (!(snp_name_gwas %in% names(gwas))) {
        stop(paste("Column", snp_name_gwas, "not found in gwas!"))
    }
    colnames(allLD)[colnames(allLD) == snp_name_ld] <- "SNP"
    colnames(allLD)[colnames(allLD) == tags_name_ld] <- "TAGS"
    colnames(allLD)[colnames(allLD) == chr_name_ld] <- "CHR"
    colnames(gwas)[colnames(gwas) == snp_name_gwas] <- "SNP"
    colnames(gwas)[colnames(gwas) == p_name_gwas] <- "P"

    if (matchBy == "SNP") {
        tags <- c(strsplit(allLD$TAGS[allLD$SNP == snp], split=tagsSplit,
                       fixed=TRUE)[[1]], snp)
    } else if (matchBy == "BP") {
        tags <- c(strsplit(allLD$TAGS[allLD$SNP == snp], split=tagsSplit,
                       fixed=TRUE)[[1]], gwas$BP[gwas$SNP == snp])
    } else {
        stop("matchBy has to be SNP or BP, but ", matchBy, "provided.")
    }

    chr <- as.numeric(allLD$CHR[allLD$SNP == snp])

    if ("NONE" %in% tags) {
        minP_ld_snp <- snp
    } else {
        if (matchBy == "SNP") {
            ld_snps <- dplyr::select_(gwas[gwas$SNP %in% tags,], ~SNP, ~P)
        } else {
            ld_snps <- dplyr::select_(gwas[gwas$CHR == chr & gwas$BP %in% tags,],
                                      ~SNP, ~P)
        }
        if (nrow(ld_snps) > 1) {
            if (is.negLog) {
                minP_ld_snp <- ld_snps$SNP[which.max(ld_snps$P)]
            } else {
                minP_ld_snp <- ld_snps$SNP[which.min(ld_snps$P)]
            }
        } else {
            minP_ld_snp <- ld_snps$SNP
        }
    }
    return(minP_ld_snp)
}



