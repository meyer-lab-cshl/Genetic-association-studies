perIndividualQC <- function(qcdir, alg, maleTh, femaleTh, imissTh, hetTh,
                            highIBDTh, externalSex=NULL, family=NULL,
                            plot=FALSE, verbose=TRUE) {
    if (verbose) {
        message("Identification of individuals with discordant sex ",
                "information")
    }
    fail_sex <- check_sex(qcdir=qcdir, alg=alg, maleTh=maleTh,
                             femaleTh=femaleTh, externalSex=externalSex,
                             plot=plot)
    if (!is.null(fail_sex$fail_sex)) {
        write.table(fail_sex$fail_sex[,1:2],
                    file=paste(qcdir,"/",alg, ".fail-sexcheck.IDs", sep=""),
                    quote=FALSE, row.names=FALSE, col.names=FALSE)
    }
    if (!is.null(fail_sex$mixup)) {
        write.table(fail_sex$mixup[,1:2],
                    file=paste(qcdir,"/",alg, ".sexcheck_mixup.IDs", sep=""),
                    quote=FALSE, row.names=FALSE, col.names=FALSE)
    }
    if (verbose) {
        message("Identification of individuals with outlying missing ",
                "genotype or heterozygosity rates")
    }
    fail_het_imiss <-
        check_heterozygosity_and_missingness(qcdir=qcdir, alg=alg,
                                             imissTh=imissTh, hetTh=hetTh,
                                             plot=plot)
    if (!is.null(fail_het_imiss$fail_imiss)) {
        write.table(fail_het_imiss$fail_imiss[,1:2],
                    file=paste(qcdir, "/", alg, ".fail-imiss.IDs", sep=""),
                    quote=FALSE, row.names=FALSE, col.names=FALSE)
    }
    if (!is.null(fail_het_imiss$fail_het)) {
        write.table(fail_het_imiss$fail_het[,1:2],
                    file=paste(qcdir, "/", alg, ".fail-het.IDs", sep=""),
                    quote=FALSE, row.names=FALSE, col.names=FALSE)
    }
    if (verbose) message("Identification of related individuals")
    fail_relatedness <- check_relatedness(qcdir=qcdir, alg=alg, famfile=family,
                                          highIBDTh=highIBDTh, plot=plot)
    if (!is.null(fail_relatedness$fail_highIBD)) {
        write.table(ID[which(ID[,1] %in% fail_highIBD[,1]),],
                    file=paste(qcdir,"/", alg,".fail_IBD.txt", sep=""),
                    row.names=FALSE, quote=FALSE, col.names=TRUE, sep="\t")
    }
    if (verbose) message("Identification of individuals of divergent ancestry")
    fail_ancestry <- check_ancestry(qcdir=qcdir, alg=alg, plot=plot)
    if (!is.null(fail_ancestry)) {
        write.table(fail_ancestry,
                    file=paste(qcdir, "/",alg,".fail-ancestry.IDs",sep=""),
                    quote=FALSE, row.names=FALSE, col.names=FALSE)
    }

    system(paste("cat ",qcdir ,"/",  alg, ".fail-*.IDs | sort | uniq >",
                 qcdir, "/", alg, ".fail.IDs", sep=""), wait=TRUE)

    fail_list <- list(missing_genotype=as.vector(fail_het_imiss$fail_imiss$IID),
                      highIBD=as.vector(fail_relatedness$IID),
                      outlying_heterozygosity=
                          as.vector(fail_het_imiss$fail_het$IID),
                      mismatched_sex=as.vector(fail_sex$fail_sex$IID),
                      diff_ancestry=as.vector(fail_ancestry$IID))
    return(fail_list)
}

#' Overview of per individual QC
#'
#' afdsgsd
overviewPerIndividualQC <- function(fail_list, plot) {
    list2counts <- function(element, all_names) {
        all_names[!(all_names %in% element)] <- 0
        all_names[all_names %in% element] <- 1
        return(as.numeric(all_names))
    }
    fail_list <- fail_list[!sapply(fail_list, is.null)]
    unique_samples_fail_all <- unique(unlist(fail_list))

    # a) overview QC fails independent of ethnicity
    fail_list_wo_ancestry <- fail_list[!names(fail_list) == "diff_ancestry"]
    unique_samples_fail_wo_ancestry <- unique(unlist(fail_list_wo_ancestry))

    fail_counts_wo_ancestry <- sapply(fail_list_wo_ancestry, list2counts,
                                      unique_samples_fail_wo_ancestry)
    rownames(fail_counts_wo_ancestry) <- unique_samples_fail_wo_ancestry


    # b) overview of QC and ancestry fails
    fail_counts_all <- sapply(list(QC_fail=unique_samples_fail_wo_ancestry,
                                   Ancestry_fail= fail_list$diff_ancestry),
                              list2counts, unique_samples_fail_all)
    rownames(fail_counts_all) <- unique_samples_fail_all

    if (plot) {
        limma::vennDiagram(fail_counts_wo_ancestry,
                           names=names(fail_list_wo_ancestry),
                           main="QC Fail overview (all ethnicities)")
        limma::vennDiagram(fail_counts_all,
                           names= names(fail_counts_all),
                           main=paste("Overlap between 'per individual QC'",
                                      "and ancestry analysis"))
    }
    nr_fail_samples = length(unique_samples_fail_all)
    return(nr_fail_samples=nr_fail_samples, fail_count_all=fail_counts_all,
           fail_counts_wo_ancestry=fail_counts_wo_ancestry)
}

perMarkerQC <- function(qcdir, alg, mafTh=0.01, macTh=20, hweTh=1e-5,
                        lmissTh=0.01, plot=FALSE, verbose=TRUE) {
    if (verbose) message("Identification of SNPs with high missingness rate")
    fail_SNP_missingness <- check_SNP_missingness(qcdir=qcdir, alg=alg,
                                                  lmissTh=lmissTh, plot=plot)
    if (verbose) message("Identification of SNPs with deviation from HWE")
    fail_HWE <- check_HWE(qcdir=qcdir, alg=alg, hweTh=hweTh, plot=plot)
    if (verbose) message("Remove markers with a low minor allele frequency")
    fail_maf <- check_maf(qcdir=qcdir, alg=alg, mafTh=mafTh, macTh=macTh,
                          plot=plot)
    fail_list <- list(SNP_missingness=fail_SNP_missingness$SNP,
                      hwe=fail_hwe$SNP, maf=fail_maf$SNP)
    return(fail_list)
}

#' Identification of individuals with discordant sex information
#'
#' Evaluates and depicts results from plink --check-sex. Reads resulting datafile
#' qcdir/alg.sexcheck and returns IIDs for samples whose SNPSEX != PEDSEX (where
#' the SNPSEX is determined by the heterozygosity rate across X-chromosomal
#' variants).
#' Mismatching SNPSEX and PEDSEX IDs can indicate plating errors, sample-mixup
#' or generally samples with poor genotyping. In the latter case, these IDs are
#' likely to fail other QC steps as well.
#' Optionally, an external file with sample IDs [ID] and sex [Sex] can be
#' provided to double check if external and PEDSEX data (often processed at
#' different centers) match. If a mismatch between PEDSEX and SNPSEX was
#' detected, by SNPSEX == Sex, PEDSEX of these individuals can optionally be
#' updated (fixMixup=TRUE).
#' @param qcdir [character] path/to/directory/with/QC/results containing
#' alg.sexcheck results as returned by plink --check-sex. Users needs writing
#' permission to qcdir.
#' @param alg [character] prefix of QC-ed plink files, ie alg.bed, alg.bim,
#' alg.fam and alg.sexcheck.
#' @param maleTh [double] Threshold of X-chromosomal heterozygosity rate for
#' males.
#' @param femaleTh [double] Threshold of X-chromosomal heterozygosity rate for
#' females.
#' @param plot [logical] Should results be depicted as plots? When choosing this
#' option, make sure you have X-forwarding/graphical interface available for
#' interactive plotting. Alternatively, enclose calls to this function in
#' pdf(outfile) check_sex(arguments) dev.off().
#' @param fixMixup [logical] Should PEDSEX of individuals with mismatch between
#' PEDSEX and Sex, with Sex==SNPSEX automatically corrected: this will directly
#' change the alg.bim/.bed/.fam files!
#' @param externalFemale [integer/character] Identifier for 'female' in
#' externalSex.
#' @param externalMale [integer/character] Identifier for 'male' in externalSex.
#' @return named list with i) fail_sex: dataframe with FID, IID, PEDSEX, SNPSEX
#' and Sex (if externalSex was provided) of individuals failing sex check and
#' 2) mixup: dataframe with FID, IID, PEDSEX, SNPSEX and Sex (if externalSex
#' was provided) of individuals whose PEDSEX != Sex and Sex == SNPSEX.
check_sex <- function(qcdir, alg, externalSex=NULL, maleTh=0.8, femaleTh=0.2,
                      plot=FALSE, fixMixup=FALSE, externalFemale="F",
                      externalMale="M") {
    if (!file.exists(qcdir,"/", alg, ".sexcheck",sep="")){
        stop("plink --check-sex results file: ", qcdir,"/", alg,
             ".sexcheck does not exist.")
    }
    sexcheck <- read.table(paste(qcdir,"/", alg, ".sexcheck",sep=""),
                           header=TRUE, stringsAsFactors=FALSE)

    names_sexcheck <- c("FID", "IID", "PEDSEX", "SNPSEX", "STATUS", "F")
    if (!all(names_sexcheck) == names(sexcheck)) {
        stop("Header of", qcdir,"/", alg, ".sexcheck is not correct. Was your
             file generated with plink --check-sex?")
    }

    if (is.null(externalSex)) {
        fail_sex <- sexcheck[sexcheck$STATUS == "PROBLEM",]
        if (nrow(fail_sex) == 0) fail_sex <- NULL
        return(list(fail_sex=fail_sex, mixup=NULL))
    }

    if (!("Sex" %in% names(externalSex))) {
        stop("Column Sex not found in externalSex!"))
    }
    if (!("ID" %in% names(externalSex))) {
        stop("Column ID not found in externalSex!"))
    }
    sexcheck_fuse <- merge(sexcheck, externalSex, by.x="IID", by.y="ID")
    sex_mismatch <-
        apply(dplyr::select(data_fuse, Sex, PEDSEX, SNPSEX), 1, function(ind) {
                  # Mismatch Sex in pheno and geno files
                  if (ind[1] == externalFemale && ind[2] %in% c(0, 1)) {
                      return(ifelse(ind[3] == 2, FALSE, TRUE))
                  }
                  if (ind[1] == externalMale && ind[2] %in% c(0, 2)) {
                      return(ifelse(ind[3] == 1, FALSE, TRUE))
                  }
                  # Mismatch Sex assigned and estimated
                  if (ind[1] == externalFemale && ind[2] == 2) {
                      return(ifelse(ind[3] == 1, TRUE, NA))
                  }
                  if (ind[1] == externalMale && ind[2] == 1) {
                      return(ifelse(ind[3] == 2, TRUE, NA))
                  }
            })
    # SNPSEX != (Sex in pheno file, PEDSEX)
    fail_sex <-
       dplyr::select(sexcheck_fuse, FID, IID, Sex, PEDSEX,
                     SNPSEX, F)[which(sex_mismatch),]
    # SNPSEX == Sex in pheno file and PEDSEX != (SNPSEX, Sex)
    mixup_geno_pheno <-
       dplyr::select(sexcheck_fuse, FID, IID, Sex, PEDSEX,
                     SNPSEX, F)[which(!sex_mismatch),]
    # Fix mismatch between PEDSEX and sex
    if (fixMixup) {
        if (nrow(mixup_geno_pheno) != 0) {
            file_mixup <- paste(qcdir, "/", alg, ".mismatched_sex_geno_pheno",
                                sep="")
            write.table(dplyr::select(mixup_geno_pheno, FID, IID, SNPSEX),
                       file=file_mixup,
                       row.names=FALSE, quote=FALSE, col.names=FALSE)
            system(paste("plink --bfile ", qcdir, "/", alg, "\ ",
                               "--update-sex ", file_mixup, "\ ",
                               "--make-bed \ ",
                               "--out ", qcdir , "/", alg, sep=""))
        } else {
            if (verbose) {
                message("All assigned genotype sexes (PEDSEX) match external",
                        " sex assignment (Sex)")
            }
            mixup_geno_pheno <- NULL
        }
    } else {
        # Append mismatch to fail_sex, ie treat as fail IDs
        if (nrow(mixup_geno_pheno) != 0) {
            fail_sex <- rbind(fail_sex, mixup_geno_pheno)
        } else {
            mixup_geno_pheno <- NULL
        }
    }
    if (plot) {
        sexcheck$color <- "#1b9e77"
        data$color[sexcheck$PEDSEX==2] <- "#d95f02"
        plot(sexcheck$PEDSEX, sexcheck$F, axes=FALSE, pch=20, main="Sex Check",
             xlim=c(0,3),
             col=sexcheck$color, xlab="Reported Gender", ylab="ChrX Inbreeding")
        axis(1, at=c(1,2), labels=c("Male","Female"))
        axis(2, at=c(maleTh, femaleTh), labels=c(maleTh,femaleTh))
        segments(0.5, maleTh, 1.5, maleTh, lty=2, col="red")
        segments(1.5, femaleTh, 2.5, femaleTh, lty=2, col="red")
        maptools::pointLabel(x=fail_sex$PEDSEX, y=fail_sex$F,
                             labels=fail_sex$IID, cex=.6,
                             offset =.2)
        box()
    }
    return(list(fail_sex=fail_sex, mixup=mixup_geno_pheno))
}

#' Identification of individuals with outlying missing genotype or
#' heterozygosity rates
#'
#' Evaluates and depicts results from plink --missing (missing genotype rates
#' per individual) and plink --het (heterozygosity rates per individual).
#' Non-systematic failures in genotyping and outlying heterozygosity (hz) rates
#' per individual are often proxies for DNA sample quality. Larger than expected 
#' heterozygosity can indicate possible DNA contamination.
#' The mean heterozygosity in Plink is computed as hz_mean = (N-O)/N, where
#' N: number of non-missing genotypes and O:observed number of homozygous
#' genotypes for a given individual.
#' Mean heterozygosity can differ between populations and SNP genotyping panels.
#' Within a population and genotyping panel, a reduced heterozygosity rate can
#' indicate inbreeding - these individuals will then be returned by
#' check_relatedness as individuals that fail the relatedness filters.
#'
#' @param qcdir [character] path/to/directory/with/QC/results containing
#' alg.imiss and alg.het results as returned by plink --missing and plink
#' --het. Users needs writing permission to qcdir.
#' @param alg [character] prefix of QC-ed plink files, ie alg.bed, alg.bim,
#' alg.fam, alg.het and alg.imiss.
#' @param imissTh [double] Threshold for acceptable missing genotype rate in any
#' individual
#' @param hetTh [double] Threshold for acceptable deviation from mean
#' heterozygosity in any individual. Expressed as multiples of standard
#' deviation of heterozygosity (het), i.e. individuals outside mean(het) +/-
#' hetTh*sd(het) will be returned as failing heterozygosity check.
#' @param plot [logical] Should results be depicted as plots? When choosing this
#' option, make sure you have X-forwarding/graphical interface available for
#' interactive plotting. Alternatively, enclose calls to this function in
#' pdf(outfile) check_heterozygosity_and_missingness(arguments) dev.off().
#' @return named [list] with i) fail_imiss dataframe containing FID (Family ID),
#' IID (Within-family ID), MISS_PHENO (Phenotype missing? (Y/N)), N_MISS (Number
#' of missing genotype call(s), not including obligatory missings), N_GENO (
#' Number of potentially valid call(s)), F_MISS (Missing call rate) of
#' individuals failing missing genotype check and ii) fail_het dataframe
#' containing FID (Family ID), IID (Within-family ID), O(HOM) (Observed number
#' of homozygotes), E(HOM) (Expected number of homozygotes), N(NM) (Number of
#' non-missing autosomal genotypes), F (Method-of-moments F coefficient
#' estimate) of individuals failing  outlying heterozygosity check.
check_heterozygosity_and_missingness <- function(qcdir, alg, imissTh=0.03,
                                                 hetTh=3, plot=FALSE) {

    if (!file.exists(qcdir,"/", alg, ".imiss",sep="")){
        stop("plink --missing output file: ", qcdir,"/", alg,
             ".imiss does not exist.")
    }
    if (!file.exists(qcdir,"/", alg, ".het",sep="")){
        stop("plink --het output file: ", qcdir,"/", alg,
             ".het does not exist.")
    }
    names_imiss <- c("FID", "IID", "MISS_PHENO", "N_MISS", "N_GENO", "F_MISS")
    imiss <- read.table(paste(qcdir, "/", alg, ".imiss", sep=""), header=TRUE,
                        as.is=TRUE)
    if (!all(names_imiss) == names(imiss)) {
        stop("Header of", qcdir,"/", alg, ".imiss is not correct. Was your
             file generated with plink --imiss?")
    }
    fail_imiss <- imiss[imiss$F_MISS > imissTh,]

    names_het <- c("FID", "IID", "O(HOM)", "E(HOM)", "N(NM)", "F")
    het <- read.table(paste(qcdir, "/", alg, ".het", sep=""), header=TRUE)
    if (!all(names_het) == names(het)) {
        stop("Header of", qcdir,"/", alg, ".het is not correct. Was your
             file generated with plink --het?")
    }
    fail_het <- het[het$F < (mean(het$F)  - hetTh*sd(het$F)) |
                     het$F > (mean(het$F) + hetTh*sd(het$F)),]

    if (plot) {
        nr_samples <- nrow(imiss)
        imiss$logF_MISS <- log10(imiss[,6])
        data <- merge(imiss, het, by="IID")
        F_minus_5sd <- mean(data$F) - (5*sd(data$F))
        F_add_5sd <- mean(data$F) + (5*sd(data$F))
        F_minus_4sd <- mean(data$F) - (4*sd(data$F))
        F_add_4sd <- mean(data$F) + (4*sd(data$F))
        F_minus_3sd <- mean(data$F) - (3*sd(data$F))
        F_add_3sd <- mean(data$F) + (3*sd(data$F))
        F_minus_2sd <- mean(data$F) - (2*sd(data$F))
        F_add_2sd <- mean(data$F) + (2*sd(data$F))
        F_minus_1sd <- mean(data$F) - (1*sd(data$F))
        F_add_1sd <- mean(data$F) + (1*sd(data$F))
        plot(data$logF_MISS,data$F, xlim=c(-4,0), pch=20,
             col=densCols(data$logF_MISS, data$F),
             main="Heterozygosity by Missingness for All samples",
             xlab="Proportion of missing SNPs",
             ylab="Heterozygosity (SD from mean)", axes=FALSE)
        box()
        fail_het_imiss <- data[which(data[,1] %in%
                                     union(fail_het[,1], fail_imiss[,1])),]
        maptools::pointLabel(x=fail_het_imiss$logF_MISS, y=fail_het_imiss$F,
                   labels=fail_het_imiss[,1], cex=0.6, offset=.2)
        axis(1, at=c(-4,-3,-2, log10(0.03), log10(0.05),-1,0),
             labels=c(0.0001, 0.001, 0.01, 0.03, 0.05, 0.01, 1))
        axis(2,
             at=c(F_minus_5sd, F_minus_4sd, F_minus_3sd,
                  F_add_3sd, F_add_4sd, F_add_5sd),
             labels=c("-5", "-4", "-3" ,"+3", "+4", "+5"), las=2)
        abline(h=F_minus_5sd, col="azure4", lty=3)
        abline(h=F_add_5sd, col="azure4", lty=3)
        abline(h=F_minus_4sd, col="azure4", lty=3)
        abline(h=F_add_4sd, col="azure4", lty=3)
        abline(h=F_minus_3sd, col="azure4", lty=3)
        abline(h=F_add_3sd, col="azure4", lty=3)
        abline(h=F_minus_2sd, col="azure4", lty=3)
        abline(h=F_add_2sd, col="azure4",lty=3)
        abline(h=F_minus_1sd, col="azure4",lty=3)
        abline(h=F_add_1sd, col="azure4",lty=3)
        abline(h=mean(data$F) - (hetTh*sd(data$F)), col="red", lty=2)
        abline(h=mean(data$F) + (hetTh*sd(data$F)), col="red", lty=2)
        abline(v=log10(0.03), col="azure4", lty=3)
        abline(v=log10(0.05), col="azure4", lty=3)
        abline(v=log10(imissTh), col="red", lty=2)
    }
    return(list(fail_imiss=fail_imiss,
                fail_het=fail_het))
}

#' Identification of related individuals
#'
#' Evaluates and depicts results of plink --genome.
#' plink --genome calculates identity by state (IBS) for each pair of individuals
#' based on the average proportion of alleles shared at genotyped SNPs. The
#' degree of recent shared ancestry, i.e. the identity by descent (IBD) can
#' be estimated from the genome-wide IBS. The proportion of IBD between two
#' individuals is returned by --genome as PI_HAT.
#' check_relatedness uses two external scripts for parsing of the alg.genome
#' files to find pairs of samples whose proportion of IBD is larger than the
#' specified highIBDTh. Subsequently, the individual with the greater genotype
#' missingness rate is selected and returned as the individual failing the
#' relatedness check. Optionally, a famfile with additional
#' relatedness information can be provided.
#'
#' @param qcdir [character] path/to/directory/with/QC/results containing
#' alg.imiss and alg.genome results as returned by plink --missing and plink
#' --genome. Users needs writing permission to qcdir.
#' @param alg [character] prefix of QC-ed plink files, ie alg.bed, alg.bim,
#' alg.fam, alg.genome and alg.imiss.
#' @param famfile [character] path/to/file/with/additional/family/information.
#' @param highIBDTh [double] Threshold for acceptable proportion of IBD between
#' pair of individuals.
#' @param plot [logical] Should results be depicted as plots? When choosing this
#' option, make sure you have X-forwarding/graphical interface available for
#' interactive plotting. Alternatively, enclose calls to this function in
#' pdf(outfile) check_relatedness(arguments) dev.off().
#' @return [data.frame] with 

check_relatedness <- function(qcdir, alg, famfile=NULL, highIBDTh=0.1875,
                              plot=FALSE) {
    if (!file.exists(qcdir,"/", alg, ".imiss",sep="")){
        stop("plink --missing output file: ", qcdir,"/", alg,
             ".imiss does not exist.")
    }
    if (!file.exists(qcdir,"/", alg, ".genome",sep="")){
        stop("plink --genome output file: ", qcdir,"/", alg,
             ".genome does not exist.")
    }
    names_imiss <- c("FID", "IID", "MISS_PHENO", "N_MISS", "N_GENO", "F_MISS")
    imiss <- read.table(paste(qcdir, "/", alg, ".imiss", sep=""), header=TRUE,
                        as.is=TRUE)
    if (!all(names_imiss) == names(imiss)) {
        stop("Header of", qcdir,"/", alg, ".imiss is not correct. Was your
             file generated with plink --imiss?")
    }
    genome <- c("FID", "IID", "MISS_PHENO", "N_MISS", "N_GENO", "F_MISS")
    genome <- read.table(paste(qcdir, "/", alg, ".genome", sep=""), header=TRUE,
                        as.is=TRUE)
    if (!all(names_genome) == names(genome)) {
        stop("Header of", qcdir,"/", alg, ".genome is not correct. Was your
             file generated with plink --genome?")
    }
    if (!is.null(famfile)) {
        system(paste("~/GWAS/analysis/genotyping/QC/highIBD_relatedness.pl ",
                     "--file ", qcdir, "/", alg," --thres ", highIBDTh,
                    " --famfile ", famfile, sep=""), wait=TRUE)
    } else {
        system(paste("~/GWAS/analysis/genotyping/QC/highIBD_relatedness_Anderson.pl  ",
                     qcdir, "/", alg, sep=""))
    }
    is.IBD =  system(paste("cat ", qcdir, "/", alg, ".fail-IBD.IDs | wc -l",
                               sep=""), intern=TRUE)
    if (is.IBD != 0 ) {
        fail_highIBD <- read.table(paste(qcdir,"/", alg, ".fail-IBD.IDs",
                                         sep=""))
    } else {
        fail_highIBD  <- NULL
    }
    if (plot) {
        hist(genome$PI_HAT, col="darkblue",
             main="Estimated IBD (PI_HAT) for All Pairs",
             xlab="Estimated pairwise IBD", breaks=c(0.0125*0:80),
             ylab="# of pairs")
        abline(v=highIBDTh, lty=2, col="red")
        genome <- genome[genome$PI_HAT > 0.05,]
        hist(genome$PI_HAT,col="darkblue",
             main="Estimated IBD (PI_HAT) for pairs >0.05",
             xlab="Estimated pairwise IBD", breaks=c(0.0125*0:80),
             ylab="# of pairs")
            abline(v=highIBDTh, lty=2, col="red")
    }
    return(fail_IBD=fail_highIBD)
}


#' Identification of individuals of divergent ancestry
#'
#' Evaluates and depicts results of plink --pca on merged genotypes from
#' individuals to be QCed and individuals of reference population of known
#' genotypes. Currently, check ancestry only supports automatic selection of
#' individuals of European descent. It uses information from principal
#' components 1 and 2 returned by plnik --pca to find the center of the European
#' reference samples (mean(PC1_europeanRef), mean(PC2_europeanRef). It then
#' computes the maximum euclidean distance (maxDist) of the European reference
#' samples from this centre. All study samples whose euclidean distance from the
#' centre falls outside the circle described by the radius r=europeanTh* maxDist
#' are considered non-European and their IDs are returned as failing the
#' ancestry check.
#'
#' @param qcdir [character] path/to/directory/with/QC/results containing
#' prefixMergedDataset.eigenvec results as returned by plink --pca.
#' User needs writing permission to qcdir.
#' @param alg [character] prefix of QC-ed plink files, ie alg.bed, alg.bim,
#' alg.fam.
#' @param prefixMergedDataset [character] Prefix of merged dataset (study and
#' reference samples) used in plink --pca, resulting in
#' prefixMergedDataset.eigenvec
#' @param europeanTh [double] Scaling factor of radius to be drawn around center
#' of reference European samples, with study samples inside this radius
#' considered to be of European descent and samples outside this radius of
#' non-European descent. The radius is computed as the maximum Euclidean distance
#' of reference Europeans to the centre of reference European samples.
#' @param refSamplesFile [character] /Path/to/File/with/reference samples. Needs
#' columns with sample identifiers [refSamplesIID] corresponding to IIDs in
#' prefixMergedDataset.eigenvec and population identifier [refSamplesPop]
#' corresponding to population IDs [refColorsPop] in refColorsfile/refColors.
#' @param refColorsFile [character, optional]
#' /Path/to/File/with/Population/Colors; if not provided and is.null(refColors)
#' ggplot default colors for qualitative scale are used if plot is TRUE.
#' @param refSamplesIID [character] Column name of reference sample IDs in
#' refSamples/refSamplesFile.
#' @param refSamplesPop [character] Column name of reference sample population
#' IDs in refSamples/refSamplesFile.
#' @param refColorsColor [character] Column name of population colors in
#' refColors/refColorsFile
#' @param refColorsPop [character] Column name of reference sample population
#' IDs in refColors/refColorsFile.
#' @param studyColor [character] Color to be used for study population if plot
#' is TRUE.
#' @param plot [logical] Should results be depicted as plots? When choosing this
#' option, make sure you have X-forwarding/graphical interface available for
#' interactive plotting. Alternatively, enclose calls to this function in
#' pdf(outfile) check_ancestry(arguments) dev.off().
#' @return [data.frame] with FID and IID of non-European individuals.
check_ancestry <- function(qcdir, alg, prefixMergedDataset, europeanTh=1.5,
                           refSamples=NULL, refColors=NULL,
                           refSamplesFile=NULL, refColorsFile=NULL,
                           refSamplesIID="IID", refSamplesPop="Pop",
                           refColorsColor="Color", refColorsPop="Pop",
                           studyColor="#2c7bb6",
                           plot=FALSE) {
    if (!file.exists(qcdir,"/", alg, ".fam",sep="")){
        stop("plink family file: ", qcdir,"/", alg, ".fam does not exist.")
    }
    samples <- data.table::fread(paste(qcdir, "/", alg, ".fam", sep=""),
                                 header=FALSE, stringsAsFactors=FALSE,
                                 data.table=FALSE)[,1:2]
    colnames(samples) <- c("FID", "IID")
    if (!file.exists(qcdir,"/", prefixMergedDataset, ".eigenvec",sep="")){
        stop("plink --pca output file: ", qcdir,"/", prefixMergedDataset,
             ".eigenvec does not exist.")
    }
    pca_data <- data.table::fread(paste(qcdir, "/", prefixMergedDataset,
                                        ".eigenvec", sep=""),
                                  stringsAsFactors=FALSE, data.table=FALSE)
    colnames(pca_data) <- c("FID", "IID", paste("PC",1:(ncol(pca_data)-2),
                                                sep=""))

    if (is.null(refSamples) && !file.exists(refSamplesFile)) {
        stop("refSamples are not specified and refSamplesFile file",
             refSamplesFile, "does not exist".)
    }
    if (!is.null(refSamplesFile)) {
        refSamples <- read.table(refSamplesFile, header=TRUE,
                                 stringsAsFactors=FALSE)
    }
    if (!is.null(refColorsFile) && !file.exists(refColorsFile)) {
        stop("refColorsFile file", refColorsFile, "does not exist".)
    }
    if (!is.null(refColorsFile)) {
        refColors <- read.table(refColorsFile, header=TRUE,
                                 stringsAsFactors=FALSE)
    }
    if (!(refSamplesIID  %in% names(refSamples))) {
        stop(paste("Column", refSamplesIID, "not found in refSamples."))
    }
    if (!(refSamplesPop %in% names(refSamples))) {
        stop(paste("Column", refSamplesPop, "not found in refSamples."))
    }
    names(refSamples)[names(refSamples) == refSampleID] <- "IID"
    names(refSamples)[names(refSamples) == refSamplePop] <- "Pop"

    if (!is.null(refColors)) {
        if (!(refColorsColor  %in% names(refColors))) {
            stop(paste("Column", refColorsColor, "not found in refColors."))
        }
        if (!(refColorsPop %in% names(refColors))) {
            stop(paste("Column", refColorsPop, "not found in refColors."))
        }
        names(refColors)[names(refColors) == refColorsColor] <- "Color"
        names(refColors)[names(refColors) == refColorsPop] <- "Pop"
        refColors <- dplyr::select(refColors, Pop, Color)
    } else {
        refColors <- data.frame(Pop=unique(as.character(refSamples$Pop)))
        refColors$Color <- 1:nrow(refColors)
    }
    refSamples <- merge(refSamples, refColors, by="Pop", all.X=TRUE)

    ## Combine pca data and population information ####
    data_all <- merge(pca_data, refSamples, by="IID", all.x=TRUE)
    data_all$Pop[is.na(data_all$Pop)] <- alg
    data_all$Color[is.na(data_all$Color)] <- samplesColor
    data_all <- data_all[order(data_all$Pop, decreasing=TRUE),]

    refColors <- rbind(refColors,c(alg, samplesColor))
    data_all$Color <- as.factor(data_all$Color)
    data_all$Pop <- factor(data_all$Pop, levels=refColors$Pop)

    ## Find mean coordinates and distances of reference Europeans ####
    all_european <- dplyr::filter(data_all, Pop %in% c("CEU", "TSI"))
    euro_pc1_mean <- mean(all_european$PC1)
    euro_pc2_mean <- mean(all_european$PC2)

    all_european$euclid_dist <- sqrt((all_european$PC1 - euro_pc1_mean)^2 +
                                     (all_european$PC2 - euro_pc2_mean)^2)

    max_euclid_dist <- max(all_european$euclid_dist)

    ## Find samples' distances to reference Europeans ####
    data_alg <- dplyr::filter(data_all, Pop == alg)
    data_alg$euclid_dist <- sqrt((data_alg$PC1 - euro_pc1_mean)^2 +
                             (data_alg$PC2 - euro_pc2_mean)^2)
    fail_ancestry <-
        dplyr::select(data_alg, FID, IID)[which(data_alg$euclid_dist >
                                                (max_euclid_dist * europeanTh))]
    if (plot) {
        p <- ggplot2::ggplot()
        p <- p + ggplot2::geom_point(data=data_all,
                                     ggplot2::aes(x=PC1, y=PC2, color=Pop)) +
            ggplot2::scale_color_manual(values=refColors$Color,
                                        name="Population") +
            ggforce::geom_circle(aes(x0=euro_pc1_mean, y0=euro_pc2_mean,
                                     r=(max_euclid_dist * europeanTh))) +
            ggplot2::theme_bw()
        print(p)
    }
    return(fail_ancestry)
}


#' Identification of SNPs with high missingness rates
#'
#' Uses plink --remove alg.fail.IDs --missing --freq to calculate the rates of
#' missing genotype calls and frequency for all variants in the individuals that
#' passed the perIndividualQC.
#'
#' @param qcdir [character] /path/to/directory/with/QC/results containing
#' alg.bim, alg.bed, alg.fam files and alg.fail.IDs containing IIDs of
#' individuals that failed QC. User needs writing permission to qcdir.
#' @param alg [character] prefix of QC-ed plink files, i.e. alg.bed, alg.bim,
#' alg.fam.
#' @param lmissTh [double] Threshold for acceptable variant missing rate across
#' samples.
#' @param plot [logical] Should results be depicted as plots? When choosing this
#' option, make sure you have X-forwarding/graphical interface available for
#' interactive plotting. Alternatively, enclose calls to this function in
#' pdf(outfile) check_SNP_missingness(arguments) dev.off().
#' @return [data.frame] with CHR (Chromosome code), SNP (Variant identifier),
#' CLST (Cluster identifier. Only present with --within/--family), N_MISS
#' (Number of missing genotype call(s), not counting obligatory missings),
#' N_CLST  (Cluster size; does not include nonmales on Ychr; Only present with
#' --within/--family), N_GENO (Number of potentially valid call(s)), F_MISS 
#' (Missing call rate) for all SNPs failing the lmissTh.
check_SNP_missingness <- function(qcdir, alg, lmissTh=0.01, plot=FALSE) {
    if (!file.exists(qcdir,"/", alg, ".fam",sep="")){
        stop("plink family file: ", qcdir,"/", alg, ".fam does not exist.")
    }
    if (!file.exists(qcdir,"/", alg, ".bim",sep="")){
        stop("plink snp file: ", qcdir,"/", alg, ".bim does not exist.")
    }
    if (!file.exists(qcdir,"/", alg, ".bed",sep="")){
        stop("plink binary file: ", qcdir,"/", alg, ".bed does not exist.")
    }
    if (!file.exists(qcdir,"/", alg, ".fail.IDs",sep="")){
        message("File with individuals that failed perIndividual QC: ",
                qcdir,"/", alg, ".fail.IDs does not exist. Continue ",
                "check_SNP_missingness for all samples in ", qcdir,"/", alg,
                ".fam")
        suffix <- ""
        system(paste("plink --bfile ", qcdir, "/", alg,
                     "--missing ",
                     "--freq ",
                     "--out ", qcdir, "/", alg, suffix, sep=""))
    } else {
        suffix <- ".no_failIDs"
        system(paste("plink --bfile ", qcdir, "/", alg,
                     " --remove ", qcdir ,"/", alg, ".fail.IDs --missing ",
                     "--freq ",
                     "--out ", qcdir, "/", alg, suffix, sep=""))
    }

    lmiss <- read.table(paste(qcdir, "/", alg, suffix, ".lmiss",sep=""),
                       as.is=TRUE, header=TRUE)
    if (plot) {
        frq <- read.table(paste(qcdir, "/", alg, suffix, ".frq", sep=""),
                      header=TRUE, as.is=TRUE)
        lmiss_frq <- merge(lmiss, frq)
        dataLowMAF <- lmiss_frq[lmiss_frq$MAF < 0.05,]
        hist(log10(dataLowMAF$F_MISS), axes=F, xlim=c(-4,0), col="burlywood",
             ylab="Number of SNPs", xlab="Proportion of missing data",
             main="Missingness for SNPs with MAF<0.05")
        axis(side=2, labels=TRUE)
        axis(side=1, labels=FALSE)
        mtext(c("0.0001","0.001", lmissTh, "0.01", "0.1", "1"), side=1,
              at=c(-4, -3, log10(lmissTh), -2, -1, 0), line=1)
        abline(v=log10(lmissTh), lty=2, col="red")
        dataHighMAF <- lmiss_frq[lmiss_frq$MAF >=0.05,]
        hist(log10(dataHighMAF$F_MISS), axes=FALSE, xlim=c(-4,0),
             col="burlywood",
             ylab="Number of SNPs", xlab="Proportion of missing data",
             main="Missingness for SNPs with MAF>=0.05")
        axis(side=2, labels=TRUE)
        axis(side=1, labels=FALSE)
        mtext(c("0.0001","0.001", lmissTh, "0.01", "0.1", "1"), side=1,
              at=c(-4, -3, log10(lmissTh), -2, -1, 0), line=1)
        abline(v=log10(lmissTh), lty=2, col="red")
    }
    fail_missingness <- lmiss[lmiss$F_MISS > lmissTh,]
    return(fail_missingness)
}

#' Identification of SNPs showing a significant deviation from Hardy-Weinberg-
#' equilibrium (HWE)
#'
#' Uses plink --remove alg.fail.IDs --hwe to calculate the observed and expected
#' heterozygote frequencies per SNP in the individuals that
#' passed the perIndividualQC and compute the deviation of the frequencies from
#' Hardy-Weinberg equilibrium (HWE) by HWE exact test.
#' If plot is TRUE, the p-values of the HWE exact test are displayed at
#' histograms, where the hweTh is used to depict the QC cut-off for SNPs.
#'
#' @param qcdir [character] /path/to/directory/with/QC/results containing
#' alg.bim, alg.bed, alg.fam files and alg.fail.IDs containing IIDs of
#' individuals that failed QC. User needs writing permission to qcdir.
#' @param alg [character] prefix of QC-ed plink files, i.e. alg.bed, alg.bim,
#' alg.fam.
#' @param hweTh [double] Significance threshold for deviation from HWE.
#' @param plot [logical] Should results be depicted as plots? When choosing this
#' option, make sure you have X-forwarding/graphical interface available for
#' interactive plotting. Alternatively, enclose calls to this function in
#' pdf(outfile) check_HWE(arguments) dev.off().
#' @return [data.frame] with CHR (Chromosome code), SNP (Variant identifier),
#' TEST (Type of test: one of {'ALL', 'AFF', 'UNAFF', 'ALL(QT)', 'ALL(NP)'}), A1
#' (Allele 1; usually minor), A2 (Allele 2; usually major), GENO ('/'-separated
#' genotype counts: A1 hom, het, A2 hom), O(HET) (Observed heterozygote frequency
#' E(HET) (Expected heterozygote frequency), P (Hardy-Weinberg equilibrium exact
#' test p-value) for all SNPs that failed the hweTh.
check_HWE <- function(qcdir, alg, hweTh=1e-5, plot=FALSE) {
    if (!file.exists(qcdir,"/", alg, ".fam",sep="")){
        stop("plink family file: ", qcdir,"/", alg, ".fam does not exist.")
    }
    if (!file.exists(qcdir,"/", alg, ".bim",sep="")){
        stop("plink snp file: ", qcdir,"/", alg, ".bim does not exist.")
    }
    if (!file.exists(qcdir,"/", alg, ".bed",sep="")){
        stop("plink binary file: ", qcdir,"/", alg, ".bed does not exist.")
    }
    if (!file.exists(qcdir,"/", alg, ".fail.IDs",sep="")){
        message("File with individuals that failed perIndividual QC: ",
                qcdir,"/", alg, ".fail.IDs does not exist. Continue ",
                "check_HWE for all samples in ", qcdir,"/", alg,
                ".fam")
        suffix <- ""
        system(paste("plink --bfile ", qcdir, "/", alg,
                     "--hardy ",
                     "--out ", qcdir, "/", alg, suffix, sep=""))
    } else {
        suffix <- ".no_failIDs"
        system(paste("plink --bfile ", qcdir, "/", alg,
                     " --remove ", qcdir, "/", alg, ".fail.IDs --hardy ",
                     "--out ", qcdir, "/", alg, suffix, sep=""))
    }
    hwe <- read.table(paste(qcdir, "/", alg, suffix, ".hwe", sep=""),
                       header=TRUE, as.is=TRUE)
    hwe <- hwe[grepl("ALL", hwe$TEST),]
    if (plot) {
        hist(-log10(hwe$P), ylab="Number of SNPs",col="palegreen",
             xlab="-log(P_hwe)", main="HWE_Pvalue for All SNPs")
        hwe <- hwe[hwe$P <0.01,]
        hist(-log10(hwe$P), ylab="Number of SNPs", col="palegreen",
             xlab="-log(P_hwe)",
             main="HWE_Pvalue for SNPs with HWE_Pvalue <0.01")
        abline(v=-log10(hweTh), lty=2, col="red")
    }
    fail_hwe <- hwe[hwe$P < hweTh,]
    return(fail_hwe)
}

#' Identification of SNPs with a very low minor allele frequency
#'
#' Uses plink --remove alg.fail.IDs --freq to calculate the minor allele
#' frequencies for all variants in the individuals that
#' passed the perIndividualQC.
#' If plot is TRUE, the minor allele frequency distributions are plotted as
#' histograms.
#'
#' @param qcdir [character] /path/to/directory/with/QC/results containing
#' alg.bim, alg.bed, alg.fam files and alg.fail.IDs containing IIDs of
#' individuals that failed QC. User needs writing permission to qcdir.
#' @param alg [character] prefix of QC-ed plink files, i.e. alg.bed, alg.bim,
#' alg.fam.
#' @param mafTh [double] Threshold for minor allele frequency cut-off.
#' @param macTh [double] Threshold for minor allele cut cut-off, if both mafTh
#' and macTh are specied, macTh is used (macTh = mafTh*2*NrSamples)
#' @param plot [logical] Should results be depicted as plots? When choosing this
#' option, make sure you have X-forwarding/graphical interface available for
#' interactive plotting. Alternatively, enclose calls to this function in
#' pdf(outfile) check_maf(arguments) dev.off().
#' @return [data.frame] with CHR (Chromosome code), SNP (Variant identifier), A1
#' (Allele 1; usually minor), A2 (Allele 2; usually major), MAF (Allele 1
#' frequency), NCHROBS (Number of allele observations) for all SNPs that failed
#' the mafTh/macTh.
check_maf <- function(qcdir, alg, mafTh=0.01, macTh=20, verbose=FALSE,
                      plot=FALSE) {
    if (!file.exists(qcdir,"/", alg, ".fam",sep="")){
        stop("plink family file: ", qcdir,"/", alg, ".fam does not exist.")
    }
    if (!file.exists(qcdir,"/", alg, ".bim",sep="")){
        stop("plink snp file: ", qcdir,"/", alg, ".bim does not exist.")
    }
    if (!file.exists(qcdir,"/", alg, ".bed",sep="")){
        stop("plink binary file: ", qcdir,"/", alg, ".bed does not exist.")
    }
    if (!file.exists(qcdir,"/", alg, ".fail.IDs",sep="")){
        message("File with individuals that failed perIndividual QC: ",
                qcdir,"/", alg, ".fail.IDs does not exist. Continue ",
                "check_maf for all samples in ", qcdir,"/", alg,
                ".fam")
        suffix <- ""
        system(paste("plink --noweb --bfile ", qcdir, "/", alg,
                     " --freq --make-bed ",
                     "--out ", qcdir, "/", alg, suffix, sep=""))
    } else {
        suffix <- ".no_failIDs"
        system(paste("plink --noweb --bfile ", qcdir, "/", alg,
                     " --remove ", qcdir ,"/", alg, ".fail.IDs --freq ",
                     "--make-bed ",
                     "--out ", qcdir, "/", alg, suffix, sep=""))
    }
    maf <- read.table(paste(qcdir,"/",alg, suffix, ".frq",sep=""),
                       header=TRUE, as.is=TRUE)
    samples <- read.table(paste(qcdir,"/",alg, suffix, ".fam",sep=""),
                          header=TRUE, as.is=TRUE)
    if (is.null(mafTh) && is.null(macTh)) {
        stop("Either mafTh or macTh need to be provided")
    }
    if(!is.null(macTh)) {
        mafTh <- macTh/(2*nrow(samples))
    }
    if (verbose) {
        if (!is.null(mafTh) && !is.null(macTh)) {
            message("Both mafTh and macTh provided, macTh=", macTh,
                    " is used (corresponds to mafTh=", mafTh, ")")
        } else if (!is.null(mafTh)) {
            message("The mafTh is ", mafTh)
        } else {
            message("The macTh is ", macTh," which corresponds to a mafTh=",
                    mafTh, ")")
        }
    }
    if (plot) {
        hist(maf$MAF, col="lightblue1", ylab="Number of SNPs", breaks=50,
             xlab="Minor Allele Frequency",
             main="Minor Allele Frequency (MAF) for All SNPs")
        abline(v=mafTh, lty=2, col="red")
    }
    failed_maf <- maf[maf$MAF < mafTh,]
    return(failed_maf)
}

#' Create QC-ed dataset
#'
#' Remove samples that fail per individual QC
#' Remove markers that fail per marker QC
clean_data <- function(qcdir, alg, MAC, highIBDfilter=FALSE, MAF=NULL) {
    # convert MAC into MAF
    MAF <- mafTh/(2*(nr_samples - nr_fail_samples))
    if (highIBDfilter == "IBD") {
        # i) with sample filtering based on all fail-IDs (including relatedness)
        # and variant filter for missing genotypes, MAF and HWE
        system(paste("plink --bfile ", qcdir, "/", alg,
                     " --remove ", qcdir,"/", alg, ".fail.IDs ",
                     "--hwe ", hweTh,
                     " --geno ", lmissTh,
                     " --maf ", MAF ,
                     " --make-bed --out ", qcdir, "/", alg, ".clean", sep=""))
    }
    # ii) with sample filtering based on fail-IDs not including failed IBD and
    # variant filtering for missing genotypes and HWE (keep related samples and
    # MAF since only part of the cohort; smart relatedness filtering when
    # combining cohorts)
    system(paste("cat ", qcdir ,"/",  alg, ".fail-ancestry.IDs  ",
                 qcdir ,"/",  alg, ".fail-het.IDs ",
                 qcdir ,"/",  alg, ".fail-sexcheck.IDs  ",
                 qcdir ,"/",  alg, ".fail-imiss.IDs | sort | uniq >",
                 qcdir, "/", alg, ".noIBDfail.IDs", sep=""), wait=TRUE)
    system(paste("plink --bfile ", qcdir, "/", alg,
                 " --remove ", qcdir,"/", alg, ".noIBDfail.IDs ",
                 "--hwe ", hweTh,
                 " --geno ", lmissTh,
                 " --make-bed --out ", qcdir, "/", alg, ".clean.related",
                 sep=""))
    return(MAF)
}

## functions ###
