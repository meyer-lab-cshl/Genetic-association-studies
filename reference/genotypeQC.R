###################################################################################################
###     Script based on "Anderson et al. (2010) Data qualtiy control in gentic case-control     ###
###     association studies." Nature protocols 5(9):1564-73                                     ###
###                                                                                             ###
###                                                                                             ###
###         by Hannah Meyer                                                                     ###
###################################################################################################

#############
# libraries #
#############
#suppressPackageStartupMessages(library(optparse))
# library("geneplotter")

library ("R.methodsS3", lib.loc ="/homes/hannah/R/x86_64-pc-linux-gnu-library/3.0")
library ("R.oo", lib.loc ="/homes/hannah/R/x86_64-pc-linux-gnu-library/3.0")
library ("R.utils", lib.loc ="/homes/hannah/R/x86_64-pc-linux-gnu-library/3.0")
library("calibrate", lib.loc ="/homes/hannah/R/x86_64-pc-linux-gnu-library/3.0")
library("sp", lib.loc= "/homes/hannah/R/x86_64-pc-linux-gnu-library/3.0")
library("maptools", lib.loc= "/homes/hannah/R/x86_64-pc-linux-gnu-library/3.0")
library("limma", lib.loc = "/homes/hannah/R/x86_64-pc-linux-gnu-library/3.1")

#############
# functions #
#############

### function parameters ###

# qcdir                  directory to QC data/raw data?
# alg                    snp calling algorithm/flie name
# sexCheckfilter         should sample's sex be checked?  yes: sexCheckfilter == "SEX"
# maleTh                 threshold for homozygosity rates expected on X-chromosome for males
# femaleTh               threshold for homozygosity rates expected on X-chromosome for females
# imissTh                threshold for missing genotypes
# hetTh                  threshold for heterozygosity rate
# C1                     threshold for PC1
# C2                     threshold for PC2
# C3                     threshold for PC3
# highIBDfilter          Should Threshold for High IBD be applied? yes: highIBDfilter == "IBD"
# highIBDTh              Threshold for high IBD 
# crypticIBDTh           Threshold unknown relatedness
# lmissTh                Threshold for missing SNPs
# hweTh                  P-Value threshold for deviation from Hardy-Weinberg-equilibrium
# mafTh                  Threshold for minor allele frequency
# ethnicityFilter        should samples be filter based on ethnicity (Default: yes (caucasians)

#example: alg="gencall";qcdir="/homes/hannah/GWAS/data/genotype/MRI_genotype/QC/gencall.alleth";sexCheckfilter="SEX";maleTh=0.8;femaleTh=0.2;imissTh=0.03;hetTh=3;C1=0.008;C2=0;C3=-0.024;highIBDfilter="IBD";highIBDTh=0.125;
#crypticIBDTh=100;lmissTh=0.01;hweTh=1e-3;mafTh=20;omnix_file="~/GWAS/data/genotype/MRI_genotype/omnix_hhtmri_20141121/gencall_qc/pipeline_summary.csv";ethnicityFilter="no";sample="~/GWAS/data/phenotype/2Dphenotype/20141203GenScan.txt"

#example: alg="gencall";qcdir="/homes/hannah/GWAS/data/genotype/MRI_genotype/QC/gencall.combined";sexCheckfilter="SEX";maleTh=0.8;femaleTh=0.2;imissTh=0.03;hetTh=3;C1=0.008;C2=0;C3=-0.024;highIBDfilter="IBD";highIBDTh=0.125;
#crypticIBDTh=100;lmissTh=0.01;hweTh=1e-3;mafTh=20;omnix_file="~/GWAS/data/genotype/MRI_genotype/omnix_hhtmri_20141121/gencall_qc/pipeline_summary.csv";ethnicityFilter="yes";sample="~/GWAS/data/phenotype/2Dphenotype/20141203GenScan.txt"

#example: alg="gencall";qcdir="/homes/hannah/GWAS/data/genotype/MRI_genotype/QC/gencall.combined.it2";sexCheckfilter="SEX";maleTh=0.8;femaleTh=0.2;imissTh=0.03;hetTh=3;C1=0.008;C2=0;C3=-0.024;highIBDfilter="IBD";highIBDTh=0.125;
#crypticIBDTh=100;lmissTh=0.01;hweTh=1e-3;mafTh=20;omnix_file="~/GWAS/data/genotype/MRI_genotype/omnix_hhtmri_20141121/gencall_qc/pipeline_summary.csv";ethnicityFilter="yes";sample="~/GWAS/data/phenotype/2Dphenotype/20160209_All_BRU_family_format.txt"
#center="sanger"

#example: alg="gencall.singapore12";qcdir="/homes/hannah/GWAS/data/genotype/MRI_genotype/QC/gencall.singapore12";sexCheckfilter="SEX";maleTh=0.8;femaleTh=0.2;imissTh=0.03;hetTh=3;C1=0.005;C2=0;C3=-0.02;highIBDfilter="IBD";highIBDTh=0.125;
#eigenStr="eigenvec";crypticIBDTh=100;lmissTh=0.01;hweTh=1e-3;mafTh=20;omnix_file="~/GWAS/data/genotype/MRI_genotype/QC/gencall.singapore12/gencall.singapore12.raw.fam";ethnicityFilter="yes";sample="~/GWAS/data/phenotype/2Dphenotype/20160209_All_BRU_family_format.txt"
#center="Singapore"

#example: alg="gencall.singapore3";qcdir="/homes/hannah/GWAS/data/genotype/MRI_genotype/QC/gencall.singapore3";sexCheckfilter="SEX";maleTh=0.8;femaleTh=0.2;imissTh=0.03;hetTh=3;C1=0.012;C2=0.01;C3=-0.02;highIBDfilter="IBD";highIBDTh=0.125;
#eigenStr="eigenvec";crypticIBDTh=100;lmissTh=0.01;hweTh=1e-3;mafTh=20;omnix_file="~/GWAS/data/genotype/MRI_genotype/QC/gencall.singapore3/gencall.singapore3.raw.fam";ethnicityFilter="yes";sample="~/GWAS/data/phenotype/2Dphenotype/All_BRU_Hannah_FamilyNo_9.2.16.txt"
#center="Singapore"


### function objection ###
list2counts <- function(element, all_names) {
    all_names[!(all_names %in% element)] <- 0
    all_names[all_names %in% element] <- 1
    return(as.numeric(all_names))
}
myQC<- function(qcdir,alg,sexCheckfilter,maleTh,femaleTh,imissTh,hetTh,C1, C2, C3, highIBDfilter,highIBDTh,crypticIBDTh,lmissTh,hweTh,mafTh,omnix_file, ethnicityFilter, sample, center, eigenStr) {

    print(c(alg,sexCheckfilter,maleTh,femaleTh,imissTh,hetTh,C1, C2, C3, highIBDfilter,highIBDTh,lmissTh,hweTh,mafTh, omnix_file,ethnicityFilter, sample, center))

    cat("Read sample description file (generated by ID_mapping.R\n")
    SampleID <- read.table(file=sample, sep="\t", header=TRUE, quote="", stringsAsFactors=FALSE, na.strings=c("NA",""))[,c("Bru.Number","Ethnicity", "Sex")]
    if (any(duplicated(SampleID$Bru.Number))) {
        SampleID <- SampleID[-which(duplicated(SampleID$Bru.Number)),]
    }
    SampleID$Bru.Number <- toupper(SampleID$Bru.Number)
    if(center == "sanger") {
        omnix <- read.table(omnix_file, header=TRUE, sep=",")[,c("supplier_name", "sample")]
        omnix$sample <- gsub(".*_", "", omnix$sample)
        omnix$supplier_name <- toupper(omnix$supplier_name)
        eigenStr="eigenvec"
    } else {
        omnix <- data.frame(sample=read.table(omnix_file, header=FALSE, stringsAsFactors=FALSE)[,1])
        omnix$sample <- gsub("[0-9]*\\_([0-9A-Z]*)\\_.*", "\\1", omnix$sample)
        omnix$supplier_name=omnix$sample
        eigenStr="eigenvec"
    }

    ID_red <- SampleID[which(SampleID[,1] %in% omnix[,1]),]
    omnix_red <- omnix[which(omnix[,1] %in% ID_red[,1]),]
    omnix_match <- omnix_red[match(ID_red[,1], omnix_red[,1]),]

    # write ethniciies
    ID <- data.frame(omnix_match[,2], ID_red)
    ID$Ethnicity <- gsub("CH", "Chinese", ID$Ethnicity)
    ID$Ethnicity <- gsub("ASC", "Asian Subcontinent", ID$Ethnicity)
    ID$Ethnicity <- gsub("AC", "Afro-caribbean", ID$Ethnicity)
    ID$Ethnicity <- gsub("AF", "African", ID$Ethnicity)
    ID$Ethnicity <- gsub("JAP", "Japanese", ID$Ethnicity)
    ID$Ethnicity <- gsub("M", "Mixed", ID$Ethnicity)
    ID$Ethnicity <- gsub("O_U", "Other/Unknown", ID$Ethnicity)
    ID$Ethnicity <- gsub("^C$", "Caucasian", ID$Ethnicity, perl=TRUE)

    names(ID) <- c("omnix","BRU", "Ethnicity", "Sex")
    ID_file <- paste(qcdir,"/", alg,"_sampleID.txt", sep="")
    write.table(ID, file=ID_file, sep="\t",quote=FALSE, row.names=FALSE)

    pdf(paste(qcdir,"/", alg,".pdf",sep=""),width=10,height=8)


    ####################################
    ###### Per-individual QC ###########
    ####################################

    cat("i) identification of individuals with discordant sex information\n")
    
    ##### i) identification of individuals with discordant sex information
    ##### --> highlight plating errors and sample mix-up #####
    ##### --> calculate homozygosity rates across all X-chromosome SNPs and compare with expected rates: male -> rate = 1, female -Y rate < .2 ####

    if (sexCheckfilter =="SEX") {
        fail_sex <- NULL
        data <- read.table(paste(qcdir,"/",alg,".sexcheck",sep=""),header=T)
        data_fuse <- merge(data, ID, by.x=1, by.y=1)
        gender_mixup <- data_fuse[c(intersect(which(data_fuse$Sex == "F"), which(data_fuse$PEDSEX ==1)),intersect(which(data_fuse$Sex =="M"), which(data_fuse$PEDSEX ==2)), which(data_fuse$PEDSEX ==0)),]
        if(nrow(gender_mixup) != 0) {
            for( i in 1:nrow(gender_mixup)) {
                if ((data_fuse$Sex[i] == "F" && data_fuse$SNPSEX[i] == 2)  || (data_fuse$Sex[i] == "M" && data_fuse$SNPSEX[i] == 1)) {
                    data_fuse[which(data_fuse[,1] %in% gender_mixup[i,1]),3] <- data_fuse[which(data_fuse[,1] %in%  gender_mixup[i,1]),4]
                    system(paste("perl -i -e 'while(<>) {chomp $_; if ($_=~m/^", as.vector(gender_mixup[i,1]),"/) {@array = split \" \", $_; $array[4] =", data_fuse[which(data_fuse[,1] %in% gender_mixup[i,1]),4],"; $_= join(\" \", @array) } print $_, \"\\n\"}' ", qcdir,"/", alg, ".fam" , sep=""))
                } else if (data_fuse$PEDSEX[i] == 0) {
                    if ((data_fuse$Sex[i] == "F" && data_fuse$SNPSEX[i] == 2)  || (data_fuse$Sex[i] == "M" && data_fuse$SNPSEX[i] == 1)) {
                        data_fuse[which(data_fuse[,1] %in% gender_mixup[i,1]),3] <- data_fuse[which(data_fuse[,1] %in%  gender_mixup[i,1]),4]
                        system(paste("perl -i -e 'while(<>) {chomp $_; if ($_=~m/^",gender_mixup[i,1],"/) {@array = split \" \", $_; $array[4] =", data_fuse[which(data_fuse[,1] %in% gender_mixup[i,1]),4],"; $_= join(\" \", @array) } print $_, \"\\n\"}' ", qcdir,"/", alg, ".fam" , sep=""))
                    } else {
                        fail_sex <- rbind(fail_sex, data_fuse[i,])
                    }
                } else {
                    fail_sex <- rbind(fail_sex, data_fuse[i,])
                }
            }
            system(paste("plink --bfile ", qcdir, "/", alg, " --check-sex --out ", qcdir, "/",  alg, "_secondsexcheck; mv ",qcdir,"/",alg, "_secondsexcheck.log ", qcdir,"/plink_log/", alg, "_secondsexcheck.log; grep 'PROBLEM' ", qcdir, "/",  alg, "_secondsexcheck.sexcheck > ", qcdir, "/",alg,"_secondsexcheck.failsex", sep=""))
            data <- read.table(paste(qcdir,"/",alg,"_secondsexcheck.sexcheck",sep=""),header=T)
            data_fuse <- merge(data, ID, by.x=1, by.y=1)
        }
        fail_sex <- rbind(fail_sex,  data_fuse[c(intersect(which(data_fuse[,3] == 1),which(data_fuse[,6]  < maleTh)), intersect(which(data_fuse[,3] ==2), which(data_fuse[,6] > femaleTh))),])
        data$color="blue"
        data$color[data$PEDSEX==2] <- "pink"
        plot(data$PEDSEX, data$F, axes=F, pch=20, main="Sex Check", xlim=c(0,3), col=data$color, xlab="Reported Gender", ylab="ChrX Inbreeding")
        axis(1,at=c(1,2),labels=c("Male","Female"))
        axis(2,at=c(maleTh,femaleTh),labels=c(maleTh,femaleTh))
        segments(0.5,maleTh,1.5,maleTh,lty=2,col="red")
        segments(1.5,femaleTh,2.5,femaleTh,lty=2,col="red")
        pointLabel(fail_sex[,3], fail_sex[,6], fail_sex[,1], cex=.6, offset =.2)
        box()
        write.table(fail_sex[,1:2], file=paste(qcdir,"/",alg, ".fail-sexcheck.IDs", sep=""), append=FALSE, quote=FALSE, row.names=FALSE, col.names=FALSE)

    }

    cat("ii) identification of individuals with outlying missing genotype or heterozygosity rates\n")

    ##### ii) identifaction of individuals with outlying missing genotype or heterozygosity rates #####
    ##### --> genotype failure and heterozygosity (hz) rate per individual are meassure for DNA sample quality ######
    ##### --> hz genotypes in individuals: excessive -> DNA contamination?, reduced -> inbreeding? #####[Ma.
    ##### --> mean hz = (N-O)/N, N: number of non-missing genotypes, O:  observed number of homozygous genotypes for a given individual ####
    ##### --> mean hz differs between populations and SNP genotyping panels ####

    imiss=read.table(paste(qcdir,"/",alg,".imiss",sep=""),header=T,as.is=T)

    nr_samples = nrow(imiss)
    fail_imiss <- imiss[imiss$F_MISS > imissTh,]
     # write table with IDs of missing genotypes
    write.table(fail_imiss[,1:2], file=paste(qcdir,"/",alg,".fail-imiss.IDs",sep=""),append=F,quote=F,row.names=F,col.names=F)
    imiss$logF_MISS = log10(imiss[,6])
    het=read.table(paste(qcdir,"/",alg,".het",sep=""),header=T)
    data = merge(imiss, het)
    F_minus_5sd = mean(data$F)-(5*sd(data$F))
    F_add_5sd = mean(data$F)+(5*sd(data$F))
    F_minus_4sd = mean(data$F)-(4*sd(data$F))
    F_add_4sd = mean(data$F)+(4*sd(data$F))
    F_minus_3sd = mean(data$F)-(3*sd(data$F))
    F_add_3sd = mean(data$F)+(3*sd(data$F))
    F_minus_2sd = mean(data$F)-(2*sd(data$F))
    F_add_2sd = mean(data$F)+(2*sd(data$F))
    F_minus_1sd = mean(data$F)-(1*sd(data$F))
    F_add_1sd = mean(data$F)+(1*sd(data$F))
    fail_het <- data[data$F <(mean(data$F) -hetTh*sd(data$F)) | data$F > (mean(data$F) + hetTh*sd(data$F)),]
    # write table of IDs with outlying heterozygosity rates 
    write.table(fail_het[,1:2], file=paste(qcdir,"/", alg,".fail-het.IDs",sep=""),append=F,quote=F,row.names=F,col.names=F)
    plot(data$logF_MISS,data$F, xlim=c(-4,0),pch=20, col=densCols(data$logF_MISS,data$F), main="Heterozygosity by Missingness for All samples", xlab="Proportion of missing SNPs", ylab="Heterozygosity (SD from mean)",axes=F)
    box()
    fail_het_imiss <- data[which(data[,1] %in% union(fail_het[,1], fail_imiss[,1])),]
    pointLabel(x=fail_het_imiss$logF_MISS,y=fail_het_imiss$F, labels=fail_het_imiss[,1], cex=0.6, offset=.2)
    axis(1,at=c(-4,-3,-2, log10(0.03), log10(0.05),-1,0),labels=c(0.0001,0.001,0.01,0.03, 0.05,0.01, 1))
    axis(2,at=c(F_minus_5sd,F_minus_4sd,F_minus_3sd,F_add_3sd,F_add_4sd,F_add_5sd),labels=c("-5","-4","-3","+3","+4","+5"),las=2)
    abline(h=F_minus_5sd, col="azure4",lty=3)
    abline(h=F_add_5sd, col="azure4",lty=3)
    abline(h=F_minus_4sd, col="azure4",lty=3)
    abline(h=F_add_4sd, col="azure4",lty=3)
    abline(h=F_minus_3sd, col="azure4",lty=3)
    abline(h=F_add_3sd, col="azure4",lty=3)
    abline(h=F_minus_2sd, col="azure4",lty=3)
    abline(h=F_add_2sd, col="azure4",lty=3)
    abline(h=F_minus_1sd, col="azure4",lty=3)
    abline(h=F_add_1sd, col="azure4",lty=3)
    abline(h=mean(data$F)-(hetTh*sd(data$F)), col="red",lty=2)
    abline(h=mean(data$F)+(hetTh*sd(data$F)),col="red",lty=2)
    abline(v=log10(0.03), col="azure4", lty=3)
    abline(v=log10(0.05), col="azure4", lty=3)
    abline(v=log10(imissTh), col="red", lty=2)


    cat("iii) identification of duplicated or related (cryptic) individuals\n")
    
    ##### iii) identification of duplicated or related (cryptic) individuals #####
    ##### sample ancestry: samples should be unrelated --> maximum relatedness between any pair less than that second-degree relative  #####
    ##### --> if present in cohort, cohort allele frequency not a fair reflection of populations allele frequency #####
    ##### --> calculate Identity by state (IBS) for each pair of individuals based on the average proportion of alleles shared at genotyped SNPs 
    ##### --> use only independent SNPs: remove regions of extended LD and prune other regions i.e. no pair of SNP w/i given window is correlated (r2 >2)
    ##### --> duplicates: IBS = 1
    ##### --> degree of recent shared ancestry = identity by descent (IBD) can be estimated from genome wide IBS (using PLINK)
    ##### --> IBD >0.98: duplicates
    ##### sample pair-wise IBD #####

    data <- read.table(paste(qcdir,"/",alg,".genome",sep=""), header=T, as.is=T)
    hist(data$PI_HAT,col="darkblue", main="Estimated IBD (PI_HAT) for All Pairs", xlab="Estimated pairwise IBD", breaks=c(0.0125*0:80), ylab="# of pairs")
    abline(v=highIBDTh,lty=2,col="red")
    data <- data[data$PI_HAT >0.05,]
    hist(data$PI_HAT,col="darkblue", main="Estimated IBD (PI_HAT) for pairs >0.05", xlab="Estimated pairwise IBD", breaks =c(0.0125*0:80), ylab="# of pairs")
    if (highIBDfilter == "IBD") {
        abline(v=highIBDTh,lty=2,col="red")
        system(paste("~/GWAS/analysis/genotyping/highIBD_relatedness.pl --file ", qcdir, "/", alg," --thres ", highIBDTh, " --famfile /homes/hannah/GWAS/data/phenotype/2Dphenotype/20160209_All_BRU_family_format.txt", sep=""), wait=TRUE)
        is.IBD =  system(paste("cat ", qcdir, "/", alg, ".fail-IBD.IDs | wc -l", sep=""), intern=TRUE)
    } else {
        is.IBD=0
    }
    if (is.IBD != 0 ) {
        fail_highIBD <- read.table(paste(qcdir,"/", alg, ".fail-IBD.IDs", sep=""))
        write.table(ID[which(ID[,1] %in% fail_highIBD[,1]),], file=paste(qcdir,"/", alg,".fail_IBD.txt", sep=""), row.names=FALSE, append=FALSE, quote=FALSE, col.names=TRUE, sep="\t")
    } else {
        fail_highIBD  <- NULL
    }

    # write table of cryptic ID
    if (crypticIBDTh != 0) {
        system(paste("awk '$10>0.1 {print $1,$2 \"\\n\" $3,$4}' ", qcdir,"/",alg, ".genome | sort | uniq -c | sort -nr -k1,1 | awk '$1 >", crypticIBDTh  ," {print $2,$3}' > ", alg, ".fail-crypticIBD.IDs", sep=""))
        is.crypticIBD =  system(paste("cat ", qcdir, "/", alg, ".fail-crypticIBD.IDs | wc -l", sep=""), intern=TRUE)
    } else {
        is.crypticIBD=0
    }
    if (is.crypticIBD != 0) {
        fail_crypticIBD <- read.table(paste(qcdir,"/", alg, ".fail-crypticIBD.IDs", sep=""))
        write.table(ID[which(ID[,2] %in% fail_crypticIBD[,1]),], file=paste(qcdir,"/", alg,".fail_crypticIBD.txt", sep=""), row.names=FALSE, append=FALSE, quote=FALSE, col.names=TRUE, sep="\t")
    } else {
        fail_crypticIBD <- NULL
    }

    cat("iv) identification of individuals of divergent ancestry\n")
    ##### iv) identification of individuals of divergent ancestry #####
    ##### --> done in PLINK --genome, followed by smartpca.pl of combined HapMapIII reference set and sample cohort

    system(paste("Rscript ~/GWAS/analysis/genotyping/pca.R --vanilla --default-packages=R.utils --qcdir=", qcdir, " --pcafile=", alg, ".HapMapIII.pruned.pca.", eigenStr, " --alg=", alg, " --C1=",C1, " --C2=",C2, " --C3=", C3, " --IDfile=", ID_file, sep=""), wait=TRUE)
    data <-read.table(file=paste(qcdir,"/", alg, ".pca", sep=""))
    
    # ancestry object contains list of sample IDs and PCA determined ethnicity
    ancestry <-  readRDS(paste(qcdir, "/", alg, ".ancestry.rds", sep=""))
    
    # samples failing european ancestry
    samples_fail_ancestry <- matrix(nc=2, rep(data[union(which(data$PC1 < C1), which(data$PC2 < C2)),1],2))
    
    ## overview QC fails     
    fail_list =  list(missing_genotype=as.vector(fail_imiss[,1]), highIBD=as.vector(fail_highIBD[,1]), crypticIBD=as.vector(fail_crypticIBD[,1]), outlying_heterozygosity=as.vector(fail_het[,1]), mismatched_sex=as.vector(fail_sex[,1]), diff_ancestry= as.vector(samples_fail_ancestry[,1]))
    fail_list = fail_list[!sapply(fail_list, is.null)]
    unique_samples_fail_all <- unique(unlist(fail_list))
   
    # a) overview QC fails independent of ethnicity
    fail_list_wo_ancestry = fail_list[!names(fail_list) == "diff_ancestry"] 
    unique_samples_fail_wo_ancestry <- unique(unlist(fail_list_wo_ancestry))
    
    fail_counts_wo_ancestry <- sapply(fail_list_wo_ancestry,list2counts, unique_samples_fail_wo_ancestry)
    rownames(fail_counts_wo_ancestry) <- unique_samples_fail_wo_ancestry              

    vennDiagram(fail_counts_wo_ancestry, names= names(fail_list_wo_ancestry), main="QC Fail overview (all ethnicities)")
    
    # b) overview of  QC and ancestry fails 
    fail_counts_all <- sapply(list(QC_fail= unique_samples_fail_wo_ancestry, Ancestry_fail= fail_list$diff_ancestry),list2counts, unique_samples_fail_all)
    rownames(fail_counts_all) <- unique_samples_fail_all              
    
    vennDiagram(fail_counts_all, names= names(fail_counts_all), main="Overlap between 'per individual QC' and ancestry analysis")


    # c) overview QC fails per ethnicity 
    fail_list_ancestries <- lapply(ancestry[!names(ancestry) == "nonmatching"], function(eth, qc_list) {
                                        qc_per_eth <- data.frame(sapply(qc_list, list2counts, eth))
                                        rownames(qc_per_eth) <-eth
                                        qc_per_eth_red <- qc_per_eth[!rowSums(qc_per_eth) == 0,]
                                       return(qc_per_eth_red)
                            }, qc_list= fail_list_wo_ancestry)

    plot_venn <- lapply(seq_along(fail_list_ancestries), function(eth) {vennDiagram(fail_list_ancestries[[eth]], names=colnames(fail_list_ancestries[[eth]]), main=paste("QC Fail overview (",names(fail_list_ancestries)[eth], ")"))})  

    write_fail_reason <- lapply(seq_along(fail_list_ancestries), function(index, ID, eth_list) {
                                if (length(eth_list[[index]])[1] !=0) {
                                    eth <- eth_list[[index]]
                                    eth_name <- names(eth_list)[index]
                                    ID_red <- ID[which(ID$omnix %in% rownames(eth)),]
                                    eth <- eth[match(ID_red$omnix, rownames(eth)),]
                                    eth <- cbind(eth, as.vector(ID_red$Ethnicity), as.vector(ID_red$Sex))
                                    eth_tech <- eth[,c("missing_genotype", "outlying_heterozygosity","mismatched_sex")]
                                    eth_tech <- eth_tech[apply(eth_tech, 1, function(row) if (sum(as.numeric(row)) == 0) { return(FALSE) } else {return(TRUE)}),]
                                    eth_tech <- cbind(rownames(eth_tech), eth_tech)
                                    eth_tech <- merge(eth_tech, ID, by=1)#[,c(1,5,2,3,4,6,7)]
                                    colnames(eth_tech)[1] <- "omnix"
                                    write.table(eth, file=paste(qcdir,"/", alg, ".fail_reason_QC_", eth_name,".txt", sep=""), sep="\t",  quote=FALSE, append=FALSE)
                                    write.table(eth_tech, file=paste(qcdir,"/", alg, ".fail_reason_QC_", eth_name,"_tech.txt", sep=""), sep="\t", row.names=FALSE,  quote=FALSE, append=FALSE)
                                }
                            }, ID= ID, eth_list =fail_list_ancestries)

    if (ethnicityFilter !=  "yes") {
        # exclude asians since too few samples
        write.table(cbind(ancestry$asian, ancestry$asian), file=paste(qcdir, "/",alg,".fail-ancestry-asian.IDs",sep=""), append=FALSE, quote=FALSE, row.names=FALSE, col.names=FALSE)
        nr_fail_samples = length(c(ancestry$asian, unique_samples_fail_wo_ancestry))
    }

    if ( ethnicityFilter == "yes") {
        write.table(samples_fail_ancestry, file=paste(qcdir, "/",alg,".fail-ancestry.IDs",sep=""), append=FALSE, quote=FALSE, row.names=FALSE, col.names=FALSE)
        nr_fail_samples = length(unique_samples_fail_all)

        if(nrow(gender_mixup) != 0) {
            gender_mixup_eth <- data[which(data[,1] %in% gender_mixup[,1]),]
            write.table(cbind(gender_mixup[,c(1,2,4,8,9,3)],gender_mixup_eth[,2:3]),  file=paste(qcdir,"/", alg, ".gender_mixup.txt", sep=""), append=FALSE, quote=FALSE, row.names=FALSE, sep="\t")
        }
        if(nrow(fail_sex) !=0) {
            fail_sex_eth <- data[which(data[,1] %in% fail_sex[,1]),]
            write.table(cbind(fail_sex[,c(1,7,3,4,9,6,8)],fail_sex_eth[,2:3]), file=paste(qcdir,"/",alg, ".failsex.txt", sep=""),sep="\t",  append=FALSE, quote=FALSE, row.names=FALSE, col.names=TRUE)
        }

    }


    ####################################
    ###### Per-marker (SNP) QC #########
    ####################################


    cat("###### Per-marker (SNP) QC #########\n")
    system(paste("cat ",qcdir ,"/",  alg, ".fail-*.IDs | sort | uniq >", qcdir, "/", alg, ".fail.IDs", sep=""), wait=TRUE)
    system(paste("plink --noweb --bfile ", qcdir, "/", alg, " --remove ", qcdir ,"/", alg, ".fail.IDs --missing --out ", qcdir, "/", alg, ".no_failID", sep=""))
    system(paste("plink --noweb --bfile ", qcdir, "/", alg, " --remove ", qcdir ,"/", alg, ".fail.IDs --freq --out ", qcdir, "/", alg, sep=""))


    cat("i) identification of SNPs with an excessive missing genotype\n")
    ##### i) identification of SNPs with an excessive missing genotype ####

    ##### SNP missingness histogram #####
    lmiss <-read.table(paste(qcdir,"/",alg,".no_failID.lmiss",sep=""), as.is=T, header=T)
    frq <- read.table(paste(qcdir,"/",alg,".frq",sep=""), header=T,as.is=T)
    data <- merge(lmiss,frq)
    dataLowMAF <- data[data$MAF <0.05,]
    hist(log10(dataLowMAF$F_MISS),axes=F,xlim=c(-4,0),col="burlywood",ylab="Number of SNPs",xlab="Proportion of missing data", main="Missingness for SNPs with MAF<0.05")
    axis(side=2,labels=T)
    axis(side=1,labels=F)
    mtext(c("0.0001","0.001",lmissTh,"0.01","0.1","1"),side=1,at=c(-4,-3,log10(lmissTh),-2,-1,0),line=1)
    abline(v=log10(lmissTh),lty=2,col="red")
    dataHighMAF <- data[data$MAF >=0.05,]
    hist(log10(dataHighMAF$F_MISS),axes=F,xlim=c(-4,0),col="burlywood",ylab="Number of SNPs",xlab="Proportion of missing data", main="Missingness for SNPs with MAF>=0.05")
    axis(side=2,labels=T)
    axis(side=1,labels=F)
    mtext(c("0.0001","0.001",lmissTh,"0.01","0.1","1"),side=1,at=c(-4,-3,log10(lmissTh),-2,-1,0),line=1)
    abline(v=log10(lmissTh),lty=2,col="red")


    cat("ii) identification of SNPs showing a significant deviation from Hardy-Weinberg-equilibrium (HWE)\n")
    ##### ii) identification of SNPs showing a significant deviation from Hardy-Weinberg-equilibrium (HWE) #####
    ##### --> only control samples should be used, as deviation from HWE may also indicate selection

    if (ethnicityFilter == 'yes') {
        ##### --> display HWE SNP P-values as histogram for all SNPs and SNPs with P-value less then 0.01 #####
        system(paste("plink --noweb --bfile ", qcdir, "/", alg, " --remove ", qcdir, "/", alg, ".fail.IDs --hardy --out ", qcdir, "/", alg, sep=""))
        data <- read.table(paste(qcdir,"/",alg,".hwe",sep=""), header=T, as.is=T)
        data <- data[grepl("ALL", data$TEST) & data$P >0.000000000001,]
        hist(-log10(data$P),ylab="Number of SNPs",col="palegreen", xlab="-log(P_hwe)", main="HWE_Pvalue for All SNPs")
        data <- data[data$P <0.01,]
        hist(-log10(data$P),ylab="Number of SNPs",col="palegreen",xlab="-log(P_hwe)", main="HWE_Pvalue for SNPs with HWE_Pvalue <0.01")
        abline(v=-log10(hweTh),lty=2,col="red")
    }
    if (ethnicityFilter == 'no') {
        ancestry_filter <- lapply(ancestry[-c(4,5)], function(x) if (length(which(x %in% unique_samples_fail_wo_ancestry)) != 0) return(x[-which(x %in% unique_samples_fail_wo_ancestry)]) else return(x))
        write.ancestry_filter <- lapply(seq_along(ancestry_filter), function(x) write.table(cbind(ancestry_filter[[x]], ancestry_filter[[x]]),paste(qcdir,"/", names(ancestry_filter)[x],"_ids.txt", sep=""), col.names=FALSE, row.names=FALSE, quote=FALSE))
        hwe.filter <- lapply(seq_along(ancestry[-c(4,5)]), function(x, ancestry) {
                             system(paste("plink --noweb --bfile ", qcdir, "/", alg, " --keep ", qcdir, "/", names(ancestry)[x], "_ids.txt --hardy --out ", qcdir, "/",names(ancestry)[x], sep=""))

                             data <- read.table(paste(qcdir,"/",names(ancestry)[x],".hwe",sep=""), header=T, as.is=T)
                             data <- data[grepl("ALL", data$TEST) & data$P >0.000000000001,]
                             removeSNP <- data[which(data$P < hweTh),]
                             hist(-log10(data$P),ylab="Number of SNPs",col="palegreen", xlab="-log(P_hwe)", main=paste("HWE_Pvalue for All SNPs (", names(ancestry)[x]," ancestry)", sep=""))
                             data <- data[data$P <0.01,]
                             hist(-log10(data$P),ylab="Number of SNPs",col="palegreen",xlab="-log(P_hwe)", main=paste("HWE_Pvalue for SNPs with HWE_Pvalue <0.01\n(", names(ancestry)[x]," ancestry)", sep=""))
                             abline(v=-log10(hweTh),lty=2,col="red")
                             return(removeSNP)
                     }, ancestry=ancestry[-c(4,5)])
       hwe.filter <- do.call(rbind, hwe.filter)
       hwe.filter.snps <- paste(unique(hwe.filter$SNP), collapse=", ")
       write.table(unique(hwe.filter$SNP), file=paste(qcdir,"/hwe.filter.snps.txt", sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE)
    }

    cat("iv) removal of markers with a very low minor allele frequency\n")
    ##### iv) removal of markers with a very low minor allele frequency #####
    ##### --> display SNP minor allele frequency as histogram ##### 

    data <- read.table(paste(qcdir,"/",alg,".frq",sep=""), header=T, as.is=T)
    hist(data$MAF,col="lightblue1",ylab="Number of SNPs",breaks=50, xlab="Minor Allele Frequency", main="Minor Allele Frequency (MAF) for All SNPs")
    abline(v=mafTh,lty=2,col="red")

    #######################################################################################
    ###### Clean Data: apply analyses steps described above to dataset (using PLINK) ######
    #######################################################################################

    cat("Clean data\n")
    # convert MAC into MAF
    MAF = mafTh/(2*(nr_samples - nr_fail_samples))
    write.table(MAF, file=paste(qcdir, "/MAFbasedonMAC", mafTh, ".txt", sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE)
    # remove SNPs that failed per-marker analysis: --remove ".fail.IDs"; calculate --freq
    if (ethnicityFilter == 'yes') {
        # generate two cleaned datasets:
        if (highIBDfilter == "IBD") {
            # i) with sample filtering based on all fail-IDs  (including relatedness) and variant filter for missing genotypes, MAF and HWE 
            system(paste("plink --noweb --bfile ", qcdir, "/", alg, " --remove ", qcdir,"/", alg, ".fail.IDs --hwe ", hweTh, " --geno ", lmissTh, " --maf ", MAF ," --make-bed --out ", qcdir, "/", alg, ".clean", sep=""))
        }
        # ii) with sample filtering based on fail-IDs not including failed IBD and variant filtering for missing genotypes and HWE (keep related samples and MAF since only part of the cohort; sophisticated relatedness filtering when combining cohorts)
        system(paste("cat ",qcdir ,"/",  alg, ".fail-ancestry.IDs  ",qcdir ,"/",  alg, ".fail-het.IDs ",qcdir ,"/",  alg, ".fail-sexcheck.IDs  ",qcdir ,"/",  alg, ".fail-imiss.IDs | sort | uniq >", qcdir, "/", alg, ".noIBDfail.IDs", sep=""), wait=TRUE)
        system(paste("plink --noweb --bfile ", qcdir, "/", alg, " --remove ", qcdir,"/", alg, ".noIBDfail.IDs --hwe ", hweTh, " --geno ", lmissTh, " --make-bed --out ", qcdir, "/", alg, ".clean.related", sep=""))
    }
    if (ethnicityFilter == 'no') {
        system(paste("plink --noweb --bfile ", qcdir, "/", alg, " --remove ", qcdir,"/", alg, ".fail.IDs --exclude ",qcdir, "/hwe.filter.snps.txt --geno ", lmissTh, " --maf ", MAF ," --make-bed --out ", qcdir, "/", alg, ".clean", sep=""))
    }
    dev.off()
}


###############
### Program ###
###############
# defaults
center="sanger"
eigenStr="evec"


args <- commandArgs(asValue=TRUE)

cat(unlist(args))

myQC(alg=args$alg, qcdir=args$qcdir, sexCheckfilter=args$sexCheckfilter, maleTh=as.numeric(args$maleTh), femaleTh=as.numeric(args$femaleTh), imissTh=as.numeric(args$imissTh), hetTh=as.numeric(args$hetTh), 
     C1=as.numeric(args$C1), C2=as.numeric(args$C2), C3=as.numeric(args$C3), highIBDfilter=args$highIBDfilter, highIBDTh=as.numeric(args$highIBDTh), crypticIBDTh=as.numeric(args$crypticIBDTh), 
     lmissTh=as.numeric(args$lmissTh), hweTh=as.numeric(args$hweTh), mafTh=as.numeric(args$mafTh), omnix_file=args$omnixfile, ethnicityFilter=args$ethnicityFilter, sample=args$sample, center=args$center, eigenStr=eigenStr)
