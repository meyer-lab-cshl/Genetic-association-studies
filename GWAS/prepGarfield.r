#' Prepares GWAS results for GARFIELD analysis
#'
#' GARFIELD uses p-values and chromosomal positions of genetic association
#' studies to do greedy LD pruning and enrichement analysis of the results based
#' on genomic annotations. In order to run GARFIELD, several datasets are needed
#' which can be downloaded from \url{https://www.ebi.ac.uk/birney-srv/GARFIELD/}.
#' In addition, the study GWAS results have to
#' supplied in a specific format and saved in the same location as other
#' required files. prepGarfield takes GWAS results and saves them in
#' GARFIELD-format to the specified directory. If this directory is different
#' from the directory where the required annotation datasets are stored, a
#' symbolic link to the annotation directory can be created that will be found
#' by \code{\link[garfield]{garfield.run}}.
#'
#' @param gwas [NrSNP x 3 (+x)] data.frame  with chromosome identifier (chr_name
#' column), chromosomol base pair position (bp_name column) and association
#' p-values (p_name column) for NrSNP genetic variants.
#' @param directory [character] full/path/to/directory where GARFIELD results
#' will be stored; do not supply path with tilde expansion as this is not
#' supported on every platform.
#' @param trait_name [character] Name of analysed trait; will be used as
#' subdirectory name within directory to save files prepared for GARFIELD.
#' @param chr_name [character] Name of chromosome identifier column. Chromosomes
#' can be specified as simple numbers i.e. e.g. 1 or '1' or as 'chr1'.
#' @param bp_name [character] Name of base pair position column. Entries have to
#' be integers.
#' @param p_name [character] Name of p-value column. Entries have to be floats.
#' @param garfielddir [character] Full/path/to/garfiled/directory; Path were
#' "annotation","maftssd","pval" and "tags" subdirectories with per chromosome
#' files of input data are stored. For more information check
#' \code{\link[garfield]{garfield.run}}.
#' @return Returns NULL; results in GARFIELD-formated files written to
#' directory/trait_name/chr* and link created from directory/trait_name to
#' garfielddir/trait_name.


prepGarfield <- function(gwas, trait_name, directory, chr_name="CHR",
                        bp_name="BP", p_name="P", garfielddir=NULL) {
    if (grepl("~", directory)) {
        stop("directory contains ~: path expansion not guaranteed on
             every platform (see path.expand{base}), please provide full file
             path to the directory")
    }
    d <- paste(directory, "/", trait_name, sep="")
    if (!dir.exists(d)) dir.create(d, recursive=TRUE)

    if (!(chr_name %in% names(gwas))) {
        stop(paste("Column", chr_name, "not found!"))
   }
    if (!(bp_name %in% names(gwas))) {
        stop(paste("Column", bp_name, "not found!"))
    }
    if (!(p_name %in% names(gwas))) {
        stop(paste("Column", p_name, "not found!"))
    }

    colnames(gwas)[colnames(gwas) == chr_name] <- "CHR"
    colnames(gwas)[colnames(gwas) == bp_name] <- "BP"
    colnames(gwas)[colnames(gwas) == p_name] <- "P"

    if (!is.numeric(gwas$BP)) stop(paste(bp_name, "column should be numeric."))
    if (!is.numeric(gwas$P)) stop(paste(p_name, "column should be numeric."))

    chromosomes <- unique(gwas$CHR)
    writeGarfield <- sapply(chromosomes, function(chr) {
               perChr <- dplyr::filter(gwas, CHR == chr)
               chr <- gsub("chr", "", chr)
               write.table(dplyr::select(perChr, BP, P),
                           paste(d, "/chr", chr, sep=""),
                           sep=" ", col.names=FALSE, row.names=FALSE)
                        })
    if (!is.null(garfielddir)) {
         if(nzchar(Sys.readlink(garfielddir), keepNA=TRUE)) {
             message("Symbolic link ", garfielddir, "/", trait_name,
                     " exists already, won't be recreated. If this is not the ",
                     "correct link, delete link and rerun.")
         } else {
            system(paste("ln -s ", d, " ", garfielddir, sep=""))
         }
     }
}
