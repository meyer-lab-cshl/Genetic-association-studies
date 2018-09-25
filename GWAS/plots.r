#' Manhattan plot for genome-wide association studies
#'
#' The results of a genome-wide association study are visualised the genomic
#' location of the tested genetic variant on the x-axis versus its  -log10(p)
#' of the association on the y-axis. Optionally, lines for suggestive and
#' genome-wide significance can be drawn.
#'
#' @param d [NrSNPs x 3 (+2)] dataframe with columns: i) mandatory columns for
#' [chr] for chromosome, [bp] for chromosome position, [p] for assocation pvalue
#' and ii) optional columns for [snp] variant identifier and [compare] for
#' analysis identifiers (if results from different analyses depicted in same
#' plot). [chr], [bp], [p] have to be numeric, [snp] is character, [compare] a
#' factor.
#' @param chr [character] column name of chromosome column.
#' @param bp [character] column name of chromosome position column.
#' @param p [character] column name of association p-value column.
#' @param snp [character] column name of (optional) variant identifier column.
#' @param compare [character] column name of (optional) compare column.
#' @param min.y [double] minimum y-value to plot.
#' @param max.y [double] maximum y-value to plot.
#' @param is.negLog [logical] indicates that [p] column is already converted to
#' -log10(p-value)
#' @param compareAnalysis [logical] should different analysis types be compared
#' and depicted in different colors? If so, provide additional factor column
#' [compare] with information about analysis group.
#' @param highlight [logical] vector of variant IDs to be highlighted on plot.
#' Variant ID's have to be present in the [snp] column.
#' colorHighlight [character] color for hightlighted variants.
#' color [character] vector of colors for manhattan plot.
#' @param genomewideline [double] y-value to draw genomewide significance line
#' at.
#' @param colorGenomewide [character] colors of genome-wide significance line.
#' @param linetypeGenomewide [integer] linetype of genome-wide significance line.
#' @param size.x.labels [integer] size of x-axis labels
#' @param size.y.labels [integer] size of y-axis labels
#' @param size.points [double] size of plotting points
#' @param raster [logical] set to use ggrastr::geom_point_rast for plotting
#' points i.e. rasterising points of plot. Recommended for large numbers of
#' values to plot; allows for saving final plot as .pdf.
#' return ggplot2 object of manhattan plot.

manhattan <- function(d, chr = "CHR", bp = "BP", p = "P", snp="SNP",
                    compare="TYPE", compareAnalysis=FALSE,
                    title=NULL, max.y="max", min.y="min", is.negLog=FALSE,
                    highlight=NULL, colorHighlight="green",
                    color=c("#67a9cf", "#016c59"), a=0.5,
                    genomewideline=-log10(5e-8), colorGenomewide="gray90",
                    linetypeGenomewide=1,
                    size.x.labels=12, size.y.labels=12, size.points=1,
                    raster=TRUE) {

    if (!(chr %in% names(d))) stop(paste("Column", chr, "not found!"))
    if (!(bp %in% names(d))) stop(paste("Column", bp, "not found!"))
    if (!(p %in% names(d))) stop(paste("Column", p, "not found!"))

    if (!is.numeric(d[[bp]])) stop(paste(bp, "column should be numeric."))
    if (!is.numeric(d[[p]])) stop(paste(p, "column should be numeric."))
    if (!is.numeric(d[[chr]])) {
        stop(paste(chr, "column should be numeric. Does your [chr] column",
                    "chromsomes in chr1 etc format? Are there 'X', 'Y',",
                    " 'MT', etc? If so, change them to numeric encoding."))
    }

    names(d)[names(d) == chr] <- "CHR"
    names(d)[names(d) == bp] <- "BP"
    names(d)[names(d) == p] <- "P"

    if (!is.null(d[[snp]])) {
        names(d)[names(d) == snp] <- "SNP"
    }
    if (!is.null(d[[compare]])) {
        names(d)[names(d) == compare] <- "TYPE"
    }

    d <- na.omit(d)

    if (!is.negLog) {
        if (any(d$P < 0 | d$P >= 1)) stop ("P-values have to be in range (0,1]")
        d  <- d[order(d$CHR, d$BP),]
        message("Pvalues are converted to negative log10 pvalues")
        d$logp <- -log10(d$P)
    } else {
        d <- d[order(d$CHR, d$BP),]
        message("log10(p values) are used to depict results")
        d$logp <- d$P
    }

    d$pos <- NA
    ticks <- NULL
    lastbase <- 0
    numchroms <- length(unique(d$CHR))

    if (numchroms == 1) {
        d$pos <- d$BP
    } else {
        for (i in unique(d$CHR)) {
            if (i == 1) {
                d[d$CHR==i, ]$pos <- d[d$CHR==i, ]$BP
            } else {
                lastbase <- lastbase + max(subset(d, CHR==i-1)$BP)
                d[d$CHR==i, ]$pos <- d[d$CHR==i, ]$BP + lastbase
            }
            ticks <- c(ticks,
                    d[d$CHR==i, ]$pos[floor(length(d[d$CHR==i, ]$pos)/2)+1])
        }
        ticklim <- c(min(d$pos),max(d$pos))
    }

    if (max.y == "max") {
        maxy <- ceiling(max(d$logp))
    } else {
        maxy <- max.y
    }
    if (min.y == "min") {
        miny <- floor(min(d$logp))
    } else {
        miny <- min.y
    }
    if (maxy < 8) {
        maxy <- 8
    }

    mycols <- rep(color, max(d$CHR))
    ylab <-expression(-log[10](italic(p)))

    if (numchroms == 1) {
        p <- ggplot2::ggplot(data=d, ggplot2::aes(x=pos, y=logp))
        if (! raster) {
            p <- p + ggplot2::geom_point(size=size.points)
        } else {
            p <- p + ggrastr::geom_point_rast(size=size.points)
        }
        p <- p + ggplot2::ylab(expression(-log[10](italic(p)))) +
            ggplot2::xlab(paste("Chromosome", unique(d$CHR),"position"))
    } else {
        p <- ggplot2::ggplot(data=d, ggplot2::aes(x=pos, y=logp))
        p <- p + ggplot2::ylab(expression(-log[10](italic(p))))
        p  <- p + ggplot2::scale_x_continuous(name="Chromosome", breaks=ticks,
                               limits=ticklim, expand=c(0.01,0.01),
                               labels=(unique(d$CHR)))
        p <- p + ggplot2::scale_y_continuous(limits = c(miny, maxy),
                                             expand=c(0.01,0.01))
    }

    if (compareAnalysis) {
        if (!raster) {
            p <- p + ggplot2::geom_point(ggplot2::aes(color=TYPE, alpha=a),
                                         size=size.points)
        } else {
            p <- p + ggrastr::geom_point_rast(ggplot2::aes(color=TYPE, alpha=a),
                                         size=size.points)
        }
        p <- p + ggplot2::scale_colour_manual(values=color)
    } else {
        if (!raster) {
            p <- p + ggplot2::geom_point(ggplot2::aes(color=as.factor(CHR)),
                                         size=size.points)
        } else {
            p <- p + ggrastr::geom_point_rast(ggplot2::aes(color=as.factor(CHR)),
                                         size=size.points)
        }
        p <- p + ggplot2::scale_colour_manual(values=mycols, guide=FALSE)
        p <- p + ggplot2::theme(legend.position="none")
    }
    if (!is.null(highlight)) {
        if (any(!(highlight %in% as.vector(d$SNP)))) {
            warning("SNPs selected for highlighting do not exist in d")
        }
        d.annotate <- d[as.numeric(substr(d$SNP,3,100)) %in% highlight, ]
        p <- p + ggplot2::geom_point(data=d.annotate, colour=I(colorHighlight),
                                         size=size.points)
    }

    if (is.null(title)) {
        p <- p + ggplot2::theme(title=title)
    }
    p <- p + ggplot2::theme_classic()
    p <- p + ggplot2::theme(
        axis.text.x=ggplot2::element_text(size=size.x.labels, colour="grey50"),
        axis.text.y=ggplot2::element_text(size=size.y.labels, colour="grey50"),
        axis.ticks=ggplot2::element_blank()
    )

    if (genomewideline) {
        p <- p + ggplot2::geom_segment(x=min(d$pos), xend=max(d$pos),
                        y=genomewideline, yend=genomewideline,
                        colour=colorGenomewide,
                        linetype=linetypeGenomewide)
    }
    p
}


#' Quantile-quantile plot for genome-wide association studies.
#'
#' The calibration of a genome-wide association study is visualised by depicting
#' the -log10 of the observed p-values versus the -log10 p-vlaues of expected
#' p-values.
#'
#' @param values vector with [NrSNPs] observed p-values [double].
#' @param ci threshold for confidence interval [double].
#' @param highlight vector of length [NrSNPs], with 0/1 indicating if a SNP
#' should be highlighted or not.
#' @param name [character] Title of highlight color scale.
#' @param size.title [int] size of plot title.
#' @param size.text [int] size of plot labels
#' @param raster [logical] set to use ggrastr::geom_point_rast for plotting
#' points i.e. rasterizing points of plot. Recommended for large numbers of
#' values to plot; allows for saving final plot as .pdf.
#' return ggplot2 object of quantile-quantile plot


qqplot <- function(pvalues, ci=0.95, is.negLog=FALSE,
                   highlight=NULL, name="", size.title=12,
                   size.text=12, raster=TRUE) {
    N  <- length(pvalues)
    if (is.negLog) {
        observed <- sort(pvalues, decreasing=TRUE)
    } else {
        observed <- -log10(sort(pvalues))
    }

    df <- data.frame(
        observed <- observed,
        expected <- -log10(1:N / N),
        clower   <- -log10(qbeta(ci,     1:N, N - 1:N + 1)),
        cupper   <- -log10(qbeta(1 - ci, 1:N, N - 1:N + 1))
    )
    if (!is.null(highlight)) {
        df$highlight <- highlight
    }
    xlabel <- expression(Expected~~-log[10](italic(p)))
    ylabel <- expression(Observed~~-log[10](italic(p)))

    p <- ggplot2::ggplot(df)
    p <- p + ggplot2::geom_ribbon(ggplot2::aes(x=expected, ymin=clower,
                                               ymax=cupper), fill="gray90") +
        ggplot2::geom_segment(ggplot2::aes(x=0, y=0, xend=max(df$expected),
                         yend=max(df$expected)), color="gray10") +
        ggplot2::xlim(0, max(df$expected)) +
        ggplot2::labs(x=xlabel, y=ylabel) +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.title=ggplot2::element_text(size=size.title),
              axis.text=ggplot2::element_text(size=size.text)
              )
    if (!is.null(highlight)) {
        if (!raster) {
            p <- p + ggplot2::geom_point(ggplot2::aes(x=expected, y=observed,
                                                      color=highlight))
        } else {
            p <- p + ggrastr::geom_point_rast(ggplot2::aes(x=expected,
                                                           y=observed,
                                                           color=highlight))
       }
       p <- p + ggplot2::scale_color_manual(values=c("#32806E","gray10"),
                                            name=name)
    } else {
        if (!raster) {
            p <- p + ggplot2::geom_point(ggplot2::aes(x=expected, y=observed),
                                         col="gray10")
        } else {
            p <- p + ggrastr::geom_point_rast(ggplot2::aes(x=expected,
                                                           y=observed),
                                              col="gray10")
       }
    }
    p
}
