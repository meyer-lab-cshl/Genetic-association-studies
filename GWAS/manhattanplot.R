manhattan <- function(d, chr = "CHR", bp = "BP", p = "P", snp="SNP",
					title=NULL, max.y="max", min.y="min", 
                    suggestiveline=0, genomewideline=-log10(5e-8), 
                    size.x.labels=12, size.y.labels=12, is.negLog10=FALSE,
                    mtvsst=FALSE, trial=FALSE, xscale=FALSE, highlight=NULL, 
                    color=c("#67a9cf", "#016c59"), a=0.5, 
					colorGenomewide = "gray90", colorSuggestive = "gray90",
					colorHighlight="green",
					linetypeGenomewide = 1, linetypeSuggestive = 2) {
    if (!(chr %in% names(d)))
        stop(paste("Column", chr, "not found!"))
    if (!(bp %in% names(d)))
        stop(paste("Column", bp, "not found!"))
    if (!(p %in% names(d)))
        stop(paste("Column", p, "not found!"))
    if (!is.numeric(d[[chr]]))
        stop(paste(chr, "column should be numeric. Do you have 'X', 'Y',", 
			" 'MT', etc? If so change to numbers and try again."))
    if (!is.numeric(d[[bp]]))
        stop(paste(bp, "column should be numeric."))
    if (!is.numeric(d[[p]]))
        stop(paste(p, "column should be numeric."))

    names(d) <- toupper(names(d))
	names(d)[names(d) == chr] <- "CHR"
	names(d)[names(d) == bp] <- "BP"
	names(d)[names(d) == p] <- "P"
    
    if (!is.null(d[[snp]])) {
		names(d)[names(d) == snp] <- "SNP"
	}

    if (!is.negLog10) {
        d = subset(d[order(d$CHR, d$BP), ], (P > 0 & P <= 1 & is.numeric(P)))
        message("Pvalues are converted to negative log10 pvalues")
        d$logp = -log10(d$P)
    } else {
        d = d[order(d$CHR, d$BP), ]
        message("log10(p values) are used to depict results")
        d$logp = d$P
    }    
    ylab=expression(-log[10](italic(p)))
    
    d <- d[d$CHR %in% 1:22, ]
    
	d <- na.omit(d)
	
    d$pos <- NA
	ticks <- NULL
	lastbase <- 0
	
	d <- d[order(d$CHR),]
	
	numchroms <- length(unique(d$CHR))
	if (numchroms == 1) {
		d$pos <- d$BP
	} else {
		for (i in unique(d$CHR)) {
			if (i == 1) {
				d[d$CHR==i, ]$pos <- d[d$CHR==i, ]$BP
			}	else {
				lastbase <- lastbase + max(subset(d,CHR==i-1)$BP)
				d[d$CHR==i, ]$pos <- d[d$CHR==i, ]$BP + lastbase
			}
			ticks <- c(ticks, 
					d[d$CHR==i, ]$pos[floor(length(d[d$CHR==i, ]$pos)/2)+1])
		}
		ticklim <- c(min(d$pos),max(d$pos))
	}
	
	mycols <- rep(color, max(d$CHR))
	
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
	
	if (numchroms == 1) {
		p <- ggplot2::ggplot(data=d, aes(x=pos, y=logp))
		p <- p + ggplot2::geom_point()
		p <- p + ggplot2::ylab(expression(-log[10](italic(p)))) + 
			ggplot2::xlab(paste("Chromosome", unique(d$CHR),"position"))
	} else {
		p <- ggplot2::ggplot(data=d, aes(x=pos, y=logp))
		p <- p + ggplot2::ylab(expression(-log[10](italic(p))))
		if (!xscale) {
			p  <- p + 
			ggplot2::scale_x_continuous(name="Chromosome", breaks=ticks, 
							   limits=ticklim, expand=c(0.01,0.01),
                               labels=(unique(d$CHR)))
		}
		if (xscale) {
			 p  <- p + 
			ggplot2::scale_x_continuous(name="Chromosome", breaks=ticks, 
							   expand=c(0,0), limits=ticklim, 
							   labels=(unique(d$CHR)))
		}
		p <- p + ggplot2::scale_y_continuous(limits = c(miny, maxy),
                                             expand=c(0.01,0.01))
	}
	
	if (trial) {
		p <- p + 
			ggplot2::geom_point(aes(color=factor(SETUP), 
								shape=factor(ANALYSIS)), alpha=a)
		p <- p + ggplot2::facet_grid(THRESHOLD ~ TYPE)
	} else if (mtvsst) {
		p <- p + ggplot2::geom_point(aes(color=TYPE, shape=MARKER), alpha = a)
		p <- p + ggplot2::scale_colour_manual(values=color)            
		p <- p + ggplot2::scale_shape_manual(values=c(20,15), guide=FALSE)
	} else {
		p <- p + ggplot2::geom_point(aes(color=as.factor(CHR)))
		p <- p + ggplot2::scale_colour_manual(values=mycols, guide=FALSE)
		p <- p + ggplot2::theme(legend.position = "none") 
	}
	if (!is.null(highlight)) {
		if (any(!(highlight %in% as.vector(d$SNP)))) {
			warning(paste("You're trying to highlight SNPs that don't", "
					exist in your results."))
		}	
		d.annotate=d[as.numeric(substr(d$SNP,3,100)) %in% highlight, ]
		p <- p + ggplot2::geom_point(data=d.annotate, colour=I(colorHighlight)) 
	}
	
    if (is.null(title)) {
	    p <- p + ggplot2::theme(title=title)
    }
	p <- p + ggplot2::theme_classic()
	p <- p + ggplot2::theme(
		axis.text.x=element_text(size=size.x.labels, colour="grey50"), 
		axis.text.y=element_text(size=size.y.labels, colour="grey50"), 
		axis.ticks=element_blank()
	)
	
	if (suggestiveline) {
		p <- p + ggplot2::geom_segment(x=min(d$pos), xend=max(d$pos),
                       	y=suggestiveline, yend=suggestiveline,
					   	colour=colorSuggestive,
					   	linetype=linetypeSuggestive)
	}
	if (genomewideline) {
		p <- p + ggplot2::geom_segment(x=min(d$pos), xend=max(d$pos),
					  	y=genomewideline, yend=genomewideline,
						colour=colorGenomewide,
						linetype=linetypeGenomewide)
	}
	p
}
