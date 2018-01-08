qqplot <- function(pvalues, ci = 0.95, highlight=NULL, name="", size.title=12,
                   size.text=12) {
    N  <- length(pvalues)
    df <- data.frame(
        observed = -log10(sort(pvalues)),
        expected = -log10(1:N / N),
        clower   = -log10(qbeta(ci,     1:N, N - 1:N + 1)),
        cupper   = -log10(qbeta(1 - ci, 1:N, N - 1:N + 1))
    )
    if(!is.null(highlight)) {
        df$highlight <- highlight
    }
    xlabel=expression(Expected~~-log[10](italic(p))) 
    ylabel= expression(Observed~~-log[10](italic(p)))
    p <- ggplot2::ggplot(df) 
    p <- p + ggplot2::geom_ribbon(aes(x=expected, ymin=clower, ymax=cupper), 
                            fill="gray90") +
        ggplot2::geom_segment(aes(x = 0, y = 0, xend = max(df$expected), 
                         yend = max(df$expected)),
                         color="gray10") +
        ggplot2::xlim(0, max(df$expected)) +
        ggplot2::labs(x=xlabel, y=ylabel) +
        ggplot2::theme_bw() +
        theme(axis.title=element_text(size=size.title),
              axis.text=element_text(size=size.text)
              )
    if(!is.null(highlight)) {
        p <- p + ggplot2::geom_point(aes(x=expected, y=observed, 
										color=highlight)) +
            ggplot2::scale_color_manual(values=c("#32806E","gray10"), 
										name=name)
        
    } else {
        p <- p + ggplot2::geom_point(aes(x=expected, y=observed), col="gray10") 
    }
    p
}

