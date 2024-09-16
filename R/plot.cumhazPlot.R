plot.cumhazPlot <- function(x, ...) {
  if (is.null(x$m)) {
    nplots = length(x$distributions)
    nro <- ifelse(nplots %in% 1:3, 1, ifelse(nplots %in% 4:6, 2, 3))
    x$m <- matrix(1:nplots, byrow = TRUE, nrow = nro)
  }
  if (!x$ggp) {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    layout(x$m)
    par(pch = 16, las = 1, mar = c(4, 4.5, 2, 1), font.lab = 3, ...)
    for (i in x$distributions) {
      if (!is.null(x$xscale[[i]])) {
        plot(x$xscale[[i]], x$yscale[[i]], col = x$colour, xlab = x$xlabs[[i]],
             ylab = x$ylabs[[i]], main = x$titl[[i]])
        lines(x$xscale[[i]], x$regline[[i]])
      }
    }
  } else {
    plolis <- vector("list", length(x$distributions))
    names(plolis) <- x$distributions
    for (i in x$distributions) {
      if (!is.null(x$xscale[[i]])) {
        plolis[[i]] <- local({
          i <- i
          p1 <- ggplot(mapping = aes(x = x$xscale[[i]], y = x$yscale[[i]])) +
            geom_point(colour = x$colour) +
            geom_segment(aes(x = x$xscale[[i]][1], y = x$regline[[i]][1],
                             xend = rev(x$xscale[[i]])[1], yend = rev(x$regline[[i]])[1])) +
            labs(title = x$titl[[i]], size = 6, fontface = "bold") +
            xlab(x$xlabs[[i]]) +
            ylab(x$ylabs[[i]])
        })
      }
    }
    grid.arrange(grobs = plolis, layout_matrix = x$m)
  }
}
