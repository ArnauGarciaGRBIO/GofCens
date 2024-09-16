plot.kmPlot <- function(x, ...) {
  if (is.null(x$m)) {
    nplots <- length(x$distributions)
    nro <- ifelse(nplots %in% 1:3, 1, ifelse(nplots %in% 4:6, 2, 3))
    x$m <- matrix(1:nplots, byrow = TRUE, nrow = nro)
  }
  if (!x$ggp) {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    layout(x$m)
    par(pch = 16, las = 1, mar = c(4, 4.5, 2, 1), font.lab = 4, ...)
    for (i in x$distributions) {
      plot(x$survKM, col = x$colour, xlab = "Time",
           ylab = expression(bolditalic(hat(S)(t))), main = x$titl[[i]])
      lines(x$sqtime, x$srvline[[i]])
    }
  } else {
    plolis <- vector("list", length(x$distributions))
    names(plolis) <- x$distributions
    for (i in x$distributions) {
      if (!is.null(x$srvline[[i]])) {
        plolis[[i]] <- local({
          i <- i
          tmpdat <- data.frame(x = x$sqtime, y = x$srvline[[i]])
          data <- data.frame(x$times, x$cens)
          names(data) <- c("times", "cens")
          p1 <- ggsurvplot(x$survKM, data = data,
                           ggtheme = theme_minimal(),
                           xlab = expression(bolditalic(Time)),
                           ylab = expression(bolditalic(hat(S)(t))),
                           censor = FALSE, legend = "none",
                           title = x$titl[[i]],
                           font.main = c(14, "bold", "black"),
                           palette = x$colour)$plot +
            geom_point(aes(tmpdat$x, tmpdat$y), size = 1, data = tmpdat) +
            geom_line(aes(tmpdat$x, tmpdat$y), data = tmpdat)

        })
      }
    }
    suppressWarnings(grid.arrange(grobs = plolis, layout_matrix = x$m))
  }
}
