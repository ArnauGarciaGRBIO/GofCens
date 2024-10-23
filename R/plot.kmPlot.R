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
      plot(x$survKM, col = c(x$colour[2], x$colour[3], x$colour[3]), xlab = "Time",
           ylab = expression(bolditalic(hat(S)(t))), main = x$titl[[i]])
      lines(x$sqtime, x$srvline[[i]], col = x$colour[1])
    }
  } else {
    plolis <- vector("list", length(x$distributions))
    names(plolis) <- x$distributions
    for (i in x$distributions) {
      if (!is.null(x$srvline[[i]])) {
        plolis[[i]] <- local({
          i <- i
          tmpdat <- data.frame(x = x$sqtime, y = x$srvline[[i]])
          km_data <- data.frame(
            time = x$survKM$time,
            surv = x$survKM$surv,
            lower = x$survKM$lower,
            upper = x$survKM$upper
          )
          p1 <- ggplot(data = km_data, aes(x = km_data$time, y = km_data$surv)) +
            geom_step(color = x$colour[2], linewidth = 1) +
            geom_ribbon(aes(ymin = km_data$lower, ymax = km_data$upper), alpha = 0.2,
                        fill = x$colour[3]) +
            labs(
              title = x$titl[[i]],
              x = expression(bolditalic(Time)),
              y = expression(bolditalic(hat(S)(t)))
            ) +
            theme(
              plot.title = element_text(hjust = 0.5),
              axis.title = element_text(size = 8),
              axis.text = element_text(size = 6)
            ) +
            geom_point(aes(tmpdat$x, tmpdat$y), size = 0.5, data = tmpdat) +
            geom_line(aes(tmpdat$x, tmpdat$y), data = tmpdat, linewidth = 0.5,
                      color = x$colour[1])
        })
      }
    }
    suppressWarnings(grid.arrange(grobs = plolis, layout_matrix = x$m))
  }
}
