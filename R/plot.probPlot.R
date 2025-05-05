plot.probPlot <- function(x, ...) {
  index <- 0
  if ("ER" %in% x$plots){
    for (i in 1:length(x$uPoint)) {
      while (x$theorQQ[i] > x$empiricF[1, index + 1]) {
        index <- index + 1
      }
      if (index != 0) {
        x$uEstim[i] <- x$empiricF[2, index]
      }
    }
  }

  howmany <- length(x$plots)
  nCol <- length(x$colour)
  if (is.null(x$m)) {
    if (howmany == 1) {
      x$m <- matrix(1)
    } else {
      if (howmany == 2) {
        x$m <- matrix(1:2, nrow = 1)
      } else {
        if (howmany == 3) {
          x$m <- matrix(c(1, 2, 3, 0), byrow = TRUE, nrow = 2)
        } else {
          x$m <- matrix(c(1, 3, 2, 4), nrow = 2)
        }
      }
    }
  }
  if (!x$ggp) {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    layout(x$m)
    par(col = 1, las = 1, mar = c(4.5, 5, 2, 1), oma = c(0, 0, 1, 0), pch = 16,
        ...)
    if ("PP" %in% x$plots) {
      plot(1 - x$survTim, x$theorPP, col = x$colour[0 %% nCol + 1],
           xlab = "Empirical cumulative distribution",
           ylab = "Theoretical cumulative distribution",
           main = "P-P plot")
      abline(0, 1)
    }
    if ("QQ" %in% x$plots) {
      plot(x$tim, x$theorQQ, col = x$colour[1 %% nCol + 1],
           xlab = "Sample quantiles",
           ylab = "Theoretical quantiles",
           main = "Q-Q plot")
      abline(0, 1)
    }
    if ("SP" %in% x$plots) {
      plot(2 / pi * asin(sqrt(1 - x$survTim)), 2 / pi * asin(sqrt(x$theorPP)),
           col = x$colour[2 %% nCol + 1],
           xlab = expression(bold(2 / pi %*% arcsin(hat(F)[n](t)^{1/2}))),
           ylab = expression(bold(2 / pi %*% arcsin(hat(F)[0](t)^{1/2}))),
           main = "SP plot")
      abline(0, 1)
    }
    if ("ER" %in% x$plots) {
      plot(x$uPoint, x$uEstim, col = x$colour[3 %% nCol + 1],
           xlab = expression(bold(hat(F)[u](t))),
           ylab = expression(bold(hat(F)[u](paste(hat(F)[0]^{-1})(hat(F)(t))))),
           main = "ER plot")
      abline(0, 1)
    }
    if (x$mtitle) {
      title(paste("Probability plots for", x$outp$Distribution, "distribution"),
            outer = TRUE, cex.main = 1.5)
    }
  } else {
    ERx <- ERy <- PPx <- PPy <- QQx <- QQy <- SPx <- SPy <- NULL
    PP <- QQ <- SP <- ER <- NULL
    if ("PP" %in% x$plots) {
      PP <- ggplot(data = data.frame(PPx = 1 - x$survTim, PPy = x$theorPP),
                   aes(x = PPx, y = PPy)) +
                  geom_point(colour = x$colour[0 %% nCol + 1]) +
                  xlab("Empirical cumulative distribution") +
                  ylab("Theoretical cumulative distribution") +
                  geom_abline(intercept = 0) +
                  annotate("text", label = "P-P plot", x = Inf, y = -Inf, hjust = 1,
                           vjust = -1, size = 6, fontface = "bold")
    }
    if ("QQ" %in% x$plots) {
      QQ <- ggplot(data = data.frame(QQx = x$tim, QQy = x$theorQQ),
                   aes(x = QQx, y = QQy)) +
                  geom_point(colour = x$colour[1 %% nCol + 1]) +
                  xlab("Sample quantiles") +
                  ylab("Theoretical quantiles") +
                  geom_abline(intercept  =  0) +
                  annotate("text", label = "Q-Q plot", x = Inf, y = -Inf, hjust = 1,
                           vjust = -1, size = 6, fontface = "bold")
    }
    if ("SP" %in% x$plots) {
      SP <- ggplot(data = data.frame(SPx = 2 / pi * asin(sqrt(1 - x$survTim)),
                                     SPy = 2 / pi * asin(sqrt(x$theorPP))),
                   aes(x = SPx, y = SPy)) +
                  geom_point(colour = x$colour[2 %% nCol + 1]) +
                  xlab(expression(bold(2 / pi %*% arcsin(hat(F)[n](t)^{1/2})))) +
                  ylab(expression(bold(2 / pi %*% arcsin(hat(F)[0](t)^{1/2})))) +
                  geom_abline(intercept = 0) +
                  annotate("text", label = "SP plot", x = Inf, y = -Inf, hjust = 1,
                           vjust = -1, size = 6, fontface = "bold")
    }
    if ("ER" %in% x$plots) {
      ER <- ggplot(data = data.frame(ERx = x$uPoint, ERy =x$uEstim), aes(x = ERx, y = ERy)) +
        geom_point(colour = x$colour[3 %% nCol + 1]) +
        xlab(expression(hat(F)[u](t))) +
        ylab(expression(hat(F)[u](paste(hat(F)[0]^{-1})(hat(F)(t))))) +
        geom_abline(intercept = 0) +
        annotate("text", label = "ER plot", x = Inf, y = -Inf, hjust = 1,
                 vjust = -1, size = 6, fontface = "bold")
    }
    plolis <- list(PP, QQ, SP, ER)[c("PP", "QQ", "SP", "ER") %in% x$plots]
    if (x$mtitle) {
      grid.arrange(grobs = plolis, layout_matrix = x$m,
                   top = textGrob(paste("Probability plots for",
                                        x$outp$Distribution, "distribution"),
                                  gp = gpar(fontsize = 15, font = 2)))
    } else {
      grid.arrange(grobs = plolis, layout_matrix = x$m)
    }
  }
}
