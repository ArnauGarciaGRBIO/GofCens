print.summary.KScens <- function(x, degs = 3, ...) {
  if (x$outp == "table") {
    cat("Distribution:", x$Distribution, "\n")
    if (!is.null(x$Hypothesis)) {
      cat("\nNull hypothesis:\n")
      header1 <- c("Parameter", "Value")
      max_col_width1 <- max(nchar(header1), nchar(names(x$Hypothesis)))
      cat(sprintf("%s | %s\n", strrep("-", max_col_width1),
                  strrep("-", max_col_width1)))
      cat(sprintf("%-*s | %-*s\n", max_col_width1, header1[1],
                  max_col_width1, header1[2]))
      cat(sprintf("%s | %s\n", strrep("-", max_col_width1),
                  strrep("-", max_col_width1)))
      for (i in 1:length(x$Hypothesis)) {
        cat(sprintf("%-*s | %-*s\n", max_col_width1, names(x$Hypothesis)[i],
                    max_col_width1, unname(x$Hypothesis)[i]))
      }
      cat(sprintf("%s | %s\n", strrep("-", max_col_width1),
                  strrep("-", max_col_width1)))
    }
    cat("\nKS Test results:\n")
    header <- c("Metric", "Value")
    max_col_width <- max(nchar(header), nchar(names(x$Test)))
    cat(sprintf("%s | %s\n", strrep("-", max_col_width),
                strrep("-", max_col_width)))
    cat(sprintf("%-*s | %-*s\n", max_col_width, header[1],
                max_col_width, header[2]))
    cat(sprintf("%s | %s\n", strrep("-", max_col_width),
                strrep("-", max_col_width)))
    for (i in 1:length(x$Test)) {
      cat(sprintf("%-*s | %-*s\n", max_col_width, names(x$Test)[i],
                  max_col_width, unname(round(x$Test, degs))[i]))
    }
    cat(sprintf("%s | %s\n", strrep("-", max_col_width),
                strrep("-", max_col_width)))
    cat("\nParameter estimates:\n")
    header1 <- c("Parameter", "Value", "s.e.")
    max_col_width1 <- max(nchar(header1), nchar(names(x$Estimates)))
    cat(sprintf("%s | %s | %s\n", strrep("-", max_col_width1),
                strrep("-", max_col_width1), strrep("-", max_col_width1)))
    cat(sprintf("%-*s | %-*s | %-*s\n", max_col_width1, header1[1],
                max_col_width1, header1[2], max_col_width1, header1[3]))
    cat(sprintf("%s | %s | %s\n", strrep("-", max_col_width1),
                strrep("-", max_col_width1),
                strrep("-", max_col_width1)))
    for (i in 1:length(x$Estimates)) {
      cat(sprintf("%-*s | %-*s | %-*s\n", max_col_width1,
                  names(x$Estimates)[i], max_col_width1,
                  unname(round(x$Estimates, degs))[i],
                  max_col_width1,
                  unname(round(x$StdErrors, degs))[i]))
    }
    cat(sprintf("%s | %s | %s\n", strrep("-", max_col_width1),
                strrep("-", max_col_width1),
                strrep("-", max_col_width1)))
    cat("\n")
    if (x$print.AIC) {
      cat( "AIC:", round(x$aic, degs), "\n")
    }
    if (x$print.BIC) {
      cat( "BIC:", round(x$bic, degs), "\n")
    }
    cat("\n")
    if (x$print.infoBoot) {
      cat( "Number of bootstrap samples:", x$BS, "\n")
      cat("\n")
    }
  } else {
    cat("Distribution:", x$Distribution, "\n")
    if (!is.null(x$Hypothesis)) {
      cat("\nNull hypothesis:\n")
      print(x$Hypothesis)
    }
    cat("\nKS Test results:\n")
    print(round(x$Test, degs))
    cat("\nParameter estimates (se):\n")
    for (i in 1:length(x$Estimates)) {
      cat(names(x$Estimates)[i], strrep(" ",
                                        2*degs+7-nchar(names(x$Estimates)[i])),
          strrep(" ", 5), sep = "")
    }
    cat("\n")
    for(i in 1:length(x$Estimates)) {
      cat(unname(round(x$Estimates, degs))[i], " ",
          "(", unname(round(x$StdErrors, degs))[i], ")", strrep(" ", 5),sep = "")
    }
    cat("\n")
    cat("\n")
    if (x$print.AIC) {
      cat( "AIC:", round(x$aic, degs), "\n")
    }
    if (x$print.BIC) {
      cat( "BIC:", round(x$bic, degs), "\n")
    }
    cat("\n")
    if (x$print.infoBoot) {
      cat( "Number of bootstrap samples:", x$BS, "\n")
      cat("\n")
    }
  }
  invisible(x)
}
