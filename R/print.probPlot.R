print.probPlot <- function(x, ...) {
  if (x$prnt) {
    cat("Distribution:", x$outp$Distribution, "\n")
    cat(" ", "\n")
    if (!all(sapply(x$params0, is.null))) {
      cat("Parameters used in probability plots:\n")
      for (dist in x$distr) {
        if (dist %in% c("gumbel", "weibull", "normal", "logistic", "lognormal",
                        "loglogistic")) {
          if ("location" %in% names(x$outp$Parameters)) {
            cat("Location:", x$outp$Parameters[1], "\n")
          }
          if ("shape" %in% names(x$outp$Parameters)) {
            cat("   Shape:", x$outp$Parameters[1], "\n")
          }
          cat("   Scale:", x$outp$Parameters[2], "\n")
        } else if (dist == "exponential") {
          cat("   Scale:",  x$outp$Parameters, "\n")
        } else {
          cat("  Shape1:", x$outp$Parameters[1], "\n")
          cat("  Shape2:", x$outp$Parameters[2], "\n")
          cat("  Domain:", x$outp$interval.domain[1], "-",
              x$outp$interval.domain[2], "\n")
        }
      }
      cat(" ", "\n")
    }
    cat("Parameter estimates:\n")
    for (dist in x$distr) {
      if (dist %in% c("gumbel", "weibull", "normal", "logistic", "lognormal",
                      "loglogistic")) {
        if ("location" %in% names(x$outp$Estimates)) {
          cat("Location (se): ", round(x$outp$Estimates[1], x$degs), " ",
              "(", round(x$outp$StdErrors[1], x$degs), ")", "\n", sep = "")
        }
        if ("shape" %in% names(x$outp$Estimates)) {
          cat("   Shape (se): ", round(x$outp$Estimates[1], x$degs), " ",
              "(", round(x$outp$StdErrors[1], x$degs), ")", "\n", sep = "")
        }
        cat("   Scale (se): ", round(x$outp$Estimates[2], x$degs), " ",
            "(", round(x$outp$StdErrors[2], x$degs), ")", "\n", sep = "")
      } else if (dist == "exponential") {
        cat("   Scale (se): ",  round(x$outp$Estimates,x$degs), " ",
            "(", round(x$outp$StdErrors, x$degs), ")", "\n", sep = "")
      } else {
        cat("  Shape1 (se): ", round(x$outp$Estimates[1], x$degs), " ",
            "(", round(x$outp$StdErrors[1], x$degs), ")", "\n", sep = "")
        cat("  Shape2 (se): ", round(x$outp$Estimates[2], x$degs), " ",
            "(", round(x$outp$StdErrors[2], x$degs), ")", "\n", sep = "")
        cat("  Domain:", round(x$outp$interval.domain[1], x$degs), "-",
            round(x$outp$interval.domain[2], x$degs), "\n")
      }
      cat("\n")
      if (x$print.AIC) {
        cat( "   AIC:", round(x$outp$aic, x$degs), "\n")
      }
      if (x$print.BIC) {
        cat( "   BIC:", round(x$outp$bic, x$degs), "\n")
      }
      cat("\n")
    }
  }
  invisible(list(Distribution = x$distr, Estimates = round(x$outp$Estimates,
                                                           x$degs)))
}
