print.cumhazPlot <- function(x, ...) {
  if (x$prnt) {
    cat("Parameter estimates\n")
    cat(" ", "\n")
    for (dist in x$distributions) {
      if (!is.null(x$params[[dist]])) {
        cat(dist, "\n", sep = "")
        if (dist %in% c("gumbel", "weibull", "normal", "logistic", "lognormal",
                        "loglogistic")) {
          if ("location" %in% names(x$params[[dist]])) {
            cat("Location (se): ", round(x$params[[dist]][1], x$degs), " " ,
                "(", round(x$se[[dist]][1], x$degs), ")", "\n", sep = "")
          }
          if ("shape" %in% names(x$params[[dist]])) {
            cat("   Shape (se): ", round(x$params[[dist]][1], x$degs), " " ,
                "(",round(x$se[[dist]][1], x$degs), ")", "\n", sep = "")
          }
          cat("   Scale (se): ", round(x$params[[dist]][2], x$degs), " " ,
              "(", round(x$se[[dist]][2], x$degs), ")", "\n", sep = "")
        } else {
          if (dist == "exponential") {
            cat("   Scale (se): ", round(x$params[[dist]], x$degs), " " ,
                "(", round(x$se[[dist]], x$degs), ")", "\n", sep = "")
          } else {
            cat("  Shape1 (se): ", round(x$params[[dist]]$parameters[1], x$degs), " " ,
                "(", round(x$se[[dist]][1], x$degs), ")", "\n", sep = "")
            cat("  Shape2 (se): ", round(x$params[[dist]]$parameters[2], x$degs), " " ,
                "(", round(x$se[[dist]][2], x$degs), ")", "\n", sep = "")
            cat("  Domain:", round(x$params[[dist]]$domain[1], x$degs), "-",
                round(x$params[[dist]]$domain[2],x$degs), "\n")
          }
        }
        if(x$print.AIC){
          cat( "   AIC:", round(x$aic[[dist]], x$degs), "\n")
        }
        if(x$print.BIC){
          cat( "   BIC:", round(x$bic[[dist]], x$degs), "\n")
        }
        cat("\n")
      }
    }
  }
  if ("beta" %in% x$distributions) {
    x$params$beta <- x$params$beta$parameters
  }
  invisible(lapply(x$params, round, x$degs))
}
