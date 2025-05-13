print.chisqcens <- function(x, ...) {
  if(!inherits(x, "chisqcens")){
    stop("Use only 'chisqcens' objects")
  }
  cat("Null hypothesis: the data follows a",
      x$Distribution, "distribution \n")
  if(!is.null(x$Hypothesis)){
    cat("\nParameters:\n")
    print(x$Hypothesis)
  }
  cat("\nChi-squared Test results:\n")
  print(round(x$Test, 3))
  cat("\n")

}
