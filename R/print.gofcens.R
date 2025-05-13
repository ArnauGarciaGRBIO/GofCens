print.gofcens <- function(x, ...) {
  if(!inherits(x, "gofcens")){
    stop("Use only 'gofcens' objects")
  }
  cat("Null hypothesis: the data follows a",
      x$Distribution, "distribution \n")
  if(!is.null(x$Hypothesis)){
    cat("\nParameters:\n")
    print(x$Hypothesis)
  }
  cat("\nTest statistics\n")
  print(round(x$Test, 3))
  cat("\np-values\n")
  print(round(x$pval, 3))
  cat("\n")
}
