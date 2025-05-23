print.KScens <- function(x, ...) {
  if (!inherits(x, "KScens")) {
    stop("Use only 'KScens' objects")
  }
  cat("Null hypothesis: the data follows a",
      x$Distribution, "distribution \n")
  if (!is.null(x$Hypothesis)) {
    cat("\nParameters:\n")
    print(x$Hypothesis)
  }
  cat("\nKS Test results:\n")
  print(round(x$Test, 3))
  cat("\n")
}
