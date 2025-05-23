print.ADcens <- function(x, ...) {
  if (!inherits(x, "ADcens")) {
    stop("Use only 'ADcens' objects")
  }
  cat("Null hypothesis: the data follows a",
      x$Distribution, "distribution \n")
  if (!is.null(x$Hypothesis)) {
    cat("\nParameters:\n")
    print(x$Hypothesis)
  }
  cat("\nAD Test results:\n")
  print(round(x$Test, 3))
  cat("\n")
}
