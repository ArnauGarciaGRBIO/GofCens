print.CvMcens <- function(x, ...) {
  if(!inherits(x, "CvMcens")){
    stop("Use only 'CvMcens' objects")
  }
  cat("Null hypothesis: the data follows a",
      x$Distribution, "distribution \n")
  if(!is.null(x$Hypothesis)){
    cat("\nParameters:\n")
    print(x$Hypothesis)
  }
  cat("\nCvM Test results:\n")
  print(round(x$Test, 3))
  cat("\n")
}
