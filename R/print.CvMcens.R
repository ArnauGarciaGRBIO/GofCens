print.CvMcens <- function(x, ...) {
  if(!inherits(x, "CvMcens")){
    stop("Use only 'CvMcens' objects")
  }
  if(x$complete) {
    str <- "complete"
  } else {
    str <- "right-censored"
  }
  cat("Testing via Cramér-von Mises test if the", str, "data follows a",
      x$Distribution, "distribution \n")
  if(!is.null(x$Hypothesis)){
    cat("\nNull hypothesis:\n")
    print(x$Hypothesis)
  }
  cat("\nCvM Test results:\n")
  print(round(x$Test, 3))
  cat("\nParameter estimates (se):\n")
  for (i in 1:length(x$Estimates)) {
    cat(names(x$Estimates)[i], strrep(" ",
                                      6+7-nchar(names(x$Estimates)[i])),
        strrep(" ", 5), sep = "")
  }
  for(i in 1:length(x$Estimates)){
    cat(unname(round(x$Estimates, 3))[i], " ",
        "(", unname(round(x$StdErrors, 3))[i], ")", strrep(" ", 5),sep = "")
  }
  cat("\n")
  cat("\n")

}
