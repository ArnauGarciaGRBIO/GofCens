summary.chisqcens <- function(object, outp = c("list", "table"),
                            print.AIC = TRUE, print.BIC = TRUE,
                            print.infoBoot = FALSE, ...){
  if(!inherits(object, "chisqcens")){
    stop("Use only 'chisqcens' objects")
  }

  outp <- match.arg(outp)
  if (!outp %in% c("list", "table")) {
    stop("Invalid value of outp. Use 'table' or 'list'.")
  }
  object$outp <- outp
  object$print.AIC <- print.AIC
  object$print.BIC <- print.BIC
  object$print.infoBoot <- print.infoBoot

  class(object) <- c("summary.chisqcens", class(object))
  object
}
