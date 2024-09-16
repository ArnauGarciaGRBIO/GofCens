cumhazPlot.formula <- function(formula, data, ...) {
  call <- match.call
  m <- match.call(expand.dots = FALSE)
  m[[1]] <- as.name("model.frame")
  m$... <- NULL
  m <- eval(m, parent.frame())
  Y <- model.extract(m, "response")
  if (!inherits(Y, "Surv")) {
    stop("The left-hand side of the formula must be a 'Surv' object.")
  }
  times <- Y[, 1]
  cens <- Y[, 2]
  cumhazPlot(times, cens, ...)
}
