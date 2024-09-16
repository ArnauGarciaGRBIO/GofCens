gofcens.default <- function(times, cens = rep(1, length(times)),
                            distr = c("exponential", "gumbel", "weibull", "normal",
                                      "lognormal", "logistic", "loglogistic", "beta"),
                            betaLimits = c(0, 1), igumb = c(10, 10),
                            BS = 999, params0 = list(shape = NULL, shape2 = NULL,
                                                     location = NULL, scale = NULL),
                            ...) {
  if (!is.numeric(times)) {
    stop("Variable times must be numeric!")
  }
  if (any(times <= 0)) {
    stop("Times must be strictly positive!")
  }
  if (any(!cens %in% 0:1)) {
    stop("Censoring status must be either 0 or 1!")
  }
  if (!is.list(params0)) {
    stop("params0 must be a list!")
  }
  distr <- match.arg(distr)
  if (distr == "beta" && any(times < betaLimits[1] | times > betaLimits[2])) {
    msg <- paste0("Times must be within limits! Try with 'betaLimits = c(",
                  pmax(0, min(times) - 1), ", ", ceiling(max(times) + 1), ")'.")
    stop(msg)
  }
  if (!all(sapply(params0, is.null))) {
    if (distr == "exponential" && is.null(params0$scale)) {
      stop("Argument 'params0' requires a value for the scale parameter.")
    }
    if (distr %in% c("weibull", "loglogistic") &&
        (is.null(params0$shape) || is.null(params0$scale))) {
      stop("Argument 'params0' requires values for the shape and scale parameters.")
    }
    if (distr %in% c("gumbel", "normal", "lognormal", "logistic") &&
        (is.null(params0$location) || is.null(params0$scale))) {
      stop("Argument 'params0' requires values for the location and scale parameters.")
    }
    if (distr == "beta" && (is.null(params0$shape) || is.null(params0$shape2))) {
      stop("Argument 'params0' requires values for both shape parameters.")
    }
  }
  dd <- data.frame(left = as.vector(times), right = ifelse(cens == 1, times, NA))
  alpha0 <- params0$shape
  gamma0 <- params0$shape2
  mu0 <- params0$location
  beta0 <- params0$scale
  alphaML <- gammaML <- muML <- betaML <- NULL
  alphaSE <- gammaSE <- muSE <- betaSE <- NULL
  aic <- bic <- NULL
  if (distr == "exponential") {
    if (!is.null(beta0)) {
      hypo <- c(scale = beta0)
    }
    paramsML <- survreg(Surv(times, cens) ~ 1, dist = "exponential")
    muu <- unname(coefficients(paramsML))
    betaML <- 1 / exp(-muu)
    betaSE <- sqrt(paramsML$var[1])*exp(muu)
    aic <- 2 - 2*paramsML$loglik[1]
    bic <- log(length(times)) - 2*paramsML$loglik[1]
  }
  if (distr == "gumbel") {
    if (!is.null(mu0) && !is.null(beta0)) {
      hypo <- c(location = mu0, scale = beta0)
    }
    paramsML <- try(suppressMessages(fitdistcens(dd, "gumbel",
                                                 start = list(alpha = igumb[1],
                                                              scale = igumb[2]))),
                    silent = TRUE)
    if (is(paramsML, "try-error")) {
      stop("Function failed to estimate the parameters.\n
            Try with other initial values.")
    }
    muML <- unname(paramsML$estimate[1])
    betaML <- unname(paramsML$estimate[2])
    muSE <- unname(paramsML$sd[1])
    betaSE <- unname(paramsML$sd[2])
    aic <- paramsML$aic
    bic <- paramsML$bic
  }
  if (distr == "weibull") {
    if (!is.null(alpha0) && !is.null(beta0)) {
      hypo <- c(shape = alpha0, scale = beta0)
    }
    paramsML <- fitdistcens(dd, "weibull")
    alphaML <- unname(paramsML$estimate[1])
    betaML <- unname(paramsML$estimate[2])
    alphaSE <- unname(paramsML$sd[1])
    betaSE <- unname(paramsML$sd[2])
    aic <- paramsML$aic
    bic <- paramsML$bic
  }
  if (distr == "normal") {
    if (!is.null(mu0) && !is.null(beta0)) {
      hypo <- c(location = mu0, scale = beta0)
    }
    paramsML <- fitdistcens(dd, "norm")
    muML <- unname(paramsML$estimate[1])
    betaML <- unname(paramsML$estimate[2])
    muSE <- unname(paramsML$sd[1])
    betaSE <- unname(paramsML$sd[2])
    aic <- paramsML$aic
    bic <- paramsML$bic
  }
  if (distr == "lognormal") {
    if (!is.null(mu0) && !is.null(beta0)) {
      hypo <- c(location = mu0, scale = beta0)
    }
    paramsML <- fitdistcens(dd, "lnorm")
    muML <- unname(paramsML$estimate[1])
    betaML <- unname(paramsML$estimate[2])
    muSE <- unname(paramsML$sd[1])
    betaSE <- unname(paramsML$sd[2])
    aic <- paramsML$aic
    bic <- paramsML$bic
  }
  if (distr == "logistic") {
    if (!is.null(mu0) && !is.null(beta0)) {
      hypo <- c(location = mu0, scale = beta0)
    }
    paramsML <- fitdistcens(dd, "logis")
    muML <- unname(paramsML$estimate[1])
    betaML <- unname(paramsML$estimate[2])
    muSE <- unname(paramsML$sd[1])
    betaSE <- unname(paramsML$sd[2])
    aic <- paramsML$aic
    bic <- paramsML$bic
  }
  if (distr == "loglogistic") {
    if (!is.null(alpha0) && !is.null(beta0)) {
      hypo <- c(shape = alpha0, scale = beta0)
    }
    paramsML <- survreg(Surv(times, cens) ~ 1, dist = "loglogistic")
    alphaML <- 1 / exp(unname(paramsML$icoef)[2])
    betaML <- exp(unname(paramsML$icoef)[1])
    alphaSE <- sqrt(paramsML$var[4])*exp(-unname(paramsML$icoef)[2])
    betaSE <- sqrt(paramsML$var[1])*exp(unname(paramsML$icoef)[1])
    aic <- 2*2 - 2*paramsML$loglik[1]
    bic <- log(length(times))*2 - 2*paramsML$loglik[1]
  }
  if (distr == "beta") {
    if (!is.null(alpha0) && !is.null(gamma0)) {
      hypo <- c(shape = alpha0, shape2 = gamma0)
    }
    aBeta <- betaLimits[1]
    bBeta <- betaLimits[2]
    paramsML <- fitdistcens((dd - aBeta) / (bBeta - aBeta), "beta")
    alphaML <- unname(paramsML$estimate[1])
    gammaML <- unname(paramsML$estimate[2])
    alphaSE <- unname(paramsML$sd[1])
    gammaSE <- unname(paramsML$sd[2])
    aic <- paramsML$aic
    bic <- paramsML$bic
  }
  KStest <- KScens(times, cens, distr, betaLimits, igumb, params0 = params0)
  KS <- as.vector(KStest$Test[1])
  KSp <- as.vector(KStest$Test[2])
  CvMtest <- CvMcens(times, cens, distr, betaLimits, igumb, BS = BS,
                     params0 = params0)
  CvM <- as.vector(CvMtest$Test[1])
  CvMp <- as.vector(CvMtest$Test[2])
  ADtest <- ADcens(times, cens, distr, betaLimits, igumb, BS = BS,
                   params0 = params0)
  AD <- as.vector(ADtest$Test[1])
  ADp <- as.vector(ADtest$Test[2])
  if (all(sapply(params0, is.null))) {
    output <- list(Distribution = distr,
                   Test = c(KS = KS, CvM = CvM, AD = AD),
                   pval = c(KS = KSp, CvM = CvMp, AD = ADp),
                   Estimates = c(shape = alphaML, shape2 = gammaML,
                                       location = muML, scale = betaML),
                   StdErrors = c(shapeSE = alphaSE, shape2SE = gammaSE,
                                 locationSE = muSE, scaleSE = betaSE),
                   aic = aic, bic = bic)
  } else {
    output <- list(Distribution = distr,
                   Hypothesis = hypo,
                   Test = c(KS = KS, CvM = CvM, AD = AD),
                   pval = c(KS = KSp, CvM = CvMp, AD = ADp),
                   Estimates = c(shape = alphaML, shape2 = gammaML,
                                       location = muML, scale = betaML),
                   StdErrors = c(shapeSE = alphaSE, shape2SE = gammaSE,
                                 locationSE = muSE, scaleSE = betaSE),
                   aic = aic, bic = bic)
  }
  class(output) <- "gofcens"
  output
}
