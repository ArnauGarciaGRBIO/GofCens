probPlot.default <- function(times, cens = rep(1, length(times)),
                             distr = c("exponential", "gumbel", "weibull", "normal",
                                       "lognormal", "logistic", "loglogistic", "beta"),
                             plots = c("PP", "QQ", "SP", "ER"),
                             colour = c("green4", "deepskyblue4", "yellow3", "mediumvioletred"),
                             mtitle = TRUE, ggp = FALSE, m = NULL, betaLimits = c(0, 1),
                             igumb = c(10, 10), prnt = TRUE, degs = 3,
                             params0 = list(shape = NULL, shape2 = NULL,
                                            location = NULL, scale = NULL),
                             print.AIC = TRUE, print.BIC = TRUE,
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
  if (!is.logical(ggp) || !is.logical(prnt)) {
    stop("ggp and prnt must be logicals!")
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
  survKM <- survfit(Surv(times, cens) ~ 1)
  tim <- summary(survKM)$time
  survTim <- summary(survKM)$surv
  uPointSurv <- survfit(Surv(tim, rep(1, length(tim))) ~ 1)
  uPoint <- 1 - uPointSurv$surv
  empiricF <- rbind(c(0, tim, Inf), c(0, uPoint, 1))
  uEstim <- rep(0, length(uPoint))
  alpha0 <- params0$shape
  gamma0 <- params0$shape2
  mu0 <- params0$location
  beta0 <- params0$scale
  alphaML <- gammaML <- muML <- betaML <- NULL
  alphaSE <- gammaSE <- muSE <- betaSE <- NULL
  aic <- bic <- NULL
  if (distr == "exponential") {
    paramsML <- survreg(Surv(times, cens) ~ 1, dist = "exponential")
    muu <- unname(coefficients(paramsML))
    betaML <- 1 / exp(-muu)
    betaSE <- sqrt(paramsML$var[1])*exp(muu)
    aic <- 2 - 2*paramsML$loglik[1]
    bic <- log(length(times)) - 2*paramsML$loglik[1]
    if (is.null(beta0)) {
      rateExp <- exp(-muu)
      outp <- list(Distribution = "exponential", Estimates = betaML,
                   StdErrors = betaSE, aic = aic, bic = bic)
    } else {
      rateExp <- 1 / beta0
      hypo <- c(scale = beta0)
      outp <- list(Distribution = "exponential", Parameters = hypo,
                   Estimates = betaML, StdErrors = betaSE, aic = aic, bic = bic)
    }
    theorPP <- pexp(tim, rateExp)
    theorQQ <- qexp(1 - survTim, rateExp)
  }
  if (distr == "gumbel") {
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
    if (is.null(mu0) || is.null(beta0)) {
      locGum <- muML
      scaleGum <- betaML
      outp <- list(Distribution = "Gumbel",
                   Estimates = c(location = muML, scale = betaML),
                   StdErrors = c(locationSE = muSE, scaleSE = betaSE),
                   aic = aic, bic = bic)
    } else {
      locGum <- mu0
      scaleGum <- beta0
      hypo <- c(location = mu0, scale = beta0)
      outp <- list(Distribution = "Gumbel", Parameters = hypo,
                   Estimates = c(location = muML, scale = betaML),
                   StdErrors = c(locationSE = muSE, scaleSE = betaSE),
                   aic = aic, bic = bic)
    }
    theorPP <- pgumbel(tim, locGum, scaleGum)
    theorQQ <- qgumbel(1 - survTim, locGum, scaleGum)
  }
  if (distr == "weibull") {
    paramsML <- fitdistcens(dd, "weibull")
    alphaML <- unname(paramsML$estimate[1])
    betaML <- unname(paramsML$estimate[2])
    alphaSE <- unname(paramsML$sd[1])
    betaSE <- unname(paramsML$sd[2])
    aic <- paramsML$aic
    bic <- paramsML$bic
    if (is.null(alpha0) || is.null(beta0)) {
      shapeWei <- alphaML
      scaleWei <- betaML
      outp <- list(Distribution = "Weibull",
                   Estimates = c(shape = alphaML, scale = betaML),
                   StdErrors = c(shapeSE = alphaSE, scaleSE = betaSE),
                   aic = aic, bic = bic)
    } else {
      shapeWei <- alpha0
      scaleWei <- beta0
      hypo <- c(shape = alpha0, scale = beta0)
      outp <- list(Distribution = "Weibull", Parameters = hypo,
                   Estimates = c(shape = alphaML, scale = betaML),
                   StdErrors = c(shapeSE = alphaSE, scaleSE = betaSE),
                   aic = aic, bic = bic)
    }
    theorPP <- pweibull(tim, shapeWei, scaleWei)
    theorQQ <- qweibull(1 - survTim, shapeWei, scaleWei)
  }
  if (distr == "normal") {
    paramsML <- fitdistcens(dd, "norm")
    muML <- unname(paramsML$estimate[1])
    betaML <- unname(paramsML$estimate[2])
    muSE <- unname(paramsML$sd[1])
    betaSE <- unname(paramsML$sd[2])
    aic <- paramsML$aic
    bic <- paramsML$bic
    if (is.null(mu0) || is.null(beta0)) {
      locNorm <- muML
      scaleNorm <- betaML
      outp <- list(Distribution = "normal",
                   Estimates = c(location = muML, scale = betaML),
                   StdErrors = c(locationSE = muSE, scaleSE = betaSE),
                   aic = aic, bic = bic)
    } else {
      locNorm <- mu0
      scaleNorm <- beta0
      hypo <- c(location = mu0, scale = beta0)
      outp <- list(Distribution = "normal", Parameters = hypo,
                   Estimates = c(location = muML, scale = betaML),
                   StdErrors = c(locationSE = muSE, scaleSE = betaSE),
                   aic = aic, bic = bic)
    }
    theorPP <- pnorm(tim, locNorm, scaleNorm)
    theorQQ <- qnorm(1 - survTim, locNorm, scaleNorm)
  }
  if (distr == "lognormal") {
    paramsML <- fitdistcens(dd, "lnorm")
    muML <- unname(paramsML$estimate[1])
    betaML <- unname(paramsML$estimate[2])
    muSE <- unname(paramsML$sd[1])
    betaSE <- unname(paramsML$sd[2])
    aic <- paramsML$aic
    bic <- paramsML$bic
    if (is.null(mu0) || is.null(beta0)) {
      locLnorm <- muML
      scaleLnorm <- betaML
      outp <- list(Distribution = "log-normal",
                   Estimates = c(location = muML, scale = betaML),
                   StdErrors = c(locationSE = muSE, scaleSE = betaSE),
                   aic = aic, bic = bic)
    } else {
      locLnorm <- mu0
      scaleLnorm <- beta0
      hypo <- c(location = mu0, scale = beta0)
      outp <- list(Distribution = "log-normal", Parameters = hypo,
                   Estimates = c(location = muML, scale = betaML),
                   StdErrors = c(locationSE = muSE, scaleSE = betaSE),
                   aic = aic, bic = bic)
    }
    theorPP <- plnorm(tim, locLnorm, scaleLnorm)
    theorQQ <- qlnorm(1 - survTim, locLnorm, scaleLnorm)
  }
  if (distr == "logistic") {
    paramsML <- fitdistcens(dd, "logis")
    muML <- unname(paramsML$estimate[1])
    betaML <- unname(paramsML$estimate[2])
    muSE <- unname(paramsML$sd[1])
    betaSE <- unname(paramsML$sd[2])
    aic <- paramsML$aic
    bic <- paramsML$bic
    if (is.null(mu0) || is.null(beta0)) {
      locLogis <- muML
      scaleLogis <- betaML
      outp <- list(Distribution = "logistic",
                   Estimates = c(location = muML, scale = betaML),
                   StdErrors = c(locationSE = muSE, scaleSE = betaSE),
                   aic = aic, bic = bic)
    } else {
      locLogis <- mu0
      scaleLogis <- beta0
      hypo <- c(location = mu0, scale = beta0)
      outp <- list(Distribution = "logistic", Parameters = hypo,
                   Estimates = c(location = muML, scale = betaML),
                   StdErrors = c(locationSE = muSE, scaleSE = betaSE),
                   aic = aic, bic = bic)
    }
    theorPP <- plogis(tim, locLogis, scaleLogis)
    theorQQ <- qlogis(1 - survTim, locLogis, scaleLogis)
  }
  if (distr == "loglogistic") {
    paramsML <- survreg(Surv(times, cens) ~ 1, dist = "loglogistic")
    alphaML <- 1 / exp(unname(paramsML$icoef)[2])
    betaML <- exp(unname(paramsML$icoef)[1])
    alphaSE <- sqrt(paramsML$var[4])*exp(-unname(paramsML$icoef)[2])
    betaSE <- sqrt(paramsML$var[1])*exp(unname(paramsML$icoef)[1])
    aic <- 2*2 - 2*paramsML$loglik[1]
    bic <- log(length(times))*2 - 2*paramsML$loglik[1]
    if (is.null(alpha0) || is.null(beta0)) {
      shapeLoglog <- alphaML
      scaleLoglog <- betaML
      outp <- list(Distribution = "log-logistic",
                   Estimates = c(shape = alphaML, scale = betaML),
                   StdErrors = c(shapeSE = alphaSE, scaleSE = betaSE),
                   aic = aic, bic = bic)
    } else {
      shapeLoglog <- alpha0
      scaleLoglog <- beta0
      hypo <- c(shape = alpha0, scale = beta0)
      outp <- list(Distribution = "log-logistic", Parameters = hypo,
                   Estimates = c(shape = alphaML, scale = betaML),
                   StdErrors = c(shapeSE = alphaSE, scaleSE = betaSE),
                   aic = aic, bic = bic)
    }
    theorPP <- pllogis(tim, shapeLoglog, scale = scaleLoglog)
    theorQQ <- qllogis(1 - survTim, shapeLoglog, scale = scaleLoglog)
  }
  if (distr == "beta") {
    aBeta <- betaLimits[1]
    bBeta <- betaLimits[2]
    paramsML <- fitdistcens((dd - aBeta) / (bBeta - aBeta), "beta")
    alphaML <- unname(paramsML$estimate[1])
    gammaML <- unname(paramsML$estimate[2])
    alphaSE <- unname(paramsML$sd[1])
    gammaSE <- unname(paramsML$sd[2])
    aic <- paramsML$aic
    bic <- paramsML$bic
    if (is.null(alpha0) || is.null(gamma0)) {
      shape1Beta <- alphaML
      shape2Beta <- gammaML
      outp <- list(Distribution = "beta",
                   Estimates = c(shape = alphaML, shape2 = gammaML),
                   StdErrors = c(shapeSE = alphaSE, shape2SE = gammaSE),
                   interval.domain = betaLimits,
                   aic = aic, bic = bic)
    } else {
      shape1Beta <- alpha0
      shape2Beta <- gamma0
      hypo <- c(shape = alpha0, shape2 = gamma0)
      outp <- list(Distribution = "beta", Parameters = hypo,
                   Estimates = c(shape = alphaML, shape2 = gammaML),
                   StdErrors = c(shapeSE = alphaSE, shape2SE = gammaSE),
                   interval.domain = betaLimits,
                   aic = aic, bic = bic)
    }
    theorPP <- pbeta((tim - aBeta)/(bBeta - aBeta), shape1Beta, shape2Beta)
    theorQQ <- qbeta((1 - survTim), shape1Beta, shape2Beta) * (bBeta - aBeta)
               + aBeta
  }
  output <- list(times = times, cens = cens, distr = distr, plots = plots,
                 colour = colour, mtitle = mtitle, ggp = ggp, m = m,
                 betaLimits = betaLimits, igumb = igumb, degs = degs, prnt = prnt,
                 params0 = params0, tim = tim, survTim = survTim, uPoint = uPoint,
                 uEstim = uEstim, empiricF = empiricF, theorPP = theorPP,
                 theorQQ = theorQQ, alphaML = alphaML, gammaML = gammaML,
                 muSE = muSE, betaSE = betaSE, alphaSE = alphaSE, gammaSE = gammaSE,
                 muSE = muSE, betaSE = betaSE, aic = aic, bic = bic,
                 print.AIC = print.AIC, print.BIC = print.BIC,
                 outp = outp)
  class(output) <- "probPlot"
  print(output)
  plot(output)
}
