chisqcens.default <- function(times, cens = rep(1, length(times)), M,
                              distr = c("exponential", "gumbel", "weibull", "normal",
                                        "lognormal", "logistic", "loglogistic", "beta"),
                              betaLimits=c(0, 1), igumb = c(10, 10),
                              BS = 999, params0 = list(shape = NULL, shape2 = NULL,
                                                       location = NULL, scale = NULL),
                              tol = 1e-04, ...) {
  if (!is.numeric(times)) {
    stop("Variable times must be numeric!")
  }
  if (any(times <= 0)) {
    stop("Times must be strictly positive!")
  }
  if (any(!cens %in% 0:1)) {
    stop("Status indicator must be either 0 or 1!")
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
  rnd <- -log(tol, 10)
  times <- round(pmax(times, tol), rnd)
  dd <- data.frame(left = as.vector(times), right = ifelse(cens == 1, times, NA))
  alpha0 <- params0$shape
  gamma0 <- params0$shape2
  mu0 <- params0$location
  beta0 <- params0$scale
  alphaML <- gammaML <- muML <- betaML <- NULL
  alphaSE <- gammaSE <- muSE <- betaSE <- NULL
  aic <- bic <- NULL
  censKM <- survfit(Surv(times, 1 - cens) ~ 1)
  survKM <- survfit(Surv(times, cens) ~ 1)
  n <- length(times)
  Morig <- M
  cb <- unique(quantile(survKM, probs = seq(0, 1, 1 / M))$quantile)
  if (anyNA(cb)) {
    cb <- cb[!is.na(cb)]
  }
  Mout <- length(cb) - 1
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
    expStat <- function(dat) {
      if (is.null(beta0)) {
        muu <- unname(coefficients(survreg(Surv(dat$times, dat$cens) ~ 1,
                                           dist = "exponential")))
        betahat <- 1 / exp(-muu)
      } else {
        betahat <- beta0
      }
      survKM <- survfit(Surv(dat$times, dat$cens) ~ 1)
      cb <- unique(quantile(survKM, probs = seq(0, 1, 1 / M))$quantile)
      if (anyNA(cb)) {
        cb <- cb[!is.na(cb)]
      }
      Mred <- length(cb) - 1
      Fhat <- c(seq(0, 1, 1 / M)[1:Mred], 1)
      cb[Mred + 1] <- Inf
      obsfreq <- n * diff(Fhat)
      expProb <- diff(pexp(cb[1:(Mred + 1)], 1 / betahat))
      v <- (obsfreq - n * expProb) / sqrt(n * expProb)
      tn <- t(v) %*% v
      return(tn)
    }
    expRnd <- function(dat, mle) {
      out <- dat
      n <- nrow(dat)
      unifn <- runif(n)
      survtimes <- round(pmax(rexp(n, mle),  tol), rnd)
      censtimes <- as.vector(quantile(censKM, unifn)$quantile)
      censtimes[is.na(censtimes)] <- Inf
      out$times <- pmin(survtimes, censtimes)
      out$cens <- as.numeric(survtimes < censtimes)
      out
    }
    beta <- ifelse(is.null(beta0), betaML, beta0)
    bts <- boot(data.frame(times, cens), expStat, R = BS, sim = "parametric",
                ran.gen = expRnd, mle = 1 / beta)
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
    gumbStat <- function(dat) {
      if (is.null(mu0) || is.null(beta0)) {
        dd <- data.frame(left = as.vector(dat$times),
                         right = ifelse(dat$cens == 1, dat$times, NA))
        paramsBSML <- fitdistcens(dd, "gumbel", start = list(alpha = muML,
                                                             scale = betaML))
        muhat <- unname(paramsBSML$estimate[1])
        betahat <- unname(paramsBSML$estimate[2])
      } else {
        muhat <- mu0
        betahat <- beta0
      }
      survKM <- survfit(Surv(dat$times, dat$cens) ~ 1)
      cb <- unique(quantile(survKM, probs = seq(0, 1, 1 / M))$quantile)
      if (anyNA(cb)) {
        cb <- cb[!is.na(cb)]
      }
      Mred <- length(cb) - 1
      Fhat <- c(seq(0, 1, 1 / M)[1:Mred], 1)
      cb[Mred + 1] <- Inf
      obsfreq <- n * diff(Fhat)
      gumbProb <- diff(pgumbel(cb[1:(M + 1)], muhat, betahat))
      gumbProb[1] <- gumbProb[1] + pgumbel(cb[1], muhat, betahat)
      v <- (obsfreq - n * gumbProb) / sqrt(n * gumbProb)
      tn <- t(v) %*% v
      return(tn)
    }
    gumbRnd <- function(dat, mle) {
      out <- dat
      n <- nrow(dat)
      unifn <- runif(n)
      survtimes <- round(pmax(rgumbel(n, mle[1], mle[2]), tol), rnd)
      censtimes <- as.vector(quantile(censKM, unifn)$quantile)
      censtimes[is.na(censtimes)] <- Inf
      out$times <- pmin(survtimes, censtimes)
      out$cens <- as.numeric(survtimes < censtimes)
      out
    }
    if (is.null(mu0) || is.null(beta0)) {
      mu <- muML
      beta <- betaML
    } else {
      mu <- mu0
      beta <- beta0
    }
    bts <- boot(data.frame(times, cens), gumbStat, R = BS, sim = "parametric",
                ran.gen = gumbRnd, mle = c(mu, beta))
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
    weiStat <- function(dat) {
      if (is.null(alpha0) || is.null(beta0)) {
        dd <- data.frame(left = as.vector(dat$times),
                         right = ifelse(dat$cens == 1, dat$times, NA))
        paramsBSML <- fitdistcens(dd, "weibull")
        alphahat <- unname(paramsBSML$estimate[1])
        betahat <- unname(paramsBSML$estimate[2])
      } else {
        alphahat <- alpha0
        betahat <- beta0
      }
      survKM <- survfit(Surv(dat$times, dat$cens) ~ 1)
      cb <- unique(quantile(survKM, probs = seq(0, 1, 1 / M))$quantile)
      if (anyNA(cb)) {
        cb <- cb[!is.na(cb)]
      }
      Mred <- length(cb) - 1
      Fhat <- c(seq(0, 1, 1 / M)[1:Mred], 1)
      cb[Mred + 1] <- Inf
      obsfreq <- n * diff(Fhat)
      weiProb <- diff(pweibull(cb[1:(Mred + 1)], alphahat, betahat))
      v <- (obsfreq - n * weiProb) / sqrt(n * weiProb)
      tn <- t(v) %*% v
      return(tn)
    }
    weiRnd <- function(dat, mle) {
      out <- dat
      n <- nrow(dat)
      unifn <- runif(n)
      survtimes <- round(pmax(rweibull(n, mle[1], mle[2]),  tol), rnd)
      censtimes <- as.vector(quantile(censKM, unifn)$quantile)
      censtimes[is.na(censtimes)] <- Inf
      out$times <- pmin(survtimes, censtimes)
      out$cens <- as.numeric(survtimes < censtimes)
      out
    }
    if (is.null(alpha0) || is.null(beta0)) {
      alpha <- alphaML
      beta <- betaML
    } else {
      alpha <- alpha0
      beta <- beta0
    }
    bts <- boot(data.frame(times, cens), weiStat, R = BS, sim = "parametric",
                ran.gen = weiRnd, mle = c(alpha, beta))
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
    normStat <- function(dat) {
      if (is.null(mu0) || is.null(beta0)) {
        dd <- data.frame(left = as.vector(dat$times),
                         right = ifelse(dat$cens == 1, dat$times, NA))
        paramsBSML <- fitdistcens(dd, "norm")
        muhat <- unname(paramsBSML$estimate[1])
        betahat <- unname(paramsBSML$estimate[2])
      } else {
        muhat <- mu0
        betahat <- beta0
      }
      survKM <- survfit(Surv(dat$times, dat$cens) ~ 1)
      cb <- unique(quantile(survKM, probs = seq(0, 1, 1 / M))$quantile)
      if (anyNA(cb)) {
        cb <- cb[!is.na(cb)]
      }
      Mred <- length(cb) - 1
      Fhat <- c(seq(0, 1, 1 / M)[1:Mred], 1)
      cb[Mred + 1] <- Inf
      obsfreq <- n * diff(Fhat)
      normProb <- diff(pnorm(cb[1:(Mred + 1)], muhat, betahat))
      normProb[1] <- normProb[1] + pnorm(cb[1], muhat, betahat)
      v <- (obsfreq - n * normProb) / sqrt(n * normProb)
      tn <- t(v) %*% v
      return(tn)
    }
    normRnd <- function(dat, mle) {
      out <- dat
      n <- nrow(dat)
      unifn <- runif(n)
      survtimes <- round(pmax(rnorm(n, mle[1], mle[2]),  tol), rnd)
      censtimes <- as.vector(quantile(censKM, unifn)$quantile)
      censtimes[is.na(censtimes)] <- Inf
      out$times <- pmin(survtimes, censtimes)
      out$cens <- as.numeric(survtimes < censtimes)
      out
    }
    if (is.null(mu0) || is.null(beta0)) {
      mu <- muML
      beta <- betaML
    } else {
      mu <- mu0
      beta <- beta0
    }
    bts <- boot(data.frame(times, cens), normStat, R = BS, sim = "parametric",
                ran.gen = normRnd, mle = c(mu, beta))
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
    lnormStat <- function(dat) {
      if (is.null(mu0) || is.null(beta0)) {
        dd <- data.frame(left = as.vector(dat$times),
                         right = ifelse(dat$cens == 1, dat$times, NA))
        paramsBSML <- fitdistcens(dd, "lnorm")
        muhat <- unname(paramsBSML$estimate[1])
        betahat <- unname(paramsBSML$estimate[2])
      } else {
        muhat <- mu0
        betahat <- beta0
      }
      survKM <- survfit(Surv(dat$times, dat$cens) ~ 1)
      cb <- unique(quantile(survKM, probs = seq(0, 1, 1 / M))$quantile)
      if (anyNA(cb)) {
        cb <- cb[!is.na(cb)]
      }
      Mred <- length(cb) - 1
      Fhat <- c(seq(0, 1, 1 / M)[1:Mred], 1)
      cb[Mred + 1] <- Inf
      obsfreq <- n * diff(Fhat)
      lnormProb <- diff(plnorm(cb[1:(Mred + 1)], muhat, betahat))
      v <- (obsfreq - n * lnormProb) / sqrt(n * lnormProb)
      tn <- t(v) %*% v
      return(tn)
    }
    lnormRnd <- function(dat, mle) {
      out <- dat
      n <- nrow(dat)
      unifn <- runif(n)
      survtimes <- round(pmax(rlnorm(n, mle[1], mle[2]), tol), rnd)
      censtimes <- as.vector(quantile(censKM, unifn)$quantile)
      censtimes[is.na(censtimes)] <- Inf
      out$times <- pmin(survtimes, censtimes)
      out$cens <- as.numeric(survtimes < censtimes)
      out
    }
    if (is.null(mu0) || is.null(beta0)) {
      mu <- muML
      beta <- betaML
    } else {
      mu <- mu0
      beta <- beta0
    }
    bts <- boot(data.frame(times, cens), lnormStat, R = BS, sim = "parametric",
                ran.gen = lnormRnd, mle = c(mu, beta))
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
    logiStat <- function(dat) {
      if (is.null(mu0) || is.null(beta0)) {
        dd <- data.frame(left = as.vector(dat$times),
                         right = ifelse(dat$cens == 1, dat$times, NA))
        paramsBSML <- fitdistcens(dd, "logis")
        muhat <- unname(paramsBSML$estimate[1])
        betahat <- unname(paramsBSML$estimate[2])
      } else {
        muhat <- mu0
        betahat <- beta0
      }
      survKM <- survfit(Surv(dat$times, dat$cens) ~ 1)
      cb <- unique(quantile(survKM, probs = seq(0, 1, 1 / M))$quantile)
      if (anyNA(cb)) {
        cb <- cb[!is.na(cb)]
      }
      Mred <- length(cb) - 1
      Fhat <- c(seq(0, 1, 1 / M)[1:Mred], 1)
      cb[Mred + 1] <- Inf
      obsfreq <- n * diff(Fhat)
      logiProb <- diff(plogis(cb[1:(Mred + 1)], muhat, betahat))
      logiProb[1] <- logiProb[1] + plogis(cb[1], muhat, betahat)
      v <- (obsfreq - n * logiProb) / sqrt(n * logiProb)
      tn <- t(v) %*% v
      return(tn)
    }
    logiRnd <- function(dat, mle) {
      out <- dat
      n <- nrow(dat)
      unifn <- runif(n)
      survtimes <- round(pmax(rlogis(n, mle[1], mle[2]),  tol), rnd)
      censtimes <- as.vector(quantile(censKM, unifn)$quantile)
      censtimes[is.na(censtimes)] <- Inf
      out$times <- pmin(survtimes, censtimes)
      out$cens <- as.numeric(survtimes < censtimes)
      out
    }
    if (is.null(mu0) || is.null(beta0)) {
      mu <- muML
      beta <- betaML
    } else {
      mu <- mu0
      beta <- beta0
    }
    bts <- boot(data.frame(times, cens), logiStat, R = BS, sim = "parametric",
                ran.gen = logiRnd, mle = c(mu, beta))
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
    llogiStat <- function(dat) {
      if (is.null(alpha0) || is.null(beta0)) {
        paramsBSML <- unname(survreg(Surv(dat$times, dat$cens) ~ 1,
                                     dist = "loglogistic")$icoef)
        alphahat <- 1 / exp(paramsBSML[2])
        betahat <- exp(paramsBSML[1])
      } else {
        alphahat <- alpha0
        betahat <- beta0
      }
      survKM <- survfit(Surv(dat$times, dat$cens) ~ 1)
      cb <- unique(quantile(survKM, probs = seq(0, 1, 1 / M))$quantile)
      if (anyNA(cb)) {
        cb <- cb[!is.na(cb)]
      }
      Mred <- length(cb) - 1
      Fhat <- c(seq(0, 1, 1 / M)[1:Mred], 1)
      cb[Mred + 1] <- Inf
      obsfreq <- n * diff(Fhat)
      llogiProb <- diff(pllogis(cb[1:(Mred + 1)], alphahat, scale = betahat))
      v <- (obsfreq - n * llogiProb) / sqrt(n * llogiProb)
      tn <- t(v) %*% v
      return(tn)
    }
    llogiRnd <- function(dat, mle) {
      out <- dat
      n <- nrow(dat)
      unifn <- runif(n)
      survtimes <- round(pmax(rllogis(n, mle[1], scale = mle[2]),  tol), rnd)
      censtimes <- as.vector(quantile(censKM, unifn)$quantile)
      censtimes[is.na(censtimes)] <- Inf
      out$times <- pmin(survtimes, censtimes)
      out$cens <- as.numeric(survtimes < censtimes)
      out
    }
    if (is.null(alpha0) || is.null(beta0)) {
      alpha <- alphaML
      beta <- betaML
    } else {
      alpha <- alpha0
      beta <- beta0
    }
    bts <- boot(data.frame(times, cens), llogiStat, R = BS, sim = "parametric",
                ran.gen = llogiRnd, mle = c(alpha, beta))
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
    betaStat <- function(dat) {
      if (is.null(alpha0) || is.null(gamma0)) {
        dd <- data.frame(left = as.vector(dat$times),
                         right = ifelse(dat$cens == 1, dat$times, NA))
        paramsBSML <- fitdistcens((dd - aBeta) / (bBeta - aBeta), "beta")
        alphahat <- unname(paramsBSML$estimate[1])
        gammahat <- unname(paramsBSML$estimate[2])
      } else {
        alphahat <- alpha0
        gammahat <- gamma0
      }
      survKM <- survfit(Surv(dat$times, dat$cens) ~ 1)
      cb <- unique(quantile(survKM, probs = seq(0, 1, 1 / M))$quantile)
      if (anyNA(cb)) {
        cb <- cb[!is.na(cb)]
      }
      Mred <- length(cb) - 1
      Fhat <- c(seq(0, 1, 1 / M)[1:Mred], 1)
      cb[Mred + 1] <- Inf
      obsfreq <- n * diff(Fhat)
      betaProb <- diff(pbeta(cb[1:(M + 1)] * (bBeta - aBeta) + aBeta, alphahat,
                             gammahat))
      v <- (obsfreq - n * betaProb) / sqrt(n * betaProb)
      tn <- t(v) %*% v
      return(tn)
    }
    betaRnd <- function(dat, mle) {
      out <- dat
      n <- nrow(dat)
      unifn <- runif(n)
      survtimes <- round(pmax(rbeta(n, alpha, gamma) * (bBeta - aBeta) + aBeta,
                              tol), rnd)
      censtimes <- as.vector(quantile(censKM, unifn)$quantile)
      censtimes[is.na(censtimes)] <- Inf
      out$times <- pmin(survtimes, censtimes)
      out$cens <- as.numeric(survtimes < censtimes)
      out
    }
    if (is.null(alpha0) || is.null(gamma0)) {
      alpha <- alphaML
      gamma <- gammaML
    } else {
      alpha <- alpha0
      gamma <- gamma0
    }
    bts <- boot(data.frame(times, cens), betaStat, R = BS, sim = "parametric",
                ran.gen = betaRnd, mle = c(alpha, gamma))
  }
  tn <- bts$t0
  pval <- (sum(bts$t[, 1] > bts$t0[1]) + 1) / (bts$R + 1)
  if (all(sapply(params0, is.null))) {
    output <- list(Distribution = distr,
                   Test = c("Statistic" = tn, "p-value" = pval),
                   Estimates = c(shape = alphaML, shape2 = gammaML,
                                       location = muML, scale = betaML),
                   StdErrors = c(shapeSE = alphaSE, shape2SE = gammaSE,
                                 locationSE = muSE, scaleSE = betaSE),
                   Cellnumbers = c("Original" = Morig, "Final" = Mout),
                   aic = aic, bic = bic,
                   BS = BS)
  } else {
    output <- list(Distribution = distr,
                   Hypothesis = hypo,
                   Test = c("Statistic" = tn, "p-value" = pval),
                   Estimates = c(shape = alphaML, shape2 = gammaML,
                                       location = muML, scale = betaML),
                   StdErrors = c(shapeSE = alphaSE, shape2SE = gammaSE,
                                 locationSE = muSE, scaleSE = betaSE),
                   Cellnumbers = c("Original" = Morig, "Final" = Mout),
                   aic = aic, bic = bic,
                   BS = BS)
  }
  class(output) <- "chisqcens"
  output
}
