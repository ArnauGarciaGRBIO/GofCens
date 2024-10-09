cumhazPlot.default <- function(times, cens = rep(1, length(times)),
                               distr = "all6", colour = 1, betaLimits = c(0, 1),
                               igumb = c(10, 10), ggp = FALSE, m = NULL,
                               prnt = TRUE, degs = 3, print.AIC = TRUE, print.BIC = TRUE, ...) {
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
  if (length(distr) == 1 && distr == "all6") {
    distributions <- c("weibull", "loglogistic", "lognormal", "gumbel",
                       "logistic", "normal")
  } else {
    distributions <- match.arg(distr, c("exponential", "weibull", "loglogistic",
                                        "lognormal", "gumbel", "logistic",
                                        "normal", "beta"), several.ok = TRUE)
  }
  if ("beta" %in% distributions && any(times < betaLimits[1] |
                                       times > betaLimits[2])) {
    warning("Beta distribution is ignored because of out-of-bounds values.
  Try with 'betaLimits = c(", pmax(0, min(times) - 1), ", ",
            ceiling(max(times) + 1), ")'.",
            immediate. = TRUE)
    cat("\n")
    distributions <- distributions[distributions != "beta"]
  }
  if (length(distributions) == 0) {
    stop("No distribution to be fitted!")
  }
  dd <- data.frame(left = as.vector(times), right = ifelse(cens == 1, times, NA))
  survNA <- survfit(Surv(times, cens) ~ 1, stype = 2, ctype = 1)
  Haz <- -log(summary(survNA)$surv)
  tim <- summary(survNA)$time
  genlis <- vector("list", length(distributions))
  names(genlis) <- distributions
  for (v in c("params", "se", "aic", "bic", "xscale", "yscale", "xlabs",
              "ylabs", "regline", "titl")) {
    assign(v, genlis)
  }
  xla = "Time"
  xlogla = "Log(Time)"
  if ("exponential" %in% distributions) {
    tryCatch({
      fitExpo <- survreg(Surv(times, cens) ~ 1, dist = "exponential")
      muu <- unname(coefficients(fitExpo))
      params$exponential <- 1 / exp(-muu)
      names(params$exponential) <- "scale"
      se$exponential <- sqrt(fitExpo$var[1])*exp(muu)
      names(se$exponential) <- "scale (se)"
      aic$exponential <- 2 - 2*fitExpo$loglik[1]
      bic$exponential <- log(length(times)) - 2*fitExpo$loglik[1]
      xscale$exponential <- tim
      yscale$exponential <- Haz
      xlabs$exponential <- xla
      ylabs$exponential <- expression(bolditalic(hat(Lambda)(Time)))
      regline$exponential <- exp(-muu) * tim
      titl$exponential <- "Exponential"
    }, error = function(e) e)
  }
  if ("gumbel" %in% distributions) {
    tryCatch({
      fitGumb <- fitdistcens(dd, "gumbel",
                             start = list(alpha = igumb[1], scale = igumb[2]))
      locGumb <- fitGumb$estimate[1]
      scaleGumb <- fitGumb$estimate[2]
      params$gumbel <- c(locGumb, scaleGumb)
      names(params$gumbel) <- c("location", "scale")
      locGumbSE <- fitGumb$sd[1]
      scaleGumbSE <- fitGumb$sd[2]
      se$gumbel <- c(locGumbSE, scaleGumbSE)
      names(se$gumbel) <- c("location (se)", "scale (se)")
      aic$gumbel <- fitGumb$aic
      bic$gumbel <- fitGumb$bic
      xscale$gumbel <- tim
      yscale$gumbel <- -log(-log(1 - exp(-Haz)))
      xlabs$gumbel <- xla
      ylabs$gumbel <- expression(bolditalic(-log(-log(1 - exp(-hat(Lambda)(Time))))))
      regline$gumbel <- (tim - locGumb) / scaleGumb
      titl$gumbel <- "Gumbel"
    }, error = function(e) e)
  }
  if ("weibull" %in% distributions) {
    tryCatch({
      fitWei <- fitdistcens(dd, "weibull")
      shapeWei <- fitWei$estimate[1]
      scaleWei <- fitWei$estimate[2]
      params$weibull <- c(shapeWei, scaleWei)
      shapeWeiSE <- fitWei$sd[1]
      scaleWeiSE <- fitWei$sd[2]
      se$weibull <- c(shapeWeiSE, scaleWeiSE)
      names(se$weibull) <- c("shape (se)", "scale (se)")
      aic$weibull <- fitWei$aic
      bic$weibull <- fitWei$bic
      xscale$weibull <- log(tim)
      yscale$weibull <- log(Haz)
      xlabs$weibull <- xlogla
      ylabs$weibull <- expression(bolditalic(log(hat(Lambda)(Time))))
      regline$weibull <- shapeWei * (-log(scaleWei) + log(tim))
      titl$weibull <- "Weibull"
    }, error = function(e) e)
  }
  if ("normal" %in% distributions) {
    tryCatch({
      fitNorm <- fitdistcens(dd, "norm")
      locNorm <- fitNorm$estimate[1]
      scaleNorm <- fitNorm$estimate[2]
      params$normal <- fitNorm$estimate
      names(params$normal) <- c("location", "scale")
      locNormSE <- fitNorm$sd[1]
      scaleNormSE <- fitNorm$sd[2]
      se$normal <- c(locNormSE, scaleNormSE)
      names(se$normal) <- c("location (se)", "scale (se)")
      aic$normal <- fitNorm$aic
      bic$normal <- fitNorm$bic
      xscale$normal <- tim
      yscale$normal <- qnorm(1 - exp(-Haz))
      xlabs$normal <- xla
      ylabs$normal <- expression(bolditalic(paste(Phi^-1,(1 - exp(-hat(Lambda)(Time))))))
      regline$normal <- (tim - locNorm) / scaleNorm
      titl$normal <- "Normal"
    }, error = function(e) e)
  }
  if ("logistic" %in% distributions) {
    tryCatch({
      fitLog <- fitdistcens(dd, "logis")
      locLogis <- fitLog$estimate[1]
      scaleLogis <- fitLog$estimate[2]
      params$logistic <- fitLog$estimate
      locLogisSE <- fitLog$sd[1]
      scaleLogisSE <- fitLog$sd[2]
      se$logistic <- fitLog$sd
      names(se$logistic) <- c("location (se)", "scale (se)")
      aic$logistic <- fitLog$aic
      bic$logistic <- fitLog$bic
      xscale$logistic <- tim
      yscale$logistic <- log(exp(Haz) - 1)
      xlabs$logistic <- xla
      ylabs$logistic <- expression(bolditalic(log(exp(hat(Lambda)(Time)) - 1)))
      regline$logistic <- (tim - locLogis) / scaleLogis
      titl$logistic <- "Logistic"
    }, error = function(e) e)
  }
  if ("lognormal" %in% distributions) {
    tryCatch({
      fitLnorm <- fitdistcens(dd, "lnorm")
      locLnorm <- fitLnorm$estimate[1]
      scaleLnorm <- fitLnorm$estimate[2]
      params$lognormal <- c(locLnorm, scaleLnorm)
      names(params$lognormal) <- c("location", "scale")
      locLnormSE <- fitLnorm$sd[1]
      scaleLnormSE <- fitLnorm$sd[2]
      se$lognormal <- c(locLnormSE, scaleLnormSE)
      names(se$lognormal) <- c("location (se)", "scale (se)")
      aic$lognormal <- fitLnorm$aic
      bic$lognormal <- fitLnorm$bic
      xscale$lognormal <- log(tim)
      yscale$lognormal <- qnorm(1 - exp(-Haz))
      xlabs$lognormal <- xlogla
      ylabs$lognormal <- expression(bolditalic(paste(Phi^-1,(1 - exp(-hat(Lambda)(Time))))))
      regline$lognormal <- (log(tim) - locLnorm)/scaleLnorm
      titl$lognormal <- "Lognormal"
    }, error = function(e) e)
  }
  if ("loglogistic" %in% distributions) {
    tryCatch({
      fitLoglog <- survreg(Surv(times, cens) ~ 1, dist = "loglogistic")
      shapeLoglogis <- 1 / exp(unname(fitLoglog$icoef)[2])
      scaleLoglogis <- exp(unname(fitLoglog$icoef)[1])
      params$loglogistic <- c(shapeLoglogis, scaleLoglogis)
      names(params$loglogistic) <- c("shape", "scale")
      shapeLoglogisSE <- sqrt(fitLoglog$var[4])*exp(-unname(fitLoglog$icoef)[2])
      scaleLoglogisSE <- sqrt(fitLoglog$var[1])*exp(unname(fitLoglog$icoef)[1])
      se$loglogistic <- c(shapeLoglogisSE, scaleLoglogisSE)
      names(se$loglogistic) <- c("shape (se)", "scale (se)")
      aic$loglogistic <- 2*2 - 2*fitLoglog$loglik[1]
      bic$loglogistic <- log(length(times))*2 - 2*fitLoglog$loglik[1]
      xscale$loglogistic <- log(tim)
      yscale$loglogistic <- log(exp(Haz) - 1)
      xlabs$loglogistic <- xlogla
      ylabs$loglogistic <- expression(bolditalic(log(exp(hat(Lambda)(Time) - 1))))
      regline$loglogistic <- shapeLoglogis * (log(tim) - log(scaleLoglogis))
      titl$loglogistic <- "Log-logistic"
    }, error = function(e) e)
  }
  if ("beta" %in% distributions) {
    tryCatch({
      aBeta <- betaLimits[1]
      bBeta <- betaLimits[2]
      fitBeta <- fitdistcens((dd - aBeta) / (bBeta - aBeta), "beta")
      shape1Beta <- fitBeta$estimate[1]
      shape2Beta <- fitBeta$estimate[2]
      params$beta <- list(parameters = c(shape1Beta, shape2Beta),
                          domain = betaLimits)
      shape1BetaSE <- fitBeta$sd[1]
      shape2BetaSE <- fitBeta$sd[2]
      se$beta <- c(shape1BetaSE, shape2BetaSE)
      names(se$beta) <- c("shape1 (se)", "shape2 (se)")
      aic$beta <- fitBeta$aic
      bic$beta <- fitBeta$bic
      xscale$beta <- tim
      yscale$beta <- qbeta(1 - exp(-Haz), shape1Beta, shape2Beta)
      xlabs$beta <- xla
      ylabs$beta <- expression(bolditalic(paste(F^-1,(1 - exp(-hat(Lambda)(Time))))))
      regline$beta <- (tim - aBeta) / (bBeta - aBeta)
      titl$beta <- "Beta"
    }, error = function(e) e)
  }
  output <- list(times = times, cens = cens, distributions = distributions,
                 params = params, se = se, aic = aic , bic = bic, xscale = xscale,
                 yscale = yscale, xlabs = xlabs, ylabs = ylabs,
                 regline = regline, titl = titl,colour = colour, ggp = ggp,
                 m = m, betaLimits = betaLimits, igumb = igumb,
                 prnt = prnt, degs = degs, print.AIC = print.AIC, print.BIC = print.BIC)
  class(output) <- "cumhazPlot"
  print(output)
  plot(output)
}
