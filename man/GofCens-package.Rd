\name{GofCens-package}
\alias{GofCens-package}
\alias{GofCens}
\docType{package}
\title{Goodness-of-Fit Methods for Right-Censored Data.}
\description{This package implements both graphical tools and goodness-of-fit
  tests for right-censored data. All functions can handle complete as well as
  right-censored data. It has implemented:
    \enumerate{
      \item Kolmogorov-Smirnov (\code{\link[GofCens]{KScens}}), Cramér-von Mises (\code{\link[GofCens]{CvMcens}}),
      and Anderson-Darling (\code{\link[GofCens]{ADcens}}) tests, which use the empirical distribution
      function for complete data and are extended for right-censored data.

      \item Generalized chi-squared-type (\code{\link[GofCens]{chisqcens}}) test, which is based
      on the squared differences between observed and expected counts using random
      cells with right-censored data.

      \item A series of graphical tools such as probability (\code{\link[GofCens]{probPlot}})
      or cumulative hazard (\code{\link[GofCens]{cumhazPlot}}) plots to guide the decision about
      the most suitable parametric model for the data.
    }
}

\details{
  All functions of the \pkg{GofCens} package can be used to check the goodness of fit of the following 8 distributions. The list shows the parametrizations of the
  survival functions.

  \enumerate{
    \item Exponential Distribution [Exp\eqn{(\beta)}] \deqn{S(t)=e^{-\frac{t}{\beta}}}
    \item Weibull Distribution [Wei(\eqn{\alpha,\,\beta})]  \deqn{S(t)=e^{-(\frac{t}{\beta})^\alpha}}
    \item Gumbel Distribution [Gum(\eqn{\mu,\,\beta})]  \deqn{S(t)=1 - e^{-e^{-\frac{t-\mu}{\beta}}}}
    \item Log-Logistic Distribution [LLogis(\eqn{\alpha, \beta})]  \deqn{S(t)=\frac{1}{1 + \left(\frac{t}{\beta}\right)^\alpha}}
    \item Logistic Distribution [Logis(\eqn{\mu,\beta})]  \deqn{S(t)=\frac{e^{-\frac{t -\mu}{\beta}}}{1 + e^{-\frac{t - \mu}{\beta}}}}
    \item Log-Normal Distribution [LN(\eqn{\mu,\beta})]  \deqn{S(t)=\int_{\frac{\log t - \mu}{\beta}}^\infty \!\frac{1}{\sqrt{2 \pi}}}
    \item Normal Distribution [N(\eqn{\mu,\beta})]  \deqn{S(t)=\int_t^\infty \! \frac{1}{\beta\sqrt{2\pi}}e^{-\frac{(x - \mu)^2}{2 \beta^2}} dx}
    \item 4-Param. Beta Distribution [Beta(\eqn{\alpha, \gamma, a, b})] \deqn{S(t)=1 - \frac{B_{(\alpha, \gamma, a, b)}(t)}{B(\alpha, \gamma)}}
  }

  The list of the parameters of the theoretical distribution can be set manually using the argument \code{params} of each function. In that case, the correspondence is: \eqn{\alpha} is the \code{shape} value, \eqn{\gamma} is the \code{shape2} value, \eqn{\mu} is the \code{location} value and \eqn{\beta} is the \code{scale} value.

  In addition, the functions \code{KScens}, \code{CvMcens}, \code{ADcens}, and \code{chisqcens} can be used with any other
  distribution for which the corresponding density (\code{dname}), distribution (\code{pname}), and random generator
  (\code{rname}) functions are defined.

  \tabular{ll}{
    Package: \tab GofCens\cr
    Type: \tab Package\cr
    Version: \tab 1.5\cr
    Date: \tab 2025-05-29\cr
    License: \tab GPL (>= 2)\cr}
}
\author{Klaus Langohr, Mireia Besalú, Matilde Francisco, Arnau Garcia, Guadalupe Gómez

  Maintainer: Klaus Langohr <klaus.langohr@upc.edu>}
