\name{chisqcens}
\alias{chisqcens}
\alias{chisqcens.default}
\alias{chisqcens.formula}
\title{General chi-squared statistics for right-censored data.}
\description{
  Function \code{chisqcens} computes the general chi-squared test statistic for
  right-censored data introduced by Kim (1993) and the respective p-value
  using bootstrapping.
}
\usage{
\method{chisqcens}{default}(times, cens = rep(1, length(times)), M,
          distr = c("exponential", "gumbel", "weibull", "normal",
                    "lognormal", "logistic", "loglogistic", "beta"),
          betaLimits=c(0, 1), igumb = c(10, 10), BS = 999,
          params0 = list(shape = NULL, shape2 = NULL,
                         location = NULL, scale = NULL), tol = 1e-04, ...)
\method{chisqcens}{formula}(formula, data,...)
}
\arguments{
  \item{times}{Numeric vector of times until the event of interest.}
  \item{cens}{Status indicator (1, exact time; 0, right-censored time). If not
              provided, all times are assumed to be exact.}
  \item{M}{Number indicating the number of cells that will be considered.}
  \item{distr}{A string specifying the name of the distribution to be studied.
               The possible distributions are the exponential (\code{"exponential"}),
               the Weibull (\code{"weibull"}), the Gumbel (\code{"gumbel"}),
               the normal (\code{"normal"}), the lognormal (\code{"lognormal"}),
               the logistic (\code{"logistic"}), the loglogistic (\code{"loglogistic"}),
               and the beta (\code{"beta"}) distribution.}
  \item{betaLimits}{Two-components vector with the lower and upper bounds of the
                    Beta distribution. This argument is only required, if the beta
                    distribution is considered.}
  \item{igumb}{Two-components vector with the initial values for the estimation
               of the Gumbel distribution parameters.}
  \item{BS}{Number of bootstrap samples.}
  \item{params0}{List specifying the parameters of the theoretical distribution.
                 By default, parameters are set to \code{NULL} and estimated with
                 the maximum likelihood method. This argument is only considered,
                 if all parameters of the studied distribution are specified.}
  \item{tol}{Precision of survival times.}
  \item{formula}{A formula with a numeric vector as response (which assumes no censoring) or \code{Surv} object.}
  \item{data}{Data frame for variables in \code{formula}.}
  \item{\dots}{Additional arguments.}
}
\details{The function implements the test introduced by Kim (1993) and returns
         the value of the test statistic.

         The cell boundaries of the test are obtained via the quantiles, which
         are based on the Kaplan-Meier estimate of the distribution function.
         In the presence of right-censored data, it is possible that not all
         quantiles are estimated, and in this case, the value of \code{M}
         provided by the user is reduced.

         The parameter estimation is acomplished with the \code{fitdistcens}
         function of the \pkg{fitdistrplus} package.

         The precision of the survival times is important mainly in
         the data generation step of the bootstrap samples.
}
\value{
\code{chisqcens} returns an object of \link[base]{class} \code{"chisqcens"}.

An object of class \code{"chisqcens"} is a list containing the following components:
  \item{\code{Distribution}}{Null distribution.}
  \item{\code{Hypothesis}}{Parameters under the null hypothesis (if \code{params0} is provided).}
  \item{\code{Test}}{Vector containing the value of the test statistic (\code{Statistic}) and the
  estimated p-value (\code{p-value}).}
  \item{\code{Estimates}}{Vector with the maximum likelihood estimates of the parameters
    of the distribution under study.}
  \item{\code{StdErrors}}{Vector containing the estimated standard errors.}
  \item{\code{CellNumber}}{Vector with two values: the original cell number
                           introduced by the user and the final cell number used.}
  \item{\code{aic}}{The Akaike information criterion.}
  \item{\code{bic}}{The so-called BIC or SBC (Schwarz Bayesian criterion).}
  \item{\code{BS}}{The number of bootstrap samples used.}
}

\references{
  J. H. Kim. \emph{Chi-Square Goodness-of-Fit Tests for Randomly Censored Data.}
  In: The Annals of Statistics, 21 (3) (1993), 1621-1639.
}
\author{K. Langohr, M. Besalú, M. Francisco, A. Garcia, G. Gómez.}

\examples{
# Complete data
set.seed(123)
chisqcens(times = rgumbel(100, 12, scale = 4), M = 8, distr = "gumbel",
          BS = 99)

\dontrun{
# Censored data
set.seed(123)
colonsamp <- colon[sample(nrow(colon), 300), ]
chisqcens(Surv(time, status) ~ 1, colonsamp, M = 6, distr = "normal")
}
}
