\name{ADcens}
\alias{ADcens}
\alias{ADcens.default}
\alias{ADcens.formula}
\title{Anderson-Darling test for complete and right-censored data}
\description{
  Function \code{ADcens} computes the Anderson-Darling test statistic and p-value for complete
  and right-censored data against eight possible distributions using bootstrapping.
}
\usage{
\method{ADcens}{default}(times, cens = rep(1, length(times)),
       distr = c("exponential", "gumbel", "weibull", "normal",
                 "lognormal", "logistic", "loglogistic", "beta"),
       betaLimits = c(0, 1), igumb = c(10, 10), BS = 999,
       params0 = list(shape = NULL, shape2 = NULL,
                      location = NULL, scale = NULL), tol = 1e-04, ...)
\method{ADcens}{formula}(formula, data, ...)
}
\arguments{
  \item{times}{Numeric vector of times until the event of interest.}
  \item{cens}{Status indicator (1, exact time; 0, right-censored time). If not
    provided, all times are assumed to be exact.}
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
\details{
  The parameter estimation is acomplished with the \code{fitdistcens}
  function of the \pkg{fitdistrplus} package.

  To avoid long computation times due to bootstrapping, an alternative
  with complete data is the function \code{ad.test} of the \pkg{goftest} package.

  The precision of the survival times is important mainly in
  the data generation step of the bootstrap samples.
}

\value{
\code{ADcens} returns an object of \link[base]{class} \code{"ADcens"}.

An object of class \code{"ADcens"} is a list containing the following components:
  \item{\code{Distribution}}{Null distribution.}
  \item{\code{Hypothesis}}{Parameters under the null hypothesis (if \code{params0} is provided).}
  \item{\code{Test}}{Vector containing the value of the Anderson-Darling statistic (\code{AD}) and the
  estimated p-value (\code{p-value}).}
  \item{\code{Estimates}}{Vector with the maximum likelihood estimates of the parameters
    of the distribution under study.}
  \item{\code{StdErrors}}{Vector containing the estimated standard errors.}
  \item{\code{aic}}{The Akaike information criterion.}
  \item{\code{bic}}{The so-called BIC or SBC (Schwarz Bayesian criterion).}
}

\references{
  G. Marsaglia and J. Marsaglia. \emph{Evaluating the Aderson-Darling Distrinution}. In: Journal os Statistical Software, Articles, 9 (2) (2004), 1-5.
}

\author{K. Langohr, M. Besalú, M. Francisco, A. Garcia, G. Gómez.}

\section{Warning }{If the amount of data is large, the execution time of the
  function can be elevated. The parameter \code{BS} can
  limit the number of random censored samples generated and
  reduce the execution time.
}

\seealso{Function \code{ad.test} (Package \pkg{goftest}) for complete data and
  function \link{gofcens} for statistics and p-value of th Kolmogorov-Smirnov,
  Cramér von-Mises and Anderson-Darling together for right-censored data.
}

\examples{
# Complete data
set.seed(123)
ADcens(times = rweibull(100, 12, scale = 4), distr = "weibull",
       BS = 199)

\dontrun{
# Censored data
set.seed(123)
colonsamp <- colon[sample(nrow(colon), 300), ]
ADcens(Surv(time, status) ~ 1, colonsamp, distr = "normal")
}
}
