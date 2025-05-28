\name{gofcens}
\alias{gofcens}
\alias{gofcens.formula}
\alias{gofcens.default}
\title{Kolmogorov-Smirnov, Cramér-von Mises, and Anderson-Darling tests for
       complete and right-censored data}
\description{
  Function \code{gofcens} computes the Kolmogorov-Smirnov, Cramér-von Mises, and
  Anderson-Darling statistics and p-values for right-censored data against eight
  possible predefined or user-specified distributions using bootstrapping. This
  function also accounts for complete data.
}
\usage{
\method{gofcens}{default}(times, cens = rep(1, length(times)),
        distr = c("exponential", "gumbel", "weibull", "normal",
                  "lognormal", "logistic", "loglogistic", "beta"),
        betaLimits = c(0, 1), igumb = c(10, 10), BS = 999,
        params0 = list(shape = NULL, shape2 = NULL,
                       location = NULL, scale = NULL, theta = NULL),
        tol = 1e-04, start = NULL, ...)
\method{gofcens}{formula}(formula, data, ...)
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
               and the beta (\code{"beta"}) distribution. In addition, if the
               character string used is \code{"name"}, every distribution for which the
               corresponding density (\code{dname}), distribution (\code{pname}),
               and random generator (\code{rname}) functions are defined, can be used.}
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
  \item{start}{A named list giving the initial values of parameters of the named
               distribution or a function of data computing initial values and
               returning a named list. This argument may be omitted (default)
               for the eight prespecified distributions. See more details in
               \link[fitdistrplus]{mledist}.}
  \item{\dots}{Additional arguments for the \code{boot} function of the
    \pkg{boot} package.}
}
\details{All p-values are calculated via bootstrapping methods. For the three
         hypothesis tests, the same data generated with the bootstrapping method are used.

         The precision of the survival times is important mainly in
         the data generation step of the bootstrap samples.

         When dealing with complete data, we recommend the use of functions
         \code{\link[stats]{ks.test}} of the \pkg{stats} package and
         \code{cvm.test} and \code{ad.test} of the \pkg{goftest} package.
}

\value{
\code{gofcens} returns an object of \link[base]{class} \code{"gofcens"}.

An object of class \code{"gofcens"} is a list containing the following components:
  \item{\code{Distribution}}{Null distribution.}
  \item{\code{Hypothesis}}{Parameters under the null hypothesis (if \code{params0} is provided).}
  \item{\code{Test}}{Vector containing the values of the Kolmogovor-Smirnov (\code{KS}), Cramér-von Mises (\code{CvM}),
        and Anderson-Darling (\code{AD}) test statistics and the estimated p-value (\code{p-value}).}
  \item{\code{Estimates}}{Vector with the maximum likelihood estimates of the parameters
    of the distribution under study.}
  \item{\code{StdErrors}}{Vector containing the estimated standard errors.}
  \item{\code{aic}}{The Akaike information criterion.}
  \item{\code{bic}}{The so-called BIC or SBC (Schwarz Bayesian criterion).}
  \item{\code{BS}}{The number of bootstrap samples used.}
}

\references{
J. A. Koziol and S. B. Green. \emph{A Cramér-von Mises statistic for randomly
censored data}. In: Biometrika, 63 (3) (1976), 465-474. URL: https://doi.org/10.1093/biomet/63.3.465

A. N. Pettitt and M. A. Stephens. \emph{Modified Cramér-von Mises statistics
for censored data.} In: Biometrika, 63 (2) (1976), 291-298. URL: https://doi.org/10.1093/biomet/63.2.291
}
\author{K. Langohr, M. Besalú, M. Francisco, A. Garcia, G. Gómez.}

\section{Warning }{If the amount of data is large, the execution time of the
                   function can be elevated. The parameter \code{BS} can
                   limit the number of random censored samples generated and
                   reduce the execution time.}

\seealso{\link[stats]{ks.test} (Package \pkg{stats}), \code{cvm.test}
         (Package \pkg{goftest}), and \code{ad.test}
         (Package \pkg{goftest}) for complete data, and \link{KScens} for the
         Kolmogorov-Smirnov test for right-censored data, which returns the p-value.
}

\examples{
\dontrun{
# Complete data
set.seed(123)
gofcens(times = rweibull(100, 12, scale = 4), distr = "weibull", BS = 499)
summary(gofcens(times = rweibull(100, 12, scale = 4), distr = "exponential"),
        outp = "table", print.infoBoot = TRUE)

# Censored data
data(colon)
set.seed(123)
colonsamp <- colon[sample(nrow(colon), 300), ]
gofcens(Surv(time, status) ~ 1, colonsamp, distr = "normal")
}
}
