\name{gofcens}
\alias{gofcens}
\alias{gofcens.default}
\alias{gofcens.formula}
\title{Kolmogorov-Smirnov, Cramér-von Mises, and Anderson-Darling statistics for
       complete and right-censored data}
\description{
  Function \code{gofcens} computes the Kolmogorov-Smirnov, Cramér-von Mises, and
  Anderson-Darling statistics ans p-values for complete and right-censored data against
  eight possible distributions using bootstrapping.
}
\usage{
\method{gofcens}{default}(times, cens = rep(1, length(times)),
        distr = c("exponential", "gumbel", "weibull", "normal",
                  "lognormal", "logistic", "loglogistic", "beta"),
        betaLimits = c(0, 1), igumb = c(10, 10), BS = 999,
        params0 = list(shape = NULL, shape2 = NULL, location = NULL,
                       scale = NULL), ...)
\method{gofcens}{formula}(formula, data,...)
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
  \item{formula}{A formula with a numeric vector as response (which assumes no censoring) or \code{Surv} object.}
  \item{data}{Data frame for variables in \code{formula}.}
  \item{\dots}{Additional arguments.}
}
\details{Fleming et al. (1980) proposed a modified Kolmogorov-Smirnov test to be
         used with right-censored data. Koziol and Green (1976) proposed a Cramér-von
         Mises statistic for randomly censored data. This function reproduces
         this test for a given survival data and a theorical distribution. In
         presence of ties, different authors provide slightly different definitions
         of the product-limit estimator, what might provide different values of
         the test statistic.

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
}

\references{
T. R. Fleming et al. \emph{Modified Kolmogorov-Smirnov test procedure with application
to arbitrarily right-censored data}. In: Biometrics 36 (1980), 607-625.

J. A. Koziol and S. B. Green. \emph{A Cramér-von Mises statistic for randomly
censored data}. In: Biometrika, 63 (3) (1976), 465-474.

A. N. Pettitt and M. A. Stephens. \emph{Modified Cramér-von Mises statistics
for censored data.} In: Biometrika, 63 (2) (1976), 291-298.
}
\author{K. Langohr, M. Besalú, M. Francisco, A. Garcia, G. Gómez.}

\section{Warning }{If the amount of data is large, the execution time of the
                   function can be elevated. The parameter \code{BS} can
                   limit the number of random censored samples generated and
                   reduce the execution time.}

\seealso{\link[stats]{ks.test} (Package \code{stats}), \code{cvm.test}
         (Package \code{goftest}), and \code{ad.test}
         (Package \code{goftest}) for complete data, and \link{KScens} for the
         Kolmogorov-Smirnov test for right-censored data, which returns the p-value.
}

\examples{
\dontrun{
# Complete data
set.seed(123)
gofcens(times = rweibull(100, 12, scale = 4), distr = "weibull", BS = 499)

# Censored data
data(colon)
set.seed(123)
colonsamp <- colon[sample(nrow(colon), 300), ]
gofcens(Surv(time, status) ~ 1, colonsamp, distr = "normal")
}
}
