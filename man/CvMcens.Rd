\name{CvMcens}
\alias{CvMcens}
\alias{CvMcens.formula}
\alias{CvMcens.default}
\title{Cramér-von Mises test for complete and right-censored data}
\description{
  Function \code{CvMcens} computes the Cramér-von Mises statistic and p-value for complete
  and right-censored data against eight possible distributions.
}
\usage{
\method{CvMcens}{default}(times, cens = rep(1, length(times)),
        distr = c("exponential", "gumbel", "weibull", "normal",
                  "lognormal", "logistic", "loglogistic", "beta"),
        betaLimits = c(0, 1), igumb = c(10, 10), BS = 999,
        params0 = list(shape = NULL, shape2 = NULL,
                       location = NULL, scale = NULL, theta = NULL),
                       tol = 1e-04, start = NULL, ...)
\method{CvMcens}{formula}(formula, data,...)
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
               corresponiding density (\code{dname}), probability (\code{pname})
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
\details{
         Koziol and Green (1976) proposed a Cramér-von
         Mises statistic for randomly censored data. This function reproduces
         this test for a given survival data and a theorical distribution. In
         presence of ties, different authors provide slightly different definitions
         of the product-limit estimator, what might provide different values of
         the test statistic.

         The parameter estimation is acomplished with the \code{fitdistcens}
         function of the \pkg{fitdistrplus} package.

         To avoid long computation times due to bootstrapping, an alternative
         with complete data is the function \code{cvm.test} of the \pkg{goftest} package.

         The precision of the survival times is important mainly in
         the data generation step of the bootstrap samples.
}
\value{
\code{CvMcens} returns an object of \link[base]{class} \code{"CvMcens"}.

An object of class \code{"CvMcens"} is a list containing the following components:
  \item{\code{Distribution}}{Null distribution.}
  \item{\code{Hypothesis}}{Parameters under the null hypothesis (if \code{params0} is provided).}
  \item{\code{Test}}{Vector containing the value of the Cramér-von Mises statistic (\code{CvM}) and the
  estimated p-value (\code{p-value}).}
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

\seealso{Function \code{cvm.test} (Package \pkg{goftest}) for complete data and
         \link{gofcens} for statistics and p-value of Kolmogorov-Smirnov,
         Cramér von-Mises and Anderson-Darling together for right-censored data.
}

\examples{
# Complete data
set.seed(123)
CvMcens(times = rweibull(100, 12, scale = 4), distr = "weibull",
        BS = 199)
summary(CvMcens(times = rweibull(100, 12, scale = 4), distr = "normal",
        BS = 99), degs = 4, print.AIC = FALSE, print.BIC = FALSE)

\dontrun{
# Censored data
set.seed(123)
colonsamp <- colon[sample(nrow(colon), 300), ]
CvMcens(Surv(time, status) ~ 1, colonsamp, distr = "normal")
}
}
