\name{KScens}
\alias{KScens}
\alias{KScens.formula}
\alias{KScens.default}
\title{Kolmogorov-Smirnov test for complete and right-censored data}
\description{
  Function \code{KScens} computes the Kolmogorov-Smirnov statistic and p-value for complete
  and right-censored data against eight possible distributions using either bootsrapping or
  a modified test.
}

\usage{
\method{KScens}{default}(times, cens = rep(1, length(times)),
       distr = c("exponential", "gumbel", "weibull", "normal",
                 "lognormal", "logistic", "loglogistic", "beta"),
       betaLimits = c(0, 1), igumb = c(10, 10), BS = 999,
       params0 = list(shape = NULL, shape2 = NULL, location = NULL,
                      scale = NULL, theta = NULL),
       tol = 1e-04, boot = TRUE, start = NULL, ...)
\method{KScens}{formula}(formula, data, ...)
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
  \item{boot}{Logical to indicate if the p-value is computed using bootstrapping or using the
              the modified Kolmogorov-Smirnov test (see details). Default is \code{TRUE}.}
  \item{start}{A named list giving the initial values of parameters of the named
               distribution or a function of data computing initial values and
               returning a named list. This argument may be omitted (default)
               for the eight prespecified distributions. See more details in
               \link[fitdistrplus]{mledist}.}
  \item{\dots}{Additional arguments for the \code{boot} function of the
    \pkg{boot} package.}
}

\details{By default the p-value is computed via bootstrapping methods.

         The parameter estimation is acomplished with the \code{fitdistcens}
         function of the \pkg{fitdistrplus} package.

         To avoid long computation times due to bootstrapping, an alternative
         with complete data is the function \code{\link[stats]{ks.test}} of the \pkg{stats} package.

         The precision of the survival times is important mainly in
         the data generation step of the bootstrap samples.

         If \code{boot = FALSE} a modified test is used to compute the p-value.
         Fleming et al. (1980) proposed a modified Kolmogorov-Smirnov test to use
         with right-censored data. This function reproduces this test for a
         given survival data and a theorical distribution. The approximation for
         the p-value is acceptable when it is smaller than 0.8 and excellent when
         it is smaller than 0.2. The output of the function follows the notation
         of Fleming et al. (1980).

         In presence of ties, different authors provide slightly different
         definitions of \eqn{\widehat{F}_n(t)}, with which other values of
         the test statistic might be obtained.}

\value{
\code{KScens} returns an object of \link[base]{class} \code{"KScens"}.

An object of class \code{"KScens"} is a list containing the following components:
  \item{\code{Distribution}}{Null distribution.}
  \item{\code{Hypothesis}}{Parameters under the null hypothesis (if \code{params0} is provided).}
  \item{\code{Test}}{Vector containing the value of the modified Kolmogorov-Smirnov statistic (\code{A}), the
  estimated p-value (\code{p-value}), the estimation of the image of the last recorded time (\code{F(y_m)}) and
  the last recorded time (\code{y_m}).}
  \item{\code{Estimates}}{Vector with the maximum likelihood estimates of the parameters
    of the distribution under study.}
  \item{\code{StdErrors}}{Vector containing the estimated standard errors.}
  \item{\code{aic}}{The Akaike information criterion.}
  \item{\code{bic}}{The so-called BIC or SBC (Schwarz Bayesian criterion).}
  \item{\code{BS}}{The number of bootstrap samples used. If the modified test is used, a 0 is returned.}
}

\references{
T. R. Fleming et al. \emph{Modified Kolmogorov-Smirnov test procedure with application
to arbitrarily right-censored data}. In: Biometrics 36 (1980), 607-625. URL: https://doi.org/10.2307/2556114
}

\author{K. Langohr, M. Besalú, M. Francisco, A. Garcia, G. Gómez.}

\seealso{Function \link[stats]{ks.test} (Package \pkg{stats}) for complete data and
         \link{gofcens} for statistics and p-value of Kolmogorov-Smirnov, Cramér von-Mises
         and Anderson-Darling together for right-censored data.
}

\examples{
# Complete data with bootstrapping
set.seed(123)
KScens(times = rweibull(100, 12, scale = 4), distr = "weibull", BS = 99)

# Censored data with bootstrapping
KScens(Surv(time, status) ~ 1, colon, distr = "norm", BS = 99)

# Censored data using the modified test
KScens(Surv(time, status) ~ 1, colon, distr = "norm", boot = FALSE)

data(nba)
summary(KScens(Surv(survtime, cens) ~ 1, nba, "logis", boot = FALSE), degs = 2)
KScens(Surv(survtime, cens) ~ 1, nba, "beta", betaLimits = c(0, 80),
       boot = FALSE)
}
