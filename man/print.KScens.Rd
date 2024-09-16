\name{print.KScens}
\alias{print.KScens}
\title{Printing method for \code{KScens} object.}
\description{Printing method for \code{KScens} object.}
\usage{
\method{print}{KScens}(x, prnt = TRUE, outp = c("list", "table"),  degs = 3, print.AIC = TRUE,
                         print.BIC = TRUE, ...)
}
\arguments{
  \item{x}{An object of class \code{KScens}.}
  \item{prnt}{Logical to indicate if the estimations of the Kolmogorov-Smirnov statistic and p-value should be printed. Default is \code{TRUE}.}
  \item{outp}{Indicator of how the output will be displayed. The possible formats are \code{list} and \code{table}.}
  \item{degs}{Integer indicating the number of decimal places of the numeric
    results of the output.}
  \item{print.AIC}{Logical to indicate if the AIC of the model should be printed. Default is \code{TRUE}}
  \item{print.BIC}{Logical to indicate if the BIC of the model should be printed. Default is \code{TRUE}}
  \item{\dots}{Additional arguments.}
}

\value{
If \code{prnt = TRUE}, a list or table (if \code{outp = "table"}) containing the following components:
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
The list is also returned invisibly.
}

\author{K. Langohr, M. Besalú, M. Francisco, G. Gómez, A. Garcia.}
\examples{
# List output
set.seed(123)
KScens(times = rweibull(1000, 12, scale = 4), distr = "weibull")

# Table output
set.seed(123)
print(KScens(times = rweibull(1000, 12, scale = 4), distr = "weibull"),
      outp = "table")
}
