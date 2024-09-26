\name{print.gofcens}
\alias{print.gofcens}
\title{Printing method for \code{gofcens} object.}
\description{Printing method for \code{gofcens} object.}
\usage{
\method{print}{gofcens}(x, prnt = TRUE, outp = c("list", "table"),  degs = 3, print.AIC = TRUE,
                         print.BIC = TRUE, print.infoBoot = FALSE,  ...)
}
\arguments{
  \item{x}{An object of class \code{gofcens}.}
  \item{prnt}{Logical to indicate if the values of the Kolmogovor-Smirnov, Cramér-von Mises,
              and Anderson-Darling test statistics along with the p-values should be printed.
              Default is \code{TRUE}.}
  \item{outp}{Indicator of how the output will be displayed. The possible formats are \code{list} and \code{table}.}
  \item{degs}{Integer indicating the number of decimal places of the numeric
    results of the output.}
  \item{print.AIC}{Logical to indicate if the AIC of the model should be printed. Default is \code{TRUE}}
  \item{print.BIC}{Logical to indicate if the BIC of the model should be printed. Default is \code{TRUE}}
  \item{print.infoBoot}{Logical to indicate if the number of bootstrap samples used should be printed. Default is \code{FALSE}}
  \item{\dots}{Additional arguments.}
}

\value{
If \code{prnt = TRUE}, a list or table (if \code{outp = "table"}) containing the following components:
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
The list is also returned invisibly.
}

\author{K. Langohr, M. Besalú, M. Francisco, A. Garcia, G. Gómez.}
\examples{
\dontrun{
# List output
set.seed(123)
gofcens(times = rweibull(100, 12, scale = 4), distr = "weibull",
        BS = 149)

# Table output
set.seed(123)
print(gofcens(times = rweibull(100, 12, scale = 4), distr = "weibull",
              BS = 149), outp = "table")
}
}
