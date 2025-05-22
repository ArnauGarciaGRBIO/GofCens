\name{Methods tests}
\alias{methods.CvMcens}
\alias{methods.ADcens}
\alias{methods.KScens}
\alias{methods.gofcens}
\alias{methods.chisqcens}
\alias{print.CvMcens}
\alias{print.ADcens}
\alias{print.KScens}
\alias{print.gofcens}
\alias{print.chisqcens}
\alias{summary.CvMcens}
\alias{summary.ADcens}
\alias{summary.KScens}
\alias{summary.gofcens}
\alias{summary.chisqcens}
\alias{print.summary.CvMcens}
\alias{print.summary.ADcens}
\alias{print.summary.KScens}
\alias{print.summary.gofcens}
\alias{print.summary.chisqcens}
\title{
Methods for \code{ADcens}, \code{chisqcens}, \code{CvMcens}, \code{gofcens}
and \code{KScens} objects.
}
\description{
Print, summary and print summary methods for \code{ADcens}, \code{chisqcens},
\code{CvMcens}, \code{gofcens} and \code{KScens} objects.
}
\usage{
\method{print}{ADcens}(x, \dots)

\method{summary}{ADcens}(object, outp = c("list", "table"),
                            print.AIC = TRUE, print.BIC = TRUE,
                            print.infoBoot = FALSE, \dots)

\method{print.summary}{ADcens}(x, degs = 3, \dots)

\method{print}{chisqcens}(x, \dots)

\method{summary}{chisqcens}(object, outp = c("list", "table"),
                            print.AIC = TRUE, print.BIC = TRUE,
                            print.infoBoot = FALSE, \dots)

\method{print.summary}{chisqcens}(x, degs = 3, \dots)

\method{print}{CvMcens}(x, \dots)

\method{summary}{CvMcens}(object, outp = c("list", "table"),
                            print.AIC = TRUE, print.BIC = TRUE,
                            print.infoBoot = FALSE, \dots)

\method{print.summary}{CvMcens}(x, degs = 3, \dots)

\method{print}{gofcens}(x, \dots)

\method{summary}{gofcens}(object, outp = c("list", "table"),
                            print.AIC = TRUE, print.BIC = TRUE,
                            print.infoBoot = FALSE, \dots)

\method{print.summary}{gofcens}(x, degs = 3, \dots)

\method{print}{KScens}(x, \dots)

\method{summary}{KScens}(object, outp = c("list", "table"),
                            print.AIC = TRUE, print.BIC = TRUE,
                            print.infoBoot = FALSE, \dots)

\method{print.summary}{KScens}(x, degs = 3, \dots)

}
\arguments{
  \item{object, x}{An object of class \code{ADcens} \code{chisqcens},
  \code{CvMcens}, \code{gofcens} or \code{KScens}.}
  \item{outp}{Indicator of how the output will be displayed. The possible formats are \code{list} and \code{table}.}
  \item{print.AIC}{Logical to indicate if the AIC of the model should be printed. Default is \code{TRUE}}
  \item{print.BIC}{Logical to indicate if the BIC of the model should be printed. Default is \code{TRUE}}
  \item{print.infoBoot}{Logical to indicate if the number of bootstrap samples used should be printed. Default is \code{FALSE}}
  \item{degs}{Integer indicating the number of decimal places of the numeric
    results of the output. By default is 3.}
  \item{\dots}{Optional arguments.}
}
\details{
\describe{
\item{\code{print()}}{Prints only essential information of the object.}
\item{\code{summary()}}{Produces result summaries, printing information in a more extensive
                        form and with different formats available.}
\item{\code{print.summary()}}{Prints the summary output.}
}
}
\value{
\describe{
\item{\code{print()}}{Basic information is returned on the screen.}
\item{\code{summary()}}{A list with the elements: \cr
\itemize{
    \item{\code{Distribution}: Null distribution.}
  \item{\code{Hypothesis}: Parameters under the null hypothesis (if \code{params0} is provided).}
  \item{\code{Test}: Vector containing the value of the Anderson-Darling statistic (\code{AD}) and the
  estimated p-value (\code{p-value}).}
  \item{\code{Estimates}: Vector with the maximum likelihood estimates of the parameters
    of the distribution under study.}
  \item{\code{StdErrors}: Vector containing the estimated standard errors.}
  \item{\code{aic}: The Akaike information criterion.}
  \item{\code{bic}: The so-called BIC or SBC (Schwarz Bayesian criterion).}
  \item{\code{BS}: The number of bootstrap samples used.}
}
}
\item{\code{print.summary()}}{\code{x}, with the invisible flag set to prevent printing. }
}
}

\author{
K. Langohr, M. Besalú, M. Francisco, A. Garcia, G. Gómez.
}

\examples{
# Generating data
set.seed(123)
survt <- round(rlnorm(300, 2, 1), 2)
censt <- round(rexp(300, 1 / 20), 2)
times <- pmin(survt, censt)
cens <- as.numeric(survt <= censt)

# Print method
set.seed(123)
CvMcens(times, cens, distr = "weibull", BS = 99)

# List output from summary method
set.seed(123)
summary(ADcens(times = rweibull(100, 12, scale = 4), distr = "weibull",
       BS = 149))

# Table output from summary method
set.seed(123)
summary(ADcens(times = rweibull(100, 12, scale = 4), distr = "weibull",
             BS = 99), outp = "table")
}

