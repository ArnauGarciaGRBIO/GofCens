\name{kmPlot}
\alias{kmPlot}
\alias{kmPlot.default}
\alias{kmPlot.formula}
\alias{kmPlot.probPlot}
\alias{kmPlot.probPlot}
\title{Plot of the Kaplen-Meier and parametric estimations}
\description{
  Function \code{kmPlot} is a function that generates a plot that combines a
  Kaplan-Meier survival curve and a parametric survival curve in the same graph.
  It is useful for comparing non-parametric survival estimates with the fitted
  parametric survival model.
}

\usage{
\method{kmPlot}{default}(times, cens = rep(1, length(times)), distr = "all6", colour = 1,
       betaLimits = c(0, 1), igumb = c(10, 10), ggp = FALSE, m = NULL,
       prnt = TRUE, degs = 3, print.AIC = TRUE, print.BIC = TRUE,...)
\method{kmPlot}{formula}(formula, data, ...)
}

\arguments{
  \item{times}{Numeric vector of times until the event of interest.}
  \item{cens}{Status indicator (1, exact time; 0, right-censored time). If not
    provided, all times are assumed to be exact.}
  \item{distr}{A string specifying the name of the distribution to be studied.
    The possible distributions are
    the Weibull (\code{"weibull"}), the Gumbel (\code{"gumbel"}),
    the normal (\code{"normal"}), the lognormal (\code{"lognormal"}),
    the logistic (\code{"logistic"}), the loglogistic (\code{"loglogistic"}),
    the exponential (\code{"exponential"})
    and the beta (\code{"beta"}) distribution. Default is \code{"all6"} and
    includes the fisrt 6 listed which are the most used distributions.}
  \item{colour}{Vector with three components indicating the colours of the
    displayed plots. The first element is for the survival curve, the second
    for the Kaplan-Meier curve, and the last one for the confidence intervals.}
  \item{betaLimits}{Two-components vector with the lower and upper bounds of the
    Beta distribution. This argument is only required, if the beta
    distribution is considered.}
  \item{igumb}{Two-components vector with the initial values for the estimation of
    the Gumbel distribution parameters.}
  \item{ggp}{Logical to use or not the \pkg{ggplot2} package to draw the plots.
             Default is \code{FALSE}.}
  \item{m}{Optional layout for the plots to be displayed.}
  \item{prnt}{Logical to indicate if the maximum likelihood estimates of the
    parameters should be printed. Default is \code{TRUE}.}
  \item{degs}{Integer indicating the number of decimal places of the numeric
              results of the output.}
  \item{formula}{A formula with a numeric vector as response (which assumes no censoring) or \code{Surv} object.}
  \item{data}{Data frame for variables in \code{formula}.}
  \item{print.AIC}{Logical to indicate if the AIC of the model should be printed. Default is \code{TRUE}}
  \item{print.BIC}{Logical to indicate if the BIC of the model should be printed. Default is \code{TRUE}}
  \item{\dots}{Optional arguments for function \code{par}, if \code{ggp = FALSE}.}
}

\details{
  The parameter estimation is acomplished with the \code{fitdistcens}
  function of the \pkg{fitdistrplus} package.
}

\value{
  If \code{prnt = TRUE}, the following output is returned:
  \item{\code{Distribution}}{Distribution under study.}
  \item{\code{Estimates}}{A list with the maximum likelihood estimates of the
                             parameters of all distributions considered.}
  \item{\code{StdErrors}}{Vector containing the estimated standard errors.}
  \item{\code{aic}}{The Akaike information criterion.}
  \item{\code{bic}}{The so-called BIC or SBC (Schwarz Bayesian criterion).}
  In addition, a list with the same contents is returned invisibly.
}

\references{
  Peterson Jr, Arthur V. \emph{Expressing the Kaplan-Meier estimator as a function of empirical subsurvival functions}.
  In: Journal of the American Statistical Association 72.360a (1977): 854-858.
}

\author{K. Langohr, M. Besalú, M. Francisco, A. Garcia, G. Gómez.}

\examples{
# Plots for complete data and default distributions
set.seed(123)
x <- rexp(1000, 0.1)
kmPlot(x)

# Plots for censored data using ggplot2
kmPlot(Surv(time, status) ~ 1, colon, distr= "lognormal", ggp = TRUE)

# Plots for censored data from three distributions
data(nba)
kmPlot(Surv(survtime, cens) ~ 1, nba, distr = c("normal", "weibull", "lognormal"),
       prnt = FALSE)
}

