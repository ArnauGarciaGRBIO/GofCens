\name{cumhazPlot}
\alias{cumhazPlot}
\title{Cumulative hazard plots to check the goodness of fit of parametric models}
\description{
  Function \code{cumhazPlot} uses the cumulative hazard plot to check if a certain distribution
  is an appropiate choice for the data.
}

\usage{
\method{cumhazPlot}{default}(times, cens = rep(1, length(times)), distr = "all6", colour = 1,
           betaLimits = c(0, 1), igumb = c(10, 10), ggp = FALSE, m = NULL,
           prnt = TRUE, degs = 3, print.AIC = TRUE, print.BIC = TRUE,...)
\method{cumhazPlot}{formula}(formula, data, ...)
}

\arguments{
  \item{times}{Numeric vector of times until the event of interest.}
  \item{cens}{Status indicator (1, exact time; 0, right-censored time). If not
              provided, all times are assumed to be exact.}
  \item{distr}{A string specifying the names of the distributions to be studied.
               The possible distributions are the exponential (\code{"exponential"}),
               the Weibull (\code{"weibull"}), the Gumbel (\code{"gumbel"}),
               the normal (\code{"normal"}), the lognormal (\code{"lognormal"}),
               the logistic (\code{"logistic"}), the loglogistic (\code{"loglogistic"}),
               and the beta (\code{"beta"}) distribution. By default, \code{distr}
               is set to \code{"all6"}, which means that the cumulative hazard
               plots are drawn for the Weibull, loglogistic, lognormal, Gumbel,
               logistic, and normal distributions.}
  \item{colour}{Colour of the points. Default colour: black.}
  \item{betaLimits}{Two-components vector with the lower and upper bounds of the
                    Beta distribution. This argument is only required, if the beta
                    distribution is considered.}
  \item{igumb}{Two-components vector with the initial values for the estimation
               of the Gumbel distribution parameters.}
  \item{ggp}{Logical to use or not the \pkg{ggplot2} package to draw the plots.
             Default is \code{FALSE}.}
  \item{m}{Optional layout for the plots to be displayed.}
  \item{prnt}{Logical to indicate if the maximum likelihood estimates of the
              parameters of all distributions considered should be printed.
              Default is \code{TRUE}.}
  \item{degs}{Integer indicating the number of decimal places of the numeric
              results of the output.}
  \item{formula}{A formula with a numeric vector as response (which assumes no censoring) or \code{Surv} object.}
  \item{data}{Data frame for variables in \code{formula}.}
  \item{print.AIC}{Logical to indicate if the AIC of the model should be printed. Default is \code{TRUE}}
  \item{print.BIC}{Logical to indicate if the BIC of the model should be printed. Default is \code{TRUE}}
  \item{\dots}{Optional arguments for function \code{par}, if \code{ggplo = FALSE}.}
}

\details{The cumulative hazard plot is based on transforming the cumulative
         hazard function \eqn{\Lambda} in such a way that it becomes linear in \eqn{t}
         or \eqn{\log(t)}. This transformation is specific for each distribution.
         The function uses the data to compute the Nelson-Aalen estimator of the
         cumulative hazard function, \eqn{\widehat{\Lambda}}, and the
         maximum likelihood estimators of the parameters of the theoretical
         distribution under study. If the distribution fits the data, the plot is
         expected to be a straight line.

         The parameter estimation is acomplished with the \code{fitdistcens}
         function of the \pkg{fitdistrplus} package.}

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

\author{K. Langohr, M. Besalú, M. Francisco, A. Garcia, G. Gómez.}

\examples{
# Complete data and default distributions
set.seed(123)
x <- rlogis(1000, 50, 5)
cumhazPlot(x, lwd = 2)

# Censored data comparing three distributions
data(nba)
cumhazPlot(Surv(survtime, cens) ~ 1, nba, distr = c("expo", "normal", "gumbel"))
}
