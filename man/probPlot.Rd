\name{probPlot}
\alias{probPlot}
\alias{probPlot.default}
\alias{probPlot.formula}
\alias{print.probPlot}
\alias{plot.probPlot}
\title{Probability plots to check the goodness of fit of parametric models}
\description{
\code{probPlot} provides four types of probability plots: P-P plot, Q-Q plot, Stabilised probability plot, and Empirically Rescaled plot to check if a certain distribution is an appropiate choice for the data.
}

\usage{
\method{probPlot}{default}(times, cens = rep(1, length(times)),
         distr = c("exponential", "gumbel", "weibull", "normal",
                   "lognormal", "logistic", "loglogistic", "beta"),
         plots = c("PP", "QQ", "SP", "ER"),
         colour = c("green4", "deepskyblue4", "yellow3",
                    "mediumvioletred"), mtitle = TRUE, ggp = FALSE,
         m = NULL, betaLimits = c(0, 1), igumb = c(10, 10),
         prnt = TRUE, degs = 3,
         params0 = list(shape = NULL, shape2 = NULL,
                        location = NULL, scale = NULL), print.AIC = TRUE,
                        print.BIC = TRUE, ...)
\method{probPlot}{formula}(formula, data, ...)
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
  \item{plots}{Vector stating the plots to be displayed. Possible choices are
               the P-P plot (\code{"PP"}), the Q-Q plot (\code{"QQ"}),
               the Stabilised Probability plot (\code{"SP"}), and the
               Empirically Rescaled plot (\code{"ER"}). By default, all four
               plots are displayed.}
  \item{colour}{Vector indicating the colours of the displayed plots.
                The vector will be recycled if its length is smaller
                than the number of plots to be displayed.}
  \item{mtitle}{Logical to add or not the title "Probability plots for a \code{distr}
                distribution" to the plot. Default is \code{TRUE}.}
  \item{ggp}{Logical to use or not the \pkg{ggplot2} package to draw the plots.
               Default is \code{FALSE}.}
  \item{m}{Optional layout for the plots to be displayed.}
  \item{betaLimits}{Two-components vector with the lower and upper bounds of the
                    Beta distribution. This argument is only required, if the beta
                    distribution is considered.}
  \item{igumb}{Two-components vector with the initial values for the estimation
               of the Gumbel distribution parameters.}
  \item{prnt}{Logical to indicate if the maximum likelihood estimates of the
              parameters should be printed. Default is \code{TRUE}.}
  \item{degs}{Integer indicating the number of decimal places of the numeric
              results of the output.}
  \item{params0}{List specifying the parameters of the theoretical distribution.
                 By default, parameters are set to \code{NULL} and estimated with
                 the maximum likelihood method. This argument is only considered,
                 if all parameters of the studied distribution are specified.}
  \item{formula}{A formula with a numeric vector as response (which assumes no censoring) or \code{Surv} object.}
  \item{data}{Data frame for variables in \code{formula}.}
  \item{print.AIC}{Logical to indicate if the AIC of the model should be printed. Default is \code{TRUE}}
  \item{print.BIC}{Logical to indicate if the BIC of the model should be printed. Default is \code{TRUE}}
  \item{\dots}{Optional arguments for function \code{par}, if \code{ggp = FALSE}.}
}
\details{By default, function \code{probPlot} draws four plots: P-P plot,
         SP plot, Q-Q plot, and EP plot. Following, a description is given for
         each plot.

         The \bold{Probability-Probability plot} (P-P plot) depicts the empirical
         distribution, \eqn{\widehat{F}(t)}, which is obtained with the Kaplan-Meier
         estimator if data are right-censored, versus the theoretical cumulative
         distribution function (cdf), \eqn{\widehat{F_0}(t)}. If the data come
         from the chosen distribution, the points of the resulting graph are
         expected to lie on the identity line.

         The \bold{Stabilised Probability plot} (SP plot), proposed by Michael (1983),
         is a transformation of the P-P plot. It stabilises the variance of the
         plotted points. If \eqn{F_0 = F} and the parameters of \eqn{F_0} are known,
         \eqn{\widehat{F_0}(t)} corresponds to the cdf of a uniform order statistic,
         and the arcsin transformation stabilises its variance. If the data come
         from distribution \eqn{F_0}, the SP plot will resemble the identity line.

         The \bold{Quartile-Quartile plot} (Q-Q plot) is similar to the P-P plot,
         but it represents the sample quantiles versus the theoretical ones,
         that is, it plots \eqn{t} versus \eqn{\widehat{F}_0^{-1}(\widehat{F}(t))}.
         Hence, if \eqn{F_0} fits the data well, the resulting plot will resemble
         the identity line.

         A drawback of the Q-Q plot is that the plotted points are not evenly spread.
         Waller and Turnbull (1992) proposed the \bold{Empirically Rescaled plot}
         (EP plot), which plots \eqn{\widehat{F}_u(t)} against
         \eqn{\widehat{F}_u(\widehat{F}_0^{-1}(\widehat{F}(t)))}, where
         \eqn{\widehat{F}_u(t)} is the empirical cdf of the points corresponding
         to the uncensored observations. Again, if \eqn{\widehat{F}_0} fits the
         data well, the ER plot will resemble the identity line.

         By default, all four probability plots are drawn and the maximum
         likelihood estimates of the parameters of the chosen parametric model
         are returned. The parameter estimation is acomplished with the
         \code{fitdistcens} function of the \pkg{fitdistrplus} package.
}

\value{
  If \code{prnt = TRUE}, the following output is returned:
   \item{\code{Distribution}}{Distribution under study.}
   \item{\code{Parameters}}{Parameters used to draw the plots (if \code{params0}
                            is provided).}
   \item{\code{Estimates}}{A list with the maximum likelihood estimates of the
                             parameters of all distributions considered.}
  \item{\code{StdErrors}}{Vector containing the estimated standard errors.}
  \item{\code{aic}}{The Akaike information criterion.}
  \item{\code{bic}}{The so-called BIC or SBC (Schwarz Bayesian criterion).}
  In addition, a list with the same contents is returned invisibly.
}

\references{
  J. R. Michael. \emph{The Stabilized Probability Plot}. In: Biometrika 70 (1) (1983),
  11-17.

  L.A. Waller and B.W. Turnbull. \emph{Probability Plotting with Censored Data}. In:
  American Statistician 46 (1) (1992), 5-12.
}

\author{K. Langohr, M. Besalú, M. Francisco, A. Garcia, G. Gómez.}

\examples{
# P-P, Q-Q, SP, and EP plots for complete data
set.seed(123)
x <- rlnorm(1000, 3, 2)
probPlot(x)
probPlot(x, distr = "lognormal")

# P-P, Q-Q, SP, and EP plots for censored data using ggplot2
probPlot(Surv(time, status) ~ 1, colon, "weibull", ggp = TRUE)

# P-P, Q-Q and SP plots for censored data and lognormal distribution
data(nba)
probPlot(Surv(survtime, cens) ~ 1, nba, "lognorm", plots = c("PP", "QQ", "SP"),
         ggp = TRUE, m = matrix(1:3, nr = 1))
}

