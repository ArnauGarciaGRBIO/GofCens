\name{nba}
\alias{nba}
\docType{data}
\title{Survival times of former NBA players.}
\description{Survival times of former NBA players after their NBA career.}
\usage{data("nba")}
\format{
  A data frame with 3962 observations on the following 3 variables.
  \describe{
    \item{\code{id}}{Player ID}
    \item{\code{survtime}}{Time (in years) from end of NBA career until either
    death or July 31, 2019.}
    \item{\code{cens}}{Death indicator (1, exact survival time; 0, right-censored
    survival time).}
  }
}
\details{
The survival times of former NBA players were analyzed by Martínez et al. (2022).
}
\source{
J. A. Martínez, K. Langohr, J. Felipo, L. Consuegra and M. Casals. \emph{Data set on
mortality of national basketball association (NBA) players.} In: Data in Brief,
45 (2022). URL: https://doi.org/10.1016/j.dib.2022.108615
}
\examples{
data(nba)
cumhazPlot(Surv(survtime, cens) ~ 1, nba)
}
\keyword{datasets}
