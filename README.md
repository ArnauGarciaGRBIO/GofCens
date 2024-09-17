# GofCens: Goodness-of-Fit Methods for Complete and Right-Censored Data <img src="man/figures/GofCens_logo1.png" align="right" alt="" width="160" />


<!-- badges: start -->
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version-last-release/GofCens)](https://cran.r-project.org/package=GofCens)
[![](https://cranlogs.r-pkg.org/badges/grand-total/GofCens)](https://cran.r-project.org/package=GofCens)
[![Download counter](https://cranlogs.r-pkg.org/badges/GofCens)](https://cran.r-project.org/package=GofCens)
<!-- badges: end -->

The **GofCens** package include the following graphical tools and goodness-of-fit tests for complete and right-censored data: 
- Kolmogorov-Smirnov, Cramér-von Mises, and Anderson-Darling tests, which use the empirical distribution function for complete data and are extended for right-censored data.
- Generalized chi-squared-type test, which is based on the squared differences between observed and expected counts using random cells with right-censored data.
- A series of graphical tools such as probability or cumulative hazard plots to guide the decision about the most suitable parametric model for the data.

## Installation
**GofCens** can be installed from [CRAN](https://cran.r-project.org/):
```{r CRAN-instalation, eval = FALSE}
install.packages("GofCens")
```


## Minimal Example
To conduct goodness-of-fit tests with right censored data we can use the `KScens()`, `CvMcens()`, `ADcens()` and `chisqcens()` functions. We illustrate this by means of the `colon` dataset:
```{r, eval = FALSE}
# Kolmogorov-Smirnov
set.seed(123)
KScens(Surv(time, status) ~ 1, colon, distr = "norm")

# Cramér-von Mises
colonsamp <- colon[sample(nrow(colon), 300), ]
CvMcens(Surv(time, status) ~ 1, colonsamp, distr = "normal")

# Anderson-Darling
ADcens(Surv(time, status) ~ 1, colonsamp, distr = "normal")

# Generalized chi-squared-type test
chisqcens(Surv(time, status) ~ 1, colonsamp, M = 6, distr = "normal")
```
The graphical tools provide nice plots via the functions `cumhazPlot()`, `kmPlot()` and `probPlot()`. See several examples using the `nba` data set:
```{r, eval = FALSE}
data(nba)
cumhazPlot(Surv(survtime, cens) ~ 1, nba, distr = c("expo", "normal", "gumbel"))
kmPlot(Surv(survtime, cens) ~ 1, nba, distr = c("normal", "weibull", "lognormal"),
       prnt = FALSE)
probPlot(Surv(survtime, cens) ~ 1, nba, "lognorm", plots = c("PP", "QQ", "SP"),
         ggp = TRUE, m = matrix(1:3, nr = 1))
``` 

