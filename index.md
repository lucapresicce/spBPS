# spBPS

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/spBPS?color=blue)](https://CRAN.R-project.org/package=spBPS)
[![Downloads](https://cranlogs.r-pkg.org/badges/spBPS?color=orange)](https://CRAN.R-project.org/package=spBPS)
[![CRAN RStudio mirror
downloads](https://cranlogs.r-pkg.org/badges/grand-total/spBPS?color=green)](https://CRAN.R-project.org/package=spBPS)

## Overview

This package provides the principal functions to perform accelerated
modeling for univariate and multivariate spatial regressions. The
package is used mostly within the novel working paper *“Bayesian
Transfer Learning for Artificially Intelligent Geospatial Systems: A
Predictive Stacking Approach” ([**Luca
Presicce**](https://lucapresicce.github.io/) and Sudipto Banerjee,
2024+)“*. To guarantee the reproducibility of scientific results, in the
[Bayesian-Transfer-Learning-for-GeoAI](https://github.com/lucapresicce/Bayesian-Transfer-Learning-for-GeoAI)
repository are also available all the scripts of code used for
simulations, data analysis, and results presented in the Manuscript and
its Supplemental material.

## Installation

If installing from CRAN, use the following.

``` r
install.packages("spBPS")
```

For a quick installation of the development version, run the following
command in R. We use the `devtools` R package to install. Then, check
for its presence on your device, otherwise install it:

``` r
if (!require(devtools)) {
  install.packages("devtools", dependencies = TRUE)
}
```

Once you have installed *devtools*, we can proceed. Let’s install the
`spBPS` package!

``` r
devtools::install_github("lucapresicce/spFFBS")
```

## Usage

Once successfully installed, load the library in R.

``` r
library(spBPS)
```

Cool! You are ready to start, now you too could perform ***fast &
feasible*** Bayesian geostatistical modeling!

## Contacts

|            |                                                                                                                                                                                                  |
|:-----------|:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------:|
| Author     |                                                      Luca Presicce (<l.presicce@campus.unimib.it>) & Sudipto Banerjee (<sudipto@ucla.edu>)                                                       |
| Maintainer |                                                                          Luca Presicce (<l.presicce@campus.unimib.it>)                                                                           |
| Reference  | [**Luca Presicce**](https://lucapresicce.github.io/) and Sudipto Banerjee (2024+) *“Bayesian Transfer Learning for Artificially Intelligent Geospatial Systems: A Predictive Stacking Approach”* |
