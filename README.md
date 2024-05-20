# Univariate and Multivariate Accelerated Spatial Modeling by Bayesian Predictive Stacking

This package provides the principal functions to perform accelerated modeling for univariate and multivariate spatial regressions. The package is used mostly within the novel working paper *"Building Artificially Intelligent Geostatistical Systems Using Bayesian Predictive Stacking" ([**Luca Presicce**](https://lucapresicce.github.io/) and Sudipto Banerjee, 2024+)"*. In order to guarantee the reproducibility of scientific results, in the [Bayesian-Transfer-Learning-and-Divide-Conquer-Models-for-Massive-Spatial-Datasets](https://github.com/lucapresicce/Bayesian-Transfer-Learning-and-Divide-Conquer-Models-for-Massive-Spatial-Datasets) repository are also available all the scripts of code used for simulations, data analysis, and results presented in the Manuscript and its Supplemental material.


--------------------------------------------------------------------------------
## Roadmap

| Folder | Description |
| :--- | :---: |
| `R` | contains funtions in R |
| `src` | contains function in Rcpp/C++ |

--------------------------------------------------------------------------------
## Guided installation
Since the package is not already available on CRAN (already submitted, and hopefully soon available), we use the `devtools` R package to install. Then, check for its presence on your device, otherwise install it:
```{r, echo = F, eval = F, collapse = TRUE}
if (!require(devtools)) {
  install.packages("devtools", dependencies = TRUE)
}
```
Once you have installed *devtools*, we can proceed. Let's install the `spBPS` package!
```{r}
devtools::install_github("lucapresicce/spBPS")
```
Cool! You are ready to start, now you too could perform **_fast & feasible_** Bayesian geostatistical modeling!

<!--
## Tutorial for usage
-->

--------------------------------------------------------------------------------
## Contacts

| | |
| :--- | :---: |
| Author | Luca Presicce (l.presicce@campus.unimib.it) & Sudipto Banerjee (sudipto@ucla.edu) |
| Maintainer | Luca Presicce (l.presicce@campus.unimib.it) |
| Reference | [**Luca Presicce**](https://lucapresicce.github.io/) and Sudipto Banerjee (2024+) *"Building Artificially Intelligent Geostatistical Systems Using Bayesian Predictive Stacking"*  |

<!--
Maintainer: l.presicce@campus.unimib.it
Reference: **Luca Presicce** and Sudipto Banerjee (2024+) *"Accelerated Meta-Kriging for massive Spatial dataset"* 
-->

