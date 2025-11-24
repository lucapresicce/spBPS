# Draw from the conditional posterior predictive for a set of unobserved covariates

Draw from the conditional posterior predictive for a set of unobserved
covariates

## Usage

``` r
r_pred_cond_MvT(data, X_u, d_u, d_us, hyperpar, poster, post)
```

## Arguments

- data:

  [list](https://rdrr.io/r/base/list.html) two elements: first named
  \\Y\\, second named \\X\\

- X_u:

  [matrix](https://rdrr.io/r/base/matrix.html) unobserved instances
  covariate matrix

- d_u:

  [matrix](https://rdrr.io/r/base/matrix.html) unobserved instances
  distance matrix

- d_us:

  [matrix](https://rdrr.io/r/base/matrix.html) cross-distance between
  unobserved and observed instances matrix

- hyperpar:

  [list](https://rdrr.io/r/base/list.html) two elemets: first named
  \\\alpha\\, second named \\\phi\\

- poster:

  [list](https://rdrr.io/r/base/list.html) output from `fit_cpp_MvT`
  function

- post:

  [list](https://rdrr.io/r/base/list.html) output from `post_draws_MvT`
  function

## Value

[list](https://rdrr.io/r/base/list.html) posterior predictive samples
