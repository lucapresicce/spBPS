# Evaluate the density of a set of unobserved response with respect to the conditional posterior predictive

Evaluate the density of a set of unobserved response with respect to the
conditional posterior predictive

## Usage

``` r
d_pred_cpp(data, X_u, Y_u, d_u, d_us, hyperpar, poster)
```

## Arguments

- data:

  [list](https://rdrr.io/r/base/list.html) two elements: first named
  \\Y\\, second named \\X\\

- X_u:

  [matrix](https://rdrr.io/r/base/matrix.html) unobserved instances
  covariate matrix

- Y_u:

  [matrix](https://rdrr.io/r/base/matrix.html) unobserved instances
  response matrix

- d_u:

  [matrix](https://rdrr.io/r/base/matrix.html) unobserved instances
  distance matrix

- d_us:

  [matrix](https://rdrr.io/r/base/matrix.html) cross-distance between
  unobserved and observed instances matrix

- hyperpar:

  [list](https://rdrr.io/r/base/list.html) two elemets: first named
  \\\delta\\, second named \\\phi\\

- poster:

  [list](https://rdrr.io/r/base/list.html) output from `fit_cpp`
  function

## Value

[vector](https://rdrr.io/r/base/vector.html) posterior predictive
density evaluations
