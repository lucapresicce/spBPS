# Draw from the joint posterior predictive for a set of unobserved covariates

Draw from the joint posterior predictive for a set of unobserved
covariates

## Usage

``` r
r_pred_joint(data, X_u, d_u, d_us, hyperpar, poster, R)
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
  \\\delta\\, second named \\\phi\\

- poster:

  [list](https://rdrr.io/r/base/list.html) output from `fit_cpp`
  function

- R:

  [integer](https://rdrr.io/r/base/integer.html) number of posterior
  predictive samples

## Value

[list](https://rdrr.io/r/base/list.html) posterior predictive samples
