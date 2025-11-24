# Perform prediction for BPS accelerated models - loop over prediction set

Perform prediction for BPS accelerated models - loop over prediction set

## Usage

``` r
spPredict_BPS(data, X_u, priors, coords, crd_u, hyperpar, W, R, J)
```

## Arguments

- data:

  [list](https://rdrr.io/r/base/list.html) two elements: first named
  \\Y\\, second named \\X\\

- X_u:

  [matrix](https://rdrr.io/r/base/matrix.html) unobserved instances
  covariate matrix

- priors:

  [list](https://rdrr.io/r/base/list.html) priors: named
  \\\mu_b\\,\\V_b\\,\\a\\,\\b\\

- coords:

  [matrix](https://rdrr.io/r/base/matrix.html) sample coordinates for X
  and Y

- crd_u:

  [matrix](https://rdrr.io/r/base/matrix.html) unboserved instances
  coordinates

- hyperpar:

  [list](https://rdrr.io/r/base/list.html) two elemets: first named
  \\\delta\\, second named \\\phi\\

- W:

  [matrix](https://rdrr.io/r/base/matrix.html) set of stacking weights

- R:

  [integer](https://rdrr.io/r/base/integer.html) number of desired
  samples

- J:

  [integer](https://rdrr.io/r/base/integer.html) number of desired
  partition of prediction set

## Value

[list](https://rdrr.io/r/base/list.html) BPS posterior predictive
samples
