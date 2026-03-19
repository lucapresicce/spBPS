# Compute the BPS spatial prediction given a set of stacking weights

Compute the BPS spatial prediction given a set of stacking weights

## Usage

``` r
BPS_pred_MvT(data, X_u, priors, coords, crd_u, hyperpar, W, R)
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
  \\\mu_B\\,\\V_r\\,\\\Psi\\,\\\nu\\

- coords:

  [matrix](https://rdrr.io/r/base/matrix.html) sample coordinates for X
  and Y

- crd_u:

  [matrix](https://rdrr.io/r/base/matrix.html) unboserved instances
  coordinates

- hyperpar:

  [list](https://rdrr.io/r/base/list.html) two elemets: first named
  \\\alpha\\, second named \\\phi\\

- W:

  [matrix](https://rdrr.io/r/base/matrix.html) set of stacking weights

- R:

  [integer](https://rdrr.io/r/base/integer.html) number of desired
  samples

## Value

[list](https://rdrr.io/r/base/list.html) BPS posterior predictive
samples
