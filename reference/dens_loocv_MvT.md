# Compute the LOOCV of the density evaluations for fixed values of the hyperparameters

Compute the LOOCV of the density evaluations for fixed values of the
hyperparameters

## Usage

``` r
dens_loocv_MvT(data, priors, coords, hyperpar)
```

## Arguments

- data:

  [list](https://rdrr.io/r/base/list.html) two elements: first named
  \\Y\\, second named \\X\\

- priors:

  [list](https://rdrr.io/r/base/list.html) priors: named
  \\\mu_B\\,\\V_r\\,\\\Psi\\,\\\nu\\

- coords:

  [matrix](https://rdrr.io/r/base/matrix.html) sample coordinates for X
  and Y

- hyperpar:

  [list](https://rdrr.io/r/base/list.html) two elemets: first named
  \\\alpha\\, second named \\\phi\\

## Value

[vector](https://rdrr.io/r/base/vector.html) posterior predictive
density evaluations
