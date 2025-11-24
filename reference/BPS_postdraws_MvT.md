# Compute the BPS posterior samples given a set of stacking weights

Compute the BPS posterior samples given a set of stacking weights

## Usage

``` r
BPS_postdraws_MvT(data, priors, coords, hyperpar, W, R, par)
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

- W:

  [matrix](https://rdrr.io/r/base/matrix.html) set of stacking weights

- R:

  [integer](https://rdrr.io/r/base/integer.html) number of desired
  samples

- par:

  if `TRUE` only \\\beta\\ and \\\Sigma\\ are sampled (\\\omega\\ is
  omitted)

## Value

[matrix](https://rdrr.io/r/base/matrix.html) BPS posterior samples
