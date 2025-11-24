# Compute the BPS posterior samples given a set of stacking weights

Compute the BPS posterior samples given a set of stacking weights

## Usage

``` r
BPS_postdraws(data, priors, coords, hyperpar, W, R)
```

## Arguments

- data:

  [list](https://rdrr.io/r/base/list.html) two elements: first named
  \\Y\\, second named \\X\\

- priors:

  [list](https://rdrr.io/r/base/list.html) priors: named
  \\\mu_b\\,\\V_b\\,\\a\\,\\b\\

- coords:

  [matrix](https://rdrr.io/r/base/matrix.html) sample coordinates for X
  and Y

- hyperpar:

  [list](https://rdrr.io/r/base/list.html) two elemets: first named
  \\\delta\\, second named \\\phi\\

- W:

  [matrix](https://rdrr.io/r/base/matrix.html) set of stacking weights

- R:

  [integer](https://rdrr.io/r/base/integer.html) number of desired
  samples

## Value

[matrix](https://rdrr.io/r/base/matrix.html) BPS posterior samples
