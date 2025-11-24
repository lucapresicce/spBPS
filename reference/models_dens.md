# Return the CV predictive density evaluations for all the model combinations

Return the CV predictive density evaluations for all the model
combinations

## Usage

``` r
models_dens(data, priors, coords, hyperpar, useKCV, K)
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

- useKCV:

  if `TRUE` K-fold cross validation is used instead of LOOCV (no
  `default`)

- K:

  [integer](https://rdrr.io/r/base/integer.html) number of folds

## Value

[matrix](https://rdrr.io/r/base/matrix.html) posterior predictive
density evaluations (each columns represent a different model)
