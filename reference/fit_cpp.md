# Compute the parameters for the posteriors distribution of \\\beta\\ and \\\Sigma\\ (i.e. updated parameters)

Compute the parameters for the posteriors distribution of \\\beta\\ and
\\\Sigma\\ (i.e. updated parameters)

## Usage

``` r
fit_cpp(data, priors, coords, hyperpar)
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

## Value

[list](https://rdrr.io/r/base/list.html) posterior update parameters
