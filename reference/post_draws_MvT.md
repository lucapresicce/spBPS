# Sample R draws from the posterior distributions

Sample R draws from the posterior distributions

## Usage

``` r
post_draws_MvT(poster, R, par, p)
```

## Arguments

- poster:

  [list](https://rdrr.io/r/base/list.html) output from `fit_cpp`
  function

- R:

  [integer](https://rdrr.io/r/base/integer.html) number of posterior
  samples

- par:

  if `TRUE` only \\\beta\\ and \\\Sigma\\ are sampled (\\\omega\\ is
  omitted)

- p:

  [integer](https://rdrr.io/r/base/integer.html) if `par = TRUE`, it
  specifies the column number of \\X\\

## Value

[list](https://rdrr.io/r/base/list.html) posterior samples
