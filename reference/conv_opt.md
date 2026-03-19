# Solver for Bayesian Predictive Stacking of Predictive densities convex optimization problem

Solver for Bayesian Predictive Stacking of Predictive densities convex
optimization problem

## Usage

``` r
conv_opt(scores)
```

## Arguments

- scores:

  [matrix](https://rdrr.io/r/base/matrix.html) \\N \times K\\ of
  expected predictive density evaluations for the K models considered

## Value

W [matrix](https://rdrr.io/r/base/matrix.html) of Bayesian Predictive
Stacking weights for the K models considered
