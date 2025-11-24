# Compute the BPS weights by convex optimization

Compute the BPS weights by convex optimization

## Usage

``` r
CVXR_opt(scores)
```

## Arguments

- scores:

  [matrix](https://rdrr.io/r/base/matrix.html) \\N \times K\\ of
  expected predictive density evaluations for the K models considered

## Value

conv_opt [function](https://rdrr.io/r/base/function.html) to perform
convex optimiazion with CVXR R package
