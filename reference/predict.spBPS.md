# Predict at new locations using a fitted spBPS model

Generates posterior predictive samples at new spatial locations using a
fitted `spBPS` object. Memory usage is controlled via `pred_batch_size`

## Usage

``` r
# S3 method for class 'spBPS'
predict(
  object,
  newdata,
  draws = 200L,
  pred_batch_size = 200L,
  n_cores = 1L,
  return_summary = FALSE,
  probs = c(0.025, 0.975),
  ...
)
```

## Arguments

- object:

  Object of class `"spBPS"` returned by
  [`spBPS()`](https://lucapresicce.github.io/spBPS/reference/spBPS.md).

- newdata:

  List with elements `X` (u x p covariate matrix) and `coords` (u x d
  coordinate matrix).

- draws:

  Integer. Number of posterior predictive draws. Default 200.

- pred_batch_size:

  Integer. Number of prediction sites processed per batch. Default 200.
  Smaller values use less RAM; larger values may be faster.

- n_cores:

  Integer. Number of OMP threads. Default 1.

- return_summary:

  Logical. If `FALSE` (default) return full draw arrays. If `TRUE`
  return summary statistics.

- probs:

  Numeric vector of length 2. Quantile probabilities for credible
  intervals when `return_summary = TRUE`. Default `c(0.025, 0.975)`.

- ...:

  Currently unused.

## Value

When `return_summary = FALSE`, a list with arrays `Wu`, `Yu`, `MY` of
dimension u x q x R. When `return_summary = TRUE`, a list with u x q
matrices `mu_Yu`, `sd_Yu`, `q_lo`, `q_hi`, `mu_Wu`.

## See also

[`spBPS`](https://lucapresicce.github.io/spBPS/reference/spBPS.md)

## Examples

``` r
# \donttest{
n <- 1000
p <- 2
q <- 1

Y <- matrix(rnorm(n * q), ncol = q)
X <- matrix(rnorm(n * p), ncol = p)
coords <- matrix(runif(n * 2), ncol = 2)

data    <- list(Y = Y, X = X)
priors  <- list(mu_B = matrix(0, nrow = p, ncol = q),
                V_r  = diag(10, p),
                Psi  = diag(1, q),
                nu   = 3)
hyperpar <- list(alpha = 0.5, phi = 1)

res <- spBPS(data, priors, coords, hyperpar, subset_size = 200)
#> 
#> ====================================================
#>          Welcome to spBPS Bayesian Engine
#> ====================================================
#> 
#> Partitioning data into K = 5 subsets ...
#> 
#> Computing local stacking weights over J = 1 models ...
#> Local weights computed in 0.1 s.
#> 
#> Computing global stacking weights over K = 5 partitions ...
#> Global weights computed in 0.0 s.
#> 
#> ====================================================
#>      spBPS pipeline completed successfully!
#> ====================================================
#>   Fitting:     0.1 s
#>   Combination: 0.0 s
#>   TOTAL:       0.1 s
#> 

u <- 500
p <- 2
q <- 1

X_u <- matrix(rnorm(u * p), ncol = p)
crd_u <- matrix(runif(u * 2), ncol = 2)

# Predictive sampling
pred <- predict(res, newdata=list(X=X_u, coords=crd_u),
                draws=200L, pred_batch_size=500L)
#> [predict.spBPS]  u = 500  draws = 200  cores = 1
#>   Completed in 0.3 s.

# Summary statistics only (memory-efficient for large u)
pred_sum <- predict(res, newdata=list(X=X_u, coords=crd_u),
                    draws=500L, return_summary=TRUE,
                    probs=c(0.025, 0.975))
#> [predict.spBPS]  u = 500  draws = 500  cores = 1
#>   Completed in 0.5 s.
# }

```
