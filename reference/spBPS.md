# Unified spatial BPS workflow (multivariate path, works for q = 1)

Orchestrates subsetting, local stacking weight estimation, global
stacking combination, and optional posterior or predictive simulation
using Double Bayesian Predictive Stacking for latent spatial regression.
Works for both multivariate outcomes and the univariate case via
`q = 1`.

## Usage

``` r
spBPS(
  data,
  priors,
  coords,
  hyperpar,
  subset_size = 500L,
  K = NULL,
  cv_folds = 5L,
  rp = 1,
  combine_method = c("bps", "pseudoBMA"),
  draws = 0L,
  newdata = NULL,
  include_latent = FALSE,
  n_cores = 1L,
  pred_batch_size = 200L
)
```

## Arguments

- data:

  List with matrices `Y` (response, n x q) and `X` (covariates, n x p).

- priors:

  List of priors for the multivariate model: `mu_B` (p x q mean matrix),
  `V_r` (p x p covariance), `Psi` (q x q scale), `nu` (degrees of
  freedom).

- coords:

  Matrix of observation coordinates (n x d).

- hyperpar:

  List with elements `alpha` and `phi` (vectors allowed); defines the
  grid of models over which stacking weights are computed.

- subset_size:

  Target subset size when `K` is not provided. Default 500.

- K:

  Optional number of subsets. When `NULL`, computed as
  `ceiling(nrow(Y) / subset_size)` and lower-bounded at 1.

- cv_folds:

  Number of folds for local cross-validation. Default 5.

- rp:

  Fraction of rows used when recomputing global stacking weights (passed
  to `BPS_combine`). Ignored when `combine_method = "pseudoBMA"`.

- combine_method:

  Global combination method: Bayesian Predictive Stacking (`"bps"`) or
  pseudo-BMA (`"pseudoBMA"`).

- draws:

  Number of posterior/predictive draws to return (0 to skip sampling).
  When positive and `newdata` is supplied, joint posterior and
  predictive draws are returned. When positive and `newdata` is `NULL`,
  only posterior draws (beta, sigma) are returned.

- newdata:

  Optional list with `X` (u x p) and `coords` (u x d) for prediction
  locations. Required when `draws > 0` and predictions are desired. For
  large `u`, use `pred_batch_size` to control memory usage.

- include_latent:

  Unused; kept for compatibility.

- n_cores:

  Number of cores for parallel computation. Controls both the local
  weight estimation (parallel over K subsets) and predictive sampling.
  Default 1.

- pred_batch_size:

  Batch size for streaming prediction. Controls the maximum number of
  prediction sites processed at once, hence the peak memory usage.
  Default 200. Rule of thumb: peak RAM (MB) ? batch_size x q x draws x 8
  / 1e6. Set `NULL` for automatic selection (min(200, u)).

## Value

Object of class `"spBPS"` ? a list with components:

- subsets:

  Partition information: Y_list, X_list, crd_list.

- weights_global:

  K-vector of global stacking weights.

- weights_local:

  K-list of local stacking weight vectors.

- epd:

  K-list of n_k x J log-density matrices.

- priors, hyperpar:

  Stored for use by
  [`predict.spBPS`](https://lucapresicce.github.io/spBPS/reference/predict.spBPS.md).

- timings:

  Named numeric vector: fitting, combination, sampling (seconds).

- posterior:

  (if draws \> 0) List with three-dimensional arrays: `beta` (p x q x
  R), `sigma` (q x q x R), `model` (R).

- predictive:

  (if draws \> 0 and newdata supplied) List with three-dimensional
  arrays: `Wu`, `Yu`, `MY` (each u x q x R).

## See also

[`predict.spBPS`](https://lucapresicce.github.io/spBPS/reference/predict.spBPS.md)
for generating predictions on new data after fitting.

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
# }
```
