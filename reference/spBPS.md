# Unified spatial BPS workflow (multivariate path, works for q = 1)

Orchestrates subsetting, local stacking weight estimation, global
stacking combination, and optional posterior or predictive simulation
using the multivariate Student-t spatial model. Works for both
multivariate outcomes and the univariate case via `q = 1`.

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
  cores = NULL
)
```

## Arguments

- data:

  List with matrices `Y` (response) and `X` (covariates).

- priors:

  List of priors for the multivariate model (`mu_B`, `V_r`, `Psi`,
  `nu`).

- coords:

  Matrix of observation coordinates.

- hyperpar:

  List with elements `alpha` and `phi` (vectors allowed).

- subset_size:

  Target subset size when `K` is not provided. Default 500.

- K:

  Optional number of subsets. When `NULL`, computed as
  `ceiling(nrow(Y) / subset_size)` and lower-bounded at 1.

- cv_folds:

  Number of folds for local cross-validation (default 5).

- rp:

  Fraction of rows used when recomputing global stacking weights (passed
  to `BPS_combine`). Ignored when `combine_method = "pseudoBMA"`.

- combine_method:

  Choose between Bayesian Predictive Stacking (`"bps"`) or pseudo-BMA
  (`"pseudoBMA"`) for combining subsets.

- draws:

  Number of joint posterior/predictive draws to return (0 to skip). When
  positive, `newdata` must be supplied because draws are obtained via
  `BPS_post_MvT` which jointly samples posterior and predictive.

- newdata:

  Optional list with `X` and `coords` for prediction locations; required
  when either draw count is positive.

- include_latent:

  Logical; if `TRUE`, posterior draws include latent processes.

- cores:

  Optional integer; when \>1 a parallel backend is registered internally
  via `doParallel::registerDoParallel(cores)` for the fit and draw
  loops. When `NULL`, the existing foreach backend (if any) is used.

## Value

List with components `subsets`, `weights_global`, `weights_local`,
`epd`, and optional `posterior` and `predictive` draws.

## Examples

``` r
# \donttest{
n <- 1000
p <- 2
q <- 1

Y <- matrix(rnorm(n*q), ncol = q)
X <- matrix(rnorm(n*p), ncol = p)
coords <- matrix(runif(n*2), ncol = 2)

data <- list(Y = Y, X = X)
priors <- list(mu_B = matrix(0, nrow = p, ncol = q),
                             V_r = diag(10, p),
                             Psi = diag(1, q),
                             nu = 3)
hyperpar <- list(alpha = 0.5, phi = 1)
subset_size <- 200

res <- spBPS(data, priors, coords, hyperpar, subset_size = subset_size)
#> 
#> ====================================================
#>          Welcome to spBPS Bayesian Engine
#> ====================================================
#> 
#> Pritioning data into K = 5 subsets ... 
#> 
#> Computing local stacking weights over J = 1 models ...
#> Local weights computed.
#> 
#> Computing global stacking weights over K = 5 partitions ...
#> Global weights computed.
#> 
#> ====================================================
#>      spBPS pipeline completed successfully!
#> ====================================================
#> 

# }
```
