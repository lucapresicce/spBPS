# Combine subset models wiht BPS

Combine subset models wiht BPS

## Usage

``` r
BPS_combine(fit_list, K, rp)
```

## Arguments

- fit_list:

  [list](https://rdrr.io/r/base/list.html) K fitted model outputs
  composed by two elements each: first named \\epd\\, second named \\W\\

- K:

  [integer](https://rdrr.io/r/base/integer.html) number of folds

- rp:

  [double](https://rdrr.io/r/base/double.html) percentage of
  observations to take into account for optimization (`default=1`)

## Value

[matrix](https://rdrr.io/r/base/matrix.html) posterior predictive
density evaluations (each columns represent a different model)

## Examples

``` r
# \donttest{
## Generate subsets of data
n <- 100
p <- 3
X <- matrix(rnorm(n*p), nrow = n, ncol = p)
Y <- matrix(rnorm(n), nrow = n, ncol = 1)
crd <- matrix(runif(n*2), nrow = n, ncol = 2)
data_part <- subset_data(data = list(Y = Y, X = X, crd = crd), K = 10)

## Select competitive set of values for hyperparameters
delta_seq <- c(0.1, 0.2, 0.3)
phi_seq <- c(3, 4, 5)

## Perform Bayesian Predictive Stacking within subsets
fit_list <- vector(length = 10, mode = "list")
for (i in 1:10) {
    Yi <- data_part$Y_list[[i]]
    Xi <- data_part$X_list[[i]]
    crd_i <- data_part$crd_list[[i]]
    p <- ncol(Xi)
    bps <- spBPS::BPS_weights(data = list(Y = Yi, X = Xi),
                               priors = list(mu_b = matrix(rep(0, p)),
                                             V_b = diag(10, p),
                                             a = 2,
                                             b = 2), coords = crd_i,
                                             hyperpar = list(delta = delta_seq,
                                                             phi = phi_seq),
                                             K = 5)
     w_hat <- bps$W
     epd <- bps$epd
     fit_list[[i]] <- list(epd, w_hat) }

## Combination weights between partitions using Bayesian Predictive Stacking
comb_bps <- BPS_combine(fit_list = fit_list, K = 10, rp = 1)
# }
```
