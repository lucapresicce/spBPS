# Compute the BPS spatial prediction given a set of stacking weights

Compute the BPS spatial prediction given a set of stacking weights

## Usage

``` r
BPS_pred(data, X_u, priors, coords, crd_u, hyperpar, W, R)
```

## Arguments

- data:

  [list](https://rdrr.io/r/base/list.html) two elements: first named
  \\Y\\, second named \\X\\

- X_u:

  [matrix](https://rdrr.io/r/base/matrix.html) unobserved instances
  covariate matrix

- priors:

  [list](https://rdrr.io/r/base/list.html) priors: named
  \\\mu_b\\,\\V_b\\,\\a\\,\\b\\

- coords:

  [matrix](https://rdrr.io/r/base/matrix.html) sample coordinates for X
  and Y

- crd_u:

  [matrix](https://rdrr.io/r/base/matrix.html) unboserved instances
  coordinates

- hyperpar:

  [list](https://rdrr.io/r/base/list.html) two elemets: first named
  \\\delta\\, second named \\\phi\\

- W:

  [matrix](https://rdrr.io/r/base/matrix.html) set of stacking weights

- R:

  [integer](https://rdrr.io/r/base/integer.html) number of desired
  samples

## Value

[list](https://rdrr.io/r/base/list.html) BPS posterior predictive
samples

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

## Select competetive set of values for hyperparameters
delta_seq <- c(0.1, 0.2, 0.3)
phi_seq <- c(3, 4, 5)

## Fit local models
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

## Model combination weights between partitions using Bayesian Predictive Stacking
comb_bps <- BPS_combine(fit_list = fit_list, K = 10, rp = 1)
Wbps <- comb_bps$W
W_list <- comb_bps$W_list

## Generate prediction points
m <- 50
X_new <- matrix(rnorm(m*p), nrow = m, ncol = p)
crd_new <- matrix(runif(m*2), nrow = m, ncol = 2)

## Perform posterior predictive sampling
R <- 250
subset_ind <- sample(1:10, R, TRUE, Wbps)
predictions <- vector(length = R, mode = "list")
for (r in 1:R) {
  ind_s <- subset_ind[r]
  Ys <- matrix(data_part$Y_list[[ind_s]])
  Xs <- data_part$X_list[[ind_s]]
  crds <- data_part$crd_list[[ind_s]]
  Ws <- W_list[[ind_s]]
  result <- spBPS::BPS_pred(data = list(Y = Ys, X = Xs), coords = crds,
                            X_u = X_new, crd_u = crd_new,
                            priors = list(mu_b = matrix(rep(0, p)),
                                          V_b = diag(10, p),
                                          a = 2,
                                          b = 2),
                                          hyperpar = list(delta = delta_seq,
                                                          phi = phi_seq),
                                          W = Ws, R = 1)

  predictions[[r]] <- result}

# }
```
