# Compute the BPS weights by convex optimization

Compute the BPS weights by convex optimization

## Usage

``` r
BPS_weights_MvT(data, priors, coords, hyperpar, K)
```

## Arguments

- data:

  [list](https://rdrr.io/r/base/list.html) two elements: first named
  \\Y\\, second named \\X\\

- priors:

  [list](https://rdrr.io/r/base/list.html) priors: named
  \\\mu_B\\,\\V_r\\,\\\Psi\\,\\\nu\\

- coords:

  [matrix](https://rdrr.io/r/base/matrix.html) sample coordinates for X
  and Y

- hyperpar:

  [list](https://rdrr.io/r/base/list.html) two elemets: first named
  \\\alpha\\, second named \\\phi\\

- K:

  [integer](https://rdrr.io/r/base/integer.html) number of folds

## Value

[matrix](https://rdrr.io/r/base/matrix.html) posterior predictive
density evaluations (each columns represent a different model)

## Examples

``` r
# \donttest{
## Generate subsets of data
n <- 100
p <- 3
q <- 2
X <- matrix(rnorm(n*p), nrow = n, ncol = p)
Y <- matrix(rnorm(n*q), nrow = n, ncol = q)
crd <- matrix(runif(n*2), nrow = n, ncol = 2)

## Select competitive set of values for hyperparameters
alfa_seq <- c(0.7, 0.8, 0.9)
phi_seq <- c(3, 4, 5)

## Perform Bayesian Predictive Stacking within subsets
bps <- spBPS::BPS_weights_MvT(data = list(Y = Y, X = X),
                              priors = list(mu_B = matrix(0, nrow = p, ncol = q),
                                            V_r = diag(10, p),
                                            Psi = diag(1, q),
                                            nu = 3), coords = crd,
                                            hyperpar = list(alpha = alfa_seq,
                                                            phi = phi_seq),
                                            K = 5)
# }
```
