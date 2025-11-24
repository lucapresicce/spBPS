# Compute the Euclidean distance matrix

Compute the Euclidean distance matrix

## Usage

``` r
arma_dist(X)
```

## Arguments

- X:

  [matrix](https://rdrr.io/r/base/matrix.html) (tipically of \\N\\
  coordindates on \\\mathbb{R}^2\\ )

## Value

[matrix](https://rdrr.io/r/base/matrix.html) distance matrix of the
elements of \\X\\

## Examples

``` r
## Compute the Distance matrix of dimension (n x n)
n <- 100
p <- 2
X <- matrix(runif(n*p), nrow = n, ncol = p)
distance.matrix <- arma_dist(X)
```
