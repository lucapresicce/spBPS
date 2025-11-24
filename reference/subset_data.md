# Function to subset data for meta-analysis

Function to subset data for meta-analysis

## Usage

``` r
subset_data(data, K)
```

## Arguments

- data:

  [list](https://rdrr.io/r/base/list.html) three elements: first named
  \\Y\\, second named \\X\\, third named \\crd\\

- K:

  [integer](https://rdrr.io/r/base/integer.html) number of desired
  subsets

## Value

[list](https://rdrr.io/r/base/list.html) subsets of data, and the set of
indexes

## Examples

``` r
## Create a list of K random subsets given a list with Y, X, and crd
n <- 100
p <- 3
q <- 2
X <- matrix(rnorm(n*p), nrow = n, ncol = p)
Y <- matrix(rnorm(n*q), nrow = n, ncol = q)
crd <- matrix(runif(n*2), nrow = n, ncol = 2)
subsets <- subset_data(data = list(Y = Y, X = X, crd = crd), K = 10)
```
