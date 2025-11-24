# Function to subset data for meta-analysis

Function to subset data for meta-analysis

## Usage

``` r
forceSymmetry_cpp(mat)
```

## Arguments

- mat:

  [matrix](https://rdrr.io/r/base/matrix.html) not-symmetric matrix

## Value

[matrix](https://rdrr.io/r/base/matrix.html) symmetric matrix (lower
triangular of `mat` is used)

## Examples

``` r
## Force matrix to be symmetric (avoiding numerical problems)
n <- 4
X <- matrix(runif(n*n), nrow = n, ncol = n)
X <- forceSymmetry_cpp(mat = X)
```
