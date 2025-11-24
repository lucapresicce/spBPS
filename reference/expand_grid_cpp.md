# Build a grid from two vector (i.e. equivalent to `expand.grid()` in `R`)

Build a grid from two vector (i.e. equivalent to
[`expand.grid()`](https://rdrr.io/r/base/expand.grid.html) in `R`)

## Usage

``` r
expand_grid_cpp(x, y)
```

## Arguments

- x:

  [vector](https://rdrr.io/r/base/vector.html) first vector of numeric
  elements

- y:

  [vector](https://rdrr.io/r/base/vector.html) second vector of numeric
  elements

## Value

[matrix](https://rdrr.io/r/base/matrix.html) expanded grid of
combinations

## Examples

``` r
## Create a matrix from all combination of vectors
x <- seq(0, 10, length.out = 100)
y <- seq(-1, 1, length.out = 20)
grid <- expand_grid_cpp(x = x, y = y)
```
