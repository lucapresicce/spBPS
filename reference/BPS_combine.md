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
