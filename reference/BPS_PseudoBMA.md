# Combine subset models wiht Pseudo-BMA

Combine subset models wiht Pseudo-BMA

## Usage

``` r
BPS_PseudoBMA(fit_list)
```

## Arguments

- fit_list:

  [list](https://rdrr.io/r/base/list.html) K fitted model outputs
  composed by two elements each: first named \\epd\\, second named \\W\\

## Value

[matrix](https://rdrr.io/r/base/matrix.html) posterior predictive
density evaluations (each columns represent a different model)
