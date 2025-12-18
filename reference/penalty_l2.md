# L2 penalty function (Ridge)

Creates a penalty function that computes the L2 norm squared (sum of
squares). Used for parameter shrinkage.

## Usage

``` r
penalty_l2(weights = NULL)
```

## Arguments

- weights:

  Optional parameter weights (default: all 1)

## Value

Penalty function

## Examples

``` r
penalty <- penalty_l2()
penalty(c(1, -2, 3))  # Returns 14
#> [1] 14

# Weighted L2
penalty <- penalty_l2(weights = c(1, 2, 1))
penalty(c(1, -2, 3))  # Returns 1^2 + (2*2)^2 + 3^2 = 26
#> [1] 26
```
