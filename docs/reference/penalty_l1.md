# L1 penalty function (LASSO)

Creates a penalty function that computes the L1 norm (sum of absolute
values). Used for sparsity-inducing regularization.

## Usage

``` r
penalty_l1(weights = NULL)
```

## Arguments

- weights:

  Optional parameter weights (default: all 1)

## Value

Penalty function

## Examples

``` r
penalty <- penalty_l1()
penalty(c(1, -2, 3))  # Returns 6
#> [1] 6

# Weighted L1
penalty <- penalty_l1(weights = c(1, 2, 1))
penalty(c(1, -2, 3))  # Returns 1*1 + 2*2 + 1*3 = 8
#> [1] 8
```
