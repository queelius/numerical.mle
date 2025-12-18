# Elastic net penalty (combination of L1 and L2)

Creates a penalty combining L1 and L2 norms. The parameter alpha
controls the balance: alpha=1 is pure LASSO, alpha=0 is pure Ridge.

## Usage

``` r
penalty_elastic_net(alpha = 0.5, weights = NULL)
```

## Arguments

- alpha:

  Balance between L1 and L2 (numeric in \[0,1\], default: 0.5)

- weights:

  Optional parameter weights (default: all 1)

## Value

Penalty function

## Examples

``` r
# Equal mix of L1 and L2
penalty <- penalty_elastic_net(alpha = 0.5)

# More L1 (more sparsity)
penalty <- penalty_elastic_net(alpha = 0.9)

# More L2 (more shrinkage)
penalty <- penalty_elastic_net(alpha = 0.1)
```
