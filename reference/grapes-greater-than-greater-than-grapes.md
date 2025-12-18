# Sequential Solver Composition

Chains two solvers sequentially. The result of the first solver becomes
the starting point for the second. This enables coarse-to-fine
strategies.

## Usage

``` r
s1 %>>% s2
```

## Arguments

- s1:

  First solver function

- s2:

  Second solver function

## Value

A new solver function that runs s1 then s2

## Examples

``` r
# Coarse-to-fine: grid search to find good region, then gradient ascent
strategy <- grid_search(n = 5) %>>% gradient_ascent()
#> Error in grid_search(n = 5): argument "lower" is missing, with no default

# Three-stage refinement
strategy <- grid_search(n = 3) %>>% gradient_ascent() %>>% newton_raphson()
#> Error in grid_search(n = 3): argument "lower" is missing, with no default
```
