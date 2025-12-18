# Parallel Solver Racing

Runs multiple solvers and returns the best result (highest
log-likelihood). Useful when unsure which method will work best for a
given problem.

## Usage

``` r
s1 %|% s2
```

## Arguments

- s1:

  First solver function

- s2:

  Second solver function

## Value

A new solver function that runs both and picks the best

## Examples

``` r
# Race gradient-based vs derivative-free
strategy <- gradient_ascent() %|% nelder_mead()

# Race multiple methods
strategy <- gradient_ascent() %|% bfgs() %|% nelder_mead()
```
