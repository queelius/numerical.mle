# Backtracking line search step

Performs a backtracking line search to find a step size that improves
the objective function.

## Usage

``` r
.backtracking_step(
  loglike,
  direction,
  theta_current,
  max_step,
  constraint,
  backtrack_ratio,
  max_iter,
  min_step,
  debug
)
```

## Arguments

- loglike:

  Log-likelihood function

- direction:

  Search direction vector

- theta_current:

  Current parameter values

- max_step:

  Maximum step size

- constraint:

  Domain constraints

- backtrack_ratio:

  Backtracking multiplier (0 \< r \< 1)

- max_iter:

  Maximum iterations for line search

- min_step:

  Minimum step size threshold

- debug:

  Print debug information

## Value

List with success (logical) and theta (new parameter values)
