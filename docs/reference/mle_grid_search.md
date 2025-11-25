# MLE via grid search

Performs exhaustive grid search over a bounded parameter space.
Optionally refines each grid point using a local solver.

## Usage

``` r
mle_grid_search(loglike, lower, upper, grid_size, refine_solver = NULL, ...)
```

## Arguments

- loglike:

  Log-likelihood function

- lower:

  Lower bounds for parameters (numeric vector)

- upper:

  Upper bounds for parameters (numeric vector)

- grid_size:

  Grid resolution. Either a single integer (same resolution per
  dimension) or a vector of integers (one per dimension).

- refine_solver:

  Optional local solver to refine each grid point. If NULL, only
  evaluates loglike at grid points without refinement.

- ...:

  Additional arguments passed to refine_solver

## Value

mle object with best solution found, including:

- theta.hat:

  Best parameter estimate

- loglike:

  Log-likelihood at best point

- grid_size:

  Grid resolution used

- n_evaluated:

  Number of grid points evaluated

## Examples

``` r
if (FALSE) { # \dontrun{
# Simple grid search without refinement
loglike <- function(theta) {
  -(theta[1]^2 + theta[2]^2)
}

result <- mle_grid_search(
  loglike = loglike,
  lower = c(-5, -5),
  upper = c(5, 5),
  grid_size = 20
)

# Grid search with local refinement
score <- function(theta) {
  -2 * theta
}

result <- mle_grid_search(
  loglike = loglike,
  lower = c(-5, -5),
  upper = c(5, 5),
  grid_size = 10,
  refine_solver = mle_gradient_ascent,
  score = score,
  config = mle_config_gradient(eta = 0.1, max_iter = 20)
)

# Different resolution per dimension
result <- mle_grid_search(
  loglike = loglike,
  lower = c(-5, -5),
  upper = c(5, 5),
  grid_size = c(20, 10)  # 20 points in dim 1, 10 in dim 2
)
} # }
```
