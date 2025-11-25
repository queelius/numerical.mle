# Internal direction-based optimizer

Core optimization algorithm used by gradient-based solvers. This is an
internal function not meant to be called directly by users.

## Usage

``` r
.mle_optimize_direction(
  loglike,
  direction_fn,
  theta0,
  config,
  constraint,
  use_linesearch
)
```

## Arguments

- loglike:

  Log-likelihood function

- direction_fn:

  Function computing search direction

- theta0:

  Initial parameters

- config:

  Configuration object (mle_config or subclass)

- constraint:

  Domain constraints (mle_constraint object)

- use_linesearch:

  Whether to use backtracking line search

## Value

mle_numerical object with optimization results
