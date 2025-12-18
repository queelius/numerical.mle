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

List with optim()-compatible format:

- par:

  Final parameter estimates

- value:

  Log-likelihood at solution

- convergence:

  0 if converged, 1 otherwise

- iterations:

  Number of iterations taken

- path:

  Optimization path (if trace=TRUE)

## Details

Returns an optim()-compatible list that can be passed to
algebraic.mle::mle_numerical().
