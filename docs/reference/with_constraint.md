# Quick constrained optimization

Convenience wrapper for constrained optimization with simplified
interface. Automatically creates constraint object from support and
projection functions.

## Usage

``` r
with_constraint(solver, support, project, ...)
```

## Arguments

- solver:

  Solver function (e.g., mle_grad, mle_nr)

- support:

  Support function (returns TRUE if theta is valid)

- project:

  Projection function (maps invalid theta to valid theta)

- ...:

  Arguments passed to solver

## Value

mle object from solver

## Examples

``` r
if (FALSE) { # \dontrun{
# Constrain parameters to be positive
result <- with_constraint(
  solver = mle_grad,
  support = function(theta) all(theta > 0),
  project = function(theta) pmax(theta, 1e-8),
  loglike = loglike,
  score = score,
  theta0 = c(1, 1)
)
} # }
```
