# Quick gradient ascent with sensible defaults

Convenience wrapper for mle_gradient_ascent with simplified interface.
Automatically enables line search for better convergence.

## Usage

``` r
mle_grad(loglike, score, theta0, use_linesearch = TRUE, ...)
```

## Arguments

- loglike:

  Log-likelihood function

- score:

  Score function (gradient)

- theta0:

  Initial parameters

- use_linesearch:

  Use backtracking line search (default: TRUE)

- ...:

  Additional config parameters passed to mle_config_linesearch or
  mle_config_gradient

## Value

mle_gradient_ascent object

## Examples

``` r
if (FALSE) { # \dontrun{
# Quick usage with defaults
result <- mle_grad(loglike, score, theta0 = c(0, 1))

# Override config parameters
result <- mle_grad(
  loglike, score, theta0 = c(0, 1),
  max_iter = 200,
  rel_tol = 1e-6
)

# Without line search
result <- mle_grad(
  loglike, score, theta0 = c(0, 1),
  use_linesearch = FALSE,
  eta = 0.1
)
} # }
```
