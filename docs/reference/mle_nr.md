# Quick Newton-Raphson with sensible defaults

Convenience wrapper for mle_newton_raphson with simplified interface.
Always uses line search for stability (recommended for Newton-Raphson).

## Usage

``` r
mle_nr(loglike, score, fisher, theta0, inverted = FALSE, ...)
```

## Arguments

- loglike:

  Log-likelihood function

- score:

  Score function

- fisher:

  Fisher information or covariance function

- theta0:

  Initial parameters

- inverted:

  Is fisher the covariance matrix? (default: FALSE)

- ...:

  Additional config parameters passed to mle_config_linesearch

## Value

mle_newton_raphson object

## Examples

``` r
if (FALSE) { # \dontrun{
# Quick usage with Fisher information matrix
result <- mle_nr(loglike, score, fisher, theta0 = c(0, 1))

# With covariance matrix (inverted FIM)
result <- mle_nr(
  loglike, score, covariance,
  theta0 = c(0, 1),
  inverted = TRUE
)

# Override config
result <- mle_nr(
  loglike, score, fisher,
  theta0 = c(0, 1),
  max_iter = 50,
  max_step = 0.5
)
} # }
```
