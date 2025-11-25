# Maximum likelihood estimation via gradient ascent

Performs gradient ascent optimization to find the MLE. This method uses
the score function (gradient of log-likelihood) to iteratively improve
parameter estimates.

## Usage

``` r
mle_gradient_ascent(
  loglike,
  score,
  theta0,
  config = mle_config_gradient(),
  constraint = mle_constraint()
)
```

## Arguments

- loglike:

  Log-likelihood function taking theta as input

- score:

  Score function (gradient of log-likelihood) taking theta as input

- theta0:

  Initial parameter guess (numeric vector)

- config:

  Configuration object (mle_config_gradient or mle_config_linesearch).
  Use mle_config_linesearch() for adaptive step sizes (recommended), or
  mle_config_gradient() for fixed step size.

- constraint:

  Optional domain constraints (mle_constraint object)

## Value

mle_numerical object with class mle_gradient_ascent containing:

- theta.hat:

  MLE estimate

- loglike:

  Log-likelihood at MLE

- score:

  Score vector at MLE (should be near zero)

- info:

  Fisher information matrix (negative Hessian)

- sigma:

  Covariance matrix (inverse of Fisher information)

- iter:

  Number of iterations

- converged:

  Convergence status

- config:

  Configuration used

- path:

  Optimization path (if trace=TRUE in config)

## Examples

``` r
if (FALSE) { # \dontrun{
# Normal distribution MLE
data <- rnorm(100, mean = 5, sd = 2)

loglike <- function(theta) {
  sum(dnorm(data, mean = theta[1], sd = theta[2], log = TRUE))
}

score <- function(theta) {
  mu <- theta[1]
  sigma <- theta[2]
  c(
    sum((data - mu) / sigma^2),
    sum((data - mu)^2 / sigma^3 - 1/sigma)
  )
}

# With line search (recommended)
result <- mle_gradient_ascent(
  loglike = loglike,
  score = score,
  theta0 = c(0, 1),
  config = mle_config_linesearch(max_step = 1.0)
)

# With fixed step size
result <- mle_gradient_ascent(
  loglike = loglike,
  score = score,
  theta0 = c(0, 1),
  config = mle_config_gradient(eta = 0.1)
)

# With constraints (positive variance only)
constraint <- mle_constraint(
  support = function(theta) theta[2] > 0,
  project = function(theta) c(theta[1], max(theta[2], 1e-8))
)

result <- mle_gradient_ascent(
  loglike = loglike,
  score = score,
  theta0 = c(0, 1),
  config = mle_config_linesearch(),
  constraint = constraint
)

# Check convergence
print(result$converged)
print(result$theta.hat)
print(result$score)  # Should be near zero
} # }
```
