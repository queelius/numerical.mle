# Maximum likelihood estimation via Newton-Raphson

Performs Newton-Raphson optimization to find the MLE. This second-order
method uses both the score (gradient) and Fisher information matrix
(Hessian) to achieve faster convergence than gradient ascent.

## Usage

``` r
mle_newton_raphson(
  loglike,
  score,
  fisher,
  theta0,
  config = mle_config_linesearch(),
  constraint = mle_constraint(),
  inverted = FALSE
)
```

## Arguments

- loglike:

  Log-likelihood function taking theta as input

- score:

  Score function (gradient of log-likelihood) taking theta as input

- fisher:

  Fisher information matrix function. Either the FIM itself or its
  inverse (covariance matrix), depending on the `inverted` parameter.

- theta0:

  Initial parameter guess (numeric vector)

- config:

  Configuration object (typically mle_config_linesearch). Newton-Raphson
  benefits from line search to ensure stability.

- constraint:

  Optional domain constraints (mle_constraint object)

- inverted:

  Logical. If TRUE, `fisher` is the covariance matrix (inverse of FIM).
  If FALSE (default), `fisher` is the FIM.

## Value

mle_numerical object with class mle_newton_raphson containing:

- theta.hat:

  MLE estimate

- loglike:

  Log-likelihood at MLE

- score:

  Score vector at MLE (should be near zero)

- info:

  Fisher information matrix

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
# Normal distribution MLE with Newton-Raphson
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

fisher <- function(theta) {
  n <- length(data)
  sigma <- theta[2]
  matrix(c(
    n / sigma^2, 0,
    0, 2*n / sigma^2
  ), nrow = 2)
}

# Standard usage with FIM
result <- mle_newton_raphson(
  loglike = loglike,
  score = score,
  fisher = fisher,
  theta0 = c(0, 1),
  config = mle_config_linesearch()
)

# Using inverted FIM (covariance matrix)
covariance <- function(theta) {
  MASS::ginv(fisher(theta))
}

result <- mle_newton_raphson(
  loglike = loglike,
  score = score,
  fisher = covariance,
  theta0 = c(0, 1),
  inverted = TRUE
)

# With constraints
constraint <- mle_constraint(
  support = function(theta) theta[2] > 0,
  project = function(theta) c(theta[1], max(theta[2], 1e-8))
)

result <- mle_newton_raphson(
  loglike = loglike,
  score = score,
  fisher = fisher,
  theta0 = c(0, 1),
  config = mle_config_linesearch(),
  constraint = constraint
)

# Faster convergence than gradient ascent
print(result$iter)  # Typically fewer iterations
print(result$score)  # Should be very close to zero
} # }
```
