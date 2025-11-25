# MLE via random restarts

Attempts to find a global maximum by running a local solver from
multiple random starting points. This helps escape local maxima and find
better solutions when the likelihood surface is multimodal.

## Usage

``` r
mle_random_restart(loglike, solver, theta0_sampler, n_trials = 100L, ...)
```

## Arguments

- loglike:

  Log-likelihood function

- solver:

  Solver function (e.g., mle_gradient_ascent, mle_newton_raphson). Must
  accept loglike, theta0, and additional arguments.

- theta0_sampler:

  Function generating random initial parameters. Called without
  arguments, must return a valid theta0 vector.

- n_trials:

  Number of random trials to perform (integer, default: 100)

- ...:

  Additional arguments passed to solver

## Value

Best mle object found across all trials, with additional attribute
n_trials indicating the number of trials performed.

## Examples

``` r
if (FALSE) { # \dontrun{
# Multimodal likelihood with multiple local maxima
loglike <- function(theta) {
  # Mixture of two peaks
  -((theta[1]-5)^2 + (theta[2]-5)^2) / 10 -
   ((theta[1]+3)^2 + (theta[2]+3)^2) / 10
}

score <- function(theta) {
  c(
    -(theta[1]-5) / 5 - (theta[1]+3) / 5,
    -(theta[2]-5) / 5 - (theta[2]+3) / 5
  )
}

# Random sampler for initial points
sampler <- function() {
  runif(2, min = -10, max = 10)
}

# Try 50 random starting points
result <- mle_random_restart(
  loglike = loglike,
  solver = mle_gradient_ascent,
  theta0_sampler = sampler,
  n_trials = 50,
  score = score,
  config = mle_config_linesearch(max_iter = 50)
)

print(result$theta.hat)
print(result$loglike)
print(result$n_trials)
} # }
```
