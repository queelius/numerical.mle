# numerical.mle

An R package providing numerical maximum likelihood estimation (MLE)
solvers with a clean, composable API.

## Installation

Install from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("queelius/numerical.mle")
```

## Features

- **Type-safe configuration objects** for optimization parameters
- **Core solvers**: Gradient ascent and Newton-Raphson
- **Meta-solvers**: Grid search and random restart for global
  optimization
- **Function transformers**: Subsampling (stochastic gradient) and
  penalties (L1/L2/elastic net)
- **Constraint support**: Domain constraints with projection functions

## Quick Start

``` r
library(numerical.mle)

# Generate sample data
set.seed(42)
x <- rnorm(100, mean = 5, sd = 2)

# Define log-likelihood and score functions
loglike <- function(theta) {
  if (theta[2] <= 0) return(-Inf)
  sum(dnorm(x, mean = theta[1], sd = theta[2], log = TRUE))
}

score <- function(theta) {
  if (theta[2] <= 0) return(c(NA, NA))
  mu <- theta[1]
  sigma <- theta[2]
  n <- length(x)
  c(sum(x - mu) / sigma^2,
    -n / sigma + sum((x - mu)^2) / sigma^3)
}

# Fit using gradient ascent with line search
result <- mle_grad(loglike, score, theta0 = c(0, 1))
print(coef(result))
#> NULL
```

## Configuration Objects

Create reusable, type-safe configurations:

``` r
# Base configuration
config <- mle_config(max_iter = 500, abs_tol = 1e-8)

# Gradient descent with fixed step size (eta)
config_grad <- mle_config_gradient(eta = 0.01, max_iter = 1000)

# Gradient descent with backtracking line search
config_ls <- mle_config_linesearch(max_step = 1.0, backtrack_ratio = 0.5)
```

## Constrained Optimization

Define domain constraints:

``` r
# Positive parameters only
constraint <- mle_constraint(
  support = function(theta) all(theta > 0),
  project = function(theta) pmax(theta, 1e-8)
)

# Use with any solver
result <- mle_gradient_ascent(
  loglike = loglike,
  score = score,
  theta0 = c(1, 1),
  config = mle_config_linesearch(),
  constraint = constraint
)
```

## Meta-Solvers

### Grid Search

``` r
result <- mle_grid_search(

  loglike = loglike,
  lower = c(-10, 0.1),
  upper = c(10, 5),
  grid_size = 20
)
```

### Random Restart

``` r
result <- mle_random_restart(

  loglike = loglike,
  solver = function(theta0) mle_grad(loglike, score, theta0),
  theta0_sampler = function() c(runif(1, -10, 10), runif(1, 0.1, 5)),
  n_trials = 10
)
```

## Function Transformers

### Stochastic Gradient (Mini-batching)

``` r
# For large datasets
loglike_stochastic <- with_subsampling(loglike, data = x, subsample_size = 32)
```

### Regularization

``` r
# L2 regularization (Ridge)
loglike_l2 <- with_penalty(loglike, penalty_l2(), lambda = 0.1)

# L1 regularization (Lasso)
loglike_l1 <- with_penalty(loglike, penalty_l1(), lambda = 0.1)

# Elastic net
loglike_elastic <- with_penalty(loglike, penalty_elastic_net(alpha = 0.5), lambda = 0.1)
```

## Documentation

Full documentation: <https://queelius.github.io/numerical.mle/>

## License

MIT
