# compositional.mle

An R package for **composable maximum likelihood estimation**. Solvers
are first-class functions that combine via sequential chaining, parallel
racing, and random restarts.

## Installation

``` r
# install.packages("devtools")
devtools::install_github("queelius/compositional.mle")
```

## Design Philosophy

Following SICP principles, the package provides:

1.  **Primitive solvers** -
    [`gradient_ascent()`](https://queelius.github.io/compositional.mle/reference/gradient_ascent.md),
    [`newton_raphson()`](https://queelius.github.io/compositional.mle/reference/newton_raphson.md),
    [`bfgs()`](https://queelius.github.io/compositional.mle/reference/bfgs.md),
    [`nelder_mead()`](https://queelius.github.io/compositional.mle/reference/nelder_mead.md),
    etc.
2.  **Composition operators** - `%>>%` (sequential), `%|%` (race),
    [`with_restarts()`](https://queelius.github.io/compositional.mle/reference/with_restarts.md)
3.  **Closure property** - Combining solvers yields a solver

## Quick Start

``` r
library(compositional.mle)

# Generate sample data
set.seed(42)
x <- rnorm(100, mean = 5, sd = 2)

# Define the problem (separate from solver strategy)
problem <- mle_problem(
  loglike = function(theta) {
    if (theta[2] <= 0) return(-Inf)
    sum(dnorm(x, theta[1], theta[2], log = TRUE))
  },
  score = function(theta) {
    mu <- theta[1]; sigma <- theta[2]; n <- length(x)
    c(sum(x - mu) / sigma^2,
      -n / sigma + sum((x - mu)^2) / sigma^3)
  },
  constraint = mle_constraint(
    support = function(theta) theta[2] > 0,
    project = function(theta) c(theta[1], max(theta[2], 1e-8))
  )
)

# Simple solve
result <- gradient_ascent()(problem, theta0 = c(0, 1))
result$theta.hat
#> [1] 5.065030 2.072274
```

## Composing Solvers

### Sequential Chaining (`%>>%`)

Chain solvers for coarse-to-fine optimization:

``` r
# Grid search -> gradient ascent -> Newton-Raphson
strategy <- grid_search(lower = c(-10, 0.5), upper = c(10, 5), n = 5) %>>%
  gradient_ascent(max_iter = 50) %>>%
  newton_raphson(max_iter = 20)

result <- strategy(problem, theta0 = c(0, 1))
result$theta.hat
#>     Var1     Var2 
#> 5.065030 2.072274
```

### Parallel Racing (`%|%`)

Race multiple methods, keep the best:

``` r
# Try multiple approaches, pick winner by log-likelihood
strategy <- gradient_ascent() %|% bfgs() %|% nelder_mead()

result <- strategy(problem, theta0 = c(0, 1))
c(result$theta.hat, loglike = result$loglike)
#>                             loglike 
#>    5.065030    2.072274 -214.758518
```

### Random Restarts

Escape local optima with multiple starting points:

``` r
strategy <- with_restarts(
  gradient_ascent(),
  n = 10,
  sampler = uniform_sampler(c(-10, 0.5), c(10, 5))
)

result <- strategy(problem, theta0 = c(0, 1))
result$theta.hat
#> [1] 5.065030 2.072274
```

## Available Solvers

| Factory                                                                                          | Method                           | Requires      |
|--------------------------------------------------------------------------------------------------|----------------------------------|---------------|
| [`gradient_ascent()`](https://queelius.github.io/compositional.mle/reference/gradient_ascent.md) | Steepest ascent with line search | score         |
| [`newton_raphson()`](https://queelius.github.io/compositional.mle/reference/newton_raphson.md)   | Second-order Newton              | score, fisher |
| [`bfgs()`](https://queelius.github.io/compositional.mle/reference/bfgs.md)                       | Quasi-Newton BFGS                | score         |
| [`lbfgsb()`](https://queelius.github.io/compositional.mle/reference/lbfgsb.md)                   | L-BFGS-B with box constraints    | score         |
| [`nelder_mead()`](https://queelius.github.io/compositional.mle/reference/nelder_mead.md)         | Simplex (derivative-free)        | \-            |
| [`grid_search()`](https://queelius.github.io/compositional.mle/reference/grid_search.md)         | Exhaustive grid                  | \-            |
| [`random_search()`](https://queelius.github.io/compositional.mle/reference/random_search.md)     | Random sampling                  | \-            |

## Function Transformers

``` r
# Stochastic gradient (mini-batching)
loglike_sgd <- with_subsampling(loglike, data = x, subsample_size = 32)

# Regularization
loglike_l2 <- with_penalty(loglike, penalty_l2(), lambda = 0.1)
loglike_l1 <- with_penalty(loglike, penalty_l1(), lambda = 0.1)
```

## Documentation

Full documentation: <https://queelius.github.io/compositional.mle/>

## License

MIT
