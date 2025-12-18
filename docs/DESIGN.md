# numerical.mle Design Document

## Philosophy

Following SICP principles: 1. **Primitive expressions** - Basic solvers
(gradient ascent, Newton-Raphson, etc.) 2. **Means of combination** -
Composition operators (`%>>%`, `%|%`, `with_restarts`) 3. **Means of
abstraction** - Solver factories, problem specification

**Key property**: Closure - combining solvers yields a solver.

## Core Abstractions

### 1. mle_problem

Encapsulates the statistical estimation problem, separate from
optimization strategy.

``` r
problem <- mle_problem(

  loglike,
  score = NULL,           # Auto-computed if NULL
  fisher = NULL,          # Auto-computed if NULL
  constraint = NULL,      # mle_constraint object
  theta_names = NULL,     # Parameter names for nice output
  n_obs = NULL            # For AIC/BIC
)
```

**Key features**: - Lazy numerical differentiation when analytic forms
not provided - Immutable - create new problems via
`update(problem, ...)` - Validates inputs on construction

### 2. Solver Functions

A solver is a function: `(problem, theta0, trace) -> mle_result`

Solver *factories* return solver functions:

``` r
# Factory pattern
gradient_ascent <- function(
  learning_rate = 1.0,
  line_search = TRUE,
  max_iter = 100,
  tol = 1e-8
) {

  # Returns a solver function
  function(problem, theta0, trace = mle_trace()) {
    # ... optimization logic ...
    # Returns mle_result
  }
}
```

This means: -
[`gradient_ascent()`](https://queelius.github.io/compositional.mle/reference/gradient_ascent.md)
returns a solver - `gradient_ascent(max_iter = 200)` returns a
configured solver - All solvers have the same signature:
`(problem, theta0, trace) -> result`

### 3. Composition Operators

**Sequential** (`%>>%`): Chain solvers, passing result as next starting
point

``` r
grid_search(n = 10) %>>% gradient_ascent() %>>% newton_raphson()
```

**Parallel/Race** (`%|%`): Run both, select best

``` r
gradient_ascent() %|% nelder_mead() %|% bfgs()
```

**Restarts**: Multiple starting points

``` r
gradient_ascent() %>% with_restarts(n = 20, sampler = ...)
```

**Conditional**:

``` r
gradient_ascent() %>% unless_converged(newton_raphson())
```

### 4. Tracing System

``` r
trace <- mle_trace(
  values = TRUE,      # Track log-likelihood
  path = TRUE,        # Track parameter values
  gradients = TRUE,   # Track gradient norms
  timing = TRUE       # Track wall-clock time
)

result <- solver(problem, theta0, trace = trace)

# Analyze
plot(result)                    # Convergence plot
optimization_path(result)       # Data frame of path
```

### 5. Results: mle_result

Extends algebraic.mle::mle_numerical with: - `$converged` - logical -
`$iterations` - count - `$trace` - optimization trace (if requested) -
`$chain` - for composed solvers, list of intermediate results -
`$solver` - which solver produced this result

## Solver Inventory

### Gradient-Based (require score)

| Factory                                                                                          | Description                               | Wraps  |
|--------------------------------------------------------------------------------------------------|-------------------------------------------|--------|
| [`gradient_ascent()`](https://queelius.github.io/compositional.mle/reference/gradient_ascent.md) | Steepest ascent with optional line search | native |
| [`bfgs()`](https://queelius.github.io/compositional.mle/reference/bfgs.md)                       | Quasi-Newton BFGS                         | optim  |
| `lbfgs()`                                                                                        | Limited-memory BFGS with box constraints  | optim  |

### Second-Order (require score + fisher)

| Factory                                                                                        | Description            | Wraps  |
|------------------------------------------------------------------------------------------------|------------------------|--------|
| [`newton_raphson()`](https://queelius.github.io/compositional.mle/reference/newton_raphson.md) | Classic Newton-Raphson | native |
| [`fisher_scoring()`](https://queelius.github.io/compositional.mle/reference/fisher_scoring.md) | Uses expected Fisher   | native |

### Derivative-Free

| Factory                                                                                      | Description         | Wraps  |
|----------------------------------------------------------------------------------------------|---------------------|--------|
| [`nelder_mead()`](https://queelius.github.io/compositional.mle/reference/nelder_mead.md)     | Simplex method      | optim  |
| [`grid_search()`](https://queelius.github.io/compositional.mle/reference/grid_search.md)     | Exhaustive grid     | native |
| [`random_search()`](https://queelius.github.io/compositional.mle/reference/random_search.md) | Random sampling     | native |
| `sim_anneal()`                                                                               | Simulated annealing | native |

### Coordinate Methods

| Factory               | Description             | Wraps  |
|-----------------------|-------------------------|--------|
| `coordinate_ascent()` | One parameter at a time | native |

## Example Usage

``` r
library(numerical.mle)

# Generate data
set.seed(42)
data <- rnorm(100, mean = 5, sd = 2)

# Define the problem
problem <- mle_problem(
  loglike = function(theta) {
    sum(dnorm(data, theta[1], theta[2], log = TRUE))
  },
  constraint = mle_constraint(
    support = function(theta) theta[2] > 0,
    project = function(theta) c(theta[1], max(theta[2], 1e-8))
  ),
  theta_names = c("mu", "sigma"),
  n_obs = length(data)
)

# Simple solve
result <- gradient_ascent()(problem, c(0, 1))

# Composed strategy: coarse to fine
strategy <-
  grid_search(lower = c(-10, 0.1), upper = c(10, 5), n = 5) %>>%
  gradient_ascent(max_iter = 50) %>>%
  newton_raphson(max_iter = 20)

result <- strategy(problem, c(0, 1))

# Robust global search
strategy <-
  gradient_ascent() %>%
  with_restarts(n = 10, sampler = function() c(runif(1, -10, 10), runif(1, 0.1, 5)))

result <- strategy(problem, c(0, 1))

# Race different methods
strategy <- gradient_ascent() %|% bfgs() %|% nelder_mead()
result <- strategy(problem, c(0, 1))

# With tracing
result <- gradient_ascent()(
  problem,
  c(0, 1),
  trace = mle_trace(path = TRUE, values = TRUE)
)
plot(result)
```

## File Organization

    R/
      problem.R           # mle_problem(), is_mle_problem(), update.mle_problem()
      solver.R            # Solver protocol, is_solver(), solve()
      compose.R           # %>>%, %|%, with_restarts(), unless_converged()
      trace.R             # mle_trace(), trace methods, plotting
      result.R            # mle_result class, print/summary methods

      solvers/
        gradient.R        # gradient_ascent()
        newton.R          # newton_raphson(), fisher_scoring()
        quasi_newton.R    # bfgs(), lbfgs()
        derivative_free.R # nelder_mead(), grid_search(), random_search()
        annealing.R       # sim_anneal()
        coordinate.R      # coordinate_ascent()

## Open Questions

1.  **Parallel execution**: Should `%|%` actually run in parallel
    (future/parallel package)?

2.  **Automatic differentiation**: Should we support autodiff packages
    for score/fisher?

3.  **Caching**: Should problem cache score/fisher evaluations?

4.  **Verbose output**: How to handle progress reporting during
    optimization?

5.  **Early stopping**: Should composed solvers support early
    termination criteria?
