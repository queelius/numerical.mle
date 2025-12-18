# BFGS Solver

Creates a solver using the BFGS quasi-Newton method via
[`optim()`](https://rdrr.io/r/stats/optim.html). BFGS approximates the
Hessian from gradient information, providing second-order-like
convergence without computing the Hessian directly.

## Usage

``` r
bfgs(max_iter = 100L, tol = 1e-08, report = 0L)
```

## Arguments

- max_iter:

  Maximum number of iterations

- tol:

  Convergence tolerance (passed to optim's reltol)

- report:

  Reporting frequency (0 = no reporting)

## Value

A solver function with signature (problem, theta0, trace) -\> mle_result

## Details

BFGS is often a good default choice: it's more robust than
Newton-Raphson (no matrix inversion issues) and faster than gradient
ascent (uses curvature information).

The solver automatically uses the score function from the problem if
available, otherwise computes gradients numerically.

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic usage
result <- bfgs()(problem, c(0, 1))

# Race BFGS against gradient ascent
strategy <- bfgs() %|% gradient_ascent()
} # }
```
