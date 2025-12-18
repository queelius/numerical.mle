# Nelder-Mead Solver (Derivative-Free)

Creates a solver using the Nelder-Mead simplex method via
[`optim()`](https://rdrr.io/r/stats/optim.html). This is a
derivative-free method useful when gradients are unavailable or
unreliable.

## Usage

``` r
nelder_mead(max_iter = 500L, tol = 1e-08)
```

## Arguments

- max_iter:

  Maximum number of iterations

- tol:

  Convergence tolerance

## Value

A solver function

## Details

Nelder-Mead doesn't use gradient information, making it robust but
potentially slower. It's useful as a fallback when gradient-based
methods fail, or for problems with non-smooth likelihoods.

## Examples

``` r
if (FALSE) { # \dontrun{
# Use when gradients are problematic
result <- nelder_mead()(problem, c(0, 1))

# Race against gradient methods
strategy <- gradient_ascent() %|% nelder_mead()
} # }
```
