# Create line search configuration

Extends gradient configuration with backtracking line search parameters.
Line search adaptively finds step sizes that improve the objective.

## Usage

``` r
mle_config_linesearch(
  max_step = 1,
  backtrack_ratio = 0.5,
  max_iter_ls = 10L,
  min_step = 1e-08,
  norm = function(x) max(abs(x)),
  max_iter = 100L,
  abs_tol = NULL,
  rel_tol = 1e-05,
  trace = FALSE,
  debug = FALSE,
  debug_freq = 1L
)
```

## Arguments

- max_step:

  Maximum step size per iteration (numeric, default: 1.0)

- backtrack_ratio:

  Backtracking multiplier, must be in (0,1) (default: 0.5)

- max_iter_ls:

  Maximum line search iterations (integer, default: 10)

- min_step:

  Minimum step size threshold (numeric, default: 1e-8)

- norm:

  Distance measure function (default: max absolute value)

- max_iter:

  Maximum iterations (integer, default: 100)

- abs_tol:

  Absolute tolerance (numeric or NULL, default: NULL to use rel_tol)

- rel_tol:

  Relative tolerance (numeric, default: 1e-5)

- trace:

  Store optimization path (logical, default: FALSE)

- debug:

  Print debug information (logical, default: FALSE)

- debug_freq:

  Debug output frequency (integer, default: 1)

## Value

An mle_config_linesearch object

## Examples

``` r
# Conservative line search
config <- mle_config_linesearch(
  max_step = 0.5,
  backtrack_ratio = 0.8
)

# Aggressive line search
config <- mle_config_linesearch(
  max_step = 10.0,
  backtrack_ratio = 0.3,
  max_iter_ls = 20
)
```
