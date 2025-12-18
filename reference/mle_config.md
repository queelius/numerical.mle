# Create optimization configuration

Creates a base configuration object for MLE optimization algorithms.
This object stores convergence criteria and debugging options.

## Usage

``` r
mle_config(
  max_iter = 100L,
  abs_tol = NULL,
  rel_tol = 1e-05,
  trace = FALSE,
  debug = FALSE,
  debug_freq = 1L
)
```

## Arguments

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

An mle_config object

## Examples

``` r
# Basic configuration
config <- mle_config(max_iter = 200, rel_tol = 1e-6)

# Configuration with tracing
config <- mle_config(trace = TRUE, debug = TRUE)
```
