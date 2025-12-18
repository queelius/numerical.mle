# Create gradient-based optimization configuration

Extends base configuration with gradient-specific parameters like
learning rate and distance metric.

## Usage

``` r
mle_config_gradient(
  eta = 1,
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

- eta:

  Learning rate / step size (numeric, default: 1.0)

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

An mle_config_gradient object

## Examples

``` r
# Basic gradient configuration
config <- mle_config_gradient(eta = 0.1, max_iter = 500)

# With custom norm (L2 norm)
config <- mle_config_gradient(
  eta = 0.01,
  norm = function(x) sqrt(sum(x^2))
)
```
