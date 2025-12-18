# Create a Trace Configuration

Specifies what information to track during optimization.

## Usage

``` r
mle_trace(
  values = FALSE,
  path = FALSE,
  gradients = FALSE,
  timing = FALSE,
  every = 1L
)
```

## Arguments

- values:

  Track log-likelihood values at each iteration

- path:

  Track parameter values at each iteration

- gradients:

  Track gradient norms at each iteration

- timing:

  Track wall-clock time

- every:

  Record every nth iteration (1 = all iterations)

## Value

An mle_trace configuration object

## Examples

``` r
# Track everything
trace <- mle_trace(values = TRUE, path = TRUE, gradients = TRUE)

# Minimal tracing (just convergence path)
trace <- mle_trace(values = TRUE)

# Sample every 10th iteration for long runs
trace <- mle_trace(values = TRUE, path = TRUE, every = 10)
```
