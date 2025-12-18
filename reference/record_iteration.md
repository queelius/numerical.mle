# Record an iteration to trace

Record an iteration to trace

## Usage

``` r
record_iteration(recorder, theta, value = NULL, gradient = NULL)
```

## Arguments

- recorder:

  Trace recorder from new_trace_recorder

- theta:

  Current parameters

- value:

  Current log-likelihood (or NULL)

- gradient:

  Current gradient (or NULL)
