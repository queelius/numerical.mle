# Create domain constraint specification

Specifies domain constraints for optimization. The support function
checks if parameters are valid, and the project function maps invalid
parameters back to valid ones.

## Usage

``` r
mle_constraint(support = function(theta) TRUE, project = function(theta) theta)
```

## Arguments

- support:

  Function testing if theta is in support (returns TRUE/FALSE)

- project:

  Function projecting theta onto support

## Value

An mle_constraint object

## Examples

``` r
# Positive parameters only
constraint <- mle_constraint(
  support = function(theta) all(theta > 0),
  project = function(theta) pmax(theta, 1e-8)
)

# Parameters in [0, 1]
constraint <- mle_constraint(
  support = function(theta) all(theta >= 0 & theta <= 1),
  project = function(theta) pmax(0, pmin(1, theta))
)

# No constraints (default)
constraint <- mle_constraint()
```
