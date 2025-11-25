# num_iterations

num_iterations

## Usage

``` r
num_iterations(x, ...)
```

## Arguments

- x:

  the \`mle\` object

- ...:

  additional arguments to pass

## Value

the number of iterations used to find the MLE

## Examples

``` r
loglike <- function(theta) {
   -sum(dnorm(theta, mean = theta[1], sd = theta[2], log = TRUE))
}
score <- function(theta) {
  -numDeriv::grad(loglike, theta)
}
sol <- mle_gradient_raphson(theta0 = theta, score = score, loglike = loglike)
#> Error in mle_gradient_raphson(theta0 = theta, score = score, loglike = loglike): could not find function "mle_gradient_raphson"
num_iterations(sol)
#> Error: object 'sol' not found
```
