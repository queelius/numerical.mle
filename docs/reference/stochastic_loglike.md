# stochastic loglikelihood constructor good for large datasets. if applied to a gradient ascent method, this will perform stochastic gradient ascent.

stochastic loglikelihood constructor good for large datasets. if applied
to a gradient ascent method, this will perform stochastic gradient
ascent.

## Usage

``` r
stochastic_loglike(log.p, obs, options)
```

## Arguments

- log.p:

  log pdf (or pmf) of the parametric model being fit to \`obs\`
  parameters. it can also just be proportional to the log pdf, since
  sometimes the normalizing constant is unknown or hard to compute.

- obs:

  a matrix, vector, or data frame of observations

- options:

  a list of options
