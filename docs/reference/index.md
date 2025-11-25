# Package index

## Configuration

Create configuration objects for optimization algorithms

- [`mle_config()`](https://queelius.github.io/numerical.mle/reference/mle_config.md)
  : Create optimization configuration
- [`mle_config_gradient()`](https://queelius.github.io/numerical.mle/reference/mle_config_gradient.md)
  : Create gradient-based optimization configuration
- [`mle_config_linesearch()`](https://queelius.github.io/numerical.mle/reference/mle_config_linesearch.md)
  : Create line search configuration
- [`mle_constraint()`](https://queelius.github.io/numerical.mle/reference/mle_constraint.md)
  : Create domain constraint specification
- [`is_mle_config()`](https://queelius.github.io/numerical.mle/reference/is_mle_config.md)
  : Check if object is an mle_config
- [`is_mle_constraint()`](https://queelius.github.io/numerical.mle/reference/is_mle_constraint.md)
  : Check if object is an mle_constraint

## Core Solvers

Main optimization algorithms

- [`mle_gradient_ascent()`](https://queelius.github.io/numerical.mle/reference/mle_gradient_ascent.md)
  : Maximum likelihood estimation via gradient ascent
- [`mle_newton_raphson()`](https://queelius.github.io/numerical.mle/reference/mle_newton_raphson.md)
  : Maximum likelihood estimation via Newton-Raphson

## Meta-Solvers

Global optimization strategies

- [`mle_grid_search()`](https://queelius.github.io/numerical.mle/reference/mle_grid_search.md)
  : MLE via grid search
- [`mle_random_restart()`](https://queelius.github.io/numerical.mle/reference/mle_random_restart.md)
  : MLE via random restarts

## Convenience Wrappers

Quick access to solvers with sensible defaults

- [`mle_grad()`](https://queelius.github.io/numerical.mle/reference/mle_grad.md)
  : Quick gradient ascent with sensible defaults
- [`mle_nr()`](https://queelius.github.io/numerical.mle/reference/mle_nr.md)
  : Quick Newton-Raphson with sensible defaults
- [`with_constraint()`](https://queelius.github.io/numerical.mle/reference/with_constraint.md)
  : Quick constrained optimization

## Function Transformers

Transform log-likelihood functions

- [`with_subsampling()`](https://queelius.github.io/numerical.mle/reference/with_subsampling.md)
  : Create stochastic log-likelihood with subsampling
- [`with_penalty()`](https://queelius.github.io/numerical.mle/reference/with_penalty.md)
  : Add penalty term to log-likelihood
- [`penalty_l1()`](https://queelius.github.io/numerical.mle/reference/penalty_l1.md)
  : L1 penalty function (LASSO)
- [`penalty_l2()`](https://queelius.github.io/numerical.mle/reference/penalty_l2.md)
  : L2 penalty function (Ridge)
- [`penalty_elastic_net()`](https://queelius.github.io/numerical.mle/reference/penalty_elastic_net.md)
  : Elastic net penalty (combination of L1 and L2)
- [`compose()`](https://queelius.github.io/numerical.mle/reference/compose.md)
  : Compose multiple function transformations

## Generic Methods

Methods for mle_numerical objects

- [`is_converged()`](https://queelius.github.io/numerical.mle/reference/is_converged.md)
  : is_converged
- [`is_mle_numerical()`](https://queelius.github.io/numerical.mle/reference/is_mle_numerical.md)
  : is_mle_numerical
- [`num_iterations()`](https://queelius.github.io/numerical.mle/reference/num_iterations.md)
  : num_iterations
- [`mle_numerical()`](https://queelius.github.io/numerical.mle/reference/mle_numerical.md)
  : mle_numerical

## Package

Package documentation

- [`numerical.mle-package`](https://queelius.github.io/numerical.mle/reference/numerical.mle.md)
  [`numerical.mle`](https://queelius.github.io/numerical.mle/reference/numerical.mle.md)
  : \`numerical.mle\`: A package for numerically solving maximum
  likelihood estimators from log-likelihood functions.
