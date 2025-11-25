test_that("mle_config creates valid configuration objects", {
  # Default configuration
  config <- mle_config()
  expect_s3_class(config, "mle_config")
  expect_equal(config$max_iter, 100)
  expect_null(config$abs_tol)
  expect_equal(config$rel_tol, 1e-5)
  expect_false(config$trace)
  expect_false(config$debug)
  expect_equal(config$debug_freq, 1)

  # Custom configuration
  config <- mle_config(
    max_iter = 500,
    abs_tol = 1e-8,
    rel_tol = 1e-6,
    trace = TRUE,
    debug = TRUE,
    debug_freq = 10
  )
  expect_equal(config$max_iter, 500)
  expect_equal(config$abs_tol, 1e-8)
  expect_equal(config$rel_tol, 1e-6)
  expect_true(config$trace)
  expect_true(config$debug)
  expect_equal(config$debug_freq, 10)
})

test_that("mle_config validates inputs", {
  # Invalid max_iter
  expect_error(mle_config(max_iter = -10))
  expect_error(mle_config(max_iter = 0))

  # Invalid tolerances
  expect_error(mle_config(abs_tol = -1e-5))
  expect_error(mle_config(rel_tol = -1e-5))
  expect_error(mle_config(rel_tol = 0))

  # Invalid logical flags
  expect_error(mle_config(trace = "yes"))
  expect_error(mle_config(debug = 1))

  # Invalid debug_freq
  expect_error(mle_config(debug_freq = 0))
  expect_error(mle_config(debug_freq = -1))
})

test_that("mle_config coerces numeric to integer", {
  config <- mle_config(max_iter = 200.5, debug_freq = 5.5)
  expect_type(config$max_iter, "integer")
  expect_equal(config$max_iter, 200)
  expect_type(config$debug_freq, "integer")
  expect_equal(config$debug_freq, 5)
})

test_that("mle_config_gradient extends mle_config", {
  config <- mle_config_gradient()
  expect_s3_class(config, c("mle_config_gradient", "mle_config"))
  expect_equal(config$eta, 1.0)
  expect_true(is.function(config$norm))

  # Custom configuration
  custom_norm <- function(x) sqrt(sum(x^2))
  config <- mle_config_gradient(
    eta = 0.1,
    norm = custom_norm,
    max_iter = 200
  )
  expect_equal(config$eta, 0.1)
  expect_identical(config$norm, custom_norm)
  expect_equal(config$max_iter, 200)
})

test_that("mle_config_gradient validates inputs", {
  # Invalid eta
  expect_error(mle_config_gradient(eta = 0))
  expect_error(mle_config_gradient(eta = -0.1))

  # Invalid norm
  expect_error(mle_config_gradient(norm = "not a function"))
  expect_error(mle_config_gradient(norm = 42))
})

test_that("mle_config_gradient norm function works", {
  config <- mle_config_gradient()

  # Default is max absolute value
  expect_equal(config$norm(c(-3, 2, 1)), 3)
  expect_equal(config$norm(c(1, 5, -2)), 5)

  # Custom L2 norm
  l2_norm <- function(x) sqrt(sum(x^2))
  config <- mle_config_gradient(norm = l2_norm)
  expect_equal(config$norm(c(3, 4)), 5)
})

test_that("mle_config_linesearch extends mle_config_gradient", {
  config <- mle_config_linesearch()
  expect_s3_class(config, c("mle_config_linesearch", "mle_config_gradient", "mle_config"))
  expect_equal(config$eta, 1.0)  # max_step becomes eta
  expect_equal(config$backtrack_ratio, 0.5)
  expect_equal(config$max_iter_ls, 10)
  expect_equal(config$min_step, 1e-8)

  # Custom configuration
  config <- mle_config_linesearch(
    max_step = 2.0,
    backtrack_ratio = 0.3,
    max_iter_ls = 20,
    min_step = 1e-10
  )
  expect_equal(config$eta, 2.0)
  expect_equal(config$backtrack_ratio, 0.3)
  expect_equal(config$max_iter_ls, 20)
  expect_equal(config$min_step, 1e-10)
})

test_that("mle_config_linesearch validates inputs", {
  # Invalid backtrack_ratio
  expect_error(mle_config_linesearch(backtrack_ratio = 0))
  expect_error(mle_config_linesearch(backtrack_ratio = 1))
  expect_error(mle_config_linesearch(backtrack_ratio = 1.5))
  expect_error(mle_config_linesearch(backtrack_ratio = -0.1))

  # Invalid max_iter_ls
  expect_error(mle_config_linesearch(max_iter_ls = 0))
  expect_error(mle_config_linesearch(max_iter_ls = -5))

  # Invalid min_step
  expect_error(mle_config_linesearch(min_step = 0))
  expect_error(mle_config_linesearch(min_step = -1e-8))
})

test_that("mle_constraint creates valid constraint objects", {
  # Default (unconstrained)
  constraint <- mle_constraint()
  expect_s3_class(constraint, "mle_constraint")
  expect_true(is.function(constraint$support))
  expect_true(is.function(constraint$project))

  # Support should accept anything by default
  expect_true(constraint$support(c(-10, 20, 30)))
  expect_true(constraint$support(c(0, 0, 0)))

  # Project should be identity by default
  theta <- c(1, -2, 3)
  expect_equal(constraint$project(theta), theta)
})

test_that("mle_constraint with custom functions works", {
  # Positive constraint
  constraint <- mle_constraint(
    support = function(theta) all(theta > 0),
    project = function(theta) pmax(theta, 1e-8)
  )

  # Test support
  expect_true(constraint$support(c(1, 2, 3)))
  expect_false(constraint$support(c(1, -2, 3)))
  expect_false(constraint$support(c(0, 1, 2)))

  # Test projection
  expect_equal(constraint$project(c(1, -2, 3)), c(1, 1e-8, 3))
  expect_equal(constraint$project(c(-1, -2, -3)), c(1e-8, 1e-8, 1e-8))
})

test_that("mle_constraint validates inputs", {
  expect_error(mle_constraint(support = "not a function"))
  expect_error(mle_constraint(project = 42))
  expect_error(mle_constraint(support = NULL))
})

test_that("mle_constraint with box constraints works", {
  # Box constraints: 0 <= theta <= 1
  constraint <- mle_constraint(
    support = function(theta) all(theta >= 0 & theta <= 1),
    project = function(theta) pmax(0, pmin(1, theta))
  )

  # Test support
  expect_true(constraint$support(c(0.5, 0.5)))
  expect_true(constraint$support(c(0, 1)))
  expect_false(constraint$support(c(-0.1, 0.5)))
  expect_false(constraint$support(c(0.5, 1.1)))

  # Test projection
  expect_equal(constraint$project(c(-0.5, 0.5, 1.5)), c(0, 0.5, 1))
  expect_equal(constraint$project(c(0.3, 2.0)), c(0.3, 1.0))
})

test_that("is_mle_config works correctly", {
  expect_true(is_mle_config(mle_config()))
  expect_true(is_mle_config(mle_config_gradient()))
  expect_true(is_mle_config(mle_config_linesearch()))

  expect_false(is_mle_config(list(max_iter = 100)))
  expect_false(is_mle_config("not a config"))
  expect_false(is_mle_config(NULL))
})

test_that("is_mle_constraint works correctly", {
  expect_true(is_mle_constraint(mle_constraint()))
  expect_false(is_mle_constraint(list(support = function(x) TRUE)))
  expect_false(is_mle_constraint("not a constraint"))
  expect_false(is_mle_constraint(NULL))
})

test_that("Configuration objects are reusable", {
  # Create a config and use it multiple times
  config <- mle_config_linesearch(max_iter = 500, max_step = 0.5)

  # Should be able to use this config for different solvers
  expect_equal(config$max_iter, 500)
  expect_equal(config$eta, 0.5)

  # Config should not be modified by use
  config_copy <- config
  expect_identical(config, config_copy)
})

test_that("Constraint objects are reusable", {
  # Create a constraint and use it multiple times
  constraint <- mle_constraint(
    support = function(theta) all(theta > 0),
    project = function(theta) pmax(theta, 0.01)
  )

  # Should work for different parameter vectors
  expect_true(constraint$support(c(1, 2)))
  expect_false(constraint$support(c(-1, 2)))
  expect_equal(constraint$project(c(-1, 2)), c(0.01, 2))
  expect_equal(constraint$project(c(0, 3)), c(0.01, 3))
})
