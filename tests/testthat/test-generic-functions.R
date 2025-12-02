# Tests for generic functions

test_that("mle_numerical constructor creates valid objects", {
  # Create a simple mle_numerical object
  result <- mle_numerical(
    theta.hat = c(0, 1),
    loglike = -10,
    score = c(0.01, 0.02),
    info = matrix(c(1, 0, 0, 1), nrow = 2),
    sigma = matrix(c(1, 0, 0, 1), nrow = 2),
    iter = 50L,
    converged = TRUE
  )

  expect_true(inherits(result, "mle_numerical"))
  expect_true(inherits(result, "mle"))
  expect_equal(result$iter, 50L)
  expect_true(result$converged)
})

test_that("mle_numerical validates inputs", {
  expect_error(
    mle_numerical(
      theta.hat = c(0, 1),
      loglike = -10,
      score = c(0.01, 0.02),
      info = matrix(c(1, 0, 0, 1), nrow = 2),
      sigma = matrix(c(1, 0, 0, 1), nrow = 2),
      iter = "not a number",
      converged = TRUE
    )
  )

  expect_error(
    mle_numerical(
      theta.hat = c(0, 1),
      loglike = -10,
      score = c(0.01, 0.02),
      info = matrix(c(1, 0, 0, 1), nrow = 2),
      sigma = matrix(c(1, 0, 0, 1), nrow = 2),
      iter = 50L,
      converged = "not logical"
    )
  )
})

test_that("is_converged works correctly", {
  converged_result <- mle_numerical(
    theta.hat = c(0, 1),
    loglike = -10,
    score = c(0.01, 0.02),
    info = matrix(c(1, 0, 0, 1), nrow = 2),
    sigma = matrix(c(1, 0, 0, 1), nrow = 2),
    iter = 50L,
    converged = TRUE
  )

  not_converged_result <- mle_numerical(
    theta.hat = c(0, 1),
    loglike = -10,
    score = c(0.01, 0.02),
    info = matrix(c(1, 0, 0, 1), nrow = 2),
    sigma = matrix(c(1, 0, 0, 1), nrow = 2),
    iter = 100L,
    converged = FALSE
  )

  expect_true(is_converged(converged_result))
  expect_false(is_converged(not_converged_result))
})

test_that("is_mle_numerical works correctly", {
  result <- mle_numerical(
    theta.hat = c(0, 1),
    loglike = -10,
    score = c(0.01, 0.02),
    info = matrix(c(1, 0, 0, 1), nrow = 2),
    sigma = matrix(c(1, 0, 0, 1), nrow = 2),
    iter = 50L,
    converged = TRUE
  )

  expect_true(is_mle_numerical(result))
  expect_false(is_mle_numerical(list(a = 1)))
  expect_false(is_mle_numerical("string"))
  expect_false(is_mle_numerical(123))
  expect_false(is_mle_numerical(NULL))
})

test_that("num_iterations works correctly", {
  result <- mle_numerical(
    theta.hat = c(0, 1),
    loglike = -10,
    score = c(0.01, 0.02),
    info = matrix(c(1, 0, 0, 1), nrow = 2),
    sigma = matrix(c(1, 0, 0, 1), nrow = 2),
    iter = 75L,
    converged = TRUE
  )

  expect_equal(num_iterations(result), 75L)
})

test_that("generic functions work with solver results", {
  # Test with actual solver output
  set.seed(123)
  x <- rnorm(100, mean = 5, sd = 2)

  loglike <- function(theta) {
    if (theta[2] <= 0) return(-Inf)
    sum(dnorm(x, mean = theta[1], sd = theta[2], log = TRUE))
  }

  score <- function(theta) {
    if (theta[2] <= 0) return(c(NA, NA))
    mu <- theta[1]
    sigma <- theta[2]
    n <- length(x)
    d_mu <- sum(x - mu) / sigma^2
    d_sigma <- -n / sigma + sum((x - mu)^2) / sigma^3
    c(d_mu, d_sigma)
  }

  result <- mle_gradient_ascent(
    loglike = loglike,
    score = score,
    theta0 = c(0, 1),
    config = mle_config_linesearch(max_iter = 200)
  )

  expect_true(is_mle_numerical(result))
  expect_true(is.logical(is_converged(result)))
  expect_true(is.numeric(num_iterations(result)))
  expect_true(num_iterations(result) > 0)
})
