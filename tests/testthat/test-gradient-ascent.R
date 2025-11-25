## Tests for gradient ascent solver

test_that("mle_gradient_ascent works on simple quadratic problem", {
  loglike <- function(theta) -(theta[1]^2 + theta[2]^2)
  score <- function(theta) -2 * theta

  result <- mle_gradient_ascent(
    loglike = loglike,
    score = score,
    theta0 = c(2, 2),
    config = mle_config_linesearch(max_iter = 100)
  )

  expect_s3_class(result, "mle_gradient_ascent")
  expect_s3_class(result, "mle_numerical")
  expect_true(abs(result$theta.hat[1]) < 0.1)
  expect_true(abs(result$theta.hat[2]) < 0.1)
})

test_that("mle_gradient_ascent respects max_iter", {
  loglike <- function(theta) -(theta[1]^2 + theta[2]^2)
  score <- function(theta) -2 * theta

  result <- mle_gradient_ascent(
    loglike = loglike,
    score = score,
    theta0 = c(10, 10),
    config = mle_config_gradient(eta = 0.001, max_iter = 5)
  )

  expect_equal(result$iter, 5)
  expect_false(result$converged)
})

test_that("mle_gradient_ascent works with constrained optimization", {
  loglike <- function(theta) -(theta[1]^2 + theta[2]^2)
  score <- function(theta) -2 * theta

  constraint <- mle_constraint(
    support = function(theta) all(theta > 0),
    project = function(theta) pmax(theta, 0.01)
  )

  result <- mle_gradient_ascent(
    loglike = loglike,
    score = score,
    theta0 = c(2, 2),
    config = mle_config_linesearch(max_iter = 100),
    constraint = constraint
  )

  # Should converge near boundary
  expect_true(all(result$theta.hat > 0))
  expect_true(all(result$theta.hat < 0.2))
})

test_that("mle_gradient_ascent validates inputs", {
  loglike <- function(theta) -sum(theta^2)
  score <- function(theta) -2 * theta

  expect_error(mle_gradient_ascent(
    loglike = "not a function",
    score = score,
    theta0 = c(1, 1)
  ))

  expect_error(mle_gradient_ascent(
    loglike = loglike,
    score = "not a function",
    theta0 = c(1, 1)
  ))
})

test_that("mle_gradient_ascent stores path when trace=TRUE", {
  loglike <- function(theta) -sum(theta^2)
  score <- function(theta) -2 * theta

  config <- mle_config_linesearch(max_iter = 20, trace = TRUE)

  result <- mle_gradient_ascent(
    loglike = loglike,
    score = score,
    theta0 = c(5, 5),
    config = config
  )

  expect_true(!is.null(result$path))
  expect_equal(ncol(result$path), 2)
  expect_true(nrow(result$path) > 0)
})

test_that("mle_grad convenience wrapper works", {
  loglike <- function(theta) -(theta[1]^2 + theta[2]^2)
  score <- function(theta) -2 * theta

  result <- mle_grad(loglike, score, theta0 = c(3, 3), max_iter = 50)

  expect_s3_class(result, "mle_gradient_ascent")
  expect_true(abs(result$theta.hat[1]) < 0.5)
  expect_true(abs(result$theta.hat[2]) < 0.5)
})

test_that("mle_gradient_ascent works on 1D problems", {
  loglike <- function(theta) -(theta - 5)^2
  score <- function(theta) -2 * (theta - 5)

  result <- mle_gradient_ascent(
    loglike = loglike,
    score = score,
    theta0 = 0,
    config = mle_config_linesearch(max_iter = 100)
  )

  expect_true(abs(result$theta.hat - 5) < 0.5)
})
