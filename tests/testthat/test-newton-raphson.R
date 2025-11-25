## Tests for Newton-Raphson solver

test_that("mle_newton_raphson works on simple quadratic problem", {
  loglike <- function(theta) -(theta[1]^2 + theta[2]^2)
  score <- function(theta) -2 * theta
  fisher <- function(theta) matrix(c(2, 0, 0, 2), nrow = 2)

  result <- mle_newton_raphson(
    loglike = loglike,
    score = score,
    fisher = fisher,
    theta0 = c(2, 2),
    config = mle_config_linesearch(max_iter = 50)
  )

  expect_s3_class(result, "mle_newton_raphson")
  expect_s3_class(result, "mle_numerical")
  expect_true(abs(result$theta.hat[1]) < 0.1)
  expect_true(abs(result$theta.hat[2]) < 0.1)
})

test_that("mle_newton_raphson works with inverted FIM (covariance)", {
  loglike <- function(theta) -(theta[1]^2 + theta[2]^2)
  score <- function(theta) -2 * theta
  covar <- function(theta) matrix(c(0.5, 0, 0, 0.5), nrow = 2)

  result <- mle_newton_raphson(
    loglike = loglike,
    score = score,
    fisher = covar,
    theta0 = c(2, 2),
    config = mle_config_linesearch(max_iter = 50),
    inverted = TRUE
  )

  expect_true(abs(result$theta.hat[1]) < 0.1)
  expect_true(abs(result$theta.hat[2]) < 0.1)
})

test_that("mle_newton_raphson validates inputs", {
  loglike <- function(theta) -sum(theta^2)
  score <- function(theta) -2 * theta
  fisher <- function(theta) diag(2, length(theta))

  expect_error(mle_newton_raphson(
    loglike = "not a function",
    score = score,
    fisher = fisher,
    theta0 = c(1, 1)
  ))

  expect_error(mle_newton_raphson(
    loglike = loglike,
    score = "not a function",
    fisher = fisher,
    theta0 = c(1, 1)
  ))

  expect_error(mle_newton_raphson(
    loglike = loglike,
    score = score,
    fisher = "not a function",
    theta0 = c(1, 1)
  ))
})

test_that("mle_newton_raphson respects max_iter", {
  loglike <- function(theta) -(theta[1]^2 + theta[2]^2)
  score <- function(theta) -2 * theta
  fisher <- function(theta) matrix(c(0.001, 0, 0, 0.001), nrow = 2)  # Very small FIM

  result <- mle_newton_raphson(
    loglike = loglike,
    score = score,
    fisher = fisher,
    theta0 = c(100, 100),
    config = mle_config_linesearch(max_iter = 3)
  )

  expect_equal(result$iter, 3)
})

test_that("mle_newton_raphson works with constrained optimization", {
  loglike <- function(theta) -(theta[1]^2 + theta[2]^2)
  score <- function(theta) -2 * theta
  fisher <- function(theta) matrix(c(2, 0, 0, 2), nrow = 2)

  constraint <- mle_constraint(
    support = function(theta) all(theta > 0),
    project = function(theta) pmax(theta, 0.01)
  )

  result <- mle_newton_raphson(
    loglike = loglike,
    score = score,
    fisher = fisher,
    theta0 = c(2, 2),
    config = mle_config_linesearch(max_iter = 50),
    constraint = constraint
  )

  expect_true(all(result$theta.hat > 0))
})

test_that("mle_nr convenience wrapper works", {
  loglike <- function(theta) -(theta[1]^2 + theta[2]^2)
  score <- function(theta) -2 * theta
  fisher <- function(theta) matrix(c(2, 0, 0, 2), nrow = 2)

  result <- mle_nr(loglike, score, fisher, theta0 = c(3, 3), max_iter = 30)

  expect_s3_class(result, "mle_newton_raphson")
  expect_true(abs(result$theta.hat[1]) < 0.5)
  expect_true(abs(result$theta.hat[2]) < 0.5)
})

test_that("mle_newton_raphson works on 1D problems", {
  loglike <- function(theta) -(theta - 5)^2
  score <- function(theta) -2 * (theta - 5)
  fisher <- function(theta) matrix(2, nrow = 1)

  result <- mle_newton_raphson(
    loglike = loglike,
    score = score,
    fisher = fisher,
    theta0 = 0,
    config = mle_config_linesearch(max_iter = 50)
  )

  expect_true(abs(result$theta.hat - 5) < 0.5)
})
