## Tests for Newton-Raphson solver

test_that("newton_raphson works on simple quadratic problem", {
  problem <- mle_problem(
    loglike = function(theta) -(theta[1]^2 + theta[2]^2),
    score = function(theta) -2 * theta,
    fisher = function(theta) matrix(c(2, 0, 0, 2), nrow = 2)
  )

  solver <- newton_raphson(max_iter = 50, line_search = TRUE)
  result <- solver(problem, c(2, 2))

  expect_s3_class(result, "mle_newton_raphson")
  expect_s3_class(result, "mle_numerical")
  expect_true(abs(result$theta.hat[1]) < 0.1)
  expect_true(abs(result$theta.hat[2]) < 0.1)
})

test_that("newton_raphson respects max_iter", {
  problem <- mle_problem(
    loglike = function(theta) -(theta[1]^2 + theta[2]^2),
    score = function(theta) -2 * theta,
    fisher = function(theta) matrix(c(0.001, 0, 0, 0.001), nrow = 2)  # Very small FIM
  )

  solver <- newton_raphson(max_iter = 3)
  result <- solver(problem, c(100, 100))

  expect_equal(result$iterations, 3)
})

test_that("newton_raphson works with constrained optimization", {
  constraint <- mle_constraint(
    support = function(theta) all(theta > 0),
    project = function(theta) pmax(theta, 0.01)
  )

  problem <- mle_problem(
    loglike = function(theta) -(theta[1]^2 + theta[2]^2),
    score = function(theta) -2 * theta,
    fisher = function(theta) matrix(c(2, 0, 0, 2), nrow = 2),
    constraint = constraint
  )

  solver <- newton_raphson(max_iter = 50, line_search = TRUE)
  result <- solver(problem, c(2, 2))

  expect_true(all(result$theta.hat > 0))
})

test_that("newton_raphson works on 1D problems", {
  problem <- mle_problem(
    loglike = function(theta) -(theta - 5)^2,
    score = function(theta) -2 * (theta - 5),
    fisher = function(theta) matrix(2, nrow = 1)
  )

  solver <- newton_raphson(max_iter = 50, line_search = TRUE)
  result <- solver(problem, 0)

  expect_true(abs(result$theta.hat - 5) < 0.5)
})

test_that("newton_raphson without line_search", {
  problem <- mle_problem(
    loglike = function(theta) -(theta[1]^2 + theta[2]^2),
    score = function(theta) -2 * theta,
    fisher = function(theta) matrix(c(2, 0, 0, 2), nrow = 2)
  )

  solver <- newton_raphson(max_iter = 50, line_search = FALSE)
  result <- solver(problem, c(1, 1))

  # Should still converge on this well-behaved problem
  expect_true(abs(result$theta.hat[1]) < 0.1)
  expect_true(abs(result$theta.hat[2]) < 0.1)
})

test_that("fisher_scoring works on normal MLE", {
  set.seed(123)
  data <- rnorm(100, mean = 5, sd = 2)

  problem <- mle_problem(
    loglike = function(theta) {
      if (theta[2] <= 0) return(-Inf)
      sum(dnorm(data, theta[1], theta[2], log = TRUE))
    },
    score = function(theta) {
      mu <- theta[1]
      sigma <- theta[2]
      n <- length(data)
      c(sum(data - mu) / sigma^2,
        -n / sigma + sum((data - mu)^2) / sigma^3)
    },
    fisher = function(theta) {
      sigma <- theta[2]
      n <- length(data)
      matrix(c(n / sigma^2, 0, 0, 2 * n / sigma^2), nrow = 2)
    },
    constraint = mle_constraint(
      support = function(theta) theta[2] > 0,
      project = function(theta) c(theta[1], max(theta[2], 1e-6))
    )
  )

  solver <- fisher_scoring(max_iter = 50)
  result <- solver(problem, c(0, 1))

  # fisher_scoring is an alias for newton_raphson
  expect_s3_class(result, "mle_newton_raphson")
  expect_true(abs(result$theta.hat[1] - mean(data)) < 0.5)
  expect_true(abs(result$theta.hat[2] - sd(data)) < 0.5)
})
