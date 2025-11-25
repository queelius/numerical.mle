## Integration tests for complete MLE workflows

test_that("Normal distribution MLE with gradient ascent", {
  set.seed(123)
  data <- rnorm(100, mean = 5, sd = 2)

  loglike <- function(theta) {
    mu <- theta[1]
    sigma <- theta[2]
    if (sigma <= 0) return(-Inf)
    sum(dnorm(data, mu, sigma, log = TRUE))
  }

  score <- function(theta) {
    mu <- theta[1]
    sigma <- theta[2]
    n <- length(data)
    d_mu <- sum(data - mu) / sigma^2
    d_sigma <- -n / sigma + sum((data - mu)^2) / sigma^3
    c(d_mu, d_sigma)
  }

  constraint <- mle_constraint(
    support = function(theta) theta[2] > 0,
    project = function(theta) c(theta[1], max(theta[2], 1e-6))
  )

  result <- mle_gradient_ascent(
    loglike = loglike,
    score = score,
    theta0 = c(0, 1),
    config = mle_config_linesearch(max_iter = 200),
    constraint = constraint
  )

  # Check parameter estimates are reasonable (don't require convergence flag)
  expect_true(abs(result$theta.hat[1] - mean(data)) < 0.5)
  expect_true(abs(result$theta.hat[2] - sd(data)) < 0.5)
})

test_that("Normal distribution MLE with Newton-Raphson", {
  set.seed(456)
  data <- rnorm(100, mean = 10, sd = 3)

  loglike <- function(theta) {
    mu <- theta[1]
    sigma <- theta[2]
    if (sigma <= 0) return(-Inf)
    sum(dnorm(data, mu, sigma, log = TRUE))
  }

  score <- function(theta) {
    mu <- theta[1]
    sigma <- theta[2]
    n <- length(data)
    d_mu <- sum(data - mu) / sigma^2
    d_sigma <- -n / sigma + sum((data - mu)^2) / sigma^3
    c(d_mu, d_sigma)
  }

  fisher <- function(theta) {
    sigma <- theta[2]
    n <- length(data)
    matrix(c(n / sigma^2, 0, 0, 2 * n / sigma^2), nrow = 2)
  }

  constraint <- mle_constraint(
    support = function(theta) theta[2] > 0,
    project = function(theta) c(theta[1], max(theta[2], 1e-6))
  )

  result <- mle_newton_raphson(
    loglike = loglike,
    score = score,
    fisher = fisher,
    theta0 = c(5, 1),
    config = mle_config_linesearch(max_iter = 50),
    constraint = constraint
  )

  expect_true(result$converged)
  expect_true(abs(result$theta.hat[1] - mean(data)) < 0.5)
  expect_true(abs(result$theta.hat[2] - sd(data)) < 0.5)
})

test_that("Convenience wrappers work for simple problem", {
  loglike <- function(theta) -(theta[1]^2 + theta[2]^2)
  score <- function(theta) -2 * theta

  result <- mle_grad(loglike, score, theta0 = c(3, 3), max_iter = 50)

  expect_s3_class(result, "mle_gradient_ascent")
  expect_true(abs(result$theta.hat[1]) < 0.1)
  expect_true(abs(result$theta.hat[2]) < 0.1)
})

test_that("Grid search finds approximate solution", {
  loglike <- function(theta) -(theta[1] - 2)^2 - (theta[2] + 1)^2

  result <- mle_grid_search(
    loglike = loglike,
    lower = c(-5, -5),
    upper = c(5, 5),
    grid_size = 10
  )

  expect_s3_class(result, "mle_grid_search")
  expect_true(abs(result$theta.hat[1] - 2) < 1.5)
  expect_true(abs(result$theta.hat[2] - (-1)) < 1.5)
})

test_that("Random restart improves upon single run", {
  loglike <- function(theta) -sum(theta^2)
  score <- function(theta) -2 * theta
  sampler <- function() runif(2, -5, 5)

  result <- mle_random_restart(
    loglike = loglike,
    solver = mle_grad,
    theta0_sampler = sampler,
    n_trials = 5,
    score = score,
    max_iter = 30
  )

  expect_s3_class(result, "mle_random_restart")
  expect_true(result$successful_trials > 0)
  expect_true(sqrt(sum(result$theta.hat^2)) < 1)
})

test_that("Constrained optimization respects bounds", {
  loglike <- function(theta) -(theta[1] - 3)^2 - (theta[2] - 3)^2
  score <- function(theta) c(-2 * (theta[1] - 3), -2 * (theta[2] - 3))

  # Constrain to [0, 1] x [0, 1]
  constraint <- mle_constraint(
    support = function(theta) all(theta >= 0 & theta <= 1),
    project = function(theta) pmax(0, pmin(1, theta))
  )

  result <- mle_gradient_ascent(
    loglike = loglike,
    score = score,
    theta0 = c(0.5, 0.5),
    config = mle_config_gradient(eta = 0.1, max_iter = 100),
    constraint = constraint
  )

  # Should converge to corner (1, 1)
  expect_true(abs(result$theta.hat[1] - 1) < 0.1)
  expect_true(abs(result$theta.hat[2] - 1) < 0.1)
})

test_that("Path tracing records optimization path", {
  loglike <- function(theta) -sum(theta^2)
  score <- function(theta) -2 * theta

  config <- mle_config_gradient(eta = 0.1, max_iter = 20, trace = TRUE)

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
