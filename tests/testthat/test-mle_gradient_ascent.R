test_that("mle_gradient_ascent finds normal distribution mean", {
  # MLE for normal distribution mean (variance known)
  set.seed(123)
  true_mean <- 5
  true_sd <- 2
  n <- 100
  data <- rnorm(n, mean = true_mean, sd = true_sd)

  # Log-likelihood
  loglike <- function(mu) {
    sum(dnorm(data, mean = mu, sd = true_sd, log = TRUE))
  }

  # Score function (gradient of log-likelihood)
  score <- function(mu) {
    sum((data - mu) / true_sd^2)
  }

  result <- mle_gradient_ascent(
    theta0 = 0,
    score = score,
    options = list(
      loglike = loglike,
      line_search = TRUE,
      max_iter = 100,
      rel_tol = 1e-5
    )
  )

  expect_s3_class(result, "mle_gradient_ascent")
  expect_s3_class(result, "mle_local_search")
  expect_s3_class(result, "mle_numerical")
  expect_true(result$converged)
  expect_equal(result$theta.hat, mean(data), tolerance = 1e-3)

  # Should have Fisher information and covariance
  expect_true(!is.null(result$info))
  expect_true(!is.null(result$sigma))
  expect_true(!is.null(result$score))
})

test_that("mle_gradient_ascent validates score function", {
  # score must be a function
  expect_error(
    mle_gradient_ascent(
      theta0 = 0,
      score = "not_a_function",
      options = list()
    ),
    "score must be a function"
  )
})

test_that("mle_gradient_ascent works with multidimensional parameters", {
  # 2D normal MLE for mean vector (covariance known)
  set.seed(456)
  n <- 100
  true_mu <- c(3, -1)
  data <- cbind(rnorm(n, true_mu[1], 1), rnorm(n, true_mu[2], 1))

  loglike <- function(mu) {
    sum(dnorm(data[,1], mean = mu[1], sd = 1, log = TRUE)) +
    sum(dnorm(data[,2], mean = mu[2], sd = 1, log = TRUE))
  }

  score <- function(mu) {
    c(sum(data[,1] - mu[1]),
      sum(data[,2] - mu[2]))
  }

  result <- mle_gradient_ascent(
    theta0 = c(0, 0),
    score = score,
    options = list(
      loglike = loglike,
      line_search = TRUE,
      max_iter = 100,
      rel_tol = 1e-5
    )
  )

  expect_true(result$converged)
  expect_length(result$theta.hat, 2)
  expect_equal(result$theta.hat[1], mean(data[,1]), tolerance = 1e-3)
  expect_equal(result$theta.hat[2], mean(data[,2]), tolerance = 1e-3)

  # Info should be a 2x2 matrix
  expect_equal(dim(result$info), c(2, 2))
  expect_equal(dim(result$sigma), c(2, 2))
})

test_that("mle_gradient_ascent works without loglike (no info/sigma)", {
  set.seed(789)
  data <- rnorm(50, mean = 2, sd = 1.5)

  score <- function(mu) {
    sum((data - mu) / (1.5^2))
  }

  result <- mle_gradient_ascent(
    theta0 = 0,
    score = score,
    options = list(
      loglike = NULL,  # No log-likelihood provided
      line_search = FALSE,
      eta = 0.01,
      max_iter = 1000,
      rel_tol = 1e-5
    )
  )

  expect_s3_class(result, "mle_gradient_ascent")
  expect_true(result$converged)

  # Without loglike, info and sigma should be NULL
  expect_true(is.null(result$info))
  expect_true(is.null(result$sigma))

  # But score should still be evaluated
  expect_true(!is.null(result$score))
})

test_that("mle_gradient_ascent finds Poisson lambda MLE", {
  # Poisson MLE: lambda = mean of data
  set.seed(321)
  true_lambda <- 4.5
  n <- 200
  data <- rpois(n, lambda = true_lambda)

  loglike <- function(lambda) {
    if (lambda <= 0) return(-Inf)
    sum(dpois(data, lambda = lambda, log = TRUE))
  }

  score <- function(lambda) {
    if (lambda <= 0) return(0)
    sum(data / lambda - 1)
  }

  sup <- function(lambda) lambda > 0

  result <- mle_gradient_ascent(
    theta0 = 1,
    score = score,
    options = list(
      loglike = loglike,
      sup = sup,
      line_search = TRUE,
      max_iter = 100,
      rel_tol = 1e-5
    )
  )

  expect_true(result$converged)
  expect_true(result$theta.hat > 0)
  expect_equal(result$theta.hat, mean(data), tolerance = 1e-3)
})

test_that("mle_gradient_ascent handles constrained optimization", {
  # Maximize log-likelihood but constrain parameter to [1, 5]
  set.seed(111)
  data <- rnorm(50, mean = 7, sd = 1)  # True mean outside constraint

  loglike <- function(mu) {
    sum(dnorm(data, mean = mu, sd = 1, log = TRUE))
  }

  score <- function(mu) {
    sum(data - mu)
  }

  # Constraint: mu in [1, 5]
  sup <- function(mu) mu >= 1 && mu <= 5
  proj <- function(mu) pmax(1, pmin(5, mu))

  result <- mle_gradient_ascent(
    theta0 = 3,
    score = score,
    options = list(
      loglike = loglike,
      sup = sup,
      proj = proj,
      line_search = TRUE,
      max_iter = 100,
      rel_tol = 1e-5
    )
  )

  # Should converge to boundary at 5
  expect_true(result$theta.hat >= 1)
  expect_true(result$theta.hat <= 5)
  expect_equal(result$theta.hat, 5, tolerance = 1e-2)
})

test_that("mle_gradient_ascent iteration count is reasonable", {
  set.seed(222)
  data <- rnorm(30, mean = 1, sd = 1)

  loglike <- function(mu) sum(dnorm(data, mean = mu, sd = 1, log = TRUE))
  score <- function(mu) sum(data - mu)

  result <- mle_gradient_ascent(
    theta0 = 0,
    score = score,
    options = list(
      loglike = loglike,
      line_search = TRUE,
      max_iter = 1000,
      rel_tol = 1e-6
    )
  )

  expect_true(result$converged)
  # Should converge in reasonable number of iterations
  expect_true(result$iter < 100)
  expect_true(result$iter > 0)
})

test_that("mle_gradient_ascent evaluates score at solution", {
  set.seed(333)
  data <- rnorm(40, mean = 2.5, sd = 1)

  loglike <- function(mu) sum(dnorm(data, mean = mu, sd = 1, log = TRUE))
  score <- function(mu) sum(data - mu)

  result <- mle_gradient_ascent(
    theta0 = 0,
    score = score,
    options = list(
      loglike = loglike,
      line_search = TRUE,
      max_iter = 100,
      rel_tol = 1e-5
    )
  )

  expect_true(!is.null(result$score))

  # At MLE, score should be close to 0
  expect_equal(abs(result$score), 0, tolerance = 1e-2)
})
