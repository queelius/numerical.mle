test_that("mle_numerical constructor creates valid object", {
  # Create a simple mle_numerical object
  theta_hat <- c(1, 2)
  loglike <- function(theta) sum(dnorm(c(1,2,3), theta[1], theta[2], log=TRUE))
  score <- function(theta) c(0, 0)
  info <- matrix(c(1, 0, 0, 1), nrow=2)
  sigma <- solve(info)

  result <- mle_numerical(
    theta.hat = theta_hat,
    loglike = loglike,
    score = score,
    info = info,
    sigma = sigma,
    iter = 10L,
    converged = TRUE
  )

  expect_s3_class(result, "mle_numerical")
  expect_equal(result$theta.hat, theta_hat)
  expect_equal(result$iter, 10L)
  expect_true(result$converged)
})

test_that("mle_numerical constructor validates inputs", {
  theta_hat <- c(1, 2)
  loglike <- function(theta) 0

  # iter must be numeric
  expect_error(
    mle_numerical(theta_hat, loglike, NULL, NULL, NULL, "not_numeric", TRUE),
    "is.numeric\\(iter\\) is not TRUE"
  )

  # converged must be logical
  expect_error(
    mle_numerical(theta_hat, loglike, NULL, NULL, NULL, 10L, "not_logical"),
    "is.logical\\(converged\\) is not TRUE"
  )
})

test_that("is_mle_numerical correctly identifies mle_numerical objects", {
  theta_hat <- c(1, 2)
  loglike <- function(theta) 0

  result <- mle_numerical(theta_hat, loglike, NULL, NULL, NULL, 10L, TRUE)

  expect_true(is_mle_numerical(result))
  expect_false(is_mle_numerical(list()))
  expect_false(is_mle_numerical(NULL))
})

test_that("is_converged returns convergence status", {
  theta_hat <- c(1, 2)
  loglike <- function(theta) 0

  result_converged <- mle_numerical(theta_hat, loglike, NULL, NULL, NULL, 10L, TRUE)
  result_not_converged <- mle_numerical(theta_hat, loglike, NULL, NULL, NULL, 100L, FALSE)

  expect_true(is_converged(result_converged))
  expect_false(is_converged(result_not_converged))
})

test_that("num_iterations returns iteration count", {
  theta_hat <- c(1, 2)
  loglike <- function(theta) 0

  result <- mle_numerical(theta_hat, loglike, NULL, NULL, NULL, 42L, TRUE)

  expect_equal(num_iterations(result), 42L)
})

test_that("stochastic_loglike creates subsampling function", {
  # Simple log density for normal distribution
  log_density <- function(x, theta) {
    dnorm(x, mean = theta[1], sd = theta[2], log = TRUE)
  }

  # Create sample data
  set.seed(123)
  data <- rnorm(100, mean = 5, sd = 2)

  # Create stochastic log-likelihood with sampling without replacement
  stoch_loglike <- stochastic_loglike(log_density, data, m = 10, replace = FALSE)

  expect_type(stoch_loglike, "closure")

  # Should be able to evaluate it
  theta <- c(5, 2)
  result <- stoch_loglike(theta)
  expect_type(result, "double")
  expect_length(result, 1)
})

test_that("stochastic_loglike validates inputs", {
  log_density <- function(x, theta) dnorm(x, theta[1], theta[2], log = TRUE)
  data <- rnorm(100)

  # m must be <= length(data)
  expect_error(
    stochastic_loglike(log_density, data, m = 200, replace = FALSE),
    "m <= length\\(data\\) is not TRUE"
  )

  # log_density must be a function
  expect_error(
    stochastic_loglike("not_a_function", data, m = 10, replace = FALSE),
    "is.function\\(log_density\\) is not TRUE"
  )
})

test_that("stochastic_loglike works with replace=TRUE", {
  log_density <- function(x, theta) {
    dnorm(x, mean = theta[1], sd = theta[2], log = TRUE)
  }

  set.seed(456)
  data <- rnorm(50, mean = 3, sd = 1.5)

  # Can sample more than data size with replacement
  stoch_loglike <- stochastic_loglike(log_density, data, m = 50, replace = TRUE)

  theta <- c(3, 1.5)
  result <- stoch_loglike(theta)
  expect_type(result, "double")
  expect_length(result, 1)
  expect_false(is.nan(result))
})
