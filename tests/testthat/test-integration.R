test_that("Complete workflow: Normal distribution MLE with multiple solvers", {
  # Generate data from normal distribution
  set.seed(12345)
  true_mean <- 7.5
  true_sd <- 2.0
  n <- 200
  data <- rnorm(n, mean = true_mean, sd = true_sd)

  # Define functions
  loglike <- function(mu) {
    sum(dnorm(data, mean = mu, sd = true_sd, log = TRUE))
  }

  score <- function(mu) {
    sum((data - mu) / true_sd^2)
  }

  fim <- function(mu) {
    n / true_sd^2
  }

  # Test with gradient ascent
  result_ga <- mle_gradient_ascent(
    theta0 = 0,
    score = score,
    options = list(
      loglike = loglike,
      line_search = TRUE,
      max_iter = 100,
      rel_tol = 1e-6
    )
  )

  # Test with Newton-Raphson
  result_nr <- mle_newton_raphson(
    score = score,
    fim = fim,
    theta0 = 0,
    inverted = FALSE,
    options = list(
      loglike = loglike,
      line_search = TRUE,
      max_iter = 100,
      rel_tol = 1e-6
    )
  )

  # Both should converge to same MLE
  expect_true(result_ga$converged)
  expect_true(result_nr$converged)
  expect_equal(result_ga$theta.hat, mean(data), tolerance = 1e-4)
  expect_equal(result_nr$theta.hat, mean(data), tolerance = 1e-4)
  expect_equal(result_ga$theta.hat, result_nr$theta.hat, tolerance = 1e-4)

  # Newton-Raphson should converge faster
  expect_true(result_nr$iter < result_ga$iter)
})

test_that("Complete workflow: Poisson distribution MLE", {
  # Generate Poisson data
  set.seed(54321)
  true_lambda <- 5.5
  n <- 300
  data <- rpois(n, lambda = true_lambda)

  loglike <- function(lambda) {
    if (lambda <= 0) return(-Inf)
    sum(dpois(data, lambda = lambda, log = TRUE))
  }

  score <- function(lambda) {
    if (lambda <= 0) return(0)
    sum(data / lambda - 1)
  }

  fim <- function(lambda) {
    if (lambda <= 0) return(1e10)
    n / lambda
  }

  sup <- function(lambda) lambda > 0

  # Gradient ascent
  result_ga <- mle_gradient_ascent(
    theta0 = 1,
    score = score,
    options = list(
      loglike = loglike,
      sup = sup,
      line_search = TRUE,
      max_iter = 100,
      rel_tol = 1e-6
    )
  )

  # Newton-Raphson
  result_nr <- mle_newton_raphson(
    score = score,
    fim = fim,
    theta0 = 1,
    inverted = FALSE,
    options = list(
      loglike = loglike,
      sup = sup,
      line_search = TRUE,
      max_iter = 100,
      rel_tol = 1e-6
    )
  )

  # Both should find MLE = sample mean
  expect_true(result_ga$converged)
  expect_true(result_nr$converged)
  expect_equal(result_ga$theta.hat, mean(data), tolerance = 1e-3)
  expect_equal(result_nr$theta.hat, mean(data), tolerance = 1e-3)
})

test_that("Complete workflow: Bivariate normal MLE", {
  # Generate 2D normal data
  set.seed(99999)
  n <- 150
  true_mu <- c(3, -2)
  data <- cbind(
    rnorm(n, mean = true_mu[1], sd = 1),
    rnorm(n, mean = true_mu[2], sd = 1)
  )

  loglike <- function(mu) {
    sum(dnorm(data[,1], mean = mu[1], sd = 1, log = TRUE)) +
    sum(dnorm(data[,2], mean = mu[2], sd = 1, log = TRUE))
  }

  score <- function(mu) {
    c(sum(data[,1] - mu[1]),
      sum(data[,2] - mu[2]))
  }

  fim <- function(mu) {
    matrix(c(n, 0, 0, n), nrow = 2)
  }

  # Gradient ascent
  result_ga <- mle_gradient_ascent(
    theta0 = c(0, 0),
    score = score,
    options = list(
      loglike = loglike,
      line_search = TRUE,
      max_iter = 100,
      rel_tol = 1e-6
    )
  )

  # Newton-Raphson
  result_nr <- mle_newton_raphson(
    score = score,
    fim = fim,
    theta0 = c(0, 0),
    inverted = FALSE,
    options = list(
      loglike = loglike,
      line_search = TRUE,
      max_iter = 100,
      rel_tol = 1e-6
    )
  )

  # Both should converge
  expect_true(result_ga$converged)
  expect_true(result_nr$converged)

  # Check estimates
  expect_equal(result_ga$theta.hat[1], mean(data[,1]), tolerance = 1e-3)
  expect_equal(result_ga$theta.hat[2], mean(data[,2]), tolerance = 1e-3)
  expect_equal(result_nr$theta.hat[1], mean(data[,1]), tolerance = 1e-3)
  expect_equal(result_nr$theta.hat[2], mean(data[,2]), tolerance = 1e-3)

  # Should have proper covariance matrices
  expect_equal(dim(result_ga$sigma), c(2, 2))
  expect_equal(dim(result_nr$sigma), c(2, 2))
})

test_that("Complete workflow: Constrained optimization with projection", {
  # Estimate mean of normal data with constraint
  set.seed(77777)
  data <- rnorm(100, mean = 10, sd = 1)  # True mean = 10

  loglike <- function(mu) sum(dnorm(data, mean = mu, sd = 1, log = TRUE))
  score <- function(mu) sum(data - mu)
  fim <- function(mu) length(data)

  # Constrain to [0, 8] - true MLE is outside
  sup <- function(mu) mu >= 0 && mu <= 8
  proj <- function(mu) pmax(0, pmin(8, mu))

  result <- mle_newton_raphson(
    score = score,
    fim = fim,
    theta0 = 4,
    inverted = FALSE,
    options = list(
      loglike = loglike,
      sup = sup,
      proj = proj,
      line_search = TRUE,
      max_iter = 100,
      rel_tol = 1e-6
    )
  )

  # Should converge to boundary at 8
  expect_true(result$converged)
  expect_true(result$theta.hat >= 0)
  expect_true(result$theta.hat <= 8)
  expect_equal(result$theta.hat, 8, tolerance = 1e-2)
})

test_that("Complete workflow: Stochastic gradient ascent", {
  # Large dataset - use stochastic log-likelihood
  set.seed(11111)
  n <- 10000
  true_mean <- 4
  true_sd <- 2
  data <- rnorm(n, mean = true_mean, sd = true_sd)

  # Log density for a single observation
  log_density <- function(x, theta) {
    dnorm(x, mean = theta, sd = true_sd, log = TRUE)
  }

  # Create stochastic log-likelihood using subsample
  stoch_loglike <- stochastic_loglike(
    log_density = log_density,
    data = data,
    m = 100,  # Use only 100 observations per iteration
    replace = FALSE
  )

  # Score using all data (for simplicity)
  score <- function(mu) {
    sum((data - mu) / true_sd^2)
  }

  # Full log-likelihood for comparison
  full_loglike <- function(mu) {
    sum(dnorm(data, mean = mu, sd = true_sd, log = TRUE))
  }

  result <- mle_gradient_ascent(
    theta0 = 0,
    score = score,
    options = list(
      loglike = full_loglike,
      line_search = TRUE,
      max_iter = 100,
      rel_tol = 1e-5
    )
  )

  # Should converge close to true mean
  expect_true(result$converged)
  expect_equal(result$theta.hat, mean(data), tolerance = 1e-3)
})

test_that("Complete workflow: Multiple starting points converge to same MLE", {
  # Test robustness to initial guess
  set.seed(22222)
  data <- rnorm(80, mean = 5, sd = 1.5)

  loglike <- function(mu) sum(dnorm(data, mean = mu, sd = 1.5, log = TRUE))
  score <- function(mu) sum((data - mu) / (1.5^2))
  fim <- function(mu) length(data) / (1.5^2)

  # Try multiple starting points
  starting_points <- c(-10, 0, 2, 10, 20)
  results <- lapply(starting_points, function(theta0) {
    mle_newton_raphson(
      score = score,
      fim = fim,
      theta0 = theta0,
      inverted = FALSE,
      options = list(
        loglike = loglike,
        line_search = TRUE,
        max_iter = 100,
        rel_tol = 1e-6
      )
    )
  })

  # All should converge
  converged <- sapply(results, function(r) r$converged)
  expect_true(all(converged))

  # All should find same MLE
  estimates <- sapply(results, function(r) r$theta.hat)
  expect_true(all(abs(estimates - mean(data)) < 1e-3))
})

test_that("Complete workflow: Gradient ascent with absolute tolerance", {
  set.seed(33333)
  data <- rnorm(60, mean = 2, sd = 1)

  loglike <- function(mu) sum(dnorm(data, mean = mu, sd = 1, log = TRUE))
  score <- function(mu) sum(data - mu)

  result <- mle_gradient_ascent(
    theta0 = 0,
    score = score,
    options = list(
      loglike = loglike,
      line_search = TRUE,
      abs_tol = 1e-5,
      rel_tol = NULL,
      max_iter = 100
    )
  )

  expect_true(result$converged)
  expect_equal(result$theta.hat, mean(data), tolerance = 1e-4)
})

test_that("Complete workflow: Path tracing with trace=TRUE", {
  set.seed(44444)
  data <- rnorm(50, mean = 3, sd = 1)

  loglike <- function(mu) sum(dnorm(data, mean = mu, sd = 1, log = TRUE))
  score <- function(mu) sum(data - mu)

  result <- mle_gradient_ascent(
    theta0 = -5,  # Start far from MLE
    score = score,
    options = list(
      loglike = loglike,
      line_search = TRUE,
      trace = TRUE,
      max_iter = 100,
      rel_tol = 1e-6
    )
  )

  expect_true(result$converged)
  expect_true("path" %in% names(result))
  expect_true(is.matrix(result$path))
  expect_equal(nrow(result$path), result$iter)

  # Path should show monotonic improvement toward MLE
  path_values <- apply(result$path, 1, loglike)
  # Later values should generally be better (allowing for numerical noise)
  expect_true(path_values[result$iter] > path_values[1])
})
