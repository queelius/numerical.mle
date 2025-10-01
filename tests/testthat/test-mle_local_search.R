test_that("mle_local_search finds maximum with gradient direction", {
  # Simple 1D normal distribution MLE
  # True MLE for mean is sample mean
  set.seed(123)
  data <- rnorm(100, mean = 5, sd = 2)

  # Log-likelihood for normal mean (sd known)
  loglike <- function(mu) {
    sum(dnorm(data, mean = mu, sd = 2, log = TRUE))
  }

  # Score function (gradient)
  score <- function(mu) {
    sum((data - mu) / 4)  # derivative of log-likelihood
  }

  # Direction function is just the score
  dir <- score

  result <- mle_local_search(
    dir = dir,
    theta0 = 0,  # Starting point
    loglike = loglike,
    options = list(
      loglike = loglike,
      line_search = TRUE,
      max_iter = 100,
      rel_tol = 1e-5
    )
  )

  expect_s3_class(result, "mle_local_search")
  expect_s3_class(result, "mle_numerical")
  expect_true(result$converged)
  expect_equal(result$theta.hat, mean(data), tolerance = 1e-3)
})

test_that("mle_local_search works without line search", {
  set.seed(456)
  data <- rnorm(50, mean = 3, sd = 1.5)

  loglike <- function(mu) {
    sum(dnorm(data, mean = mu, sd = 1.5, log = TRUE))
  }

  score <- function(mu) {
    sum((data - mu) / (1.5^2))
  }

  result <- mle_local_search(
    dir = score,
    theta0 = 0,
    loglike = loglike,
    options = list(
      line_search = FALSE,
      eta = 0.01,  # Small fixed step size
      max_iter = 1000,
      rel_tol = 1e-5
    )
  )

  expect_s3_class(result, "mle_local_search")
  expect_true(result$converged)
  expect_equal(result$theta.hat, mean(data), tolerance = 1e-2)
})

test_that("mle_local_search validates initial guess in support", {
  loglike <- function(x) -x^2
  dir <- function(x) -2*x

  # Support constraint: x must be positive
  sup <- function(x) x > 0

  # Initial guess outside support should error
  expect_error(
    mle_local_search(
      dir = dir,
      theta0 = -5,
      options = list(sup = sup)
    ),
    "Initial guess.*not in support"
  )
})

test_that("mle_local_search requires loglike for line search", {
  dir <- function(x) 1

  # line_search = TRUE but no loglike should error
  expect_error(
    mle_local_search(
      dir = dir,
      theta0 = 0,
      options = list(
        line_search = TRUE,
        loglike = NULL
      )
    ),
    "Line search requires log-likelihood"
  )
})

test_that("mle_local_search respects max_iter", {
  # Direction that always points away from convergence
  dir <- function(x) 1

  result <- mle_local_search(
    dir = dir,
    theta0 = 0,
    options = list(
      line_search = FALSE,
      eta = 0.1,
      max_iter = 10,
      rel_tol = 1e-10  # Very tight tolerance, won't converge
    )
  )

  expect_false(result$converged)
  expect_true(result$iter <= 10)
})

test_that("mle_local_search with projection function", {
  # Maximize -x^2, but constrain x to [-1, 1]
  loglike <- function(x) -x^2
  score <- function(x) -2*x

  # Projection onto [-1, 1]
  proj <- function(x) pmax(-1, pmin(1, x))
  sup <- function(x) x >= -1 && x <= 1

  result <- mle_local_search(
    dir = score,
    theta0 = 0.5,
    loglike = loglike,
    options = list(
      loglike = loglike,
      proj = proj,
      sup = sup,
      line_search = FALSE,
      eta = 0.5,
      max_iter = 100,
      rel_tol = 1e-5
    )
  )

  expect_true(result$theta.hat >= -1)
  expect_true(result$theta.hat <= 1)
  expect_equal(result$theta.hat, 0, tolerance = 1e-3)
})

test_that("mle_local_search stores path when trace=TRUE", {
  set.seed(789)
  data <- rnorm(30, mean = 2, sd = 1)

  loglike <- function(mu) sum(dnorm(data, mean = mu, sd = 1, log = TRUE))
  score <- function(mu) sum(data - mu)

  result <- mle_local_search(
    dir = score,
    theta0 = 0,
    loglike = loglike,
    options = list(
      loglike = loglike,
      line_search = TRUE,
      trace = TRUE,
      max_iter = 50,
      rel_tol = 1e-5
    )
  )

  expect_true("path" %in% names(result))
  expect_true(is.matrix(result$path))
  expect_equal(nrow(result$path), result$iter)
})

test_that("mle_local_search works with multidimensional parameters", {
  # 2D normal distribution MLE for mean vector
  set.seed(321)
  n <- 100
  data <- cbind(rnorm(n, 2, 1), rnorm(n, 3, 1))

  loglike <- function(mu) {
    sum(dnorm(data[,1], mean = mu[1], sd = 1, log = TRUE)) +
    sum(dnorm(data[,2], mean = mu[2], sd = 1, log = TRUE))
  }

  score <- function(mu) {
    c(sum(data[,1] - mu[1]), sum(data[,2] - mu[2]))
  }

  result <- mle_local_search(
    dir = score,
    theta0 = c(0, 0),
    loglike = loglike,
    options = list(
      loglike = loglike,
      line_search = TRUE,
      max_iter = 100,
      rel_tol = 1e-5
    )
  )

  expect_s3_class(result, "mle_local_search")
  expect_length(result$theta.hat, 2)
  expect_equal(result$theta.hat[1], mean(data[,1]), tolerance = 1e-3)
  expect_equal(result$theta.hat[2], mean(data[,2]), tolerance = 1e-3)
})

test_that("mle_local_search handles absolute tolerance", {
  set.seed(111)
  data <- rnorm(50, mean = 4, sd = 1)

  loglike <- function(mu) sum(dnorm(data, mean = mu, sd = 1, log = TRUE))
  score <- function(mu) sum(data - mu)

  result <- mle_local_search(
    dir = score,
    theta0 = 0,
    loglike = loglike,
    options = list(
      loglike = loglike,
      line_search = TRUE,
      abs_tol = 1e-4,
      rel_tol = NULL,  # Use absolute tolerance instead
      max_iter = 100
    )
  )

  expect_true(result$converged)
  expect_equal(result$theta.hat, mean(data), tolerance = 1e-3)
})

test_that("mle_local_search uses custom norm function", {
  set.seed(222)
  data <- rnorm(30, mean = 1, sd = 1)

  loglike <- function(mu) sum(dnorm(data, mean = mu, sd = 1, log = TRUE))
  score <- function(mu) sum(data - mu)

  # L2 norm instead of default infinity norm
  l2_norm <- function(x) sqrt(sum(x^2))

  result <- mle_local_search(
    dir = score,
    theta0 = 5,
    loglike = loglike,
    options = list(
      loglike = loglike,
      line_search = TRUE,
      norm = l2_norm,
      rel_tol = 1e-5,
      max_iter = 100
    )
  )

  expect_true(result$converged)
  expect_equal(result$theta.hat, mean(data), tolerance = 1e-3)
})
