test_that("mle_newton_raphson finds normal distribution mean", {
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

  # Score function (gradient)
  score <- function(mu) {
    sum((data - mu) / true_sd^2)
  }

  # Fisher Information Matrix (scalar in 1D case)
  fim <- function(mu) {
    n / true_sd^2
  }

  result <- mle_newton_raphson(
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

  expect_s3_class(result, "mle_newton_raphson")
  expect_s3_class(result, "mle_local_search")
  expect_s3_class(result, "mle_numerical")
  expect_true(result$converged)
  expect_equal(result$theta.hat, mean(data), tolerance = 1e-4)

  # Should have score, sigma, and info
  expect_true(!is.null(result$score))
  expect_true(!is.null(result$sigma))
  expect_true(!is.null(result$info))
})

test_that("mle_newton_raphson validates inputs", {
  # score must be a function
  expect_error(
    mle_newton_raphson(
      score = "not_a_function",
      fim = function(x) 1,
      theta0 = 0
    ),
    "score.*must be a function"
  )

  # fim must be a function
  expect_error(
    mle_newton_raphson(
      score = function(x) x,
      fim = "not_a_function",
      theta0 = 0
    ),
    "fim.*must be a function"
  )

  # inverted must be logical
  expect_error(
    mle_newton_raphson(
      score = function(x) x,
      fim = function(x) 1,
      theta0 = 0,
      inverted = "not_logical"
    ),
    "inverted must be a logical"
  )
})

test_that("mle_newton_raphson works with inverted=TRUE (covariance)", {
  set.seed(456)
  data <- rnorm(50, mean = 3, sd = 1.5)

  score <- function(mu) {
    sum((data - mu) / (1.5^2))
  }

  # Provide covariance matrix instead of FIM
  covar <- function(mu) {
    1.5^2 / length(data)
  }

  result <- mle_newton_raphson(
    score = score,
    fim = covar,
    theta0 = 0,
    inverted = TRUE,  # fim is actually covariance
    options = list(
      line_search = FALSE,
      eta = 1,
      max_iter = 100,
      rel_tol = 1e-6
    )
  )

  expect_s3_class(result, "mle_newton_raphson")
  expect_true(result$converged)
  expect_equal(result$theta.hat, mean(data), tolerance = 1e-3)
})

test_that("mle_newton_raphson works with multidimensional parameters", {
  # 2D normal MLE
  set.seed(789)
  n <- 100
  true_mu <- c(2, -3)
  data <- cbind(rnorm(n, true_mu[1], 1), rnorm(n, true_mu[2], 1))

  loglike <- function(mu) {
    sum(dnorm(data[,1], mean = mu[1], sd = 1, log = TRUE)) +
    sum(dnorm(data[,2], mean = mu[2], sd = 1, log = TRUE))
  }

  score <- function(mu) {
    c(sum(data[,1] - mu[1]),
      sum(data[,2] - mu[2]))
  }

  # Fisher Information Matrix (diagonal for independent components)
  fim <- function(mu) {
    matrix(c(n, 0, 0, n), nrow = 2)
  }

  result <- mle_newton_raphson(
    score = score,
    fim = fim,
    theta0 = c(0, 0),
    inverted = FALSE,
    options = list(
      loglike = loglike,
      line_search = TRUE,
      max_iter = 50,
      rel_tol = 1e-6
    )
  )

  expect_true(result$converged)
  expect_length(result$theta.hat, 2)
  expect_equal(result$theta.hat[1], mean(data[,1]), tolerance = 1e-4)
  expect_equal(result$theta.hat[2], mean(data[,2]), tolerance = 1e-4)

  # Check dimensions
  expect_equal(dim(result$sigma), c(2, 2))
  expect_equal(dim(result$info), c(2, 2))
  expect_length(result$score, 2)
})

test_that("mle_newton_raphson converges faster than gradient ascent", {
  # Newton-Raphson should converge in fewer iterations due to second-order info
  set.seed(321)
  data <- rnorm(100, mean = 4, sd = 2)

  loglike <- function(mu) sum(dnorm(data, mean = mu, sd = 2, log = TRUE))
  score <- function(mu) sum((data - mu) / 4)
  fim <- function(mu) length(data) / 4

  result <- mle_newton_raphson(
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

  expect_true(result$converged)
  # Newton-Raphson should converge very quickly for this quadratic problem
  expect_true(result$iter < 20)
})

test_that("mle_newton_raphson handles Poisson distribution", {
  # Poisson MLE with Newton-Raphson
  set.seed(111)
  true_lambda <- 3.5
  n <- 150
  data <- rpois(n, lambda = true_lambda)

  loglike <- function(lambda) {
    if (lambda <= 0) return(-Inf)
    sum(dpois(data, lambda = lambda, log = TRUE))
  }

  score <- function(lambda) {
    if (lambda <= 0) return(0)
    sum(data / lambda - 1)
  }

  # Fisher information for Poisson
  fim <- function(lambda) {
    if (lambda <= 0) return(1e10)
    n / lambda
  }

  sup <- function(lambda) lambda > 0

  result <- mle_newton_raphson(
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

  expect_true(result$converged)
  expect_true(result$theta.hat > 0)
  expect_equal(result$theta.hat, mean(data), tolerance = 1e-4)
})

test_that("mle_newton_raphson handles constrained optimization", {
  # Constrain parameter to bounded region
  set.seed(222)
  data <- rnorm(50, mean = 8, sd = 1)  # True mean outside constraint

  loglike <- function(mu) sum(dnorm(data, mean = mu, sd = 1, log = TRUE))
  score <- function(mu) sum(data - mu)
  fim <- function(mu) length(data)

  # Constraint: mu in [2, 6]
  sup <- function(mu) mu >= 2 && mu <= 6
  proj <- function(mu) pmax(2, pmin(6, mu))

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

  # Should converge to boundary at 6
  expect_true(result$theta.hat >= 2)
  expect_true(result$theta.hat <= 6)
  expect_equal(result$theta.hat, 6, tolerance = 1e-2)
})

test_that("mle_newton_raphson score is near zero at MLE", {
  set.seed(333)
  data <- rnorm(60, mean = 1.5, sd = 1)

  loglike <- function(mu) sum(dnorm(data, mean = mu, sd = 1, log = TRUE))
  score <- function(mu) sum(data - mu)
  fim <- function(mu) length(data)

  result <- mle_newton_raphson(
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

  expect_true(result$converged)

  # At MLE, score should be approximately 0
  expect_equal(abs(result$score), 0, tolerance = 1e-3)
})

test_that("mle_newton_raphson works without line search", {
  set.seed(444)
  data <- rnorm(40, mean = 2, sd = 1)

  score <- function(mu) sum(data - mu)
  fim <- function(mu) length(data)

  result <- mle_newton_raphson(
    score = score,
    fim = fim,
    theta0 = 0,
    inverted = FALSE,
    options = list(
      line_search = FALSE,
      eta = 1,  # Full Newton step
      max_iter = 100,
      rel_tol = 1e-6
    )
  )

  expect_s3_class(result, "mle_newton_raphson")
  expect_true(result$converged)
  expect_equal(result$theta.hat, mean(data), tolerance = 1e-4)
})
