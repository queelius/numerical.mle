## Grid Search Tests

test_that("mle_grid_search works without refinement", {
  loglike <- function(theta) -(theta[1]^2 + theta[2]^2)

  result <- mle_grid_search(
    loglike = loglike,
    lower = c(-2, -2),
    upper = c(2, 2),
    grid_size = 5
  )

  expect_s3_class(result, "mle_grid_search")
  expect_length(result$theta.hat, 2)
  expect_true(all(result$theta.hat >= -2 & result$theta.hat <= 2))
  expect_true(abs(result$theta.hat[1]) <= 1)  # Should be near 0
  expect_true(abs(result$theta.hat[2]) <= 1)
  expect_true(result$n_evaluated > 0)
})

test_that("mle_grid_search validates inputs", {
  loglike <- function(theta) -sum(theta^2)

  # Invalid loglike
  expect_error(
    mle_grid_search(
      loglike = "not a function",
      lower = c(-1, -1),
      upper = c(1, 1),
      grid_size = 5
    )
  )

  # Mismatched dimensions
  expect_error(
    mle_grid_search(
      loglike = loglike,
      lower = c(-1, -1),
      upper = c(1, 1, 1),
      grid_size = 5
    )
  )

  # Invalid bounds
  expect_error(
    mle_grid_search(
      loglike = loglike,
      lower = c(1, 1),
      upper = c(-1, -1),
      grid_size = 5
    )
  )
})

test_that("mle_grid_search with 1D problem", {
  loglike <- function(theta) -(theta - 3)^2

  result <- mle_grid_search(
    loglike = loglike,
    lower = 0,
    upper = 6,
    grid_size = 10
  )

  expect_length(result$theta.hat, 1)
  expect_true(abs(result$theta.hat - 3) < 1)
})

test_that("mle_grid_search records metadata", {
  loglike <- function(theta) -sum(theta^2)

  result <- mle_grid_search(
    loglike = loglike,
    lower = c(-1, -1),
    upper = c(1, 1),
    grid_size = 3
  )

  expect_equal(result$grid_size, 3)
  expect_equal(result$n_evaluated, 9)
  expect_true(is.finite(result$loglike))
})

## Random Restart Tests

test_that("mle_random_restart basic functionality", {
  loglike <- function(theta) -sum(theta^2)
  score <- function(theta) -2 * theta
  sampler <- function() runif(2, -5, 5)

  result <- mle_random_restart(
    loglike = loglike,
    solver = mle_grad,
    theta0_sampler = sampler,
    n_trials = 5,
    score = score,
    max_iter = 20
  )

  expect_s3_class(result, "mle_random_restart")
  expect_equal(result$n_trials, 5)
  expect_true(result$successful_trials <= 5)
  expect_true(result$successful_trials > 0)
})

test_that("mle_random_restart validates inputs", {
  loglike <- function(theta) -sum(theta^2)
  solver <- mle_grad
  sampler <- function() runif(2, -5, 5)

  # Invalid loglike
  expect_error(
    mle_random_restart(
      loglike = "not a function",
      solver = solver,
      theta0_sampler = sampler,
      n_trials = 5
    )
  )

  # Invalid solver
  expect_error(
    mle_random_restart(
      loglike = loglike,
      solver = "not a function",
      theta0_sampler = sampler,
      n_trials = 5
    )
  )

  # Invalid sampler
  expect_error(
    mle_random_restart(
      loglike = loglike,
      solver = solver,
      theta0_sampler = "not a function",
      n_trials = 5
    )
  )

  # Invalid n_trials
  expect_error(
    mle_random_restart(
      loglike = loglike,
      solver = solver,
      theta0_sampler = sampler,
      n_trials = 0
    )
  )
})

test_that("mle_random_restart works with Newton-Raphson", {
  loglike <- function(theta) -sum(theta^2)
  score <- function(theta) -2 * theta
  fisher <- function(theta) diag(2, length(theta))
  sampler <- function() runif(2, -5, 5)

  result <- mle_random_restart(
    loglike = loglike,
    solver = mle_nr,
    theta0_sampler = sampler,
    n_trials = 3,
    score = score,
    fisher = fisher,
    max_iter = 10
  )

  expect_s3_class(result, "mle_random_restart")
  expect_true(result$successful_trials > 0)
})

test_that("mle_random_restart stores best solution", {
  loglike <- function(theta) -(theta[1] - 3)^2 - (theta[2] + 2)^2
  score <- function(theta) c(-2 * (theta[1] - 3), -2 * (theta[2] + 2))
  sampler <- function() runif(2, -10, 10)

  result <- mle_random_restart(
    loglike = loglike,
    solver = mle_grad,
    theta0_sampler = sampler,
    n_trials = 10,
    score = score,
    max_iter = 30
  )

  # Should find optimum near (3, -2)
  expect_true(abs(result$theta.hat[1] - 3) < 1)
  expect_true(abs(result$theta.hat[2] - (-2)) < 1)
})
