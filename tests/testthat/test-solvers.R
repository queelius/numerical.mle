# Tests for solver factories

# Helper to create standard normal MLE problem
make_normal_problem <- function(n = 100, true_mu = 5, true_sigma = 2) {
  set.seed(42)
  data <- rnorm(n, mean = true_mu, sd = true_sigma)

  mle_problem(
    loglike = function(theta) {
      if (theta[2] <= 0) return(-Inf)
      sum(dnorm(data, theta[1], theta[2], log = TRUE))
    },
    score = function(theta) {
      mu <- theta[1]
      sigma <- theta[2]
      n <- length(data)
      d_mu <- sum(data - mu) / sigma^2
      d_sigma <- -n / sigma + sum((data - mu)^2) / sigma^3
      c(d_mu, d_sigma)
    },
    constraint = mle_constraint(
      support = function(theta) theta[2] > 0,
      project = function(theta) c(theta[1], max(theta[2], 1e-8))
    ),
    theta_names = c("mu", "sigma"),
    n_obs = n
  )
}

test_that("gradient_ascent solves normal MLE", {
  problem <- make_normal_problem()
  solver <- gradient_ascent(max_iter = 200)

  result <- solver(problem, c(0, 1))

  expect_true(is_mle_numerical(result))
  expect_equal(result$theta.hat[1], 5, tolerance = 0.5)
  expect_equal(result$theta.hat[2], 2, tolerance = 0.5)
})

test_that("newton_raphson solves normal MLE", {
  problem <- make_normal_problem()
  solver <- newton_raphson(max_iter = 50)

  result <- solver(problem, c(4, 1.5))

  expect_true(is_mle_numerical(result))
  expect_equal(result$theta.hat[1], 5, tolerance = 0.5)
  expect_equal(result$theta.hat[2], 2, tolerance = 0.5)
})

test_that("bfgs solves normal MLE", {
  problem <- make_normal_problem()
  solver <- bfgs(max_iter = 100)

  # BFGS needs reasonable starting point to avoid overshooting
  result <- solver(problem, c(4, 1.5))

  expect_true(is_mle_numerical(result))
  expect_equal(result$theta.hat[1], 5, tolerance = 0.5)
  expect_equal(result$theta.hat[2], 2, tolerance = 0.5)
})

test_that("nelder_mead solves normal MLE", {
  problem <- make_normal_problem()
  solver <- nelder_mead(max_iter = 500)

  result <- solver(problem, c(0, 1))

  expect_true(is_mle_numerical(result))
  expect_equal(result$theta.hat[1], 5, tolerance = 0.5)
  expect_equal(result$theta.hat[2], 2, tolerance = 0.5)
})

test_that("grid_search finds reasonable starting point", {
  problem <- make_normal_problem()
  solver <- grid_search(lower = c(0, 0.5), upper = c(10, 4), n = 10)

  result <- solver(problem, c(0, 1))

  expect_true(is_mle_numerical(result))
  # Grid search won't be exact but should be in ballpark
  expect_true(abs(result$theta.hat[1] - 5) < 3)
  expect_true(abs(result$theta.hat[2] - 2) < 2)
})

test_that("random_search finds reasonable starting point", {
  problem <- make_normal_problem()
  sampler <- uniform_sampler(c(0, 0.5), c(10, 4))
  solver <- random_search(sampler, n = 100)

  result <- solver(problem, c(0, 1))

  expect_true(is_mle_numerical(result))
})
