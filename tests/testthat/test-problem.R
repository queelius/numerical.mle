# Tests for mle_problem

test_that("mle_problem creates valid object", {
  problem <- mle_problem(
    loglike = function(theta) -sum(theta^2)
  )

  expect_true(is_mle_problem(problem))
  expect_true(is.function(problem$loglike))
  expect_true(is.function(problem$constraint$support))
})

test_that("mle_problem with constraint", {
  problem <- mle_problem(
    loglike = function(theta) -sum(theta^2),
    constraint = mle_constraint(
      support = function(theta) all(theta > 0),
      project = function(theta) pmax(theta, 1e-8)
    )
  )

  expect_true(problem$constraint$support(c(1, 2)))
  expect_false(problem$constraint$support(c(-1, 2)))
  expect_equal(problem$constraint$project(c(-1, 2)), c(1e-8, 2))
})

test_that("get_score returns analytic or numerical", {
  # With analytic score
  problem1 <- mle_problem(
    loglike = function(theta) -sum(theta^2),
    score = function(theta) -2 * theta
  )

  score1 <- get_score(problem1)
  expect_equal(score1(c(1, 2)), c(-2, -4))

  # Without analytic score (numerical)
  problem2 <- mle_problem(
    loglike = function(theta) -sum(theta^2)
  )

  score2 <- get_score(problem2)
  expect_equal(score2(c(1, 2)), c(-2, -4), tolerance = 1e-5)
})

test_that("get_fisher returns analytic or numerical", {
  problem <- mle_problem(
    loglike = function(theta) -0.5 * sum(theta^2)
  )

  fisher <- get_fisher(problem)
  result <- fisher(c(1, 2))

  # Hessian of -0.5*sum(theta^2) is -I, so Fisher is I
  expect_equal(diag(result), c(1, 1), tolerance = 1e-5)
})

test_that("update.mle_problem works", {
  problem1 <- mle_problem(
    loglike = function(theta) -sum(theta^2),
    n_obs = 100
  )

  problem2 <- update(problem1, n_obs = 200)

  expect_equal(problem1$n_obs, 100)
  expect_equal(problem2$n_obs, 200)
})
