## Tests for function transformers

test_that("with_subsampling creates valid transformed function", {
  loglike_base <- function(theta, data) {
    sum(dnorm(data, mean = theta, sd = 1, log = TRUE))
  }

  data <- rnorm(100)
  loglike_stoch <- with_subsampling(loglike_base, data = data, subsample_size = 10)

  expect_s3_class(loglike_stoch, "loglike_subsampled")
  expect_true(is.function(loglike_stoch))

  # Should work when called
  val <- loglike_stoch(0)
  expect_true(is.finite(val))
})

test_that("with_subsampling validates inputs", {
  loglike <- function(theta, data) sum(data)

  expect_error(with_subsampling("not a function", data = 1:10, subsample_size = 5))
  expect_error(with_subsampling(loglike, data = 1:10, subsample_size = 0))
  expect_error(with_subsampling(loglike, data = 1:10, subsample_size = 20))
})

test_that("penalty_l1 computes L1 norm correctly", {
  penalty <- penalty_l1()

  expect_equal(penalty(c(1, -2, 3)), 6)
  expect_equal(penalty(c(0, 0, 0)), 0)
  expect_equal(penalty(c(-5)), 5)
})

test_that("penalty_l1 with weights works", {
  weights <- c(1, 2, 1)
  penalty <- penalty_l1(weights)

  expect_equal(penalty(c(1, -2, 3)), 1*1 + 2*2 + 1*3)
})

test_that("penalty_l2 computes L2 norm squared correctly", {
  penalty <- penalty_l2()

  expect_equal(penalty(c(1, -2, 3)), 1 + 4 + 9)
  expect_equal(penalty(c(0, 0, 0)), 0)
  expect_equal(penalty(c(3, 4)), 9 + 16)
})

test_that("penalty_l2 with weights works", {
  weights <- c(1, 2, 1)
  penalty <- penalty_l2(weights)

  expect_equal(penalty(c(1, -2, 3)), 1^2 + (2*(-2))^2 + 3^2)
})

test_that("penalty_elastic_net combines L1 and L2", {
  # Pure L1 (alpha = 1)
  penalty_l1_pure <- penalty_elastic_net(alpha = 1)
  penalty_l1_ref <- penalty_l1()
  expect_equal(penalty_l1_pure(c(1, -2, 3)), penalty_l1_ref(c(1, -2, 3)))

  # Pure L2 (alpha = 0)
  penalty_l2_pure <- penalty_elastic_net(alpha = 0)
  penalty_l2_ref <- penalty_l2()
  expect_equal(penalty_l2_pure(c(1, -2, 3)), penalty_l2_ref(c(1, -2, 3)))
})

test_that("penalty_elastic_net validates alpha", {
  expect_error(penalty_elastic_net(alpha = -0.1))
  expect_error(penalty_elastic_net(alpha = 1.5))
})

test_that("with_penalty creates penalized log-likelihood", {
  loglike_base <- function(theta) -sum(theta^2)
  lambda <- 0.1

  loglike_pen <- with_penalty(loglike_base, penalty_l2(), lambda = lambda)

  expect_s3_class(loglike_pen, "loglike_penalized")
  expect_true(is.function(loglike_pen))

  # Value should be loglike - lambda * penalty
  theta <- c(1, 2, 3)
  expected <- loglike_base(theta) - lambda * penalty_l2()(theta)
  expect_equal(loglike_pen(theta), expected)
})

test_that("with_penalty validates inputs", {
  loglike <- function(theta) -sum(theta^2)

  expect_error(with_penalty("not a function", penalty_l2()))
  expect_error(with_penalty(loglike, "not a function"))
  expect_error(with_penalty(loglike, penalty_l2(), lambda = -1))
})
