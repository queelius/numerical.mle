test_that("clip_step limits step size correctly", {
  # Step within max_norm - should not be clipped
  step <- c(0.5, 0.3)
  result <- clip_step(step, max_norm = 1)
  expect_equal(result, step)

  # Step exceeding max_norm - should be clipped
  step <- c(3, 4)  # norm = 5
  result <- clip_step(step, max_norm = 1)
  expect_equal(sqrt(sum(result^2)), 1, tolerance = 1e-10)
  # Direction should be preserved
  expect_equal(result / sqrt(sum(result^2)), step / sqrt(sum(step^2)), tolerance = 1e-10)
})

test_that("clip_step handles edge cases", {
  # Zero step
  step <- c(0, 0)
  result <- clip_step(step, max_norm = 1)
  expect_equal(result, step)

  # Single element
  step <- 2
  result <- clip_step(step, max_norm = 1)
  expect_equal(abs(result), 1)

  # Very large step
  step <- c(1000, 1000)
  result <- clip_step(step, max_norm = 0.1)
  expect_equal(sqrt(sum(result^2)), 0.1, tolerance = 1e-10)
})

test_that("backtracking_line_search finds better point", {
  # Simple quadratic function with maximum at x = 2
  f <- function(x) -(x - 2)^2 + 10

  # Starting from x = 0, direction towards maximum
  x0 <- 0
  dir <- 1  # positive direction

  result <- backtracking_line_search(
    f = f,
    dir = dir,
    x0 = x0,
    max_step = 1,
    fix = NULL,
    sup = function(x) TRUE,
    debug = FALSE,
    max_iter = 100,
    min_eta = 1e-8,
    r = 0.5
  )

  expect_true(result$found_better)
  expect_true(result$max > f(x0))
  expect_true(result$argmax > x0)
})

test_that("backtracking_line_search respects support constraint", {
  # Function defined only for positive values
  f <- function(x) {
    if (x <= 0) return(-Inf)
    log(x)
  }

  # Starting from positive position, moving in negative direction
  x0 <- 5
  dir <- -1

  sup <- function(x) x > 0

  result <- backtracking_line_search(
    f = f,
    dir = dir,
    x0 = x0,
    max_step = 10,  # Would go negative without support constraint
    fix = NULL,
    sup = sup,
    debug = FALSE,
    max_iter = 100,
    min_eta = 1e-8,
    r = 0.5
  )

  # Should find a point, but it should remain positive
  if (result$found_better) {
    expect_true(result$argmax > 0)
  }
})

test_that("backtracking_line_search uses projection function", {
  # Function with constrained domain
  f <- function(x) -(x - 3)^2 + 5

  x0 <- 0
  dir <- 1

  # Projection to keep x in [0, 2]
  proj <- function(x) {
    pmax(0, pmin(2, x))
  }

  sup <- function(x) x >= 0 && x <= 2

  result <- backtracking_line_search(
    f = f,
    dir = dir,
    x0 = x0,
    max_step = 5,
    fix = proj,
    sup = sup,
    debug = FALSE,
    max_iter = 100,
    min_eta = 1e-8,
    r = 0.5
  )

  expect_true(result$argmax >= 0)
  expect_true(result$argmax <= 2)
})

test_that("backtracking_line_search handles no improvement case", {
  # Constant function - no improvement possible
  f <- function(x) 0

  x0 <- 0
  dir <- 1

  result <- backtracking_line_search(
    f = f,
    dir = dir,
    x0 = x0,
    max_step = 1,
    fix = NULL,
    sup = function(x) TRUE,
    debug = FALSE,
    max_iter = 10,
    min_eta = 1e-8,
    r = 0.5
  )

  # May or may not find "better" (equal), but should return valid result
  expect_type(result, "list")
  expect_true("found_better" %in% names(result))
  expect_true("argmax" %in% names(result))
  expect_true("max" %in% names(result))
})

test_that("grad_descent minimizes simple function", {
  # Simple quadratic to minimize: (x-3)^2
  f <- function(x) (x - 3)^2
  df <- function(x) 2 * (x - 3)

  x0 <- 0
  result <- grad_descent(
    f = f,
    x0 = x0,
    df = df,
    sup = function(x) TRUE,
    eps = 1e-6,
    lr = 0.1,
    debug = FALSE,
    max_iter = 1000
  )

  expect_true(result$converged)
  expect_equal(result$param, 3, tolerance = 1e-4)
})

test_that("grad_descent handles multidimensional optimization", {
  # Minimize sum of squares: (x-2)^2 + (y-3)^2
  f <- function(x) (x[1] - 2)^2 + (x[2] - 3)^2
  df <- function(x) c(2*(x[1] - 2), 2*(x[2] - 3))

  x0 <- c(0, 0)
  result <- grad_descent(
    f = f,
    x0 = x0,
    df = df,
    sup = function(x) TRUE,
    eps = 1e-6,
    lr = 0.1,
    debug = FALSE,
    max_iter = 1000
  )

  expect_true(result$converged)
  expect_equal(result$param[1], 2, tolerance = 1e-4)
  expect_equal(result$param[2], 3, tolerance = 1e-4)
})

test_that("grad_descent respects support constraints", {
  # Minimize (x-10)^2, but constrain x to [0, 5]
  f <- function(x) (x - 10)^2
  df <- function(x) 2 * (x - 10)
  sup <- function(x) x >= 0 && x <= 5

  x0 <- 2
  result <- grad_descent(
    f = f,
    x0 = x0,
    df = df,
    sup = sup,
    eps = 1e-6,
    lr = 0.1,
    debug = FALSE,
    max_iter = 1000
  )

  # Should converge to boundary at x=5
  expect_true(result$param >= 0)
  expect_true(result$param <= 5)
})

test_that("grad_descent handles max_iter limit", {
  # Function that would take many iterations
  f <- function(x) (x - 1000)^2
  df <- function(x) 2 * (x - 1000)

  x0 <- 0
  result <- grad_descent(
    f = f,
    x0 = x0,
    df = df,
    sup = function(x) TRUE,
    eps = 1e-10,
    lr = 0.01,  # Small learning rate
    debug = FALSE,
    max_iter = 10  # Very few iterations
  )

  expect_false(result$converged)
  expect_true(result$iter <= 10)
})
