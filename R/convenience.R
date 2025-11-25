#' Quick gradient ascent with sensible defaults
#'
#' Convenience wrapper for mle_gradient_ascent with simplified interface.
#' Automatically enables line search for better convergence.
#'
#' @param loglike Log-likelihood function
#' @param score Score function (gradient)
#' @param theta0 Initial parameters
#' @param use_linesearch Use backtracking line search (default: TRUE)
#' @param ... Additional config parameters passed to mle_config_linesearch
#'   or mle_config_gradient
#' @return mle_gradient_ascent object
#' @examples
#' \dontrun{
#' # Quick usage with defaults
#' result <- mle_grad(loglike, score, theta0 = c(0, 1))
#'
#' # Override config parameters
#' result <- mle_grad(
#'   loglike, score, theta0 = c(0, 1),
#'   max_iter = 200,
#'   rel_tol = 1e-6
#' )
#'
#' # Without line search
#' result <- mle_grad(
#'   loglike, score, theta0 = c(0, 1),
#'   use_linesearch = FALSE,
#'   eta = 0.1
#' )
#' }
#' @export
mle_grad <- function(
  loglike,
  score,
  theta0,
  use_linesearch = TRUE,
  ...
) {
  config <- if (use_linesearch) {
    mle_config_linesearch(...)
  } else {
    mle_config_gradient(...)
  }

  mle_gradient_ascent(
    loglike = loglike,
    score = score,
    theta0 = theta0,
    config = config
  )
}

#' Quick Newton-Raphson with sensible defaults
#'
#' Convenience wrapper for mle_newton_raphson with simplified interface.
#' Always uses line search for stability (recommended for Newton-Raphson).
#'
#' @param loglike Log-likelihood function
#' @param score Score function
#' @param fisher Fisher information or covariance function
#' @param theta0 Initial parameters
#' @param inverted Is fisher the covariance matrix? (default: FALSE)
#' @param ... Additional config parameters passed to mle_config_linesearch
#' @return mle_newton_raphson object
#' @examples
#' \dontrun{
#' # Quick usage with Fisher information matrix
#' result <- mle_nr(loglike, score, fisher, theta0 = c(0, 1))
#'
#' # With covariance matrix (inverted FIM)
#' result <- mle_nr(
#'   loglike, score, covariance,
#'   theta0 = c(0, 1),
#'   inverted = TRUE
#' )
#'
#' # Override config
#' result <- mle_nr(
#'   loglike, score, fisher,
#'   theta0 = c(0, 1),
#'   max_iter = 50,
#'   max_step = 0.5
#' )
#' }
#' @export
mle_nr <- function(
  loglike,
  score,
  fisher,
  theta0,
  inverted = FALSE,
  ...
) {
  mle_newton_raphson(
    loglike = loglike,
    score = score,
    fisher = fisher,
    theta0 = theta0,
    config = mle_config_linesearch(...),
    inverted = inverted
  )
}

#' Quick constrained optimization
#'
#' Convenience wrapper for constrained optimization with simplified interface.
#' Automatically creates constraint object from support and projection functions.
#'
#' @param solver Solver function (e.g., mle_grad, mle_nr)
#' @param support Support function (returns TRUE if theta is valid)
#' @param project Projection function (maps invalid theta to valid theta)
#' @param ... Arguments passed to solver
#' @return mle object from solver
#' @examples
#' \dontrun{
#' # Constrain parameters to be positive
#' result <- with_constraint(
#'   solver = mle_grad,
#'   support = function(theta) all(theta > 0),
#'   project = function(theta) pmax(theta, 1e-8),
#'   loglike = loglike,
#'   score = score,
#'   theta0 = c(1, 1)
#' )
#' }
#' @export
with_constraint <- function(
  solver,
  support,
  project,
  ...
) {
  constraint <- mle_constraint(support = support, project = project)
  solver(..., constraint = constraint)
}
