#' MLE via random restarts
#'
#' Attempts to find a global maximum by running a local solver from multiple
#' random starting points. This helps escape local maxima and find better
#' solutions when the likelihood surface is multimodal.
#'
#' @param loglike Log-likelihood function
#' @param solver Solver function (e.g., mle_gradient_ascent, mle_newton_raphson).
#'   Must accept loglike, theta0, and additional arguments.
#' @param theta0_sampler Function generating random initial parameters.
#'   Called without arguments, must return a valid theta0 vector.
#' @param n_trials Number of random trials to perform (integer, default: 100)
#' @param ... Additional arguments passed to solver
#' @return Best mle object found across all trials, with additional attribute
#'   n_trials indicating the number of trials performed.
#' @examples
#' \dontrun{
#' # Multimodal likelihood with multiple local maxima
#' loglike <- function(theta) {
#'   # Mixture of two peaks
#'   -((theta[1]-5)^2 + (theta[2]-5)^2) / 10 -
#'    ((theta[1]+3)^2 + (theta[2]+3)^2) / 10
#' }
#'
#' score <- function(theta) {
#'   c(
#'     -(theta[1]-5) / 5 - (theta[1]+3) / 5,
#'     -(theta[2]-5) / 5 - (theta[2]+3) / 5
#'   )
#' }
#'
#' # Random sampler for initial points
#' sampler <- function() {
#'   runif(2, min = -10, max = 10)
#' }
#'
#' # Try 50 random starting points
#' result <- mle_random_restart(
#'   loglike = loglike,
#'   solver = mle_gradient_ascent,
#'   theta0_sampler = sampler,
#'   n_trials = 50,
#'   score = score,
#'   config = mle_config_linesearch(max_iter = 50)
#' )
#'
#' print(result$theta.hat)
#' print(result$loglike)
#' print(result$n_trials)
#' }
#' @export
mle_random_restart <- function(
  loglike,
  solver,
  theta0_sampler,
  n_trials = 100L,
  ...
) {
  # Coerce to integer if numeric
  if (is.numeric(n_trials) && !is.integer(n_trials)) {
    n_trials <- as.integer(n_trials)
  }

  stopifnot(
    is.function(loglike),
    is.function(solver),
    is.function(theta0_sampler),
    is.integer(n_trials), n_trials > 0
  )

  best_result <- NULL
  best_loglike <- -Inf
  successful_trials <- 0

  for (trial in seq_len(n_trials)) {
    theta0 <- theta0_sampler()

    result <- tryCatch(
      solver(loglike = loglike, theta0 = theta0, ...),
      error = function(e) {
        warning(sprintf("Trial %d failed: %s", trial, e$message))
        NULL
      }
    )

    if (!is.null(result)) {
      successful_trials <- successful_trials + 1

      # Get log-likelihood from result
      current_loglike <- if (!is.null(result$loglike)) {
        result$loglike
      } else {
        loglike(result$theta.hat)
      }

      if (!is.null(current_loglike) && !is.na(current_loglike) &&
          !is.nan(current_loglike) && current_loglike > best_loglike) {
        best_result <- result
        best_loglike <- current_loglike
      }
    }
  }

  if (is.null(best_result)) {
    stop(sprintf("All %d random restart trials failed", n_trials))
  }

  # Add meta information
  class(best_result) <- c("mle_random_restart", class(best_result))
  best_result$n_trials <- n_trials
  best_result$successful_trials <- successful_trials
  best_result
}
