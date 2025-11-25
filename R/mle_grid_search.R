#' MLE via grid search
#'
#' Performs exhaustive grid search over a bounded parameter space.
#' Optionally refines each grid point using a local solver.
#'
#' @param loglike Log-likelihood function
#' @param lower Lower bounds for parameters (numeric vector)
#' @param upper Upper bounds for parameters (numeric vector)
#' @param grid_size Grid resolution. Either a single integer (same resolution
#'   per dimension) or a vector of integers (one per dimension).
#' @param refine_solver Optional local solver to refine each grid point.
#'   If NULL, only evaluates loglike at grid points without refinement.
#' @param ... Additional arguments passed to refine_solver
#' @return mle object with best solution found, including:
#'   \item{theta.hat}{Best parameter estimate}
#'   \item{loglike}{Log-likelihood at best point}
#'   \item{grid_size}{Grid resolution used}
#'   \item{n_evaluated}{Number of grid points evaluated}
#' @examples
#' \dontrun{
#' # Simple grid search without refinement
#' loglike <- function(theta) {
#'   -(theta[1]^2 + theta[2]^2)
#' }
#'
#' result <- mle_grid_search(
#'   loglike = loglike,
#'   lower = c(-5, -5),
#'   upper = c(5, 5),
#'   grid_size = 20
#' )
#'
#' # Grid search with local refinement
#' score <- function(theta) {
#'   -2 * theta
#' }
#'
#' result <- mle_grid_search(
#'   loglike = loglike,
#'   lower = c(-5, -5),
#'   upper = c(5, 5),
#'   grid_size = 10,
#'   refine_solver = mle_gradient_ascent,
#'   score = score,
#'   config = mle_config_gradient(eta = 0.1, max_iter = 20)
#' )
#'
#' # Different resolution per dimension
#' result <- mle_grid_search(
#'   loglike = loglike,
#'   lower = c(-5, -5),
#'   upper = c(5, 5),
#'   grid_size = c(20, 10)  # 20 points in dim 1, 10 in dim 2
#' )
#' }
#' @export
mle_grid_search <- function(
  loglike,
  lower,
  upper,
  grid_size,
  refine_solver = NULL,
  ...
) {
  stopifnot(
    is.function(loglike),
    is.numeric(lower),
    is.numeric(upper),
    length(lower) == length(upper),
    all(lower <= upper),
    is.numeric(grid_size), all(grid_size > 0)
  )

  if (!is.null(refine_solver)) {
    stopifnot(is.function(refine_solver))
  }

  # Generate grid
  grid <- .generate_grid(lower, upper, grid_size)

  best_result <- NULL
  best_loglike <- -Inf
  n_evaluated <- 0

  for (i in seq_len(nrow(grid))) {
    theta0 <- grid[i, ]

    if (!is.null(refine_solver)) {
      # Use solver to refine grid point
      result <- tryCatch(
        refine_solver(loglike = loglike, theta0 = theta0, ...),
        error = function(e) NULL
      )
    } else {
      # Just evaluate at grid point
      ll <- tryCatch(
        loglike(theta0),
        error = function(e) -Inf
      )

      result <- if (!is.na(ll) && !is.nan(ll) && is.finite(ll)) {
        algebraic.mle::mle(
          theta.hat = theta0,
          loglike = ll,
          superclasses = "mle_grid"
        )
      } else {
        NULL
      }
    }

    if (!is.null(result)) {
      n_evaluated <- n_evaluated + 1

      # Get log-likelihood from result
      current_loglike <- if (!is.null(result$loglike)) {
        result$loglike
      } else {
        tryCatch(loglike(result$theta.hat), error = function(e) -Inf)
      }

      if (is.finite(current_loglike) && current_loglike > best_loglike) {
        best_result <- result
        best_loglike <- current_loglike
      }
    }
  }

  if (is.null(best_result)) {
    stop("Grid search found no valid solution")
  }

  # Add meta information
  class(best_result) <- c("mle_grid_search", class(best_result))
  best_result$grid_size <- grid_size
  best_result$n_evaluated <- n_evaluated
  best_result$n_grid_points <- nrow(grid)
  best_result
}