#' Internal direction-based optimizer
#'
#' Core optimization algorithm used by gradient-based solvers.
#' This is an internal function not meant to be called directly by users.
#'
#' @param loglike Log-likelihood function
#' @param direction_fn Function computing search direction
#' @param theta0 Initial parameters
#' @param config Configuration object (mle_config or subclass)
#' @param constraint Domain constraints (mle_constraint object)
#' @param use_linesearch Whether to use backtracking line search
#' @return mle_numerical object with optimization results
#' @keywords internal
.mle_optimize_direction <- function(
  loglike,
  direction_fn,
  theta0,
  config,
  constraint,
  use_linesearch
) {
  # Build convergence checker
  check_converged <- .make_convergence_checker(config)

  converged <- FALSE
  path <- if (config$trace) {
    matrix(nrow = config$max_iter, ncol = length(theta0))
  } else {
    NULL
  }

  theta_current <- theta0
  iter <- 1L

  repeat {
    # Compute search direction
    direction <- direction_fn(theta_current)

    # Debug output
    if (config$debug && iter %% config$debug_freq == 0) {
      .print_iteration(iter, theta_current, direction, loglike)
    }

    # Take step
    if (use_linesearch) {
      step_result <- .backtracking_step(
        loglike = loglike,
        direction = direction,
        theta_current = theta_current,
        max_step = config$eta,
        constraint = constraint,
        backtrack_ratio = config$backtrack_ratio,
        max_iter = config$max_iter_ls,
        min_step = config$min_step,
        debug = config$debug
      )

      if (!step_result$success) {
        if (config$debug) {
          cat("Line search failed to find improvement\n")
        }
        # If line search fails, check if we're at an optimum (small gradient)
        dir_norm <- sqrt(sum(direction^2))
        grad_tol <- if (!is.null(config$abs_tol)) config$abs_tol else sqrt(config$rel_tol)
        if (dir_norm < grad_tol) {
          converged <- TRUE
          if (config$debug) {
            cat(sprintf("Converged: gradient norm %.6f < tol %.6f\n", dir_norm, grad_tol))
          }
        }
        break
      }

      theta_next <- step_result$theta
    } else {
      # Simple gradient step
      theta_next <- theta_current + config$eta * direction

      # Project if necessary
      if (!constraint$support(theta_next)) {
        if (config$debug) {
          cat("Projecting theta onto support\n")
        }
        theta_next <- constraint$project(theta_next)
      }
    }

    # Store path
    if (config$trace) {
      path[iter, ] <- as.vector(theta_next)
    }

    # Check convergence
    if (check_converged(theta_next, theta_current, config)) {
      converged <- TRUE
      theta_current <- theta_next
      break
    }

    # Check iteration limit
    if (iter >= config$max_iter) {
      theta_current <- theta_next
      break
    }

    theta_current <- theta_next
    iter <- iter + 1L
  }

  # Build result
  loglike_value <- if (!is.null(loglike)) loglike(theta_current) else NULL

  result <- algebraic.mle::mle(
    theta.hat = theta_current,
    loglike = loglike_value,
    superclasses = c("mle_local_search", "mle_numerical")
  )

  result$iter <- iter
  result$converged <- converged
  result$config <- config

  if (config$trace) {
    result$path <- path[1:iter, , drop = FALSE]
  }

  result
}

#' Create convergence checker function
#'
#' @param config Configuration object
#' @return Function that checks convergence
#' @keywords internal
.make_convergence_checker <- function(config) {
  if (!is.null(config$abs_tol)) {
    function(theta_new, theta_old, config) {
      config$norm(theta_new - theta_old) < config$abs_tol
    }
  } else {
    function(theta_new, theta_old, config) {
      norm_diff <- config$norm(theta_new - theta_old)
      norm_new <- config$norm(theta_new)
      if (norm_new == 0) {
        return(norm_diff < config$rel_tol)
      }
      norm_diff < config$rel_tol * norm_new
    }
  }
}

#' Print iteration information
#'
#' @param iter Iteration number
#' @param theta Current parameters
#' @param direction Search direction
#' @param loglike Log-likelihood function (can be NULL)
#' @keywords internal
.print_iteration <- function(iter, theta, direction, loglike) {
  if (is.null(loglike)) {
    cat(sprintf("| %-4d | Direction: (%s) | Theta: %s |\n",
                iter, toString(round(direction, 3)), toString(round(theta, 3))))
  } else {
    cat(sprintf("| %-8d | Direction: (%s) | Theta: %s | LogLik: %.5f |\n",
                iter, toString(round(direction, 3)), toString(round(theta, 3)),
                loglike(theta)))
  }
}

#' Backtracking line search step
#'
#' Performs a backtracking line search to find a step size that improves
#' the objective function.
#'
#' @param loglike Log-likelihood function
#' @param direction Search direction vector
#' @param theta_current Current parameter values
#' @param max_step Maximum step size
#' @param constraint Domain constraints
#' @param backtrack_ratio Backtracking multiplier (0 < r < 1)
#' @param max_iter Maximum iterations for line search
#' @param min_step Minimum step size threshold
#' @param debug Print debug information
#' @return List with success (logical) and theta (new parameter values)
#' @keywords internal
.backtracking_step <- function(
  loglike,
  direction,
  theta_current,
  max_step,
  constraint,
  backtrack_ratio,
  max_iter,
  min_step,
  debug
) {
  norm_fn <- function(x) sqrt(sum(x^2))
  dir_norm <- norm_fn(direction)

  if (dir_norm == 0) {
    if (debug) {
      cat("Direction has zero norm\n")
    }
    return(list(success = FALSE, theta = theta_current))
  }

  loglike_current <- loglike(theta_current)
  step_size <- max_step / dir_norm
  iter <- 0L

  while (step_size >= min_step && iter < max_iter) {
    theta_candidate <- theta_current + step_size * direction

    if (debug && iter == 0) {
      cat(sprintf("  Line search: initial step_size = %.6f\n", step_size))
    }

    # Handle support violations
    if (!constraint$support(theta_candidate)) {
      theta_candidate <- constraint$project(theta_candidate)
      # Recalculate step size based on projection
      step_size <- norm_fn(theta_candidate - theta_current)
      if (debug) {
        cat(sprintf("  Projected to support, new step_size = %.6f\n", step_size))
      }
    }

    if (constraint$support(theta_candidate)) {
      loglike_candidate <- loglike(theta_candidate)

      if (!is.nan(loglike_candidate) && !is.na(loglike_candidate) &&
          loglike_candidate > loglike_current) {
        if (debug) {
          cat(sprintf("  Line search succeeded: step_size = %.6f, improvement = %.6f\n",
                      step_size, loglike_candidate - loglike_current))
        }
        return(list(success = TRUE, theta = theta_candidate))
      }
    }

    step_size <- backtrack_ratio * step_size
    iter <- iter + 1L
  }

  if (debug) {
    cat(sprintf("  Line search failed after %d iterations\n", iter))
  }
  list(success = FALSE, theta = theta_current)
}

#' Generate grid of parameter values
#'
#' @param lower Lower bounds
#' @param upper Upper bounds
#' @param grid_size Grid resolution per dimension
#' @return Matrix where each row is a parameter vector
#' @keywords internal
.generate_grid <- function(lower, upper, grid_size) {
  # Handle scalar grid_size
  if (length(grid_size) == 1) {
    grid_size <- rep(grid_size, length(lower))
  }

  # Create sequences for each dimension
  grids <- lapply(seq_along(lower), function(i) {
    seq(lower[i], upper[i], length.out = grid_size[i])
  })

  # Create full grid
  grid_list <- expand.grid(grids)
  as.matrix(grid_list)
}
