#' BFGS Solver
#'
#' Creates a solver using the BFGS quasi-Newton method via \code{optim()}.
#' BFGS approximates the Hessian from gradient information, providing
#' second-order-like convergence without computing the Hessian directly.
#'
#' @param max_iter Maximum number of iterations
#' @param tol Convergence tolerance (passed to optim's reltol)
#' @param report Reporting frequency (0 = no reporting)
#' @return A solver function with signature (problem, theta0, trace) -> mle_result
#'
#' @details
#' BFGS is often a good default choice: it's more robust than Newton-Raphson
#' (no matrix inversion issues) and faster than gradient ascent (uses
#' curvature information).
#'
#' The solver automatically uses the score function from the problem if
#' available, otherwise computes gradients numerically.
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' result <- bfgs()(problem, c(0, 1))
#'
#' # Race BFGS against gradient ascent
#' strategy <- bfgs() %|% gradient_ascent()
#' }
#'
#' @export
bfgs <- function(max_iter = 100L, tol = 1e-8, report = 0L) {
  max_iter <- as.integer(max_iter)
  report <- as.integer(report)

  function(problem, theta0, trace = mle_trace()) {
    stopifnot(is_mle_problem(problem))
    stopifnot(is.numeric(theta0))

    loglike <- problem$loglike
    score_fn <- get_score(problem)
    constraint <- problem$constraint

    # Check initial point
    if (!constraint$support(theta0)) {
      theta0 <- constraint$project(theta0)
    }

    # optim minimizes, so negate for maximization
    fn <- function(theta) {
      if (!constraint$support(theta)) return(Inf)
      -loglike(theta)
    }

    gr <- function(theta) {
      if (!constraint$support(theta)) return(rep(NA, length(theta)))
      -score_fn(theta)
    }

    # Run optim
    result <- optim(
      par = theta0,
      fn = fn,
      gr = gr,
      method = "BFGS",
      control = list(
        maxit = max_iter,
        reltol = tol,
        trace = report,
        REPORT = if (report > 0) report else 1
      ),
      hessian = TRUE
    )

    # Convert to mle_numerical
    sol <- list(
      par = result$par,
      value = -result$value,  # un-negate
      convergence = result$convergence,
      hessian = -result$hessian  # un-negate
    )

    mle_result <- algebraic.mle::mle_numerical(
      sol = sol,
      superclasses = "mle_bfgs"
    )

    mle_result$iterations <- result$counts["function"]
    mle_result$solver <- "bfgs"
    mle_result$optim_result <- result

    mle_result
  }
}

#' L-BFGS-B Solver (Box Constrained)
#'
#' Creates a solver using L-BFGS-B, a limited-memory BFGS variant that
#' supports box constraints (lower and upper bounds on parameters).
#'
#' @param lower Lower bounds on parameters (can be -Inf)
#' @param upper Upper bounds on parameters (can be Inf)
#' @param max_iter Maximum number of iterations
#' @param tol Convergence tolerance
#' @return A solver function
#'
#' @details
#' Unlike the constraint system in mle_problem (which uses projection),
#' L-BFGS-B handles box constraints natively within the algorithm.
#' Use this when you have simple bound constraints.
#'
#' @examples
#' \dontrun{
#' # Positive parameters only
#' solver <- lbfgsb(lower = c(-Inf, 0), upper = c(Inf, Inf))
#' result <- solver(problem, c(0, 1))
#' }
#'
#' @export
lbfgsb <- function(lower = -Inf, upper = Inf, max_iter = 100L, tol = 1e-8) {
  max_iter <- as.integer(max_iter)

  function(problem, theta0, trace = mle_trace()) {
    stopifnot(is_mle_problem(problem))
    stopifnot(is.numeric(theta0))

    loglike <- problem$loglike
    score_fn <- get_score(problem)

    # Expand bounds if scalar
    if (length(lower) == 1) lower <- rep(lower, length(theta0))
    if (length(upper) == 1) upper <- rep(upper, length(theta0))

    # optim minimizes
    fn <- function(theta) -loglike(theta)
    gr <- function(theta) -score_fn(theta)

    result <- optim(
      par = theta0,
      fn = fn,
      gr = gr,
      method = "L-BFGS-B",
      lower = lower,
      upper = upper,
      control = list(maxit = max_iter, factr = tol / .Machine$double.eps),
      hessian = TRUE
    )

    sol <- list(
      par = result$par,
      value = -result$value,
      convergence = result$convergence,
      hessian = -result$hessian
    )

    mle_result <- algebraic.mle::mle_numerical(
      sol = sol,
      superclasses = "mle_lbfgsb"
    )

    mle_result$iterations <- result$counts["function"]
    mle_result$solver <- "lbfgsb"
    mle_result$optim_result <- result

    mle_result
  }
}

#' Nelder-Mead Solver (Derivative-Free)
#'
#' Creates a solver using the Nelder-Mead simplex method via \code{optim()}.
#' This is a derivative-free method useful when gradients are unavailable
#' or unreliable.
#'
#' @param max_iter Maximum number of iterations
#' @param tol Convergence tolerance
#' @return A solver function
#'
#' @details
#' Nelder-Mead doesn't use gradient information, making it robust but
#' potentially slower. It's useful as a fallback when gradient-based
#' methods fail, or for problems with non-smooth likelihoods.
#'
#' @examples
#' \dontrun{
#' # Use when gradients are problematic
#' result <- nelder_mead()(problem, c(0, 1))
#'
#' # Race against gradient methods
#' strategy <- gradient_ascent() %|% nelder_mead()
#' }
#'
#' @export
nelder_mead <- function(max_iter = 500L, tol = 1e-8) {
  max_iter <- as.integer(max_iter)

  function(problem, theta0, trace = mle_trace()) {
    stopifnot(is_mle_problem(problem))
    stopifnot(is.numeric(theta0))

    loglike <- problem$loglike
    constraint <- problem$constraint

    # Check initial point
    if (!constraint$support(theta0)) {
      theta0 <- constraint$project(theta0)
    }

    # optim minimizes
    fn <- function(theta) {
      if (!constraint$support(theta)) return(Inf)
      -loglike(theta)
    }

    result <- optim(
      par = theta0,
      fn = fn,
      method = "Nelder-Mead",
      control = list(maxit = max_iter, reltol = tol),
      hessian = TRUE
    )

    sol <- list(
      par = result$par,
      value = -result$value,
      convergence = result$convergence,
      hessian = -result$hessian
    )

    mle_result <- algebraic.mle::mle_numerical(
      sol = sol,
      superclasses = "mle_nelder_mead"
    )

    mle_result$iterations <- result$counts["function"]
    mle_result$solver <- "nelder_mead"
    mle_result$optim_result <- result

    mle_result
  }
}
