#' Create an MLE Problem Specification
#'
#' Encapsulates a maximum likelihood estimation problem, separating the
#' statistical specification from the optimization strategy.
#'
#' @param loglike Log-likelihood function taking parameter vector theta
#' @param score Score function (gradient of log-likelihood). If NULL,
#'   computed numerically via numDeriv::grad when needed.
#' @param fisher Fisher information matrix function. If NULL, computed
#'   numerically via numDeriv::hessian when needed.
#' @param constraint Domain constraints as mle_constraint object
#' @param theta_names Character vector of parameter names for nice output
#' @param n_obs Number of observations (for AIC/BIC computation)
#' @return An mle_problem object
#'
#' @details
#' The problem object provides lazy evaluation of derivatives. If you don't
#' provide analytic score or fisher functions, they will be computed
#' numerically when first requested and cached.
#'
#' @examples
#' # With analytic derivatives
#' problem <- mle_problem(
#'   loglike = function(theta) sum(dnorm(data, theta[1], theta[2], log = TRUE)),
#'   score = function(theta) {
#'     c(sum(data - theta[1]) / theta[2]^2,
#'       -length(data)/theta[2] + sum((data - theta[1])^2) / theta[2]^3)
#'   },
#'   constraint = mle_constraint(
#'     support = function(theta) theta[2] > 0,
#'     project = function(theta) c(theta[1], max(theta[2], 1e-8))
#'   ),
#'   theta_names = c("mu", "sigma")
#' )
#'
#' # Without analytic derivatives (computed numerically)
#' problem <- mle_problem(
#'   loglike = function(theta) sum(dnorm(data, theta[1], theta[2], log = TRUE)),
#'   constraint = mle_constraint(
#'     support = function(theta) theta[2] > 0
#'   )
#' )
#'
#' @export
mle_problem <- function(
  loglike,
  score = NULL,
  fisher = NULL,
  constraint = NULL,

  theta_names = NULL,
  n_obs = NULL
) {
  stopifnot(is.function(loglike))
  if (!is.null(score)) stopifnot(is.function(score))
  if (!is.null(fisher)) stopifnot(is.function(fisher))
  if (!is.null(constraint)) stopifnot(inherits(constraint, "mle_constraint"))
  if (!is.null(theta_names)) stopifnot(is.character(theta_names))
  if (!is.null(n_obs)) stopifnot(is.numeric(n_obs), n_obs > 0)

  # Default constraint: no constraints

if (is.null(constraint)) {
    constraint <- mle_constraint()
  }

  structure(
    list(
      loglike = loglike,
      .score = score,
      .fisher = fisher,
      constraint = constraint,
      theta_names = theta_names,
      n_obs = n_obs,
      .cache = new.env(parent = emptyenv())
    ),
    class = "mle_problem"
  )
}

#' @export
print.mle_problem <- function(x, ...) {
  cat("MLE Problem\n")
  cat("  Parameters:", if (!is.null(x$theta_names)) paste(x$theta_names, collapse = ", ") else "unnamed", "\n")
  cat("  Score:", if (!is.null(x$.score)) "analytic" else "numerical", "\n")
  cat("  Fisher:", if (!is.null(x$.fisher)) "analytic" else "numerical", "\n")
  cat("  Constraints:", if (!identical(x$constraint, mle_constraint())) "yes" else "none", "\n")
  if (!is.null(x$n_obs)) cat("  Observations:", x$n_obs, "\n")
  invisible(x)
}

#' Get score function from problem
#'
#' Returns the score (gradient) function, computing numerically if not provided.
#'
#' @param problem An mle_problem object
#' @return Score function
#' @export
get_score <- function(problem) {
  if (!is.null(problem$.score)) {
    problem$.score
  } else {
    # Return numerical score function
    function(theta) {
      numDeriv::grad(problem$loglike, theta)
    }
  }
}

#' Get Fisher information function from problem
#'
#' Returns the Fisher information matrix function, computing numerically if not provided.
#'
#' @param problem An mle_problem object
#' @return Fisher information function
#' @export
get_fisher <- function(problem) {
  if (!is.null(problem$.fisher)) {
    problem$.fisher
  } else {
    # Return numerical Fisher (negative Hessian)
    function(theta) {
      -numDeriv::hessian(problem$loglike, theta)
    }
  }
}

#' Check if object is an mle_problem
#'
#' @param x Object to test
#' @return Logical
#' @export
is_mle_problem <- function(x) {
  inherits(x, "mle_problem")
}

#' Update an mle_problem
#'
#' Create a new problem with some fields updated.
#'
#' @param object An mle_problem
#' @param ... Named arguments to update
#' @return New mle_problem
#' @export
update.mle_problem <- function(object, ...) {
  args <- list(...)
  current <- list(
    loglike = object$loglike,
    score = object$.score,
    fisher = object$.fisher,
    constraint = object$constraint,
    theta_names = object$theta_names,
    n_obs = object$n_obs
  )
  current[names(args)] <- args
  do.call(mle_problem, current)
}
