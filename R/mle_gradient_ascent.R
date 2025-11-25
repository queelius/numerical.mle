#' Maximum likelihood estimation via gradient ascent
#'
#' Performs gradient ascent optimization to find the MLE. This method
#' uses the score function (gradient of log-likelihood) to iteratively
#' improve parameter estimates.
#'
#' @param loglike Log-likelihood function taking theta as input
#' @param score Score function (gradient of log-likelihood) taking theta as input
#' @param theta0 Initial parameter guess (numeric vector)
#' @param config Configuration object (mle_config_gradient or mle_config_linesearch).
#'   Use mle_config_linesearch() for adaptive step sizes (recommended),
#'   or mle_config_gradient() for fixed step size.
#' @param constraint Optional domain constraints (mle_constraint object)
#' @return mle_numerical object with class mle_gradient_ascent containing:
#'   \item{theta.hat}{MLE estimate}
#'   \item{loglike}{Log-likelihood at MLE}
#'   \item{score}{Score vector at MLE (should be near zero)}
#'   \item{info}{Fisher information matrix (negative Hessian)}
#'   \item{sigma}{Covariance matrix (inverse of Fisher information)}
#'   \item{iter}{Number of iterations}
#'   \item{converged}{Convergence status}
#'   \item{config}{Configuration used}
#'   \item{path}{Optimization path (if trace=TRUE in config)}
#' @examples
#' \dontrun{
#' # Normal distribution MLE
#' data <- rnorm(100, mean = 5, sd = 2)
#'
#' loglike <- function(theta) {
#'   sum(dnorm(data, mean = theta[1], sd = theta[2], log = TRUE))
#' }
#'
#' score <- function(theta) {
#'   mu <- theta[1]
#'   sigma <- theta[2]
#'   c(
#'     sum((data - mu) / sigma^2),
#'     sum((data - mu)^2 / sigma^3 - 1/sigma)
#'   )
#' }
#'
#' # With line search (recommended)
#' result <- mle_gradient_ascent(
#'   loglike = loglike,
#'   score = score,
#'   theta0 = c(0, 1),
#'   config = mle_config_linesearch(max_step = 1.0)
#' )
#'
#' # With fixed step size
#' result <- mle_gradient_ascent(
#'   loglike = loglike,
#'   score = score,
#'   theta0 = c(0, 1),
#'   config = mle_config_gradient(eta = 0.1)
#' )
#'
#' # With constraints (positive variance only)
#' constraint <- mle_constraint(
#'   support = function(theta) theta[2] > 0,
#'   project = function(theta) c(theta[1], max(theta[2], 1e-8))
#' )
#'
#' result <- mle_gradient_ascent(
#'   loglike = loglike,
#'   score = score,
#'   theta0 = c(0, 1),
#'   config = mle_config_linesearch(),
#'   constraint = constraint
#' )
#'
#' # Check convergence
#' print(result$converged)
#' print(result$theta.hat)
#' print(result$score)  # Should be near zero
#' }
#' @importFrom MASS ginv
#' @importFrom numDeriv hessian
#' @export
mle_gradient_ascent <- function(
  loglike,
  score,
  theta0,
  config = mle_config_gradient(),
  constraint = mle_constraint()
) {
  # Validation
  stopifnot(
    is.function(loglike),
    is.function(score),
    is.numeric(theta0),
    inherits(config, "mle_config"),
    inherits(constraint, "mle_constraint")
  )

  if (!constraint$support(theta0)) {
    stop("Initial guess `theta0` not in support")
  }

  # Determine if we use line search based on config class
  use_linesearch <- inherits(config, "mle_config_linesearch")

  # Call internal solver with direction function = score
  result <- .mle_optimize_direction(
    loglike = loglike,
    direction_fn = score,
    theta0 = theta0,
    config = config,
    constraint = constraint,
    use_linesearch = use_linesearch
  )

  # Augment with method-specific information
  result$score <- score(result$theta.hat)

  if (!is.null(loglike)) {
    # Compute Fisher information as negative Hessian
    result$info <- tryCatch(
      -numDeriv::hessian(loglike, result$theta.hat),
      error = function(e) NULL
    )

    if (!is.null(result$info)) {
      # Compute covariance as inverse of Fisher information
      result$sigma <- tryCatch(
        MASS::ginv(result$info),
        error = function(e) NULL
      )
    }
  }

  class(result) <- c("mle_gradient_ascent", class(result))
  result
}
