#' Maximum likelihood estimation via Newton-Raphson
#'
#' Performs Newton-Raphson optimization to find the MLE. This second-order method
#' uses both the score (gradient) and Fisher information matrix (Hessian) to
#' achieve faster convergence than gradient ascent.
#'
#' @param loglike Log-likelihood function taking theta as input
#' @param score Score function (gradient of log-likelihood) taking theta as input
#' @param fisher Fisher information matrix function. Either the FIM itself or
#'   its inverse (covariance matrix), depending on the \code{inverted} parameter.
#' @param theta0 Initial parameter guess (numeric vector)
#' @param config Configuration object (typically mle_config_linesearch).
#'   Newton-Raphson benefits from line search to ensure stability.
#' @param constraint Optional domain constraints (mle_constraint object)
#' @param inverted Logical. If TRUE, \code{fisher} is the covariance matrix
#'   (inverse of FIM). If FALSE (default), \code{fisher} is the FIM.
#' @return mle_numerical object with class mle_newton_raphson containing:
#'   \item{theta.hat}{MLE estimate}
#'   \item{loglike}{Log-likelihood at MLE}
#'   \item{score}{Score vector at MLE (should be near zero)}
#'   \item{info}{Fisher information matrix}
#'   \item{sigma}{Covariance matrix (inverse of Fisher information)}
#'   \item{iter}{Number of iterations}
#'   \item{converged}{Convergence status}
#'   \item{config}{Configuration used}
#'   \item{path}{Optimization path (if trace=TRUE in config)}
#' @examples
#' \dontrun{
#' # Normal distribution MLE with Newton-Raphson
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
#' fisher <- function(theta) {
#'   n <- length(data)
#'   sigma <- theta[2]
#'   matrix(c(
#'     n / sigma^2, 0,
#'     0, 2*n / sigma^2
#'   ), nrow = 2)
#' }
#'
#' # Standard usage with FIM
#' result <- mle_newton_raphson(
#'   loglike = loglike,
#'   score = score,
#'   fisher = fisher,
#'   theta0 = c(0, 1),
#'   config = mle_config_linesearch()
#' )
#'
#' # Using inverted FIM (covariance matrix)
#' covariance <- function(theta) {
#'   MASS::ginv(fisher(theta))
#' }
#'
#' result <- mle_newton_raphson(
#'   loglike = loglike,
#'   score = score,
#'   fisher = covariance,
#'   theta0 = c(0, 1),
#'   inverted = TRUE
#' )
#'
#' # With constraints
#' constraint <- mle_constraint(
#'   support = function(theta) theta[2] > 0,
#'   project = function(theta) c(theta[1], max(theta[2], 1e-8))
#' )
#'
#' result <- mle_newton_raphson(
#'   loglike = loglike,
#'   score = score,
#'   fisher = fisher,
#'   theta0 = c(0, 1),
#'   config = mle_config_linesearch(),
#'   constraint = constraint
#' )
#'
#' # Faster convergence than gradient ascent
#' print(result$iter)  # Typically fewer iterations
#' print(result$score)  # Should be very close to zero
#' }
#' @importFrom MASS ginv
#' @export
mle_newton_raphson <- function(
  loglike,
  score,
  fisher,
  theta0,
  config = mle_config_linesearch(),
  constraint = mle_constraint(),
  inverted = FALSE
) {
  # Validation
  stopifnot(
    is.function(loglike),
    is.function(score),
    is.function(fisher),
    is.numeric(theta0),
    inherits(config, "mle_config"),
    inherits(constraint, "mle_constraint"),
    is.logical(inverted)
  )

  if (!constraint$support(theta0)) {
    stop("Initial guess `theta0` not in support")
  }

  # Build covariance and FIM functions based on inverted flag
  if (inverted) {
    covar_fn <- fisher
    fim_fn <- function(x) MASS::ginv(covar_fn(x))
  } else {
    fim_fn <- fisher
    covar_fn <- function(x) MASS::ginv(fim_fn(x))
  }

  # Newton direction: covariance * score (H^{-1} * grad)
  direction_fn <- function(theta) {
    as.vector(covar_fn(theta) %*% score(theta))
  }

  # Newton-Raphson typically benefits from line search
  use_linesearch <- inherits(config, "mle_config_linesearch")

  result <- .mle_optimize_direction(
    loglike = loglike,
    direction_fn = direction_fn,
    theta0 = theta0,
    config = config,
    constraint = constraint,
    use_linesearch = use_linesearch
  )

  # Augment with method-specific information
  result$score <- score(result$theta.hat)
  result$info <- fim_fn(result$theta.hat)
  result$sigma <- covar_fn(result$theta.hat)

  class(result) <- c("mle_newton_raphson", class(result))
  result
}
