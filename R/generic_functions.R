#' is_converged
#'
#' Function to determine whether a `mle_numerical` object has converged.
#'
#' @param x the `mle` object
#' @param ... additional arguments to pass
#' @export
is_converged <- function(x, ...) {
    UseMethod("is_converged")
}

#' @export
is_converged.mle_numerical <- function(x, ...) {
    x$converged
}

#' is_mle_numerical
#' 
#' Function to determine whether an object `x` is of type `mle_numerical`.
#'
#' @param x the `mle` object
#' @param ... additional arguments to pass
#' @export
is_mle_numerical <- function(x) {
    inherits(x, "mle_numerical")
}

#' num_iterations
#' @param x the `mle` object
#' @param ... additional arguments to pass
#' @return the number of iterations used to find the MLE
#' @examples
#' loglike <- function(theta) {
#'    -sum(dnorm(theta, mean = theta[1], sd = theta[2], log = TRUE))
#' }
#' score <- function(theta) {
#'   -numDeriv::grad(loglike, theta)
#' }
#' sol <- mle_gradient_raphson(theta0 = theta, score = score, loglike = loglike)
#' num_iterations(sol)
#' @export
num_iterations <- function(x, ...) {
    UseMethod("num_iterations")
}

#' @export
num_iterations.mle_numerical <- function(x, ...) {
    x$iter
}

#' mle_numerical
#' 
#' a constructor for the mle_numerical class.
#' @importFrom algebraic.mle mle
#' @export
mle_numerical <- function(theta.hat, loglike, score, info, sigma, iter, converged) {
    stopifnot(is.numeric(iter))
    stopifnot(is.logical(converged))

    sol <- mle(
        theta.hat = theta.hat,
        loglike = loglike,
        score = score,
        info = info,
        sigma = sigma,
        superclasses = c("mle_numerical"))

    sol$iter <- iter
    sol$converged <- converged
    sol
}

#' stochastic loglikelihood constructor
#' good for large datasets. if applied to a gradient ascent method, this
#' will perform stochastic gradient ascent.
#' @param log.p log pdf (or pmf) of the parametric model being fit to `obs`
#' parameters. it can also just be proportional to the log pdf, since sometimes
#' the normalizing constant is unknown or hard to compute.
#' @param obs a matrix, vector, or data frame of observations
#' @param options a list of options
#' @export
stochastic_loglike <- function(log_density, data, m, replace = FALSE)
{
    stopifnot(m <= length(data))
    stopifnot(is.function(log_density))

    function(theta) {
        sum(log_density(sample(x = data, size = m, replace = replace), theta))
    }
}