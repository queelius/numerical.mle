#' mle_newton_raphson
#' 
#' Find an MLE using the Newton-Raphson method.
#'
#' @param score score Function, the gradient of log-likelihood
#' @param fim FIM function
#' @param theta0 initial guess of MLE
#' @param inverted logical, if TRUE `fim` is covariance instead of FIM
#' @inherit mle_local_search options
#'
#' @return an object of class `mle_newton_raphson`, which is an `mle` object
#' @importFrom MASS ginv
#' @export
mle_newton_raphson <- function(
    score,
    fim,
    theta0,
    inverted = FALSE,
    options = list())
{
    if (!is.function(score)) {
        stop("`score` must be a function, the gradient of a log-likelihood")
    } else if (!is.function(fim)) {
        stop("`fim` must be a function, a FIM or covariance matrix (inverted)")
    } else if (!is.logical(inverted)) {
        stop("inverted must be a logical, if TRUE `fim` is covariance instead of FIM") # nolint: line_length_linter.
    }
    covar <- NULL
    if (inverted) {
        covar <- fim
        fim <- function(x) ginv(covar(x))
    } else {
        covar <- function(x) ginv(fim(x))
    }
    dir <- function(x) covar(x) %*% score(x)

    sol <- mle_local_search(theta0 = theta0, dir = dir, options = options)
    class(sol) <- c("mle_newton_raphson", class(sol))
    sol$score <- score(sol$theta.hat)
    sol$sigma <- covar(sol$theta.hat)
    sol$info <- fim(sol$theta)
    sol
}
