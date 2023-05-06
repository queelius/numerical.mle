#' mle_gradient_ascent
#' 
#' MLE method using gradient ascent
#'
#' @param theta0 initial guess of MLE
#' @param score score function, the gradient of a log-likelihood
#' @inherit mle_local_search options
#' 
#' @return an object of class `mle_gradient_ascent`, which is an `mle` object
#' @importFrom MASS ginv
#' @importFrom numDeriv hessian
#' @export
mle_gradient_ascent <- function(theta0, score, options) {
    if (!is.function(score)) {
        stop("score must be a function, the gradient of a log-likelihood")
    }
    sol <- mle_local_search(theta0 = theta0, dir = score, options = options)
    class(sol) <- c("mle_gradient_ascent", class(sol))
    if (options$loglike) {
        sol$info <- -hessian(options$loglike, sol$theta.hat)
        sol$sigma <- ginv(sol$info)
    }
    sol$score <- score(sol$theta.hat)
    sol
}
