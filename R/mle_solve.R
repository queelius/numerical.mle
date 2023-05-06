#' mle_solve
#' 
#' Uses various functions and algorithms to find the MLE.
#' 
#' @importFrom numDeriv grad
#' 
#' @param theta0 numeric vector, initial guess for MLE
mle_solve <- function(options, theta0, method = c("optim"), ...) {

    for (m in method)
    {
        mle <- NULL
        if (method == "newton_raphson") {

            if (is.null(options$score)) {
                options$score <- function(x) {
                    numDeriv::grad(options$loglike, x)
                }
            }
            if (is.null(options$info)) {
                options$fim <- function(x) {
                    -numDeriv::hessian(options$loglike, x)
                }
            }

            mle <- mle_newton_raphson(
                theta0,
                options$score,
                options$info,
                options)

        } else if (method == "optim") {

            sol <- optim(
                par = theta0,
                fn = options$loglike,
                gr = options$score,
                hessian = TRUE, ...)

            mle <- mle_optim(sol)

        } else if (method == "gradient_ascent") {

            if (is.null(options$score)) {
                options$score <- function(x) numDeriv::grad(options$loglike, x)
            }

            mle <- mle_gradient_ascent(
                theta0,
                options$score,
                options)

        } else if (method == "sim_anneal") {

            mle <- mle_sim_anneal(
                theta0 = theta0,
                options = options,
                loglike = options$loglike)
    
        } else {
            stop("unknown method")
        }

        theta0 <- point(mle)
    }

    mle
}

