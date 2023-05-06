


#' mle_random_search
#' 
#' MLE method using random search
#' @export
mle_random_search <- function(rtheta, loglike, options)
{
    if (!is.function(loglike)) {
        stop("loglike must be a function, the log-likelihood")
    }

    defaults <- list(
        sup = function(x) TRUE,
        loglike = NULL,
        mle_solver = NULL,
        rtheta = NULL,
        trials = 100
    )
    options <- modifyList(defaults, options)

    mle <- options$mle_solver(rtheta(), options)
    for (i in 1:options$trials) {
        theta0 <- options$rtheta()
        sol <- options$mle_solver(theta0, options)
        if (is.null(sol)) {
            next
        }
        if (loglike(sol) > loglike(mle)) {
            mle <- sol
        }
    }

    mle
}