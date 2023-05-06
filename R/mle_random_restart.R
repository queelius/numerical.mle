
#' mle_random_restart
#' 
#' MLE method using random restart
#' 
#' An interesting idea is to use random restarts to find the MLE. This is
#' similar to the idea of simulated annealing, but instead of using a
#' temperature parameter, we just try a bunch of random starting and then try
#' to find the MLE from each of those starting points using `mle_solver`. The
#' best MLE is returned.
#' 
#' A reasonable `rtheta0` is a multivariate normal distribution with a wide
#' variance. This is because we want to explore the entire parameter space.
#' 
#' If we already have a candidate MLE `x`, we can use that as the sampling distribution of
#' the MLE (MVN with mean `point(x)` and variance-covariance `vcov(x)`). This
#' will generate random starting points that are close to the candidate MLE, which
#' may help us find a better MLE.
#' 
#' @param rtheta0 initial guess of MLE distribution, defaults to MVN with wide variance
#' @importFrom algebraic.mle loglike
#' @export
mle_random_restart <- function(
    rtheta0,
    options)
{
    defaults <- list(
        sup = function(x) TRUE,
        trials = 100,
        loglike = NULL,
        mle_solver = NULL,
        rtheta0 = NULL)

    options <- modifyList(defaults, options)

    mle <- options$mle_solver(theta0, options)
    for (i in 1:options$trials) {
        theta0 <- options$rtheta0()
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
