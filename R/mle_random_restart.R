
#' MLE method using random restart
#' 
#' An interesting idea is to use random restarts to find the MLE. This is
#' similar to the idea of simulated annealing, but instead of using a
#' temperature parameter, we just try a bunch of random starting points and then
#' find a local maximum from each of those starting points using `mle_solver`.
#' The best solution is returned.
#' 
#' @param rtheta0 initial guess sampler
#' @param mle_solver function to find MLE
#' @param ntrials number of random restarts to try
#' @param ... additional arguments to pass to `mle_solver`
#' @importFrom algebraic.mle loglike
#' @export
mle_random_restart <- function(
    rtheta0,
    mle_solver,
    ntrials = 100L,    
    ...)
{
    best_sol <- NULL
    for (i in 1:ntrials) {
        sol <- mle_solver(rtheta0(), ...)
        if (is.null(sol)) {
            next
        }
        if (is.null(best_sol) || loglik_val(sol) > loglik_val(best_sol)) {
            best_sol <- sol
        }
    }
    best_sol
}
