#' mle_optim
#'
#' This function takes the output of `optim` and turns it into an `mle` object.
#' 
#' @param sol the output of `optim`
#' @return a `numerical_mle` object, specialized for `optim` (stats package)
#' solutions.
#' @importFrom MASS ginv
#' @export
mle_optim <- function(sol) {
    sigma <- NULL
    info <- NULL
    if (!is.null(sol$hessian)) {
        info <- -sol$hessian
        sigma <- ginv(info)
    }
    mle.sol <- mle(
        theta.hat=sol$par,
        loglike=sol$value,
        score=NULL,
        sigma=sigma,
        info=info,
        obs=NULL,
        nobs=NULL,
        superclasses=c("mle_optim","numerical_mle"))
    mle.sol$converged <- sol$convergence == 0
    mle.sol$optim_data <- sol
    mle.sol
}
