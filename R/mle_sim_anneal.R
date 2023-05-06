#' mle_sim_anneal
#'
#' This function takes the output of `sim_anneal` and turns it into an `mle` object.
#' 
#' @importFrom algebraic.mle mle
#' 
#' @export
mle_sim_anneal <- function(theta0, loglike = NULL, options = list()) {
    if (!is.null(loglike)) {
        options$loglike <- loglike
    }
    sol <- sim_anneal(x0 = theta0, obj_fn = options$loglike, options = options)
    mle.sol <- mle(
        theta.hat = sol$argmax,
        loglike = sol$max,
        score = NULL,
        sigma = NULL,
        info = NULL,
        obs = NULL,
        nobs = NULL,
        superclasses = c("mle_sim_anneal", "numerical_mle"))
    mle.sol$converged <- TRUE
    mle.sol
}



#' sim_anneal
#'
#' This function implements the simulated annealing algorithm,
#' which is a global optimization algorithm that is useful for
#' finding a good starting point for a local optimization algorithm.

#' We do not return this as an MLE object because, to be a good
#' estimate of the MLE, the gradient of `f` evaluated
#' at its solution should be close to zero, assuming the MLE
#' is interior to the domain of `f`. However, since this algorithm
#' is not guided by gradient information, it is not sensitive to
#' the gradient of `f` and instead only seeks to maximize `f`.
#' 
#' @param obj_fn Objective function to maximize, default is NULL
#'               (must be specified in options)
#' @param x0 Initial guess, default is NULL (must be specified in options)
#' @param options List of optional arguments
#' @param ... Additional arguments that may be passed to `options$neigh`
#' @describeIn sim_anneal options
#' @field t_init Initial temperature
#' @field t_end Final temperature
#' @field alpha Cooling factor
#' @field iter_per_temp Number of iterations per temperature
#' @field max_iter Maximum number of iterations, used instead of t_end
#'        if not NULL, defaults to NULL
#' @field debug If TRUE, print debugging information to the console
#' @field trace If TRUE, track the history of positions and values
#' @field sup Support function, returns TRUE if x is in the domain of f
#' @field neigh Neighborhood function, returns a random neighbor of x
#' @field debug_freq Frequency of debug output, defaults to 10
#' @field obj_fn Objective function to maximize, if not specified in
#'        formal parameter `obj_fn`
#' @field x0 Initial guess, if not specified in formal parameter `x0`
#' @return list with best solution (argmax) and its corresponding
#'         objective function value (max), and optionally path
#' @importFrom stats runif
#' @export
sim_anneal <- function(x0 = NULL, obj_fn = NULL, options = list(), ...) {

    defaults <- list(
        t_init = 100,
        t_end = 1e-6,
        alpha = 0.95,
        iter_per_temp = 100,
        max_iter = 1e20,
        debug = FALSE,
        trace = FALSE,
        sup = function(x) TRUE,
        neigh = function(x) x + rnorm(length(x)),
        debug_freq = 10)

    if (!is.null(obj_fn)) {
        options$obj_fn <- obj_fn
    }
    if (!is.null(x0)) {
        options$x0 <- x0
    }
    x0 <- options$x0
    f <- options$obj_fn

    options <- modifyList(defaults, options)

    stopifnot(options$t_init > 0,
              options$alpha > 0, options$alpha < 1,
              options$iter_per_temp > 0,
              (is.null(options$max_iter) || options$max_iter > 0),
              is.function(options$neigh),
              is.function(options$sup),
              is.logical(options$debug),
              is.logical(options$trace))

    if (!is.function(f)) {
        stop("obj_fn must be a function (to maximize)")
    }
    if (!options$sup(x0)) {
        stop("initial guess x0 not in support")
    }

    argmax <- x0
    fmax <- f(argmax)
    fx0 <- fmax
    t <- options$t_init
    iter <- 0L
    path <- matrix(nrow=0,ncol=length(x0))

    while (t > options$t_end) {
        iter <- iter + 1L
        x <- options$neigh(x0, ...)

        if (iter %% options$debug_freq == 0) {
            cat(sprintf("| %d | (%s) -> (%s) |\n",
                iter,
                paste(round(x0,4), collapse=", "),
                paste(round(x0,4), collapse=", ")))
        }
        if (!options$sup(x)) {
            if (options$debug) {
                cat("(", x, ") not in support, skipping\n")
            }
            next
        }
        fx <- f(x)
        if (is.nan(fx)) {
            if (options$debug) {
                cat("obj_fn at (", x, ") is NaN, skipping\n")
            }
            next
        }

        if (fx > fmax) {
            if (options$debug) {
                cat(sprintf("Found obj_fn = %.4f at solution = (%s)\n",
                    fx, paste(round(x,4), collapse=", ")))
            }
            argmax <- x
            fmax <- fx
        }

        if (exp((fx - fx0) / t) > runif(1)) {
            if (options$debug) {
                cat(sprintf("| %d | %.4f | (%s) |\n",
                iter, fx, paste(round(x,4), collapse=", ")))
            }
            x0 <- x
            fx0 <- fx
            if (options$trace) {
                path <- rbind(path, x)
            }
        }
        if (iter %% options$iter_per_temp == 0) {
            t <- t * options$alpha
        }
    }

    sol <- list(argmax = argmax, max = fmax, options = options)
    if (options$trace) {
        sol$path <- path
    }
    sol
}
