#' mle_local_search
#'
#' Performs a local search to find the MLE, assuming the MLE is an
#' interior point of the support and that an initial guess `theta0`
#' that is near the MLE is provided. Use a global search method
#' like `sim_anneal` to find a good initial guess.
#'
#' @param theta0 numeric, initial guess
#' @param dir function, promising direction function
#' @param options list, options for the local search, see function description.
#' @describeIn mle_local_search options
#' @field sup function, domain of support for log-likelihood
#' @field eta numeric, learning rate, defaults to 1
#' @field max_iter integer, maximum number of iterations, defaults to 1000
#' @field max_iter_ls integer, maximum number of iterations for the
#'        line search, defaults to 1000
#' @field abs_tol numeric, tolerance for convergence, defaults to NULL
#'        (use rel_tol instead)
#' @field rel_tol numeric, relative tolerance for convergence, defaults to 1e-5
#' @field r numeric, backtracking line search parameter, defaults to 0.5
#' @field proj function, projection function to enforce domain of support
#' @field norm function, distance measure for convergence checks, defaults
#'        to the to the infinity norm.
#' @field debug logical, output debugging information if TRUE; default FALSE
#' @field trace logical, if TRUE store the path of the search in the `path`
#'        attribute of the output; default FALSE
#' @field line_search logical, if TRUE, perform a line search; default TRUE
#'        in this case, learning rate `eta` refers to the maximum step size
#'        that can be taken per iteration.
#' @field debug_freq integer, frequency of debug output, defaults to 1
#' 
#' @return an `mle` object with additional attributes `iter` and `converged`
#'         and optionally `path` if `trace` is TRUE.
#' @export
mle_local_search <- function(
    dir,
    theta0,
    loglike = NULL,
    options = list())
{
    defaults <- list(
        loglike = NULL,
        abs_tol = NULL,
        rel_tol = 1e-5,
        proj = function(theta) theta,
        sup = function(theta) TRUE,
        norm = function(dx) max(abs(dx)),
        eta = 1,
        r = 0.5,
        max_iter = 100L,
        max_iter_ls = 10L,
        debug = FALSE,
        trace = FALSE,
        line_search = FALSE,
        debug_freq = 1L
    )
    options <- modifyList(defaults, options)
    stopifnot(
        is.function(dir),
        options$eta > 0,
        options$r > 0, options$r < 1,
        options$max_iter > 0,
        options$max_iter_ls > 0,
        is.logical(options$debug),
        is.logical(options$trace),
        (is.null(options$loglike) || is.function(options$loglike)),
        is.function(options$proj),
        is.function(options$norm),
        is.function(options$sup))

    if (options$line_search && is.null(options$loglike)) {
        stop("Line search requires log-likelihood function (options$loglike)")
    }

    if (!options$sup(theta0)) {
        stop("Initial guess `theta0` not in support")
    }

    check_converged <- NULL
    if (is.null(options$abs_tol)) {
        stopifnot(!is.null(options$rel_tol), options$rel_tol > 0)
        check_converged <- function(x1,x0) {
            options$norm(x1-x0) < options$rel_tol * options$norm(x1)
        }
    } else {
        stopifnot(options$abs_tol > 0)
        check_converged <- function(x1,x0) {
            options$norm(x1-x0) < options$abs_tol
        }
    }

    converged <- FALSE
    path <- NULL
    if (options$trace) {
        path <- matrix(nrow=options$max_iter, ncol=length(theta0))
    }

    iter <- 1L
    repeat {
        d <- dir(theta0)
        if (options$debug && iter %% options$debug_freq == 0) {
            if (is.null(options$loglike)) {
                cat(sprintf("| %-4d | Direction: (%s) | Solution: %s |\n",
                    iter, toString(round(d,3)), toString(round(theta0,3))))
            } else {
                cat(sprintf("| %-8d | Direction: (%s) | Solution: %s | Log-likelihood: %.5f |\n", # nolint: line_length_linter.
                    iter, toString(round(d,3)), toString(round(theta0,3)),
                    options$loglike(theta0),5))
            }
        }

        theta1 <- NULL
        if (options$line_search) {
            # we use backtracking for an approximate line search
            res <- backtracking_line_search(
                f = options$loglike,
                dir = d,
                x0 = theta0,
                max_step = options$eta,
                sup = options$sup,
                fix = options$proj,
                debug = options$debug,
                r = options$r,
                max_iter = options$max_iter_ls)

            if (!res$found_better) {
                if (options$debug) {
                    cat("Backtracking failed to find a better solution\n")
                }
                break
            }
            theta1 <- res$argmax            
        } else {
            theta1 <- theta0 + options$eta * d
            if (!options$sup(theta1)) {
                if (options$debug) {
                    cat("Projecting `theta1` onto support\n")
                }
                theta1 <- options$proj(theta1)
            }
        }

        if (options$trace) {
            path[iter,] <- as.vector(theta1)
        }

        good_enuf <- check_converged(theta1, theta0)
        theta0 <- theta1
        if (good_enuf) {
            converged <- TRUE
            break
        } else if (iter == options$max_iter) {
            break
        }

        iter <- iter + 1L
    }

    if (!converged && options$debug) {
        cat("Failed to converge\n")
    }

    loglike_value <- NULL
    if (!is.null(options$loglike)) {
        loglike_value <- options$loglike(theta0)
    }

    sol <- mle(theta.hat = theta0, loglike = loglike_value,
        superclasses = c("mle_local_search", "mle_numerical"))
    sol$iter <- iter
    sol$converged <- converged
    sol$options <- options

    if (options$trace) {
        sol$path <- path[1:iter, , drop = FALSE]
    }
    sol
}
