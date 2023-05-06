


#' mle_grid_search
#' 
#' MLE method using grid search
#' @export

mle_grid_search <- function(loglike, options) {

    if (!is.function(loglike)) {
        stop("loglike must be a function, the log-likelihood")
    }

    defaults <- list(
        sup = function(x) TRUE,
        loglike = NULL,
        mle_solver = NULL,
        lower <- NULL,
        upper <- NULL
    )
    options <- modifyList(defaults, options)

    # subdivide parameter support into a grid each of size `grid_size` (for
    # along each dimension)
    grid <- subdivide_region(options$lower, options$upper, options$grid_size)

    for (i in 1:nrow(grid)) {
        theta0 <- grid[i, ]
        sol <- options$mle_solver(theta0, options)
        if (is.null(sol)) {
            next
        }
        if (loglike(sol) > loglike(mle)) {
            mle <- sol
        }
    }
}


#' subdivide_region
#' @param lower lower bounds of parameter support (vector)
#' @param upper upper bounds of parameter support (vector)
#' @param grid_size maximum size of each dimension of each hypercube
#'                  that makes up the grid
subdivide_region <- function(lower, upper, grid_size) {

    if (length(lower) != length(upper)) {
        stop("lower and upper must be the same length")
    }

    if (any(lower > upper)) {
        stop("lower must be less than (or equal to) upper")
    }

    if (any(grid_size <= 0)) {
        stop("grid_size must be greater than 0")
    }

    # number of dimensions
    n <- length(lower)

    # number of hypercubes in each dimension
    n_hypercubes <- ceiling((upper - lower) / grid_size)

    # number of hypercubes in total
    n_total <- prod(n_hypercubes)

    # initialize grid
    grid <- matrix(NA, nrow = n_total, ncol = n)

    # initialize hypercube index
    hypercube <- rep(1, n)

    # initialize grid index
    i <- 1

    # loop over hypercubes
    while (hypercube[1] <= n_hypercubes[1]) {

        # loop over dimensions
        for (j in 1:n) {

            # set grid value
            grid[i, j] <- lower[j] + (hypercube[j] - 1) * grid_size[j]

            # increment hypercube index
            hypercube[j] <- hypercube[j] + 1

            # if hypercube index exceeds number of hypercubes in dimension
            if (hypercube[j] > n_hypercubes[j]) {

                # reset hypercube index
                hypercube[j] <- 1
            } else {

                # break out of loop over dimensions
                break
            }
        }

        # increment grid index
        i <- i + 1
    }

    # return grid
    grid
}