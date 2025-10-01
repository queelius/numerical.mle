# Mock algebraic.mle::mle function for testing
# This is a minimal implementation to allow tests to run without the full algebraic.mle package

if (!require("algebraic.mle", quietly = TRUE)) {
  # Create a simple mle constructor mock
  mle <- function(theta.hat, loglike = NULL, score = NULL, info = NULL,
                  sigma = NULL, obs = NULL, nobs = NULL, superclasses = NULL) {
    result <- list(
      theta.hat = theta.hat,
      loglike = loglike,
      score = score,
      info = info,
      sigma = sigma,
      obs = obs,
      nobs = nobs
    )

    classes <- c(superclasses, "mle")
    class(result) <- classes

    return(result)
  }

  # Export the mock function
  assign("mle", mle, envir = .GlobalEnv)
}
