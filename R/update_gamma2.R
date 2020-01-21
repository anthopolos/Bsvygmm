#' @title Update variance \code{gamma2} of stratum level intercepts in the latent class membership model.
#'
#' @description This uses an inverse gamma prior or a uniform prior to update the variance of the stratum level intercepts in the latent class membership model.
#' @param deltaStratum A \code{K-1} column matrix of stratum level intercepts.
#' @param prior.shape A scalar for the prior shape, if an inverse gamma prior is used.
#' @param prior.scale A scalar for the prior scale, if an inverse gamma prior is used.
#' @param hierVar A string equal to \code{Unif} or \code{IG} to indicate a uniform prior distribution or inverse gamma prior distribution.
#' @export
#' @return A \code{K-1} vector of variances.
update_gamma2 <- function(deltaStratum, prior.shape, prior.scale, hierVar) {

  K <- ncol(as.matrix(deltaStratum)) + 1

  if (K > 2) {

    values <- rep(NA, (K - 1))

    for (k in 1:(K - 1)) {

      # Scale matrix
      llik <- crossprod(deltaStratum[ , k], deltaStratum[ , k])

      # Degrees of freedom
      df <- nrow(deltaStratum)

      if (hierVar == "IG") {
        # Posterior scale
        post_scale <-  llik / 2 + prior.scale # q x 1
        # Posterior shape
        post_shape <- df / 2 + prior.shape # A scalar
      } else if (hierVar == "Unif") {
        post_scale <-  llik / 2
        # Posterior shape
        post_shape <- (df - 1) / 2  # A scalar
      }

      values[k] <- 1 / rgamma(1, post_shape, post_scale)
    }

  } else if (K == 2) {

    # Scale matrix
    llik <- crossprod(deltaStratum, deltaStratum)

    # Degrees of freedom
    df <- length(deltaStratum)

    if (hierVar == "IG") {
      # Posterior scale
      post_scale <-  llik / 2 + prior.scale # q x 1
      # Posterior shape
      post_shape <- df / 2 + prior.shape # A scalar
    } else if (hierVar == "Unif") {
      post_scale <-  llik / 2
      # Posterior shape
      post_shape <- (df - 1) / 2  # A scalar
    }

    values <- 1 / rgamma(1, post_shape, post_scale)

  }

  return(values)

}
