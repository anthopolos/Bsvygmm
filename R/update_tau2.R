#' @title Update variance \code{tau2} of stratum level intercepts in the latent class membership model.
#'
#' @description This uses an inverse gamma prior or a uniform prior to update the variance of the area segment (aka cluster) level intercepts in the latent class membership model.
#' @param u A \code{K-1} column matrix of area segment level intercepts.
#' @param prior.shape A scalar for the prior shape, if an inverse gamma prior is used.
#' @param prior.scale A scalar for the prior scale, if an inverse gamma prior is used.
#' @param hierVar A string equal to \code{Unif} or \code{IG} to indicate a uniform prior distribution or inverse gamma prior distribution.
#' @export
#' @return A \code{K-1} vector of variances.
update_tau2 <- function(u, prior.shape, prior.scale, hierVar) {

  K <- ncol(as.matrix(u)) + 1

  if (K > 2) {

    values <- rep(NA, (K - 1))

    for (k in 1:(K - 1)) {
      # Scale matrix
      llik <- crossprod(u[ , k], u[ , k])
      # Degrees of freedom
      df <- nrow(u)

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
    llik <- crossprod(u, u)

    # Degrees of freedom
    df <- length(u)

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
