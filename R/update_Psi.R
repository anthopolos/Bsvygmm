#' @title Update stratum level variance in the longitudinal outcomes model.
#'
#' @description An inverge gamma or uniform prior is used to update the variance of the stratum level intercepts.
#' @param C Latent class assignments for each subject.
#' @param zeta  A \code{K} column matrix of intercepts at the stratum level.
#' @param prior.scale A scalar for the prior scale, if an inverse gamma prior is used.
#' @param prior.shape A scalar for the prior shape, if an inverse gamma prior is used.
#' @param hierVar A string equal to \code{Unif} or \code{IG} to indicate a uniform or inverse gamma prior respectively.
#' @export
#' @return A \code{K} length vector of latent class specific variances.
update_Psi <- function(C, zeta, prior.scale, prior.shape, hierVar){

  # Number of classes
  K <- length(table(C))

  # Storage
  llik <- rep(NA, K)
  df <- rep(NA, K)

  # Latent class specific variance
  values <- rep(NA, K)

  for (k in 1:K) {

    # Cluster level random effects
    zetak  <- zeta[ , k]

    # Scale matrix
    llik <- crossprod(zetak, zetak)

    # Degrees of freedom
    df <- length(zetak)

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

  return(values)

}


