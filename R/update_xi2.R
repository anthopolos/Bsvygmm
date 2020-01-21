#' @title Update variance \code{xi2} of area segment level spatial intercepts in the latent class membership model.
#'
#' @description This uses an inverse gamma prior or a uniform prior to update the variance of the area segment (aka cluster) level spatial intercepts in the latent class membership model.
#' @param nu A \code{K-1} column matrix of area segment level spatial intercepts.
#' @param ADJ A binary symmetric adjacency matrix. Equal to \code{NULL} if spatial correlations are not included in the latent class membership model.
#' @param prior.shape A scalar for the prior shape, if an inverse gamma prior is used.
#' @param prior.scale A scalar for the prior scale, if an inverse gamma prior is used.
#' @param hierVar A string equal to \code{Unif} or \code{IG} to indicate a uniform prior distribution or inverse gamma prior distribution.
#' @export
#' @return A \code{K-1} vector of variances.
update_xi2 <- function(nu, ADJ, prior.shape, prior.scale, hierVar) {

  K <- ncol(as.matrix(nu)) + 1
  nJ <- dim(ADJ)[1]
  m <- rowSums(ADJ)

  if (K > 2) {

    values <- rep(NA, (K - 1))

    for (k in 1:(K - 1)) {

      ### Posterior scale
      # A vector of conditional means for nu for class k
      mu <- apply(ADJ == 1, 1, function(x){mean(nu[x, k])})
      llik <- sum(m*(nu[ , k]^2 - 2*nu[ , k]*mu + mu^2))
      df <- nJ

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

  }

  if (K == 2) {

    ### Posterior scale
    # A vector of conditional means for nu for class k
    mu <- apply(ADJ == 1, 1, function(x){mean(nu[x])})

    llik <- sum(m*(nu^2 - 2*nu*mu + mu^2))
    df <- nJ

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
