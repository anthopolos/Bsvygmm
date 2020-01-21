#' @title Update variance-covariance matrix of the subject level effects.
#'
#' @description An inverse Wishart prior is used to update the variance-covariance of the subject level effects.
#'
#' @param C Latent class assignments for each subject.
#' @param q A scalar for the number of subject effects.
#' @param b  A \code{q} column matrix of subject level effects.
#' @param prior.df A scalar for the prior degrees of freedom.
#' @param prior.scale A \code{q} by \code{q} prior scale matrix.
#' @export
#' @return Latent class specific variance-covariances.
update_Phi <- function(C, q, b, prior.df, prior.scale){

  # Number of classes
  K <- length(table(C))

  # Storage
  llik <- array(NA, dim = c(q, q, K))
  df <- rep(NA, K)

  for (k in 1:K) {

    # Subject level information
    ind_sub <- which(C == k)

    #Subject level random effects
    bk <- b[ind_sub, ] # Here we choose class k

    # Scale matrix
    llik[ , , k] <- crossprod(bk, bk)

    # Degrees of freedom
    df[k] <- length(ind_sub)

  }

  # Unstructured and latent class specific
  values <- array(NA, dim = c(q, q, K))

  for (k in 1:K) {

    # Posterior degrees of freedom
    post_df <- df[k] + prior.df
    # Posterior scale matrix
    post_scale <- llik[ , , k] + prior.scale # q x q

    values[ , , k] <- MCMCpack::riwish(v = post_df, S = post_scale)

  }

  return(values)

}

