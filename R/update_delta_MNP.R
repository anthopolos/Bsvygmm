#' @title Update regression coefficients \code{delta} in the multinomial probit model.
#'
#' @description This uses a normal prior distribution to update the regression coefficients in the multinomial probit model. The regression coefficients for \code{k=1} are fixed to zero.
#' @param Z A \code{K-1} column matrix of latent variable \code{Z}.
#' @param delta A \code{K-1} column matrix of regression coefficients in the multinomial probit model with \code{k=1} fixed to zero.
#' @param alpha A \code{K-1} column matrix of regression coefficients for the B-spline.
#' @param deltaStratum A \code{K-1} column matrix of stratum specific effects.
#' @param u A \code{K-1} column matrix of area segment specific effects for modeling correlations among subjects in the same area segment. Equals \code{NULL} if not desired.
#' @param nu A \code{K-1} column matrix of spatial random effects. Equals \code{NULL} if not desired.
#' @param W An \code{s} column design matrix of covariates, including a column of one's for intercept.
#' @param B An \code{R} column design matrix of basis functions for the B-spline.
#' @param clusterIDSub Area segment or cluster identifier for each subject.
#' @param stratumIDSub Stratum identifier for each subject.
#' @param prior.mu A vector of prior means supplied through \code{priors}.
#' @param prior.Sigma A covariance matrix supplied through \code{priors}.
#' @param clusterIDSub Area segment or cluster identifier for each subject.
#' @param stratumIDSub Stratum identifier for each subject.
#' @param spline Logical. If \code{TRUE}, a B-spline is included for the basis functions in argument \code{B} for the latent class membership model.
#' @export
#' @return A \code{K-1} column matrix of regression coefficients in the multinomial probit model with \code{k=1} fixed to zero.
update_delta_MNP <- function(Z, delta, alpha, deltaStratum, u, nu, W, B, prior.mu, prior.Sigma, clusterIDSub, stratumIDSub, spline) {

  # Number of classes
  K <- ncol(delta) + 1

  # Number of subjects
  n <- nrow(Z)


  if (!is.null(nu)) {
    nJ <- nrow(nu)
    #Only expand nu for clusterIDSubs that appear in data
    clusterIDSubf <- factor(clusterIDSub, levels = seq(1, nJ, by = 1))
  }

  # Update delta with proper conjugate, note the error distribution on latent variable is N(X\beta, 1) so posterior variance of delta does not get updated
  #See Chib 1993
  post_var <- solve(t(W) %*% W + solve(prior.Sigma))

  mu <- matrix(NA, nrow = n, ncol = (K - 1))

  for (k in 1:(K - 1)) {

    # 1. Spatially structured heterogeneity at the cluster level only
    if (!is.null(nu) & is.null(u)) {
      if (spline) {
        mu[ , k] <- Z[ , k] - B %*% alpha[ , k] - nu[clusterIDSubf, k] - deltaStratum[factor(stratumIDSub), k]
      } else {
        mu[ , k] <-  Z[ , k] - nu[clusterIDSubf, k] - deltaStratum[factor(stratumIDSub), k]
      }

    }

    # 2. Unstructured heterogeneity at the cluster level only
    if (is.null(nu) & !is.null(u)) {
      if (spline) {
        mu[ , k] <- Z[ , k] - B %*% alpha[ , k] - u[factor(clusterIDSub), k] - deltaStratum[factor(stratumIDSub), k]
      } else {
        mu[ , k] <- Z[ , k] - u[factor(clusterIDSub), k] - deltaStratum[factor(stratumIDSub), k]
      }

    }

    # 3. Spatially structured and unstructured heterogeneity at the cluster level only
    if (!is.null(nu) & !is.null(u)) {
      if (spline) {
        mu[ , k] <-  Z[ , k] - B %*% alpha[ , k] - u[factor(clusterIDSub), k] - deltaStratum[factor(stratumIDSub), k] - nu[clusterIDSubf, k]
      } else {
        mu[ , k] <- Z[ , k] - u[factor(clusterIDSub), k] - deltaStratum[factor(stratumIDSub), k] - nu[clusterIDSubf, k]
      }

    }


    post_mean <- post_var %*% (t(W) %*% mu[ , k] + solve(prior.Sigma) %*% prior.mu)

    delta[ , k] <- mnormt::rmnorm(1, mean = post_mean, varcov = post_var)

  } # End of k loop

  return(delta)

}
