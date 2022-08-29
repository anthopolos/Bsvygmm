#' @title Update latent variable \code{Z} in multinomial probit model.
#'
#' @description The latent variable in the multinomial probit model is updated using a truncated normal distribution. The updating algorithm allows for either or both area segment intercepts for correlations among subjects in the same area segment and spatial correlations among neighboring area segments. The number of latent classes is assumed to be greater than 2.
#' @param Z \code{K-1} column matrix of latent variables Z.
#' @param C Latent class assignments for each subject.
#' @param delta \code{K-1} column matrix of latent class-specific regression coefficients.
#' @param alpha \code{K-1} column matrix of regression coefficients for the B-spline.
#' @param deltaStratum A \code{K-1} column matrix of stratum specific effects. Equals \code{NULL} if not desired.
#' @param u A \code{K-1} column matrix of area segment specific effects for modeling correlations among subjects in the same area segment. Equals \code{NULL} if not desired.
#' @param nu A \code{K-1} column matrix of spatial random effects. Equals \code{NULL} if not desired.
#' @param W An \code{s} column design matrix of covariates, including a column of one's for intercept.
#' @param B An \code{R} column design matrix of basis functions for the B-spline.
#' @param clusterIDSub Area segment or cluster identifier for each subject.
#' @param stratumIDSub Stratum identifier for each subject.
#' @param spline Logical. If \code{TRUE}, a B-spline is included for the basis functions in argument \code{B} for the latent class membership model.
#' @export
#' @return A \code{K - 1} matrix of latent variable \code{Z}.
update_Z_MNP <- function(Z, C, delta, alpha, deltaStratum, u, nu, W, B, clusterIDSub, stratumIDSub, spline) {

  # Number of subjects
  n <- length(C)

  # Number of classes
  K <- ncol(delta) + 1

  # Transform C to multinomial binary matrix
  CMN <- matrix(0, nrow = n, ncol = K)
  for (i in 1:n) {
    for (k in 1:K) {
      CMN[i , C[i]] <- 1
    }
  }

  ### Lower, upper matrices for K - 1 columns
  lower <- upper <- matrix(NA, nrow = n, ncol = (K - 1))

  # The first class (k = 1) is the reference class
  for (i in 1:n) {
    for (k in 2:K) {
      if (CMN[i, k] == 1) {
        lower[i, (k - 1)] <- max(c(0, Z[i, -(k - 1)]))
        upper[i, (k - 1)] <- Inf
      } else if (CMN[i, k] == 0) {
        lower[i, (k - 1)] <- -Inf
        upper[i, (k - 1)] <- max(c(0, Z[i, -(k - 1)]))
      }
    }
  }


  if (!is.null(nu)) {
    nJ <- nrow(nu)
    #Only expand nu for clusterIDs that appear in data
    clusterIDSubf <- factor(clusterIDSub, levels = seq(1, nJ, by = 1))
  }


  mu <- matrix(NA, nrow = n, ncol = (K - 1))

  for (k in 1:(K - 1)) {

    # 1. Spatially structured heterogeneity at the cluster level only
    if (!is.null(nu) & is.null(u)) {
      if (spline) {
        mu[ , k] <- W %*% delta[ , k] + B %*% alpha[ , k] + nu[clusterIDSubf, k]
      } else {
        mu[ , k] <-  W %*% delta[ , k] + nu[clusterIDSubf, k]
      }

    }

    # 2. Unstructured heterogeneity at the cluster level only
    if (is.null(nu) & !is.null(u)) {
      if (spline) {
        mu[ , k] <- W %*% delta[ , k] + B %*% alpha[ , k] + u[factor(clusterIDSub), k]
      } else {
        mu[ , k] <- W %*% delta[ , k] + u[factor(clusterIDSub), k]
      }

    }


    # 3. Spatially structured and unstructured heterogeneity at the cluster level
    if (!is.null(nu) & !is.null(u)) {
      if (spline) {
        mu[ , k] <-  W %*% delta[ , k] + B %*% alpha[ , k] + u[factor(clusterIDSub), k] + nu[clusterIDSubf, k]
      } else {
        mu[ , k] <- W %*% delta[ , k] + u[factor(clusterIDSub), k] + nu[clusterIDSubf, k]
      }

    }

    # 4. No spatially structured and unstructured heterogeneity at the cluster level
    if (is.null(nu) & is.null(u)) {
      if (spline) {
        mu[ , k] <-  W %*% delta[ , k] + B %*% alpha[ , k]
      } else {
        mu[ , k] <- W %*% delta[ , k]
      }

    }

    # Add stratum level random effect if desired
    if (!is.null(deltaStratum)) {
      mu[ , k] <- mu[ , k] + deltaStratum[factor(stratumIDSub), k]
    }

    Z[ , k] <- truncnorm::rtruncnorm(n = n, mean = mu[ , k], sd = 1, a = lower[ , k], b = upper[ , k])

  } # End k loop


  return(Z)

}
