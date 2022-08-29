#' @title Update spatial area segment (aka cluster) level effects \code{nu} in the latent class membership model.
#'
#' @description The area segment (aka cluster) level effects that model spatial correlations among neighboring area segments are updated using an intrinsic conditional autoregressive prior distribution..
#' @param Z \code{K-1} column matrix of latent variables Z.
#' @param delta \code{K-1} column matrix of latent class-specific regression coefficients.
#' @param alpha \code{K-1} column matrix of regression coefficients for the B-spline.
#' @param deltaStratum A \code{K-1} column matrix of stratum specific effects. Equals \code{NULL} if not desired.
#' @param u A \code{K-1} column matrix of area segment specific effects for modeling correlations among subjects in the same area segment. Equals \code{NULL} if not desired.
#' @param nu A \code{K-1} column matrix of spatial random effects.
#' @param xi2 A \code{K-1} vector of variances for the arae segment level intercepts.
#' @param W An \code{s} column design matrix of covariates, including a column of one's for intercept.
#' @param B An \code{R} column design matrix of basis functions for the B-spline.
#' @param ADJ A binary symmetric adjacency matrix. Equal to \code{NULL} if spatial correlations are not included in the latent class membership model.
#' @param clusterIDSub Area segment or cluster identifier for each subject.
#' @param stratumIDSub Stratum identifier for each subject.
#' @param spline Logical. If \code{TRUE}, a B-spline is included for the basis functions in argument \code{B} for the latent class membership model.
#' @export
#' @return A \code{K - 1} matrix of spatial intercepts at the area segment level.
update_nu <- function(Z, delta, alpha, deltaStratum, u, nu, xi2, W, B, ADJ, clusterIDSub, stratumIDSub, spline) {

  K <- ncol(as.matrix(delta)) + 1

  ### Spatial information
  # Number of spatial units in the study area
  nJ <- nrow(ADJ)

  # Redefine clusterIDSub as a factor that can miss certain levels
  #In this way we can get a 0 llik contribution to spatial units in which no subjects live
  clusterIDSubf <- factor(clusterIDSub, levels = seq(1, nJ, by = 1))

  # Number of neighbors for each spatial unit
  m <- rowSums(ADJ)

  # Number of subjects in each cluster, defining clusterIDSub as a factor with the complete set of spatial units as levels
  nis <- as.vector(table(clusterIDSubf))

  if (K > 2) {
    ### Store values
    values <- matrix(NA, nrow = nJ, ncol = (K - 1))

    ### Loop over classes
    for (k in 1:(K - 1)) {

      # Choose among mean components, including cluster level unstructured heterogeneity and spline
      if (!is.null(u)) {
        if (spline) {
          mu <- Z[ , k] - W %*% delta[ , k] - B %*% alpha[ , k] - u[factor(clusterIDSub), k]
        } else {
          mu <- Z[ , k] - W %*% delta[ , k] - u[factor(clusterIDSub), k]
        }

      }

      if (is.null(u)) {
        if (spline) {
          mu <- Z[ , k] - W %*% delta[ , k] - B %*% alpha[ , k]
        } else {
          mu <- Z[ , k] - W %*% delta[ , k]
        }
      }

      # Stratum level effect if desired
      if (!is.null(deltaStratum)) {
        mu <- mu - deltaStratum[factor(stratumIDSub), k]
      }

      # Prior mean for all nu j
      mu_nu <- apply(ADJ == 1, 1, function(x){mean(nu[x, k])})

      # Loop over spatial units
      for (j in 1:nJ) {
        # Posterior variance
        post_var <- solve(m[j] / xi2[k] + nis[j])

        # Posterior mean likelihood and prior contributions, variance of sampling distribution fixed to 1
        lik <- sum(mu[clusterIDSubf == j])

        prior <- m[j] / xi2[k] * mu_nu[j]
        post_mean <- post_var * (lik + prior)

        values[j , k] <- rnorm(1, mean = post_mean, sd = sqrt(post_var))
      }

      # Center on the fly
      values[ , k] <- values[ , k] - mean(values[ , k])

    }

  }

  if (K == 2) {

    ### Store values
    values <- rep(NA, nJ)

    ### Choose among mean components, including cluster level unstructured heterogeneity and spline
    if (!is.null(u)) {
      if (spline) {
        mu <- Z - W %*% delta - B %*% alpha - u[factor(clusterIDSub)]
      } else {
        mu <- Z - W %*% delta - u[factor(clusterIDSub)]
      }

    }

    if (is.null(u)) {
      if (spline) {
        mu <- Z - W %*% delta - B %*% alpha
      } else {
        mu <- Z - W %*% delta
      }
    }

    if (!is.null(deltaStratum)) {
      mu <- mu - deltaStratum[factor(stratumIDSub)]
    }

    # Prior mean for all nu j
    mu_nu <- apply(ADJ == 1, 1, function(x){mean(nu[x])})

    # Loop over spatial units
    for (j in 1:nJ) {
      # Posterior variance
      post_var <- solve(m[j] / xi2 + nis[j])

      # Posterior mean likelihood and prior contributions, variance of sampling distribution fixed to 1
      lik <- sum(mu[clusterIDSubf == j])
      prior <- m[j] / xi2 * mu_nu[j]
      post_mean <- post_var * (lik + prior)

      values[j] <- rnorm(1, mean = post_mean, sd = sqrt(post_var))
    }

    # Center on the fly
    values <- values - mean(values)

  }

  return(values)

}

