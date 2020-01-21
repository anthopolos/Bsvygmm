#' @title Update stratum level intercepts in the multinomial probit model.
#'
#' @description This uses a normal prior distribution with mean 0 to update the stratum level intercepts in the multinomial probit model. The number of latent classes is assumed to be greater than 2.
#' @param Z \code{K-1} column matrix of latent variables Z.
#' @param delta \code{K-1} column matrix of latent class-specific regression coefficients.
#' @param alpha \code{K-1} column matrix of regression coefficients for the B-spline.
#' @param u A \code{K-1} column matrix of area segment specific effects for modeling correlations among subjects in the same area segment. Equals \code{NULL} if not desired.
#' @param nu A \code{K-1} column matrix of spatial random effects. Equals \code{NULL} if not desired.
#' @param gamma2 A \code{K-1} vector of variances for the stratum level intercepts.
#' @param W An \code{s} column design matrix of covariates, including a column of one's for intercept.
#' @param B An \code{R} column design matrix of basis functions for the B-spline.
#' @param clusterIDSub Area segment or cluster identifier for each subject.
#' @param stratumIDSub Stratum identifier for each subject.
#' @param spline Logical. If \code{TRUE}, a B-spline is included for the basis functions in argument \code{B} for the latent class membership model.
#' @export
#' @return A \code{K - 1} matrix of stratum specific intercepts \code{delta_{sk}}.
update_deltaStratum_MNP <- function(Z, delta, alpha, u, nu, gamma2, W, B, clusterIDSub, stratumIDSub, spline) {

  K <- ncol(Z) + 1

  M <- length(unique(stratumIDSub))

  # Number of clusters in each stratum
  nj <- as.vector(table(stratumIDSub))

  if (!is.null(nu)) {
    # Redefine clusterIDSub as a factor that can miss certain levels
    nJ <- nrow(nu)
    clusterIDSubf <- factor(clusterIDSub, levels = seq(1, nJ, by = 1))
  }

  ### Store values
  values <- matrix(NA, nrow = M, ncol = (K - 1))

  for (k in 1:(K - 1)) {

    ### Posterior variance
    post_var <- solve(diag(nj + as.vector(solve(gamma2[k])), M))

    # 1. Unstructured heterogeneity at the cluster level only
    if (is.null(nu) & !is.null(u)) {
      if (spline) {
        llikSub <- Z[ , k] - W %*% delta[ , k] - B %*%  alpha[ , k] - u[factor(clusterIDSub), k]
      } else {
        llikSub <- Z[ , k] - W %*% delta[ , k] - u[factor(clusterIDSub), k]
      }
    }

    # 2. Spatially structured heterogeneity at the cluster level only
    if (!is.null(nu) & is.null(u)) {
      if (spline) {
        llikSub <- Z[ , k] - W %*% delta[ , k] - B %*%  alpha[ , k] - nu[clusterIDSubf, k]
      } else {
        llikSub <- Z[ , k] - W %*% delta[ , k] - nu[clusterIDSubf, k]
      }
    }

    # 3. Spatially structured and unstructured heterogeneity at the cluster level only
    if (!is.null(nu) & !is.null(u)) {
      if (spline) {
        llikSub <- Z[ , k] - W %*% delta[ , k] - B %*% alpha[ , k] - u[factor(clusterIDSub), k] - nu[clusterIDSubf, k]
      } else {
        llikSub <- Z[ , k] - W %*% delta[ , k] - u[factor(clusterIDSub), k] - nu[clusterIDSubf, k]
      }
    }

    # Posterior mean calculation
    #Sum lik contribution to stratum level
    lik <- as.vector(tapply(llikSub, stratumIDSub, sum))

    ### Posterior mean
    post_mean <- post_var %*% as.matrix(lik)

    values[ , k] <- rnorm(M, mean = post_mean, sd = sqrt(diag(post_var)))

  } # End k loop

  return(values)

}
