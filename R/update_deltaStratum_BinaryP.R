#' @title Update stratum level intercepts in the probit model.
#'
#' @description This uses a normal prior distribution with mean 0 to update the stratum level intercepts in the probit model.
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
update_deltaStratum_BinaryP <- function(Z, delta, alpha, u, nu, gamma2, W, B, clusterIDSub, stratumIDSub, spline) {

  M <- length(unique(stratumIDSub))

  if (!is.null(nu)) {
    # Redefine clusterIDSub as a factor that can miss certain levels
    nJ <- length(nu)
    clusterIDSubf <- factor(clusterIDSub, levels = seq(1, nJ, by = 1))
  }

  ### Posterior variance calculation
  # Number of clusters in each stratum
  nj <- as.vector(table(stratumIDSub))
  post_var <- solve(diag(nj + as.vector(solve(gamma2)), M))
  # 1. Unstructured heterogeneity at the cluster level only
  if (is.null(nu) & !is.null(u)) {
    if (spline) {
      llikSub <- Z - W %*% delta - B %*%  alpha - u[factor(clusterIDSub)]
    } else {
      llikSub <- Z - W %*% delta - u[factor(clusterIDSub)]
    }
  }

  # 2. Spatially structured heterogeneity at the cluster level only
  if (!is.null(nu) & is.null(u)) {
    if (spline) {
      llikSub <- Z - W %*% delta - B %*%  alpha - nu[clusterIDSubf]
    } else {
      llikSub <- Z - W %*% delta - nu[clusterIDSubf]
    }
  }

  # 3. Spatially structured and unstructured heterogeneity at the cluster level only
  if (!is.null(nu) & !is.null(u)) {
    if (spline) {
      llikSub <- Z - W %*% delta - B %*% alpha - u[factor(clusterIDSub)] - nu[clusterIDSubf]
    } else {
      llikSub <- Z - W %*% delta - u[factor(clusterIDSub)] - nu[clusterIDSubf]
    }
  }

  ### Posterior mean
  #Sum lik contribution to stratum level
  lik <- as.vector(tapply(llikSub, stratumIDSub, sum))
  post_mean <- post_var %*% as.matrix(lik)

  # Draw random effects
  values <- rnorm(M, mean = post_mean, sd = sqrt(diag(post_var)))

  return(values)

}
