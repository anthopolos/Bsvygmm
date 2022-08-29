#' @title Update latent variable \code{Z} in the probit model.
#'
#' @description The latent variable in the probit model is updated using a truncated normal distribution. The updating algorithm allows for either or both area segment intercepts for correlations among subjects in the same area segment and spatial correlations among neighboring area segments.
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
update_Z_BinaryP <- function(Z, C, delta, alpha, deltaStratum, u, nu, W, B, clusterIDSub, stratumIDSub, spline) {

  n <- length(C)

  values <- rep(NA, n)

  if (!is.null(nu)) {
    nJ <- length(nu)
    #Only expand nu for clusterIDSubs that appear in data
    clusterIDSubf <- factor(clusterIDSub, levels = seq(1, nJ, by = 1))
  }

  # 1. Spatially structured heterogeneity at the cluster level only
  if (!is.null(nu) & is.null(u)) {
    if (spline) {
      mu <- W %*% delta + B %*% alpha + nu[clusterIDSubf]
    } else {
      mu <- W %*% delta + nu[clusterIDSubf]
    }
  }


  # 2. Unstructured heterogeneity at the cluster level only
  if (is.null(nu) & !is.null(u)) {
    if (spline) {
      mu <- W %*% delta + B %*% alpha + u[factor(clusterIDSub)]
    } else {
      mu <- W %*% delta + u[factor(clusterIDSub)]
    }
  }

  # 3. Spatially structured and unstructured heterogeneity at the cluster level

  if (!is.null(nu) & !is.null(u)) {
    if (spline) {
      mu <- W %*% delta + B %*% alpha + u[factor(clusterIDSub)] + nu[clusterIDSubf]
    } else {
      mu <- W %*% delta + u[factor(clusterIDSub)] +  nu[clusterIDSubf]
    }
  }

  # 4. No spatially structured and unstructured heterogeneity at the cluster level
  if (is.null(nu) & is.null(u)) {
    if (spline) {
      mu <- W %*% delta + B %*% alpha
    } else {
      mu <- W %*% delta
    }
  }


  # Add stratum level random effect if desired
  if (!is.null(deltaStratum)) {
    mu <- mu + deltaStratum[factor(stratumIDSub)]
  }


  values[C == 1] <- truncnorm::rtruncnorm(n = n, mean = mu, sd = 1, a = -Inf, b = 0)[C == 1]
  values[C == 2] <- truncnorm::rtruncnorm(n = n, mean = mu, sd = 1, a = 0, b = Inf)[C == 2]

  return(values)

}
