#' @title Update regression coefficients \code{alpha} for the B-spline in the probit model.
#'
#' @description This uses a normal prior distribution to update the regression coefficients for the B-spline in the probit model. The mean of the prior distribution is fixed to 0. The regression coefficients for \code{k=1} are fixed to zero.
#' @param Z A \code{K-1} column matrix of latent variable \code{Z}.
#' @param delta A \code{K-1} column matrix of regression coefficients in the probit model with \code{k=1} fixed to zero.
#' @param deltaStratum A \code{K-1} column matrix of stratum specific effects. Equals \code{NULL} if not desired.
#' @param u A \code{K-1} column matrix of area segment specific effects for modeling correlations among subjects in the same area segment. Equals \code{NULL} if not desired.
#' @param nu A \code{K-1} column matrix of spatial random effects. Equals \code{NULL} if not desired.
#' @param W An \code{s} column design matrix of covariates, including a column of one's for intercept.
#' @param B An \code{R} column design matrix of basis functions for the B-spline.
#' @param clusterIDSub Area segment or cluster identifier for each subject.
#' @param stratumIDSub Stratum identifier for each subject.
#' @param prior.Sigma A covariance matrix supplied through \code{priors}.
#' @param clusterIDSub Area segment or cluster identifier for each subject.
#' @param stratumIDSub Stratum identifier for each subject.
#' @export
#' @return A \code{K-1} column matrix of regression coefficients for the B-spline with \code{k=1} fixed to zero.
update_alpha_BSpline_BinaryP <- function(Z, delta, deltaStratum, u, nu, W, B, prior.Sigma, clusterIDSub, stratumIDSub) {

  K <- ncol(as.matrix(delta)) + 1
  nB <- ncol(B)

  if (!is.null(nu)) {
    nJ <- length(nu)
    #Only expand nu for clusterIDSubs that appear in data
    clusterIDSubf <- factor(clusterIDSub, levels = seq(1, nJ, by = 1))
  }

  post_var <- solve(t(B) %*% B + solve(prior.Sigma))

  # 1. Spatially structured heterogeneity at the cluster level only
  if (!is.null(nu) & is.null(u)) {
    mu <-  Z - W %*% delta - nu[clusterIDSubf]
  }

  # 2. Unstructured heterogeneity at the cluster level only
  if (is.null(nu) & !is.null(u)) {
    mu <-  Z - W %*% delta - u[factor(clusterIDSub)]
  }

  # 3. Spatially structured and unstructured heterogeneity at the cluster level only
  if (!is.null(nu) & !is.null(u)) {
    mu <-  Z - W %*% delta - u[factor(clusterIDSub)] - nu[clusterIDSubf]
  }

  # 4. No spatially structured and unstructured heterogeneity at the cluster level only
  if (is.null(nu) & is.null(u)) {
    mu <-  Z - W %*% delta
  }


  # Account for stratum level intercept
  if (!is.null(deltaStratum)) {
    mu <- mu - deltaStratum[factor(stratumIDSub)]
  }


  post_mean <- post_var %*% (t(B) %*% mu)

  values <- matrix(mnormt::rmnorm(1, mean = post_mean, varcov = post_var), nrow = nB, ncol = K - 1)

  return(values)

}

