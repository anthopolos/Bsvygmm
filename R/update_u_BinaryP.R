#' @title Update area segment (aka cluster) level intercepts in the probit model.
#'
#' @description This uses a normal prior distribution with mean 0 to update the area segment (aka cluster) level intercepts in the probit model. These intercepts model correlation among subjects in the same area segment.
#' @param Z A \code{K-1} column matrix of latent variable \code{Z}.
#' @param delta A \code{K-1} column matrix of regression coefficients in the multinomial probit model with \code{k=1} fixed to zero.
#' @param alpha A \code{K-1} column matrix of regression coefficients for the B-spline.
#' @param deltaStratum A \code{K-1} column matrix of stratum specific effects.
#' @param nu A \code{K-1} column matrix of spatial random effects. Equals \code{NULL} if not desired.
#' @param W An \code{s} column design matrix of covariates, including a column of one's for intercept.
#' @param B An \code{R} column design matrix of basis functions for the B-spline.
#' @param clusterIDSub Area segment or cluster identifier for each subject.
#' @param stratumIDSub Stratum identifier for each subject.
#' @param spline Logical. If \code{TRUE}, a B-spline is included for the basis functions in argument \code{B} for the latent class membership model.
#' @export
#' @return A \code{K - 1} matrix of area segment specific intercepts \code{u}.
update_u_BinaryP <- function(Z, delta, alpha, deltaStratum, nu, tau2, W, B, clusterIDSub, stratumIDSub, spline) {

  J <- length(unique(clusterIDSub))

  # Redefine clusterIDSub as a factor that can miss certain levels
  if (!is.null(nu)) {
    nJ <- length(nu)
    clusterIDSubf <- factor(clusterIDSub, levels = seq(1, nJ, by = 1))
  }

  ### Compute the posterior variance
  #Variance in likelihood model fixed to 1
  # Number of subjects in cluster j
  nj <- as.vector(table(clusterIDSub))
  post_var <- solve(diag(nj + as.vector(solve(tau2)), J))

  # Modeling choices regarding spatial random effect
  if (!is.null(nu)) {
    if (spline) {
      llikSub <- Z - W %*% delta - B %*% alpha - deltaStratum[factor(stratumIDSub)] - nu[clusterIDSubf]
    } else {
      llikSub <- Z - W %*% delta - deltaStratum[factor(stratumIDSub)] - nu[clusterIDSubf]
    }
  }

  if (is.null(nu)) {
    if (spline) {
      llikSub <- Z - W %*% delta - B %*%  alpha - deltaStratum[factor(stratumIDSub)]
    } else {
      llikSub <- Z - W %*% delta - deltaStratum[factor(stratumIDSub)]
    }
  }

  # Posterior mean calculation
  #Sum lik contribution to cluster level
  lik <- as.vector(tapply(llikSub, clusterIDSub, sum))

  #Prior contribution cancels due to mean 0
  post_mean <- post_var %*% as.matrix(lik)

  # Draw random effects
  values <- rnorm(J, mean = post_mean, sd = sqrt(diag(post_var)))

  return(values)

}
