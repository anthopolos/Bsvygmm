#' @title Update area segment (aka cluster) level intercepts in the multinomial probit model.
#'
#' @description This uses a normal prior distribution with mean 0 to update the area segment (aka cluster) level intercepts in the multinomial probit model. These intercepts model correlation among subjects in the same area segment.
#' @param Z A \code{K-1} column matrix of latent variable \code{Z}.
#' @param delta A \code{K-1} column matrix of regression coefficients in the multinomial probit model with \code{k=1} fixed to zero.
#' @param alpha A \code{K-1} column matrix of regression coefficients for the B-spline.
#' @param deltaStratum A \code{K-1} column matrix of stratum specific effects. Equals \code{NULL} if not desired.
#' @param nu A \code{K-1} column matrix of spatial random effects. Equals \code{NULL} if not desired.
#' @param tau2 A \code{K-1} vector of variance terms for the independent random effects \code{u}. Equals \code{NULL} if not desired.
#' @param W An \code{s} column design matrix of covariates, including a column of one's for intercept.
#' @param B An \code{R} column design matrix of basis functions for the B-spline.
#' @param clusterIDSub Area segment or cluster identifier for each subject.
#' @param stratumIDSub Stratum identifier for each subject.
#' @param spline Logical. If \code{TRUE}, a B-spline is included for the basis functions in argument \code{B} for the latent class membership model.
#' @export
#' @return A \code{K - 1} matrix of area segment specific intercepts \code{u}.
update_u_MNP <- function(Z, delta, alpha, deltaStratum, nu, tau2, W, B, clusterIDSub, stratumIDSub, spline) {

  K <- ncol(delta) + 1
  J <- length(unique(clusterIDSub))

  if (!is.null(nu)) {
    # Redefine clusterIDSub as a factor that can miss certain levels
    nJ <- nrow(nu)
    clusterIDSubf <- factor(clusterIDSub, levels = seq(1, nJ, by = 1))
  }

  nj <- as.vector(table(clusterIDSub))

  ### Store values
  values <- matrix(NA, nrow = J, ncol = (K - 1))

  for (k in 1:(K - 1)) {

    ### Posterior variance
    post_var <- solve(diag(nj + as.vector(solve(tau2[k])), J))

    # With spatial random effects
    if (!is.null(nu)) {
      if (spline) {
        llikSub <- Z[ , k] - W %*% delta[ , k] - B %*% alpha[ , k] - nu[clusterIDSubf, k]
      } else {
        llikSub <- Z[ , k] - W %*% delta[ , k] - nu[clusterIDSubf, k]
      }
    }

    # Without spatial random effects
    if (is.null(nu)) {
      if (spline) {
        llikSub <- Z[ , k] - W %*% delta[ , k] - B %*%  alpha[ , k]
      } else {
        llikSub <- Z[ , k] - W %*% delta[ , k]
      }
    }

    # Account for stratum level random effects
    if (!is.null(deltaStratum)) {
      llikSub <- llikSub - deltaStratum[factor(stratumIDSub), k]
    }

    # Posterior mean calculation
    #Sum lik contribution to cluster level
    lik <- as.vector(tapply(llikSub, clusterIDSub, sum))

    # Prior contribution cancels because of zero mean
    ### Posterior mean
    post_mean <- post_var %*% as.matrix(lik)

    values[ , k] <- rnorm(J, mean = post_mean, sd = sqrt(diag(post_var)))

  } # End k loop

  return(values)

}

