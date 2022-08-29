#' @title Get the effective sample size.
#'
#' @description The reference for the approach is in Jones, R. H. (2011). Bayesian information criterion for longitudinal and clustered data. Statistics in Medicine, 30(25), 3050-3056.
#'
#' @param sigma2 Current iteration of observation level variance.
#' @param Phi Current iteration of between subject variance-covariance.
#' @param Omega Current iteration of between area segment level variance. Equals \code{NULL} if not desired.
#' @param Psi Current iteration of between stratum level variance.  Equals \code{NULL} if not desired.
#' @param priorPik Current iteration of prior probability of class membership.
#' @param Y Longitudinal outcomes.
#' @param Vr Design matrix for random effects.
#' @param subjectID Subject ID for each observation.
#' @return Effective sample size at the current iteration.
get_ESS <- function(sigma2, Phi, Omega, Psi, priorPik, Y, Vr, subjectID) {

  q <- ncol(as.matrix(Vr))
  K <- ncol(priorPik)

  subjectIDSub <- unique(subjectID)

  ### What is the sample size for each subject based on marginalizing over latent class and random effects
  ESS_sub <- sapply(unique(subjectID), function(x) {

    Yi <- Y[which(subjectID == x)]
    Vri <- matrix(Vr[which(subjectID == x), ], nrow = length(Yi), ncol = q)

    # Prior probability vector for x^th subject
    priorPiki <- matrix(priorPik[which(subjectIDSub == x), ], nrow = K, ncol = 1)

    # Marginal variance computation
    Phii <- Reduce("+", sapply(1:K, function(x) {

      priorPiki[x] * Phi[ , , x]

    }, simplify = FALSE))

    sigma2i <- sum(sigma2 * priorPiki)

    margRi <- (Vri %*% Phii %*% t(Vri)) + diag(sigma2i, nrow = length(Yi), ncol = length(Yi))

    # Include cluster level variance
    if(!is.null(Omega)) {
      Omegai <- sum(Omega * priorPiki)
      margRi <- margRi + matrix(rep(Omegai, length(Yi) * length(Yi)), nrow = length(Yi), ncol = length(Yi))
    }

    # Include stratum level variance
    if(!is.null(Psi)) {
      Psii <- sum(Psi * priorPiki)
      margRi <- margRi + matrix(rep(Psii, length(Yi) * length(Yi)), nrow = length(Yi), ncol = length(Yi))
    }

    sum(solve(cov2cor(margRi)))

  })

  sum(ESS_sub)

}
