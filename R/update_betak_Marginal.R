#' @title Update the regression coefficients in the longitudinal outcomes model.
#'
#' @description A normal prior is used to update the regression coefficients in the longitudinal outcomes model. A partially marginalized sampler is used in which \code{zeta}, \code{rho}, and \code{b} are integrated out.
#' @param C Latent class assignments for each subject.
#' @param sigma2 Observation level variance.
#' @param Phi Variance-covariance of the subject specific effects.
#' @param Omega Variance at the area segment or cluster level.  Equals \code{NULL} if not desired.
#' @param Psi Variance at the stratum level.  Equals \code{NULL} if not desired.
#' @param Y Longitudinal measurements.
#' @param Vf A design matrix for fixed effects in the longitudinal outcome model.
#' @param Vr A design matrix for the subject-level random effects in the longitudinal outcome model. \code{Vr} is a subset or equal to \code{Vf}.
#' @param prior.mu Prior mean.
#' @param prior.Sigma Prior variance.
#' @param subjectID Subject identifier for each observation.
#' @export
#' @return A \code{K} column matrix of latent class specific regression coefficients.
update_betak_Marginal <- function(C, sigma2, Phi, Omega, Psi,  prior.mu, prior.Sigma, Y, Vf, Vr, subjectID) {

  K <- length(table(C))
  p <- ncol(as.matrix(Vf))
  q <- ncol(as.matrix(Vr))

  C_expand <- C[factor(subjectID)]

  values <- matrix(NA, nrow = p, ncol = K)

  for (k in 1:K) {

    ### Observation level
    ind_obs <- which(C_expand == k)
    subjectIDk <- subjectID[ind_obs]
    Yk <- Y[ind_obs]
    Vfk <- Vf[ind_obs, ]
    Vrk <- Vr[ind_obs, ]

    Phik <- Phi[ , , k]

    # Cluster level information included or not
    if (!is.null(Omega)) {
      Omegak <- Omega[k]
    } else {
      Omegak <- NULL
    }

    # Stratum level information included or not
    if (!is.null(Psi)) {
      Psik <- Psi[k]
    } else {
      Psik <- NULL
    }

    sigma2k <- sigma2[k]

    ### Compute posterior variance
    llikContr <- sapply(unique(subjectIDk), function(x) {

      # Observations to subject x
      Yi <- Yk[which(subjectIDk == x)]
      Vri <- matrix(as.matrix(Vrk)[which(subjectIDk == x), ], nrow = length(Yi), ncol = q)
      Vfi <- matrix(Vfk[which(subjectIDk == x), ], nrow = length(Yi), ncol = p)

      # Marginal variance for each subject based on class k depends on whether cluster and stratum level random effects were included
      Rik <- (Vri %*% Phik %*% t(Vri)) + diag(sigma2k, nrow = length(Yi), ncol = length(Yi))

      # Cluster level
      if (!is.null(Omegak)) {
        Rik <- Rik + matrix(rep(Omegak, length(Yi) * length(Yi)), nrow = length(Yi), ncol = length(Yi))
      }

      # Stratum level
      if (!is.null(Psik)) {
        Rik <- Rik + matrix(rep(Psik, length(Yi) * length(Yi)), nrow = length(Yi), ncol = length(Yi))
      }

      # Likelihood contributions in marginal model
      list(t(Vfi) %*% solve(Rik)  %*% Vfi, t(Yi) %*% solve(Rik) %*% Vfi)

    }, simplify = FALSE)

    ### Sum elements from element 1 of nested list to compute posterior variance
    post_var <- solve(Reduce("+", lapply(llikContr, function(x) x[[1]])) + solve(prior.Sigma))
    ### Posterior mean with element 2
    post_mean <- post_var %*% (t(Reduce("+", lapply(llikContr, function(x) x[[2]]))) + solve(prior.Sigma) %*% prior.mu)

    values[ , k] <- mnormt::rmnorm(1, mean = post_mean, varcov = post_var)

  }

  return(values)

}
