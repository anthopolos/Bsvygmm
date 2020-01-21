#' @title Get the observed data likelihood after integrating out the random effects.
#'
#' @description The reference for the approach is in Celeux, G., Forbes, F., Robert, C. P., & Titterington, D. M. (2006). Deviance Information Criteria for Missing Data Models. Bayesian Analysis, 1(4), 651-674.
#'
#' @param priorPik Current iteration of prior probability of class membership.
#' @param beta Current iteration of fixed effects.
#' @param sigma2 Current iteration of observation level variance.
#' @param Phi Current iteration of between subject variance-covariance.
#' @param Omega Current iteration of between area segment level variance.
#' @param Psi Current iteration of between stratum level variance.
#' @param Y Longitudinal outcomes.
#' @param Vf Design matrix for fixed effects.
#' @param Vr Design matrix for random effects.
#' @param subjectID Subject ID for each observation.
#' @return Observed data likelihood at the current iteration.
get_observed_llik <- function(priorPik, beta, sigma2, Phi, Omega, Psi, Y, Vf, Vr, subjectID) {

  n <- length(unique(subjectID))
  K <- dim(priorPik)[2]
  q <- ncol(Vr)

  # Marginal mean
  mu <- Vf %*% beta

  ### Storage for likelihood contribution of each subject to each class
  llik_y_pik <- matrix(NA, nrow = n, ncol = K)

  ### K loop for likelihood contribution by latent class
  for (k in 1:K) {

    muk <- mu[ , k]
    sigma2k <- sigma2[k]

    Phik <- Phi[ , , k]
    Omegak <- Omega[k]
    Psik <- Psi[k]

    llik_y_pikTemp <- sapply(unique(subjectID), function(x) {

      Yi <- Y[which(subjectID == x)]
      muki <- muk[which(subjectID == x)]
      Vri <- matrix(Vr[which(subjectID == x), ], nrow = length(Yi), ncol = q)
      Rki <- ((Vri %*% Phik %*% t(Vri)) + diag(sigma2k, nrow = length(Yi), ncol = length(Yi)) + matrix(rep(Omegak, length(Yi) * length(Yi)), nrow = length(Yi), ncol = length(Yi)) + matrix(rep(Psik, length(Yi) * length(Yi)), nrow = length(Yi), ncol = length(Yi)))
      mvtnorm::dmvnorm(Yi, mean = muki, sigma = Rki, log = FALSE)
    })

    llik_y_pik[ , k] <- llik_y_pikTemp

  }

  # Log likelihood calculation at iteration $l$:
  # \log f(y)^l = & \log \prod_{i = 1}^n \sum_{k = 1}^K \p_{ik}^l f(y_i \ | \ \Theta_k^l)
  # = & \sum_{i = 1}^n \log \sum_{k = 1}^K \p_{ik}^l f(y_i \ | \ \Theta_k^l)
  # For iteration $l$, calculate \log \sum_{k = 1}^K \p_{ik}^l f(y_i \ | \ \Theta_k^l)
  #Store_llikTemp will be 1 by n
  store_llikTemp <- apply(priorPik * llik_y_pik, 1, get_log_sum_exponential)

  ### Log likelihood estimate
  return(sum(store_llikTemp))

}
