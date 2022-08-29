#' @title Update latent class assignments.
#'
#' @description To update the latent class assignments, a partially marginalized sampler is used in which \code{zeta}, \code{rho}, and \code{b} are integrated out.
#'
#' @param priorPik Current iteration of prior probability of class membership.
#' @param beta A \code{K} column matrix of regression coefficients associated with \code{Vf}.
#' @param Phi Variance-covariance of the subject specific effects.
#' @param Omega Variance at the area segment or cluster level. Equals \code{NULL} if not desired.
#' @param Psi Variance at the stratum level. Equals \code{NULL} if not desired.
#' @param sigma2 Observation level variance.
#' @param Y Longitudinal measurements.
#' @param Vf A design matrix for fixed effects in the longitudinal outcome model.
#' @param Vr A design matrix for the subject-level random effects in the longitudinal outcome model. \code{Vr} is a subset or equal to \code{Vf}.
#' @param subjectID Subject identifier for each observation.
#' @export
#' @return A list containing a vector of latent class assignments and a \code{K} column matrix of posterior probabilities of belonging to each class.
update_C_Marginal <- function(priorPik, beta, Phi, Omega, Psi, sigma2, Y, Vf, Vr, subjectID) {

  ### Function to handle numerical underflow in posterior probability that subject i belongs to class k, log sum of exponentials trick as described in Stan manual
  LnToRaw <- function(a){

    f <- function(a){
      A <- max(a) + log(sum(exp(a - max(a))))
      exp(a - A)
    }

    t(apply(a, 1, f))

  }

  ### Number of classes
  K <- ncol(as.matrix(priorPik))

  ### Number of subjects
  n <- length(unique(subjectID))

  ### Likelihood contribution for subject i to each class k
  q <- ncol(as.matrix(Vr))
  mu <- Vf %*% beta

  llik_y_pik <- matrix(NA, nrow = n, ncol = K)

  for (k in 1:K) {

    muk <- mu[ , k]
    sigma2k <- sigma2[k]

    Phik <- Phi[ , , k]

    if (!is.null(Omega)) {
      Omegak <- Omega[k]
    } else {
      Omegak <- NULL
    }

    if (!is.null(Psi)) {
      Psik <- Psi[k]
    } else {
      Psik <- NULL
    }

    # Compute subject level likelihood contribution
    # Log likelihood contribution for the subject level equation

    # Variance for time point 1: \sigma^2_{b_{1}} + \sigma^2_{\rho_{1}} + \sigma^2_{\epsilon}
    # Covariance between time point 1 and 2: Cov(b_{i1}, b_{i2}) + \sigma^2_{\rho_{1}}
    llik_y_pikTemp <- sapply(unique(subjectID), function(x) {

      Yi <- Y[which(subjectID == x)]
      muki <- muk[which(subjectID == x)]
      Vri <- matrix(as.matrix(Vr)[which(subjectID == x), ], nrow = length(Yi), ncol = q)

      # Compute the variance covariance of Y_i depending on whether Omega and Psi are included in the longitudinal model
      vcovYi <- (Vri %*% Phik %*% t(Vri)) + diag(sigma2k, nrow = length(Yi), ncol = length(Yi))

      if (!is.null(Omegak)) {
        vcovYi <- vcovYi + matrix(rep(Omegak, length(Yi) * length(Yi)), nrow = length(Yi), ncol = length(Yi))
      }

      if (!is.null(Psik)) {
        vcovYi <- vcovYi + matrix(rep(Psik, length(Yi) * length(Yi)), nrow = length(Yi), ncol = length(Yi))
      }

      # Density evaluated at Yi
      mvtnorm::dmvnorm(Yi, mean = muki, sigma = vcovYi, log = TRUE)
    })

    llik_y_pik[ , k] <- llik_y_pikTemp

  }

  # Numerator
  # priorPik has class K from undifferenced parameterization in position 1 so re-ordering here is unnecessary
  pik_num <- log(priorPik) + llik_y_pik
  # Posterior probability of subject i belonging in each class
  pik <- LnToRaw(pik_num)

  pik <- as.matrix(unname(pik))

  C <- Hmisc::rMultinom(pik, 1)

  list(C = C, pik = pik)

}
