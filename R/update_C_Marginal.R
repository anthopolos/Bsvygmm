#' @title Update latent class assignments.
#'
#' @description To update the latent class assignments, a partially marginalized sampler is used in which \code{zeta}, \code{rho}, and \code{b} are integrated out.
#'
#' @param K Number of latent classes.
#' @param Z Latent draws from the latent class membership model.
#' @param beta A \code{K} column matrix of regression coefficients associated with \code{Vf}.
#' @param Phi Variance-covariance of the subject specific effects.
#' @param Omega Variance at the area segment or cluster level.
#' @param Psi Varianec at the stratum level.
#' @param sigma2 Observation level variance.
#' @param Y Longitudinal measurements.
#' @param Vf A design matrix for fixed effects in the longitudinal outcome model.
#' @param Vr A design matrix for the subject-level random effects in the longitudinal outcome model. \code{Vr} is a subset or equal to \code{Vf}.
#' @param subjectID Subject identifier for each observation.
#' @export
#' @return A list containing a vector of latent class assignments and a \code{K} column matrix of posterior probabilities of belonging to each class.
update_C_Marginal <- function(K, Z, beta, Phi, Omega, Psi, sigma2, Y, Vf, Vr, subjectID) {

  ### Function to handle numerical underflow in posterior probability that subject i belongs to class k, log sum of exponentials trick as described in Stan manual
  LnToRaw <- function(a){

    f <- function(a){
      A <- max(a) + log(sum(exp(a - max(a))))
      exp(a - A)
    }

    t(apply(a, 1, f))

  }

  ### Number of subjects
  n <- length(unique(subjectID))

  ### Prior probability of class assignment
  priorPik <- matrix(NA, nrow = n, ncol = K)

  if (K == 2) {
    # K = 2:

    # First column is reference
    priorPik[ , K] <- pnorm(Z)
    priorPik[ , 1] <- 1 - priorPik[ , K]

  } else {

    # K > 2:

    # Construct original matrix R from Z based parameterization
    #We have been working with the differenced parameterization such that the coefficients represent a difference from baseline class K, where after a simple reparameterization, class K is switched to class 0 (the reference group moves from the Kth column to the first column)
    #If we assume the original parameterization is represented by R, then
    #Z_{i1} = R_{i1} - R_{iK}; Z_{i2} = R_{i2} - R_{iK}, and so forth
    #R_{iK} \sim N(0, 1) since \delta_K and \theta_K are fixed to 0
    #R is an n x K matrix, where the K^th column is that with \delta_k fixed to 0 (the K^th column is the reference column)
    R <- matrix(NA, nrow = n, ncol = K)
    R[ , K] <- rnorm(n, mean = 0, sd = 1)
    for (k in 1:(K - 1)) {
      R[ , k] <- Z[, k] + R[ , K]
    }

    # Covariance matrix
    D <- matrix(1, nrow = (K - 1), ncol = (K - 1))
    diag(D) <- 2

    # Values where distribution is evaluated
    vals <- matrix(0, nrow = n, ncol = (K - 1))
    for (k in 1:K) {
      mur_diffk <- R[ , -k] - R[ , k]
      priorPik[ , k] <- mnormt::pmnorm(vals, mean = mur_diffk, varcov = D)
    }

    priorPik <- as.matrix(unname(priorPik))

  }

  ### Likelihood contribution for subject i to each class k
  p <- ncol(as.matrix(Vf))
  mu <- Vf %*% beta

  llik_y_pik <- matrix(NA, nrow = n, ncol = K)

  for (k in 1:K) {

    muk <- mu[ , k]
    sigma2k <- sigma2[k]

    Phik <- Phi[ , , k]
    Omegak <- Omega[k]
    Psik <- Psi[k]


    # Compute subject level likelihood contribution
    # Log likelihood contribution for the subject level equation

    # Variance for time point 1: \sigma^2_{b_{1}} + \sigma^2_{\rho_{1}} + \sigma^2_{\epsilon}
    # Covariance between time point 1 and 2: Cov(b_{i1}, b_{i2}) + \sigma^2_{\rho_{1}}
    llik_y_pikTemp <- sapply(unique(subjectID), function(x) {

      Yi <- Y[which(subjectID == x)]
      muki <- muk[which(subjectID == x)]
      Vfi <- matrix(Vf[which(subjectID == x), ], nrow = length(Yi), ncol = p)
      mvtnorm::dmvnorm(Yi, mean = muki, sigma = ((Vfi %*% Phik %*% t(Vfi)) + diag(sigma2k, nrow = length(Yi), ncol = length(Yi)) + matrix(rep(Omegak, length(Yi) * length(Yi)), nrow = length(Yi), ncol = length(Yi)) + matrix(rep(Psik, length(Yi) * length(Yi)), nrow = length(Yi), ncol = length(Yi))), log = TRUE)
    })

    llik_y_pik[ , k] <- llik_y_pikTemp

  }

  # Reordering is necessary K > 2 but not for K = 2
  if (K > 2) {
    # Reorder latent classes to be on original parameterization
    #Because the McCulloch and Rossi parameterization effectively switches latent class K to class 0, and because the prior pik are computed on the original un-differenced scale, we need to make sure that priorPik and llik_y_pik are consistently ordered.
    #In llik_y_pik, the the reference group in first column is returned to the last column to be consistent with the original parameterization
    llik_y_pikReorder <- cbind(llik_y_pik[, 2:K], llik_y_pik[, 1])

    # Numerator
    pik_num <- log(priorPik) + llik_y_pikReorder

    # Posterior probability of subject i belonging in each class
    pikReorder <- LnToRaw(pik_num)

    # Draw C from multinomial, first returning pik to the formulation in McCulloch and Rossi that is the basis for latent normal bounds updating. This means that column K moves to column 1 to be the reference class. Reordering pik will funnel through to C.
    pik <- cbind(pikReorder[ , K], pikReorder[ , 1:(K - 1)])
  } else {

    # K = 2
    # Numerator
    pik_num <- log(priorPik) + llik_y_pik

    # Posterior probability of subject i belonging in each class
    pik <- LnToRaw(pik_num)

  }

  pik <- as.matrix(unname(pik))
  C <- Hmisc::rMultinom(pik, 1)

  list(C = C, pik = pik)

}
