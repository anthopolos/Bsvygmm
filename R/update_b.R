#' @title Update the subject level effects in the longitudinal model.
#'
#' @description This uses a normal prior distribution with mean 0 to update the subject level effects in the longitudinal outcomes model.
#'
#' @param C Latent class assignments for each subject.
#' @param rho  A \code{K} column matrix of area segment level intercepts.
#' @param zeta  A \code{K} column matrix of stratum level intercepts.
#' @param beta A \code{K} column matrix of regression coefficients associated with \code{Vf}.
#' @param Phi A variance-covariance matrix for the subject level effects.
#' @param sigma2 A \code{K} length vector of observation level variances.
#' @param Y Longitudinal measurements.
#' @param Vf A design matrix for fixed effects in the longitudinal outcome model.
#' @param Vr A design matrix for the subject-level random effects in the longitudinal outcome model. \code{Vr} is a subset or equal to \code{Vf}.
#' @param subjectID Subject identifier for each observation.
#' @param clusterIDObs Area segment or cluster identifier for each observation.
#' @param stratumIDObs Stratum identifier for each observation.
#' @export
#' @return A \code{q} column matrix of subject level intercepts.
update_b <- function(C, rho, zeta, beta, Phi, sigma2, Y, Vf, Vr, subjectID, clusterIDObs, stratumIDObs) {

  K <- length(table(C))
  n <- length(unique(subjectID))
  p <- ncol(as.matrix(Vf))
  q <- ncol(as.matrix(Vr))

  # Expand C
  C_expand <- C[factor(subjectID)]

  values <- matrix(NA, nrow = n, ncol = q)

  for (k in 1:K) {

    ### Observation level data
    ind_obs <- which(C_expand == k)
    subjectIDk <- subjectID[ind_obs]
    clusterIDObsk <- clusterIDObs[ind_obs]
    stratumIDObsk <- stratumIDObs[ind_obs]
    Yk <- Y[ind_obs]
    Vrk <- Vr[ind_obs, ]
    Vfk <- Vf[ind_obs, ]

    ### Subject level data
    ind_sub <- which(C == k)
    nk <- length(ind_sub)

    ### Cluster level information
    ind_cluster <- match(sort(unique(clusterIDObsk)), sort(unique(clusterIDObs)))
    ### Stratum level information
    ind_stratum <- match(sort(unique(stratumIDObsk)), sort(unique(stratumIDObs)))

    ### Unknown parameters
    sigma2k <- sigma2[k]
    betak <- beta[ , k]
    rhok  <- rho[ind_cluster, k]
    zetak  <- zeta[ind_stratum, k]

    rhokObs <- rhok[factor(clusterIDObsk)]
    zetakObs <- zetak[factor(stratumIDObsk)]

    Phik <- Phi[ , , k]

    bk <- sapply(unique(subjectIDk), function(x) {

      # Data for subject x in latent class k
      Yi <- Yk[which(subjectIDk == x)]
      Vri <- matrix(Vrk[which(subjectIDk == x), ], nrow = length(Yi), ncol = q)
      Vfi <- matrix(Vfk[which(subjectIDk == x), ], nrow = length(Yi), ncol = p)
      rhoi <- rhokObs[which(subjectIDk == x)]
      zetai <- zetakObs[which(subjectIDk == x)]

      # Posterior variance
      post_var <- solve(t(Vri) %*% Vri / sigma2k + solve(Phik))

      # Posterior mean
      post_mean <- post_var %*% t(t(Yi) %*% Vri - t(betak) %*% t(Vfi) %*% Vri - t(rhoi) %*% Vri - t(zetai) %*% Vri) / sigma2k

      bi <- mnormt::rmnorm(1, mean = post_mean, varcov = post_var)

      return(bi)

    } )

    values[ind_sub, ] <- t(bk)

  }

  return(values)

}
