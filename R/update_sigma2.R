#' @title Update the observation level variance in the longitudinal outcomes model.
#'
#' @description An inverge gamma or uniform prior is used to update the variance at the observation level.
#' @param C Latent class assignments for each subject.
#' @param zeta  A \code{K} column matrix of intercepts at the stratum level.
#' @param rho  A \code{K} column matrix of intercepts at the area segment level.
#' @param b  A \code{q} column matrix of subject specific effects.
#' @param beta A \code{K} column matrix of regression coefficients associated with \code{Vf}.
#' @param Y Longitudinal measurements.
#' @param Vf A design matrix for fixed effects in the longitudinal outcome model.
#' @param Vr A design matrix for the subject-level random effects in the longitudinal outcome model. \code{Vr} is a subset or equal to \code{Vf}.
#' @param prior.shape A scalar for the prior shape, if an inverse gamma prior is used.
#' @param prior.scale A scalar for the prior scale, if an inverse gamma prior is used.
#' @param subjectID Subject identifier for each observation.
#' @param clusterIDObs Area segment or cluster identifier for each observation.
#' @param stratumIDObs Stratum identifier for each observation.
#' @param hierVar A string equal to \code{Unif} or \code{IG} to indicate a uniform or inverse gamma prior respectively.
#' @export
#' @return A \code{K} length vector of latent class specific variances.
update_sigma2 <- function(C, zeta, rho, b, beta, Y, Vr, Vf, prior.shape, prior.scale, subjectID, clusterIDObs, stratumIDObs, hierVar){

  K <- length(table(C))

  C_expand <- C[factor(subjectID)]

  values <- rep(NA, K)

  for (k in 1:K) {

    ### Observation level data
    ind_obs <- which(C_expand == k)
    Yk <- Y[ind_obs]
    Vrk <- Vr[ind_obs, ]
    Vfk <- Vf[ind_obs, ]
    subjectIDk <- subjectID[ind_obs]
    clusterIDObsk <- clusterIDObs[ind_obs]
    stratumIDObsk <- stratumIDObs[ind_obs]

    ### Cluster level information
    ind_cluster <- match(sort(unique(clusterIDObsk)), sort(unique(clusterIDObs)))

    ### Stratum level information
    ind_stratum <- match(sort(unique(stratumIDObsk)), sort(unique(stratumIDObs)))

    ### Unknown paramters
    betak <- beta[ , k]
    rhok  <- rho[ind_cluster, k]
    zetak  <- zeta[ind_stratum, k]

    rhokObs <- rhok[factor(clusterIDObsk)]
    zetakObs <- zetak[factor(stratumIDObsk)]

    # Construct N x n*q matrix of X to facilitate matrix algebra
    #[1 1 0 0 0 / 1 2 0 0 0 / 1 3 0 0 0], for q=2 columns 1 and 2 are subject 1's intercept and random slope for time
    sp_Vr_sub <- lapply(split(as.data.frame(Vrk), subjectIDk, drop = TRUE), as.matrix)
    convert_Vr_sub <- as.matrix(Matrix::bdiag(sp_Vr_sub)) # N x nq

    ### Subject level information
    ind_sub <- which(C == k)

    bk <- b[ind_sub, ] # Here we choose class k
    muy <- Vfk %*% betak + convert_Vr_sub %*% c(t(bk)) + rhokObs + zetakObs

    llik <- crossprod(Yk - muy, Yk - muy)
    df <- length(Yk)

    if (hierVar == "IG") {
      post_scale <- llik / 2 + prior.scale
      post_shape <- df / 2 + prior.shape
    } else if (hierVar == "Unif") {
      post_scale <-  llik / 2
      # Posterior shape
      post_shape <- (df - 1) / 2  # A scalar
    }

    values[k] <- 1 / rgamma(1, post_shape, post_scale)

  }

  return(values)

}
