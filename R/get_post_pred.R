#' @title Get posterior predictive draws at each iteration.
#'
#' @description The reference for the approach is from Gelman BDA page 530-31. To conduct model checking, we take posterior predictive draws of \code{Y} conditional on \code{C}.
#'
#' @param K Number of latent classes.
#' @param C Categorical variable for the latent class of each subject.
#' @param b A \code{q} column matrix of subject-level effects.
#' @param rho A \code{K} column matrix of area segment level intercepts. Equals \code{NULL} if not desired.
#' @param zeta A \code{K} column matrix of stratum level intercepts. Equals \code{NULL} if not desired.
#' @param beta Current iteration of fixed effects.
#' @param sigma2 Current iteration of observation level variance.
#' @param Vf Design matrix for fixed effects.
#' @param Vr Design matrix for random effects.
#' @param subjectID Subject identifier for each observation.
#' @param clusterIDObs Area segment or cluster identifier for each observation.
#' @param stratumIDObs Stratum identifier for each observation.
#' @export
#' @return Draws from the posterior predictive distribution at the current iteration.
get_postPred <- function(K, C, b, rho, zeta, beta, sigma2, Vr, Vf, subjectID, clusterIDObs, stratumIDObs) {

    C_expand <- C[factor(subjectID)]

    # At each iteration store Ydraws from each latent class
    Ydraw <- rep(NA, length(subjectID))

    for (k in 1:K) {

      ### Subset information based on the C
      #Observation level information
      ind_obs <- which(C_expand == k)
      subjectIDk <- subjectID[ind_obs]
      clusterIDObsk <- clusterIDObs[ind_obs]
      stratumIDObsk <- stratumIDObs[ind_obs]
      Vfk <- Vf[ind_obs, ]
      Vrk <- Vr[ind_obs, ]

      #Subject level information
      ind_sub <- which(C == k)

      #Cluster level information
      ind_cluster <- match(sort(unique(clusterIDObsk)), sort(unique(clusterIDObs)))

      #Stratum level information
      ind_stratum <- match(sort(unique(stratumIDObsk)), sort(unique(stratumIDObs)))

      # Unknown parameters that do not need to be redrawn
      bk <- b[ind_sub, ]
      sigma2k <- sigma2[k]
      betak <- beta[ , k]

      sp_Vr_sub <- lapply(split(as.data.frame(Vrk), subjectIDk, drop = TRUE), as.matrix)
      convert_Vr_sub <- as.matrix(Matrix::bdiag(sp_Vr_sub))

      # Calculate linear predictor and draw Y for class k
      lp <- as.vector(Vfk %*% betak + as.vector(convert_Vr_sub %*% c(t(bk))))

      # Include stratum level intercepts or not
      #Expand rho and zeta to observation levels for class k
      if (!is.null(zeta)) {
        zetak  <- zeta[ind_stratum, k]
        zetakObs <- zetak[factor(stratumIDObsk)]
        lp <- lp + zetakObs
      }

      # Include cluster level intercepts or not
      if (!is.null(rho)) {
        rhok  <- rho[ind_cluster, k]
        rhokObs <- rhok[factor(clusterIDObsk)]
        lp <- lp + rhokObs
      }


      ### Draw Ypred for this iteration l
      Ydraw[ind_obs] <- lp + rnorm(length(lp), 0, sqrt(sigma2k))

    }

  return(Ydraw)

}
