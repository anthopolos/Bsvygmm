#' @title Get posterior predictive draws at each iteration.
#'
#' @description The reference for the approach is from algorithm 3.7 on page 91 of Fruhwirth-Schnatter, S. (2006) Finite Mixture and Markov Switching Models. Springer Science & Business Media, New York. The latent classes are redrawn using the probabilities of latent class membership obtained from the latent class membership model.
#'
#' @param K Number of latent classes.
#' @param beta Current iteration of fixed effects.
#' @param sigma2 Current iteration of observation level variance.
#' @param priorPik Current iteration of prior probability of class membership.
#' @param b A \code{q} column matrix of subject-level effects.
#' @param rho A \code{K} column matrix of area segment level intercepts.
#' @param zeta A \code{K} column matrix of stratum level intercepts.
#' @param Vf Design matrix for fixed effects.
#' @param Vr Design matrix for random effects.
#' @param subjectID Subject identifier for each observation.
#' @param clusterIDObs Area segment or cluster identifier for each observation.
#' @param stratumIDObs Stratum identifier for each observation.
#' @return Draws from the posterior predictive distribution at the current iteration.
get_postPred <- function(K, beta, sigma2, priorPik, b, rho, zeta, Vr, Vf, subjectID, clusterIDObs, stratumIDObs) {

    Cdraw <- Hmisc::rMultinom(priorPik, 1)

    Cdraw_expand <- Cdraw[factor(subjectID)]

    # At each iteration store Ydraws from each latent class
    Ydraw <- rep(NA, length(subjectID))

    for (k in 1:K) {

      ### Subset information based on the Cdraw
      #Observation level information
      ind_obs <- which(Cdraw_expand == k)
      subjectIDk <- subjectID[ind_obs]
      clusterIDObsk <- clusterIDObs[ind_obs]
      stratumIDObsk <- stratumIDObs[ind_obs]
      Vfk <- Vf[ind_obs, ]
      Vrk <- Vr[ind_obs, ]

      #Subject level information
      ind_sub <- which(Cdraw == k)

      #Cluster level information
      ind_cluster <- match(sort(unique(clusterIDObsk)), sort(unique(clusterIDObs)))

      #Stratum level information
      ind_stratum <- match(sort(unique(stratumIDObsk)), sort(unique(stratumIDObs)))

      # Unknown parameters that do not need to be redrawn
      bk <- b[ind_sub, ]
      sigma2k <- sigma2[k]
      betak <- beta[ , k]
      rhok  <- rho[ind_cluster, k]
      zetak  <- zeta[ind_stratum, k]

      sp_Vr_sub <- lapply(split(as.data.frame(Vrk), subjectIDk, drop = TRUE), as.matrix)
      convert_Vr_sub <- as.matrix(Matrix::bdiag(sp_Vr_sub))

      # Calculate linear predictor and draw Y for class k
      #Expand rho and zeta to observation levels for class k
      rhokObs <- rhok[factor(clusterIDObsk)]
      zetakObs <- zetak[factor(stratumIDObsk)]

      ### Draw Ypred for this iteration l
      lp <- as.vector(Vfk %*% betak + as.vector(convert_Vr_sub %*% c(t(bk))) + rhokObs + zetakObs)
      Ydraw[ind_obs] <- lp + rnorm(length(lp), 0, sqrt(sigma2k))

    }

  return(Ydraw)

}
