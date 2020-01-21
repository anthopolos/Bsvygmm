#' @title Get discrepancy measured computed at observed and replicated data.
#'
#' @param K Number of latent classes.
#' @param C Categorical variable for the latent class of each subject.
#' @param b A \code{q} column matrix of subject-level effects.
#' @param rho A \code{K} column matrix of area segment level intercepts.
#' @param zeta A \code{K} column matrix of stratum level intercepts.
#' @param beta Current iteration of fixed effects.
#' @param sigma2 Current iteration of observation level variance.
#' @param priorPik Current iteration of prior probability of class membership.
#' @param Y Vector of longitudinal outcomes.
#' @param Vf Design matrix for fixed effects.
#' @param Vr Design matrix for random effects.
#' @param subjectID Subject identifier for each observation.
#' @param clusterIDObs Area segment or cluster identifier for each observation.
#' @param stratumIDObs Stratum identifier for each observation.
#' @return Draws of the discrepancy measure computed at the observed data and replicated data.
get_discrepancy <- function(K, C, b, rho, zeta, beta, sigma2, Y, Vr, Vf, subjectID, clusterIDObs, stratumIDObs) {

  C_expand <- C[factor(subjectID)]

  ### TObs, Rep, and RepCdraw to sum across classes
  TObs <- 0
  TRep <- 0

  for (k in 1:K) {

    ### Process MCMC sample l for each latent class
    # Subject level information
    ind_sub <- which(C == k)

    # Observation level information
    ind_obs <- which(C_expand == k)
    Yk <- Y[ind_obs]
    Vrk <- Vr[ind_obs, ]
    Vfk <- Vf[ind_obs, ]
    subjectIDk <- subjectID[ind_obs]
    clusterIDObsk <- clusterIDObs[ind_obs]
    stratumIDObsk <- stratumIDObs[ind_obs]

    # Cluster information
    ind_cluster <- match(sort(unique(clusterIDObsk)), sort(unique(clusterIDObs)))

    # Stratum information
    ind_stratum <- match(sort(unique(stratumIDObsk)), sort(unique(stratumIDObs)))

    # Unknown paramters
    bk <- b[ind_sub, ]
    sigma2k <- sigma2[k]
    betak <- beta[ , k]
    rhok  <- rho[ind_cluster, k]
    zetak  <- zeta[ind_stratum, k]

    rhokObs <- rhok[factor(clusterIDObsk)]
    zetakObs <- zetak[factor(stratumIDObsk)]

    sp_Vr_sub <- lapply(split(as.data.frame(Vrk), subjectIDk, drop = TRUE), as.matrix)
    convert_Vr_sub <- as.matrix(Matrix::bdiag(sp_Vr_sub)) # N x nq

    ### Compute discrepancy measure based on observed
    # Mean squared error page 167, Gelman BDA. Also used in Neelon 2011b (but reported with error)
    lp <- Vfk %*% betak + as.vector(convert_Vr_sub %*% c(t(bk))) + rhokObs + zetakObs

    TObsk <- (Yk - lp)^2 / sigma2k
    TObs <- TObs + sum(TObsk)

    ### Compute discrepancy measure based on Y^rep conditioning on latent class
    Ydrawk <- lp + rnorm(length(lp), 0, sqrt(sigma2k))
    TRepk <- (Ydrawk - lp)^2 / sigma2k

    TRep <- TRep + sum(TRepk)


  } # Closes K loop

  store_T <- c(TObs, TRep)

  return(store_T)

}

