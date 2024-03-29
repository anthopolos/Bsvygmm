#' @title Get posterior predictive draws of discrepancy the weighted mean squared error at each iteration.
#'
#' @description See references for the discrepancy measure in manuscript. The algorithm used is based on redrawing the complete data. This is consistent with the function \code{get_post_pred_redraw()}.
#'
#' @export
#' @return Draws of the discrepancy measure computed at the observed data and replicated data.
get_discrepancy <- function(K, C, b, sigma2, beta, zeta, rho, Phi, Omega, Psi, Z, delta, alpha, deltaStratum, u, nu, gamma2, tau2, xi2, priorPik,  W, B, ADJ, clusterIDSub, stratumIDSub, spline,  subjectID, clusterIDObs, stratumIDObs, Y, Vf, Vr, LRModelType, LCClusterModelType, LCStratumModelType) {

  #---------------- Complete data discrepancy measure
  C_expand <- C[factor(subjectID)]

  ### TObs
  TObs <- 0

  for (k in 1:K) {

    ### Subset information based on the Cdraw
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
    bkObs <- as.vector(convert_Vr_sub %*% c(t(bk)))
    lp <- as.vector(Vfk %*% betak + bkObs)

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

    ### Compute discrepancy measure
    TObsk <- (Y[ind_obs] - lp)^2 / sigma2k
    TObs <- TObs + sum(TObsk)

  } # Ends K-loop for complete data discrepancy measure


  #----------------- Replicate complete data
  ### Latent class membership model re-draw latent variables
  # Stratum level heterogeniety
  if (LCStratumModelType == "Random") {
    if (K > 2) {
      lambda_rep <- update_deltaStratum_MNP(Z = Z, delta = delta, alpha = alpha, u = u, nu = nu, gamma2 = gamma2, W = W, B = B, clusterIDSub = clusterIDSub, stratumIDSub = stratumIDSub, spline = spline)
    } else {
      lambda_rep <- update_deltaStratum_BinaryP(Z = Z, delta = delta, alpha = alpha, u = u, nu = nu, gamma2 = gamma2, W = W, B = B, clusterIDSub = clusterIDSub, stratumIDSub = stratumIDSub, spline = spline)
    }
  } else {
    lambda_rep <- NULL
  }

  # Unstructure random intercepts and variance terms for cluster and stratum levels, hierarchically centered
  if (LCClusterModelType %in% c("Both", "Unstr")) {
    if (K > 2) {
      u_rep <- update_u_MNP(Z = Z, delta = delta, alpha = alpha, deltaStratum = lambda_rep, nu = nu, tau2 = tau2, W = W, B = B, clusterIDSub = clusterIDSub, stratumIDSub = stratumIDSub, spline = spline)
    } else {
      u_rep <- update_u_BinaryP(Z = Z, delta = delta, alpha = alpha, deltaStratum = lambda_rep, nu = nu, tau2 = tau2, W = W, B = B, clusterIDSub = clusterIDSub, stratumIDSub = stratumIDSub, spline = spline)
    }
  } else {
    u_rep <- NULL
  }


  # Spatial random intercept and variance term
  if (LCClusterModelType %in% c("Both", "Str")) {
    nu_rep <- update_nu(Z = Z, delta = delta, alpha = alpha, deltaStratum = lambda_rep, u = u_rep, nu = nu, xi2 = xi2, W = W, B = B, ADJ = ADJ, clusterIDSub = clusterIDSub, stratumIDSub = stratumIDSub, spline = spline)
  } else {
    nu_rep <- NULL
  }

  # Latent variable Z and regression coefficients delta depending on K > 2 or not
  if (K > 2) {
    Z_rep <- update_Z_MNP(Z = Z, C = C, delta = delta, alpha = alpha, deltaStratum = lambda_rep, u = u_rep, nu = nu_rep, W = W, B = B, clusterIDSub = clusterIDSub, stratumIDSub = stratumIDSub, spline = spline)
  } else {
    Z_rep <- update_Z_BinaryP(Z = Z, C = C, delta = delta, alpha = alpha, deltaStratum = lambda_rep, u = u_rep, nu = nu_rep, W = W, B = B, clusterIDSub = clusterIDSub, stratumIDSub = stratumIDSub, spline = spline)
  }

  # Update priorPik
  priorPik_rep <- update_priorPik(Z = Z_rep)


  ### Re-draw latent class memberships
  postClass_rep <- update_C_Marginal(priorPik = priorPik_rep, beta = beta, Phi = Phi, Omega = Omega, Psi = Psi, sigma2 = sigma2, Y = Y, Vf = Vf, Vr = Vr, subjectID = subjectID)
  C_rep <- postClass_rep[["C"]]


  ### Longitudinal outcomes model redraw latent heterogeneity
  # Stratum level random effects
  if (LRModelType %in% c("Both", "Stratum")) {
    zeta_rep <- update_zeta(C = C_rep, b = b, rho = rho, beta = beta, sigma2 = sigma2, Psi = Psi, Y = Y, Vf = Vf, Vr = Vr, subjectID = subjectID, clusterIDObs = clusterIDObs, stratumIDObs = stratumIDObs)
  } else {
    zeta_rep <- NULL
  }

  # Cluster level random effects
  if (LRModelType %in% c("Both", "Cluster")) {
    rho_rep <- update_rho(C = C_rep, b = b, zeta = zeta_rep, beta = beta, sigma2 = sigma2, Omega = Omega, Y = Y, Vf = Vf, Vr = Vr, subjectID = subjectID, clusterIDObs = clusterIDObs, stratumIDObs = stratumIDObs)
  } else {
    rho_rep <- NULL
  }

  # Subject level random effects
  b_rep <- update_b(C = C_rep, rho = rho_rep, zeta = zeta_rep, beta = beta, Phi = Phi, sigma2 = sigma2, Y = Y, Vf = Vf, Vr = Vr, subjectID = subjectID, clusterIDObs = clusterIDObs, stratumIDObs = stratumIDObs)


  #--------------------- Replicated complete data discrepancy
  C_rep_expand <- C_rep[factor(subjectID)]

  # At each iteration store Ydraws from each latent class
  TRep <- 0

  for (k in 1:K) {

    ### Subset information based on the Cdraw
    #Observation level information
    ind_obs <- which(C_rep_expand == k)
    subjectIDk <- subjectID[ind_obs]
    clusterIDObsk <- clusterIDObs[ind_obs]
    stratumIDObsk <- stratumIDObs[ind_obs]
    Vfk <- Vf[ind_obs, ]
    Vrk <- Vr[ind_obs, ]

    #Subject level information
    ind_sub <- which(C_rep == k)

    #Cluster level information
    ind_cluster <- match(sort(unique(clusterIDObsk)), sort(unique(clusterIDObs)))

    #Stratum level information
    ind_stratum <- match(sort(unique(stratumIDObsk)), sort(unique(stratumIDObs)))

    # Unknown parameters that do not need to be redrawn
    b_repk <- b_rep[ind_sub, ]
    sigma2k <- sigma2[k]
    betak <- beta[ , k]

    sp_Vr_sub <- lapply(split(as.data.frame(Vrk), subjectIDk, drop = TRUE), as.matrix)
    convert_Vr_sub <- as.matrix(Matrix::bdiag(sp_Vr_sub))

    # Calculate linear predictor and draw Y for class k
    b_repkObs <- as.vector(convert_Vr_sub %*% c(t(b_repk)))
    lp <- as.vector(Vfk %*% betak + b_repkObs)

    # Include stratum level intercepts or not
    #Expand rho and zeta to observation levels for class k
    if (!is.null(zeta_rep)) {
      zeta_repk  <- zeta_rep[ind_stratum, k]
      zeta_repkObs <- zeta_repk[factor(stratumIDObsk)]
      lp <- lp + zeta_repkObs
    }

    # Include cluster level intercepts or not
    if (!is.null(rho_rep)) {
      rho_repk  <- rho_rep[ind_cluster, k]
      rho_repkObs <- rho_repk[factor(clusterIDObsk)]
      lp <- lp + rho_repkObs
    }
    ### Draw Y_rep for this iteration l
    Y_repk <- lp + rnorm(length(lp), 0, sqrt(sigma2k))

    ### Compute residual rep at iteration l
    TRepk <- (Y_repk - lp)^2 / sigma2k
    TRep <- TRep + sum(TRepk)

  } # Closes K loop for replicated complete discrepancy measure

  ### Store discrepancy measures
  store_T <- c(TObs, TRep)

  return(store_T)

}

