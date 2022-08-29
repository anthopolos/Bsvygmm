#' @title Get posterior predictive draws and level 1 and level 2 residuals at each iteration.
#'
#' @description This approach is based on the following references: Wang, Chen-Pin, C. Hendricks Brown, and Karen Bandeen-Roche.
#' 1. "Residual diagnostics for growth mixture models: Examining the impact of a preventive intervention on multiple trajectories of aggressive behavior." Journal of the American Statistical Association 100.471 (2005): 1054-1076.
#' 2. Gelman, Andrew, et al. "Multiple imputation for model checking: completed‐data plots with missing and latent data." Biometrics 61.1 (2005): 74-85.
#'
#' @export
#' @return Draws of residuals at levels 1 and 2 for complete and replicated data, along with longitudinal outcome and latent class membership at the current MCMC iteration. Returns a data.frame with these values in addition to subject ID, cluster ID and stratum ID.
get_post_pred_redraw <- function(K, C, b, sigma2, beta, zeta, rho, Phi, Omega, Psi, Z, delta, alpha, deltaStratum, u, nu, gamma2, tau2, xi2, priorPik,  W, B, ADJ, clusterIDSub, stratumIDSub, spline,  subjectID, clusterIDObs, stratumIDObs, Y, Vf, Vr, LRModelType, LCClusterModelType, LCStratumModelType) {

  #---------------- Complete data residuals
  ### Compute residuals under replication of complete data
  C_expand <- C[factor(subjectID)]

  # At each iteration store Ydraws from each latent class
  res_level_1 <- rep(NA, length(subjectID))
  res_level_2 <- rep(NA, length(subjectID))

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

    ### Compute residual rep at iteration l
    res_level_1[ind_obs] <- Y[ind_obs] - lp
    res_level_2[ind_obs] <- as.vector(bkObs)

  }


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


  #--------------------- Replicated complete data residuals
  C_rep_expand <- C_rep[factor(subjectID)]

  # At each iteration store Ydraws from each latent class
  Y_rep <- rep(NA, length(subjectID))
  res_level_1_rep <- rep(NA, length(subjectID))
  res_level_2_rep <- rep(NA, length(subjectID))

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
    Y_rep[ind_obs] <- lp + rnorm(length(lp), 0, sqrt(sigma2k))

    ### Compute residual rep at iteration l
    res_level_1_rep[ind_obs] <- Y_rep[ind_obs] - lp
    res_level_2_rep[ind_obs] <- as.vector(b_repkObs)

  }

  ### Print to file complete data and replicated complete data draw
  res <- data.frame(rl1 = res_level_1, rl1_rep = res_level_1_rep, rl2 = res_level_2, rl2_rep = res_level_2_rep, C = C_expand, C_rep = C_rep_expand, Y = Y, Y_rep = Y_rep, subjectID = subjectID, clusterID = clusterIDObs, stratumID = stratumIDObs)
  res <- cbind(Vf, res)
  colnames(res)[1:dim(Vf)[2]] <- paste0("time", 1:dim(Vf)[2])

  return(res)

}
