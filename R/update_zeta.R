#' @title Update the stratum level intercepts in the longitudinal model.
#'
#' @description This uses a normal prior distribution with mean 0 to update the stratum level intercepts in the longitudinal outcomes model.
#' @param C Latent class assignments for each subject.
#' @param b  A \code{q} column matrix of subject specific effects.
#' @param rho  A \code{K} column matrix of area segment (aka cluster) level intercepts.
#' @param beta A \code{K} column matrix of regression coefficients associated with \code{Vf}.
#' @param sigma2 A \code{K} length vector of observation level variances.
#' @param Psi A \code{K} length vector of stratum level variances.
#' @param Y Longitudinal measurements.
#' @param Vf A design matrix for fixed effects in the longitudinal outcome model.
#' @param Vr A design matrix for the subject-level random effects in the longitudinal outcome model. \code{Vr} is a subset or equal to \code{Vf}.
#' @param subjectID Subject identifier for each observation.
#' @param clusterIDObs Area segment or cluster identifier for each observation.
#' @param stratumIDObs Stratum identifier for each observation.
#' @export
#' @return A \code{K} column matrix of stratum level intercepts.
update_zeta <- function(C, b, rho, beta, sigma2, Psi, Y, Vf, Vr, subjectID, clusterIDObs, stratumIDObs){

  K <- length(table(C))
  M <- length(unique(stratumIDObs))
  C_expand <- C[factor(subjectID)]

  # Storage
  values <- matrix(NA, nrow = M, ncol = K)

  for (k in 1:K) {

    ### Observation level information including cluster identifier of each observation
    ind_obs <- which(C_expand == k)
    subjectIDk <- subjectID[ind_obs]
    clusterIDObsk <- clusterIDObs[ind_obs]
    stratumIDObsk <- stratumIDObs[ind_obs]
    Yk <- Y[ind_obs]
    Vrk <- Vr[ind_obs, ]
    Vfk <- Vf[ind_obs, ]

    # Construct N x n*q matrix of X to facilitate matrix algebra
    #[1 1 0 0 0 / 1 2 0 0 0 / 1 3 0 0 0], for q=2 columns 1 and 2 are subject 1's intercept and random slope for time
    sp_Vr_sub <- lapply(split(as.data.frame(Vrk), subjectIDk, drop = TRUE), as.matrix)
    convert_Vr_sub <- as.matrix(Matrix::bdiag(sp_Vr_sub)) # N x nq
    sigma2k <- sigma2[k]

    ### Subjects level information
    ind_sub <- which(C == k)
    bk <- b[ind_sub, ] # Here we choose class k
    betak <- beta[ , k]

    ### Cluster level information
    ind_cluster <- match(sort(unique(clusterIDObsk)), sort(unique(clusterIDObs)))
    rhok  <- rho[ind_cluster, k]
    rhokObs <- rhok[factor(clusterIDObsk)]

    Psik <- Psi[k]

    ### Compute posterior vairance
    #Nj is the number of *observations in a stratum of class k
    stratumIDObskf <- factor(stratumIDObsk, levels = levels(factor(stratumIDObs)))
    Nj <- as.vector(table(stratumIDObskf))
    post_var <- solve(diag(Nj * as.vector(solve(sigma2k)) + as.vector(solve(Psik)), nrow = M, ncol = M))

    ### Likelihood contribution
    llikObs <- (Yk -  Vfk %*% betak - as.vector(convert_Vr_sub %*% c(t(bk))) - rhokObs) / sigma2k
    # Sum likelihood over observations in each cluster tapply will order by factor level
    llik <- as.vector(tapply(llikObs, stratumIDObskf, sum))
    # Assign stratum that do not appear in class k and therefore do not have a likelihood contribution a value of zero
    llik[is.na(llik)] <- 0

    ### Prior information (cancels because blus \sim N(0, \Omega_k)
    prior <- 0 / Psik

    ### Posterior mean
    post_mean <- post_var %*% as.matrix(llik + prior)

    values[ , k] <- rnorm(M, mean = post_mean, sd = sqrt(diag(post_var)))

  }

  return(values)

}


