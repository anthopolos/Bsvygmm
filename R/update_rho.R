#' @title Update the area segment (aka cluster) level intercepts in the longitudinal model.
#'
#' @description This uses a normal prior distribution with mean 0 to update the arae segment level intercepts in the longitudinal outcomes model.
#'
#' @param C Latent class assignments for each subject.
#' @param b  A \code{q} column matrix of subject specific effects.
#' @param zeta  A \code{K} column matrix of stratum level intercepts.
#' @param beta A \code{K} column matrix of regression coefficients associated with \code{Vf}.
#' @param sigma2 A \code{K} length vector of observation level variances.
#' @param Omega A \code{K} length vector of area segment level variances.
#' @param Y Longitudinal measurements.
#' @param Vf A design matrix for fixed effects in the longitudinal outcome model.
#' @param Vr A design matrix for the subject-level random effects in the longitudinal outcome model. \code{Vr} is a subset or equal to \code{Vf}.
#' @param subjectID Subject identifier for each observation.
#' @param clusterIDObs Area segment or cluster identifier for each observation.
#' @param stratumIDObs Stratum identifier for each observation.
#' @export
#' @return A \code{K} column matrix of area segment level intercepts.
update_rho <- function(C, b, zeta, beta, sigma2, Omega, Y, Vf, Vr, subjectID, clusterIDObs, stratumIDObs){

  K <- length(table(C))
  J <- length(unique(clusterIDObs))
  C_expand <- C[factor(subjectID)]

  # Storage
  values <- matrix(NA, nrow = J, ncol = K)

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

    ### Stratum level information
    ind_stratum <- match(sort(unique(stratumIDObsk)), sort(unique(stratumIDObs)))
    zetak  <- zeta[ind_stratum, k]
    zetakObs <- zetak[factor(stratumIDObsk)]

    Omegak <- Omega[k]

    ### Compute posterior vairance
    #Nj is the number of *observations in cluster j of class k
    clusterIDObskf <- factor(clusterIDObsk, levels = levels(factor(clusterIDObs)))
    Nj <- as.vector(table(clusterIDObskf))
    post_var <- solve(diag(Nj * as.vector(solve(sigma2k)) + as.vector(solve(Omegak)), nrow = J, ncol = J))

    ### Likelihood contribution
    llikObs <- (Yk -  Vfk %*% betak - as.vector(convert_Vr_sub %*% c(t(bk))) - zetakObs) / sigma2k
    # Sum likelihood over observations in each cluster tapply will order by factor level
    llik <- as.vector(tapply(llikObs, clusterIDObskf, sum))
    # Assign clusters that do not appear in class k and therefore do not have a likelihood contribution a value of zero
    llik[is.na(llik)] <- 0

    ### Prior information (cancels because blus \sim N(0, \Omega_k)
    prior <- 0 / Omegak

    ### Posterior mean
    post_mean <- post_var %*% as.matrix(llik + prior)

    values[ , k] <- rnorm(J, mean = post_mean, sd = sqrt(diag(post_var)))

  }

  return(values)

}


