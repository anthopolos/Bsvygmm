#' @title Simulate data for growth mixture model for complex survey data.
#'
#' @description Data are simulated according to the growth mixture model for complex survey data.
#'
#' @param K Number of latent classes.
#' @param delta A \code{K-1} column matrix of regression coefficients for the latent class membership model.
#' @param gamma2 A \code{K-1} vector of stratum level variances for the latent class membership model.
#' @param tau2 A \code{K-1} vector of area segment level variances for the latent class membership model.
#' @param xi2 A \code{K-1} vector of area segment level spatial variances for the latent class membership model.
#' @param beta A \code{K} column matrix of regression coefficients associated with \code{Vf}.
#' @param Phi Variance-covariance of the subject specific effects.
#' @param Omega Variance at the area segment or cluster level.
#' @param Psi Varianec at the stratum level.
#' @param sigma2 Observation level variance.
#' @param W A design matrix of an intercept and covariates for the latent class membership model.
#' @param ADJ A binary symmetric adjacency matrix. Equal to \code{NULL} if spatial correlations are not included in the latent class membership model.
#' @param Vf A design matrix for fixed effects in the longitudinal outcome model.
#' @param Vr A design matrix for the subject-level random effects in the longitudinal outcome model. \code{Vr} is a subset or equal to \code{Vf}.
#' @param subjectID Subject identifier for each observation.
#' @param clusterIDObs Area segment or cluster identifier for each observation.
#' @param stratumIDObs Stratum identifier for each observation.
#' @param clusterIDSub Area segment or cluster identifier for each subject.
#' @param stratumIDSub Stratum identifier for each subject.
#' @param modelType A string. For the latent class memberhsip model, if \code{modelType} is \code{Unstr}, then the latent class membership model includes only cluster level random effects that account for within-cluster correlations. If \code{Str}, then the latent class membership model includes only cluster level random effects that account for spatial correlations among clusters. If \code{Both}, then both types of random effects are included.
#' @export
#' @return A list with simulated longitudinal otucomes and latent class memberships.
simdat <- function(K, delta, gamma2, tau2, xi2, beta, Phi, Omega, Psi, sigma2, W, ADJ, Vf, Vr, subjectID, clusterIDObs, stratumIDObs, clusterIDSub, stratumIDSub, modelType){

  #----------------------------- Generate latent class membership
  ### Spatial information for multinomial or binary probit class membership model
  if (modelType %in% c("Both", "Str")) {

    # Degree of spatial autocorrelation
    rho <- 1
    # Number of spatial units in the study area, independent of whether there are subjects found in them
    nJ <- nrow(ADJ)
    # Number of neighbors for each spatial unit
    m <- rowSums(ADJ)

    # Q matrix, add small constant to make Q non singular
    Qgen <- diag(m) - rho*(ADJ == 1) + diag(0.01, nJ)

    # Draw nus for each class k from their joint distribution
    if (K > 2) {
      nu <- matrix(NA, nrow = nJ, ncol = (K - 1))

      for (k in 1:(K - 1)) {
        # Covariance of nus
        cov.nu <- xi2[k]*solve(Qgen)

        # Generate spatial random effects from their joint distribution
        nu[ , k] <- mnormt::rmnorm(1, varcov = cov.nu)
        nu[ , k] <- nu[ , k] - mean(nu[ , k])
      }

    } else if (K == 2) {
      cov.nu <- xi2*solve(Qgen)

      # Generate spatial random effects from their joint distribution
      nu <- mnormt::rmnorm(1, varcov = cov.nu)
      nu <- nu - mean(nu)

    }

  } else {

    nu <- NULL

  }


  ### Background for latent class membership model
  n <- nrow(W)
  J <- length(unique(clusterIDSub))
  M <- length(unique(stratumIDSub))

  if (modelType %in% c("Both", "Str")) {
    clusterIDSubf <- factor(clusterIDSub, levels = seq(1, nJ, by = 1))
  }

  # Generate C for K > 2
  if (K > 2) {

    # Define stratum level random effects
    deltaStratum <- matrix(NA, nrow = M, ncol = (K - 1))
    for (k in 1:(K - 1)) {
      deltaStratum[ , k] <- rnorm(M, mean = 0, sd = sqrt(gamma2[k]))
    }

    # Cluster level unstructured heterogeniety or not?
    if (modelType %in% c("Both", "Unstr")) {
      u <- matrix(NA, nrow = J, ncol = (K - 1))
      for (k in 1:(K - 1)) {
        u[ , k] <- rnorm(J, mean = 0, sd = sqrt(tau2[k]))
      }
    } else {
      u <- NULL
    }


    # 1. Spatially structured heterogeneity at the cluster level only
    if (!is.null(nu) & is.null(u)) {
      Z <- mnormt::rmnorm(n, varcov = diag(K - 1)) + W %*% delta + deltaStratum[factor(stratumIDSub), ] + nu[clusterIDSubf, ]
    }

    # 2. Unstructured heterogeneity at the cluster level only
    if (is.null(nu) & !is.null(u)) {
        Z <- mnormt::rmnorm(n, varcov = diag(K - 1)) + W %*% delta + u[factor(clusterIDSub), ] + deltaStratum[factor(stratumIDSub), ]
    }

    # 3. Spatially structured and unstructured heterogeneity at the cluster level
    if (!is.null(nu) & !is.null(u)) {
      Z <- mnormt::rmnorm(n, varcov = diag(K - 1)) + W %*% delta + u[factor(clusterIDSub), ] + deltaStratum[factor(stratumIDSub), ] + nu[clusterIDSubf, ]
    }

    # Generate categorical variable based on max value of Z
    C <- rep(NA, n)
    for (i in 1:n) {
      MAX <- max(Z[i, ]) < 0  #All Z[i, ] < 0
      if (MAX) {
        #if the max is less than 0, then set C = 0
        C[i] <- 0
      } else {
        #otherwise, set C equal to the index of the maximum
        C[i] <- which.max(Z[i, ])
      }
    }
    C <- C + 1

  } # End of for K > 2

  # Generate data for K == 2
  if (K == 2) {

    # Define stratum level random effects
    deltaStratum <- rnorm(M, mean = 0, sd = sqrt(gamma2))

    # Cluster level unstructured heterogeniety or not?
    if (modelType %in% c("Both", "Unstr")) {
      u <- rnorm(J, mean = 0, sd = sqrt(tau2))
    } else {
      u <- NULL
    }

    # 1. Spatially structured heterogeneity at the cluster level only
    if (!is.null(nu) & is.null(u)) {
      Z <- mnormt::rmnorm(n, varcov = diag(K - 1)) + W %*% delta + deltaStratum[factor(stratumIDSub)] + nu[clusterIDSubf]
    }


    # 2. Unstructured heterogeneity at the cluster level only
    if (is.null(nu) & !is.null(u)) {
      Z <- mnormt::rmnorm(n, varcov = diag(K - 1)) + W %*% delta + u[factor(clusterIDSub)] + deltaStratum[factor(stratumIDSub)]
    }

    # 3. Spatially structured and unstructured heterogeneity at the cluster level
    if (!is.null(nu) & !is.null(u)) {
      Z <- mnormt::rmnorm(n, varcov = diag(K - 1)) + W %*% delta + u[factor(clusterIDSub)] + deltaStratum[factor(stratumIDSub)] + nu[clusterIDSubf]
    }


    C <- ifelse(Z > 0, 1, 0)
    C <- C + 1

  } # End of for K = 2


  #------------- Longitudinal outcomes generation
  # Expand C for data generation at the observation level
  C_expand <- C[factor(subjectID)]

  N <- length(subjectID)
  Y <- rep(NA, N)
  q <- ncol(Vr)

  for (k in 1:K) {

    ### Observation level information
    ind_obs <- which(C_expand == k)
    subjectIDk <- subjectID[ind_obs]
    clusterIDObsk <- clusterIDObs[ind_obs]
    stratumIDObsk <- stratumIDObs[ind_obs]
    Vfk <- Vf[ind_obs, ]
    Vrk <- Vr[ind_obs, ]

    ### Subject level information
    ind_sub <- which(C == k)

    ### Cluster level information
    ind_cluster <- match(sort(unique(clusterIDObsk)), sort(unique(clusterIDObs)))

    ### Stratum level information
    ind_stratum <- match(sort(unique(stratumIDObsk)), sort(unique(stratumIDObs)))

    ### Generate stratum level random effects for class k
    zetak <- rnorm(length(ind_stratum), mean = 0, sd = sqrt(Psi[k]))

    ### Generate cluster level random effects for class k
    rhok <- rnorm(length(ind_cluster), mean = 0, sd = sqrt(Omega[k]))

    ### Resultant bk is q by nk
    bk <- sapply(1:length(ind_sub), function(x){
        mnormt::rmnorm(1, mean = rep(0, q), varcov = Phi[ , , k])})

    ### Draw Y
    # Stack random effects i.e., a n*q x 1 vector of \eta_11 \eta_12 ... subject 1 random intercept, subject 1 random slope
    stackb <- c(bk)
    # Convert Z to (Nk) x (nk*q) i.e., total number of measurments sum(n_i) by number of unique subjects*number of random effects
    sp_n_sub <- lapply(split(as.data.frame(Vrk), subjectIDk, drop = TRUE), as.matrix)
    convert_n_sub <- as.matrix(Matrix::bdiag(sp_n_sub))

    # Expand rho to observation levels
    rhokObs <- rhok[factor(clusterIDObsk)]

    # Expand zeta to observation levels
    zetakObs <- zetak[factor(stratumIDObsk)]

    # Calculate linear predictor and draw Y
    lp <- Vfk %*% beta[ , k] + as.vector(convert_n_sub %*% stackb) + rhokObs + zetakObs
    Ytemp <- lp + rnorm(length(as.vector(lp)), 0, sqrt(sigma2[k]))

    Y[ind_obs] <- Ytemp

  }

  list(Y = Y, C = C)

}
