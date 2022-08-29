#' @title Simulate data for growth mixture model for complex survey data.
#'
#' @description Data are simulated according to the growth mixture model for complex survey data.
#'
#' @param K Number of latent classes.
#' @param data A data frame in long format for the analysis variables in the growth mixture model.
#' @param forms A named list with 3 elements named \code{W}, \code{Xf}, and \code{Xr}. \code{W} must be a formula for the latent class membership model (e.g., \code{ ~ female + NHW}). \code{Xf} specifies the design matrix for the fixed effects. It may equal a character string \code{dummy} if the user desires a design matrix with dummy indicators for each time point. Otherwise, \code{Xf} will equal a formula (e.g., \code{~ time}) where an intercept will be generated in the call.\code{Xr} is the design matrix for the random effects. It may equal the character string \code{same} if \code{Xr} is the same as \code{Xf}. Otherwise, \code{Xr} will equal a formula.
#' @param delta A \code{K-1} column matrix of regression coefficients for the latent class membership model.
#' @param gamma2 A \code{K-1} vector of stratum level variances for the latent class membership model.
#' @param tau2 A \code{K-1} vector of area segment level variances for the latent class membership model.
#' @param xi2 A \code{K-1} vector of area segment level spatial variances for the latent class membership model.
#' @param beta A \code{K} column matrix of regression coefficients associated with \code{Vf}.
#' @param Phi Variance-covariance of the subject specific effects.
#' @param Omega Variance at the area segment or cluster level. Equals \code{NULL} if not desired.
#' @param Psi Variance at the stratum level. Equals \code{NULL} if not desired.
#' @param sigma2 Observation level variance.
#' @param ADJ A binary symmetric adjacency matrix. Equal to \code{NULL} if spatial correlations are not included in the latent class membership model.
#' @param subjectID A string for the name of the subject identifier in \code{data}.
#' @param time A string for the name of the time variable in the dataset.
#' @param clusterID A string for the name of the area segment or cluster identifier in \code{data}.
#' @param stratumID A string for the name of the stratum identifier in \code{data}.
#' @param LCStratumModelType A string. For the latent class membership model, if \code{LCStratumModelType} is \code{Random}, then the latent class membership model includes stratum level random effects. If \code{Fixed}, then stratum level random effects are not included in the latent class membership model. Stratum fixed effects may be included in the design matrix \code{W}.
#' @param LCClusterModelType A string. For the latent class memberhsip model, if \code{LCModelType} is \code{Unstr}, then the latent class membership model includes only cluster level random effects that account for within-cluster correlations. If \code{Str}, then the latent class membership model includes only cluster level random effects that account for spatial correlations among clusters. If \code{Both}, then both types of random effects are included. If \code{None}, then no random effects are included.
#' @param LRModelType A string. For the longitudinal response model, if \code{LRModelType} is \code{Stratum}, then the longitudinal response model includes only stratum level random effects that account for within-stratum correlations. If \code{Cluster}, then the longitudinal response model includes only cluster level random effects that account for within-cluster correlations. If \code{Both}, then both types of random effects are included. If neither are desired specify as \code{None}.
#' @export
#' @return A list with simulated longitudinal otucomes and latent class memberships.
simdat <- function(K, data, forms, delta, gamma2, tau2, xi2, beta, Phi, Omega, Psi, sigma2, ADJ, subjectID, time, clusterID, stratumID, LCStratumModelType, LCClusterModelType, LRModelType){

  #----- Background
  data$subjectID <- data[ , subjectID]
  data$clusterID <- data[ , clusterID]
  data$stratumID <- data[ , stratumID]
  data$time <- data[ , time]

  if (forms[["Xf"]] == "dummy" & !is.factor(data$time)) {
    data$time <- factor(data$time)
  }


  subjectID <- data$subjectID
  clusterIDObs <- data$clusterID
  stratumIDObs <- data$stratumID

  # Latent class membership model design matrix
  var_namess <- all.vars(forms[["W"]])
  datas <- aggregate(data[ , c("subjectID", "clusterID", "stratumID", var_namess)], by = list(data$subjectID), FUN = tail, n = 1)
  W <- model.matrix(forms[["W"]], data = datas)

  # Design matrix for fixed effects in longitudinal outcomes model
  if (forms[["Xf"]] == "dummy") {
    Vf <- data.matrix(dummy::dummy(data.frame(data[ , "time"]), int = TRUE))
  } else {
    Vf <- model.matrix(forms[["Xf"]], data = data)
  }

  # Design matrix for fixed effects in longitudinal outcomes model
  if (forms[["Xr"]] == "same") {
    Vr <- Vf
  } else {
    Vr <- model.matrix(forms[["Xr"]], data = data)
  }

  clusterIDSub <- datas$clusterID
  stratumIDSub <- datas$stratumID

  #----------------------------- Generate latent class membership
  ### Spatial information for multinomial or binary probit class membership model
  if (LCClusterModelType %in% c("Both", "Str")) {

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

  if (LCClusterModelType %in% c("Both", "Str")) {
    clusterIDSubf <- factor(clusterIDSub, levels = seq(1, nJ, by = 1))
  }

  # Generate C for K > 2
  if (K > 2) {

    # K-1 by K-1 covariance matrix
    VarCovP <- matrix(1, nrow = (K - 1), ncol = (K - 1))
    diag(VarCovP) <- 2

    # Errors for each z_i1,..,z_{iK-1}
    err <- do.call("rbind", sapply(1:n, function(x) {res <- mnormt::rmnorm(1, mean = rep(0, (K-1)), varcov = VarCovP)
    return(res) }, simplify = FALSE))


    # Define stratum level random effects
    if (LCStratumModelType == "Random") {
      deltaStratum <- matrix(NA, nrow = M, ncol = (K - 1))
      for (k in 1:(K - 1)) {
        deltaStratum[ , k] <- rnorm(M, mean = 0, sd = sqrt(gamma2[k]))
      }
    } else {
      deltaStratum <- NULL
    }

    # Cluster level unstructured heterogeniety or not?
    if (LCClusterModelType %in% c("Both", "Unstr")) {
      u <- matrix(NA, nrow = J, ncol = (K - 1))
      for (k in 1:(K - 1)) {
        u[ , k] <- rnorm(J, mean = 0, sd = sqrt(tau2[k]))
      }
    } else {
      u <- NULL
    }


    # 1. Spatially structured heterogeneity at the cluster level only
    if (!is.null(nu) & is.null(u)) {
      Z <- err + W %*% delta  + nu[clusterIDSubf, ]
    }

    # 2. Unstructured heterogeneity at the cluster level only
    if (is.null(nu) & !is.null(u)) {
      Z <- err + W %*% delta + u[factor(clusterIDSub), ]
    }

    # 3. Spatially structured and unstructured heterogeneity at the cluster level
    if (!is.null(nu) & !is.null(u)) {
       Z <- err + W %*% delta + u[factor(clusterIDSub), ] + nu[clusterIDSubf, ]
    }

    # 4. No spatially structured and unstructured heterogeneity at the cluster level
    if (is.null(nu) & is.null(u)) {
       Z <-  err + W %*% delta
    }

    # Stratum level intercept if desired
    if (!is.null(deltaStratum)) {
      Z <- Z + deltaStratum[factor(stratumIDSub), ]
    }

    # Generate categorical variable based on max value of Z
    # Reference group is column 1
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
    if (LCStratumModelType == "Random") {
      deltaStratum <- rnorm(M, mean = 0, sd = sqrt(gamma2))
    } else {
      deltaStratum <- NULL
    }

    # Cluster level unstructured heterogeniety or not?
    if (LCClusterModelType %in% c("Both", "Unstr")) {
      u <- rnorm(J, mean = 0, sd = sqrt(tau2))
    } else {
      u <- NULL
    }

    # 1. Spatially structured heterogeneity at the cluster level only
    if (!is.null(nu) & is.null(u)) {
      Z <- mnormt::rmnorm(n, varcov = diag(K - 1)) + W %*% delta + nu[clusterIDSubf]
    }

    # 2. Unstructured heterogeneity at the cluster level only
    if (is.null(nu) & !is.null(u)) {
      Z <- mnormt::rmnorm(n, varcov = diag(K - 1)) + W %*% delta + u[factor(clusterIDSub)]
    }

    # 3. Spatially structured and unstructured heterogeneity at the cluster level
    if (!is.null(nu) & !is.null(u)) {
      Z <- mnormt::rmnorm(n, varcov = diag(K - 1)) + W %*% delta + u[factor(clusterIDSub)] + nu[clusterIDSubf]
    }

    # 4. No spatially structured and unstructured heterogeneity at the cluster level
    if (is.null(nu) & is.null(u)) {
      Z <- mnormt::rmnorm(n, varcov = diag(K - 1)) + W %*% delta
    }


    # Stratum level intercept if desired
    if (!is.null(deltaStratum)) {
      Z <- Z + deltaStratum[factor(stratumIDSub)]
    }

    C <- ifelse(Z > 0, 1, 0)
    C <- C + 1

  } # End of for K = 2


  #------------- Longitudinal outcomes generation
  # Expand C for data generation at the observation level
  C_expand <- C[factor(subjectID)]

  N <- length(subjectID)
  data$Y <- rep(NA, N)
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

    ### Resultant bk is q by nk
    bk <- sapply(1:length(ind_sub), function(x){
        mnormt::rmnorm(1, mean = rep(0, q), varcov = Phi[ , , k])})

    ### Draw Y
    # Stack random effects i.e., a n*q x 1 vector of \eta_11 \eta_12 ... subject 1 random intercept, subject 1 random slope
    stackb <- c(bk)
    # Convert Z to (Nk) x (nk*q) i.e., total number of measurments sum(n_i) by number of unique subjects*number of random effects
    sp_n_sub <- lapply(split(as.data.frame(Vrk), subjectIDk, drop = TRUE), as.matrix)
    convert_n_sub <- as.matrix(Matrix::bdiag(sp_n_sub))

    # Calculate linear predictor and draw Y
    lp <- Vfk %*% beta[ , k] + as.vector(convert_n_sub %*% stackb)

    ### Include cluster and stratum level effects for each class or not
    # Generate stratum level random effects for class k
    if (!is.null(Psi) & LRModelType %in% c("Stratum", "Both")) {
      zetak <- rnorm(length(ind_stratum), mean = 0, sd = sqrt(Psi[k]))
      # Expand zeta to observation levels
      zetakObs <- zetak[factor(stratumIDObsk)]
      lp <- lp + zetakObs
    }

    # Generate cluster level random effects for class k
    if (!is.null(Omega) & LRModelType %in% c("Cluster", "Both")) {
      rhok <- rnorm(length(ind_cluster), mean = 0, sd = sqrt(Omega[k]))
      # Expand rho to observation levels
      rhokObs <- rhok[factor(clusterIDObsk)]
      lp <- lp + rhokObs
    }

    Ytemp <- lp + rnorm(length(as.vector(lp)), 0, sqrt(sigma2[k]))

    data$Y[ind_obs] <- Ytemp

  }

  list(data = data, C = C)

}
