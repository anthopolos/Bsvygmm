#' @title Bayesian growth mixture model for complex survey data.
#'
#' @description A Bayesian Growth Mixture Model (GMM) for complex survey data is fit using a Gibbs sampler with only closed-form full conditionals.
#'
#' The model was constructed for data with longitudinal measurements nested in subjects, subjects nested in area segments (also called clusters), and area segments nested in strata.
#'
#' A multinomial independence probit model is use to predict each subject's latent class membership. The multinomial independence probit model uses subjects as the unit of analysis and includes options for stratum and area segment random effects to model unobserved variations at the stratum and area segment levels. For the area segment random effects, the model includes options for modeling either or both correlations among subjects within the same area segment and spatial correlations among neighboring area segments. Correlations among subjects within the same area segment are modeled using independent random effects assumed to follow a normal distribution. Spatial correlations among neighboring area segments are modeled with an intrinsic conditional autoregessive prior distribution. An option is included to use B-splines for flexibly modeling the relationship between the latent class membership probabilities and continuously measured features of the complex sampel design (e.g., the size variable used in probability proportional to size sampling or nonresponse weighting adjustment).
#'
#' In the longitudinal outcomes model, latent class-specific longitudinal trajectories are estimated. Between-subject heterogeneity is modeled with subject-specific random intercepts and slopes assumed to follow a multivariate normal distribution. There are options to include area segment and stratum-specific independent random effects modeled using latent-class specific variance parameters. The random effects are assumed to follow a normal distribution. An option is included as to whether to include both area segment and stratum-specific random effects; one or the other; or neither. The observation-level variances are also latent class-specific.
#' @param K A scalar for the number of assumed unobserved, latent classes.
#' @param data A data.frame with the variables to be used in the latent class and longitudinal outcomes models.
#' @param forms A named list with 3 elements named \code{W}, \code{Xf}, and \code{Xr}. \code{W} must be a formula for the latent class membership model (e.g., \code{ ~ female + NHW}). \code{Xf} specifies the design matrix for the fixed effects. It may equal a character string \code{dummy} if the user desires a design matrix with dummy indicators for each time point. Otherwise, \code{Xf} will equal a formula (e.g., \code{~ time}) where an intercept will be generated in the call.\code{Xr} is the design matrix for the random effects. It may equal the character string \code{same} if \code{Xr} is the same as \code{Xf}. Otherwise, \code{Xr} will equal a formula.
#' @param ADJ A binary symmetric adjacency matrix. Equal to \code{NULL} if spatial correlations are not included in the latent class membership model.
#' @param Y Character for the name of the variable with the longitudinal measurements in \code{data}.
#' @param time Character for the name of the variable with the interview wave in integer values in \code{data}.
#' @param subjectID Character for the subject identifier for each observation in \code{data}.
#' @param clusterID Character for area segment or cluster identifier for each observation in \code{data}.
#' @param stratumID Character for stratum identifier for each observation in \code{data}.
#' @param spline Logical. If \code{TRUE}, a B-spline is included for the basis functions in argument \code{B} for the latent class membership model.
#' @param B A design matrix of basis functions for the B-spline in the latent class membership model.
#' @param LCStratumModelType A string. For the latent class membership model, if \code{LCStratumModelType} is \code{Random}, then the latent class membership model includes stratum level random effects. If \code{Fixed}, then stratum level random effects are not included in the latent class membership model. Stratum fixed effects may be included in the design matrix for the fixed effects \code{W}.
#' @param LCClusterModelType A string. For the latent class membership model, if \code{LCClusterModelType} is \code{Unstr}, then the latent class membership model includes only cluster level random effects that account for within-cluster correlations. If \code{Str}, then the latent class membership model includes only cluster level random effects that account for spatial correlations among clusters. If \code{Both}, then both types of random effects are included. If \code{None}, then neither random effects are included.
#' @param LRModelType A string. For the longitudinal response model, if \code{LRModelType} is \code{Stratum}, then the longitudinal response model includes only stratum level random effects that account for within-stratum correlations. If \code{Cluster}, then the longitudinal response model includes only cluster level random effects that account for within-cluster correlations. If \code{Both}, then both types of random effects are included. If neither stratum or cluster random effects are desired, specify as \code{None}.
#' @param priors A nested list of prior parameters. All regression coefficients are assigned a normal prior distribution. All scalar variance terms are assigned an inverse gamma prior distribution or a uniform prior distribution, according to user choice. The variance-covariance matrix for the subject level random effects is assigned an inverse-Wishart distribution. If any parameters are not needed in the model, the element in the list can be set to NULL. The prior distributions are not allowed to vary by latent class. Details are provided in the \code{Details} section.
#' @param hierVar A list with two string elements to indicate, in order, the type of prior distribution for the hierarchical variances for the random effects and for the observation level data variance. Each element may equal either \code{Unif} or \code{IG} for a uniform prior distribution or inverse gamma prior distribution, respectively.
#' @param inits A list of initial values. If any parameters are not needed in the model, the element in the list can be set to NULL. The initial values are allowed to vary by latent class. Details are provided in the \code{Details} section.
#' @param n.samples Number of MCMC iterations.
#' @param burn Number of MCMC samples to discard.
#' @param monitor Logical. If \code{TRUE}, then plots of the current iteration for the first four \code{beta} will appear, after burn-in.
#' @param update A scalar for the interval at which to print the current iteration. If monitor is \code{TRUE}, then trace plots will also be updated at this interval.
#' @param n.reps MCMC updating threshold at which to start printing selected posterior simulations to file.
#' @param writeSamples Logical. If \code{TRUE}, then MCMC samples of all of the random effects, draws from the posterior predictive distribution, and the discrepancy measure used to evaluate overall model fit are written to text files, after burn-in. The text files for each of the random effects are accompanied by corresponding text files with column name information.
#'
#' @details Details are given on the specification of \code{priors} and \code{inits}.
#'
#' \describe{
#' \item{A.	The nested list \code{priors} must contain the information on the prior distributions in the order as given below.}{
#' \enumerate{
#' \item Multivariate normal prior distribution on the regression coefficients in the latent class membership model. The list entry is of the form, for example, \code{list(rep(0, s), diag(10, s))}, where \code{s} is the number of regression coefficients. \code{rep(0, s)} assigns a prior mean of 0 to all regression coefficients. A diagonal variance-covariance matrix with variance 10 is assigned.
#'
#' \item	An inverse gamma distribution on the variance of the independent random effects in the latent class membership model at the stratum level. The form is, for example, \code{list(0.1, 0.1)}, where the first element is the shape parameter and the second parameter is the scale parameter. Set this prior equal to \code{NULL} if \code{hierVar[[1]] = Unif}, in which case a uniform prior distribution is used, or if stratum level random effects are not desired.
#'
#' \item An inverse gamma distribution on the variance of the spatial random effects in the latent class membership model at the cluster level. The form is the same as in 2. Set this prior equal to \code{NULL} if \code{hierVar[[1]] = Unif}, in which case a uniform prior distribution is used. If spatial random effects are not desired, then use \code{list(NULL)}.
#'
#' \item An inverse gamma distribution on the variance of the independent random effects in the latent class membership model at the cluster level. The form is the same as in 2. Set this prior equal to \code{NULL} if \code{hierVar[[1]] = Unif}, in which case a uniform prior distribution is used. If independent random effects are not desired, then use \code{list(NULL)}.
#'
#' \item If \code{spline = TRUE}, a multivariate normal prior distribution on the regression coefficients for the B-splines in the latent class membership model. The prior mean is fixed to 0. The user must specify the prior variance-covariance. The list entry is of the form, for example, \code{list(diag(10, r))}, where \code{r} is the number of B-spline coefficients. If \code{spline = FALSE}, then use \code{list(NULL)}.
#'
#' \item Multivariate normal prior distribution on the regression coefficients in the longitudinal outcomes model. The specification is the same as in 1.
#'
#' \item An inverse gamma distribution on the variance of the independent random effects in the longitudinal outcomes model at the stratum level. The specification is the same as in 2.
#'
#' \item An inverse gamma distribution on the variance of the independent random effects in the longitudinal outcomes model at the cluster level. The specification is the same as in 4.
#'
#' \item Inverse-Wishart prior distribution on the variance-covariance matrix of the subject level random effects in the longitudinal outcomes model. The form is, for example, \code{list(q + 2, diag(0.25, q))}, where \code{q} is the number of random effects. The first element is the prior degrees of freedom, and the second element is the prior scale matrix of dimension \code{q} by \code{q}.
#'
#' \item An inverse gamma distribution on the variance of the observation-level error in the longitudinal outcomes model. The specification is the same as in 2, except that this prior is equal to \code{NULL} if \code{hierVar[[2]] = Unif}, in which case a uniform prior distribution is used.}
#'
#' An example of the completed list object for specification of the prior distributions is: \code{priors <- list(list(rep(0, s), diag(1, s)), list(.1, .1), list(.1, .1), list(.1, .1), list(diag(10, r)), list(rep(0, p), diag(10, p)), list(.1, .1), list(.1, .1), list((q + 2), diag(0.25, q)), list(.1, .1))}
#' }
#' \item{B.	The list \code{inits} must contain the information on initial values for each of the parameters, in the same order as in \code{priors}.}{
#' \enumerate{
#' \item Initial values for the regression coefficients in the latent class membership model. The list entry is of the form, for example, \code{matrix(rep(0, s * (K - 1)), nrow = s, ncol = (K - 1))}, where \code{s} is the number of regression coefficients, \code{K} is the assumed number of latent classes.
#'
#' \item Initial values for the variances of the independent random effects in the latent class membership model at the stratum level. The list entry is of the form, for example, \code{rep(0.2, K - 1)}. If stratum level random effects are not desired, then use \code{list(NULL)}.
#'
#' \item Initial values for the spatial variances. The specification is the same as in 2. If spatial random effects are not desired, then use \code{list(NULL)}.
#'
#' \item Initial values for the variances of the independent random effects in the latent class membership model at the cluster level. The specification is the same as in 2. If independent random effects are not desired, then use \code{list(NULL)}.
#'
#' \item If \code{spline = TRUE}, initial values for the regression coefficients for the B-splines in the latent class membership model. The form is the same as 1. For example, enter \code{matrix(0, nrow = r, ncol = (K - 1))}. If \code{spline = FALSE}, then use \code{list(NULL)}.
#'
#' \item Initial values for the regression coefficients in the longitudinal outcomes model. The object is a \code{p} by \code{K} matrix, where \code{p} is the number of regression coefficients and \code{K} is the number of latent classes. An example is \code{matrix(rnorm(p*K), nrow = p, ncol = K)}.
#'
#' \item Initial values for the variances of the independent random effects in the longitudinal outcomes model at the stratum level. The object is a \code{K} length vector. An example is \code{rep(0.1, K)}. If stratum level random effects are not desired, then use \code{list(NULL)}.
#'
#' \item Initial values for the variances of the independent random effects in the longitudinal outcomes model at the cluster level. The form is the same as in 7. If cluster level random effects are not desired, then use \code{list(NULL)}.
#'
#' \item Initial values for the variance-covariance matrix of the subject level random effects in the longitudinal outcomes model. The object is a \code{q} by \code{q} by \code{K} array, where \code{q} is the number of random effects and \code{K} is the number of assumed latent class. An example is \code{array(diag(0.5, q), dim = c(q, q, K))}.
#'
#' \item Initial values for the variances of the observation-level errors in the longitudinal outcomes model. The form is the same as in 7.}
#'
#' An example of the completed list object for specification of the initial values is:
#' \code{inits <- list(matrix(rep(0, s * (K - 1)), nrow = s, ncol
#' = (K - 1)), rep(0.2, K - 1), rep(0.2, K - 1), rep(0.2, K - 1), matrix(0, nrow = r, ncol = (K - 1)), matrix(rnorm(p * K), nrow = p, ncol = K), rep(0.1, K), rep(0.1, K), array(diag(0.5, q), dim = c(q, q, K)), rep(1, K))}
#' }
#' }
#'
#' @export
#' @return Model summaries are printed to the console, including posterior means and 95\% credible intervals, posterior latent class assignment, and model comparison statistics. A label switching diagnostic using Stephen's method from the \code{label.switching} package in \code{R} is printed.
#'
#' A list of matrices of stored MCMC samples (\code{n.samples-burn}) for unknown parameters, including: \enumerate{
#' \item \code{store_delta}: Regression coefficients from the latent class membership model.
#' \item \code{store_gamma2}: Variance of the stratum-level random effects in the latent class membership model.
#' \item \code{store_xi2}: Variance of the cluster-level spatial random effects in the latent class membership model.
#' \item \code{store_tau2}: Variance of the cluster-level independent random effects in the latent class membership model.
#' \item \code{store_alpha}: Regression coefficients for the B-spline in the latent class membership model.
#' \item \code{store_beta}: Regression coefficients for the longitudinal outcomes model.
#' \item \code{store_psi2}: Variance of the stratum-level random effects in the longitudinal outcomes model.
#' \item \code{store_omega2}: Variance of the cluster-level random effects in the longitudinal outcomes model.
#' \item \code{store_phi}: Variance-covariance of the subject-level random effects in the longitudinal outcomes model.
#' \item \code{store_sigma2}: Observation-level variance in the longitudinal outcomes model.
#' \item \code{store_pi}: Posterior probabilities of latent class assignment.
#' \item \code{store_priorPi}: Latent class membership probabilities.
#' }
#'
#' If \code{writeSamples = TRUE}, then samples for the random effects, draws from the posterior predictive distribution, and discrepancy measure are written to separate text files after burn-in. Corresponding text files with the column names are also written.
#'
#' For the latent class membership model, the stratum-level random effects are stored in \code{store_lambda.txt}; the spatial cluster-level random effects are stored in \code{store_nu.txt}; and the independent cluster-level random effects the subject-level random effects are stored in \code{store_u.txt}.
#'
#' For the longitudinal outcomes model, the stratum-level random effects are stored in \code{store_zeta.txt}; the cluster-level random effects are stored in \code{store_rho.txt}; the subject-level random effects are stored in \code{store_b.txt}.
#'
#'Draws of the longitudinal outcome conditional on latent class from the posterior predictive distribution are stored in \code{store_Ydraw.txt}. This file is of dimension \code{n.samples - burn} by \code{N}, where \code{N} is the number of observations (i.e., \code{length(Y)}).
#'
#'Draws of the longitudinal outcome from the posterior predictive distribution based on additionally re-drawing the latent class membership indicators are stored in \code{store_redraw.txt} with the iteration number appended. This file is of dimension \code{n.samples - n.reps} by \code{N}, where \code{N} is the number of observations (i.e., \code{length(Y)}). This file will also contain draws of residuals for the complete and replicated data, in addition to latent class membership indicator at each iteration. SubjectID, clusterID, and stratum ID will also be appended. See \code{?get_post_pred_redraw()}.
#'
#'Samples of a measure of discrepancy are stored in \code{store_T.txt}. This file is of dimension \code{n.samples - burn} by \code{2}, where the first column is the discrepancy measure computed using the observed data, and the second column is the discrepancy measure using the replicated data. See \code{?get_discrepancy_plot}. The discrepancy measure is computed by re-drawing the latent class indicators, in addition to the longitudinal outcomes.
#'
#'@examples
#'
#'#-------- Load data from Bsvygmm
#' #Data
#' data(data)
#' #Adjacency matrix
#' data(ADJ)
#' LCStratumModelType <- "Random"
#' LCClusterModelType <- "Both"
#' LRModelType <- "Both"
#'
#' #Number of latent classes
#' K <- 2
#'
#' #Outcome
#' Y <- "Y"
#' time <- "time"
#' subjectID <- "subjectID"
#' clusterID <- "clusterID"
#' stratumID <- "stratumID"
#'
#' #Latent class membership model, design matrix for fixed effects
#' form1 <- ~ x1
#' s <- length(all.vars(form1)) + 1
#' #Longitudinal outcomes model, design matrix for fixed effects; equal to "dummy" or some formula
#' form2 <- "dummy"
#' #Longitudinal outcomes model, design matrix for random effect; equal to "same" or some formula
#' form3 <- "same"
#' #Store formulas in named list
#' forms <- list(W = form1, Xf = form2, Xr = form3)
#'
#' #Number of time points and random effects is 3, indicate for prior and init specification below
#' q <- p <- 3
#'
#' spline <- FALSE
#' B <- NULL
#'
#' #Prior for hierarchical and observation level variance
#' hierVar <- list("IG", "IG")
#'
#' #Priors
#' priors <- list(list(rep(0, s), diag(1, s)), list(.1, .1), list(2, 1),
#'      list(.1, .1), list(NULL), list(rep(0, p), diag(10, p)),
#'      list(.1, .1), list(.1, .1), list((q + 2), diag(0.25, q)), list(.1, .1))
#'
#' #Initial values
#' inits <- list(matrix(rep(0, s * (K - 1)), nrow = s, ncol = (K - 1)), rep(0.2, K - 1),
#'      rep(0.2, K - 1), rep(0.2, K - 1), NULL, matrix(rnorm(p * K), nrow = p, ncol = K),
#'      rep(0.1, K), rep(0.1, K), array(diag(0.5, q),
#'      dim = c(q, q, K)), rep(1, K))
#'
#' #MCMC algorithm
#' writeSamples <- TRUE
#' n.samples <- 2000
#' burn <- 1000
#' update <- 10
#' n.reps <- 500
#' monitor <- TRUE
#' res <- Bsvygmm(K = K, data = data, forms = forms, ADJ = ADJ, Y = Y, time = time, subjectID = subjectID, clusterID = clusterID, stratumID = stratumID, spline = spline, B = B, LCStratumModelType = LCStratumModelType, LCClusterModelType = LCClusterModelType, LRModelType = LRModelType, priors = priors, hierVar = hierVar, inits = inits, n.samples = n.samples, burn = burn, monitor = monitor, update = update, n.reps = n.reps, writeSamples = writeSamples)
#'
#' #Get Bayesian posterior predictive p-value, see ?get_discrepancy_plot
#' #Set working directory to where discrepancy samples are written
#' store_T <- read.table("store_T.txt", header = FALSE, sep = ",")
#' get_discrepancy_plot(store_T)
Bsvygmm <- function(K, data, forms, ADJ, Y, time, subjectID, clusterID, stratumID,  spline, B, LCStratumModelType, LCClusterModelType, LRModelType, priors, hierVar, inits, n.samples, burn, monitor, update, n.reps, writeSamples){

  #----------------------- Background
  data$subjectID <- data[ , subjectID]
  data$clusterID <- data[ , clusterID]
  data$stratumID <- data[ , stratumID]
  data$time <- data[ , time]

  if (forms[["Xf"]] == "dummy" & !is.factor(data$time)) {
    data$time <- factor(data$time)
  }

  ### Identifiers
  subjectID <- data$subjectID
  clusterIDObs <- data$clusterID
  stratumIDObs <- data$stratumID

  ### Longitudinal outcomes model
  Y <- data[ , Y]

  ### Design matrices
  # Latent class membership model design matrix
  var_namess <- all.vars(forms[["W"]])
  datas <- aggregate(data[ , c("subjectID", "clusterID", "stratumID", var_namess)], by = list(data$subjectID), FUN = tail, n = 1)
  W <- model.matrix(forms[["W"]], data = datas)
  s <- ncol(W)
  n <- nrow(W)

  clusterIDSub <- datas$clusterID
  stratumIDSub <- datas$stratumID

  # Design matrix for fixed effects in longitudinal outcomes model
  if (forms[["Xf"]] == "dummy") {
    Vf <- data.matrix(dummy::dummy(data.frame(data[ , "time"]), int = TRUE))
    p <- ncol(Vf)
    colnames(Vf) <- paste0("time", 1:p)
  } else {
    Vf <- model.matrix(forms[["Xf"]], data = data)
    p <- ncol(Vf)
  }

  # Design matrix for fixed effects in longitudinal outcomes model
  if (forms[["Xr"]] == "same") {
    Vr <- Vf
    q <- ncol(Vr)
    colnames(Vr) <- paste0("time", 1:q)
  } else {
    Vr <- model.matrix(forms[["Xr"]], data = data)
    q <- ncol(Vr)
  }


  # Number of unique clusters in the study sample
  J <- length(unique(clusterIDSub))
  # Number of clusters in the study region, regardless of whether observations appear in them
  if (LCClusterModelType %in% c("Both", "Str")) {
    nJ <- nrow(ADJ)
  }

  # Number of unique strata in the longitudinal data
  M <- length(unique(stratumIDSub))


  if (spline) {
    #! Elaborate with spline information
    nB <- ncol(B)
  }

  ### Prior distributions
  # Multinomial probit or binary probit model for latent class assignment
  prior_mu_delta <- priors[[1]][[1]]
  prior_Sigma_delta <- priors[[1]][[2]]

  # Stratum level random effect in the latent class model and observation level variance in longitudinal class model
  # hierVar[[1]] is for the hierarchical variance of the random effects
  # hierVar[[2]] is for the hierarchical variance of the observation level data y

  if (hierVar[[1]] == "IG" & LCStratumModelType == "Random") {
     prior_shape_gamma2 <- priors[[2]][[1]]
     prior_scale_gamma2 <- priors[[2]][[2]]
  } else {
     prior_shape_gamma2 <- NULL
     prior_scale_gamma2 <- NULL
  }


  if (hierVar[[1]] == "IG" & LCClusterModelType %in% c("Both", "Str")) {
    prior_shape_xi2 <- priors[[3]][[1]]
    prior_scale_xi2 <- priors[[3]][[2]]
  } else {
    prior_shape_xi2 <- NULL
    prior_scale_xi2 <- NULL
  }


  if (hierVar[[1]] == "IG" & LCClusterModelType %in% c("Both", "Unstr")) {
    prior_shape_tau2 <- priors[[4]][[1]]
    prior_scale_tau2 <- priors[[4]][[2]]
  } else {
    prior_shape_tau2 <- NULL
    prior_scale_tau2 <- NULL
  }

  # For spline
  if (spline) {
    prior_Sigma_alpha <- priors[[5]][[1]]
  }

  # Longitudinal outcomes
  #Prior for coefficients that vary by latent class
  prior_mu_beta <- priors[[6]][[1]]
  prior_Sigma_beta <- priors[[6]][[2]]

  #Prior for variance of stratum level random effects
  if (hierVar[[1]] == "IG" & LRModelType %in% c("Both", "Stratum")) {
    prior_shape_psi2 <- priors[[7]][[1]]
    prior_scale_psi2 <- priors[[7]][[2]]
  } else {
    prior_shape_psi2 <- NULL
    prior_scale_psi2 <- NULL
  }

  #Prior for variance of cluster level random effects
  if (hierVar[[1]] == "IG" & LRModelType %in% c("Both", "Cluster")) {
    prior_shape_omega2 <- priors[[8]][[1]]
    prior_scale_omega2 <- priors[[8]][[2]]
  } else {
    prior_shape_omega2 <- NULL
    prior_scale_omega2 <- NULL
  }


  #Prior for covariance for subject level, common prior across latent classes even if class specific covariance is desired
  prior_df_phi <- priors[[9]][[1]]
  prior_scale_phi <- priors[[9]][[2]]

  if (hierVar[[2]] == "IG") {
    prior_shape_sigma2 <- priors[[10]][[1]]
    prior_scale_sigma2 <- priors[[10]][[2]]
  } else {
    prior_shape_sigma2 <- NULL
    prior_scale_sigma2 <- NULL
  }


  ### Storage of MCMC samples
  # Latent class membership model storage
  store_delta <- matrix(NA, nrow = (n.samples - burn), ncol = ((K - 1)*s), dimnames = list(paste("Iteration", (burn + 1):n.samples, sep = "_"), paste(paste("Class", rep(2:K, each = length(colnames(W))), sep = ""), rep(colnames(W), times = (K - 1)), sep = "_")))

  if (spline) {
    store_alpha <- matrix(NA, nrow = (n.samples - burn), ncol = ((K - 1) * nB), dimnames = list(paste("Iteration", (burn + 1):n.samples, sep = "_"), paste(paste("Class", rep(2:K, each = nB), sep = ""), rep(paste("alpha", 1:nB, sep = ""), times = (K - 1)), sep = "_")))
    } else {
      store_alpha <- NULL
  }

  if (LCStratumModelType == "Random") {
    store_gamma2 <- matrix(NA, nrow = (n.samples - burn), ncol = (K - 1), dimnames = list(paste("Iteration", (burn + 1):n.samples, sep = "_"), paste(paste("Class", 2:K, sep = ""), "gamma2", sep = "_")))
  } else {
    store_gamma2 <- NULL
  }

  if (LCClusterModelType %in% c("Both", "Unstr")) {
    store_tau2 <- matrix(NA, nrow = (n.samples - burn), ncol = (K - 1), dimnames = list(paste("Iteration", (burn + 1):n.samples, sep = "_"), paste(paste("Class", 2:K, sep = ""), "tau2", sep = "_")))
  } else {
    store_tau2 <- NULL
  }

  if (LCClusterModelType %in% c("Both", "Str")) {
    store_nu <- matrix(NA, nrow = (n.samples - burn), ncol = nJ*(K - 1), dimnames = list(paste("Iteration", (burn + 1):n.samples, sep = "_"), paste(paste("Class", rep(2:K, each = nJ), sep = ""), paste("nu", rep(1:nJ, times = (K - 1)), sep = ""), sep = "_")))
    store_xi2 <- matrix(NA, nrow = (n.samples - burn), ncol = (K - 1), dimnames = list(paste("Iteration", (burn + 1):n.samples, sep = "_"), paste(paste("Class", 2:K, sep = ""), "xi2", sep = "_")))
  } else {
    store_nu <- NULL
    store_xi2 <- NULL
  }

  store_pi <- matrix(NA, nrow = (n.samples - burn), ncol = (K*n), dimnames = list(paste("Iteration", (burn + 1):n.samples, sep = "_"), paste(paste0("Class", rep(1:K, each = n)), "Subject", 1:n, sep = "_")))
  store_priorPi <- matrix(NA, nrow = (n.samples - burn), ncol = (K*n), dimnames = list(paste("Iteration", (burn + 1):n.samples, sep = "_"), paste(paste0("Class", rep(1:K, each = n)), "Subject", 1:n, sep = "_")))

  store_beta <- matrix(NA, nrow = (n.samples - burn), ncol = (p*K), dimnames = list(paste("Iteration", (burn + 1):n.samples, sep = "_"), paste(paste("Class", rep(1:K, each = p), sep = ""), rep(colnames(Vf), K), sep = "_")))

  store_sigma2 <- matrix(NA, nrow = (n.samples - burn), ncol = K, dimnames = list(paste("Iteration", (burn + 1):n.samples, sep = "_"), paste(paste("Class", 1:K, sep = ""), "sigma2", sep = "_")))

  store_phi <- matrix(NA, nrow = (n.samples - burn), ncol = q*q*K, dimnames = list(paste("Iteration", (burn + 1):n.samples, sep = "_"),  paste(paste("Class", rep(1:K, each = q*q), sep = ""), rep(paste("phi", rep(1:q, times = q), rep(1:q, each = q), sep = "") , K), sep = "_")))

  if (LRModelType %in% c("Both", "Cluster")) {
    store_omega2 <- matrix(NA, nrow = (n.samples - burn), ncol = K, dimnames = list(paste("Iteration", (burn + 1):n.samples, sep = "_"),  paste(paste("Class", rep(1:K, each = 1), sep = ""), rep("omega2", K), sep = "_")))
  } else {
    store_omega2 <- NULL
  }

  if (LRModelType %in% c("Both", "Stratum")) {
    store_psi2 <- matrix(NA, nrow = (n.samples - burn), ncol = K, dimnames = list(paste("Iteration", (burn + 1):n.samples, sep = "_"),  paste(paste("Class", rep(1:K, each = 1), sep = ""), rep("psi2", K), sep = "_")))
  } else {
    store_psi2 <- NULL
  }

  # Model comparison pieces storage
  #Entropy
  store_entropy <- rep(NA, (n.samples - burn - 1))
  #Observed data likelihood for each subject
  store_observed_llik_sub <- matrix(NA, nrow = (n.samples - burn - 1), ncol = n)
  #Observed data likelihood
  store_observed_llik <- rep(NA, (n.samples - burn - 1))
  #Observed data posterior
  store_observed_posterior <- rep(NA, (n.samples - burn - 1))
  #Effective sample size
  store_ESS <- rep(NA, (n.samples - burn - 1))

  ### Initialization
  # Multinomial probit or binary probit model
  C <- sample(1:K, n, prob = rep(1 / K, K), replace = TRUE)
  delta <- inits[[1]]

  #Stratum level unstructured heterogeneity
  if (LCStratumModelType == "Random") {
    gamma2 <- inits[[2]]
    lambda <- matrix(0, nrow = M, ncol = K - 1)
  } else {
    gamma2 <- NULL
    lambda <- NULL
  }

  #Specifying cluster level heterogeneity type
  #Spatial random effects or not
  if (LCClusterModelType %in% c("Both", "Str")) {
    xi2 <- inits[[3]]
    nu <- matrix(0, nrow = nJ, ncol = K - 1)
    clusterIDSubf <- factor(clusterIDSub, levels = seq(1, nJ, by = 1))
  } else {
    xi2 <- NULL
    nu <- NULL
  }

  #Unstructured heterogeneity or not
  if (LCClusterModelType %in% c("Both", "Unstr")) {
    u <- matrix(0, nrow = J, ncol = K - 1)
    tau2 <- inits[[4]]
  } else {
    u <- NULL
    tau2 <- NULL
  }

  #Spline or not
  if (spline) {
    alpha <- inits[[5]]
  } else {
    alpha <- NULL
    B <- NULL
  }

  # Initialization of Z
  if (K > 2) {
    # K-1 by K-1 covariance matrix
    VarCovP <- matrix(1, nrow = (K - 1), ncol = (K - 1))
    diag(VarCovP) <- 2
    err <- do.call("rbind", sapply(1:n, function(x) {res <- mnormt::rmnorm(1, mean = rep(0, (K-1)), varcov = VarCovP)
  return(res) }, simplify = FALSE))
  } else if (K == 2) {
    err <- mnormt::rmnorm(n, varcov = diag(K - 1))
  }



  #All model types assuming stratum level heterogeneity
  #Initialization for stratum level heterogeniety and spatially structured heterogeneity at the cluster level only
  if (LCClusterModelType == "Str" & LCStratumModelType == "Random") {

    if (K > 2) {
      if (spline) {
        Z <- err + W %*% delta + B %*% alpha + lambda[factor(stratumIDSub), ] + nu[clusterIDSubf, ]
      } else {
        Z <- err + W %*% delta + lambda[factor(stratumIDSub), ] + nu[clusterIDSubf, ]
      }
    } else if (K == 2) {
      if (spline) {
        Z <- err + W %*% delta + B %*% alpha + lambda[factor(stratumIDSub)] + nu[clusterIDSubf]
      } else {
        Z <- err + W %*% delta + lambda[factor(stratumIDSub)] + nu[clusterIDSubf]
      }
    }

  }


  #Initialization for stratum level heterogeniety and unstructured heterogeneity at the cluster level only
  if (LCClusterModelType == "Unstr"  & LCStratumModelType == "Random") {

    if (K > 2) {
      if (spline) {
        Z <- err + W %*% delta + B %*% alpha + u[factor(clusterIDSub), ] + lambda[factor(stratumIDSub), ]
      } else {
        Z <- err + W %*% delta + u[factor(clusterIDSub), ] + lambda[factor(stratumIDSub), ]
      }
    } else if (K == 2) {
      if (spline) {
        Z <- err + W %*% delta + B %*% alpha + u[factor(clusterIDSub)] + lambda[factor(stratumIDSub)]
      } else {
        Z <- err + W %*% delta + u[factor(clusterIDSub)] + lambda[factor(stratumIDSub)]
      }
    }

  }


  #Initialization for stratum level heterogeniety and spatially structured and unstructured heterogeneity at the cluster level only
  if (LCClusterModelType == "Both" & LCStratumModelType == "Random") {
    if (K > 2) {
      if (spline) {
        Z <- err + W %*% delta + B %*% alpha + nu[clusterIDSubf, ] +  u[factor(clusterIDSub), ] + lambda[factor(stratumIDSub), ]
      } else {
        Z <- err + W %*% delta + nu[clusterIDSubf, ] + u[factor(clusterIDSub), ] + lambda[factor(stratumIDSub), ]
      }
    } else if (K == 2) {
      if (spline) {
        Z <- err + W %*% delta + B %*% alpha + nu[clusterIDSubf] + u[factor(clusterIDSub)] + lambda[factor(stratumIDSub)]
      } else {
        Z <- err + W %*% delta + nu[clusterIDSubf] + u[factor(clusterIDSub)] + lambda[factor(stratumIDSub)]
      }
    }

  }

  #Initialization for stratum level heterogeniety and no spatially structured and unstructured heterogeneity at the cluster level
  if (LCClusterModelType == "None" & LCStratumModelType == "Random") {
    if (K > 2) {
      if (spline) {
        Z <- err + W %*% delta + B %*% alpha + lambda[factor(stratumIDSub), ]
      } else {
        Z <- err + W %*% delta + lambda[factor(stratumIDSub), ]
      }
    } else if (K == 2) {
      if (spline) {
        Z <- err + W %*% delta + B %*% alpha + lambda[factor(stratumIDSub)]
      } else {
        Z <- err + W %*% delta + lambda[factor(stratumIDSub)]
      }
    }

  }




  #All model types assuming no stratum level heterogeneity
  #Initialization for spatially structured heterogeneity at the cluster level only
  if (LCClusterModelType == "Str" & LCStratumModelType == "Fixed") {

    if (K > 2) {
      if (spline) {
        Z <- err + W %*% delta + B %*% alpha + nu[clusterIDSubf, ]
      } else {
        Z <- err + W %*% delta + nu[clusterIDSubf, ]
      }
    } else if (K == 2) {
      if (spline) {
        Z <- err + W %*% delta + B %*% alpha + nu[clusterIDSubf]
      } else {
        Z <- err + W %*% delta + nu[clusterIDSubf]
      }
    }

  }


  #Initialization for unstructured heterogeneity at the cluster level only
  if (LCClusterModelType == "Unstr"  & LCStratumModelType == "Fixed") {

    if (K > 2) {
      if (spline) {
        Z <- err + W %*% delta + B %*% alpha + u[factor(clusterIDSub), ]
      } else {
        Z <- err + W %*% delta + u[factor(clusterIDSub), ]
      }
    } else if (K == 2) {
      if (spline) {
        Z <- err + W %*% delta + B %*% alpha + u[factor(clusterIDSub)]
      } else {
        Z <- err + W %*% delta + u[factor(clusterIDSub)]
      }
    }

  }


  #Initialization for spatially structured and unstructured heterogeneity at the cluster level only
  if (LCClusterModelType == "Both" & LCStratumModelType == "Fixed") {
    if (K > 2) {
      if (spline) {
        Z <- err + W %*% delta + B %*% alpha + nu[clusterIDSubf, ] +  u[factor(clusterIDSub), ]
      } else {
        Z <- err + W %*% delta + nu[clusterIDSubf, ] + u[factor(clusterIDSub), ]
      }
    } else if (K == 2) {
      if (spline) {
        Z <- err + W %*% delta + B %*% alpha + nu[clusterIDSubf] + u[factor(clusterIDSub)]
      } else {
        Z <- err + W %*% delta + nu[clusterIDSubf] + u[factor(clusterIDSub)]
      }
    }

  }


  #Initialization for no spatially structured and unstructured heterogeneity at the cluster level
  if (LCClusterModelType == "None" & LCStratumModelType == "Fixed") {
    if (K > 2) {
      if (spline) {
        Z <-  err + W %*% delta + B %*% alpha
      } else {
        Z <- err + W %*% delta
      }

    } else if (K == 2) {
      if (spline) {
        Z <- err + W %*% delta + B %*% alpha
        } else {
        Z <- err + W %*% delta
      }
    }

  }


  # Longitudinal outcomes
  beta <- inits[[6]]

  if (LRModelType %in% c("Both", "Stratum")) {

    psi2 <- inits[[7]]
    zeta <- matrix(NA, nrow = M, ncol = K)
    C_expand <- C[factor(subjectID)]

    for (k in 1:K) {
      #Subject level information
      ind_obs <- which(C_expand == k)
      stratumIDObsk <- stratumIDObs[ind_obs]

      #Stratum level information
      ind_stratum <- match(sort(unique(stratumIDObsk)), sort(unique(stratumIDObs)))
      zeta[ind_stratum, k] <- rnorm(length(ind_stratum), mean = 0, sd = sqrt(psi2[k]))
    }

  } else {
    psi2 <- NULL
    zeta <- NULL
  }


  if (LRModelType %in% c("Both", "Cluster")) {

    omega2 <- inits[[8]]
    rho <- matrix(NA, nrow = J, ncol = K)
    C_expand <- C[factor(subjectID)]

    for (k in 1:K) {
      #Subject level information
      ind_obs <- which(C_expand == k)
      clusterIDObsk <- clusterIDObs[ind_obs]

      #Cluster level information
      ind_cluster <- match(sort(unique(clusterIDObsk)), sort(unique(clusterIDObs)))

      rho[ind_cluster, k] <- rnorm(length(ind_cluster), mean = 0, sd = sqrt(omega2[k]))

    }

  } else {

    omega2 <- NULL
    rho <- NULL

  }

  phi <- inits[[9]]
  b <- matrix(NA, nrow = n, ncol = q)

  for (k in 1:K) {

    ind_sub <- which(C == k)

    b[ind_sub, ] <- t(sapply(1:length(ind_sub), function(x){
        mnormt::rmnorm(1, mean = rep(0, q), varcov = phi[ , , k])}))
  }

  sigma2 <- inits[[10]]

  ### Start time for loop
  ptm <- proc.time()

  for (i in 1:n.samples) {

    # Update parameters in latent class membership model
    #Latent variable Z and regression coefficients delta depending on K > 2 or not
    if (K > 2) {
      Z <- update_Z_MNP(Z = Z, C = C, delta = delta, alpha = alpha, deltaStratum = lambda, u = u, nu = nu, W = W, B = B, clusterIDSub = clusterIDSub, stratumIDSub = stratumIDSub, spline = spline)
      priorPik <- update_priorPik(Z = Z)
      delta <- update_delta_MNP(Z = Z, delta = delta, alpha = alpha, deltaStratum = lambda, u = u, nu = nu, W = W, B = B, prior.mu = prior_mu_delta, prior.Sigma = prior_Sigma_delta, clusterIDSub = clusterIDSub, stratumIDSub = stratumIDSub, spline = spline)
    } else {
      Z <- update_Z_BinaryP(Z = Z, C = C, delta = delta, alpha = alpha, deltaStratum = lambda, u = u, nu = nu, W = W, B = B, clusterIDSub = clusterIDSub, stratumIDSub = stratumIDSub, spline = spline)
      priorPik <- update_priorPik(Z = Z)
      delta <- update_delta_BinaryP(Z = Z, delta = delta, alpha = alpha, deltaStratum = lambda, u = u, nu = nu, W = W, B = B, prior.mu = prior_mu_delta, prior.Sigma = prior_Sigma_delta, clusterIDSub = clusterIDSub, stratumIDSub = stratumIDSub, spline = spline)
    }

    # Stratum level heterogeniety
    if (LCStratumModelType == "Random") {
      if (K > 2) {
        lambda <- update_deltaStratum_MNP(Z = Z, delta = delta, alpha = alpha, u = u, nu = nu, gamma2 = gamma2, W = W, B = B, clusterIDSub = clusterIDSub, stratumIDSub = stratumIDSub, spline = spline)
      } else {
        lambda <- update_deltaStratum_BinaryP(Z = Z, delta = delta, alpha = alpha, u = u, nu = nu, gamma2 = gamma2, W = W, B = B, clusterIDSub = clusterIDSub, stratumIDSub = stratumIDSub, spline = spline)
      }
      gamma2 <- update_gamma2(deltaStratum = lambda, prior.shape = prior_shape_gamma2, prior.scale = prior_scale_gamma2, hierVar = hierVar[[1]])
    }


    #Spline coefficients
    if (spline) {
      if (K > 2) {
        alpha <- update_alpha_BSpline_MNP(Z = Z, delta = delta, deltaStratum = lambda,  u = u, nu = nu, W = W, B = B, prior.Sigma = prior_Sigma_alpha, clusterIDSub = clusterIDSub, stratumIDSub = stratumIDSub)
      } else {
        alpha <- update_alpha_BSpline_BinaryP(Z = Z, delta = delta, deltaStratum = lambda, u = u, nu = nu, W = W, B = B, prior.Sigma = prior_Sigma_alpha, clusterIDSub = clusterIDSub, stratumIDSub = stratumIDSub)
      }
    }


    #Unstructure random intercepts and variance terms for cluster and stratum levels, hierarchically centered
    if (LCClusterModelType %in% c("Both", "Unstr")) {
      if (K > 2) {
        u <- update_u_MNP(Z = Z, delta = delta, alpha = alpha, deltaStratum = lambda, nu = nu, tau2 = tau2, W = W, B = B, clusterIDSub = clusterIDSub, stratumIDSub = stratumIDSub, spline = spline)
      } else {
        u <- update_u_BinaryP(Z = Z, delta = delta, alpha = alpha, deltaStratum = lambda, nu = nu, tau2 = tau2, W = W, B = B, clusterIDSub = clusterIDSub, stratumIDSub = stratumIDSub, spline = spline)
      }
      tau2 <- update_tau2(u = u, prior.shape = prior_shape_tau2, prior.scale = prior_scale_tau2, hierVar = hierVar[[1]])
    }


    #Spatial random intercept and variance term
    if (LCClusterModelType %in% c("Both", "Str")) {
      nu <- update_nu(Z = Z, delta = delta, alpha = alpha, deltaStratum = lambda, u = u, nu = nu, xi2 = xi2, W = W, B = B, ADJ = ADJ, clusterIDSub = clusterIDSub, stratumIDSub = stratumIDSub, spline = spline)
      xi2 <- update_xi2(nu = nu, ADJ = ADJ, prior.shape = prior_shape_xi2, prior.scale = prior_scale_xi2, hierVar = hierVar[[1]])
    }

    # Update longitudinal outcomes parameters
    #Stratum level random effects
    if (LRModelType %in% c("Both", "Stratum")) {
      zeta <- update_zeta(C = C, b = b, rho = rho, beta = beta, sigma2 = sigma2, Psi = psi2, Y = Y, Vf = Vf, Vr = Vr, subjectID = subjectID, clusterIDObs = clusterIDObs, stratumIDObs = stratumIDObs)
      psi2 <- update_Psi(C = C, zeta = zeta, prior.scale = prior_scale_psi2, prior.shape = prior_shape_psi2, hierVar = hierVar[[1]])
    }

    #Cluster level random effects
    if (LRModelType %in% c("Both", "Cluster")) {
      rho <- update_rho(C = C, b = b, zeta = zeta, beta = beta, sigma2 = sigma2, Omega = omega2, Y = Y, Vf = Vf, Vr = Vr, subjectID = subjectID, clusterIDObs = clusterIDObs, stratumIDObs = stratumIDObs)
      omega2 <- update_Omega(C = C, rho = rho, prior.scale = prior_scale_omega2, prior.shape = prior_shape_omega2, hierVar = hierVar[[1]])
    }

    #Subject level random effects
    b <- update_b(C = C, rho = rho, zeta = zeta, beta = beta, Phi = phi, sigma2 = sigma2, Y = Y, Vf = Vf, Vr = Vr, subjectID = subjectID, clusterIDObs = clusterIDObs, stratumIDObs = stratumIDObs)
    #Class specific regression coefficients (at the minimum intercept terms)
    beta <- update_betak_Marginal(C = C, sigma2 = sigma2, Phi = phi, Omega = omega2, Psi = psi2, prior.mu = prior_mu_beta, prior.Sigma = prior_Sigma_beta, Y = Y, Vf = Vf, Vr = Vr, subjectID = subjectID)

    #Subject level random effect covariance matrix
    phi <- update_Phi(C = C, q = q, b = b, prior.df = prior_df_phi, prior.scale = prior_scale_phi)

    #Level 1 parameters
    sigma2 <- update_sigma2(C = C, zeta = zeta, rho = rho, b = b, beta = beta, Y = Y, Vr = Vr, Vf = Vf, prior.shape = prior_shape_sigma2, prior.scale = prior_scale_sigma2, subjectID = subjectID, clusterIDObs = clusterIDObs, stratumIDObs = stratumIDObs, hierVar = hierVar[[2]])

    # Model comparison, current iteration i of above parameters consistent with iteration i - l of latent class indicators, placing this step here ensures this
    # Use burn + 1 so we start to save model comparison statistics after sample of latent class indicators is > burn.
    # Therefore storage will have i-burn-1 rows
    if (i > (burn + 1)) {

      #Entropy
      if (sum(pik < 1E-323) > 0) {
        pik[which(pik < 1E-323)] <- rep(1E-323, sum(pik < 1E-323))
      }
      # Entropy, used in DIC4 and ICL-BIC
      store_entropy[(i - burn - 1)] <- (-1) * sum(apply(pik * log(pik), 1, sum))

      # Observed data likelihood for each subject f(y_i) = \sum_k \pi_{ik} f(y_i | \beta_k)
      #Used in DIC3, DIC4, BIC and LPML
      #Subject specific, not on log scale
      store_observed_llik_sub[(i - burn - 1), ] <- get_observed_llik_sub(priorPik = priorPik, beta = beta, sigma2 = sigma2, Phi = phi, Omega = omega2, Psi = psi2, Y = Y, Vf = Vf, Vr = Vr, subjectID = subjectID)
      # Compute log observed data likelihood contribution as log( \prod_{i = 1}^n f(y_i) ) = \sum_{i=1}^n log f(y_i)
      store_observed_llik[(i - burn - 1)] <- sum(log(store_observed_llik_sub[(i - burn - 1), ]))

      # Observed data posterior, used in DIC4
      store_observed_posterior[(i - burn - 1)] <-  get_observed_posterior(llik = store_observed_llik[(i - burn - 1)], beta =  beta, Phi = phi, sigma2 = sigma2, Omega = omega2, Psi = psi2, delta = delta, alpha = alpha, u = u, deltaStratum = lambda, nu = nu, tau2 = tau2, gamma2 = gamma2, xi2 = xi2, priorPik = priorPik, ADJ = ADJ, priors = priors, hierVar = hierVar, LCStratumModelType = LCStratumModelType, LCClusterModelType = LCClusterModelType, LRModelType = LRModelType, spline = spline)

      # Effective sample size, used in BIC and ICL-BIC
      store_ESS[(i - burn - 1)] <- get_ESS(sigma2 = sigma2, Phi = phi, Omega = omega2, Psi = psi2, priorPik = priorPik, Y = Y, Vr = Vr, subjectID = subjectID)

    }

    if (i > burn) {

      # Draw Y from the posterior predictive distribution conditional on C
      Ydraw <- get_postPred(K = K, C = C, b = b, rho = rho, zeta = zeta, beta = beta, sigma2 = sigma2, Vr = Vr, Vf = Vf, subjectID = subjectID, clusterIDObs = clusterIDObs, stratumIDObs = stratumIDObs)

      # Discrepancy measure at observed and replicated discrepancy measure, this will be written to file; This is based on re-drawing latent class.
      # Observed measure is in the first column, replicated measure is in the second column.
      # Placed here separately with i > burn because outputted below with other parameters
       store_T <- get_discrepancy(K = K, C = C, b = b, sigma2 = sigma2, beta = beta, zeta = zeta, rho = rho, Phi = phi, Omega = omega2, Psi = psi2, Z = Z, delta = delta, alpha = alpha, deltaStratum = lambda, u = u, nu = nu, gamma2 = gamma2, tau2 = tau2, xi2 = xi2, priorPik = priorPik,  W = W, B = B, ADJ = ADJ, clusterIDSub = clusterIDSub, stratumIDSub = stratumIDSub, spline = spline,  subjectID = subjectID, clusterIDObs = clusterIDObs, stratumIDObs = stratumIDObs, Y = Y, Vf = Vf, Vr = Vr, LRModelType = LRModelType, LCClusterModelType = LCClusterModelType, LCStratumModelType = LCStratumModelType)


      if (i > (n.samples - n.reps)) {

        # Based on redrawing latent class; includes residuals and redraws of Y. Placed here to save only a selected number of simulations from n.reps

        store_redraw <- get_post_pred_redraw(K = K, C = C, b = b, sigma2 = sigma2, beta = beta, zeta = zeta, rho = rho, Phi = phi, Omega = omega2, Psi = psi2, Z = Z, delta = delta, alpha = alpha, deltaStratum = lambda, u = u, nu = nu, gamma2 = gamma2, tau2 = tau2, xi2 = xi2, priorPik = priorPik,  W = W, B = B, ADJ = ADJ, clusterIDSub = clusterIDSub, stratumIDSub = stratumIDSub, spline = spline,  subjectID = subjectID, clusterIDObs = clusterIDObs, stratumIDObs = stratumIDObs, Y = Y, Vf = Vf, Vr = Vr, LRModelType = LRModelType, LCClusterModelType = LCClusterModelType, LCStratumModelType = LCStratumModelType)
      }

    }

    # Posterior class assignment
    #We draw a new C: all parameters in the next iteration will be consistent with the new C
    postClass <- update_C_Marginal(priorPik = priorPik, beta = beta, Phi = phi, Omega = omega2, Psi = psi2, sigma2 = sigma2, Y = Y, Vf = Vf, Vr = Vr, subjectID = subjectID)

    C <- postClass[["C"]]
    pik <- postClass[["pik"]]


    if (i %% update == 0) {
      cat("Iteration:", i, "\n")
      cat("Class size:", table(C), "\n")
    }

    ### Store samples
    if (i > burn) {
      # Multinomial model
      store_delta[i - burn, ] <- c(delta)
      store_pi[i - burn, ] <- c(pik)
      store_priorPi[i - burn, ] <- c(priorPik)
      if (LCStratumModelType == "Random") {
        store_gamma2[i - burn, ] <- gamma2
      }
      if (LCClusterModelType %in% c("Both", "Unstr")) {
        store_tau2[i - burn, ] <- tau2
      }
      if (spline) {
        store_alpha[i - burn, ] <- c(alpha)
      }
      if (LCClusterModelType %in% c("Both", "Str")) {
        store_xi2[i - burn, ] <- xi2
      }

      # Longitudinal model
      if (LRModelType %in% c("Both", "Stratum")) {
        store_psi2[i - burn, ] <- psi2
      }
      if (LRModelType %in% c("Both", "Cluster")) {
        store_omega2[i - burn, ] <- omega2
      }
      store_beta[i - burn, ] <- c(beta)
      store_phi[i - burn, ] <- c(phi)
      store_sigma2[i - burn, ] <- sigma2


      # Additional samples to write to file
      if (writeSamples == TRUE) {

        #Samples of stratum level random effects in latent class membership model
        if (LCStratumModelType == "Random") {
          write.table(t(c(lambda)), file = "store_lambda.txt", sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        }
        #Samples of cluster level independent random intercepts from the class membership model
        if (LCClusterModelType %in% c("Both", "Unstr")) {
          write.table(t(c(u)), file = "store_u.txt", sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        }
        #Samples of cluster level spatial random intercepts from the class membership model
        if (LCClusterModelType %in% c("Both", "Str")) {
          write.table(t(c(nu)), file = "store_nu.txt", sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        }

        #Samples for b, rho, and zeta stored to recreate subject level trajectories
        if (LRModelType %in% c("Both", "Stratum")) {
          write.table(t(c(zeta)), file = "store_zeta.txt", sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        }
        if (LRModelType %in% c("Both", "Cluster")) {
          write.table(t(c(rho)), file = "store_rho.txt", sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        }
        #Stacked RE b_10, b_11, b_12 and so forth to match column names
        write.table(t(c(t(b))), file = "store_b.txt", sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)

        #Draws from posterior predictive distribution, redrawing latent class
        #!Ydraw <- get_postPred(K = K, beta = beta, sigma2 = sigma2, priorPik = priorPik, b = b, rho = rho, zeta = zeta, Vr = Vr, Vf = Vf, subjectID = subjectID, clusterIDObs = clusterIDObs, stratumIDObs = stratumIDObs)
        write.table(t(c(Ydraw)), file = "store_Ydraw.txt", sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)

        if (i > (n.samples - n.reps)) {
          write.table(store_redraw, file = paste0("store_redraw", "_", "iteration", i, ".txt"), sep = ",", row.names = FALSE, col.names = TRUE, append = FALSE)
        }

        #Record discrepancy measure T based on observed dataa and replicated data, must be computed above since C at iteration i-1 is consistent with parameters at i
        write.table(t(c(store_T)), file = "store_T.txt", sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)

      }

      ### Monitor selected parameters

      if (monitor == TRUE & ((i - burn) %% update == 0)) {
        par(mfrow = c(2, 2))
        plot(store_beta[1:(i - burn), 1], type = "l")
        plot(store_beta[1:(i - burn), 2], type = "l")
        plot(store_beta[1:(i - burn), 3], type = "l")
        plot(store_sigma2[1:(i - burn), 1], type = "l")
      }
    }
  }  # End of MCMC loop


  ### Register total time
  tottime <- proc.time() - ptm
  cat("Total minutes elapsed:", tottime/60, "\n", "\n")


  ### Write column names for writeSamples files

  if (writeSamples == TRUE) {
    ### Write column names to file for written to file samples
    if (LCStratumModelType == "Random") {
      colNameslambda <- paste(paste("Class", rep(2:K, each = M), sep = ""), paste("lambda", rep(1:M, times = (K - 1)), sep = ""), sep = "_")
      write.table(t(colNameslambda), file = "colNameslambda.txt", sep = ",", row.names = FALSE, col.names = FALSE, append = FALSE)
    }
    if (LCClusterModelType %in% c("Both", "Unstr")) {
      colNamesu <- paste(paste("Class", rep(2:K, each = J), sep = ""), paste("u", rep(1:J, times = (K - 1)), sep = ""), sep = "_")
      write.table(t(colNamesu), file = "colNamesu.txt", sep = ",", row.names = FALSE, col.names = FALSE, append = FALSE)
    }
    if (LCClusterModelType %in% c("Both", "Str")) {
      colNamesnu <- paste(paste("Class", rep(2:K, each = nJ), sep = ""), paste("nu", rep(1:nJ, times = (K - 1)), sep = ""), sep = "_")
      write.table(t(colNamesnu), file = "colNamesnu.txt", sep = ",", row.names = FALSE, col.names = FALSE, append = FALSE)
    }

    # Random effects
    if (LRModelType %in% c("Both", "Stratum")) {
      colNamesZeta <- paste(paste("Class", rep(1:K, each = M), sep = ""), paste("zeta", rep(1:M, times = K), sep = ""), rep(paste("RE", rep(1, each = M), sep = ""), times = K), sep = "_")
      write.table(t(colNamesZeta), file = "colNamesZeta.txt", sep = ",", row.names = FALSE, col.names = FALSE, append = FALSE)
    }
    if (LRModelType %in% c("Both", "Cluster")) {
      colNamesRho <- paste(paste("Class", rep(1:K, each = J), sep = ""), paste("rho", rep(1:J, times = K), sep = ""), rep(paste("RE", rep(1, each = J), sep = ""), times = K), sep = "_")
      write.table(t(colNamesRho), file = "colNamesRho.txt", sep = ",", row.names = FALSE, col.names = FALSE, append = FALSE)
    }
    colNamesb <- paste(paste("Subject", rep(1:n, each = q), sep = ""), rep(paste("RE", 1:q, sep = ""), times = n), sep = "_")
    write.table(t(colNamesb), file = "colNamesb.txt", sep = ",", row.names = FALSE, col.names = FALSE, append = FALSE)

  }

  ### Summary of model output

  # Model comparison summary
  #DIC4 calculation, with approximation using DIC2
  #Log likelihood evaluated at the posterior mode
  modeIndex <- get_posterior_mode(store_observed_posterior)
  llikYMode <- store_observed_llik[modeIndex]

  #DIC2 calculation following Celeux et al. 2006
  Dbar_2 <- mean(store_observed_llik)
  Dtilde_2 <- llikYMode
  pD_2 <- -2 * Dbar_2 + 2 * Dtilde_2
  DIC_2 <- -4 * Dbar_2 + 2 * Dtilde_2

  #DIC4 calculation, based on approximation in Celeux et al. 2006
  entropy_bar <- mean(store_entropy)
  DIC_4 <- DIC_2 + 2 * entropy_bar

  #DIC3 calculation
  #Sum of the log of the predictive density estimate per subject obtained by averaging over MCMC iterations
  Dbar_3 <- mean(store_observed_llik)
  Dtilde_3 <-  sum(log(apply(store_observed_llik_sub, 2, mean)))
  pD_3 <- -2 * Dbar_3 + 2 * Dtilde_3
  DIC_3 <- -4 * Dbar_3 + 2 * Dtilde_3

  #BIC calculation
  #Compute model dimension
  # K*p number of fixed effects \beta_k
  # K number of \sigma_2_k
  # s number of variables in W
  dK <- K * p + K + (K - 1) * s

  dSub <- (q * (q + 1) / 2) * K
  dK <- dK + dSub

  if (LRModelType %in% c("Stratum")) {
    dStrata <- K
    dK <- dK + dStrata
  } else if (LRModelType %in% c("Cluster")) {
    dClus <- K
    dK <- dK + dClus
  } else if (LRModelType %in% c("Both")) {
    dStrata <- K
    dClus <- K
    dK <- dK + dStrata + dClus
  }


  if (LCClusterModelType == "Unstr") {
    dK <- dK + (K - 1)
  } else if (LCClusterModelType == "Str") {
    dK <- dK + (K - 1)
  } else if (LCClusterModelType == "Both") {
    dK <- dK + 2 * (K - 1)
  }

  if (LCStratumModelType == "Random") {
    dK <- dK + (K - 1)
  }

  #Average effective sample size
  ESS_bar <- mean(store_ESS)

  ### BIC 1 is calculated using the log likelihood after integrating out the random effects and length of subjectID
  # Approximate the maximum likelihood estimator by taking the parameter values that maximize the log of the observed data likelihood over MCMC samples. This is per Fruhwirth-Schnatter, S., & Pyne, S. (2010). Bayesian Inference for Finite Mixtures of Univariate and Multivariate Skew-Normal and Skew-t Distributions. Biostatistics, 11(2), 317–336. https://doi.org/10.1093/biostatistics/kxp062
  # Per Jones (2011) SAS uses the number of clusters or subjects, i.e., n for the sample size of independent observations
  BIC_1 <- -2 * max(store_observed_llik) + dK * log(n)

  # BIC 2 is calculated using the log likelihood after integrating out the random effects and using the effective sample size (Jones 2011)
  BIC_2 <- -2 * max(store_observed_llik) + dK * log(ESS_bar)

  #ICL BIC_k
  ICL_BIC <- BIC_1 + 2 * store_entropy[which.max(store_observed_llik)]


  ### Add calculations for conditional predictive ordinate (see Gelfand and Dey 1994) and log pseudo marginal likleihood
  # n-length vector of CPOs
  #store_observed_llik_sub is \sum_{k = 1}^K \pi_{ik} f(y_i | \beta_k)
  CPO <- apply(store_observed_llik_sub, 2, function(x) {
    1 / mean(1 / x)
  })

  # Largest LPML indicates the better fitting model
  #Sum of log CPO across i = 1...n
  LPML <- sum(log(CPO))


  # Table of model comparison
  table_model_comparison <- matrix(c(BIC_1, BIC_2, ICL_BIC, DIC_3, DIC_4, LPML), nrow = 1)
  colnames(table_model_comparison) <- c("BIC1", "BIC2", "ICL-BIC", "DIC3", "DIC4", "LPML")
  rownames(table_model_comparison) <- c("value")

  # Table of effective number of parameters
  #NA for DIC4 due to approximation method from Celeux based on DIC2
  table_pD <- matrix(c(Dbar_3, Dtilde_3, pD_3, NA, NA, NA), nrow = 3, ncol = 2, byrow = FALSE)
  colnames(table_pD) <- c("DIC3", "DIC4")
  rownames(table_pD) <- c("Dbar", "Dtilde", "pD")


  # MCMC samples to return in a named list
  S <- list(store_delta = store_delta, store_gamma2 = store_gamma2, store_xi2 = store_xi2, store_tau2 = store_tau2, store_alpha = store_alpha, store_beta = store_beta, store_psi2 = store_psi2, store_omega2 = store_omega2, store_phi = store_phi, store_sigma2 = store_sigma2, store_pi = store_pi, store_priorPi = store_priorPi)


  # Posterior means and 95% Credible intervals after burn-in
  # Latent class model
  if (LCClusterModelType == "Both") {
    if (spline == TRUE) {
      objs <- c("store_delta", "store_xi2", "store_tau2", "store_alpha")
    } else {
      objs <- c("store_delta", "store_xi2", "store_tau2")
    }
  }

  if (LCClusterModelType == "Str") {
    if (spline == TRUE) {
      objs <- c("store_delta", "store_xi2",  "store_alpha")
    } else {
      objs <- c("store_delta", "store_xi2")
    }
  }

  if (LCClusterModelType == "Unstr") {
    if (spline == TRUE) {
      objs <- c("store_delta", "store_tau2",  "store_alpha")
    } else {
      objs <- c("store_delta", "store_tau2")
    }
  }

  if (LCClusterModelType == "None") {
    if (spline == TRUE) {
      objs <- c("store_delta", "store_alpha")
    } else {
      objs <- c("store_delta")
    }
  }


  if (LCStratumModelType == "Random") {
    objs <- c(objs, "store_gamma2")
  }

  #Append longitudinal outcomes model
  if (LRModelType == "Both") {
     objs <- c(objs, c("store_beta", "store_psi2", "store_omega2", "store_phi", "store_sigma2"))
  } else if (LRModelType == "Stratum") {
    objs <- c(objs, c("store_beta", "store_psi2", "store_phi", "store_sigma2"))
  } else if (LRModelType == "Cluster") {
    objs <- c(objs, c("store_beta", "store_omega2", "store_phi", "store_sigma2"))
  } else {
    objs <- c(objs, c("store_beta", "store_phi", "store_sigma2"))
  }

  temp <- lapply(S[objs], function(X){
    apply(X, 2, function(x){c(round(mean(x), digits = 4), round(quantile(x, probs = c(0.025, 0.975)), digits = 4))})
    })

  table_posterior_summary <- t(do.call("cbind", temp))
  colnames(table_posterior_summary) <- c("Post. Mean", "2.5 %", "97.5%")


  # Posterior classification
  postPi <- matrix(c(apply(S[["store_pi"]], 2, mean)), nrow = n, ncol = K, byrow = FALSE)
  postPredC <- apply(postPi, 1, which.max)

  #No. subjects in class k with post prob > 0.95, .90, .8
  sub_gr95 <- sub_gr90 <- sub_gr80 <-  rep(NA, K)
  for (k in 1:K) {
    temp <- postPi[postPredC == k, k]
    sub_gr95[k] <- sum(temp >= .95)
    sub_gr90[k] <- sum(temp >= .90)
    sub_gr80[k] <- sum(temp >= .80)
  }

  #Mean and median posterior probability by latent class
  meanPostPi <- medianPostPi <-  rep(NA, K)
  for (k in 1:K) {
    temp <- postPi[postPredC == k, k]
    meanPostPi[k] <- round(mean(temp), digits = 2)
    medianPostPi[k] <- round(median(temp), digits = 2)
  }

  table_latent_class_unit_assignment <- rbind(table(postPredC), sub_gr95, sub_gr90, sub_gr80, meanPostPi, medianPostPi)
  rownames(table_latent_class_unit_assignment) <- c("Predicted class size", "No. subjects with probability at least 0.95",
                                                    "No. subjects with probability at least 0.90",
                                                    "No. subjects with probability at least 0.80",
                                                    "Mean probability", "Median probability")
  colnames(table_latent_class_unit_assignment) <- c(paste("Class", 1:K))

  # Evaluation of label switching using Stephen's method
  temp <- array(S[["store_pi"]], dim = c(dim(store_pi)[1], n, K))
  perm_res <- label.switching::stephens(temp)

  #Number of rows not equal to the identity
  identity <- matrix(rep(c(1:K)), nrow = dim(store_pi)[1], ncol = K, byrow = TRUE)

  if (isTRUE(all.equal(perm_res[[1]], identity))) {

    cat("No evidence of label switching problem using Stephen's method from the label.switching package", "\n", "\n")

    } else {

    cat("Warning: Stephen's method from the label.switching package has detected label switching", "\n", "\n")

    }


  # Print output to screen
  cat("Background information", "\n")
  cat("Number of subjects:", n, "\n")
  cat("Number of observations:", length(subjectID), "\n")
  cat("Number of latent classes:", K, "\n")
  if (LCClusterModelType %in% c("Both", "Str")) {
    cat("Clusters in study area:", nJ, "\n")
  }
  cat("Number of area segments (clusters):", J, "\n", "\n")

  cat("Posterior latent class assignment:", "\n")
  print(table_latent_class_unit_assignment)

  cat("\n", "Reference class in latent class membership model:", 1, "\n")
  cat("Posterior means and 95% credible intervals:", "\n")
  print(table_posterior_summary)

  cat("\n", "Key to table of posterior means and 95% credible intervals:", "\n")
  if (LCStratumModelType == "Random") {cat("gamma2:", "variance of stratum-level random effects in latent class membership model", "\n")}
  if (LCClusterModelType %in% c("Both", "Str")) {cat("xi2:", "variance of cluster-level spatial random effects in latent class membership model", "\n")}
  if (LCClusterModelType %in% c("Both", "Unstr")) {cat("tau2:", "variance of cluster-level independent random effects in latent class membership model", "\n")}
  if (LRModelType %in% c("Both", "Stratum")) {cat("psi2:", "variance of stratum-level random effects in longitudinal outcomes model", "\n")}
  if (LRModelType %in% c("Both", "Cluster")) {cat("omega2:", "variance of cluster-level random effects in longitudinal outcomes model", "\n")}
  cat("phi:", "elements of variance-covariance of subject-level random effects in longitudinal outcomes model, indexed by row, then column", "\n")
  cat("sigma2:", "variance of observation-level in longitudinal outcomes model", "\n")

  cat("\n", "Model comparison statistics:", "\n")
  print(table_model_comparison)
  print(table_pD)

  return(S)

}
