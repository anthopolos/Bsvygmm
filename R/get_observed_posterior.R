#' @title Get the log observed data posterior.
#'
#' @description The reference for the approach is in Fruhwirth-Schnatter, S. and Pyne, S. (2010) Bayesian Inference for Finite Mixtures of Univariate and Multivariate Skew-Normal and Skew-t Distributions. Biostatistics, 11, 317-336. The observed data posterior uses the observed data likelihood after integrating out the random effects.

#' @param llik Current iteration of the log of the observed data likelihood.
#' @param beta Current iteration of fixed effects.
#' @param Phi Current iteration of between subject variance-covariance.
#' @param sigma2 Current iteration of observation level variance.
#' @param Omega Current iteration of between area segment level variance.
#' @param Psi Current iteration of between stratum level variance.
#' @param delta Current iteration of regression coefficients in latent class membership model.
#' @param alpha Current iteration of spline regression coefficients in latent class membership model.
#' @param u Current iteration of area segment unstructured random effects in latent class membership model.
#' @param deltaStratum Current iteration of stratum random effects in latent class membership model.
#' @param nu Current iteration of area segment structured random effects in latent class membership model.
#' @param tau2 Current iteration of variance of area segment unstructured random effects in latent class membership model.
#' @param gamma2 Current iteration of variance of stratum random effects in latent class membership model.
#' @param xi2 Current iteration of variance of area segment structured random effects in latent class membership model.
#' @param priorPik Current iteration of prior probability of class membership.
#' @param ADJ Adjacency matrix used in computing prior mean of area segment structured random effects.
#' @param priors List of prior distributions from model fitting.
#' @param hierVar List of priors for variance at observation level and hierarchical variances.
#' @param LCStratumModelType A string. For the latent class membership model, if \code{LCStratumModelType} is \code{Random}, then the latent class membership model includes stratum level random effects. If \code{Fixed}, then stratum level random effects are not included in the latent class membership model. Stratum fixed effects may be included in the design matrix \code{W}.
#' @param LCClusterModelType String for model type in the latent class membership model.
#' @param LRModelType String for model type in the longitudinal response model.
#' @param spline Logical. If \code{TRUE}, a spline is included in the latent class membership model.
#' @return Observed data likelihood at the current iteration.
get_observed_posterior <- function(llik, beta, Phi, sigma2, Omega, Psi, delta, alpha, u, deltaStratum, nu, tau2, gamma2, xi2, priorPik, ADJ, priors, hierVar, LCStratumModelType, LCClusterModelType, LRModelType, spline) {

  # Observed data posterior as a product of the observed data likelhood and the prior information, computed on the log scale
  # =& \bigg( \prod_{i = 1}^n \sum_{k = 1}^K \pi_{sjik} f(y_{sji} \ | \ \beta_k, \sigma^2_k, \Phi_k, \psi^2_k, \omega^2_k, v^f_{sji}, v^r_{sji}) \bigg) \prod_{k = 1}^K p(\beta_k) p(\sigma^2_k), p(\Phi_k) p(\omega^2_k) p(\psi^2_k) \prod_{k = 1}^{K-1} p(\delta_k) p(\delta0Cluster_k) p(\delta0Stratum_k)
  # =& \sum_{i = 1}^n \log \bigg( \sum_{k = 1}^K \pi_{sjik} f(y_{sji} \ | \ \beta_k, \sigma^2_k, \Phi_k, \psi^2_k, \omega^2_k, v^f_{sji}, v^r_{sji}) \bigg) +  \sum_{k = 1}^K \log p(\beta_k) p(\sigma^2_k), p(\Phi_k) p(\omega^2_k) p(\psi^2_k) + \sum_{k = 1}^{K-1} \log p(\delta_k) p(\delta0Cluster_k) p(\delta0Stratum_k)

  # Background for computation
  q <- dim(Phi)[1]
  K <- dim(priorPik)[2]

  if (LCStratumModelType == "Random") {
    M <- nrow(as.matrix(deltaStratum))
  }

  if (LCClusterModelType == "Unstr") {
    J <- nrow(as.matrix(u))
  } else if (LCClusterModelType == "Str") {
    nJ <- nrow(ADJ)
    ADJ <- ADJ == 1
  } else if (LCClusterModelType == "Both") {
    J <- nrow(as.matrix(u))
    nJ <- nrow(ADJ)
    ADJ <- ADJ == 1
  }

  ### Prior contributions to observed data posterior
  prior_beta <- prior_sigma2 <- prior_Omega <- prior_Psi <- prior_Phi <- rep(NA, K)
  prior_delta <- prior_gamma2 <- rep(NA, K - 1)

  if (LCStratumModelType == "Random") {
    prior_deltaStratum <- rep(NA, K - 1)
  }

  if (LCClusterModelType == "Unstr") {
    prior_u <- prior_tau2 <- rep(NA, K - 1)
  } else if (LCClusterModelType == "Str") {
    prior_nu <- prior_xi2 <- rep(NA, K - 1)
  } else if (LCClusterModelType == "Both") {
    prior_nu <- prior_xi2 <-  prior_u <- prior_tau2 <- rep(NA, K - 1)
  }

  if (spline == TRUE) {
    prior_alpha <- rep(NA, K - 1)
  }

  ### K loop for prior contributions by latent class in the longitudinal model
  for (k in 1:K) {

    Phik <- Phi[ , , k]

    # Prior contributions in longitudinal outcomes model
    prior_beta[k] <- mvtnorm::dmvnorm(beta[ , k], mean = priors[[6]][[1]], sigma = priors[[6]][[2]], log = FALSE)
    if (hierVar[[2]] == "IG") {
      prior_sigma2[k] <- 1 / dgamma(sigma2[k], priors[[10]][[1]], priors[[10]][[2]], log = FALSE)
    } else if (hierVar[[2]] == "Unif") {
      prior_sigma2[k] <- 1 / sqrt(sigma2[k])
    }

    if (LRModelType == "Both") {

      Omegak <- Omega[k]
      Psik <- Psi[k]

      if (hierVar[[1]] == "IG") {
        prior_Omega[k] <- 1 / dgamma(Omegak, priors[[8]][[1]], priors[[8]][[2]], log = FALSE)
        prior_Psi[k] <- 1 / dgamma(Psik, priors[[7]][[1]], priors[[7]][[2]], log = FALSE)
      } else if (hierVar[[1]] == "Unif") {
        prior_Omega[k] <- 1 / sqrt(Omegak)
        prior_Psi[k] <- 1 / sqrt(Psik)
      }

    } # End of Both condition


    if (LRModelType == "Stratum") {

      Psik <- Psi[k]

      if (hierVar[[1]] == "IG") {
        prior_Psi[k] <- 1 / dgamma(Psik, priors[[7]][[1]], priors[[7]][[2]], log = FALSE)
      } else if (hierVar[[1]] == "Unif") {
        prior_Psi[k] <- 1 / sqrt(Psik)
      }

    } # End of Stratum condition


    if (LRModelType == "Cluster") {

      Omegak <- Omega[k]

      if (hierVar[[1]] == "IG") {
        prior_Omega[k] <- 1 / dgamma(Omegak, priors[[8]][[1]], priors[[8]][[2]], log = FALSE)
      } else if (hierVar[[1]] == "Unif") {
        prior_Omega[k] <- 1 / sqrt(Omegak)
      }

    } # End of Cluster condition


    prior_Phi[k] <- MCMCpack::diwish(Phik, priors[[9]][[1]], priors[[9]][[2]])

  }

  # Prior contributions in latent class membership model
  for (k in 2:K) {
    prior_delta[k - 1] <- mvtnorm::dmvnorm(as.matrix(delta)[ , (k - 1)], mean = priors[[1]][[1]], sigma = priors[[1]][[2]], log = FALSE)

    if (LCStratumModelType == "Random") {
      prior_deltaStratum[k - 1] <- prod(dnorm(as.matrix(deltaStratum)[ , k - 1], mean = 0, sd = sqrt(gamma2[k - 1]), log = FALSE))
      if (hierVar[[1]] == "IG") {
        prior_gamma2[k - 1] <- 1 / dgamma(gamma2[k - 1], priors[[2]][[1]], priors[[2]][[2]], log = FALSE)
      } else if (hierVar[[1]] == "Unif") {
        prior_gamma2[k - 1] <- 1 / sqrt(gamma2[k - 1])
      }
    } # End LC Stratum condition

    if (spline == TRUE) {
      prior_alpha[k - 1] <- mvtnorm::dmvnorm(as.matrix(alpha)[ , (k - 1)], mean = rep(0, dim(as.matrix(alpha))[1]), sigma = priors[[5]][[1]], log = FALSE)
    }

    if (LCClusterModelType == "Unstr") {
      prior_u[k - 1] <- prod(dnorm(as.matrix(u)[ , k - 1], mean = 0, sd = sqrt(tau2[k - 1]), log = FALSE))
      if (hierVar[[1]] == "IG") {
        prior_tau2[k - 1] <- 1 / dgamma(tau2[k - 1], priors[[4]][[1]], priors[[4]][[2]], log = FALSE)
      } else if (hierVar[[1]] == "Unif") {
        prior_tau2[k - 1] <- 1 / sqrt(tau2[k - 1])
      }
    } # End of unstr condition

    if (LCClusterModelType == "Str") {
      nuk <- as.matrix(nu)[ , k - 1]
      munuk <- sapply(1:nrow(ADJ), function(x){mean(nuk[ADJ[x, ]])}, simplify = TRUE)
      prior_nu[k - 1] <- prod(dnorm(nuk, mean = munuk, sd = sqrt(xi2[k - 1]), log = FALSE))

      if (hierVar[[1]] == "IG") {
        prior_xi2[k - 1] <- 1 / dgamma(xi2[k - 1], priors[[3]][[1]], priors[[3]][[2]], log = FALSE)
      } else if (hierVar[[1]] == "Unif") {
        prior_xi2[k - 1] <- 1 / sqrt(xi2[k - 1])
      }
    } # End of Str condition

    if (LCClusterModelType == "Both") {
      prior_u[k - 1] <- prod(dnorm(as.matrix(u)[ , k - 1], mean = 0, sd = sqrt(tau2[k - 1]), log = FALSE))
      nuk <- as.matrix(nu)[ , k - 1]
      munuk <- sapply(1:nrow(ADJ), function(x){mean(nuk[ADJ[x, ]])}, simplify = TRUE)
      prior_nu[k - 1] <- prod(dnorm(nuk, mean = munuk, sd = sqrt(xi2[k - 1]), log = FALSE))

      if (hierVar[[1]] == "IG") {
        prior_tau2[k - 1] <- 1 / dgamma(tau2[k - 1], priors[[4]][[1]], priors[[4]][[2]], log = FALSE)
        prior_xi2[k - 1] <- 1 / dgamma(xi2[k - 1], priors[[3]][[1]], priors[[3]][[2]], log = FALSE)
      } else if (hierVar[[1]] == "Unif") {
        prior_tau2[k - 1] <- 1 / sqrt(tau2[k - 1])
        prior_xi2[k - 1] <- 1 / sqrt(xi2[k - 1])
      }

    } # End both condition


  }

  observed_posterior <- llik + sum(log(prior_beta), log(prior_Phi), log(prior_sigma2), log(prior_delta))

  if (spline == TRUE) {
    observed_posterior <- observed_posterior + sum(log(prior_alpha))
  }

  if (LCStratumModelType == "Random") {
    observed_posterior <- observed_posterior + sum(log(prior_deltaStratum), log(prior_gamma2))
  }

  if (LCClusterModelType == "Unstr") {
    observed_posterior <- observed_posterior + sum(log(prior_u), log(prior_tau2))
  } else if (LCClusterModelType == "Str") {
    observed_posterior <- observed_posterior + sum(log(prior_nu), log(prior_xi2))
  } else if (LCClusterModelType == "Both") {
    observed_posterior <- observed_posterior + sum(log(prior_u), log(prior_tau2), log(prior_nu), log(prior_xi2))
  }


  if (LRModelType == "Both") {
    observed_posterior <- observed_posterior + sum(log(prior_Omega), log(prior_Psi))
  } else if (LRModelType == "Stratum") {
    observed_posterior <- observed_posterior + sum(log(prior_Psi))
  } else if (LRModelType == "Cluster") {
    observed_posterior <- observed_posterior + sum(log(prior_Omega))
  }


  return(observed_posterior)

}
