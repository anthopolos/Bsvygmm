#' @title Get the index for the mode of the density.
#'
#' @description The reference for the approach is in Fruhwirth-Schnatter, S., & Pyne, S. (2010). Bayesian Inference for Finite Mixtures of Univariate and Multivariate Skew-Normal and Skew-t Distributions. Biostatistics, 11(2), 317-336.
#'
#' @param samples Samples.
#' @return An index for the mode.
get_posterior_mode <- function(samples) {

  # Posterior mode using observed data posterior distribution
  modeEst <- density(samples)$x[which.max(density(samples)$y)]
  modeIndex <- which(abs(samples - (modeEst)) == min(abs(samples - modeEst)))

  return(modeIndex)

  #plot(density(samples))
  #abline(v = density(samples)$x[which.max(density(samples)$y)], col = "red")

}
