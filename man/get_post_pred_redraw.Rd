% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_post_pred_redraw.R
\name{get_post_pred_redraw}
\alias{get_post_pred_redraw}
\title{Get posterior predictive draws and level 1 and level 2 residuals at each iteration.}
\usage{
get_post_pred_redraw(
  K,
  C,
  b,
  sigma2,
  beta,
  zeta,
  rho,
  Phi,
  Omega,
  Psi,
  Z,
  delta,
  alpha,
  deltaStratum,
  u,
  nu,
  gamma2,
  tau2,
  xi2,
  priorPik,
  W,
  B,
  ADJ,
  clusterIDSub,
  stratumIDSub,
  spline,
  subjectID,
  clusterIDObs,
  stratumIDObs,
  Y,
  Vf,
  Vr,
  LRModelType,
  LCClusterModelType,
  LCStratumModelType
)
}
\value{
Draws of residuals at levels 1 and 2 for complete and replicated data, along with longitudinal outcome and latent class membership at the current MCMC iteration. Returns a data.frame with these values in addition to subject ID, cluster ID and stratum ID.
}
\description{
This approach is based on the following references: Wang, Chen-Pin, C. Hendricks Brown, and Karen Bandeen-Roche.
1. "Residual diagnostics for growth mixture models: Examining the impact of a preventive intervention on multiple trajectories of aggressive behavior." Journal of the American Statistical Association 100.471 (2005): 1054-1076.
2. Gelman, Andrew, et al. "Multiple imputation for model checking: completed‐data plots with missing and latent data." Biometrics 61.1 (2005): 74-85.
}
