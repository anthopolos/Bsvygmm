% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_post_pred.R
\name{get_postPred}
\alias{get_postPred}
\title{Get posterior predictive draws at each iteration.}
\usage{
get_postPred(
  K,
  C,
  b,
  rho,
  zeta,
  beta,
  sigma2,
  Vr,
  Vf,
  subjectID,
  clusterIDObs,
  stratumIDObs
)
}
\arguments{
\item{K}{Number of latent classes.}

\item{C}{Categorical variable for the latent class of each subject.}

\item{b}{A \code{q} column matrix of subject-level effects.}

\item{rho}{A \code{K} column matrix of area segment level intercepts. Equals \code{NULL} if not desired.}

\item{zeta}{A \code{K} column matrix of stratum level intercepts. Equals \code{NULL} if not desired.}

\item{beta}{Current iteration of fixed effects.}

\item{sigma2}{Current iteration of observation level variance.}

\item{Vr}{Design matrix for random effects.}

\item{Vf}{Design matrix for fixed effects.}

\item{subjectID}{Subject identifier for each observation.}

\item{clusterIDObs}{Area segment or cluster identifier for each observation.}

\item{stratumIDObs}{Stratum identifier for each observation.}
}
\value{
Draws from the posterior predictive distribution at the current iteration.
}
\description{
The reference for the approach is from Gelman BDA page 530-31. To conduct model checking, we take posterior predictive draws of \code{Y} conditional on \code{C}.
}
