% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/update_b.R
\name{update_b}
\alias{update_b}
\title{Update the subject level effects in the longitudinal model.}
\usage{
update_b(
  C,
  rho,
  zeta,
  beta,
  Phi,
  sigma2,
  Y,
  Vf,
  Vr,
  subjectID,
  clusterIDObs,
  stratumIDObs
)
}
\arguments{
\item{C}{Latent class assignments for each subject.}

\item{rho}{A \code{K} column matrix of area segment level intercepts. Equals \code{NULL} if not desired.}

\item{zeta}{A \code{K} column matrix of stratum level intercepts. Equals \code{NULL} if not desired.}

\item{beta}{A \code{K} column matrix of regression coefficients associated with \code{Vf}.}

\item{Phi}{A variance-covariance matrix for the subject level effects.}

\item{sigma2}{A \code{K} length vector of observation level variances.}

\item{Y}{Longitudinal measurements.}

\item{Vf}{A design matrix for fixed effects in the longitudinal outcome model.}

\item{Vr}{A design matrix for the subject-level random effects in the longitudinal outcome model. \code{Vr} is a subset or equal to \code{Vf}.}

\item{subjectID}{Subject identifier for each observation.}

\item{clusterIDObs}{Area segment or cluster identifier for each observation.}

\item{stratumIDObs}{Stratum identifier for each observation.}
}
\value{
A \code{q} column matrix of subject level intercepts.
}
\description{
This uses a normal prior distribution with mean 0 to update the subject level effects in the longitudinal outcomes model.
}
