% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adj.R
\docType{data}
\name{ADJ}
\alias{ADJ}
\title{Spatial data for fitting a 2-class model with spatial and independent random effects to predict latent class membership}
\format{
A symmetric matrix with 50 rows and 50 columns.
}
\usage{
ADJ
}
\description{
An symmetric adjacency matrix containing the neighbors information, where 1 indicates a row and column are neighbors to each other and 0 otherwise. The diagonal is 0. The number of spatial units is assumed to equal the number of clusters in the longitudinal dataset.
}
\keyword{datasets}
