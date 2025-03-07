% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/optimization_functions.R
\name{bootstrap_glm}
\alias{bootstrap_glm}
\title{Random starting values for HME}
\usage{
bootstrap_glm(n = 2, ...)
}
\arguments{
\item{n}{The number of random parameter vectors to produce}

\item{...}{Values to pass to glm()}
}
\description{
Estimates a GLM model and then uses the fitted parameters and covariance
matrix to randomly sample n new parameter vectors.
}
