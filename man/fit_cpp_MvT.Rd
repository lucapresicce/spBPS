% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{fit_cpp_MvT}
\alias{fit_cpp_MvT}
\title{Compute the parameters for the posteriors distribution of \eqn{\beta} and \eqn{\Sigma} (i.e. updated parameters)}
\usage{
fit_cpp_MvT(data, priors, coords, hyperpar)
}
\arguments{
\item{data}{\link{list} two elements: first named \eqn{Y}, second named \eqn{X}}

\item{priors}{\link{list} priors: named \eqn{\mu_B},\eqn{V_r},\eqn{\Psi},\eqn{\nu}}

\item{coords}{\link{matrix} sample coordinates for X and Y}

\item{hyperpar}{\link{list} two elemets: first named \eqn{\alpha}, second named \eqn{\phi}}
}
\value{
\link{list} posterior update parameters
}
\description{
Compute the parameters for the posteriors distribution of \eqn{\beta} and \eqn{\Sigma} (i.e. updated parameters)
}
