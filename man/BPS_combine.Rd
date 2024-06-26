% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{BPS_combine}
\alias{BPS_combine}
\title{Combine subset models wiht BPS}
\usage{
BPS_combine(fit_list, K, rp)
}
\arguments{
\item{fit_list}{\link{list} K fitted model outputs composed by two elements each: first named \eqn{epd}, second named \eqn{W}}

\item{K}{\link{integer} number of folds}

\item{rp}{\link{double} percentage of observations to take into account for optimization (\code{default=1})}
}
\value{
\link{matrix} posterior predictive density evaluations (each columns represent a different model)
}
\description{
Combine subset models wiht BPS
}
