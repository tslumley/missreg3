\name{MEtaProspModInf}
\alias{MEtaProspModInf}
\title{Prospective Model Information function for models with M linear predictors.}
\description{
A sub-function called by \code{ML2Inf} to supply values and its derivatives for 
the first part of the profile loglikelihood regarding to the model of interest 
using the discrete partition version.
}
\usage{
MEtaProspModInf(theta,nderivs=2,y,x,wts=1,modelfn,off.set=0, ...)
}
\arguments{
  \item{theta}{Vector of the parameter values.}
  \item{nderivs}{Number of derivatives to be calculated, ranged from 0 (loglikelihood only)
to 2 (information matrix).}
  \item{y}{The response of interest, can be either a vector or matrix.}
  \item{x}{A 3-dimensional array (R*C*M) specifying the covariates values, 
with R the number of observations, C the length of \code{theta} and M the number 
of linear predictors used.}
  \item{wts}{An optional vector of weights (\code{n_i}) to be used in the fitting process. The
default is 1.}
  \item{modelfn}{A class of sub-functions called by \code{MEtaProspModInf} to calculate the values 
and their derivatives with respect to the linear predictor (\code{eta}'s) of X for the model of 
interest \code{f(Y|X; theta)}.}
  \item{off.set}{The offset provided in a matrix form (R*M) with R the number of 
observations and M the number of linear predictors used.}
  \item{\dots}{Further arguments passed to or from related functions.}
}
\details{
This sub-function is used to implement prospective regression models with a fixed number of
M linear predictors. It calculates the value and its derivatives for the first part of 
the profile loglikelihood in the form of \code{l*(theta,Q)} within each s-stratum \cr\cr
          \code{sum_{A(s)}{n_i^(s)*log{f(y_i^{(s)}|x_i^{(s)};theta)}} }, \cr\cr
with respect to \code{theta} through the M linear predictors (m=1,...,M), \cr\cr 
          \code{eta_{im} = o_{im}+x_{i(m)}^T*theta}
\cr \cr 
See "Description of the \code{missreg} Library" for all details.
}
\value{
  A list with the following components
  \item{loglk}{Log-likelihood obtained from the current \code{theta} values}
  \item{score}{Score vector obtained from the curent \code{theta} values when \code{nderivs>=1}; 
\code{NULL} otherwise.}
   \item{inf}{Observed information matrix obtained from the current \code{theta} values
when \code{nderivs=2}; \code{NULL} otherwise.}
}
\references{
Description of the \code{missreg} Library, Wild and Jiang, 2007.}
\author{Chris Wild, Yannan Jiang}
\seealso{\code{\link{ML2Inf}}}

\keyword{ Statistics|models|regression }
