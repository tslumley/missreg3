\name{MLInf}
\alias{MLInf}
\title{Likelihood calculation function for the discrete partition version with stratification }
\description{
An outer function of \code{ML2Inf} to provide the value, score vector and information matrix
at \code{theta} for the so-called profile loglikelihood \code{l_P(theta)} of the form 
\code{l*(theta, Q)} with strtified two-phase sampled data. It reduces to an unstratified 
approach when nstrata=1. 
}
\usage{
MLInf(theta, nderivs = 2, ProspModInf, StratModInf, x, y, 
      Aposn, Acounts, Bposn, Bcounts, rmat, Qmat, 
      xStrat = rep(1, dim(x)[1]), extra = NULL, 
      off.set = matrix(0, dim(x)[1], dim(x)[3]), 
      control.inner = mlefn.control.inner(...), ...)
}
\arguments{
  \item{theta}{Starting values for parameter \code{theta} in the regression model. 
See \code{\link{ML2Inf}} for details. }
  \item{nderivs}{Number of derivatives to be calculated. }
  \item{ProspModInf}{See \code{\link{ML2Inf}}. }
  \item{StratModInf}{See \code{\link{ML2Inf}}. }
  \item{x}{See \code{\link{ML2Inf}}. }
  \item{y}{See \code{\link{ML2Inf}}. }
  \item{Aposn}{See \code{\link{ML2Inf}}. }
  \item{Acounts}{See \code{\link{ML2Inf}}. }
  \item{Bposn}{See \code{\link{ML2Inf}}. }
  \item{Bcounts}{See \code{\link{ML2Inf}}. }
  \item{rmat}{The \code{r_k^(s)} provided in a matrix form (K*S) with S 
the number of strata and K number of distinct \code{h}-values observed.}
  \item{Qmat}{The \code{Q_k^(s)} provided in a matrix form (K*S). }
  \item{xStrat}{A vector of values 1 to S specifying the stratum membership of 
each observation.}
  \item{extra}{Provides Qmat from last iteration as starting values for next 
inner iterative loop in \code{mlefn} function call.}
  \item{off.set}{See \code{\link{ML2Inf}}. }
  \item{control.inner}{See \code{\link{ML2Inf}}. }
  \item{\dots}{Further arguments passed to or from related functions. }
}
\details{
This is the direct function called by \code{mlefn} to calculate the value, 
score vector and observed information matrix at \code{theta} for the so-called
profile loglikelihood \code{l_P(theta)} using the discrete partition version. 
It calls the inner function \code{ML2Inf} to evaluate 
\code{l*^(s)(theta,Q^(s))} within each s-stratum.  }
\value{
  A list with the following components.
  \item{loglk}{Log-likelihood obtained from the current \code{theta} values.}
  \item{score}{Score vector obtained from the current \code{theta} values when 
\code{nderivs>=1}; \code{NULL} otherwise.}
  \item{inf}{Observed information matrix obtained from the current \code{theta}
values when \code{nderivs=2}; \code{NULL} otherwise.}
  \item{extra}{A list providing updated Qmat values from previous iteration.}
}
\references{Description of the \code{missreg} Library, Wild and Jiang, 2007. }
\author{Chris Wild, Yannan Jiang}
\seealso{\code{\link{ML2Inf}}; \code{\link{mlefn}} }

\keyword{ Statistics|models|regression }

