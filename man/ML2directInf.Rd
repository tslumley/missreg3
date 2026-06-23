\name{ML2directInf}
\alias{ML2directInf}
\title{Core likelihood calculation function for the direct approach }
\description{
A sub-function called by \code{MLdirectInf} to provide the value, score vector 
and information matrix at \code{theta} for the so-called profile loglikelihood 
\code{l_P(theta)} of the form \code{l(theta, delta)} within each s-stratum
with stratified two-phase sampled data. It reduces to an unstratified approach when
nstrata=1. 
}
\usage{
ML2directInf(theta, nderivs = 2, modelfn, hmodelfn, x, y, Aposn, 
	     Acounts, Bposn, Bcounts, hvalue, Cmult, delta, 
	     off.set = matrix(0, dim(x)[1], dim(x)[3]), inxStrat, 
	     control.inner = mlefn.control.inner(...), ...)
}
\arguments{
  \item{theta}{Starting values for parameters \code{theta} in the regression model.}
  \item{nderivs}{Number of derivatives to be calculated. Either 0 (loglikelihood
value only), 1 (also return score vector), or 2 (also return information matrix). }
  \item{modelfn}{A class of sub-functions called by \code{ML2directInf} to supply 
values and their derivatives with respect to the \code{eta}'s (M linear predictors 
with respect to \code{theta}) for the model of interest \code{f(Y|X;theta)}.  }
  \item{hmodelfn}{A class of sub-functions called by \code{ML2directInf} to supply
values and their derivatives with respect to the \code{eta}'s for 
\code{pr(h_k|x_j;theta)} under the same class of models. } 
  \item{x}{A 3-dimensional array (R*C*M) specifying the covariates values,
with R the number of observations, C the length of \code{theta} and M the number of 
linear predictors used. }
  \item{y}{The response of interest, can be eitehr a vector or matrix. }
  \item{Aposn}{A vector specifying the positions of those observations with the set 
of complete (x, y)-values from s-stratum. }
  \item{Acounts}{A vector specifying the frequency of each observation (\code{n_i})
with the set of complete (x,y)-values from s-stratum. } 
  \item{Bposn}{A vector specifying the positions of those observations with the 
x-values observed in s-stratum; \code{NULL} in prospective sampling.}
  \item{Bcounts}{A vector specifying the frequency of each observation (\code{m_j})
with the x-values observed in s-stratum; \code{NULL} in prospective sampling.}
  \item{hvalue}{The \code{h_k^(s)} in the loglikelihood. }
  \item{Cmult}{The \code{r_k^(s)} in the loglikelihood. }
  \item{delta}{The \code{delta_j^(s)} in the loglikelihood. }
  \item{off.set}{The offset provided in a matrix form (R*M) with R the number
of observations and M the number of linear predictors used. }
  \item{inxStrat}{See \code{\link{ML2Inf}}. }
  \item{control.inner}{Specify control parameters for inner iterations nested within
the \code{mlefn} function call. See \code{\link{mlefn}} for details.}
  \item{\dots}{Further arguments passed to or from related functions. }
}
\details{
This is the core function in the direct approach to calculate the value, score vector
and observed information matrix at \code{theta} for the profile loglikelihood \code{l_P(theta)}
of the form \cr
\code{l^(s)(theta,delta^(s))} within each s-stratum. \cr\cr
It is an inner function called by \code{MLdirectInf}.
}
\value{
  A list with the following components.
  \item{loglk}{Log-likelihood value obtained from the current \code{theta} values.}
  \item{score}{Score vector obtained from the current \code{theta} values when \code{nderivs>=1};
\code{NULL} otherwise.}
  \item{inf}{Observed information matrix obtained from the current \code{theta} values when
\code{nderivs=2}; \code{NULL} otherwise.}
  \item{delta}{A vector of length J providing values for \code{delta_j^(s)} either as its entry values
or updated from the inner iterative process when \code{Bposn} is not \code{NULL}.}
}
\references{Description of the \code{missreg} Library, Wild and Jiang, 2007. }
\author{Chris Wild, Yannan Jiang }
\seealso{\code{\link{MLdirectInf}}; \code{\link{ML2Inf}}  }
\keyword{Statistics|models|regression }
