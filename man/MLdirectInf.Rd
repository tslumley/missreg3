\name{MLdirectInf}
\alias{MLdirectInf}
\title{Likelihood calculation function for the direct approach with stratification }
\description{
An outer function of \code{ML2directInf} to provide the value, score vector and information
matrix at \code{theta} for the profile loglikelihood \code{l_P(theta)} of the form
\code{l(theta,delta)} with stratified two-phase sampled data. It reduces to an unstratified 
approach when nstrata=1.
}
\usage{
MLdirectInf(theta, nderivs = 2, deltamat = NULL, modelfn, 
	    hmodelfn, x, y, xStrat, Aposn, Acounts, Bposn, 
	    Bcounts, hvalue, Cmult, hxStrat, 
	    off.set = matrix(0, dim(x)[1], dim(x)[3]), extra = NULL, 
	    control.inner = mlefn.control.inner(...), ...)
}
\arguments{
  \item{theta}{Starting values for parameter \code{theta} in the regression model. }
  \item{nderivs}{Number of derivatives to be calculated. See \code{\link{ML2directInf}}
for details. }
  \item{deltamat}{The \code{delta_j^(s)} provided in a matrix form (J*S) 
with S the number of strata and J the number of distinct x-values observed.}
  \item{modelfn}{See \code{\link{ML2directInf}}. }
  \item{hmodelfn}{See \code{\link{ML2directInf}}. }
  \item{x}{See \code{\link{ML2directInf}}. }
  \item{y}{See \code{\link{ML2directInf}}. }
  \item{xStrat}{A vector of values 1 to S specifying the stratum membership of each
observation. }
  \item{Aposn}{A vector specifying the positions of those observations with the set
of complete (x,y)-values. }
  \item{Acounts}{A vector specifying the frequency of each observation (\code{n_i}) 
with the set of complete (x,y)-values. }
  \item{Bposn}{A vector specifying the positions of those observations with the x-values
observed; \code{NULL} in prospective sampling.}
  \item{Bcounts}{A vector specifying the frequency of each observation (\code{m_j}) with
the x-values observed; \code{NULL} in prospective sampling. }
  \item{hvalue}{The \code{h_k} in the loglikelihood. }
  \item{Cmult}{The \code{r_k} in the loglikelihood. }
  \item{hxStrat}{A vector of value 1 to S specifying the stratum membership of
each \code{hvalue}. }
  \item{off.set}{See \code{\link{ML2directInf}}. }
  \item{extra}{Provides \code{deltamat} from last iteration as starting values for next
inner iterative loop in \code{mlefn} function call. }
  \item{control.inner}{See \code{\link{ML2directInf}}. }
  \item{\dots}{Further arguments passed to or from related functions. }
}
\details{
This is the direct function called by \code{mlefn} to calculate the value, score vector
and observed information matrix at \code{theta} for the so-called profile loglikelihood 
\code{l_P(theta)} using the direct approach. It calls the inner function \code{ML2directInf} 
to evaluate \code{l^(s)(theta,delta^(s))} within each s-stratum.
}
\value{
  A list with the following components.
  \item{loglk}{Log-likelihood obtained from the current \code{theta} values.}
  \item{score}{Score vector obtained from the current \code{theta} values when \code{nderivs>=1};
\code{NULL} otherwise.}
  \item{inf}{Observed information matrix obtained from the current \code{theta} values when 
\code{nderivs=2}; \code{NULL} otherwise.}
  \item{extra}{A list providing updated \code{deltamat} values from previous iteration.}
}
\references{Description of the \code{missreg} Library, Wild and Jiang, 2007. }
\author{Chris Wild, Yannan Jiang }
\seealso{\code{\link{ML2directInf}}; \code{\link{mlefn}} }
\keyword{Statistics|models|regression}
