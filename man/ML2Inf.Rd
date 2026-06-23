\name{ML2Inf}
\alias{ML2Inf}
\title{ Core likelihood calculation function for the discrete partition version }
\description{
  A sub-function called by \code{MLInf} to provide the value, score vector and
information matrix at \code{theta} for the so-called profile loglikelihood
\code{l*(theta, Q)} within each stratum with stratified two-phase
sampled data. It reduces to an unstratified approach when nstrata=1.
}
\usage{
ML2Inf(theta, nderivs = 2, ProspModInf, StratModInf, x, y, Aposn, 
       Acounts, Bposn, Bcounts, rvec, Qs, usage = "thetaonly", 
       thetaparts = 0, paruse = "auto", inxStrat, 
       off.set = matrix(0, dim(x)[1], dim(x)[3]), 
       control.inner = mlefn.control.inner(...), ...)
}
\arguments{
  \item{theta}{ Starting values for parameters \code{theta} in the regression
model unless \code{usage} below is \code{"combined"}, in which case \code{theta}
should contain starting values for the \code{theta} parameters followed by starting
values fo the \code{rho} or \code{xi} parameters (which will be stripped off on entry
to the function). }
  \item{nderivs}{Number of derivatives to be calculated. Either 0 (loglikelihood
value only), 1 (also return score vector), or 2 (also return information matrix). }
  \item{ProspModInf}{A class of sub-functions called by \code{ML2Inf} to supply
values and their derivatives for the first (A) part of the loglikelihood regarding
to the model of interest.  }
  \item{StratModInf}{A class of sub-functions called by \code{ML2Inf} to supply
values and their derivatives for the third (B) part of the loglikelihood 
regarding to the \code{h}-distribution at the first phase. }
  \item{x}{See \code{\link{MEtaProspModInf}}.  }
  \item{y}{See \code{\link{MEtaProspModInf}}.  }
  \item{Aposn}{A vector specifying the positions of those observations contributed to 
the A part of the loglikelihood in the data matrix. }
  \item{Acounts}{A vector specifying the frequency of each observation (\code{n_i}) 
contributed to the A part of the loglikelihood.  }
  \item{Bposn}{A vector specifying the positions of those observations contributed to
the B part of the loglikelihood in the data matrix.  }
  \item{Bcounts}{A vector specifying the frequency of each observation (\code{m_j}) 
contributed to the B part of the loglikelihood.  }
  \item{rvec}{The \code{r_k^(s)} in the loglikelihood.  }
  \item{Qs}{The \code{Q_k^(s)} in the loglikelihood.  }
  \item{usage}{Work with and report results for the following three conditions:
(1) \code{"thetaonly"} (profile other parameters); (2) \code{"combined"} (fit both 
\code{theta} and \code{rho/xi} simultaneously); and (3) \code{"Qfixed"} (fix other parameters).  }
  \item{thetaparts}{A vector of length 2 specifying the number of \code{theta} and \code{rho/xi} 
parameters as appropriate; Used only if \code{usage="combined"}.  }
  \item{paruse}{The choice of using either \code{rho} or \code{xi} parameters as follows: 
(1) \code{"rhos"}; (2) \code{"xis"}; or (3) \code{"auto"} (function to choose the \code{rho's}
if no more than one \code{r_k^(s)=0} and the \code{xi's} otherwise). Any other string is treated 
as the last option.  }
  \item{inxStrat}{Optional to enable printing a diagnostic when \code{ML2Inf} fails
and has been called from \code{MLInf}. }
  \item{off.set}{See \code{\link{MEtaProspModInf}}.  }
  \item{control.inner}{Specify control parameters for inner iterations nested within
the \code{mlefn} function call. }
  \item{\dots}{Further arguments passed to or from related functions.  }
}
\details{
This is the core function at the distrete partition version to calculate the value, 
score vector and observed information matrix at \code{theta} for the so-called profile
loglikelihood \code{l*^(s)(theta,Q^(s))}. It is an inner function called by \code{MLInf}
to supply values of \code{l*} and its derivatives within each s-stratum. 
\cr \cr
See "Description of the \code{missreg} Library" for more details.
}
\value{
  A list with the following components.
  \item{loglk }{Log-likelihood obtained from the current \code{theta} values.}
  \item{score }{Score vector obtained from the current \code{theta} values when 
\code{nderivs>=1}; \code{NULL} otherwise.}
  \item{inf }{Observed information matrix obainted from the current \code{theta}
values when \code{nderivs=2}; \code{NULL} otherwise.}
  \item{Qs}{A vector of length K providing values for \code{Q_k^(s)} either as 
its entry values or updated from the inner iterative process when \code{r_k^(s)} is
not equal to 0 and \code{usage="thetaonly"}.}
}
\references{Description of the \code{missreg} Library, Wild and Jiang, 2007. }
\author{Chris Wild, Yannan Jiang }
\seealso{\code{\link{MLInf}}; \code{\link{MEtaProspModInf}} }
\keyword{ Statistics|models|regression }
