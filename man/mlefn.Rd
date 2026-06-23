\name{mlefn}
\alias{mlefn}
\alias{mlefn.control}
\alias{mlefn.control.inner}
\title{A modified Newton maximiser}
\description{
A function to maximise, minimise or find stationary values for a (general) function.  It was
originally written to maximize a loglikelihood function so that is the language that is
employed.\cr \cr
\code{mlefn.control} and \code{mlefn.control.inner} supply parameter values to control 
the iterative process and reporting. They differ only in their defaults.
}
\usage{
mlefn(theta, loglkfn, control=mlefn.control(...), ...)

mlefn.control(messg="", niter=20, tol=1e-08, guide="uphill",
      print.progress=2, max.eigenrat=0.05, n.earlyit=0,
      constrain="no", fixed=NULL, Aconstrain=NULL, 
      cconstrain=NULL, ...)

mlefn.control.inner(messg="Inner:", niter=20, tol=1e-08, 
      guide="auto", print.progress=0, max.eigenrat=0.05, 
      n.earlyit=0, constrain="no", fixed=NULL, 
      Aconstrain=NULL, cconstrain=NULL, ...)
}
\arguments{
  \item{theta}{Starting values for the parameters of the loglikelhood function.}
  \item{loglkfn}{An inner function to compute the loglikelihood and its
     derivatives. The values returned by this function must be a list with
\code{loglk}, \code{score} and \code{inf}.  }
  \item{messg}{A labelling string to be printed as a part of warnings etc. 
Useful with nested calls to \code{mlefn}.}
  \item{niter}{Maximum number of iterations used. The default is 20.}
  \item{tol}{Level of tolerance for checking the convergence.}
  \item{guide}{Specification of the direction for convergence with the following choices:
      \describe{
        \item{\code{"uphill"}}{-- seek a maximum;}
        \item{\code{"downhill"}}{-- seek a minimum;}
        \item{\code{"no"}}{-- straight Newton approach without using loglikelihood 
		values to guide the search.}
        \item{\code{"auto"}}{-- only used in \code{mlefn.control.inner} when the inner
		function \code{loglkfn} requires a call to \code{mlefn} itself and we 
		want it to determine and supply a legitmate value for that (inner) call.}
      }
  }
  \item{print.progress}{A numeric value used to control the printing of error messages (if any);
0 should be used if no printing is required.}
  \item{max.eigenrat}{An argument used in the inner function \code{greenstadt.step}
to control the eigenvalues of the information matrix. This is the Greenstadt modification described
in page 601 of Seber and Wild (1989). }
  \item{n.earlyit}{Number of iterations to be treated as \code{"early"}; The default is 0.}
  \item{constrain}{Specification of constrain on the parameter estimates with the following choices:
      \describe{
        \item{\code{"no"}}{-- no constrains;}
        \item{\code{"fix"}}{-- fix some of the parameters at their starting values;}
      }
  }
  \item{fixed}{A vector specifying the parameters to be fixed, indicated by their orders in \code{theta}.
Used only if \code{constrain}=\code{"fix"}. }
  \item{Aconstrain}{an \code{I} matrix with number of rows equal to the number of "fixed" parameters. 
Used only if \code{constrain}=\code{"fix"}. }
  \item{cconstrain}{A vector specifying the values of those "fixed" parameters. \cr 
    Used only if \code{constrain}=\code{"fix"}. }
  \item{control}{to pass control options}
  \item{\dots}{Further arguments passed to or from related functions.}
}
\details{
This is the base function to maximise, minimise or find stationary values for \code{theta} 
using the provided \code{loglkfn} function. All semi-parametric maximum likelihood approaches
we have proposed in \code{missreg} library require this function to obtain maximum likelihood
estimates of parameters. See "Description of the \code{missreg} Library" for more details. 
}
\value{
  A list with the following components.
  \item{theta}{Updated parameter estimates at this iteration.}
  \item{loglk}{Log-likelihood obtained from the current \code{theta} values.}
  \item{score}{Score vector obtained from the current \code{theta} values.}
  \item{inf}{Observed information matrix obtained from the curent \code{theta} values.}
  \item{constrscore}{Constrained score vector if \code{constrain="fix"}; otherwise \code{NULL}.}
  \item{constrinf}{Constrained observed information matrix if \code{constrain="fix"}; otherwise \code{NULL}.}
  \item{counter}{Number of iterations performed.}
  \item{error}{A numeric value indicating the types of errors during iterations; 
a value of 0 indicates no error.} 
}
\references{
Nonlinear Regression, Seber and Wild,  1989. Wiley: New York. \cr
Description of the \code{missreg} Library, Wild and Jiang, 2007. 
}
\author{Chris Wild, Yannan Jiang}
\keyword{ Statistics|models|regression }
