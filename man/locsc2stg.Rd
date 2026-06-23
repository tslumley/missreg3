\name{locsc2stg}
\alias{locsc2stg}
\title{Linear regression with location-scale model for two-phase sampled data.}
\description{
Fits location-scale model of the form \code{Y = eta1 + exp(eta2)*error} to
data with a single continuous Y-variable and two-phase missingness structure, 
using the linear predictors \code{eta1} and \code{eta2} for specification of
the location and scale respectively.   
}
\usage{
locsc2stg(formula1, formula2, yCuts=NULL, weights=NULL, 
	  xstrata=NULL, data=list(), obstype.name="obstype",  
	  method="direct", fit=TRUE, errdistrn="logistic", 
	  errmodpars=6, xs.includes=FALSE, compactX=FALSE, 
	  compactY=TRUE, straty.maxnvals=20, start=NULL, 
	  Qstart=NULL, deltastart=NULL, int.rescale=TRUE,
          control=mlefn.control(...), 
          control.inner=mlefn.control.inner(...), ...)
}
\arguments{
  \item{formula1}{A symbolic description of the location model to be fitted 
		(eta1).}
  \item{formula2}{A symbolic description of the log-scale model to be fitted 
		(eta2). \code{~1} will fit a constant.}
  \item{yCuts}{Cutpoints used to define Y-strata. Critical when 
		\code{method="ycutmeth"}. Also required when \code{method="direct"} 
		but the starting values are not provided (See Details for more 
		descriptions).}
  \item{weights}{An optional vector of weights to be used in the fitting process.
		Should be \code{NULL} or a numeric vector.}
  \item{xstrata}{Specify names of the stratification variables to be used, e.g. 
		\code{"vname"} or \code{c("vname1","vname2",...)}. Strata are
		defined by cross-classification of all levels.}
  \item{obstype.name}{Name of the variable specifying labels for observations 
         	by sampling and variable type: \code{"uncond"}, \code{"retro"}, 
		\code{"xonly"}, \code{"y|x"} or \code{"strata"}.}
  \item{data}{A data frame containing all the variables required for analysis, 
		including those for \code{xstrata} and \code{obstype.name}.} 
  \item{method}{Two methods are implemented: \code{"ycutmeth"} and  \code{"direct"}
		(see Details for more descriptions).}
  \item{fit}{If \code{FALSE}, only stratum report will be generated without model fitting.}  
  \item{errdistrn}{A specification for the error distribution. Three choices are provided:
		standard logistic (\code{"logistic"}), standard normal (\code{"normal"})
		or student's-t distribution (\code{"t"}). 
		The default is \code{"logistic"}.}
  \item{errmodpars}{Set parameter values for the error distribution. The default is 6
		for student's-t distribution.}
  \item{xs.includes}{\code{TRUE} if weights specified for observations labelled as 
		\code{"strata"} include those observed at the second phase
	 	(i.e. \code{"retro"} or \code{"uncond"} observations).} 
  \item{compactX}{If \code{TRUE}, compress X matrix to distinct values with counts 
		before model fitting. This is only applicable to \code{method="direct"}.}
  \item{compactY}{If TRUE, limit the Y-values observed at the first phase
    		(\code{obstype="strata"}) to limited numbers of equally spaced 
		possible values. This is only applicable to \code{method="direct"}. }
  \item{straty.maxnvals}{If \code{compactY=TRUE}, specify the number of equally 
		spaced possible values spanning the range of Y observed as \code{"strata"}. 
		The default is 20.}
  \item{start}{Starting values for the regression parameters. Can be compusory if the
		program cannot produce a valid starting value at some situations.}
  \item{Qstart}{An optional starting matrix for Pr(Ystratum=i|Xstratum=j).
		Can be compulsory if the program cannot produce a valid starting
		value at some situations.}
  \item{deltastart}{An optional starting matrix for Pr(X=xk|Xstratum=j). This
		is only applicable to \code{method="direct"}.}
  \item{int.rescale}{If \code{TRUE}, all X-variables will be standardised first 
		before fitted in the model.}
  \item{control}{Specify control parameters for the iterations in \code{mlefn} call.
          	See \code{\link{mlefn}} for details.}
  \item{control.inner}{Specify control parameters for inner iterations nested 
		within \code{mlefn} call. See \code{\link{mlefn}} 
		for details.}
  \item{\dots}{Further arguments passed to or from related functions.}
}
\details{
This function fits location-scale models to continuous Y using different error 
distributions with various types of observations collected at different two-phase 
sampling schemes. More detailed descriptions of this function can be found in 
"Description of the \code{missreg} Library" (Wild and Jiang).

Two methods are implemented with either Y being categorical 
(\code{"ycutmeth"}) or at a continuous scale (\code{"direct"}). The argument 
\code{yCuts} is critical to the first approach but only required for the second 
approach when a starting value is needed. If \code{yCuts} is a vector, it defines 
the Y-strata with intervals (-infty, \code{yCuts}, infty). If \code{yCuts} 
is a matrix, the number of columns indicates the number of strata used and you
can define different cutpoints for each stratum. If you want to have differing 
numbers of cutpoints for different X-strata, pad out the bottom of any column 
that is not full with NAs.}
\value{
  \item{missReport}{Matrix containing information on deleted records with missing 
		observations.}
  \item{StrReport}{Cross tabulation of counts for different levels of \code{obstype}
		and Y-values by X-strata.}
  \item{xStrReport}{Cross tabulation of counts for \code{obstype} by X-strata 
		when \code{obstype="xonly"}.}
  \item{key}{Specify detailed classification for each of the X-strata.}
  \item{yCutsKey}{Specify the cutoff intervals for defined Y-strata within 
		each X-stratum.}
  \item{fit}{\code{TRUE} or \code{FALSE} as its argument.}
  \item{error}{The error messages returned by \code{mlefn} call. Non-zero values
		indicate an unsuccessful fit.}
  \item{coefficients}{The coefficients matrix with estimates, standard errors,
		z values and associated p-values.}
  \item{loglk}{Log-likelihood returned from final \code{mlefn} call.}
  \item{score}{Score vector returned from final \code{mlefn} call.}
  \item{inf}{Observed information matrix returned from final \code{mlefn} call.}
  \item{fitted}{The fitted values of Y obtained from the model.}
  \item{cov}{The asymptotic covariance matrix (inverse of the informnation matrix).}
  \item{cor}{The asymptotic correlation matrix.}
  \item{Qmat}{The estimated Pr(Ystratum=i|Xstratum=j) from the last iteration.}
  \item{deltamat}{The estimated \code{delta} matrix from the last iteration. \cr This
		is only applicable to \code{method="direct"}.}
}
\references{Description of the \code{missreg} Library, Wild and Jiang, 2007.}
\author{Chris Wild, Yannan Jiang}
\note{
The function \code{summary.locsc2stg} provides a complete summary of
the regression results including the Wald tests and a regression panel.  
All related output functions (\code{print.locsc2stg}, 
\code{summary.locsc2stg} and \code{print.summary.locsc2stg}) don't have 
help files provided at the moment.}

\seealso{ \code{\link{bin2stg}}}

\examples{
data(lowbirth.ls)
lowbirth.ls$sex.age <- interaction(lowbirth.ls$sex,lowbirth.ls$gest)
yCuts <- matrix(c(2550,2650,2740,2840,2900,3010,3030,3140),nrow=1)

z1 <- locsc2stg(birthwt ~ gest + mumht + bmi+ethnicdb+hyper+smoke, ~1,
		yCuts=yCuts, xstrata=c("sex.age"), data=lowbirth.ls, 
 		obstype.name=c("instudy"), xs.includes=FALSE, 
		method="ycutmeth")
summary(z1)

z2 <- locsc2stg(birthwt ~ gest + mumht + bmi+ethnicdb+hyper+smoke,~1,
               	xstrata=c("sex.age"),data=lowbirth.ls, 
	       	obstype.name=c("instudy"), xs.includes=FALSE, 
		method="direct", start=z1$coefficients, compactX=TRUE,
               	compactY=TRUE, straty.maxnvals=20)
summary(z2)

z2 <- locsc2stg(birthwt ~ gest + mumht + bmi+ethnicdb+hyper+smoke,~1,
		yCuts=yCuts, xstrata=c("sex.age"), data=lowbirth.ls, 
		obstype.name=c("instudy"), xs.includes=FALSE, 
		method="direct", start=z1$coefficients, Qstart=z1$Qmat, 
		compactX=TRUE, compactY=TRUE, straty.maxnvals=20, 
	       	control.inner=mlefn.control.inner(n.earlyit=3))
summary(z2)

z3 <- locsc2stg(birthwt ~ gest + mumht + bmi+ethnicdb+hyper+smoke,~1,
               	xstrata=c("sex.age"),data=lowbirth.ls, 
		obstype.name=c("instudy"), xs.includes=FALSE, 
		method="direct", start=z2$coefficients, 
		deltastart=z2$deltamat, compactX=TRUE,
          	compactY=TRUE, straty.maxnvals=100)
summary(z3)
}
\keyword{ Statistics|models|regression }
