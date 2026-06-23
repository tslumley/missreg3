\name{bin2stg}
\alias{bin2stg}
\title{Binary regression for two-phase sampled data }
\description{
Fits binary regression models to data with the two-phase missingness structure. 
This class includes stratified case-control data.
}
\usage{
bin2stg(formula, weights = NULL, xstrata = NULL, 
        obstype.name = "obstype", data, fit = TRUE, 
        xs.includes = FALSE, linkname = "logit", 
        start = NULL, Qstart = NULL, int.rescale = TRUE, 
        off.set = NULL, control = mlefn.control(...), 
        control.inner = mlefn.control.inner(...), ...)
}
\arguments{
  \item{formula}{A symbolic description of the model to be fitted.
		If there is only one non-NA level of the response variable presented
		in the data, that level is treated as "failure" (control).}
  \item{weights}{An optional vector of weights to be used in the fitting process.
		Should be \code{NULL} or a numeric vector.}
  \item{xstrata}{Specify names of the stratification variables to be used, e.g.
		\code{"vname"} or \code{c("vname1","vname2",...)}.
		Strata are defined by cross-classification of all levels.}
  \item{obstype.name}{Name of the variable specifying labels for observations by sampling
		and variable type: \code{"uncond"}, \code{"retro"}, \code{"xonly"},
		\code{"y|x"} or \code{"strata"}.}
  \item{data}{A data frame containing all the variables required for analysis,
		including those for \code{xstrata} and \code{obstype.name}.}
  \item{fit}{If \code{FALSE}, only stratum report will be generated without model fitting.
		\cr \cr This is useful in providing a data check, or finding internal
		ordering of the \code{xstrata} so that \code{yCuts} can be specified consistently
		with this ordering.}
  \item{xs.includes}{\code{TRUE} if \code{weights} specified for observations labelled as
		\code{"strata"} include those observed at the second phase (i.e. \code{"retro"}
		or \code{"uncond"} observations). }
  \item{linkname}{A specification for the model link function. Three choices are provided:
		\code{"logit"}, \code{"probit"} or \code{"cloglog"}. The default is \code{"logit"}.}
  \item{start}{Starting values for the regression parameters. Can be compusory if the program 
		cannot produce a valid starting value at some situations.
		\cr \cr When only part of the starting parameters are provided, names of these
		parameters will be used (if specified) to match the design matrix.
		Zeros will be used as starting values for all other parameters. This is useful
		when an updated fit is considered. }
  \item{Qstart}{An optional starting matrix for Pr(Y=i|Xstratum=j). The first row should
		be related to the successes (cases) and the second to the failures (controls).
		Can be compulsory if the program cannot produce a valid starting value at
		some situations.}
  \item{int.rescale}{If \code{TRUE}, all X variables will be standardised first before fitted
		in the model. }
  \item{off.set}{Specify an \code{a priori} known component to be included in the predictors.
		Should be \code{NULL} or a numeric vector.}
  \item{control}{Specify control parameters for the iterations in \code{mlefn} call. 
		See \code{\link{mlefn}} for details.}
  \item{control.inner}{Specify control parameters for inner iterations nested within 
		\code{mlefn} call. See \code{\link{mlefn}} for details.}
  \item{\dots}{Further arguments passed to or from related functions.}
}
\details{
This function fits binary regression models using several links with various types 
of observations collected at different two-phase sampling schemes. More detailed 
descriptions of the function and its applications can be found in "Description 
of the \code{missreg} Library" (Wild and Jiang).
}
\value{
  \item{missReport}{Matrix containing information on deleted records with missing observations.}
  \item{StrReport}{Cross tabulation of counts for different levels of \code{obstype} and Y-values 
		   by X-strata.}
  \item{xStrReport}{Cross tabulation of counts for \code{obstype} by X-strata when
		    \code{obstype="xonly"}.}
  \item{key}{Specify detailed classification for each of the X-strata.}
  \item{yKey}{Specify the Y variable and its level that the model is constructed for.}
  \item{fit}{\code{TRUE} or \code{FALSE} as its argument.}
  \item{error}{The error messages returned by \code{mlefn} call. Non-zero values indicate
		an unsuccessful fit.}
  \item{coefficients}{The coefficients matrix with estimates, standard errors, z values
		and associated p-values.}
  \item{loglk}{Log-likelihood returned from final \code{mlefn} call.}
  \item{score}{Score vector returned from final \code{mlefn} call.}
  \item{inf}{Observed information matrix returned from final \code{mlefn} call.}
  \item{fitted}{The fitted values of Y obtained by transforming the linear predictors 
		by the inverse of the link function.}
  \item{cov}{The asymptotic covariance matrix (inverse of the informnation matrix).}
  \item{cor}{The asymptotic correlation matrix.}
  \item{Qmat}{The estimated Pr(Y=i|Xstratum=j) from the last iteration.}
}
\references{Description of the \code{missreg} Library, Wild and Jiang, 2007.}
\author{Chris Wild, Yannan Jiang}
\note{The function \code{summary.bin2stg} provides a complete summary of 
the regression results including the Wald tests and a regression panel.
All related output functions (\code{print.bin2stg}, \code{summary.bin2stg}
and \code{print.summary.bin2stg}) don't have help files provided at the moment.}
\examples{
data(leprosy1)
leprosy1$age.trans <- 100 * (leprosy1$age + 7.5)^-2
z1 <- bin2stg(leprosy ~ age.trans + scar, data=leprosy1, weights=counts,
              xstrata="age", xs.includes=TRUE)
summary(z1)

data(leprosy2)
leprosy2$age.trans <- 100 * (leprosy2$age + 7.5)^-2
z2 <- bin2stg(cbind(case,control) ~ age.trans + scar, data=leprosy2,
	      xstrata="age", xs.includes=TRUE) 
summary(z2)

data(leprosy3)
leprosy3$age.trans <- 100 * (leprosy3$age + 7.5)^-2
z3 <- bin2stg(leprosy ~ age.trans + scar, data=leprosy3, weights=counts,
	      xs.includes=TRUE)
summary(z3)

data(wilms.sub)
z4 <- bin2stg(cbind(case,control) ~ stage*hist, xstrata=c("stage","inst"), 
              xs.includes=TRUE, data=wilms.sub)
summary(z4)

data(trawl)
attach(trawl)
# 265 out of 787 fish in fine net have length over 35  (caught37=NA)
# 353 out of 738 fish in test net have length over 35  (caught37=1)
# So 738 were caught from (estimate) 353*787/265 that entered
#est. pr(caught) assuming all fish over len=35 are caught
phat <- 738 / (787*353/265)  
                                         
z5 <- bin2stg(caught37 ~ I(length-35), weights=count, data=trawl,
          start=c(log(phat/(1-phat)),0), Qstart=matrix(c(phat,1-phat)))
summary(z5)

data(lowbirth.bin)
z6 <- bin2stg(sgagp~mumht+bmi+I(bmi^2) + ethnicdb + factor(occ)+ hyper + smoke,
          weights=counts, xstrata=c("ethnicdb","smokedb"),
          obstype.name=c("instudy"), data=lowbirth.bin, xs.includes=FALSE)
summary(z6)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Statistics|models|regression}
