\name{linbin2stg}
\alias{linbin2stg}
\title{Estimate binary-logistic parameters and odds ratios using linear regression with a single continuous Y-variable for two-phase sampled data. }
\description{
Fit location-scale model of the form \code{Y = eta + sigma*error} to data with a single continuous Y-variable and two-phase missingness structure, and convert to binary-logistic parameters and odds-ratio estimates with appropriate cut-points of Y. 
}
\usage{
linbin2stg(formula1, yCuts, lower.tail = TRUE, weights = NULL, 
	xstrata = NULL, data = list(), obstype.name = "obstype", 
	fit = TRUE, xs.includes = FALSE, compactX = FALSE, 
	start = NULL, Qstart = NULL, deltastart = NULL, 
	int.rescale = TRUE, control = mlefn.control(...), 
	control.inner = mlefn.control.inner(...), ...)
}
\arguments{
  \item{formula1}{A symbolic description of the location model to be fitted, i.e. eta.}
  \item{yCuts}{Cutpoint(s) used to define the binary Y-variable for logistic regression. Can be a matrix form (1*S) with S the number of xstrata.}
  \item{lower.tail}{If TRUE, define the cases being \code{{Y <= yCuts}}.}
  \item{weights}{An optional vector of weights to be used in the fitting process. Should be \code{NULL} or a numeric vector.}
  \item{xstrata}{Specify names of the stratification variables to be used, e.g. \code{"vname"} or \code{c("vname1","vname2",...)}. Strata are defined by cross-classification of all levels.}
  \item{data}{A data frame containing all the variables required for analysis, including those for \code{xstrata} and \code{obstype.name}.}
  \item{obstype.name}{Name of the variable specifying labels for observations by sampling and variable type: \code{"uncond"}, \code{"retro"}, \code{"xonly"}, \code{"y|x"} or \code{"strata"}.}
  \item{fit}{If \code{FALSE}, only stratum report will be generated without model fitting.}
  \item{xs.includes}{\code{TRUE} if weights specified for observations labelled as \code{"strata"} include those observed at the second phase (i.e. \code{"retro"} or \code{"uncond"} observations).}
  \item{compactX}{If \code{TRUE}, compress X matrix to distinct values with counts before model fitting. }
  \item{start}{Starting values for the regression parameters. Can be compusory if the program cannot produce a valid starting value at some situations.}
  \item{Qstart}{An optional starting matrix for Pr(Ystratum=i|Xstratum=j). Can be compulsory if the program cannot produce a valid starting value at some situations.}
  \item{deltastart}{An optional starting matrix for Pr(X=xk|Xstratum=j).}
  \item{int.rescale}{If \code{TRUE}, all X-variables will be standardised first before fitted in the model.}
  \item{control}{Specify control parameters for the iterations in \code{mlefn} call. See \code{\link{mlefn}} for details.}
  \item{control.inner}{Specify control parameters for inner iterations nested within \code{mlefn} call. See \code{\link{mlefn}} for details.}
  \item{\dots}{Further arguments passed to or from related functions.}
}
\details{
This function is a simple application of \code{locsc2stg} fitting linear regression models with a continuous Y using logistic error distribution. The results are then converted to much more efficient inferences about the same odds-ratio parameters being estimated by the logistic regression with the dichotomized binary outcome (case-control).  \cr\cr
More detailed descriptions of this function can be found in 
"Description of the \code{missreg} Library" (Wild and Jiang).
}
\value{
  \item{missReport}{Matrix containing information on deleted records with missing observations.}
  \item{StrReport}{Cross tabulation of counts for different levels of \code{obstype} and Y-values by X-strata.}
  \item{xStrReport}{Cross tabulation of counts for \code{obstype} by X-strata when \code{obstype="xonly"}.}
  \item{key}{Specify detailed classification for each of the X-strata.}
  \item{yCutsKey}{Specify the cutoff intervals for defined Y-strata within each X-stratum.}
  \item{fit}{\code{TRUE} or \code{FALSE} as its argument.}
  \item{error}{The error messages returned by \code{mlefn} call. Non-zero values indicate an unsuccessful fit.}
  \item{coefficients}{Linear regression coefficients.}
  \item{loglk}{Log-likelihood returned from final \code{mlefn} call.}
  \item{score}{Score vector returned from final \code{mlefn} call.}
  \item{inf}{Observed information matrix returned from final \code{mlefn} call.}
  \item{fitted}{The fitted values of Y obtained from the model.}
  \item{cov}{The asymptotic covariance matrix (inverse of the informnation matrix) of linear parameter estimates.}
  \item{cor}{The asymptotic correlation matrix of linear parameter estimates.}
  \item{bcoefficients}{Binary regression coefficients converted from linear parameters.}
  \item{bcov}{The asymptotic variance of binary parameter estimates.}
}
\references{Description of the \code{missreg} Library, Wild and Jiang, 2007.}
\author{Chris Wild, Yannan Jiang}
\note{The function \code{summary.linbin2stg} provides a complete summary of the regression results including the Wald tests and a regression panel for linear coefficients, a regression panel for binary coefficients, and associated odds-ratio estimates and confidence intervals.  All related output functions (\code{print.linbin2stg}, \code{summary.linbin2stg} and \code{print.summary.linbin2stg}) don't have 
help files provided at the moment. \cr\cr
Also note that the intercept of binary coefficients will not be available when more than one cut-point of Y is used, e.g. different for each x-stratum.}
\seealso{ \code{\link{locsc2stg}} }
\examples{
data(lowbirth.ls)
lowbirth.ls$sex.age <- interaction(lowbirth.ls$sex,lowbirth.ls$gest)
yCuts <- matrix(c(2550,2650,2740,2840,2900,3010,3030,3140),nrow=1)
yCut1 <- mean(yCuts)
 
### Multiple yCuts;
z1 <- linbin2stg(birthwt~gest+mumht+bmi+ethnicdb+hyper+smoke, 
                  yCuts=yCuts, xstrata=c("sex.age"), data=lowbirth.ls, 
                  obstype.name=c("instudy"), xs.includes=FALSE)

summary(z1)

### Single yCut;
z2 <- linbin2stg(birthwt~gest+mumht+bmi+ethnicdb+hyper+smoke, 
                  yCuts=yCut1, xstrata=c("sex.age"), data=lowbirth.ls, 
                  obstype.name=c("instudy"), xs.includes=FALSE)

summary(z2) 
}
\keyword{ Statistics|models|regression }

