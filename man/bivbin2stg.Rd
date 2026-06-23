\name{bivbin2stg}
\alias{bivbin2stg}
\title{Bivariate binary regression for two-phase sampled data}
\description{
Fits bivariate binary regression models to data with two correlated binary
Y-variables and two-phase missingness structure. 
}
\usage{
bivbin2stg(formula1, formula2, formula3, weights = NULL, 
	   xstrata = NULL, obstype.name = "obstype", data, 
	   fit = TRUE, xs.includes = FALSE, y1samp = TRUE, 
	   method = "palmgren", start = NULL, Qstart = NULL, 
           off.set = NULL, control = mlefn.control(...), 
           control.inner = mlefn.control.inner(...), ...)
}
\arguments{
  \item{formula1}{A symbolic description of the model to be fitted for Y1, 
		the binary response defining the case-control status of subjects.
                When the \code{spml2} method is considered, it provides model formula
                for Y1|Y2 where Y2 is another binary response of interest 
		normally observed at the second phase.}
  \item{formula2}{A symbolic description of the model to be fitted for Y2, 
		the second binary response of interest correlated with Y1.}  
  \item{formula3}{A symbolic description of the model to be fitted quantifying 
		the association between Y1 and Y2. \code{~1} will fit a constant model.
                This model is not required when the \code{spml2} method is considered.}
  \item{weights}{An optional vector of weights to be used in the fitting process. Should
		be \code{NULL} or a numeric vector.}
  \item{xstrata}{Specify names of the stratification variables to be used, e.g. 
		\code{"vname"} or \code{c("vname1","vname2",...)}. Strata are defined
		by cross-classification of all levels.}
  \item{obstype.name}{Name of the variables specifying labels for observations by sampling 
		and variable type: \code{"uncond"},\code{"retro"},\code{"y|x"},
                \code{"xonly"}, or \code{"strata"}.}
  \item{data}{A data frame containing all the variables required for analysis, including
		those for \code{xstrata} and \code{obstype.name}.}
  \item{fit}{If \code{FALSE}, only stratum report will be generated without model fitting.}
  \item{xs.includes}{\code{TRUE} if weights specified for observations labelled as \code{"strata"} 
		include those observed at the second phase (i.e. \code{"retro"} or 
		\code{"uncond"} observations).}
  \item{y1samp}{\code{TRUE} if Y-strata are defined by the case-control information of Y1 
		in the population.  
                \code{FALSE} if Y-strata are defined by both Y1 and Y2 with either 
		"all controls" (Y1=0 and Y2=0) or not (Y1=1 or Y2=1).}
  \item{method}{Four methods are implemented: \code{"palmgren"}, \code{"bahadur"}, 
                \code{"copula"} and \code{"spml2"} (see Details for more descriptions). 
		Note that the last method is not available when \code{y1samp=FALSE}.}
  \item{start}{Starting values for the regression parameters in Y1-model, Y2-model 
		and the association model when applicable.}
  \item{Qstart}{An optional starting matrix for Pr(Ystratum=i | Xstratum=j). 
		Can be compulsory if the program cannot produce a valid starting value at 
		some situations.}
  \item{off.set}{Specify an \code{a priori} known component to be included in the predictors.
		Should be \code{NULL} or a numeric vector.}
  \item{control}{Specify control parameters for inner iterations nested within \code{mlefn} call.
                See \code{\link{mlefn}} for details.}
  \item{control.inner}{Specify control parameters for inner iterations nested within \code{mlefn}
                       call. See \code{\link{mlefn}} for details.}
  \item{\dots}{Further arguments passed to or from related functions.}
}
\details{
This function fits bivariate binary regresison to two correlated binary outcomes Y1 and 
Y2 using several models with various types observations collected at differnt two-phase
sampling schemes. 
\cr\cr
The joint distribution of Y1 and Y2 given X can be modelled using the marginal 
distributions of Pr(Y1|X) and Pr(Y2|X) along with an association model between Y1
and Y2. Currently implemented models for this approach include the Palmgren, Bahadur 
and Copula models. When we are only interested in Pr(Y2|X), another semiparametric 
approach (called \code{spml2} method) can be used in terms of a conditional factorisation 
Pr(Y1|Y2, X)*Pr(Y2|X) both treated parametrically. 
\cr\cr
More detailed descriptions of this function can be found in "Description of the 
\code{missreg} Library" (Wild and Jiang). 
}
\value{
  \item{missReport}{Matrix containing information on deleted records with missing observations.}
  \item{StrReport}{Cross tabulation of counts for different levels of \code{obstype} and 
           	Y-values by X-strata.}
  \item{xStrReport}{Cross tabulation of counts for \code{obstype} by X-strata when 
		\code{obstype="xonly"}.}
  \item{key}{Specify detailed classification for each of the X-strata.}
  \item{ykey}{Specify the Y-variables that the model is being constructed for.}
  \item{fit}{\code{TRUE} or \code{FALSE} as its argument.}
  \item{error}{The error messages returned by \code{mlefn} call. Non-zero values indicate
		an unsuccessful fit.}
  \item{coefficients}{The coefficients matrix with estimates, standard errors, z values
		and associated p-values. Will report separately for each marginal model 
		used.}
  \item{loglk}{Log-likelihood returned from final \code{mlefn} call.}
  \item{score}{Score vector returned from final \code{mlefn} call.}
  \item{inf}{Observed information matrix returned from final \code{mlefn} call.}
  \item{fitted.Y2}{The fitted values of Y2 obtained by transforming the linear predictors by
		the inverse of the link function. Note that all methods we have implemented 
		evaluate Pr(Y2|X) which is normally the model of interest.}
  \item{cov}{The asymptotic covariance matrix (inverse of the information matrix).}
  \item{cor}{The asymptotic correlation matrix.}
  \item{Qmat}{The estimated Pr(Ystratum=i|Xstratum=j) from the last iteration.}
  ...
}
\references{Description of the \code{missreg} Library, Wild and Jiang, 2007.}
\author{Chris Wild, Yannan Jiang}
\note{The function \code{summary.bivbin2stg} provides a complete summary of the regression
	results including the Wald tests and a regression model. All related output functions
	\code{print.bivbin2stg}, \code{summary.bivbin2stg} and \code{print.summary.bivbin2stg}
	don't have help files provided at the moment.}

\seealso{\code{\link{bin2stg}}}
\examples{
### SAMPLING ON CASE-CONTROL INFORMATION OF Y1 ONLY ###
data(cotdeath)
z1 <- bivbin2stg(y1~x, y2~x, ~x, weights=wts, data=cotdeath, 
                  xs.includes=TRUE, method="palmgren")
summary(z1)

z2 <- bivbin2stg(y1~x, y2~x, ~x, weights=wts, data=cotdeath, 
                  xs.includes=TRUE, method="bahadur")
summary(z2)
 
z3 <- bivbin2stg(y1~x, y2~x, ~x, weights=wts, data=cotdeath, 
                  xs.includes=TRUE, method="copula")
summary(z3)

z4 <- bivbin2stg(y1~x*y2, y2~x, weights=wts, data=cotdeath, 
                  xs.includes=TRUE, method="spml2")
summary(z4)

data(infarct)
z5 <- bivbin2stg(sgagp~ethnic+smoked+hyper+mumwt+mumwtc2+agepreg,
		anyinf~smoked+hyper+age1st, ~age1st, weights=count,
		xstrata=c("sex", "gest"), obstype.name="instudy",
		data=infarct, xs.includes=TRUE, method="palmgren")
summary(z5)

### SAMPLING ON ALL CONTROLS (Y1=0 AND Y2=0) OR NOT ###
data(dat00)
z6 <- bivbin2stg(y1~x, y2~x, ~x, weights=wts, data=dat00, y1samp=FALSE, 
                 xstrata="v", xs.includes=FALSE, method="palmgren")
summary(z6)
}
\keyword{Statistics|models|regression}

