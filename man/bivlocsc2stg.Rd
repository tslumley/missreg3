\name{bivlocsc2stg}
\alias{bivlocsc2stg}
\title{Bivariate binary-linear regression for two-phase sampled data}
\description{
Fits bivariate binary-linear regression models to data with two associated response
variables, binary Y1 and continuous Y2, and two-phase missingness structure.}
\usage{
bivlocsc2stg(formula1, formula2, formula3, weights = NULL, 
	     xstrata = NULL, data, obstype.name = "obstype", 
	     fit = TRUE, xs.includes = FALSE, off.set = NULL, 
	     errdistrn = "normal", errmodpars = 6, start = NULL, 
	     Qstart = NULL, control = mlefn.control(...), 
             control.inner = mlefn.control.inner(...), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula1}{A symbolic description of the model to be fitted for Y1|Y2, 
where Y1 is the binary response defining the case-control status of subjects
and Y2 is a continuous response of interest observed at the second phase.}
  \item{formula2}{A symbolic description of the location model to be fitted for Y2.}
  \item{formula3}{A symbolic description of the log-scale model to be fitted for Y2.}
  \item{weights}{An optional vector of weights to be used in the fitting process. 
Should be \code{NULL} or a numeric vector.}
  \item{xstrata}{Specify names of the stratification variables to be used, e.g.
\code{"vname"} or \code{c("vname1","vnames2",...)}. 
Strata are defined by cross-classification of all levels.}
  \item{obstype.name}{Name of the variable specifying labels for observations 
by sampling and variable type: \code{"uncond"},\code{"retro"},\code{"y|x"},
\code{"xonly"}, or \code{"strata"}.}
\item{data}{The data}
  \item{fit}{If \code{FALSE}, only stratum report will be generated without model fitting.}
  \item{xs.includes}{\code{TRUE} if weights specified for observations labelled as 
\code{"strata"} include those observed at the second phase (i.e. \code{"uncond"} or
\code{"retro"} observations).}
  \item{off.set}{Specify an \code{a} \code{priori} known component to be included in 
the predictors. Should be \code{NULL} or a numeric vector.}
  \item{errdistrn}{A specification for the erro distribution. Three choices are provided:
standard logistic (\code{"logistic"}), standard normal (\code{"normal"}) or student's-t
distribution (\code{"t"}). The default is \code{"logistic"}.}
  \item{errmodpars}{Set parameter values for the error distribution. The default is 6 for
student's-t distribution.}
  \item{start}{Starting values for the regression parameters.}
  \item{Qstart}{An optional starting matrix for Pr(Ystratum=i|Xstratum=j).}
  \item{control}{Specify control parameters for the iterations in \code{mlefn} call. 
See \code{\link{mlefn}} for details.}
  \item{control.inner}{Specify control parameters for inner iterations nested within
\code{mlefn} call. See \code{\link{mlefn}} for details.}
  \item{\dots}{Further arguments passed to or from related functions.}
}
\details{
  This function extends the application of \code{SPML2} method when Y2, the second response of interest 
  associated with Y1, is a continuous variable and ideal to be analysed under the location-scale model.
  In particular, we use a logistic regression model for Y1|Y2 as in \code{bivbin2stg} when the 
  \code{SPML2} method is applied, but a linear regression model for Y2 itself. Although the function 
  allows for different error distributions (\code{"logistic"}, \code{"normal"}, and \code{"t"} are 
  implemented so far), only the \code{normal} is assumed in the strata function and should be used at this stage.}

\value{
  \item{missReport}{Matrix containing information on deleted records with missing observations.}
  \item{StrReport}{Cross tabulation of counts for different levels of \code{obstype} and Y-values by
X-strata.}
  \item{xStrReport}{Cross tabulation of counts for \code{obstype} by X-strata when \code{obstype="xonly"}.}
  \item{key}{Specify detailed classification for each of the X-strata.}
  \item{ykey}{Vector containing the names of the Y-variables and the names of the 
              level of Ys the model is being constructed for. The sequence is as
              (name of Y1, name of the level at Y1=1, name of Y2).}
  \item{fit}{\code{TRUE} or \code{FALSE} as its argument.}
  \item{error}{The error messages returned by \code{mlefn} call. Non-zero values indicate an unsuccessful fit.}
  \item{coefficients}{The coefficients matrix with estimates, standard errors, z-values and associated
p-values.}
  \item{loglk}{Log-likelihood returned from final \code{mlefn} call.}
  \item{score}{Score vector returned from final \code{mlefn} call.}
  \item{inf}{Observed information matrix returned from final \code{mlefn} call.}
  \item{fitted}{The fitted values of Y2 obtained from the model.}
  \item{cov}{The asymptotic covariance matrix (inverse of the information matrix).}
  \item{cor}{The asymptotic correlation matrix.}
  \item{Qmat}{The estimated Pr(Ystratum)=i|Xstratum=j) from the last iteration.}
}
\references{Description of the \code{missreg} Library, Wild and Jiang, 2007.}
\author{Chris Wild, Yannan Jiang}
\note{The function \code{summary.bivlocsc2stg} gives a complete summary of the regression results 
including the Wald tests and a regression panel. All related output functions (\code{print.bivlocsc2stg},
 \code{summary.bivlocsc2stg} and \code{print.summary.bivlocsc2stg}) don't have help files provided
at the moment.}

\seealso{\code{\link{locsc2stg}}; \code{\link{bivbin2stg}}}
\examples{
## Data Generation ##
N <- 5000
x <- rnorm(N)
eps <- rnorm(N) 

theta2 <- c(0.5,1,0)
y2 <- theta2[1]+theta2[2]*x+exp(theta2[3])*eps

theta1 <- c(-3,-0.5,1,0.5)
eta1 <- theta1[1]+theta1[2]*y2+theta1[3]*x+theta1[4]*y2*x
p1 <- plogis(eta1)
y1 <- 1*(runif(N)<p1)

xcut <- c(-30,-1,0,1,30)
xstrata <- as.numeric(cut(x,xcut))

indca <- (1:N)[y1==1]
indct <- sample((1:N)[y1==0],length(indca))
ind <- sort(c(indca,indct))
rest <- (1:N)[-ind]
obstype <- rep("retro",N)
obstype[rest] <- "strata"
y2[rest] <- NA; x[rest] <- NA
dat <- data.frame(y1,y2,x,xstrata,obstype)

## Proportion of cases in population (about 0.1) ##
prca <- length(indca)/N
prca

## Model fit ##
z <- bivlocsc2stg(y1~y2*x,y2~x,~1,xstrata="xstrata",data=dat,xs.includes=FALSE)
summary(z)
}
\keyword{Statistics|models|regression}
