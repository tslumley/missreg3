\name{rclusbin2}
\alias{rclusbin2}
\docType{data}
\title{
Random intercept model for clustered binary data following case-control sampling.
}
\description{
Fits random intercept models to clustered binary data after case and control sampling, wherein interest is in the relationship between a binary response (Y) that is related to the sampling variable (Z).
}
\usage{
rclusbin2(formula1, formula2, weights=NULL, ClusInd.name=NULL, IntraClus.name=NULL, yname, xstrata=NULL, ystrata, 
obstype.name="obstype", data, NMat=NULL, xs.includes=FALSE, MaxInClus=NULL, rmsingletons=FALSE, retrosamp=TRUE, 
nzval0=20, fit=TRUE, devcheck=FALSE, linkname="logit", start=NULL, Qstart=NULL, sigma=NULL, paruse="xis", 
control=mlefn.control(...), control.inner=mlefn.control.inner(...), ...)
}
\arguments{
  \item{formula1}{A symbolic description of the random intercept model to be fitted, i.e. the model of interest.}
  \item{formula2}{A symbolic description of the auxiliary model to be fitted, between the sampling (case-control) variable and the binary response of interest.}
  \item{weights}{An optional vector of weights to be used in the fitting process. Should be \code{NULL} or the name of a numeric vector in the data frame. It provides weights at the individual level when there are clusters of size greater than one. \cr When all clusters are of size one, it provides weights at cluster=individual level.}
  \item{ClusInd.name}{Name of a vector in the data frame specifying cluster membership. Can be \code{NULL} if all clusters are of size one.}
  \item{IntraClus.name}{Name of a vector in the data frame specifying intra-cluster sequence of individual subjects in a cluster. The one with the smallest i.d. is treated as the proband who were originally sampled into a study.}
  \item{yname}{Name of the binary response variable of interest in the data frame. Must be specified.}
  \item{xstrata}{Specify names of the stratification variables to be used, e.g. \code{"vname"} or \code{c("vname1","vname2",...)}. Strata are defined by cross-classifiction of all levels. \cr This function only deals with the situation when clusters are defined within xstrata.}
  \item{ystrata}{Specify name of the variable defining the case and control strata. }
  \item{obstype.name}{Name of the variable specifying labels for observations by sampling and variable type: \code{"uncond"}, \code{"retro"}, \code{"xonly"}, \code{"y|x"} or \code{"strata"}. }
  \item{data}{A data frame containing all the variables required for analysis, including those for \code{xstrata}, \code{ystrata} and 
		\code{obstype.name}.}
  \item{NMat}{Population counts in a matrix form with rows and columns corresponding to case-control strata and X-strata respectively. Should not be provided when there is any observation of the type \code{"strata"}.}
  \item{xs.includes}{\code{TRUE} if \code{weights} specified for obervations labelled as "strata" include those observed at the second phase (i.e. \code{"retro"} or \code{"uncond"} observations).}
  \item{MaxInClus}{A value specifying the maximum number of individuals allowed in a cluster. Set to \code{NULL} if there is no limit.}
  \item{rmsingletons}{If \code{TRUE}, remove clusters of size one.}
  \item{retrosamp}{As the default, must be \code{TRUE} here.}
  \item{nzval0}{Number of points to calculate the zeros and weights needed for Gauss-Hermite quadrature. The default is 20.}
  \item{fit}{If \code{FALSE}, only stratum report will be generated without model fitting.}
  \item{devcheck}{If \code{TRUE}, check the first and second derivatives. The default should be \code{FALSE}.}
  \item{linkname}{A specification for the model link function. Three choices are provide: \code{"logit"}, \code{"probit"} or "cloglog". The default is \code{"logit"}. }
  \item{start}{Starting values for all regression parameters. }
  \item{Qstart}{An optional starting matrix for Pr(Zstratum=i|Xstratum=j). Can be compulsory if the program cannot produce a valid starting value at some situations. }
  \item{sigma}{An optional starting value for \code{sigma}. The default value (when set to \code{NULL}) is 0.5.}
  \item{paruse}{As the default, must be \code{"xis"} here.}
  \item{control}{Specify control parameters for the iterations in \code{mlefn} call.}
  \item{control.inner}{Specify control parameters for inner iterations nested within \code{mlefn} call. }
  \item{\dots}{Further arguments passed to or from related function. }
}
\details{
To be added.
}
\source{
http://www.stat.auckland.ac.nz/~wild
}
\references{
Longitudinal Studies of Binary Response Data Following Case-Control and Stratified Case-Control Sampling: Design and Analysis, Schildcrout and Rathouz, BIOMETRICS 2009.
}
\seealso{ \code{\link{ghq}},  \code{\link{rclusbin}}}
\examples{
data(adhd)
head(adhd)

adhd$obstype <- rep("retro", dim(adhd)[1])
adhd$probandS <- 2 - adhd$proband #as 1/2 for case/control
adhd$sexF <- adhd$sex-1 #as 1/0 for female/male
adhd$wave1 <- ifelse(adhd$wave==1, 1, 0)
adhd$wave2 <- ifelse(adhd$wave==2, 1, 0)

adhd1 <- adhd[adhd$wave==1,]
z0 <- glm(proband ~ adhd, family=binomial, data=adhd1)
z0$coefficients

nMat <- ftable(adhd1$sex~adhd1$probandS)  # 1=male; 2=female;
nMat

### Samping ratios for boys/girls (Schildcrout & Rathouz)
pi_ctF <- 1/22.6
pi_ctM <- 1/22.4
NMat <- cbind(c(113, 96/pi_ctM), c(25, 21/pi_ctF))

z <- rclusbin2(adhd ~ wave1+wave2+wave+sexF+african+other+wave*sexF+wave*african, proband.1~adhd.1, ClusInd.name="id", IntraClus.name="wave", 
yname="adhd", ystrata="probandS", xstrata="sex", data=adhd, NMat=NMat, nzval0=40, control=mlefn.control(niter=100))

summary(z)
}
\keyword{Statistics|cluster|regression}

