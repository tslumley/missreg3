\name{rclusbin}
\alias{rclusbin}
\docType{data}
\title{
Random intercept model for clustered binary data
}
\description{
Fits random intercept models to clustered binary data with the two-phase missingness structure.
}
\usage{
rclusbin(formula, data, weights=NULL, ClusInd=NULL, IntraClus=NULL, xstrata=NULL, ystrata=NULL, 
obstype.name="obstype", NMat=NULL, xs.includes=FALSE, MaxInClus=NULL, rmsingletons=FALSE, retrosamp="proband", 
gamma=NULL, nzval0=20, fit=TRUE, devcheck=FALSE, linkname="logit", start=NULL, Qstart=NULL, sigma=NULL, 
control=mlefn.control(...), control.inner=mlefn.control.inner(...), ...)
}
\arguments{
  \item{formula}{A symbolic description of the model to be fitted.}
  \item{data}{A data frame containing all the variables required for analysis, including those for \code{xstrata}, \code{ystrata} and 
		\code{obstype.name}.}
  \item{weights}{An optional vector of weights to be used in the fitting process. Should be \code{NULL} or the name of a numeric vector in the data frame. It provides weights at the individual level when there are clusters of size greater than one. \cr When all clusters are of size one, it provides weights at cluster=individual level.}
  \item{ClusInd}{Name of a vector in the data frame specifying cluster membership. Can be \code{NULL} if all clusters are of size one.}
  \item{IntraClus}{Name of a vector in the data frame specifying intra-cluster sequence of individual subjects in a cluster. The one with the smallest i.d. is treated as the proband who were originally sampled into a study.}
  \item{xstrata}{Specify names of the stratification variables to be used,e.g. \code{"vname"} or \code{c("vname1","vname2",...)}. Strata are defined by cross-classifiction of all levels. \cr This function only deals with the situation when clusters are defined within xstrata.}
  \item{ystrata}{Specify name of the variable defing the Y-strata. Compulsory when gamma probabilities are used (see Details for more 
		descriptions).}
  \item{obstype.name}{Name of the variable specifying labels for observations by sampling and variable type: \code{"uncond"}, \code{"retro"}, \code{"xonly"}, \code{"y|x"} or \code{"strata"}. }
  \item{NMat}{Population counts in a matrix form with rows and columns corresponding to Y-strata and X-strata respectively. Should not be provided when there is any observation of the type \code{"strata"}.}
  \item{xs.includes}{\code{TRUE} if \code{weights} specified for obervations labelled as "strata" include those observed at the second phase (i.e. \code{"retro"} or \code{"uncond"} observations).}
  \item{MaxInClus}{A value specifying the maximum number of individuals allowed in a cluster. Set to \code{NULL} if there is no limit.}
  \item{rmsingletons}{If \code{TRUE}, remove clusters of size one.}
  \item{retrosamp}{Three restrospective sampling schemes can be applied based on the Y-status of all subjects in the same cluster: \code{"proband"}, \code{"allcontrol"} and \code{"gamma"} (see Details for more descriptions). The default is \code{"proband"}.}
  \item{gamma}{A vector of length 2 specifying the probabilities that individuals belong to \code{Y=1} based on their cluster status (see Details for more descriptions). }
  \item{nzval0}{Number of points to calculate the zeros and weights needed for Gauss-Hermite quadrature. The default is 20.}
  \item{fit}{If \code{FALSE}, only stratum report will be generated without model fitting.}
  \item{devcheck}{If \code{TRUE}, check the first and second derivatives. The default should be \code{FALSE}.}
  \item{linkname}{A specification for the model link function. Three choices are provide: \code{"logit"}, \code{"probit"} or "cloglog". The default is \code{"logit"}. }
  \item{start}{Starting values for the regression parameters. The first \code{p}-coefficients are parameters for X-variables. The last parameter is for the random intercept term and normally denoted as \code{w=log(sigma)}. \cr The program cannot provide starting values for all data strctures so will force you to use this whenever it is necessary. }
  \item{Qstart}{An optional starting matrix for Pr(Ystratum=i|Xstratum=j). Can be compulsory if the program cannot produce a valid starting value at some situations. }
  \item{sigma}{An optional starting value for \code{sigma}. The default value (when set to \code{NULL}) is 0.5.}
  \item{control}{Specify control parameters for the iterations in \code{mlefn} call.}
  \item{control.inner}{Specify control parameters for inner iterations nested within \code{mlefn} call. }
  \item{\dots}{Further arguments passed to or from related function. }
}
\details{
This function fits binary regression models with a random intercept of the form \code{a_i=e^{w*eps_i}} where \code{w=log(sigma)} and \code{eps_i} is standard normal for each cluster, along with a linear predictor \code{eta_{ij}=x_{ij}^T*beta} for the subject \code{j} in the \code{i^th} cluster. 
\cr \cr
 The function can be applied to both prospective and retrospective data with various types of observations collected at different two-phase sampling schemes. Three retrospective samplings are considered with the Y-strata defined as: \cr
(1) the case-control status of the proband only (\code{"proband"}); 
(2) the case-cotnrol status of all members in the same cluster (\code{"allcontrol"}). If any one of the members are cases, the cluster belongs to Y-strata=1 and otherwise Y-strata=0; \cr
(3) the case-control status of all members in the same cluster plus the gamma probabilities (\code{"gamma"}). The conditional probability of Y-strata=1 depends on \code{sum_j{Y_j}=1} (with \code{gamma_1} probability) or \code{sum_j{Y_j}>1} (with \code{gamma_2} probability). \cr Here \code{Y_j} indicates case-control status (1 for a case and 0 for a control) of the \code{j^{th}} individual in a cluster.
}
\source{
http://www.stat.auckland.ac.nz/~wild
}
\references{
Description of the \code{missreg} Library, Wild and Jiang, 2007.
}
\seealso{ \code{\link[glmmML]{ghq}},  \code{\link{rclusbin2}}}
\examples{
data(brainpairs)
brainpairs$obstype <- rep("retro", dim(brainpairs)[1])
z2 <- rclusbin(bt ~ ep + ca, ClusInd="id", IntraClus="relid", data=brainpairs)
summary(z2)

data(rdat00)
z3 <- rclusbin(y~x, ClusInd="cluster", data=rdat00, retrosamp="allcontrol")
summary(z3)
}
\keyword{Statistics|cluster|regression}

