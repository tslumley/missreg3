\name{leprosy1}
\alias{leprosy1}
\docType{data}
\title{The Leprosy data}
\description{
  The leprosy data set was described in Scott and Wild (1997, 2001) and 
used as an example of standard two-phase case-control sampled data.  
}
\usage{data(leprosy1)}
\format{
  A data frame with 42 observations on the following 5 variables.
  \describe{
    \item{\code{leprosy}}{a factor with levels \code{no} \code{yes}}
    \item{\code{age}}{a numeric vector indicating the mid-point of six 5-year age groups}
    \item{\code{scar}}{a factor with levels \code{no} \code{yes}}
    \item{\code{counts}}{a numeric vector indicating number of subjects with each observation}
    \item{\code{obstype}}{a factor with levels \code{retro} and \code{strata}, compulsory for function call}
  }
}
\details{
  The leprosy data were obtained by sampling from the results of a population
cross-sectional study of people under 35 in Northern Malawi and represented as 
a three-way contingency table in Clayton and Hills (1993). Those people with leprosy
were defined to be cases and the rest to be controls. The data were first categorized
into six 5-year age sampling strata and the numbers of cases and controls falling into
each stratum were observed. All cases have been chosen with an equal-sized control
group subsampled from the control population within each age stratum. The potential 
risk factor that indicates the presence or absence of a BCG vaccination scar was then
observed.
\cr \cr
The data are represented in three different formats in \code{leprosy1}, 
\code{leprosy2} and \code{leprosy3}. See "Description of the \code{missreg}
Library" for more details.   
}
\references{
Description of the \code{missreg} Library, Wild and Jiang, 2007
}
\examples{
data(leprosy1)
}
\keyword{datasets}
