\name{lowbirth.bin}
\alias{lowbirth.bin}
\non_function{}
\title{The Low Birthweight data}
\usage{data(lowbirth.bin)}
\description{
A subset of the data collected in the Auckland Birthweight Collaborative (ABC) Study.
}
\format{
  A data frame with 1148 observations on the following 18 variables.
  \describe{
    \item{sgagp}{A factor defining the case (sga) and control (aga) status of 
the baby.}
    \item{sex}{1=female or 2=male}
    \item{instudy}{A factor with levels \code{retro} and \code{strata}, 
used as the \code{obstype} variable in function call}
    \item{htstrat}{A factor with levels as class intervals of mother's height}
    \item{ethnicdb}{A factor with levels \code{A} (Asian), \code{E} (Euro.), 
\code{M} (Maori) or \code{P} (Pacifican)}
    \item{mstrat}{Marital status of the mother}
    \item{occ}{Mother's ccupational group, 1 to 3 (3 is highest)}
    \item{mumht}{Height of the mother in cm}
    \item{mumwt}{Weight of the mother in kg}
    \item{bmi}{Body mass index of the mother}
    \item{smoke}{Smoking status of the mother prior to pregnancy}
    \item{smokedb}{Smoking variable from database}
    \item{age1st}{Mother's age at first pregnancy}
    \item{hyper}{Any hypertension (1=yes, 0=no)}
    \item{edstratdb}{Mother's educational level}
    \item{eductrm}{Mother's age left school}
    \item{mstratdb}{As for \code{mstrat} with some levels combined}
    \item{counts}{Number of subjects with each observation (frequency)}
  }
}
\references{
Description of the \code{missreg} Library, Wild and Jiang, 2007.
}
\details{
The ABC study was conducted in 1995-1997 in order to find potential risk factors 
associated with small-for-gestational-age babies in New Zealand. It was a 
population-based case-control study with the cases being those babies with
their birthweights equal to or below the sex-specific 10th percentile for 
gestational age in the population.
\cr \cr 
The \code{lowbirh.bin} is a semi-random subset of the original data.}
\keyword{datasets}






