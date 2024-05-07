\name{ProjectionCPE}
\title{Projecton-based Approach to Estimate the CPE from Nested Survival Models.}
\alias{ProjectionCPE}
\description{
Using a projection-based approach, the concordance probability estimate (CPE) for nested Cox proportional hazards regression models is computed. Under this approach, the proportional hazards assumption is satisfied for the model containing both the standard covariates and the new covariates under investigation. 
}
\usage{
ProjectionCPE(Time, Event,StandardMarkers,NewMarkers, tau, Block=TRUE)
 
}
\arguments{
  \item{Time}{Survival Time}
  \item{Event}{Event Status}
  \item{StandardMarkers}{A vector or matrix of standard covariates.}
  \item{NewMarkers}{A vector or matrix of new covariates.}
  \item{tau}{The follow-up duration for a given study. This should be no larger than the maximum observed failure time (uncensored).}
  \item{Block}{Due to the complexity of the calculation, this is an approximation that calculates the CPE using random blocks of 50 patients. Recommended for data sets larger than 150 individuals.} 
}

\value{an object with the following elements
  \item{projCPE}{test statistic using the projection approach.}
  }

\details{
  This function implements the method described in Heller and Devlin (under review).
}

\examples{

SM0 = rnorm(300, 0,1)
SM1 = rnorm(300, 0,1)
NM0 = rnorm(300, 0,1)
NM1 = rnorm(300, 0,1)
  
S1 = exp(-1*(0.408*SM0 +0.15*SM1 + 0.684*NM0+0.15*NM1))*rweibull(300,scale = 1 , shape=1) 
C1 = runif(300, 0,1.58)
  
Time = pmin(S1, C1)
Event = 1*(S1 < C1)
 
StandardMarkers = cbind(SM0,SM1)
NewMarkers = cbind(NM0 ,NM1 )
tau = max(Time[Event==1])

ProjectionCPE(Time, Event,StandardMarkers,NewMarkers,  tau)

}
\references{
  Heller, Glenn, and Sean M. Devlin. Measuring the Impact of New Risk Factors Within Survival Models. \emph{Under Review}.
}
