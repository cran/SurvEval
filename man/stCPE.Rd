\name{stCPE}
\title{Estimate the Durability of a Baseline Risk Score or the CPE for a Fixed Follow-up Duration.}
\alias{stCPE}
\description{
This implements the methodology developed by Devlin, Gonen, and Heller (2020) to estimate the durability of a baseline risk score estimated under a Cox proportional hazards model. The same methodology can be used to estimate the concordance probability estimate (CPE) up to a fixed follow-up duration in a study. 
}
\usage{
stCPE(Time, Event, Markers, starttime, tau)
 
}
\arguments{
  \item{Time}{Survival Time.}
  \item{Event}{Event Status.}
  \item{Markers}{A vector or matrix of covariates.}
  \item{starttime}{The start time when evaluating the durability of a risk score. This should be set to 0 when evaluating the CPE up to a study's follow-up duration.}
  \item{tau}{The end time under evaluation.}
}

\value{an object with the following element
  \item{stCPE.estimate}{which is the concordance probability estimate.}
  }

\details{
  This function implements the method described in Devlin, Gonen, and Heller (2020).
}

\examples{

SM0 = rnorm(300, 0,1)
SM1 = rnorm(300, 0,1)
NM0 = rnorm(300, 0,1)
NM1 = rnorm(300, 0,1)
  
S1 <- exp(-1*(0.408*SM0 +0.15*SM1 + 0.684*NM0+0.15*NM1))*rweibull(300,scale = 1 , shape=1) 
C1 <- runif(300, 0,1.58)
  
Time <- pmin(S1, C1)
Event <- 1*(S1 < C1)
 
Markers <- cbind(SM0,SM1,NM0 ,NM1 )
tau <- max(Time[Event==1])

stCPE(Time, Event, Markers, 0, tau)

}
\references{
Devlin, Sean M., Mithat Gonen, and Glenn Heller. Measuring the temporal prognostic utility of a baseline risk score. Lifetime data analysis (2020).
}
