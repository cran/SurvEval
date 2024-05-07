\name{SurvEval-package}
\alias{SurvEval-package}
\alias{SurvEval}
\docType{package}
\title{
  Methods for the Evaluation of Survival Models.
}
\description{
This implements various tools to evaluate Cox proportional hazards regression models. 
}
\details{
Currently, this implements the methodology outlined in the following two papers. The first is Devlin, Gonen, and Heller (2020) that evaluates the discriminatory power of a baseline-derived risk score over time. This same method is used to evaluate the concordance probability estimate (CPE) up until the follow-up duration of a given study. 

The second is Heller and Devlin (under review) that develops a projection-based approach to estimate the CPE from nested survival models. 


}
\author{
Sean Devlin and Glenn Heller.

Maintainer: Sean Devlin <devlins@mskcc.org>
}
\references{
Devlin, Sean M., Mithat Gonen, and Glenn Heller. Measuring the temporal prognostic utility of a baseline risk score. Lifetime data analysis 26.4 (2020): 856-871.

Heller, Glenn, and Sean M. Devlin. Measuring the Impact of New Risk Factors Within Survival Models. Under Review.  
}
\keyword{ package }
\keyword{ CPE }
\keyword{ Survival }
