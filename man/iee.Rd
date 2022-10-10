\name{iee}
\alias{iee}
\title{INDEPENDENT ESTIMATING EQUATIONS FOR BINARY AND COUNT REGRESSION}
\description{
Independent estimating equations for binary and count regression.}
\usage{
iee(xdat,ydat,margmodel,link) 
}
\arguments{
  \item{xdat}{The \eqn{d\times p} matrix of covariates, where \eqn{d} is the length of the time-series and \eqn{p} is the number of covariates including
  the unit first column  to account for the intercept (except for ordinal regression where there is no intercept).}
  \item{ydat}{The \eqn{d}-dimensional vector of  dicrete time series reponse, where \eqn{d} is the length of the series.}
  \item{margmodel}{
  Indicates the marginal model.
  Choices are \dQuote{poisson} for Poisson, \dQuote{bernoulli} for Bernoulli, 
  and  \dQuote{nb1} , \dQuote{nb2} for the NB1 and NB2 parametrization 
  of negative binomial in Cameron and Trivedi (1998).}
  \item{link}{The link function.
  Choices are \dQuote{log} for the log link function, \dQuote{logit} for 
  the logit link function, and \dQuote{probit} for 
  the probit link function.}
  }
\details{The univariate parameters  are
estimated from the sum of univariate marginal log-likelihoods.

}
\value{A list containing the following components:
\item{coef}{The vector with the estimated regression parameters.}
\item{gam}{The vector with the estimated parameters that  are not regression parameters. This is NULL for Poisson and binary regression.} }
\references{
Cameron, A. C. and Trivedi, P. K. (1998)
\emph{Regression Analysis of Count Data.}
Cambridge: Cambridge University Press.
}
\author{
Aristidis K. Nikoloulopoulos \email{A.Nikoloulopoulos@uea.ac.uk}
}



\examples{
################################################################################
#                      NB2 regression for count time-series data
################################################################################
################################################################################
#                      read and set up data set
################################################################################
data(polio)
ydat <-polio
d=length(ydat)
tvec=1:length(ydat)
tvec1=tvec-73
xdat <- cbind(1, tvec1/1000, cos(2 * pi * tvec1 / 12), sin(2 * pi * tvec1 / 12),
                      cos(2 * pi * tvec1 / 6), sin(2 * pi * tvec1 / 6)) 
################################################################################
#                      select the marginal model
################################################################################
margmodel="nb2"

i.est<-iee(xdat,ydat,margmodel)
cat("\niest: IEE estimates\n")
print(c(i.est$reg,i.est$gam))
}
\keyword{regression}



