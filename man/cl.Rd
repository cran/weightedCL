\name{cl}
\alias{cl}
\alias{cl.ord}
\title{COMPOSITE LIKELIHOOD ESTIMATION FOR MVN COPULA}
\description{
Composite likelihood estimation for MVN copula.
}
\usage{
cl(p,q,b,gam,xdat,ydat,margmodel,link)
cl.ord(p,q,b,gam,xdat,ydat,link)
}
\arguments{
  \item{p}{The order of the autoregressive component.}
  \item{q}{The order of the moving average component.}
  \item{b}{The regression coefficients.}
  \item{gam}{The uinivariate parameters that are not regression coefficients. That is the parameter \eqn{\gamma} of negative binomial distribution or  the \eqn{q}-dimensional vector of the univariate cutpoints of ordinal model. \eqn{\gamma} is NULL for Poisson and binary regression.}
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
  the probit link function. 
}
}
\details{The
 composite likelihood method in Zhao and Joe (2005). The
univariate parameters  are estimated
from the sum of univariate marginal log-likelihoods and then the dependence
parameters are estimated from the sum of bivariate marginal log-likelihoods
with the univariate parameters fixed from the first step. 

Note that \code{cl.ord} is a variant of the code for ordinal (probit and logistic) regression.
}
\value{
A list containing the following components:
 \item{minimum}{The negative value of the sum of bivariate marginal log-likelihoods at
CL1 estimates.}
 \item{estimate}{The composite likelihood estimates.}
 \item{gradient}{The gradient at the estimated minimum of CL1.}
 \item{code}{An integer indicating why the optimization process terminated,
 same as in \code{\link{nlm}}.}
}
\references{
Zhao, Y. and Joe, H. (2005)
Composite likelihood estimation in multivariate data analysis.
\emph{The Canadian Journal of Statistics},  \bold{33}, 335--356.
}
\author{
Aristidis K. Nikoloulopoulos \email{A.Nikoloulopoulos@uea.ac.uk}
}
\seealso{
  \code{\link{wcl}}
  \code{\link{iee}}
}
\examples{
\donttest{
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
################################################################################
#                      select the  ARMA structure
################################################################################
p=2;q=1
################################################################################
#                      perform CL1 estimation
################################################################################
i.est<-iee(xdat,ydat,margmodel)
cat("\niest: IEE estimates\n")
print(c(i.est$reg,i.est$gam))
est.rho<-cl(p=p,q=q,b=i.est$reg,gam=i.est$gam,
             xdat,ydat,margmodel,link)
cat("\nest.rho: CL estimates\n")
print(est.rho$e)
################################################################################
#                         Ordinal time-series regression 
################################################################################
################################################################################
#                      read and set up data set
################################################################################
data(sleep)
ydat=sleep$sleep
bydat=oydat=ydat
bydat[ydat==4]=0
bydat[ydat<4]=1
oydat[ydat==4]=1
oydat[ydat<4]=2
oydat[ydat==2]=3
oydat[ydat==3]=4

x1=sleep$heartrate
x2=sleep$temperature
z1=(x1-mean(x1))/sd(x1)
z2=(x2-mean(x2))/sd(x2)
xdat=cbind(z1,z2)
################################################################################
#                      select the link
################################################################################
link="probit"
################################################################################
#                      select the ARMA structure
################################################################################
p=1;q=0
################################################################################
#                      perform CL1 estimation
################################################################################
i.est<-iee.ord(xdat,oydat,link)
cat("\niest: IEE estimates\n")
print(c(i.est$reg,i.est$gam))
est.rho<-cl.ord(p=p,q=q,b=i.est$reg,gam=i.est$gam,
             xdat,oydat,link)
cat("\nest.rho: CL estimates\n")
print(est.rho$e)
}
}


\keyword{regression}

