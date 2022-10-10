\name{wcl}
\alias{wcl}
\alias{wcl.ord}
\title{SOLVING THE WEIGHTED COMPOSITE LIKELIHOOD ESTIMATING EQUATIONS WITH INPUTS THE WEIGHT MATRICES AND  DATA}
\description{Solving the weighted composite likelihood estimating equations with inputs the weight matrices
  and data.}
\usage{
wcl(start,WtScMat,xdat,ydat,margmodel,link)
wcl.ord(start,WtScMat,xdat,ydat,link)
}
\arguments{
 \item{start}{A starting value of the vector of regression
 and not regression parameters. The composite likelihood estimates of regression and not regression parameters 
 is a good starting value.}
 \item{WtScMat}{A list containing the following components.
  omega: the matrix  \eqn{\boldsymbol{\Omega}^{(1)}};
  delta: the matrix \eqn{\boldsymbol{\Delta}^{(1)}};
  X: the matrix \eqn{\mathbf{X}}.}
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
\details{Obtain  estimates \eqn{{\hat{\mathbf{a}}}} of the univariate parameters
solving the weighted composite likelihood estimating equations.

Note that \code{wcl.ord} is a variant of the code for ordinal (probit and logistic) regression.
}
\value{A list containing the following components:
\item{root}{The weighted composite likelihood estimates.}
\item{f.root}{The value of the weighted composite likelihood estimating equations evaluated at the root.}
\item{iter}{The number of iterations used.}
\item{estim.precis}{The estimated precision for root.}
}
\references{
Nikoloulopoulos, A.K. (2022)
Efficient and feasible inference for high-dimensional normal copula regression models.
 \emph{Arxiv e-prints}, <arXiv:2203.04619>. \url{https://arxiv.org/abs/2203.04619}.
}
\author{
Aristidis K. Nikoloulopoulos \email{A.Nikoloulopoulos@uea.ac.uk}
}
\seealso{
  \code{\link{weightMat}},
  \code{\link{godambe}}
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
#                      obtain the weight matrices
################################################################################
WtScMat<-weightMat(b=i.est$reg,gam=i.est$gam,rh=est.rho$e,
                   p=p,q=q,xdat,margmodel)
################################################################################
#           obtain the weighted composite likelihood estimates
################################################################################  
est<-wcl(start=c(i.est$reg,i.est$gam),WtScMat,xdat,ydat,
              margmodel,link)
cat("est=parameter estimates\n")
print(est$r)                                  
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
################################################################################
#                      obtain the weight matrices
################################################################################
WtScMat<-weightMat.ord(b=i.est$reg,gam=i.est$gam,rh=est.rho$e,
                   p=p,q=q,xdat,link)
################################################################################
#           obtain the weighted composite likelihood estimates
################################################################################  
est<-wcl.ord(start=c(i.est$reg,i.est$gam),WtScMat,
                  xdat,oydat,link)
cat("est=parameter estimates\n")
print(est$r)               
}
}

\keyword{regression}
