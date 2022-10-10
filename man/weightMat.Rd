\name{weightMat}
\alias{weightMat}
\alias{weightMat.ord}
\title{WEIGHT MATRICES FOR THE WEIGHTED COMPOSITE LIKELIHOOD ESTIMATING EQUATIONS}
\description{Weight matrices for the weighted composite likelhood estimating equations.}
\usage{
weightMat(b,gam,rh,p,q,xdat,margmodel,link)
weightMat.ord(b,gam,rh,p,q,xdat,link)
}
\arguments{
 \item{b}{The regression coefficients.}
  \item{gam}{The uinivariate parameters that are not regression coefficients. That is the parameter \eqn{\gamma} of negative binomial distribution or  the \eqn{q}-dimensional vector of the univariate cutpoints of ordinal model. \eqn{\gamma} is NULL for Poisson and binary regression.}
\item{rh}{The vector of autregressive and moving average parameters in high-dimensional normal copula regression models with an ARMA(\eqn{p,q}) correlation matrix.}
 \item{p}{The order of the autoregressive component.}
  \item{q}{The order of the moving average component.}
  \item{xdat}{The \eqn{d\times p} matrix of covariates, where \eqn{d} is the length of the time-series and \eqn{p} is the number of covariates including
  the unit first column  to account for the intercept (except for ordinal regression where there is no intercept).}
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
\details{The  matrices that form the weight matrices \eqn{\mathbf{W}^{(1)}}  of the weighted composite likelihood estimating equations in Nikoloulopoulos et al. (2022).
Note that \code{weightMat.ord} is a variant of the code for ordinal (probit and logistic) regression.
}
\value{A list containing the following components:
\item{omega}{ The  \eqn{\boldsymbol{\Omega}^{(1)}} matrix.}
\item{delta}{ The  \eqn{\boldsymbol{\Delta}^{(1)}} matrix.}
\item{X}{ The \eqn{\mathbf{X}} matrix.}
}
\seealso{
  \code{\link{wcl}},
  \code{\link{godambe}},
}
\references{
Nikoloulopoulos, A.K. (2022)
Efficient and feasible inference for high-dimensional normal copula regression models.
 \emph{Arxiv e-prints}, <arXiv:2203.04619>. \url{https://arxiv.org/abs/2203.04619}.
}
\author{
Aristidis K. Nikoloulopoulos \email{A.Nikoloulopoulos@uea.ac.uk}
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
cat("\nest.rho: CL1 estimates\n")
print(est.rho$e)
WtScMat<-weightMat.ord(b=i.est$reg,gam=i.est$gam,rh=est.rho$e,xdat,ydat,id,tvec,corstr,link)
################################################################################
#                      obtain the weight matrices
################################################################################
WtScMat<-weightMat.ord(b=i.est$reg,gam=i.est$gam,rh=est.rho$e,
                   p=p,q=q,xdat,link)
}
}




