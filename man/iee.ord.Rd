\name{iee.ord}
\alias{iee.ord}
\title{
Maximum Likelihood for Ordinal Model 
}
\description{
Maximum Likelihood for Ordinal Probit and Logit: 
Newton-Raphson minimization of negative log-likelihood. 
}
\usage{
iee.ord(x,y,link,iprint=0,maxiter=20,toler=1.e-6)
}
\arguments{
\item{x}{
   vector or matrix of explanatory variables. Each row corresponds
   to an observation and each column to a variable.  The number of rows of x
   should equal the number of data values in y, and there should be fewer 
   columns than rows.  Missing values are not allowed.
}
\item{y}{
   numeric vector containing the ordinal response. 
   The values must be in the range 1,2,..., number of categories.
   Missing values are not allowed. 
}
\item{link}{The link function.Choices are  \dQuote{logit} for 
  the logit link function, and \dQuote{probit} for 
  the probit link function.}
\item{iprint}{
  logical indicator, default is FALSE, for whether the iterations for 
  numerical maximum likelihood should be printed.
}
\item{maxiter}{
  maximum number of Newton-Raphson iterations, default = 20.
}
\item{toler}{
  tolerance for convergence in Newton-Raphson iterations, default = 1.e-6.
}
}
\value{
  list of MLE of parameters and their associated standard 
  errors, in the order cutpt1,...,cutpt(number of categ-1),b1,...b(number of covariates).
\item{negloglik}{
  value of negative log-likelihood, evaluated at MLE
}
\item{gam}{
  MLE of ordered cutpoint parameters
}
\item{reg}{
  MLE of regression parameters
}
\item{cov}{
  estimated covariance matrix of the parameters
}
}
\details{
 
  The ordinal probit model is similar to the ordinal logit model.
  The parameter estimate of ordinal logit are roughly 1.8 to 2 times
  those of ordinal probit. 
}
\references{
  Anderson, J.A. and Pemberton, J.D. (1985). The grouped continuous model for
    multivariate ordered categorical variables and covariate adjustment.  
    \emph{Biometrics},  \bold{41}, 875--885.

}
\examples{
################################################################################
#                         Ordinal regression 
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
i.est<-iee.ord(xdat,ydat,link)
cat("\niest: IEE estimates\n")
print(c(i.est$reg,i.est$gam))
}

\keyword{regression}
