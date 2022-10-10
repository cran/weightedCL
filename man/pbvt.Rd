\name{pbvt}
\alias{pbvt}
\title{
BIVARIATE NORMAL AND STUDENT CDFs WITH VECTORIZED INPUTS  }
\description{
Bivariate normal and Student cdfs with vectorized inputs
}
\usage{
pbvt(z1,z2,param,icheck=FALSE)
}
\arguments{
\item{z1}{scalar or vector of reals }
\item{z2}{scalar or vector of reals }
\item{param}{vector of length 2, or matrix with 2 columns;
vectors and number of rows of matrix cannot be different if larger than 1;
for param, first column is rho, second column is df.}
\item{icheck}{TRUE if checks are made for proper inputs, default of FALSE}
}
\value{
cdf value(s)
}

\references{
Joe H (2014) \emph{CopulaModel:  Dependence Modeling with Copulas}. Software for book: \emph{Dependence Modeling with Copulas}, Chapman & Hall/CRC, 2014.
}

\examples{
cat("\n pbvt rho changing\n")
z1=.3; z2=.4; rho=seq(-.9,.9,.1); nu=2
param=cbind(rho,rep(nu,length(rho)))
out1=pbvt(z1,z2,param)
print(cbind(rho,out1))
cat("\n pbvt z1 changing\n")
z1=seq(-2,2,.4)
z2=.4; rho=.5; nu=2
out2=pbvt(z1,z2,c(rho,nu))
print(cbind(z1,out2))
}

\keyword{distribution}
