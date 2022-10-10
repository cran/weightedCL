\name{margmodel}
\alias{margmodel.ord}
\alias{dmargmodel}
\alias{pmargmodel}
\alias{dmargmodel.ord}
\alias{pmargmodel.ord}
\title{DENSITY AND CDF OF THE UNIVARIATE MARGINAL DISTRIBUTION}
\description{Density and cdf of the univariate marginal distribution.}
\usage{
dmargmodel(y,mu,gam,invgam,margmodel)
pmargmodel(y,mu,gam,invgam,margmodel)
dmargmodel.ord(y,mu,gam,link)
pmargmodel.ord(y,mu,gam,link)
}
\arguments{
  \item{y}{Vector of (non-negative integer) quantiles.}
  \item{mu}{The parameter \eqn{\mu} of the univariate distribution.}
  \item{gam}{The parameter(s) \eqn{\gamma} that are not regression parameters. \eqn{\gamma} is NULL for Poisson and Bernoulli distribution.}
  \item{invgam}{The  inverse of parameter \eqn{\gamma} of negative binomial
  distribution.}
  \item{margmodel}{
  Indicates the marginal model.
  Choices are \dQuote{poisson} for Poisson, \dQuote{bernoulli} for Bernoulli, 
  and  \dQuote{nb1} , \dQuote{nb2} for the NB1 and NB2 parametrization 
  of negative binomial in Cameron and Trivedi (1998). See details.}
  \item{link}{The link function. 
  Choices are  \dQuote{logit} for 
  the logit link function, and \dQuote{probit} for 
  the probit link function. }
  }
\details{Negative binomial distribution
NB\eqn{(\tau,\xi)} allows for overdispersion
and its probability mass function (pmf) is  given by
  \deqn{ f(y;\tau,\xi)=\frac{\Gamma(\tau+y)}{\Gamma(\tau)\; y!}
  \frac{\xi^y}{(1+\xi)^{\tau + y}},\quad \begin{matrix} y=0,1,2,\ldots, \\
                                          \tau>0,\; \xi>0,\end{matrix}  }
with mean \eqn{\mu=\tau\,\xi=\exp(\beta^T x)} and variance \eqn{\tau\,\xi\,(1+\xi)}.


Cameron and Trivedi (1998) present the NBk parametrization where
\eqn{\tau=\mu^{2-k}\gamma^{-1}} and \eqn{\xi=\mu^{k-1}\gamma}, \eqn{1\le k\le 2}.
In this function we use the NB1 parametrization
\eqn{(\tau=\mu\gamma^{-1},\; \xi=\gamma)}, and the NB2 parametrization
\eqn{(\tau=\gamma^{-1},\; \xi=\mu\gamma)}; the latter
is the same as in Lawless (1987).

\code{margmodel.ord} is a variant of the code for ordinal (probit and logistic) model. In this case, the response \eqn{Y} is assumed to have density
\deqn{f_1(y;\nu,\gamma)=F(\alpha_{y}+\nu)-F(\alpha_{y-1}+\nu),}
where \eqn{\nu=x\beta} is a function of \eqn{x}
and the \eqn{p}-dimensional regression vector \eqn{\beta}, and \eqn{\gamma=(\alpha_1,\ldots,\alpha_{K-1})} is the $q$-dimensional vector of the univariate cutpoints (\eqn{q=K-1}). Note that \eqn{F} normal leads to the probit model and \eqn{F} logistic
leads to the cumulative logit model for ordinal response.

}
\value{The density and cdf of the univariate distribution.}
\references{
Cameron, A. C. and Trivedi, P. K. (1998)
\emph{Regression Analysis of Count Data.}
Cambridge: Cambridge University Press.

Lawless, J. F. (1987)
Negative binomial and mixed Poisson regression.
\emph{The Canadian Journal
of Statistics}, \bold{15}, 209--225.


}

\examples{
y<-3
gam<-2.5
invgam<-1/2.5
mu<-0.5
margmodel<-"nb2"
dmargmodel(y,mu,gam,invgam,margmodel)
pmargmodel(y,mu,gam,invgam,margmodel)
link="probit"
dmargmodel.ord(y,mu,gam,link)
pmargmodel.ord(y,mu,gam,link)

}
\keyword{distribution}
