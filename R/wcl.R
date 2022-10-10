

# the mean values of the univariate marginal distribution 
# corresonding to the used link function 
# input:
# x the matix of the covariates 
# b the vector with the regression coefficients
# link has three options: 1. "log", 2. "logit". 3. "probit"
# output:
# the mean values of the univariate marginal distribution 
linked.mu<-function(x,b,link)
{ if(link=="log")
    { mu<-exp(x %*% b)
    }
    else
    { if(link=="logit")
      { expnu<-exp(x %*% b)
        mu<-expnu/(1+expnu)
      }
    else
    { # link=probit
      mu<-pnorm(x %*% b) }
    }
   mu
}

# Density  of the univariate marginal distribution
# input:
# y the vector of (non-negative integer) quantiles.
# mu the mean parameter of the univariate distribution.
# gam the parameter gamma of  the negative binomial distribution.
# invgam the inverse of parameter  gamma of negative binomial distribution.
# margmodel indicates the marginal model. Choices are ?poisson? for Poisson, ?bernoulli? for
# Bernoulli, and ?nb1? , ?nb2? for the NB1 and NB2 parametrization of negative
# binomial in Cameron and Trivedi (1998).
# output:
# the density of the univariate marginal  distribution
dmargmodel<-function(y,mu,gam,invgam,margmodel)
{ if(margmodel=="poisson")
    { dpois(y,mu)
    }
    else
    { if(margmodel=="bernoulli")
      { dbinom(y,size=1,prob=mu) }
    else
    { if(margmodel=="nb1")
      { dnbinom(y,prob=1/(1+gam),size=mu*invgam) }
    else
    { # margmodel=="nb2"
      dnbinom(y,size=invgam,mu=mu) }}}
}

# CDF  of the univariate marginal distribution
# input:
# y the vector of (non-negative integer) quantiles.
# mu the mean parameter of the univariate distribution.
# gam the parameter gamma of  the negative binomial distribution.
# invgam the inverse of parameter  gamma of negative binomial distribution.
# margmodel indicates the marginal model. Choices are ?poisson? for Poisson, ?bernoulli? for
# Bernoulli, and ?nb1? , ?nb2? for the NB1 and NB2 parametrization of negative
# binomial in Cameron and Trivedi (1998).
# output:
# the cdf of the univariate marginal  distribution
pmargmodel<-function(y,mu,gam,invgam,margmodel)
{ if(margmodel=="poisson")
    { ppois(y,mu)
    }
    else
    { if(margmodel=="bernoulli")
      { pbinom(y,size=1,prob=mu) }
    else
    { if(margmodel=="nb1")
      { pnbinom(y,prob=1/(1+gam),size=mu*invgam) }
    else
    { # margmodel=="nb2"
      pnbinom(y,size=invgam,mu=mu) }}}
}

# quantile  of the univariate marginal distribution
# input:
# y the vector of probabilities
# mu the mean parameter of the univariate distribution.
# gam the parameter gamma of  the negative binomial distribution.
# invgam the inverse of parameter  gamma of negative binomial distribution.
# margmodel indicates the marginal model. Choices are ?poisson? for Poisson, ?bernoulli? for
# Bernoulli, and ?nb1? , ?nb2? for the NB1 and NB2 parametrization of negative
# binomial in Cameron and Trivedi (1998).
# output:
# the quantile of the univariate marginal  distribution
qmargmodel<-function(y,mu,gam,invgam,margmodel)
{ if(margmodel=="poisson")
    { qpois(y,mu)
    }
    else
    { if(margmodel=="bernoulli")
      { qbinom(y,size=1,prob=mu) }
    else
    { if(margmodel=="nb1")
      { qnbinom(y,prob=1/(1+gam),size=mu*invgam) }
    else
    { # margmodel=="nb2"
      qnbinom(y,size=invgam,mu=mu) }}}
}

# negative univariate logikelihood assuming independence within clusters
# input:
# param the vector of regression and not regression parameters
# xdat the matrix of covariates (use the constant 1 for the first covariate)
# ydat the vector with the response
# margmodel indicates the marginal model. Choices are ?poisson? for Poisson, 
# ?bernoulli? for Bernoulli, and ?nb1? , ?nb2? for the NB1 and NB2 
# parametrization of negative binomial in Cameron and Trivedi (1998)
# link is the link function. Choices are ?log? for the log link function, 
# ?logit? for the logit link function, and ?probit? for the probit link function.
# output:
# negative univariate logikelihood assuming independence within clusters
marglik<-function(param,xdat,ydat,margmodel,link)
{ p<-dim(xdat)[2]
  b<-param[1:p]
  if(margmodel=="nb1" | margmodel=="nb2")
  { gam<-param[p+1]
    invgam<-1/gam
  }
  #else  gam<-invgam<-0
  mu<-linked.mu(as.matrix(xdat),b,link)
  -sum(log(dmargmodel(ydat,mu,gam,invgam,margmodel)))
}
 
 
# Independent estimating equations for binary, Poisson or 
# negative binomial regression.
# xdat the matrix of covariates (use the constant 1 for the first covariate)
# ydat the vector with the response
# margmodel indicates the marginal model. Choices are ?poisson? for Poisson, 
# ?bernoulli? for Bernoulli, and ?nb1? , ?nb2? for the NB1 and NB2 
# parametrization of negative binomial in Cameron and Trivedi (1998)
# link is the link function. Choices are ?log? for the log link function, 
# ?logit? for the logit link function, and ?probit? for the probit link function.
# output: A list containing the following components:
# coef the vector with the ML estimated regression parameters
# gam the ML estimate of gamma parameter
iee<-function(xdat,ydat,margmodel,link="log")
{ #if(margmodel=="bernoulli")  family=binomial else family=poisson
  if(margmodel=="nb1" | margmodel=="nb2" | margmodel=="poisson") link="log"
  if(margmodel=="bernoulli" & link!="probit") link="logit"
  if(margmodel=="poisson")
  { uni<-glm(ydat ~ xdat[,-1],family =poisson(link="log"))
    res<-as.vector(uni$coef)
    list(reg=res)
  } else {
  if(margmodel=="bernoulli")
  { if(link=="probit") 
    { uni<-glm(ydat ~ xdat[,-1],family =binomial(link="probit"))
    } else {
    uni<-glm(ydat ~ xdat[,-1],family =binomial(link="logit")) }     
    res<-as.vector(uni$coef)
    list(reg=res)
  } else
  { p<-dim(xdat)[2]
    uni<-nlm(marglik,c(rep(0,p),1),margmodel=margmodel,
    link=link,xdat=xdat,ydat=ydat,iterlim=1000)
    res1<-uni$e[1:p]
    res2<-uni$e[p+1]
    list(reg=res1,gam=res2) }}
}



# Bivariate composite likelihood for multivariate normal copula with Poisson, 
# binary, or negative binomial regression.
# input:
# r the vector of normal copula parameters
# b the regression coefficients
# gam the gamma parameter
# xdat the matrix of covariates (use the constant 1 for the first covariate)
# ydat the vector with the response
# margmodel indicates the marginal model. Choices are ?poisson? for Poisson, 
# ?bernoulli? for Bernoulli, and ?nb1? , ?nb2? for the NB1 and NB2 
# parametrization of negative binomial in Cameron and Trivedi (1998)
# link is the link function. Choices are ?log? for the log link function, 
# ?logit? for the logit link function, and ?probit? for the probit link function.
# output:
# negative bivariate composite likelihood for multivariate normal copula 
# with Poisson, binary, or negative binomial regression.
bcl<-function(rh,p,q,b,gam,xdat,ydat,margmodel,link="log")
{ phi=rh[1:p]
  if(q>0) { th=rh[(p+1):(p+q)] } else {th=NULL}
  c1=any(Mod(polyroot(c(1, -phi))) < 1.01)
  c2=any(Mod(polyroot(c(1, th))) < 1.01)
  if(c1 || c2) return(1e10)
  s<-0
  if(margmodel=="nb1" | margmodel=="nb2" | margmodel=="poisson") link="log"
  if(margmodel=="bernoulli" & link!="probit") link="logit" 
  if(margmodel=="nb1" | margmodel=="nb2") invgam<-1/gam else invgam<-NULL
  d<-length(ydat)
  rmat<-toeplitz(ARMAacf(ar=phi, ma=th, lag.max=d-1))            
  mu<-linked.mu(xdat,b,link)
  vlow<-pmargmodel(ydat-1,mu,gam,invgam,margmodel)
  tem<-dmargmodel(ydat,mu,gam,invgam,margmodel)
  vupp<-vlow+tem
  zlow=qnorm(vlow)
  zupp=qnorm(vupp)
  bivpairs=combn(1:d, 2)
  d2=ncol(bivpairs)
  zmat=matrix(NA,d2,5)
  for(j in 1:d2)
  { k1=bivpairs[,j][1]
    k2=bivpairs[,j][2]
    zmat[j,]=c(zlow[c(k1,k2)],zupp[c(k1,k2)],rmat[k1,k2])  
  }
  zmat[zmat==-Inf]=-10
  #prob=pbivnorm(zmat[,3],zmat[,4],zmat[,5])+
  #  pbivnorm(zmat[,1],zmat[,2],zmat[,5])-
  #  pbivnorm(zmat[,3],zmat[,2],zmat[,5])-
  #  pbivnorm(zmat[,1],zmat[,4],zmat[,5])
  prob=pbvt(zmat[,3],zmat[,4],cbind(zmat[,5],100))+
      pbvt(zmat[,1],zmat[,2],cbind(zmat[,5],100))-
      pbvt(zmat[,3],zmat[,2],cbind(zmat[,5],100))-
      pbvt(zmat[,1],zmat[,4],cbind(zmat[,5],100))
  -sum(log(prob))
} 


# optimization routine for composite likelihood for MVN copula
# b the regression coefficients
# gam the gamma parameter
# xdat the matrix of covariates (use the constant 1 for the first covariate)
# ydat the vector with the response
# margmodel indicates the marginal model. Choices are ?poisson? for Poisson, 
# ?bernoulli? for Bernoulli, and ?nb1? , ?nb2? for the NB1 and NB2 
# parametrization of negative binomial in Cameron and Trivedi (1998)
# link is the link function. Choices are ?log? for the log link function, 
# ?logit? for the logit link function, and ?probit? for the probit link function.
# output: A list containing the following components:
# minimum the value of the estimated minimum of CL1 for MVN copula
# estimate the CL1 estimates
# gradient the gradient at the estimated minimum of CL1
# code an integer indicating why the optimization process terminated, see nlm.
cl<-function(p,q,b,gam,xdat,ydat,margmodel,link="log")
{ if(margmodel=="nb1" | margmodel=="nb2" | margmodel=="poisson") link="log"
  if(margmodel=="bernoulli" & link!="probit") link="logit"
  mod <- glm(ydat~xdat[,-1], family="poisson")
  temp <- arima(resid(mod) , order = c(p,0,q))
  start<- temp$coef[1:(p+q)]
  nlm(bcl,start,p,q,b,gam,xdat,ydat,margmodel,link,print.level = 2)
}




# derivative of the marginal loglikelihood with respect to nu
# input:
# mu the mean parameter
# gam the gamma parameter
# invgam the inverse of gamma parameter
# ub the truncation value
# margmodel indicates the marginal model. Choices are ?poisson? for Poisson,
# ?bernoulli? for Bernoulli, and ?nb1? , ?nb2? for the NB1 and NB2
# parametrization of negative binomial in Cameron and Trivedi (1998)
# link is the link function. Choices are ?log? for the log link function,
# ?logit? for the logit link function, and ?probit? for the probit link function.
# output:
# the vector with the derivatives of the margmodel loglikelihood with respect to nu
derlik.nu<-function(mu,gam,invgam,ub,margmodel,link)
{ if(link=="probit")
  { if(mu==1 & is.finite(mu)){mu<-0.9999}
    nu<-qnorm(mu)
    (0:1-mu)/mu/(1-mu)*dnorm(nu)
  }
  else
  { if(margmodel=="nb1")
    { j<-0:(ub-1)
      s<-c(0,cumsum(1/(mu+gam*j)))
      (s-invgam*log(1+gam))*mu
    }
  else
  { if(margmodel=="nb2")
    { pr<-1/(mu*gam+1)
      (0:ub-mu)*pr
    }
  else { 0:ub-mu }}}
}





# derivative of the marginal loglikelihood with respect to nu
# input:
# mu the mean parameter
# gam the gamma parameter
# invgam the inverse of gamma parameter
# y  the value of a non-negative integer quantile
# margmodel indicates the marginal model. Choices are ?poisson? for Poisson,
# ?bernoulli? for Bernoulli, and ?nb1? , ?nb2? for the NB1 and NB2
# parametrization of negative binomial in Cameron and Trivedi (1998)
# link is the link function. Choices are ?log? for the log link function,
# ?logit? for the logit link function, and ?probit? for the probit link function.
# output:
# the derivative of the margmodel loglikelihood with respect to nu
iderlik.nu<-function(mu,gam,invgam,y,margmodel,link)
{ if(link=="probit")
  { if(mu==1 & is.finite(mu)){mu<-0.9999}
    nu<-qnorm(mu)
    (y-mu)/mu/(1-mu)*dnorm(nu)
  }
  else
  { if(margmodel=="nb1")
    { s<-0
      if(y>0)
      { j<-0:(y-1)
        s<-sum(1/(mu+gam*j))
      }
      (s-invgam*log(1+gam))*mu
    }
  else
  { if(margmodel=="nb2")
    { pr<-1/(mu*gam+1)
      (y-mu)*pr
    }
  else { y-mu }}}
}

# derivative of the NB loglikelihood with respect to gamma
# input:
# gam the gamma parameter
# invgam the inverse of gamma parameter
# mu the mean parameter
# ub the truncation value
# margmodel indicates the marginal model. Choices are  ?nb1? , ?nb2? for
# the NB1 and NB2  parametrization of negative binomial in
# Cameron and Trivedi (1998)
# output:
# the vector with the derivatives of the NB loglikelihood with respect to gamma
derlik.gam<-function(mu,gam,invgam,ub,margmodel)
{ j<-0:(ub-1)
  if(margmodel=="nb1")
  { s<-c(0,cumsum(j/(mu+gam*j)))
    s+invgam*invgam*mu*log(1+gam)-(0:ub+invgam*mu)/(1+gam)
  }
  else
  { #if(margmodel=="nb2")
    pr<-1/(mu*gam+1)
    s<-c(0,cumsum(j/(1+j*gam)))
    s-log(pr)/(gam*gam)-(0:ub+invgam)*mu*pr
  }
}

# derivative of the NB loglikelihood with respect to gamma
# input:
# gam the gamma parameter
# invgam the inverse of gamma parameter
# mu the mean parameter
# y  the value of a non-negative integer quantile
# margmodel indicates the marginal model. Choices are  ?nb1? , ?nb2? for
# the NB1 and NB2  parametrization of negative binomial in
# Cameron and Trivedi (1998)
# output:
# the derivative of the NB loglikelihood with respect to gamma
iderlik.gam<-function(mu,gam,invgam,y,margmodel)
{ s<-0
  if(margmodel=="nb1")
  { if(y>0)
  { j<-0:(y-1)
    s<-sum(j/(mu+gam*j))
  }
  s+invgam*invgam*mu*log(1+gam)-(y+invgam*mu)/(1+gam)
  }
  else
  { #if(margmodel=="nb2")
    if(y>0)
  { j<-0:(y-1)
    s<-sum(j/(1+gam*j))
  }
  pr<-1/(mu*gam+1)
  s-log(pr)/(gam*gam)-(y+invgam)*mu*pr
  }
}



# minus expectation of the second derivative of the marginal loglikelihood
# with resect to nu
# input:
# mu the mean parameter
# gam the gamma parameter
# invgam the inverse of gamma parameter
# u the univariate cdfs
# ub the truncation value
# margmodel indicates the marginal model. Choices are  ?nb1? , ?nb2? for
# the NB1 and NB2  parametrization of negative binomial in
# Cameron and Trivedi (1998)
# output:
# the vector with the minus expectations of the margmodel loglikelihood
# with respect to nu
fisher.nu<-function(mu,gam,invgam,u,ub,margmodel,link)
{ if(link=="log" & margmodel=="poisson")
    { mu }
    else
    { if(link=="logit")
      { mu*(1-mu) }
    else
    { if(margmodel=="nb1")
      { j<-0:ub
        s1<-sum(1/(mu+j*gam)/(mu+j*gam)*(1-u))
        s2<-sum(1/(mu+j*gam)*(1-u))
        (mu*s1-s2+invgam*log(1+gam))*mu
      }
    else
    { if(margmodel=="nb2")
      { pr<-1/(mu*gam+1)
        mu*pr
      }
    else
    { # link=="probit"
      if(mu==1 & is.finite(mu)){mu<-0.9999}
      nu<-qnorm(mu)
      1/mu/(1-mu)*dnorm(nu)^2}
    }}}
}




# minus expectation of the second derivative of the marginal NB loglikelihood
# with resect to gamma
# input:
# mu the mean parameter
# gam the gamma parameter
# invgam the inverse of gamma parameter
# u the univariate cdfs
# ub the truncation value
# margmodel indicates the marginal model. Choices are  ?nb1? , ?nb2? for
# the NB1 and NB2  parametrization of negative binomial in
# Cameron and Trivedi (1998)
# output:
# the vector with the minus expectations of the margmodel loglikelihood
# with respect to gamma
fisher.gam<-function(mu,gam,invgam,u,ub,margmodel)
{ j<-0:ub
  if(margmodel=="nb1")
  { pr<-1/(gam+1)
    s<-sum((j/(mu+j*gam))^2*(1-u))
    s+2*invgam*invgam*invgam*log(1+gam)*mu-2*invgam*invgam*mu*pr-mu/gam*pr
  }
  else
  { #if(margmodel=="nb2")
    s<-sum((invgam+j)^(-2)*(1-u))
    invgam^4*(s-gam*mu/(mu+invgam))
  }
}






# minus expectation of the second derivative of the marginal loglikelihood
# with resect to nu and gamma
# input:
# mu the mean parameter
# gam the gamma parameter
# invgam the inverse of gamma parameter
# u the univariate cdfs
# ub the truncation value
# margmodel indicates the marginal model. Choices are  ?nb1? , ?nb2? for
# the NB1 and NB2  parametrization of negative binomial in
# Cameron and Trivedi (1998)
# output:
# the vector with the minus expectations of the NB loglikelihood
# with respect to nu and gamma
fisher.nu.gam<-function(mu,gam,invgam,u,ub,margmodel)
{ if(margmodel=="nb1")
  { pr<-1/(gam+1)
    j<-0:ub
    s<-sum(j/(mu+j*gam)/(mu+j*gam)*(1-u))
    (s-invgam*invgam*log(1+gam)+invgam*pr)*mu
  }
  else {0}
}

# Calculating the truncation value for the univariate distribution
# input:
# mu the mean parameter of the univariate distribution.
# gam the parameter gamma of  the negative binomial distribution.
# margmodel indicates the marginal model. Choices are ?poisson? for Poisson, ?bernoulli? for
# Bernoulli, and ?nb1? , ?nb2? for the NB1 and NB2 parametrization of negative
# binomial in Cameron and Trivedi (1998).
# output:
# the truncation value--upper bound
truncation<-function(mu,gam,margmodel)
{ if(margmodel=="poisson")
    { ub<-round(max(10,mu+7*sqrt(mu),na.rm=T))
    }
    else
    { if(margmodel=="bernoulli") ub<-1
    else
    { if(margmodel=="nb1")
      { pr<-1/(gam+1)
        v<-mu/pr
        ub<-round(max(10,mu+10*sqrt(v),na.rm=T))
       }
    else
    { pr<-1/(mu*gam+1)
      v<-mu/pr
      ub<-round(max(10,mu+7*sqrt(v),na.rm=T))
    }}}
  ub
}




# covariance matrix of the scores Omega_i
# input:
# scnu the matrix of the score functions with respect to nu
# scgam the matrix of the score functions with respect to gam
# index the bivariate pair
# pmf the matrix of rectangle probabilities
# margmodel indicates the marginal model. Choices are ?poisson? for Poisson, ?bernoulli? for
# Bernoulli, and ?nb1? , ?nb2? for the NB1 and NB2 parametrization of negative
# binomial in Cameron and Trivedi (1998).
scoreCov<-function(scnu,scgam,pmf,index,margmodel)
{ j1<-index[1]
  j2<-index[2]
  cov11<-t(scnu[,j1])%*%pmf%*%scnu[,j2]
  if(margmodel=="bernoulli" | margmodel=="poisson")
  { cov11 }
  else
  { cov12<-t(scnu[,j1])%*%pmf%*%scgam[,j2]
    cov21<-t(scgam[,j1])%*%pmf%*%scnu[,j2]
    cov22<-t(scgam[,j1])%*%pmf%*%scgam[,j2]
    matrix(c(cov11,cov12,cov21,cov22),2,2)
  }
}



# weight matrix fixed at values from the CL1 estimator
# input:
# b the regression coefficients
# gam the gamma parameter
# xdat the matrix of covariates (use the constant 1 for the first covariate)
# ydat the vector with the response
# rh the vector with CL1 estimates
# margmodel indicates the marginal model. Choices are ?poisson? for Poisson, 
# ?bernoulli? for Bernoulli, and ?nb1? , ?nb2? for the NB1 and NB2 
# parametrization of negative binomial in Cameron and Trivedi (1998)
# link is the link function. Choices are ?log? for the log link function, 
# ?logit? for the logit link function, and ?probit? for the probit link function.
# output: A list containing the following components:
# omega the array with the Omega matrices
# delta the array with the Delta matrices
# X the array with the X matrices
weightMat<-function(b,gam,rh,p,q,xdat,margmodel,link="log")
{ phi=rh[1:p]
  if(q>0){ th=rh[(p+1):(p+q)] } else { th=NULL }
  if(margmodel=="nb1" | margmodel=="nb2" | margmodel=="poisson") link="log"
  if(margmodel=="bernoulli" & link!="probit") link="logit"
  if(margmodel=="nb1" | margmodel=="nb2") invgam<-1/gam else invgam<-NULL
  d<-nrow(xdat)
  bivpairs=t(combn(1:d, 2))
  d2=nrow(bivpairs)
  rmat=toeplitz(ARMAacf(ar=phi, ma=th, lag.max=d-1))
  qq<-length(gam)
  pp<-ncol(xdat)
  omega<-matrix(NA,d*(1+qq),d*(1+qq))
  X<-matrix(NA,pp+qq,d*(1+qq))
  delta<-matrix(NA,d*(1+qq),d*(1+qq))
  dom<-d*(1+qq)
  pos<-seq(1,dom-1,by=2)    #not used for binary
  mu<-linked.mu(xdat,b,link)
  ub<-truncation(mu,gam,margmodel)
  du<-scnu<-scgam<-matrix(NA,1+ub,d)
  for(j in 1:d)
    { du[,j]<-dmargmodel(0:ub,mu[j],gam,invgam,margmodel)
      scnu[,j]<-derlik.nu(mu[j],gam,invgam,ub,margmodel,link)
      if(margmodel=="nb1" | margmodel=="nb2")
      { scgam[,j]<-derlik.gam(mu[j],gam,invgam,ub,margmodel) }
    }
    u<-apply(du,2,cumsum)
    z<-qnorm(u)
    z[is.nan(z)]<-10 
    z[z==Inf]=10 
    z[z>4&margmodel=="bernoulli"]<-7
    z<-rbind(-10,z)
    x<-NULL
    if(margmodel=="bernoulli" | margmodel=="poisson")
    { diagonal<-rep(NA,d)
    } else {
    diagonal<-array(NA,c(2,2,d)) }
    for(j in 1:d)
    { f1<-fisher.nu(mu[j],gam,invgam,u[,j],ub,margmodel,link)
      if(margmodel=="nb1" | margmodel=="nb2")
      { temp<-cbind(xdat[j,],0)
        x<-cbind(x,rbind(temp,c(0,1)))
        f2<-fisher.gam(mu[j],gam,invgam,u[,j],ub,margmodel)
        f3<-fisher.nu.gam(mu[j],gam,invgam,u[,j],ub,margmodel)
        diagonal[,,j]<-matrix(c(f1,f3,f3,f2),2,2)
      }
      else
      { temp<-xdat[j,]
        x<-cbind(x,temp)
        diagonal[j]<-f1
      }
    }
    if(margmodel=="bernoulli" | margmodel=="poisson")
    { delta<-diag(diagonal)
      off<-rep(NA,d2)
    } else {
      delta<-matrix(0,dom,dom)
      minus<-0
      for(j in pos)
      { delta[j:(j+1),j:(j+1)]<-diagonal[,,j-minus]
        minus<-minus+1
      }
      off<-array(NA,c(2,2,d2))
    }
    
    for(k in 1:d2)
    { print(k)
      k1<-bivpairs[k,][1]
      k2<-bivpairs[k,][2]
      x1=z[,k1]; x2=z[,k2]
      xmat=meshgrid(x1,x2)
      cdf=t(pbvt(xmat$x,xmat$y,c(rmat[k1,k2],1000)))
      cdf1=apply(cdf,2,diff)
      pmf=apply(t(cdf1),2,diff)
      pmf=t(pmf)
      if(margmodel=="bernoulli" | margmodel=="poisson")
      {off[k]<-scoreCov(scnu,scgam,pmf,bivpairs[k,],margmodel)}
      else {off[,,k]<-scoreCov(scnu,scgam,pmf,bivpairs[k,],margmodel)}
     }
    omega<-delta
    if(margmodel=="bernoulli" | margmodel=="poisson")
    { for(j in 1:d2)
      { omega[bivpairs[j,1],bivpairs[j,2]]<-off[j]
        omega[bivpairs[j,2],bivpairs[j,1]]<-off[j]
      }}
    else
    { ch1<-0
      ch2<-0
      for(j in 1:(d-1))
      { for(r in pos[-(1:j)])
        { omega[(j+ch1):(j+1+ch1),r:(r+1)]<-off[,,(j+ch2-ch1)]
          omega[r:(r+1),(j+ch1):(j+1+ch1)]<-t(off[,,(j+ch2-ch1)])
          ch2<-ch2+1
        }
      ch1<-ch1+1
    }}
    list(omega=omega,X=x,delta=delta)
}


# the wcl estimating equations
# input:
# param the vector of regression and not regression parameters
# xdat the matrix of covariates (use the constant 1 for the first covariate)
# ydat the vector with the response
# margmodel indicates the marginal model. Choices are ?poisson? for Poisson, 
# ?bernoulli? for Bernoulli, and ?nb1? , ?nb2? for the NB1 and NB2 
# parametrization of negative binomial in Cameron and Trivedi (1998)
# link is the link function. Choices are ?log? for the log link function, 
# ?logit? for the logit link function, and ?probit? for the probit link function.
# WtScMat is a list containing the following components:
# omega the array with the Omega matrices
# delta the array with the Delta matrices
# X the array with the X matrices
# output
# the weigted scores equations
bwcl<-function(param,WtScMat,xdat,ydat,margmodel,link="log")
{ if(margmodel=="nb1" | margmodel=="nb2" | margmodel=="poisson") link="log"
  if(margmodel=="bernoulli" & link!="probit") link="logit"
  d<-length(ydat)
  p<-ncol(xdat)
  b<-param[1:p]
  if(p<length(param)) {gam<-param[p+1]; invgam<-1/gam }
  mu<-linked.mu(xdat,b,link)
  ub<-truncation(mu,gam,margmodel)
  scnu<-scgam<-matrix(NA,ub+1,d)
  for(j in 1:d)
  { scnu[,j]<-derlik.nu(mu[j],gam,invgam,ub,margmodel,link)
    if(margmodel=="nb1" | margmodel=="nb2")
      { scgam[,j]<-derlik.gam(mu[j],gam,invgam,ub,margmodel) }
  }
  sc<-NULL
    for(j in 1:d)
    { if(ydat[j]>ub)
      { scnui<-iderlik.nu(mu[j],gam,invgam,ydat[j],margmodel,link)
        if(margmodel=="nb1" | margmodel=="nb2")
        { scgami<-iderlik.gam(mu[j],gam,invgam,ydat[j],margmodel)
        } else {
        scgami<-NULL}
      }
      else {
      scnui<-scnu[ydat[j]+1,j]
      if(margmodel=="nb1" | margmodel=="nb2")
      { scgami<-scgam[ydat[j]+1,j]
      } else {
        scgami<-NULL}}
        sc<-c(sc,c(scnui,scgami))
    }
    X<-WtScMat$X
    delta<-WtScMat$delta
    omega<-WtScMat$omega
    g<-X%*%t(delta)%*%solve(omega,sc)
    g
}


# solving the wcl estimating equations
# input:
# start the starting values (IEE estimates) for the vector of
# regression and not regression parameters
# xdat the matrix of covariates (use the constant 1 for the first covariate)
# ydat the vector with the response
# margmodel indicates the marginal model. Choices are ?poisson? for Poisson, 
# ?bernoulli? for Bernoulli, and ?nb1? , ?nb2? for the NB1 and NB2 
# parametrization of negative binomial in Cameron and Trivedi (1998)
# link is the link function. Choices are ?log? for the log link function, 
# ?logit? for the logit link function, and ?probit? for the probit link function.
# WtScMat is a list containing the following components:
# omega the array with the Omega matrices
# delta the array with the Delta matrices
# X the array with the X matrices
# output:
# the wcl estimates
wcl<-function(start,WtScMat,xdat,ydat,margmodel,link="log")
{ multiroot(f=bwcl,start,atol=1e-4,rtol=1e-4,ctol=1e-4,
  WtScMat=WtScMat,xdat=xdat,ydat=ydat,margmodel=margmodel,link=link) 
}



godambe=function(b,gam,rh,p,q,xdat,margmodel,link="log")
{ WtScMat<-weightMat(b,gam,rh,p,q,xdat,margmodel,link)
  omega= WtScMat$omega
  delta= WtScMat$delta
  X= WtScMat$X
  psi<-delta%*%t(X)
  hess<-t(psi)%*%solve(omega)%*%psi
  solve(hess)
}

 
