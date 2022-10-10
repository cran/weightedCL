

# Density  of the univariate marginal distribution
# input:
# y the vector of (non-negative integer) quantiles.
# mu the mean parameter of the univariate distribution.
# gam the parameter gamma of  the negative binomial distribution.
# output:
# the density of the univariate marginal  distribution
dmargmodel.ord<-function(y,mu,gam,link)
{ cuts<-c(-10,gam,10)
lb=cuts[y]+mu
ub=cuts[y+1]+mu
if(link=="probit") res<-pnorm(ub)-pnorm(lb) else res<-plogis(ub)-plogis(lb)
res[y<1]<-0
res
}

# CDF  of the univariate marginal distribution
# input:
# y the vector of (non-negative integer) quantiles.
# mu the mean parameter of the univariate distribution.
# gam the parameter gamma of  the negative binomial distribution.
# output:
# the cdf of the univariate marginal  distribution
pmargmodel.ord<-function(y,mu,gam,link)
{ cuts<-c(-10,gam,10)
  ub=cuts[y+1]+mu # for mprobit
  #ub=cuts[y+1]-mu # for polr
  if(link=="probit") res<-pnorm(ub) else res<-plogis(ub)
  res[y<1]<-0
  res
}



# Bivariate composite likelihood for multivariate normal copula with ordinal regression.
# input:
# r the vector of normal copula parameters
# b the regression coefficients
# gam the gamma parameter
# xdat the matrix of covariates (use the constant 1 for the first covariate)
# ydat the vector with the response
# link is the link function. Choices are ?log? for the log link function, 
# ?logit? for the logit link function, and ?probit? for the probit link function.
# output:
# negative bivariate composite likelihood for multivariate normal copula 
# with ordinal regression.
bcl.ord<-function(rh,p,q,b,gam,xdat,ydat,link)
{ phi=rh[1:p]
  if(q>0) { th=rh[(p+1):(p+q)] } else {th=NULL}
  c1=any(Mod(polyroot(c(1, -phi))) < 1.01)
  c2=any(Mod(polyroot(c(1, th))) < 1.01)
  if(c1 || c2) return(1e10)
  d=length(ydat)
  rmat<-toeplitz(ARMAacf(ar=phi, ma=th, lag.max=d-1))            
  mu<-ordreg.mu(xdat,b)
  vlow<-pmargmodel.ord(ydat-1,mu,gam,link)
  tem<-dmargmodel.ord(ydat,mu,gam,link)
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
  zmat[zmat==Inf]=10
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
# link is the link function. Choices are  
# ?logit? for the logit link function, and ?probit? for the probit link function.
# output: A list containing the following components:
# minimum the value of the estimated minimum of CL1 for MVN copula
# estimate the CL1 estimates
# gradient the gradient at the estimated minimum of CL1
# code an integer indicating why the optimization process terminated, see nlm.
cl.ord<-function(p,q,b,gam,xdat,ydat,link)
{ if(link=="logit") link="logistic"
  t1=polr(as.factor(ydat)~xdat,method=link)
  #t2 <- arima(residuals.polr(t1)[,1] , order = c(p,0,q))
  t2 <- arima(resids(t1), order = c(p,0,q))
  start<- t2$coef[1:(p+q)]
  nlm(bcl.ord,start,p,q,b,gam,xdat,ydat,link,
      print.level=2)
}


# derivative of the ordinal loglikelihood with respect to gamma
# input:
# gam the gamma parameter
# mu the mean parameter
# output:
# the vector with the derivatives of the NB loglikelihood with respect to gamma
derlik.gam.ord<-function(mu,gam,u,link)
{ K<-length(gam)+1
  k<-1:K
  cuts<-c(-10,gam,10)
  lb=cuts[k]+mu
  ub=cuts[k+1]+mu
  if(link=="probit") { dlatent=dnorm; platent=pnorm } else { dlatent=dlogis; platent=plogis }
  dlatentub<-dlatent(ub)
  dlatentlb<-dlatent(lb)
  den<-platent(ub)-platent(lb)
  res<-rep(NA,K)
  for(i in 1:K)
  { if(u==i)
  { res[i]=dlatentub[i]/den[i] }
  else
  { if(u==i-1)
  { res[i]=-dlatentlb[i]/den[i] }
    else {res[i]=0}}
}
res
}

# derivative of the NB loglikelihood with respect to gamma
# input:
# gam the gamma parameter
# mu the mean parameter
# y  the value of a non-negative integer quantile
# output:
# the derivative of the NB loglikelihood with respect to gamma
iderlik.gam.ord<-function(mu,gam,y,u,link)
{ cuts<-c(-10,gam,10)
  lb=cuts[y]+mu
  ub=cuts[y+1]+mu
  if(link=="probit") { dlatent=dnorm; platent=pnorm } else { dlatent=dlogis; platent=plogis }
  den<-platent(ub)-platent(lb)
  if(u==y) dlatent(ub)/den
  else if(u==y-1) -dlatent(lb)/den  else 0
}


der.dnorm<-function(x)
{ -x*dnorm(x) }

der.dlogis<-function(x)
{ expx=exp(x)
  expx*(1-expx)/(1+expx)^3 
}

# minus expectation of the second derivative of the marginal ordinal loglikelihood
# with resect to gamma
# input:
# mu the mean parameter
# gam the gamma parameter
# u the univariate cdfs
# output:
# the vector with the minus expectations of the margmodel loglikelihood
# with respect to gamma
fisher.gam.ord<-function(mu,gam,u,v,link)
{ cuts<-c(-10,gam,10)
  K<-length(gam)+1
  k<-1:K
  lb=cuts[k]+mu
  ub=cuts[k+1]+mu
  if(link=="probit") { dlatent=dnorm; platent=pnorm; der.dlatent=der.dnorm } else { 
    dlatent=dlogis; platent=plogis; der.dlatent=der.dlogis }
  den<-platent(ub)-platent(lb)
  dlatentub<-dlatent(ub)
  dlatentlb<-dlatent(lb)
  der.dlatentub<-der.dlatent(ub)
  der.dlatentlb<-der.dlatent(lb)
  h<-rep(NA,K)
  for(k in 1:K)
  { if(u==k & v==k)
  { num1<-der.dlatentub[k]
    num2<-dlatentub[k]
    tem<-num2/den[k]
    h[k]<--num1/den[k]+tem*tem  
  }
  else
  { if((u==k & v==k-1) | (u==k-1 & v==k))
  { h[k]<--dlatentub[k]*dlatentlb[k]/den[k]/den[k]  
  }
    else
    { if(u==k-1 & v==k-1)
    { num1<-der.dlatentlb[k]
    num2<-dlatentlb[k]
    tem<-num2/den[k]
    h[k]<-num1/den[k]+tem*tem  
    } else h[k]<-0}}
 }
 sum(h*den) 
}

# the mean values of the univariate marginal distribution 
# corresonding to the used link function 
# input:
# x the matix of the covariates 
# b the vector with the regression coefficients
# output:
# the mean values of the univariate marginal distribution 
ordreg.mu<-function(x,b)
{ if(length(b)!=1) mu<-x %*% b else mu<-x*b }



# Calculating the number of categories
# input:
# gam the cutpoints
# output:
# the number of categories
noCategories<-function(gam)
{ length(gam)+1 }

# covariance matrix of the scores Omega_i
# input:
# scgam the array of the score functions with respect to gam
# index the bivariate pair
# pmf the matrix of rectangle probabilities
scoreCov.ord<-function(scgam,pmf,index)
{ j1<-index[1]
  j2<-index[2]
  q<-dim(scgam)[2]
  cov22<-matrix(NA,q,q)
  for(i1 in 1:q)
  { for(i2 in 1:q)
  { cov22[i1,i2]<-t(scgam[,i1,j1])%*%pmf%*%scgam[,i2,j2] }
  }
  cov22
}

# weight matrix fixed at values from the CL1 estimator
# input:
# b the regression coefficients
# gam the gamma parameter
# xdat the matrix of covariates (use the constant 1 for the first covariate)
# ydat the vector with the response
# rh the vector with CL1 estimates
# link is the link function. Choices are 
# ?logit? for the logit link function, and ?probit? for the probit link function.
# output: A list containing the following components:
# omega the array with the Omega matrices
# delta the array with the Delta matrices
# X the array with the X matrices
weightMat.ord<-function(b,gam,rh,p,q,xdat,link)
{ phi=rh[1:p]
  if(q>0){ th=rh[(p+1):(p+q)] } else { th=NULL }
  d=nrow(xdat)
  rmat=toeplitz(ARMAacf(ar=phi, ma=th, lag.max=d-1))
  qq<-length(gam)
  if(is.matrix(xdat))
  { dim<-dim(xdat)
    pp<-dim[2]
  } else {pp<-1}
  bivpairs=t(combn(1:d, 2))
  d2=nrow(bivpairs)
  omega<-matrix(NA,d*qq,d*qq)
  X<-matrix(NA,pp+qq,d*qq)
  delta<-matrix(NA,d*qq,d*qq)
  dom<-d*qq
  if(qq>1) pos<-seq(1,dom-1,by=qq) else pos<-seq(1,dom)    
  mu<-ordreg.mu(xdat,b)
  ub<-noCategories(gam)
  du<-matrix(NA,ub,d)
  scgam<-array(NA,c(ub,ub-1,d))
  for(j in 1:d)
  { du[,j]<-dmargmodel.ord(1:ub,mu[j],gam,link)
    for(k in 1:(ub-1))
    { scgam[,k,j]<-derlik.gam.ord(mu[j],gam,k,link) } 
  }
  u<-apply(du,2,cumsum)
  z<-qnorm(u)
  z[is.nan(z)]<-7
  z[z>10]<-7
  z<-rbind(-7,z)
  x<-NULL
  diagonal<-array(NA,c(qq,qq,d))   
  for(j in 1:d)
  { if(is.matrix(xdat)) 
  { temp1<-matrix(rep(xdat[j,],each=qq),qq) }  else { 
    temp1<-matrix(rep(xdat[j],each=qq),qq) }
    x<-cbind(x,t(cbind(temp1,diag(qq))))
    fisher<-matrix(NA,ub-1,ub-1)
    for(k1 in 1:(ub-1))
    { for(k2 in 1:(ub-1))
      { fisher[k1,k2]<-fisher.gam.ord(mu[j],gam,k1,k2,link)
      }
    }
    diagonal[,,j]<-fisher
  }
  delta<-matrix(0,dom,dom)
  minus<-0
  for(j in pos)
  { delta[j:(j+qq-1),j:(j+qq-1)]<-diagonal[,,j-minus]
       if(qq>1) minus<-minus+qq-1 else minus<-0
  }
  off<-array(NA,c(qq,qq,d2))
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
    off[,,k]<-scoreCov.ord(scgam,pmf,bivpairs[k,])
   }
   omega<-delta
   ch1<-0
   ch2<-0
   ch3<-0
   for(j in 1:(d-1))
   { for(r in pos[-(1:j)])
    { #print(c(j,r))
     omega[(1+(j-1)*qq):(j*qq),r:(r+qq-1)]<-off[,,(j+ch2-ch1)]
     omega[r:(r+qq-1),(1+(j-1)*qq):(j*qq)]<-t(off[,,(j+ch2-ch1)])
     ch2<-ch2+1
    }
    ch1<-ch1+1
   }
  list(omega=omega,X=x,delta=delta)
}


# the weigted scores equations
# input:
# param the vector of regression and not regression parameters
# xdat the matrix of covariates 
# ydat the vector with the response
# link is the link function. Choices are  
# ?logit? for the logit link function, and ?probit? for the probit link function.
# WtScMat is a list containing the following components:
# omega the array with the Omega matrices
# delta the array with the Delta matrices
# X the array with the X matrices
# output
# the wcl estimating equations
bwcl.ord<-function(param,WtScMat,xdat,ydat,link)
{ d<-length(ydat)
  if(is.matrix(xdat))
  { p<-ncol(xdat)
  } else {p<-1}
  b<-param[1:p]
  q<-length(unique(ydat))-1
  gam<-param[(p+1):(p+q)]
  mu<-ordreg.mu(xdat,b)
  ub<-noCategories(gam)
  scgam<-array(NA,c(ub,ub-1,d))
  for(j in 1:d)
  { for(k in 1:(ub-1))
  { scgam[,k,j]<-derlik.gam.ord(mu[j],gam,k,link) }
  }        
  sc<-NULL
  for(j in 1:d)
  { scgami<-NULL
      for(k in 1:(ub-1))
      { scgami<-c(scgami,scgam[ydat[j],k,j])}
      sc<-c(sc,scgami)
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
# link is the link function. Choices are 
# ?logit? for the logit link function, and ?probit? for the probit link function.
# WtScMat is a list containing the following components:
# omega the array with the Omega matrices
# delta the array with the Delta matrices
# X the array with the X matrices
# output:
# the wcl estimates
wcl.ord<-function(start,WtScMat,xdat,ydat,link)
{ multiroot(f=bwcl.ord,start,atol=1e-4,rtol=1e-4,ctol=1e-4,
            WtScMat=WtScMat,xdat=xdat,ydat=ydat,link=link) 
}


godambe.ord=function(b,gam,rh,p,q,xdat,link)
{ WtScMat<-weightMat.ord(b,gam,rh,p,q,xdat,link)
  omega= WtScMat$omega
  delta= WtScMat$delta
  X= WtScMat$X
  psi<-delta%*%t(X)
  hess<-t(psi)%*%solve(omega)%*%psi
  solve(hess)
}
