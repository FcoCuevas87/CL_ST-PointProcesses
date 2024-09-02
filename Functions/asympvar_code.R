#Sensitivity
#Births

sensi.birth.k=function(zeta.k,gradlogzeta.k,rho0){
  
  #implement estimate in Section 4.1
  
  #zeta.k and gradlogzeta.k are evaluated at data points as well as at dummy points
  #gradlogzeta.k is (n+m)xp where (n+m) is total number of points and p is number of parameters
  #gradlogzeta.k is here gradient of log zeta.
  
  temp=rho0*zeta.k/(zeta.k+rho0)^2
  
  #sensitivity is -1 times derivative of estm. fct.
  
  return(t(gradlogzeta.k)%*%diag(temp)%*%gradlogzeta.k)
  
}



sensi.birth.k.jy=function(zeta.k,gradlogzeta.k,is.datapoint,rho0){
  #browser()
  #implement estimate in Section 4.1
  #zeta.k and gradlogzeta.k are evaluated at data points as well as at dummy points
  #gradlogzeta.k is (n+m)xp where (n+m) is total number of points and p is number of parameters
  #gradlogzeta.k is here gradient of log zeta.
  temp=rho0/(zeta.k+rho0)
  phi=is.datapoint
  phi[!is.datapoint]=(zeta.k/rho0)[!is.datapoint]
  #sensitivity is -1 times derivative of estm. fct.
  return(t(gradlogzeta.k)%*%diag(temp*phi/2)%*%gradlogzeta.k)
  
}

e.birth.var.k.ours=function(zeta.k,gradlogzeta.k,rho0,is.datapoint,Bk,Y,bw){
    #first term in expression is sensitivity
    first.term=sensi.birth.k(zeta.k,gradlogzeta.k,rho0)
    #is.datapoint identifies datapoints part of zeta.k
    temp=rho0/(zeta.k+rho0)
    phi= is.datapoint
    phi[!is.datapoint]=-(zeta.k/rho0)[!is.datapoint]

    temp2=sweep(gradlogzeta.k,MARGIN=1,STATS=temp*phi,FUN="*")
    kernelmatrix=matrix(kernel(pairdist(superimpose(Bk,Y)),bw),ncol=Bk$n+Y$n)
    diag(kernelmatrix)=0
    return(first.term+t(temp2[is.datapoint,])%*%kernelmatrix[is.datapoint,is.datapoint]%*%temp2[is.datapoint,]+t(temp2[is.datapoint,])%*%kernelmatrix[is.datapoint,!is.datapoint]%*%temp2[!is.datapoint,]/2+t(temp2[!is.datapoint,])%*%kernelmatrix[!is.datapoint,is.datapoint]%*%temp2[is.datapoint,]/2)

}


e.birth.var.k=function(zeta.k,gradlogzeta.k,rho0,is.datapoint,Bk,Y,bw){
  #first term in expression is sensitivity
  first.term=sensi.birth.k(zeta.k,gradlogzeta.k,rho0)
  
  #is.datapoint identifies datapoints part of zeta.k
  temp=rho0/(zeta.k+rho0)
  temp1=temp[is.datapoint]
  temp2=sweep(gradlogzeta.k[is.datapoint,],MARGIN=1,STATS=temp1,FUN="*")
  temp=zeta.k/(zeta.k+rho0)
  temp3=temp[!is.datapoint]
  temp4=sweep(gradlogzeta.k[!is.datapoint,],MARGIN=1,STATS=temp3,FUN="*")

  kernelmatrixB=matrix(kernel(pairdist(Bk),bw),ncol=Bk$n)
  diag(kernelmatrixB)=0
  kernelmatrixY=matrix(kernel(pairdist(Y),bw),ncol=Y$n)
  diag(kernelmatrixY)=0
  
  return(first.term+t(temp2)%*%kernelmatrixB%*%temp2-t(temp4)%*%kernelmatrixY%*%temp4)
  
}



e.birth.var.k.jy=function(zeta.k,gradlogzeta.k,rho0,is.datapoint,Bk,Y,bw){
  #first term in expression is sensitivity
  first.term=sensi.birth.k(zeta.k,gradlogzeta.k,rho0)
  #is.datapoint identifies datapoints part of zeta.k
  temp=rho0/(zeta.k+rho0)
  phi=as.numeric(is.datapoint) 
  phi[!is.datapoint]=-(zeta.k/rho0)[!is.datapoint]
  temp2=sweep(gradlogzeta.k,MARGIN=1,STATS=temp*phi,FUN="*")
  kernelmatrix=matrix(kernel(pairdist(superimpose(Bk,Y)),bw),ncol=Bk$n+Y$n)
  diag(kernelmatrix)=0
  return(first.term+t(temp2)%*%kernelmatrix%*%temp2)
}

e.birth.var.k.ours=function(zeta.k,gradlogzeta.k,rho0,is.datapoint,Bk,Y,bw){
  #first term in expression is sensitivity
  first.term=sensi.birth.k(zeta.k,gradlogzeta.k,rho0)
  #is.datapoint identifies datapoints part of zeta.k
  temp=rho0/(zeta.k+rho0)
  phi=as.numeric(is.datapoint)
  phi[!is.datapoint]=-(zeta.k/rho0)[!is.datapoint]
  temp2=sweep(gradlogzeta.k,MARGIN=1,STATS=temp*phi,FUN="*")
  
  kernelmatrix=matrix(kernel(pairdist(superimpose(Bk,Y)),bw),ncol=Bk$n+Y$n)
  diag(kernelmatrix)=0
  
  return(first.term+t(temp2[is.datapoint,])%*%kernelmatrix[is.datapoint,is.datapoint]%*%temp2[is.datapoint,]+
           t(temp2[is.datapoint,])%*%kernelmatrix[is.datapoint,!is.datapoint]%*%temp2[!is.datapoint,]/2+
           t(temp2[!is.datapoint,])%*%kernelmatrix[!is.datapoint,is.datapoint]%*%temp2[is.datapoint,]/2)
  
}

#To verify
e.birth.var.k.alt=function(zeta.k,gradzeta.k,rho0,gB,gY){
  #first term in expression is sensitivity
  first.term=sensi.birth.k(zeta.k,gradzeta.k,rho0)
  
  temp=quadweights*rho0/(zeta.k+rho0)
  nqp=quadpoints$n
  pairdistances=c(pairdist(quadpints))
  gminus2=matrix(gB(pairdistances)+gY(pairdistances) - 2, nqp,nqp)
  
  temp2=sweep(gradzeta.k,MARGIN=1,STATS=temp,FUN="*")
  
  return(first.term+t(temp2)%*%gminus2%*%temp2)
}


#Deaths
sensi.death.k=function(grad.eta.k,eta.k){
  #browser()
  expeta=exp(eta.k)
  p=expeta/(1+expeta)
  temp=p^2/expeta
  
  return(t(grad.eta.k)%*%diag(temp)%*%grad.eta.k)
}

#Another try
sensi.death.k.alt2=function(grad.eta.k,eta.k,Ik){
  
  expeta=exp(eta.k)
  p=expeta/(1+expeta)
  temp=p*Ik/expeta#has mean p^2/expeta. Should be positively correlated with variance estimator.
  
  return(t(grad.eta.k)%*%diag(temp)%*%grad.eta.k)
}

#Another try
e.death.var.k.alt1=function(grad.eta.k,covmatrix){
  
  
  return(t(grad.eta.k)%*%covmatrix%*%grad.eta.k)
}

#Another try
e.death.var.k.alt3=function(grad.eta.k,variogram,Xkminus1){
  dists=pairdist(Xkminus1)
  covmatrixI=matrix(0,Xkminus1$n,Xkminus1$n)
  
  m=length(variogram$u)
  varians=mean(variogram$v[ceiling(0.75*m):m])
  
  #look up
  for (l in 1:Xkminus1$n)
    for (k in l:Xkminus1$n){
      if (dists[l,k]<variogram$max.dist)
        covmatrixI[l,k]=covmatrixI[k,l]=varians-variogram$v[abs(dists[l,k]-variogram$u)==min(abs(dists[l,k]-variogram$u))]
      else
        covmatrixI[l,k]=covmatrixI[k,l]=0
    }
  
  return(t(grad.eta.k)%*%covmatrixI%*%grad.eta.k)
}


#To use
e.death.var.k.alt2=function(grad.eta.k,eta.k,bw,Ik,Xkminus1){
  #grad.eta.k has i'th row given by gradient at i'th point. An n x p matrix.
  #Xkminus1 is point pattern of previous generation points
  #This estimator is somehow similar to White's covariance estimator in case of heteroscedasticity
  #browser()
  expeta=exp(eta.k)
  p=expeta/(1+expeta)
  #downweight or remove contribution from distant pairs using kernel. I think kernel should be function taking values from 0 to 1. E.g. kernel(d)=1[d < bw] or kernel(d)=exp(-(d/bw)^2). Using positive definite function is beneficial since then we return pos definite matrix.
  kernelmatrix=matrix(kernel(pairdist(Xkminus1),bw),ncol=Xkminus1$n)
  residual=Ik - p

  temp=sweep(grad.eta.k,MARGIN=1,STATS=residual,FUN="*")
  return(t(temp)%*%kernelmatrix%*%temp)
}


#independent DEATHS
e.death.var.k.independent=function(grad.eta.k,eta.k,Xkminus1){
  #case where deaths are conditionally independent
  
  expeta=exp(eta.k)
  p=expeta/(1+expeta)
  varI=p*(1-p)
  return(t(grad.eta.k)%*%diag(varI)%*%grad.eta.k)
}
#Compare above two procedures on simulated data.
#We can use vgm function in sp package for estimating covariogram

e.cross.var=function(eta.k,grad.logeta.k,zeta.k,grad.logzeta.k,Ik,Xkminus1,Bk,bw){
  #grad.logeta.k is covariate matrix for deaths in case of logistic regression
  #similarly, grad.logzeta.k is covariate matrix for births in case of log linear intensity function.
  #browser()
  expeta=exp(eta.k)
  p=expeta/(1+expeta)
  resi=Ik-p
  temp1=sweep(grad.logeta.k,MARGIN=1,STATS=resi,FUN="*")#matrix with Xkminus1$n rows
  temp2=sweep(grad.logzeta.k,MARGIN=1,STATS=rho0/(zeta.k+rho0),FUN="*")#matrix with Bk$n rows
  kernelmatrix=matrix(kernel(crossdist(Xkminus1,Bk),bw), nrow = Xkminus1$n, ncol = Bk$n)#matrix with Xkminus1$n rows and Bk$n columns
  
  return(t(temp1)%*%kernelmatrix%*%temp2)
  
}


