e.birth.var.k.ours.bw =function(zeta.k,gradlogzeta.k,rho0,is.datapoint,Bk,Y,bw){
  #first term in expression is sensitivity
  first.term=sensi.birth.k(zeta.k,gradlogzeta.k,rho0)
  #is.datapoint identifies datapoints part of zeta.k
  temp=rho0/(zeta.k+rho0)
  phi= is.datapoint
  phi[!is.datapoint]=-(zeta.k/rho0)[!is.datapoint]
  
  temp2=sweep(gradlogzeta.k,MARGIN=1,STATS=temp*phi,FUN="*")
  pd.obj <- pairdist(superimpose(Bk,Y))
  e.b.v <- list()
  for(i in 1:length(bw)){
    k0.mat <- kernel(pd.obj,bw[i])
    kernelmatrix=matrix(k0.mat,ncol=Bk$n+Y$n)
    diag(kernelmatrix)=0
    
    e.b.v[[i]] <- first.term+
                  t(temp2[is.datapoint,])%*%kernelmatrix[is.datapoint,is.datapoint]%*%temp2[is.datapoint,]+
                  t(temp2[is.datapoint,])%*%kernelmatrix[is.datapoint,!is.datapoint]%*%temp2[!is.datapoint,]/2+
                  t(temp2[!is.datapoint,])%*%kernelmatrix[!is.datapoint,is.datapoint]%*%temp2[is.datapoint,]/2
  }
  return(e.b.v)
  
}

e.birth.var.k.jy.bw = function(zeta.k,gradlogzeta.k,rho0,is.datapoint,Bk,Y,bw){
  #first term in expression is sensitivity
  first.term=sensi.birth.k(zeta.k,gradlogzeta.k,rho0)
  #is.datapoint identifies datapoints part of zeta.k
  temp=rho0/(zeta.k+rho0)
  phi=as.numeric(is.datapoint) 
  phi[!is.datapoint]=-(zeta.k/rho0)[!is.datapoint]
  temp2=sweep(gradlogzeta.k,MARGIN=1,STATS=temp*phi,FUN="*")
  e.b.v <- list()
  pd.obs <- pairdist(superimpose(Bk,Y))
  for(i in 1:length(bw)){
    k0 <- kernel(pd.obs,bw[i])
    kernelmatrix=matrix(k0,ncol=Bk$n+Y$n)
    diag(kernelmatrix)=0
    e.b.v[[i]] <- first.term+t(temp2)%*%kernelmatrix%*%temp2
  }
  return(e.b.v)
}

e.birth.var.k.bw=function(zeta.k,gradlogzeta.k,rho0,is.datapoint,Bk,Y,bw){
  #first term in expression is sensitivity
  first.term=sensi.birth.k(zeta.k,gradlogzeta.k,rho0)
  
  #is.datapoint identifies datapoints part of zeta.k
  temp=rho0/(zeta.k+rho0)
  temp1=temp[is.datapoint]
  temp2=sweep(gradlogzeta.k[is.datapoint,],MARGIN=1,STATS=temp1,FUN="*")
  temp=zeta.k/(zeta.k+rho0)
  temp3=temp[!is.datapoint]
  temp4=sweep(gradlogzeta.k[!is.datapoint,],MARGIN=1,STATS=temp3,FUN="*")
  pd.obj.1 <- pairdist(Bk)
  pd.obj.2 <- pairdist(Y)
  e.b.v <- list()
  for(i in 1:length(bw)){
    kernelmatrixB=matrix(kernel(pd.obj.1,bw[i]),ncol=Bk$n)
    diag(kernelmatrixB)=0
    kernelmatrixY=matrix(kernel(pd.obj.2,bw[i]),ncol=Y$n)
    diag(kernelmatrixY)=0
    e.b.v[[i]] <- first.term+
                  t(temp2)%*%kernelmatrixB%*%temp2-
                  t(temp4)%*%kernelmatrixY%*%temp4
  }
  
  return(e.b.v)
  
}

e.death.var.k.alt2.bw=function(grad.eta.k,eta.k,bw,Ik,Xkminus1){
  #grad.eta.k has i'th row given by gradient at i'th point. An n x p matrix.
  #Xkminus1 is point pattern of previous generation points
  #This estimator is somehow similar to White's covariance estimator in case of heteroscedasticity
  #browser()
  expeta=exp(eta.k)
  p=expeta/(1+expeta)
  #downweight or remove contribution from distant pairs using kernel. I think kernel should be function taking values from 0 to 1. E.g. kernel(d)=1[d < bw] or kernel(d)=exp(-(d/bw)^2). Using positive definite function is beneficial since then we return pos definite matrix.
  residual=Ik - p
  temp=sweep(grad.eta.k,MARGIN=1,STATS=residual,FUN="*")
  pd.obj <- pairdist(Xkminus1)
  e.v.d <- list()
  for(i in 1:length(bw)){
    kernelmatrix=matrix(kernel(pd.obj,bw[i]),ncol=Xkminus1$n)
    e.v.d[[i]] <- t(temp)%*%kernelmatrix%*%temp
  }
  return(e.v.d)
}


e.cross.var.bw=function(eta.k,grad.logeta.k,zeta.k,grad.logzeta.k,Ik,Xkminus1,Bk,bw){
  #grad.logeta.k is covariate matrix for deaths in case of logistic regression
  #similarly, grad.logzeta.k is covariate matrix for births in case of log linear intensity function.
  #browser()
  expeta=exp(eta.k)
  p=expeta/(1+expeta)
  resi=Ik-p
  temp1=sweep(grad.logeta.k,MARGIN=1,STATS=resi,FUN="*")#matrix with Xkminus1$n rows
  temp2=sweep(grad.logzeta.k,MARGIN=1,STATS=rho0/(zeta.k+rho0),FUN="*")#matrix with Bk$n rows
  cp.obj <- crossdist(Xkminus1,Bk)
  e.v.d <- list()
  for(i in 1:length(bw)){
    kernelmatrix=matrix(kernel(cp.obj,bw[i]), nrow = Xkminus1$n, ncol = Bk$n)#matrix with Xkminus1$n rows and Bk$n columns
    e.v.d[[i]] <- t(temp1)%*%kernelmatrix%*%temp2
  }
  return(e.v.d)
  
}
