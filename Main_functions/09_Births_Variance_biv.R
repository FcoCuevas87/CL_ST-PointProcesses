#04_Births_Variance.R

require("pbmcapply")
require("spatstat")
setwd("/data/fcp-X79-03/Spatio_temporal_PP/stBCIcode/SimulationStudy")
setting <- readRDS( paste(experiment, "/SimSetting_",N0,"_",prefix,"_all.RDS", sep = ""))
W.sim <- setting$window
prefix <- setting$prefix
N0 <- setting$N0
bw.seq <- setting$bw.seq
Time.steps <- setting$Time.steps

source("../asympvar_code.R")
source("../asympvar_code_bw.R")
source("CE_Method.R")
source("01_Simulation_Covariates.R")

kernel <- function(h,bw){
  return(as.numeric(h<bw))#simple uniform kernel used here.
}

result <- readRDS(paste(experiment, "/BDH_PP_",N0,"_",prefix,"_",setting$Kappa_fun,"_all_biv.RDS", sep = ""))
Births.fit <- readRDS(paste(experiment, "/Births_fit_",N0,"_",prefix,"_",setting$Kappa_fun,"_all_biv.RDS", sep = ""))

Covariance.estimate.births <- function(fit.B, Births, Births.dummy, Z.births, rho0, bw.seq, Time.steps){
  
  eval.covariate <- function(i, births){
    n0 <- length(Z.births[[i]])
    if(births){
      eval.list <- matrix(0, ncol = n0, nrow = Births[[i]]$n)
    }else{
      eval.list <- matrix(0, ncol = n0, nrow = Births.dummy[[i]]$n)
    }
    
    for(jj in 1:n0){
      if(births){
        eval.list[,jj] <- Z.births[[i]][[jj]][ Births[[i]] ]
      }else{
        eval.list[,jj] <- Z.births[[i]][[jj]][ Births.dummy[[i]] ]
      }
    }
    return(eval.list)
  }
  
  Sens.Var <- lapply(as.list(1:Time.steps), function(i){
    Births.cova <- eval.covariate(i, births = TRUE)
    Dummy.cova <- eval.covariate(i, births = FALSE)
    
    gradlogzeta.k  <- rbind(cbind(1,Births.cova), cbind(1,Dummy.cova))
    zeta.k <- exp( c(gradlogzeta.k%*%fit.B) )
    
    is.datapoint <- c(rep(TRUE, Births[[i]]$n), rep(FALSE, Births.dummy[[i]]$n) )
    
    sensiv.B.jy <- sensi.birth.k.jy(zeta.k, gradlogzeta.k, is.datapoint,rho0)
    
    #www <- pcfinhom(Bk, lambda = zeta.k[is.datapoint])
    #h0 <- www$r[ which.min( abs(www$iso-1) ) ]
    var.B.jy <- simplify2array(e.birth.var.k.jy.bw(zeta.k = zeta.k, 
                                                   gradlogzeta.k = gradlogzeta.k, 
                                                   rho0 = rho0, 
                                                   is.datapoint = is.datapoint, 
                                                   Births[[i]], 
                                                   Births.dummy[[i]], 
                                                   bw = bw.seq) )
    return(list(Sens = sensiv.B.jy, Var = var.B.jy))
  })
  
  sensiv.B <- lapply(as.list(1:Time.steps), function(i) Sens.Var[[i]]$Sens)
  var.B <- lapply(Sens.Var, function(x) x$Var)
  
  senssum.B <- Reduce("+", sensiv.B)
  varsum.B <- Reduce("+", var.B)
  variance <- lapply(as.list(1:length(bw.seq)), function(i) solve( senssum.B%*%solve(varsum.B[,,i])%*%senssum.B ) )
  return(variance)
}

require("pbmcapply")
Births.cova.1 <- pbmclapply(as.list(1:250), function(x){
  Z.all <- list()
  for(i in 1:Time.steps){
    Z.all[[i]] <- c(Z.deaths[i], Z.births[i], result[[x]]$Kappa1[i], result[[x]]$Kappa2[i])
  }
  now <- Sys.time()
  Cova.estim <- Covariance.estimate.births(fit.B = Births.fit[[x]]$coef,
                                           Births = result[[x]]$Births, 
                                           Births.dummy = Births.fit[[x]]$dummy, 
                                           Z.births = Z.all,
                                           rho0 = Births.fit[[x]]$rho0,
                                           bw.seq = bw.seq,
                                           Time.steps = Time.steps)
  later <- Sys.time()
  time.vec <- later - now
  attr(Cova.estim, "time") <- time.vec
  return(Cova.estim)
  }, mc.cores = 16)

result <- result[-c(1:250)]
Births.fit <- Births.fit[-c(1:250)]

Births.cova.2 <- pbmclapply(as.list(1:250), function(x){
  Z.all <- list()
  for(i in 1:Time.steps){
    Z.all[[i]] <- c(Z.deaths[i], Z.births[i], result[[x]]$Kappa1[i], result[[x]]$Kappa2[i])
  }
  now <- Sys.time()
  Cova.estim <- Covariance.estimate.births(fit.B = Births.fit[[x]]$coef,
                                           Births = result[[x]]$Births, 
                                           Births.dummy = Births.fit[[x]]$dummy, 
                                           Z.births = Z.all,
                                           rho0 = Births.fit[[x]]$rho0,
                                           bw.seq = bw.seq,
                                           Time.steps = Time.steps)
  later <- Sys.time()
  time.vec <- later - now
  attr(Cova.estim, "time") <- time.vec
  return(Cova.estim)
  }, mc.cores = 16)

result <- result[-c(1:250)]
Births.fit <- Births.fit[-c(1:250)]

Births.cova.3 <- pbmclapply(as.list(1:250), function(x){
  Z.all <- list()
  for(i in 1:Time.steps){
    Z.all[[i]] <- c(Z.deaths[i], Z.births[i], result[[x]]$Kappa1[i], result[[x]]$Kappa2[i])
  }
  now <- Sys.time()
  Cova.estim <- Covariance.estimate.births(fit.B = Births.fit[[x]]$coef,
                             Births = result[[x]]$Births, 
                             Births.dummy = Births.fit[[x]]$dummy, 
                             Z.births = Z.all,
                             rho0 = Births.fit[[x]]$rho0,
                             bw.seq = bw.seq,
                             Time.steps = Time.steps)
  later <- Sys.time()
  time.vec <- later - now
  attr(Cova.estim, "time") <- time.vec
  return(Cova.estim)
  }, mc.cores = 16)

result <- result[-c(1:250)]
Births.fit <- Births.fit[-c(1:250)]

Births.cova.4 <- pbmclapply(as.list(1:250), function(x){
  Z.all <- list()
  for(i in 1:Time.steps){
    Z.all[[i]] <- c(Z.deaths[i], Z.births[i], result[[x]]$Kappa1[i], result[[x]]$Kappa2[i])
  }
  now <- Sys.time()
  Cova.estim <- Covariance.estimate.births(fit.B = Births.fit[[x]]$coef,
                             Births = result[[x]]$Births, 
                             Births.dummy = Births.fit[[x]]$dummy, 
                             Z.births = Z.all,
                             rho0 = Births.fit[[x]]$rho0,
                             bw.seq = bw.seq,
                             Time.steps = Time.steps)
  later <- Sys.time()
  time.vec <- later - now
  attr(Cova.estim, "time") <- time.vec
  return(Cova.estim)
  }, mc.cores = 16)

Births.cova <- c(Births.cova.1,
                 Births.cova.2,
                 Births.cova.3,
                 Births.cova.4)

saveRDS(Births.cova,  paste(experiment, "/Births_Cova_N0_",N0,"_",prefix,"_",setting$Kappa_fun,"_all_biv.RDS", sep = "") )

Time.val <- lapply( Births.cova, function(x) attr(x, "time"))
saveRDS(Time.val,  paste(experiment, "/Births_Cova_N0_",N0,"_",prefix,"_",setting$Kappa_fun,"_Time.RDS", sep = "") )
