#03_Deaths_Variance.R
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

result <- readRDS( paste(experiment, "/BDH_PP_",N0,"_",prefix,"_",setting$Kappa_fun,"_all_biv.RDS", sep = "") )
Deaths.coef <- readRDS( paste(experiment, "/Deaths_fit_",N0,"_",prefix,"_",setting$Kappa_fun,"_all_biv.RDS", sep = "") )

kernel <- function(h,bw){
  return(as.numeric(h<bw))#simple uniform kernel used here.
}

Covariance.estimate.deaths <- function(fit.D, Survs, Deaths, Z.deaths, 
                                       Z.vec.d = NULL, Z.vec.s = NULL, bw.seq, Time.steps){
  
  indexing.s <- unlist( lapply(Survs, function(x) x$n) ) 
  indexing.d <- unlist( lapply(Deaths, function(x) x$n) )
  
  indexing.s <- cumsum(indexing.s)
  indexing.d <- cumsum(indexing.d)
  
  indexing.d <- cbind( c(1, indexing.d[1:(Time.steps-1)] + 1), indexing.d )
  indexing.s <- cbind( c(1, indexing.s[1:(Time.steps-1)] + 1), indexing.s )
  
  eval.covariate <- function(i, deaths){
    n0 <- length(Z.deaths[[i]])
    if(deaths){
      eval.list <- matrix(0, ncol = n0, nrow = Deaths[[i]]$n)
    }else{
      eval.list <- matrix(0, ncol = n0, nrow = Survs[[i]]$n)
    }
    
    for(jj in 1:n0){
      if(deaths){
        eval.list[,jj] <- Z.deaths[[i]][[jj]][ Deaths[[i]] ]
      }else{
        eval.list[,jj] <- Z.deaths[[i]][[jj]][ Survs[[i]] ]
      }
    }
    return(eval.list)
  }
  
  Sens.Var <- lapply(as.list(1:Time.steps), function(i){
    ### Deaths
    sensiv.D <- list()
    var.D <- list()
    Xkminus1 <- superimpose(Survs[[i]], Deaths[[i]])
    
    Deaths.cova.vec <- eval.covariate(i, deaths = TRUE)
    Survivals.cova.vec <- eval.covariate(i, deaths = FALSE)
    
    cova.D <- rbind(Survivals.cova.vec, Deaths.cova.vec)
    
    if( !is.null(Z.vec.d) & !is.null(Z.vec.s) ) cova.D <- cbind(cova.D, c(Z.vec.s[indexing.s[i,1]:indexing.s[i,2]], 
                                                                          Z.vec.d[indexing.d[i,1]:indexing.d[i,2]]) )
    
    #Flag value
    Flag.D <- rep( c(0,1), c(Survs[[i]]$n, Deaths[[i]]$n))
    grad.etha.k <- cbind(1, cova.D)
    estimated.etha.k <- c(fit.D%*%t(grad.etha.k)) #use the true as well
    
    sensiv.D <- sensi.death.k(grad.eta.k = grad.etha.k, eta.k = estimated.etha.k)
    
    var.D.list <- e.death.var.k.alt2.bw(grad.eta.k = grad.etha.k,
                                        eta.k = estimated.etha.k,
                                        bw = bw.seq,
                                        Ik = Flag.D,
                                        Xkminus1 = Xkminus1)
    
    var.D <- simplify2array(var.D.list)
    
    return(list(Sens = sensiv.D, Var = var.D))
    #cat(i,"\n")
  })
  
  sensiv.D <- lapply(as.list(1:Time.steps), function(i) Sens.Var[[i]]$Sens)
  var.D <- lapply(Sens.Var, function(x) x$Var)
  
  senssum.D <- Reduce("+", sensiv.D)
  varsum.D <- Reduce("+", var.D)
  variance <- lapply(as.list(1:length(bw.seq)), function(i) solve( senssum.D%*%solve(varsum.D[,,i])%*%senssum.D ) )
  return(variance)
}

# x <- result[[1]]
# Z.all <- list()
# Death.all <- NULL
# Surv.all <- NULL
# for(kk in 1:Time.steps){
#   Z.all[[kk]] <- c(Z.deaths[kk],Z.births[kk])
#   X.tmp <- superimpose(x$Death[[kk]], x$Surv[[kk]])
#   fun.val <- knn.val(X.tmp)
#   Death.all <- c(Death.all, fun.val[1:x$Death[[kk]]$n])
#   Surv.all <- c(Surv.all, fun.val[(x$Death[[kk]]$n + 1):X.tmp$n])
# }

Deaths.cova.list.1 <- pbmclapply(as.list(1:250),function(x){
  Z.all <- list()
  Death.all <- NULL
  Surv.all <- NULL
  for(i in 1:Time.steps){
    Z.all[[i]] <- c(Z.deaths[i],Z.births[i], result[[x]]$Kappa3[i], result[[x]]$Kappa4[i])
    X.tmp <- superimpose(result[[x]]$Death[[i]], result[[x]]$Surv[[i]])
    #fun.val <- knn.val(X.tmp)
    #Death.all <- c(Death.all, fun.val[1:result[[x]]$Death[[i]]$n])
    #Surv.all <- c(Surv.all, fun.val[(result[[x]]$Death[[i]]$n + 1):X.tmp$n])
  }
  now <- Sys.time()
  Cova.estim <- Covariance.estimate.deaths(fit.D = Deaths.coef[x,], 
                             Survs = result[[x]]$Surv, 
                             Deaths = result[[x]]$Death, 
                             Z.deaths = Z.all,
                             Z.vec.d = Death.all,
                             Z.vec.s = Surv.all,
                             bw.seq = bw.seq,
                             Time.steps = Time.steps)
  later <- Sys.time()
  time.vec <- later - now
  attr(Cova.estim, "time") <- time.vec
  return(Cova.estim)
}, mc.cores = 30)

result <- result[-c(1:250)]
Deaths.coef <- Deaths.coef[-c(1:250),]

Deaths.cova.list.2 <- pbmclapply(as.list(1:250),function(x){
  Z.all <- list()
  Death.all <- NULL
  Surv.all <- NULL
  for(i in 1:Time.steps){
    Z.all[[i]] <- c(Z.deaths[i],Z.births[i], result[[x]]$Kappa3[i], result[[x]]$Kappa4[i])
    X.tmp <- superimpose(result[[x]]$Death[[i]], result[[x]]$Surv[[i]])
    #fun.val <- knn.val(X.tmp)
    #Death.all <- c(Death.all, fun.val[1:result[[x]]$Death[[i]]$n])
    #Surv.all <- c(Surv.all, fun.val[(result[[x]]$Death[[i]]$n + 1):X.tmp$n])
  }
  now <- Sys.time()
  Cova.estim <- Covariance.estimate.deaths(fit.D = Deaths.coef[x,], 
                                           Survs = result[[x]]$Surv, 
                                           Deaths = result[[x]]$Death, 
                                           Z.deaths = Z.all,
                                           Z.vec.d = Death.all,
                                           Z.vec.s = Surv.all,
                                           bw.seq = bw.seq,
                                           Time.steps = Time.steps)
  later <- Sys.time()
  time.vec <- later - now
  attr(Cova.estim, "time") <- time.vec
  return(Cova.estim)
}, mc.cores = 30)

result <- result[-c(1:250)]
Deaths.coef <- Deaths.coef[-c(1:250),]

Deaths.cova.list.3 <- pbmclapply(as.list(1:250),function(x){
  Z.all <- list()
  Death.all <- NULL
  Surv.all <- NULL
  for(i in 1:Time.steps){
    Z.all[[i]] <- c(Z.deaths[i],Z.births[i], result[[x]]$Kappa3[i], result[[x]]$Kappa4[i])
    X.tmp <- superimpose(result[[x]]$Death[[i]], result[[x]]$Surv[[i]])
    #fun.val <- knn.val(X.tmp)
    #Death.all <- c(Death.all, fun.val[1:result[[x]]$Death[[i]]$n])
    #Surv.all <- c(Surv.all, fun.val[(result[[x]]$Death[[i]]$n + 1):X.tmp$n])
  }
  now <- Sys.time()
  Cova.estim <- Covariance.estimate.deaths(fit.D = Deaths.coef[x,], 
                                           Survs = result[[x]]$Surv, 
                                           Deaths = result[[x]]$Death, 
                                           Z.deaths = Z.all,
                                           Z.vec.d = Death.all,
                                           Z.vec.s = Surv.all,
                                           bw.seq = bw.seq,
                                           Time.steps = Time.steps)
  later <- Sys.time()
  time.vec <- later - now
  attr(Cova.estim, "time") <- time.vec
  return(Cova.estim)
}, mc.cores = 30)

result <- result[-c(1:250)]
Deaths.coef <- Deaths.coef[-c(1:250),]

Deaths.cova.list.4 <- pbmclapply(as.list(1:250),function(x){
  Z.all <- list()
  Death.all <- NULL
  Surv.all <- NULL
  for(i in 1:Time.steps){
    Z.all[[i]] <- c(Z.deaths[i],Z.births[i], result[[x]]$Kappa3[i], result[[x]]$Kappa4[i])
    X.tmp <- superimpose(result[[x]]$Death[[i]], result[[x]]$Surv[[i]])
    #fun.val <- knn.val(X.tmp)
    #Death.all <- c(Death.all, fun.val[1:result[[x]]$Death[[i]]$n])
    #Surv.all <- c(Surv.all, fun.val[(result[[x]]$Death[[i]]$n + 1):X.tmp$n])
  }
  now <- Sys.time()
  Cova.estim <- Covariance.estimate.deaths(fit.D = Deaths.coef[x,], 
                                           Survs = result[[x]]$Surv, 
                                           Deaths = result[[x]]$Death, 
                                           Z.deaths = Z.all,
                                           Z.vec.d = Death.all,
                                           Z.vec.s = Surv.all,
                                           bw.seq = bw.seq,
                                           Time.steps = Time.steps)
  later <- Sys.time()
  time.vec <- later - now
  attr(Cova.estim, "time") <- time.vec
  return(Cova.estim)
}, mc.cores = 30)


Deaths.cova.list <- c(Deaths.cova.list.1,
                      Deaths.cova.list.2,
                      Deaths.cova.list.3,
                      Deaths.cova.list.4)

saveRDS(Deaths.cova.list,  paste( experiment,"/Deaths_Cova_N0_",N0,"_",prefix,"_",setting$Kappa_fun,"_all_biv.RDS", sep = "") )

Time.val <- lapply( Deaths.cova.list, function(x) attr(x, "time"))
saveRDS( do.call("c", Time.val),  paste( experiment,"/Deaths_Cova_N0_",N0,"_",prefix,"_",setting$Kappa_fun,"_Time.RDS", sep = "") ) 

        