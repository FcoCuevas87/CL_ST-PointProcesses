#01_Fitting_model
require("spatstat")
require("pbmcapply")

setwd("/data/fcp-X79-03/Spatio_temporal_PP/stBCIcode/SimulationStudy")
source("CE_Method.R")

settings <- readRDS( paste(experiment,"/SimSetting_",N0,"_",prefix,"_all.RDS", sep = ""))
W.sim <- settings$window
Time.steps <- settings$Time.steps
N0 <- settings$N0
source("01_Simulation_Covariates.R")
result <- readRDS( paste(experiment,"/BDH_PP_",N0,"_",prefix,"_",settings$Kappa_fun,"_all_biv.RDS", sep = "") )

#   Fitting the deaths via logistic regression (clean the code)
#(Death = x$Death, Surv = x$Surv, Z.deaths = Z.all, Time.steps)
fitting_deaths <- function(Death, Surv, Z.deaths, Time.steps, year = FALSE){
  
  eval.covariate <- function(i, deaths, year = FALSE){
    n0 <- length(Z.deaths[[i]])
    if(deaths){
      eval.list <- matrix(0, ncol = n0, nrow = Death[[i]]$n)
    }else{
      eval.list <- matrix(0, ncol = n0, nrow = Surv[[i]]$n)
    }
    for(jj in 1:n0){
      if(deaths){
        eval.list[,jj] <- Z.deaths[[i]][[jj]][ Death[[i]] ]
      }else{
        eval.list[,jj] <- Z.deaths[[i]][[jj]][ Surv[[i]] ]
      }
    }
    if(year) eval.list <- cbind(eval.list, i)
    return(eval.list)
  }  
  Deaths.cova <- lapply( as.list(1:Time.steps), eval.covariate, deaths = TRUE, year)
  Survivals.cova <- lapply( as.list(1:Time.steps), eval.covariate, deaths = FALSE, year)
  Cova.Deaths <- do.call("rbind", Deaths.cova)
  Cova.Survivals <- do.call("rbind", Survivals.cova)
  Flag.vec.D <- c(rep(1,nrow(Cova.Deaths)), rep(0,nrow(Cova.Survivals)))
  values.D <- factor( Flag.vec.D )
  Cova.D <- rbind(Cova.Deaths, Cova.Survivals)

  if(year){
    Year.factor <- factor(Cova.D[,ncol(Cova.D)])
    Cova.D <- Cova.D[,-ncol(Cova.D)]
    fit.D <- glm(values.D  ~ Cova.D + Year.factor + 0, family=binomial)
  }else{
    fit.D <- glm(values.D  ~ Cova.D, family=binomial)
  }
  return(fit.D)
}

fitting_births <- function(Births, Z.births, Time.steps){
  # Fitting the Births via logistic regression
  
  rho0 <- 4*max( do.call("c", lapply(Births, function(x) x$n)) )/spatstat.geom::area(Births[[1]]$window)
  Dummy.sim <- lapply( as.list(1:Time.steps), function(i){
    W0 <- Births[[i]]$window
    X0 <- rpoispp(rho0, win = W0)
    return(X0)
  })
  
  eval.covariate <- function(i, births){
    n0 <- length(Z.births[[i]])
    if(births){
      eval.list <- matrix(0, ncol = n0, nrow = Births[[i]]$n)
    }else{
      eval.list <- matrix(0, ncol = n0, nrow = Dummy.sim[[i]]$n)
    }
    
    for(jj in 1:n0){
      if(births){
        eval.list[,jj] <- Z.births[[i]][[jj]][ Births[[i]] ]
      }else{
        eval.list[,jj] <- Z.births[[i]][[jj]][ Dummy.sim[[i]] ]
      }
    }
    return(eval.list)
  }
  
  Births.val <- lapply(as.list(1:Time.steps), eval.covariate, births = TRUE) 
  Births.val <- do.call("rbind", Births.val)
  
  Dummy.val <- lapply(as.list(1:Time.steps), eval.covariate, births = FALSE) 
  Dummy.val <- do.call("rbind", Dummy.val)
  
  Cova.B <- rbind(Births.val, Dummy.val)
  
  Flag.vec.B <- rep(c(1,0), c(nrow(Births.val), nrow(Dummy.val)))
  values.B <- factor( Flag.vec.B )
  
  fit.B <- glm(values.B  ~ Cova.B + offset(rep( -log(rho0), length(Flag.vec.B) )), family=binomial)
  return( list(coef = fit.B$coefficients, dummy = Dummy.sim, rho0 = rho0) )
}

require("parallel")

Deaths.coef <- pbmclapply(result, function(x){
  Z.all <- list()
  Death.all <- NULL
  Surv.all <- NULL
  for(kk in 1:Time.steps){
    Z.all[[kk]] <- c(Z.deaths[kk],Z.births[kk], x$Kappa3[kk], x$Kappa4[kk])
    X.tmp <- superimpose(x$Death[[kk]], x$Surv[[kk]])
    #fun.val <- knn.val(X.tmp)
    #Death.all <- c(Death.all, fun.val[1:x$Death[[kk]]$n])
    #Surv.all <- c(Surv.all, fun.val[(x$Death[[kk]]$n + 1):X.tmp$n])
  }
  now <- Sys.time()
  coef <- fitting_deaths(x$Death, x$Surv, Z.all, Time.steps, year = FALSE)$coefficients
  later <- Sys.time()
  time.vec <- later - now
  attr(coef, "time") <- time.vec
  return(coef)
  #fitting_deaths(x$Death, x$Surv, Z.all, NULL, Time.steps)
}, mc.cores = 36)

time.list.death <- lapply(Deaths.coef, function(x) attr(x, "time") )
Deaths.coef <- do.call("rbind", Deaths.coef)
#hist(Deaths.coef[,1]);abline(v = -0.25, col = "blue", lwd = 2); abline(v = mean(Deaths.coef[,1]), col = "red", lwd = 2)
#hist(Deaths.coef[,2]);abline(v = 0.25, col = "blue", lwd = 2); abline(v = mean(Deaths.coef[,2]), col = "red", lwd = 2)
#hist(Deaths.coef[,3]);abline(v = 0, col = "blue", lwd = 2); abline(v = mean(Deaths.coef[,3]), col = "red", lwd = 2)
#hist(Deaths.coef[,4]);abline(v = 0, col = "blue", lwd = 2); abline(v = mean(Deaths.coef[,4]), col = "red", lwd = 2)
#hist(Deaths.coef[,5]);abline(v = -0.25, col = "blue", lwd = 2); abline(v = mean(Deaths.coef[,5]), col = "red", lwd = 2)
cor(Deaths.coef)

saveRDS(Deaths.coef,  paste(experiment,"/Deaths_fit_",N0,"_",prefix,"_",settings$Kappa_fun,"_all_biv.RDS", sep = "")  )

Births.fit <- pbmclapply(result, function(x){
  Z.all <- NULL
  for(kk in 1:Time.steps){
    Z.all[[kk]] <- c(Z.deaths[kk],Z.births[kk], x$Kappa1[kk], x$Kappa2[kk])
  }
  now <- Sys.time()
  coef <- fitting_births(x$Births, Z.all, Time.steps)
  later <- Sys.time()
  time.vec <- later - now
  attr(coef, "time") <- time.vec
  return(coef)
}, mc.cores = 25)


time.list.births <- lapply(Births.fit, function(x) attr(x, "time") )
Births.coef <- do.call("rbind", lapply(Births.fit, function(x) x$coef ) )
births.mean <- colMeans(Births.coef)
saveRDS(Births.fit,  paste(experiment,"/Births_fit_",N0,"_",prefix,"_",settings$Kappa_fun,"_all_biv.RDS", sep = "")  )

#hist(Births.coef[,1]); abline(v = settings$beta.vec[1], col = "blue", lwd = 2); abline(v = births.mean[1], col = "red", lwd = 2)
#hist(Births.coef[,2]); abline(v = 0, col = "blue", lwd = 2); abline(v = births.mean[2], col = "red", lwd = 2)
#hist(Births.coef[,3]); abline(v = settings$beta.vec[2], col = "blue", lwd = 2); abline(v = births.mean[3], col = "red", lwd = 2)
#hist(Births.coef[,4]); abline(v = settings$beta.vec[3], col = "blue", lwd = 2); abline(v = births.mean[4], col = "red", lwd = 2)
#hist(Births.coef[,5]); abline(v = -2, col = "blue", lwd = 2); abline(v = births.mean[5], col = "red", lwd = 2)
#cor(Births.coef)


death.coef.time <- do.call("c", time.list.death)
births.coef.time <- do.call("c", time.list.births)
saveRDS(cbind(death.coef.time, births.coef.time), paste(experiment,"/Both_fit_",N0,"_",prefix,"_",settings$Kappa_fun,"_time_estim.RDS", sep = ""))

