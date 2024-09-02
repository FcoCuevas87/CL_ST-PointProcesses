
### Multivariate LGCP 
library("spatstat")
library("GeoModels")
source("CE_Method.R") #Circulant embedding with GeoModels
source("Kappa_fun.R")

#SubWindow
if(my.win == "Window1"){ W.sim <- owin(xrange = c(-2.5,500), yrange = c(-2.5, 250)); prefix <- "W1"}
if(my.win == "Window2"){ W.sim <- owin(xrange = c(-2.5,1000), yrange = c(-2.5, 500)); prefix <- "W2"}

######## Reading covariates ########
source("01_Simulation_Covariates.R")

##########################################
### DHB with LGCP only a point pattern ###
##########################################

param.B <- list(mean = 0, smooth = 1.75, scale = 4, sill  = 1, nugget = 0)
param.D <- list(mean = 0, smooth = 0.5, scale = 7, sill  = 1, nugget = 0)

seq(0,1000, l = 1000) -> x0
plot(x0, GeoCorrFct(x0, corrmodel = "Matern", param = param.B)$corr, type = "l")
x0[ which.min( abs(GeoCorrFct(x0, corrmodel = "Matern", param = param.B)$corr - 0.05)) ]
x0[ which.min( abs(GeoCorrFct(x0, corrmodel = "Matern", param = param.D)$corr - 0.05)) ]

kappa.bw <- 6.0
kappa.bw.death <- 15.0

#seq(0,1000, l = 1000) -> x0
#plot(x0, GeoCorrFct(x0, corrmodel = "Matern", param = param.D), type = "l")
#x0[ which.min( abs(GeoCorrFct(x0, corrmodel = "Matern", param = param.D) - 0.05)) ]

beta1b <- 0.1
beta.kappa <- 0.1
if(my.win == "Window2"){
  beta0b <- log(N0/integral( exp(beta1b*Z.births[[1]][W.sim])) ) #Initial number of points
}

if(my.win == "Window1") beta0b <- readRDS( paste(experiment, "/SimSetting_",N0,"_W2_all.RDS", sep = ""))$beta.vec[1]
nsim <- 1000

t0 <- sqrt( (W.sim$xrange[2] - W.sim$xrange[1])^2 + (W.sim$yrange[2] - W.sim$yrange[1])^2 )
dev.off()
library("foreach")
library("parallel")
library("doParallel")
etha.vec <- c(-0.25, +0.25)

saveRDS( list(beta.vec = c(beta0b, beta1b, beta.kappa), 
              etha.vec = etha.vec, 
              N0 = N0,
              param.B = param.B,
              param.D = param.D,
              kappa.bw = kappa.bw,
              bw.seq = bw.seq,
              Time.steps = Time.steps,
              window = W.sim,
              prefix = prefix,
              Kappa_fun = Kappa_select),  
         paste(experiment,"/SimSetting_",N0,"_",prefix,"_all.RDS", sep = "") )

n_cores <- 35
cl <- makePSOCKcluster( n_cores )
registerDoParallel(cl)

result <- foreach(k0 = 1:nsim)%dopar%{
  require("spatstat")
  
  #Initial value
  sim.result.B.1 <- GRF.CE.Sim(nsim = 2, W.sim$xrange, W.sim$xrange, M = N1, N = N1, corrmodel = "Matern", param = param.B, mean.val = 0) 
  sim.result.B.2 <- GRF.CE.Sim(nsim = 2, W.sim$xrange, W.sim$xrange, M = N1, N = N1, corrmodel = "Matern", param = param.B, mean.val = 0) 
  
  RF.all.1 <- im( matrix(sim.result.B.1$X[,1] , ncol = N1, nrow = N1), xrange = W.sim$xrange, yrange = W.sim$xrange)
  RF.all.2 <- im( matrix(sim.result.B.2$X[,1] , ncol = N1, nrow = N1), xrange = W.sim$xrange, yrange = W.sim$xrange)
  
  expand.grid <- expand.grid(seq(W.sim$xrange[1], W.sim$xrange[2], l = M), seq(W.sim$yrange[1], W.sim$yrange[2], l = N))
  RF.ini.vec.1 <- RF.all.1[expand.grid]
  RF.ini.vec.2 <- RF.all.2[expand.grid]
  
  RF.ini.1 <- im( matrix(RF.ini.vec.1 , ncol = N, nrow = M), xrange = W.sim$xrange, yrange = W.sim$yrange)
  RF.ini.2 <- im( matrix(RF.ini.vec.2 , ncol = N, nrow = M), xrange = W.sim$xrange, yrange = W.sim$yrange)
  
  X0.t.1 <- rpoispp( exp(beta1b*Z.births[[1]] + beta0b + RF.ini.1 - 0.5*param.B$sill) )
  X0.t.1 <- X0.t.1[W.sim]
  X0.t.2 <- rpoispp( exp(beta1b*Z.births[[1]] + beta0b + RF.ini.2 - 0.5*param.B$sill) )
  X0.t.2 <- X0.t.2[W.sim]
  
  kappa.list.1 <- kappa.list.2 <- list() #kappa.list.d <- 
  kappa.list.1[[1]] <- Kappa_fun(X = X0.t.1, bw = kappa.bw, ox.seq = x, oy.seq = y, Kappa_select = "Kappa_2")
  kappa.list.2[[1]] <- Kappa_fun(X = X0.t.2, bw = kappa.bw, ox.seq = x, oy.seq = y, Kappa_select = "Kappa_2")
  #kappa.list.d[[1]] <- Kappa_fun(X = X0.t, bw = kappa.bw, ox.seq = x, oy.seq = y, Kappa_select = Kappa_select)
  
  X1.t <- X2.t <- list()
  
  RF.all.1 <- im( matrix(sim.result.B.1$X[,2] , ncol = N1, nrow = N1), xrange = W.sim$xrange, yrange = W.sim$xrange)
  RF.ini.vec.1 <- RF.all.1[expand.grid]
  RF.ini.1 <- im( matrix(RF.ini.vec.1 , ncol = N, nrow = M), xrange = W.sim$xrange, yrange = W.sim$yrange)

  RF.all.2 <- im( matrix(sim.result.B.2$X[,2] , ncol = N1, nrow = N1), xrange = W.sim$xrange, yrange = W.sim$xrange)
  RF.ini.vec.2 <- RF.all.2[expand.grid]
  RF.ini.2 <- im( matrix(RF.ini.vec.2 , ncol = N, nrow = M), xrange = W.sim$xrange, yrange = W.sim$yrange)
  
  #X1.t[[1]] <- rpoispp( exp(beta1b*Z.births[[1]] + beta0b + RF.ini[W.sim] - 0.5*param.B$sill) )
  X1.t[[1]] <- rpoispp( exp(beta1b*Z.births[[1]] + beta0b + beta.kappa*kappa.list.1[[1]] + beta.kappa*kappa.list.2[[1]] 
                            + RF.ini.1[W.sim] - 0.5*param.B$sill) )
  X1.t[[1]] <- X1.t[[1]][W.sim]
  
  X2.t[[1]] <- rpoispp( exp(beta1b*Z.births[[1]] + beta0b + beta.kappa*kappa.list.2[[1]] + beta.kappa*kappa.list.1[[1]]
                            + RF.ini.2[W.sim] - 0.5*param.B$sill) )
  X2.t[[1]] <- X2.t[[1]][W.sim]
  #Simulating Noise
  sim.result.B.1 <- GRF.CE.Sim(nsim = Time.steps, W.sim$xrange, W.sim$xrange, M = N1, N = N1, corrmodel = "Matern", param = param.B, mean.val = 0) 
  sim.result.B.2 <- GRF.CE.Sim(nsim = Time.steps, W.sim$xrange, W.sim$xrange, M = N1, N = N1, corrmodel = "Matern", param = param.B, mean.val = 0) 
  
  Z.noise.B.1 <- Z.noise.B.2 <- list()
  Z.noise.B.1 <- lapply(1:Time.steps, function(ll){ 
    RF.test <- im( matrix(sim.result.B.1$X[,ll] , ncol = N1, nrow = N1), xrange = W.sim$xrange, yrange = W.sim$xrange)
    RF.ini.vec <- RF.test[expand.grid]
    RF.test <- im( matrix(RF.ini.vec , ncol = N, nrow = M), xrange = W.sim$xrange, yrange = W.sim$yrange)
    return(RF.test)
  })
  Z.noise.B.2 <- lapply(1:Time.steps, function(ll){ 
    RF.test <- im( matrix(sim.result.B.2$X[,ll] , ncol = N1, nrow = N1), xrange = W.sim$xrange, yrange = W.sim$xrange)
    RF.ini.vec <- RF.test[expand.grid]
    RF.test <- im( matrix(RF.ini.vec , ncol = N, nrow = M), xrange = W.sim$xrange, yrange = W.sim$yrange)
    return(RF.test)
  })
  
  sim.result.D.1 <- GRF.CE.Sim(nsim = Time.steps, W.sim$xrange, W.sim$xrange, M = N1, N = N1, corrmodel = "Matern", param = param.D, mean.val = 0) 
  sim.result.D.2 <- GRF.CE.Sim(nsim = Time.steps, W.sim$xrange, W.sim$xrange, M = N1, N = N1, corrmodel = "Matern", param = param.D, mean.val = 0) 
  Z.noise.D.1 <- Z.noise.D.2 <- list()
  
  Z.noise.D.1 <- lapply(1:Time.steps, function(ll){ 
    RF.test <- im( matrix(sim.result.D.1$X[,ll] , ncol = N1, nrow = N1), xrange = W.sim$xrange, yrange = W.sim$xrange)
    RF.ini.vec <- RF.test[expand.grid]    
    RF.test <- im( matrix(RF.ini.vec , ncol = N, nrow = M), xrange = W.sim$xrange, yrange = W.sim$yrange)
    return(RF.test)
  })
  Z.noise.D.2 <- lapply(1:Time.steps, function(ll){ 
    RF.test <- im( matrix(sim.result.D.2$X[,ll] , ncol = N1, nrow = N1), xrange = W.sim$xrange, yrange = W.sim$xrange)
    RF.ini.vec <- RF.test[expand.grid]    
    RF.test <- im( matrix(RF.ini.vec , ncol = N, nrow = M), xrange = W.sim$xrange, yrange = W.sim$yrange)
    return(RF.test)
  })
  
  success <- TRUE
  while(success){
    Surv.1 <- Surv.2 <- list()
    Death.1 <- Death.2  <- list() 
    Births.1 <- Births.2 <- list() 
    kappa.list.1 <- kappa.list.2 <- list()
    kappa.list.3 <- kappa.list.4 <- list()
    for(i in 1:Time.steps){
      X.previous.1 <- X1.t[[i]]
      X.previous.2 <- X2.t[[i]]
      
      #stochastic covariable
      kappa.list.1[[i]] <- Kappa_fun(X = X1.t[[i]], bw = kappa.bw, ox.seq = x, oy.seq = y, Kappa_select = "Kappa_2")
      kappa.list.2[[i]] <- Kappa_fun(X = X2.t[[i]], bw = kappa.bw, ox.seq = x, oy.seq = y, Kappa_select = "Kappa_2")
      
      X1.mark <- X1.t[[i]]
      X2.mark <- X2.t[[i]]
      marks(X1.mark) <- rep(1,X1.t[[i]]$n)
      marks(X2.mark) <- rep(1,X2.t[[i]]$n)
      kappa.list.3[[i]] <- Kappa_fun_marks_sim(X = X1.mark, psi = kappa.bw.death, ox.seq = x, oy.seq = y)
      kappa.list.4[[i]] <- Kappa_fun_marks_sim(X = X2.mark, psi = kappa.bw.death, ox.seq = x, oy.seq = y)    
      #Deaths model
      etha1 <- Z.deaths[[i]][ X.previous.1 ]
      kappa.1 <- kappa.list.3[[i]][ X.previous.1 ]
      
      etha2 <- Z.deaths[[i]][ X.previous.2 ]
      kappa.2 <- kappa.list.4[[i]][ X.previous.2 ]
      
      p.mat.1 <- pnorm(Z.noise.D.1[[i]]$v) #+ 5*Z.deaths[[i-1]] + 2.5
      p.mat.2 <- pnorm(Z.noise.D.2[[i]]$v)
      
      Latent.1 <- as.im( log(p.mat.1/(1-p.mat.1)), W = W.sim)
      Latent.var.1 <- Latent.1[ X.previous.1 ]
      flag.1 <- as.numeric( I(Latent.var.1 <=  (etha.vec[2]*etha1 - 0.25*kappa.1 + etha.vec[1])) )
      
      Latent.2 <- as.im( log(p.mat.2/(1-p.mat.2)), W = W.sim)
      Latent.var.2 <- Latent.1[ X.previous.2 ]
      flag.2 <- as.numeric( I(Latent.var.2 <=  (etha.vec[2]*etha2 - 0.25*kappa.2 + etha.vec[1])) )

      marks(X.previous.1) <- as.factor(flag.1)
      marks(X.previous.2) <- as.factor(flag.2)
      splitted <- split(X.previous.1)
      Surv.1[[i]] <- splitted$'0'   ### Heritage ###
      Death.1[[i]] <- splitted$'1'  ### Deaths ###
      
      splitted <- split(X.previous.2)
      Surv.2[[i]] <- splitted$'0'   ### Heritage ###
      Death.2[[i]] <- splitted$'1'  ### Deaths ###
      #Births model

      Births.1[[i]] <- try( rpoispp( exp(beta1b*Z.births[[i]] +  beta.kappa*kappa.list.1[[i]] - 2*kappa.list.2[[i]] + beta0b +  Z.noise.B.1[[i]] - 0.5*param.B$sill )), silent = TRUE)  ### Births ###
      Births.2[[i]] <- try( rpoispp( exp(beta1b*Z.births[[i]] +  beta.kappa*kappa.list.2[[i]] - 2*kappa.list.1[[i]] + beta0b +  Z.noise.B.2[[i]] - 0.5*param.B$sill )), silent = TRUE)
      #Births[[i]] <- try( rpoispp( exp(beta1b*Z.births[[i]] + beta0b +  Z.noise.B[[i]] - 0.5*param.B$sill )), silent = TRUE)  ### Births ###
      Births.1[[i]] <- Births.1[[i]][W.sim]
      Births.2[[i]] <- Births.2[[i]][W.sim]
      if( inherits(Births.1[[i]], "try-error")){
        i <- 1
        cat("too many points!")
        break
      }
      if( log(Births.1[[i]]$n) >= 10){
        i <- 1
        cat("too many points 2!")
        break
      }
      X1.t[[i+1]] <- superimpose(Births.1[[i]], Surv.1[[i]])
      X2.t[[i+1]] <- superimpose(Births.2[[i]], Surv.2[[i]])
      #cat(i, ",", "B", Births[[i]]$n, "D", Death[[i]]$n, " n:", X1.t[[i]]$n, "\n")
    }
    
    if(i == Time.steps){
      success <- FALSE
    }
    cat(i, "\n")
  }
  
  N.seq <- do.call("c", lapply(X1.t, function(X) X$n))
  #N.seq.b <- do.call("c", lapply(Births, function(X) X$n))
  #N.seq.s <- do.call("c", lapply(Surv, function(X) X$n))
  #N.seq.d <- do.call("c", lapply(Death, function(X) X$n))
  
  return(list(Births = Births.1, Surv = Surv.1, Death = Death.1, Kappa1 = kappa.list.1, Kappa2 = kappa.list.2, Kappa3 = kappa.list.3, Kappa4 = kappa.list.4))
  }
stopCluster(cl)

saveRDS(result, paste(experiment, "/BDH_PP_",N0,"_",prefix,"_",Kappa_select,"_all_biv.RDS", sep = "") )

#unlist( lapply(Surv, function(x) x$n) )/(unlist( lapply(Surv, function(x) x$n) ) + unlist( lapply(Death, function(x) x$n) ))
#plot( X1.t[[10]], pch = 16)
#points( X2.t[[10]], pch = 16, col = "red")

