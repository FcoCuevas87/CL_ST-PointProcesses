rm(list = ls())
load("Setting_W.Rdat")

library("fmesher")
library("INLA")
inlafitnew <- function(B, # list of point patterns of births
                       D, # list of point patterns of deaths
                       S, # list of point patterns of survivors
                       covars, # list of spatial covariates
                       st.covars, # list of spatiotemporal influence covariates
                       sp="cappfr", # selected species
                       dimyx=c(50, 100), # dimension of discretization for the birth model
                       num.threads=2 # number of CPU core for multi-core processing
)
{
  # common observation window
  W <- B[[1]]$window
  # number of time instances
  K <- length(B)
  if (length(D) != K | length(S) != K)
    stop("all lists of point patterns must have the same length")
  
  # prepare the data
  # discretization for the birth model
  xx <- seq(W$xrange[1], W$xrange[2], l=dimyx[2])
  yy <- seq(W$yrange[2], W$yrange[1], l=dimyx[1])
  xx <- rep(xx, each=dimyx[1])
  yy <- rep(yy, dimyx[2])
  ddat <- bdat <- NULL
  for (k in 1:K){
    # deaths
    DSk <- superimpose(D[[k]], S[[k]]) # superimpose deaths and survivors
    ddat <- rbind(ddat, data.frame(k=k, 
                                   status=c(rep(1, D[[k]]$n) , rep(0, S[[k]]$n)),
                                   x=DSk$x, 
                                   y=DSk$y,
                                   stB1=st.covars$stB[k, 1][[1]][DSk],
                                   stB2=st.covars$stB[k, 2][[1]][DSk],
                                   stD1=st.covars$stD[k, 1][[1]][DSk],
                                   stD2=st.covars$stD[k, 2][[1]][DSk]))
    
    # births 
    zz <- quadratcount(B[[k]], nx=dimyx[2], ny=dimyx[1])
    bdat <- rbind(bdat, data.frame(k=k, counts=c(zz), x=xx, y=yy,
                                   stB1=st.covars$stB[k, 1][[1]][list(x = xx, y = yy)],
                                   stB2=st.covars$stB[k, 2][[1]][list(x = xx, y = yy)],
                                   stD1=st.covars$stD[k, 1][[1]][list(x = xx, y = yy)],
                                   stD2=st.covars$stD[k, 2][[1]][list(x = xx, y = yy)]))
  }
  
  bdat <- cbind(bdat, lapply(covars, function(o){ o[ list(x = bdat$x, y = bdat$y) ] }) )
  ddat <- cbind(ddat, lapply(covars, function(o){ o[ list(x = ddat$x, y = ddat$y) ] }) )
  bdat$time <- as.factor(bdat$k)
  ddat$time <- as.factor(ddat$k)
  
  # death model
  now.deaths <- Sys.time()
  dmesh <- fm_mesh_2d_inla(loc=cbind(ddat$x, ddat$y))
  plot(dmesh)
  dspde <- inla.spde2.matern(mesh=dmesh, alpha=2, constr=TRUE)
  dform <- as.formula(paste(
    "as.integer(status) ~ 1 + f(dmesh$idx$loc, model=dspde, 
                   group=k, control.group=list(model='ar1')) + ",
    paste(names(covars), collapse = "+"), 
    "+ stD1 + stD2"))
  dfit <- inla(dform, data=ddat, family="binomial", num.threads=num.threads)
  d <- dfit$summary.fixed[, c(1:3, 5)]
  
  # estimated spatiotemporal random effect
  dmpj <- inla.mesh.projector(dmesh,  dims=rev(dimyx))
  dm <- list()
  for (k in 1:K){
    dm[[k]] <- inla.mesh.project(dmpj, 
                                 dfit$summary.random$`dmesh|S|idx|S|loc`$mean[(k - 1) * dmesh$n + 1:dmesh$n])
    dm[[k]] <- as.im(list(x=dmpj$x, y=dmpj$y, z=dm[[k]]), W=W)
  }
  
  later.deaths <- Sys.time()
  
  # birth model
  now.births <- Sys.time()
  bmesh <- fm_mesh_2d_inla(loc=cbind(bdat$x, bdat$y))
  plot(bmesh)
  bspde <- inla.spde2.matern(mesh=bmesh, alpha=2, constr=TRUE)
  bform <- as.formula(paste(
    "counts ~ 1 + f(bmesh$idx$loc, model=bspde, 
                   group=k, control.group=list(model='ar1')) + ",
    paste(names(covars), collapse = "+"), 
    "+ stB1 + stB2"))
  bfit <- inla(bform, data=bdat, family="poisson", num.threads=num.threads)
  b <- bfit$summary.fixed[, c(1:3, 5)]
  
  # estimated spatiotemporal random effect
  bmpj <- inla.mesh.projector(bmesh,  dims=rev(dimyx))
  bm <- list()
  for (k in 1:K)
  {
    bm[[k]] <- inla.mesh.project(bmpj, 
                                 bfit$summary.random$`bmesh|S|idx|S|loc`$mean[(k - 1) * bmesh$n + 1:bmesh$n])
    bm[[k]] <- as.im(list(x=bmpj$x, y=bmpj$y, z=bm[[k]]), W=W)
  }
  later.births <- Sys.time()
  
  return(list(beta=rbind(births = b, deaths = d), dfit=dfit, bfit=bfit, 
              btheta=inla.spde2.result(bfit, "bmesh|S|idx|S|loc", bspde),
              dtheta=inla.spde2.result(dfit, "dmesh|S|idx|S|loc", dspde),
              dm=dm, bm=bm, time.births = later.births - now.births,
              time.deaths = later.deaths - now.deaths))
}

#04_Births_Variance.R

require("pbmcapply")
require("spatstat")
setwd("/data/fcp-X79-03/Spatio_temporal_PP/stBCIcode/SimulationStudy")
settings <- readRDS( paste(experiment,"/SimSetting_",N0,"_",prefix,"_all.RDS", sep = ""))
W.sim <- settings$window
Time.steps <- settings$Time.steps
N0 <- settings$N0
source("01_Simulation_Covariates.R")
result <- readRDS( paste(experiment,"/BDH_PP_",N0,"_",prefix,"_",settings$Kappa_fun,"_all_biv.RDS", sep = "") )

result <- result[1:100]

time.vec <- matrix(0, nrow = 100, ncol = 2)
inla.result.list <- list()
i <- 1
for(i in 1:100){
  now <- Sys.time()
  inla.result <- inlafitnew(B = result[[i]]$Births, # list of point patterns of births
                         D = result[[i]]$Death, # list of point patterns of deaths
                         S = result[[i]]$Surv, # list of point patterns of survivors
                         covars = list(births = Z.births[[1]], 
                                       deaths = Z.deaths[[1]]), # list of spatial covariates
                         st.covars = list(stB = cbind(result[[i]]$Kappa1, result[[i]]$Kappa2), 
                                          stD = cbind(result[[i]]$Kappa3, result[[i]]$Kappa4)), # list of spatiotemporal influence covariates
                         dimyx= c(50, 100), # dimension of discretization for the birth model
                         num.threads=50) # number of CPU core for multi-core processing
  later <- Sys.time()
  time.vec[i,1] <- inla.result$time.births
  time.vec[i,2] <- inla.result$time.deaths
  inla.result.list[[i]] <- inla.result$beta
  cat(i, "\n")
}


saveRDS(inla.result.list, "inla_results.RDS")
saveRDS(time.vec, "time_vec.RDS")
