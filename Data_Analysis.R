rm(list = ls())
library("spatstat")
library("glmnet")
library("FNN")

#=====================================
#set seed for reproducibility
#=====================================
set.seed(123456)


#=========================================================
#source required functions
#=========================================================
source("Functions/asympvar_code.R")
source("Functions/asympvar_code_bw.R")
source("Main_Functions/Kappa_fun.R")
source("Main_Functions/Data_Analysis_prepdata.R")
source("Main_Functions/Data_Analysis_estm_cov.R")

#=======================================================
#load data
#=======================================================
load('Data/bci.spptable.rdata')
bci.spptable[bci.spptable$sp == "hybapr",]$Latin
# selected species
Xs <- BCIextract(c("cappfr"))
Xs.c <- BCIextract_complement(c("cappfr"))
table(Xs.c[[1]]$marks)


# ==============================================================
# read BCI covariate data
# ==============================================================
# BCI covariates
load("Data/bci.covars.Rda")
# adding scaled x and y coordinates as covariates
bci.covars$x <- as.im(function(x, y){ (x - 500) / 1000 },
                      W=Xs[[1]]$window, dimyx=c(101, 201))
bci.covars$y <- as.im(function(x, y){ (y - 250) / 500 },
                      W=Xs[[1]]$window, dimyx=c(101, 201))
par(mfrow=c(4, 4), mar=c(0.1, 0.1, 1.2, 2))
cvnames <- names(bci.covars)
for (i in 1:length(cvnames)) plot(bci.covars[[i]], main=cvnames[i])

# ==============================================================
# point patterns of recruits (births) and mortalities (deaths)
# ==============================================================
Time.steps <- 7
for(i in 1:(Time.steps+1)){
  flag <- is.na(Xs[[i]]$marks["dbh"])
  Xs[[i]]$marks["dbh"][flag] <- mean(Xs[[i]]$marks[,2], na.rm = T)
}

BD <- BDpatterns(Xs)
# number of recruits (births) in each census by species
sapply(BD$B, function(o){ table(o$marks$sp) })
# number of mortalities (deaths) in each census by species
sapply(BD$D, function(o){ table(o$marks$sp) })
# number of existing trees in each census by species
sapply(BD$E, function(o){ table(o$marks$sp) })
# missing dbh of recruits in each census by species
sapply(BD$B, function(o){ by(o$marks$dbh, o$marks$sp,
                           function(v){ mean(is.na(v)) }) })
# missing dbh of mortalities in each census by species
sapply(BD$D, function(o){ by(o$marks$dbh, o$marks$sp,
                             function(v){ mean(is.na(v)) }) })
# missing dbh of existing trees in each census by species
sapply(BD$E, function(o){ by(o$marks$dbh, o$marks$sp,
                             function(v){ mean(is.na(v)) }) })

sapply(BD$B, function(o){ table(o$marks$sp) })
sapply(BD$D, function(o){ table(o$marks$sp) })

sapply(BD$E, function(o){ table(o$marks$sp) })+
sapply(BD$B, function(o){ table(o$marks$sp) })

#================================================================
#sequence of truncation distances
#================================================================
bw.seq <- seq(5,155, l = 16)
length(BD$B)
length(BD$D)
length(BD$E)


#====================================================================
#Filtering data
#====================================================================
Time.steps <- 7
Z.deaths <- Z.births <- lapply(1:Time.steps, function(x) bci.covars[-c(4,5,14,15)] )

table.list <- lapply(Xs.c, function(x) data.frame(sort( table( marks(x)[,1] ) )) )

merge.tmp <- table.list[[1]]
for(i in 2:6){
  merge.tmp <- merge(merge.tmp, table.list[[i]], by = "Var1", all = TRUE, suffixes = c(paste(i-1,"_var"), paste(i,"_var")) )
}
matplot( t(merge.tmp[,-1]), type = "l")

na.flag <- is.na(rowMeans(merge.tmp[,-1]))
filter1 <- merge.tmp[!na.flag,]
flag.size <- (rowMeans(filter1[,-1]) <= 500)
filter2 <- filter1[!flag.size,]
species.names <- as.character(filter2[,1])

stats.mean <- matrix(0,length(species.names), ncol = Time.steps)
stats.sd <- matrix(0,length(species.names), ncol = Time.steps)
for(i in 1:Time.steps){
  for(j in 1:length(species.names)){
    X.temp <- Xs.c[[i]][marks(Xs.c[[i]])[,1] == species.names[j]]
    dbh <- marks(X.temp)[,2]
    stats.mean[j,i] <- mean(dbh, na.rm = TRUE)
    stats.sd[j,i] <- sd(dbh, na.rm = TRUE)
  }
}

apply(stats.mean, 1, quantile)
filter3 <- rowSums( stats.mean >= 1 ) != 0

species.names.final <- as.character(filter2[filter3,1])

#====================================================================
#Imputing missing values
#====================================================================

# For the Xs
for(i in 1:Time.steps){
  flag <- is.na(Xs[[i]]$marks["dbh"])
  Xs[[i]]$marks["dbh"][flag] <- mean(Xs[[i]]$marks[,2], na.rm = T)
}



# For the Xs.c
for(i in 1:Time.steps){
  result <- aggregate(dbh ~ sp, data = Xs.c[[i]]$marks, FUN = mean)
  for(j in 1:nrow(result)){
    flag.na <- is.na(Xs.c[[i]]$marks["dbh"])
    flag.sp <- Xs.c[[i]]$marks["sp"] == result[j,1]
    Xs.c[[i]]$marks[flag.sp & flag.na,] <- result[j,2]
  }
  cat(i, "\n")
}

#====================================================================
#Stochastic covariate influence of capparis species on capparis (births only)
#====================================================================

Kappa.list.0 <- list()
for(i in 1:Time.steps){
  Kappa_res <- Kappa_fun_marks_min(Xs[[i]], psi = 0.25,
                               ox.seq = seq(0, 1000, l = 201),
                               oy.seq = seq(0,  500, l = 101))
  Kappa.list.0[[i]] <- Kappa_res
  cat(i,"\n")
}


#====================================================================
#Influence of other species on capparis recruits and capparis deaths
#====================================================================
Kappa.list.2 <- list()
for(i in 1:Time.steps){
  Kappa_res <- as.im(matrix(0, ncol = 200, nrow = 100, byrow = T),
                     xrange = c(0,1000),
                     yrange = c(0, 500))
  for(j in 1:length(species.names.final)){
    X.temp <- Xs.c[[i]][marks(Xs.c[[i]])[,1] == species.names.final[j]]
    na.flag <- is.na(marks(Xs.c[[i]])[,2])
    marks(Xs.c[[i]])[na.flag,2] <- mean(marks(Xs.c[[i]])[!na.flag,2])
    psi <- 5

    Kappa.fun.tmp <- Kappa_fun_marks(X.temp, psi = psi,
                             ox.seq = seq(0, 1000, l = 201),
                             oy.seq = seq(0,  500, l = 101))

    #dbh <- marks(X.temp)[,2]
    Kappa_res <- Kappa_res + Kappa.fun.tmp
    cat(j, "\n")
  }
  Kappa.list.2[[i]] <- Kappa_res/length(species.names.final)
  #plot(Kappa.list.2[[i]], main = i)
  cat(i, "\n")
}

saveRDS(Kappa.list.0, "Data/Kappa_births.RDS")
saveRDS(Kappa.list.2, "Data/Kappa_deaths.RDS")


Kappa.list.0 <- readRDS("Data/Kappa_births.RDS")
Kappa.list.2 = readRDS("Data/Kappa_deaths.RDS")

#stacking the covariates in a list
for(i in 1:Time.steps){Z.births[[i]]$Kappa.nc <- Kappa.list.2[[i]]; Z.births[[i]]$Kappa.caparis <- Kappa.list.0[[i]]}
for(i in 1:Time.steps){Z.deaths[[i]]$Kappa.nc <- Kappa.list.2[[i]];}


#====================================================================
#Stochastic covariate influence of capparis species on capparis (deaths only)
#====================================================================

knn.val <- function(X){
  knn.vals <- pairdist(X)
  www <- apply( knn.vals, 1, function(x) sum( exp(-(x[x != 0]/10)^2 )) )
  return(www)
}

Death.all <- NULL
Surv.all <- NULL
for(kk in 1:Time.steps){
  X.tmp <- superimpose(BD$D[[kk]], BD$E[[kk]])
  fun.val <- knn.val(X.tmp)
  knn.vals <- pairdist(X.tmp)
  www <- apply( knn.vals, 1, function(x) sum( exp(-(x[x != 0]/10)^2 )*X.tmp$marks[which(x != 0),2]) )
  Death.all <- c(Death.all, www[1:BD$D[[kk]]$n])
  Surv.all <- c(Surv.all, www[(BD$D[[kk]]$n + 1):X.tmp$n])
}

Z.deaths.vec <- c(Surv.all, Death.all) #Vector of variables
#====================================================================
#estimation of parameters and covariance matrices
#====================================================================

Now.deaths <-Sys.time()
param.deaths <- fitting_deaths(Death = BD$D,
                               Surv = BD$E,
                               Z.deaths = Z.deaths,
                               Z.vector = c(Death.all, Surv.all),
                               Time.steps = Time.steps,
                               year = TRUE)
kernel <- function(h,bw){
  return(as.numeric(h<bw))#simple uniform kernel used here.
}

Cova.estim.deaths <- Covariance.estimate.deaths(fit.D = param.deaths$coefficients,
                                                Survs = BD$E,
                                                Deaths = BD$D,
                                                Z.deaths = Z.deaths,
                                                Z.vector = Z.deaths.vec,
                                                bw.seq = bw.seq,
                                                Time.steps = Time.steps,
                                                year = TRUE)
Later.deaths <-Sys.time()

Now.births <-Sys.time()
param.births <- fitting_births(BD$B, Z.births, Time.steps, year = TRUE)

Cova.estim.births <- Covariance.estimate.births(fit.B = param.births$coef,
                                                Births = BD$B,
                                                Births.dummy = param.births$dummy,
                                                Z.births = Z.births,
                                                rho0 = param.births$rho0,
                                                bw.seq = bw.seq,
                                                Time.steps = Time.steps,
                                                year = TRUE)
Later.births <-Sys.time()

round(cbind(param.deaths$coefficients,param.births$coef), 4)

Later <- Sys.time()


#==============================================================
#p-values and confidence intervals
#==============================================================
flag.D <- NULL
flag.R <- NULL
pval.D = pval.R = NULL
for(i in 1:16){
  Low.lim.Deaths <- param.deaths$coefficients - sqrt(diag(Cova.estim.deaths[[i]]))*qnorm(0.975)
  Upp.lim.Deaths <- param.deaths$coefficients + sqrt(diag(Cova.estim.deaths[[i]]))*qnorm(0.975)

  Low.lim.Recruits <- param.births$coef - sqrt(diag(Cova.estim.births[[i]]))*qnorm(0.975)
  Upp.lim.Recruits <- param.births$coef + sqrt(diag(Cova.estim.births[[i]]))*qnorm(0.975)

  z.Deaths=param.deaths$coefficients/sqrt(diag(Cova.estim.deaths[[i]]))
  z.Recruits=param.births$coef/sqrt(diag(Cova.estim.births[[i]]))

  pval.D=cbind(pval.D,2*(1-pnorm(abs(z.Deaths))))
  pval.R=cbind(pval.R,2*(1-pnorm(abs(z.Recruits))))

  flag.D <- cbind(flag.D, (Low.lim.Deaths <= 0 & 0 <= Upp.lim.Deaths) )
  flag.R <- cbind(flag.R, (Low.lim.Recruits <= 0 & 0 <= Upp.lim.Recruits) )
}


#==========================================================================
#INLA.fit
#==========================================================================
library("fmesher")
library("INLA")
library("pbmcapply")


B = BD$B
D = BD$D
S = BD$E

Kappa0 <- Kappa1 <- list()
for(k0 in 1:7){
  Kappa0[[k0]] <- Z.births[[k0]][12]$Kappa.nc
  Kappa1[[k0]] <- Z.births[[k0]][13]$Kappa.caparis
}

covars <- Z.births[[1]][1:11]

st.covars = list( stB = cbind(Kappa.nc = Kappa0, Kappa.caparis = Kappa1),
                  stD = cbind(Kappa.nc = Kappa0, Kappa.caparis = Kappa1))
num.threads <- 50

inlafit_data <- function(B, # list of point patterns of births
                       D, # list of point patterns of deaths
                       S, # list of point patterns of survivors
                       covars, # list of spatial covariates
                       st.covars, # list of spatiotemporal influence covariates
                       sp="cappfr", # selected species
                       dimyx=c(100, 200), # dimension of discretization for the birth model
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
    "as.integer(status) ~ time + f(dmesh$idx$loc, model=dspde,
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
    "counts ~ time + f(bmesh$idx$loc, model=bspde,
                   group=k, control.group=list(model='ar1')) + ",
    paste(names(covars), collapse = "+"),
    "+ stB1 + stB2"))
  bfit <- inla(bform, data=bdat, family="poisson", num.threads=num.threads)
  b <- bfit$summary.fixed[, c(1:3, 5)]

  # estimated spatiotemporal random effect
  bmpj <- inla.mesh.projector(bmesh,  dims=rev(dimyx))
  bm <- list()
  for(k in 1:K)
  {
    bm[[k]] <- inla.mesh.project(bmpj,
                                 bfit$summary.random$`bmesh|S|idx|S|loc`$mean[(k - 1) * bmesh$n + 1:bmesh$n])
    bm[[k]] <- as.im(list(x=bmpj$x, y=bmpj$y, z=bm[[k]]), W=W)
  }
  later.births <- Sys.time()

  return(list(beta=cbind(births = b, deaths = d), dfit=dfit, bfit=bfit,
              btheta=inla.spde2.result(bfit, "bmesh|S|idx|S|loc", bspde),
              dtheta=inla.spde2.result(dfit, "dmesh|S|idx|S|loc", dspde),
              dm=dm, bm=bm, time.births = later.births - now.births,
              time.deaths = later.deaths - now.deaths))
}

inla.fit.beta.0res <- cbind(births = b, deaths = d)

inla.fit.beta <- inlafit_data(B, # list of point patterns of births
                         D, # list of point patterns of deaths
                         S, # list of point patterns of survivors
                         covars, # list of spatial covariates
                         st.covars, # list of spatiotemporal influence covariates
                         sp="cappfr", # selected species
                         dimyx=c(100, 200), # dimension of discretization for the birth model
                         num.threads=50)



round(INLA_data_fit_results$inla.fit.beta.0res["births.mean"], 2)
round(INLA_data_fit_results$inla.fit.beta$beta["births.mean"], 2)
round(INLA_data_fit_results$inla.fit.beta$beta["deaths.mean"], 2)


INLA_data_fit_results$inla.fit.beta$time.deaths
INLA_data_fit_results$inla.fit.beta$time.births


save.image("BCI_data_results.RDAT")

