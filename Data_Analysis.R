rm(list = ls())
library("spatstat")
library("glmnet")

# ==============================================================
# read BCI data
# ==============================================================
# extract selected species point patterns from the BCI data
BCIextract <- function(species){
  # path to data directory
  path <- '/data/fcp-X79-03/Spatio_temporal_PP/Data/'
  # observation window
  W <- owin(xrange=c(0, 1000), yrange=c(0, 500), 
            unitname=c("meter", "meters"))
  # output list of point patterns
  X <- vector('list', length=8)
  for (i in 1:8)
  {
    tmp <- get(load(paste0(path, "bci.tree", i, ".rdata", sep="")))
    ok <- (!is.na(tmp$gx)) & (!is.na(tmp$gy)) 
    ok <- ok & (tmp$status == "A") & (tmp$sp %in% species)
    xx <- tmp$gx[ok]
    yy <- tmp$gy[ok]
    mm <- data.frame(sp=tmp$sp[ok], dbh=tmp$dbh[ok])
    X[[i]] <- ppp(xx, yy, window=W, marks=mm)
  }
  names(X) <- c(1983, 1985, 1990, 1995, 2000, 2005, 2010, 2015)
  return(X)
}

BCIextract_complement <- function(species){
  # path to data directory
  path <- '/data/fcp-X79-03/Spatio_temporal_PP/Data/'
  # observation window
  W <- owin(xrange=c(0, 1000), yrange=c(0, 500), 
            unitname=c("meter", "meters"))
  # output list of point patterns
  X <- vector('list', length=8)
  for (i in 1:8)
  {
    tmp <- get(load(paste0(path, "bci.tree", i, ".rdata", sep="")))
    ok <- (!is.na(tmp$gx)) & (!is.na(tmp$gy)) 
    ok <- ok & (tmp$status == "A") & !(tmp$sp %in% species)
    xx <- tmp$gx[ok]
    yy <- tmp$gy[ok]
    mm <- data.frame(sp=tmp$sp[ok], dbh=tmp$dbh[ok])
    X[[i]] <- ppp(xx, yy, window=W, marks=mm)
  }
  names(X) <- c(1983, 1985, 1990, 1995, 2000, 2005, 2010, 2015)
  return(X)
}

load('/data/fcp-X79-03/Spatio_temporal_PP/Data/bci.spptable.rdata')
bci.spptable[bci.spptable$sp == "hybapr",]$Latin
# selected species
Xs <- BCIextract(c("cappfr"))
Xs.c <- BCIextract_complement(c("cappfr"))
table(Xs.c[[1]]$marks)

# number of alive trees in each census by species
sapply(Xs, function(o){ table(o$marks$sp) })
# missing dbh in each census by species
sapply(Xs, function(o){ by(o$marks$dbh, o$marks$sp, 
                           function(v){ mean(is.na(v)) }) })
# median of dbh in each census by species
sapply(Xs, function(o){ by(o$marks$dbh, o$marks$sp, median, na.rm=TRUE) })

library("ggplot2")
library("ggpubr")
library("ggsci")
# plot of number of alive trees  in each census by species
# g1 <- ggplot(data=reshape2::melt(sapply(Xs, function(o){ table(o$marks$sp) }),
#                            varnames=c("species", "census")),
#        mapping=aes(x=census, y=value,  col=species)) + 
#   geom_line() + geom_point() + 
#   labs(y="number of alive trees") + scale_color_jama()
# # plot of missing dbh  in each census by species
# g2 <- ggplot(data=reshape2::melt(sapply(Xs, function(o){ by(o$marks$dbh, o$marks$sp, 
#                                                       function(v){ mean(is.na(v)) }) }), 
#                            varnames=c("species", "census")),
#        mapping=aes(x=census, y=value, col=species)) + 
#   geom_line() + geom_point() + 
#   labs(y="proportion of missing dbh") + scale_color_jama()
# # box plot of the base 10 logarithm of dbh in each census by species
# g3 <- ggplot(data=rbind(data.frame(Xs[[1]]$marks, census=1983),
#                   data.frame(Xs[[2]]$marks, census=1985),
#                   data.frame(Xs[[3]]$marks, census=1990),
#                   data.frame(Xs[[4]]$marks, census=1995),
#                   data.frame(Xs[[5]]$marks, census=2000),
#                   data.frame(Xs[[6]]$marks, census=2005),
#                   data.frame(Xs[[7]]$marks, census=2010),
#                   data.frame(Xs[[8]]$marks, census=2015)),
#        mapping=aes(x=factor(census), y=log10(dbh), fill=sp)) + 
#   geom_boxplot() + labs(x="census") + 
#   guides(fill=guide_legend(title="species")) + scale_fill_jama() 
# 
# ggarrange(g1 + theme(axis.text.x = element_blank(),
#                      axis.ticks.x = element_blank(),
#                      axis.title.x = element_blank(),
#                      plot.margin = margin(b=1)), 
#           g2 + theme(plot.margin = margin(t=1)), 
#           g3, nrow=3, ncol=1)
#dev.copy2pdf(file="~/Desktop/bci.pdf", width=8, height=6)

# ==============================================================
# read BCI covariate data
# ==============================================================
# BCI covariates
load("/data/fcp-X79-03/Spatio_temporal_PP/Data/bci.covars.Rda")
# adding scaled x and y coordinates as covariates
bci.covars$x <- as.im(function(x, y){ (x - 500) / 1000 }, 
                      W=Xs[[1]]$window, dimyx=c(101, 201))
bci.covars$y <- as.im(function(x, y){ (y - 250) / 500 }, 
                      W=Xs[[1]]$window, dimyx=c(101, 201))
par(mfrow=c(4, 4), mar=c(0.1, 0.1, 1.2, 2))
cvnames <- names(bci.covars)
for (i in 1:length(cvnames)) plot(bci.covars[[i]], main=cvnames[i])
#dev.copy2eps(file="~/Desktop/bcicovars.pdf", width=8, height=4)

# ==============================================================
# point patterns of recruits (births) and mortalities (deaths)
# ==============================================================
Time.steps <- 7
for(i in 1:(Time.steps+1)){
  flag <- is.na(Xs[[i]]$marks["dbh"])
  Xs[[i]]$marks["dbh"][flag] <- mean(Xs[[i]]$marks[,2], na.rm = T)
}


BDpatterns <- function(Xs){
  # number of time instances
  K <- length(Xs)
  D <- B <- E <- vector("list", length=K - 1)
  for (k in 1:(K - 1))
  {
    dd <- crossdist.default(Xs[[k]]$x, Xs[[k]]$y,
                            Xs[[k + 1]]$x, Xs[[k + 1]]$y, 
                            squared=TRUE)
    dd <- (dd == 0)
    print(table(rowSums(dd)))
    death <- (rowSums(dd) == 0)
    D[[k]] <- Xs[[k]][which(death)]
    print(table(colSums(dd)))
    birth <- (colSums(dd) == 0)
    B[[k]] <- Xs[[k + 1]][which(birth)]
    E[[k]] <- Xs[[k + 1]][which(!birth)]
  }
  names(E) <- c("1985", "1990", "1995", 
                "2000", "2005", "2010", "2015")
  names(B) <- names(D) <- names(E)
  return(list(B=B, D=D, E=E))
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

# plotting births, deaths and existing trees in each census
# plotfun <- function(k){
#   gdat <- rbind(data.frame(x=BD$B[[k]]$x, y=BD$B[[k]]$y, 
#                            BD$B[[k]]$marks, type="birth"),
#                 data.frame(x=BD$D[[k]]$x, y=BD$D[[k]]$y, 
#                            BD$D[[k]]$marks, type="death"),
#                 data.frame(x=BD$E[[k]]$x, y=BD$E[[k]]$y, 
#                            BD$E[[k]]$marks, type="exist"))
#   ggplot(data=gdat, mapping=aes(x=x, y=y, col=type)) + 
#     geom_rect(aes(xmin=0, xmax=1000, ymin=0, ymax=500), 
#               fill=NA, color="grey50", linetype=1, size=0.2) + 
#     geom_point(aes(size=dbh, shape=sp)) + 
#     scale_shape_manual(values=c(0, 1, 2)) + 
#     scale_color_manual(values=alpha(c("blue", "red", "green"), c(1, 1, 0.4))) +
#     scale_x_continuous(limits = c(0, 1000), expand = c(0, 0)) +
#     scale_y_continuous(limits = c(0, 500), expand = c(0, 0)) + 
#     ggtitle(paste("Census", names(BD$B)[k])) + 
#     guides(shape=guide_legend(title="species")) +
#     coord_fixed() + theme_void()
# }
# 
# g1 <- plotfun(1)
# g2 <- plotfun(2)
# g3 <- plotfun(3)
# g4 <- plotfun(4)
# g5 <- plotfun(5)
# g6 <- plotfun(6)
# g7 <- plotfun(7)
# 
# ggarrange(g1, g2, g3, g4, g5, g6, g7, 
#           nrow=4, ncol=2, common.legend = TRUE, legend="right")
#ggsave("~/Desktop/BDE.pdf", width = 8, height = 8.2)
setwd("/data/fcp-X79-03/Spatio_temporal_PP/stBCIcode/SimulationStudy/")

source("../asympvar_code.R")
source("../asympvar_code_bw.R")
source("Kappa_fun.R")

bw.seq <- seq(5,155, l = 16)
length(BD$B)
length(BD$D)
length(BD$E)

countings <- cbind(
  unlist( lapply(BD$B, function(x) x$n) ),
  unlist( lapply(BD$D, function(x) x$n) ),
  unlist( lapply(BD$E, function(x) x$n) )
)

countings[1:7,3] + countings[1:7,1]

countings[1:6,3] + countings[2:7,1] -  countings[2:7,2]

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

####################
#### Filling NA ####
####################

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

Kappa.list.0 <- list()
for(i in 1:Time.steps){
  Kappa_res <- Kappa_fun_marks_min(Xs[[i]], psi = 0.25,
                               ox.seq = seq(0, 1000, l = 201), 
                               oy.seq = seq(0,  500, l = 101))
  Kappa.list.0[[i]] <- Kappa_res
  cat(i,"\n")
}

Kappa.list.1 <- list()
for(i in 1:Time.steps){
  Kappa_res <- Kappa_fun_marks(Xs[[i]], psi = 2.5,
                    ox.seq = seq(0, 1000, l = 201), 
                    oy.seq = seq(0,  500, l = 101))
  Kappa.list.1[[i]] <- Kappa_res
  cat(i,"\n")
}

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

for(i in 1:Time.steps){Z.births[[i]]$Kappa.nc <- Kappa.list.2[[i]]; Z.births[[i]]$Kappa.caparis <- Kappa.list.0[[i]]}
for(i in 1:Time.steps){Z.deaths[[i]]$Kappa.nc <- Kappa.list.2[[i]];}

#saveRDS(Kappa.list.0, "Kappa_births.RDS")
#saveRDS(Kappa.list.2, "Kappa_deaths.RDS")

Kappa.list.0 <- readRDS("Kappa_births.RDS")


knn.val <- function(X){
  knn.vals <- pairdist(X)
  www <- apply( knn.vals, 1, function(x) sum( exp(-(x[x != 0]/10)^2 )) )
  return(www)
}

Death.all <- NULL
Surv.all <- NULL
for(kk in 1:Time.steps){
  #Z.all[[kk]] <- c(Z.deaths[kk],Z.births[kk])
  X.tmp <- superimpose(BD$D[[kk]], BD$E[[kk]])
  fun.val <- knn.val(X.tmp)
  knn.vals <- pairdist(X.tmp)
  www <- apply( knn.vals, 1, function(x) sum( exp(-(x[x != 0]/10)^2 )*X.tmp$marks[which(x != 0),2]) )
  Death.all <- c(Death.all, www[1:BD$D[[kk]]$n])
  Surv.all <- c(Surv.all, www[(BD$D[[kk]]$n + 1):X.tmp$n])
}

fitting_deaths <- function(Death, Surv, Z.deaths, Z.vector, Time.steps, year = FALSE){
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
  if(is.null(Z.vector)){
    Cova.D <- rbind(Cova.Deaths, Cova.Survivals)
  }else{
    Cova.D <- rbind(Cova.Deaths, Cova.Survivals)
    Cova.D <- cbind(Cova.D[,-ncol(Cova.D)], Z.vector, Cova.D[,ncol(Cova.D)])
  }

  if(year){
    Year.factor <- factor(Cova.D[,ncol(Cova.D)])
    Cova.D <- Cova.D[,-ncol(Cova.D)]
    fit.D <- glm(values.D  ~ Cova.D + Year.factor + 0, family=binomial)
  }else{
    fit.D <- glm(values.D  ~ Cova.D, family=binomial)
  }
  return(fit.D)
}

fitting_births <- function(Births, Z.births, Time.steps, year = FALSE){
  # Fitting the Births via logistic regression
  rho0 <- 4*max( do.call("c", lapply(Births, function(x) x$n)) )/spatstat.geom::area(Births[[1]]$window)
  Dummy.sim <- lapply( as.list(1:Time.steps), function(i){
    W0 <- Births[[i]]$window
    X0 <- rpoispp(rho0, win = W0)
    return(X0)
  })
  eval.covariate <- function(i, births, year){
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
    if(year) eval.list <- cbind(eval.list, i)
    return(eval.list)
  }
  Births.val <- lapply(as.list(1:Time.steps), eval.covariate, births = TRUE, year)
  Births.val <- do.call("rbind", Births.val)
  Dummy.val <- lapply(as.list(1:Time.steps), eval.covariate, births = FALSE, year)
  Dummy.val <- do.call("rbind", Dummy.val)
  Cova.B <- rbind(Births.val, Dummy.val)
  Flag.vec.B <- rep(c(1,0), c(nrow(Births.val), nrow(Dummy.val)))
  values.B <- factor( Flag.vec.B )
  if(year){
    Year.factor <- factor(Cova.B[,ncol(Cova.B)])
    Cova.B <- Cova.B[,-ncol(Cova.B)]
    fit.B <- glm(values.B  ~ Cova.B + Year.factor + 0 + offset(rep( -log(rho0), length(Flag.vec.B) )), family=binomial)
  }else{
    fit.B <- glm(values.B  ~ Cova.B + offset(rep( -log(rho0), length(Flag.vec.B) )), family=binomial)
  }
  return( list(coef = fit.B$coefficients, dummy = Dummy.sim, rho0 = rho0, fitting.obj = fit.B, flag = values.B) )
}

#BD$D <- lapply(BD$D, function(x) unmark(x))
#BD$E <- lapply(BD$E, function(x) unmark(x))
#BD$B <- lapply(BD$B, function(x) unmark(x))


Covariance.estimate.deaths <- function(fit.D, Survs, Deaths, 
                                       Z.deaths, Z.vector, 
                                       bw.seq, Time.steps, year = FALSE){
  #browser()
  eval.covariate <- function(i, deaths, year = FALSE){
    n0 <- length(Z.deaths[[i]])
    if(deaths){
      eval.list <- matrix(0, ncol = length(fit.D), nrow = Deaths[[i]]$n)
    }else{
      eval.list <- matrix(0, ncol = length(fit.D), nrow = Survs[[i]]$n)
    }
    for(jj in 1:n0){
      if(deaths){
        eval.list[,jj] <- Z.deaths[[i]][[jj]][ Deaths[[i]] ]
      }else{
        eval.list[,jj] <- Z.deaths[[i]][[jj]][ Survs[[i]] ]
      }
    }
    
    if(!is.null(Z.vector)){
      X.tmp <- superimpose(Deaths[[i]], Survs[[i]])
      knn.vals <- pairdist(X.tmp)
      www <- apply( knn.vals, 1, function(x) sum( exp(-(x[x != 0]/10)^2 )*X.tmp$marks[which(x != 0),2]) )
      Death.all <- www[1:BD$D[[kk]]$n]
      Surv.all <- www[(BD$D[[kk]]$n + 1):X.tmp$n]
      if(deaths){
        eval.list[,(n0+1)] <- www[1:Deaths[[i]]$n]
      }else{
        eval.list[,(n0+1)] <- www[(Deaths[[i]]$n + 1):X.tmp$n]
      }
      
    }
    if(year) eval.list[,(n0 + i + !is.null(Z.vector))] <- 1
    return(eval.list)
  }
  
  Sens.Var <- lapply(as.list(1:Time.steps), function(i){
    ### Deaths
    sensiv.D <- list()
    var.D <- list()
    Xkminus1 <- superimpose(Survs[[i]], Deaths[[i]], check = F)
    Deaths.cova.vec <- eval.covariate(i, deaths = TRUE, year)
    Survivals.cova.vec <- eval.covariate(i, deaths = FALSE, year)
    cova.D <- rbind(Survivals.cova.vec, Deaths.cova.vec)
    #Flag value
    Flag.D <- rep( c(0,1), c(Survs[[i]]$n, Deaths[[i]]$n))
    grad.etha.k <- cbind(cova.D)
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

Covariance.estimate.births <- function(fit.B, Births, Births.dummy, Z.births, rho0, bw.seq, Time.steps, year = FALSE){
  eval.covariate <- function(i, births, year){
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
    if(year) eval.list <- cbind(eval.list, i)
    return(eval.list)
  }
  Sens.Var <- lapply(as.list(1:Time.steps), function(i){
    Births.cova <- eval.covariate(i, births = TRUE, year)
    Dummy.cova <- eval.covariate(i, births = FALSE, year)
    if(!year){
      gradlogzeta.k  <- rbind(cbind(1,Births.cova), cbind(1,Dummy.cova))
      zeta.k <- exp( c(gradlogzeta.k%*%fit.B) )
    }else{
      zero.mat.births <- matrix(0, nrow = nrow(Births.cova), ncol = Time.steps )
      zero.mat.dummy <- matrix(0, nrow = nrow(Dummy.cova), ncol = Time.steps )
      zero.mat.births[,i] <- 1
      zero.mat.dummy[,i] <- 1
      gradlogzeta.k <- rbind(cbind(Births.cova[,-ncol(Births.cova)], zero.mat.births), cbind(Dummy.cova[,-ncol(Dummy.cova)], zero.mat.dummy))
      zeta.k <- exp( c(gradlogzeta.k%*%fit.B) )
    }
    is.datapoint <- c(rep(TRUE, Births[[i]]$n), rep(FALSE, Births.dummy[[i]]$n) )
    sensiv.B.jy <- sensi.birth.k.jy(zeta.k, gradlogzeta.k, is.datapoint,rho0)
    #www <- pcfinhom(Bk, lambda = zeta.k[is.datapoint])
    #h0 <- www$r[ which.min( abs(www$iso-1) ) ]
    Births.un <- unmark(Births[[i]])
    Births.dummy.un <- unmark(Births.dummy[[i]])
    var.B.jy <- simplify2array(e.birth.var.k.jy.bw(zeta.k = zeta.k,
                                                   gradlogzeta.k = gradlogzeta.k,
                                                   rho0 = rho0,
                                                   is.datapoint = is.datapoint,
                                                   Bk = Births.un,
                                                   Y = Births.dummy.un,
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
                                                Z.vector = c(Surv.all, Death.all),
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

round(cbind(param.deaths$coefficients,
            param.births$coef), 4)

Later <- Sys.time()

flag.D <- NULL
flag.R <- NULL
for(i in 1:16){
  Low.lim.Deaths <- param.deaths$coefficients - sqrt(diag(Cova.estim.deaths[[i]]))*qnorm(0.975)
  Upp.lim.Deaths <- param.deaths$coefficients + sqrt(diag(Cova.estim.deaths[[i]]))*qnorm(0.975)

  Low.lim.Recruits <- param.births$coef - sqrt(diag(Cova.estim.births[[i]]))*qnorm(0.975)
  Upp.lim.Recruits <- param.births$coef + sqrt(diag(Cova.estim.births[[i]]))*qnorm(0.975)

  flag.D <- cbind(flag.D, (Low.lim.Deaths <= 0 & 0 <= Upp.lim.Deaths) )
  flag.R <- cbind(flag.R, (Low.lim.Recruits <= 0 & 0 <= Upp.lim.Recruits) )
}

#matplot(bw.seq, cbind(colSums(flag.D),
#                      colSums(flag.R)) )

#matplot(bw.seq, t(flag.D), type = "l")
#matplot(bw.seq, t(flag.R), type = "l")

#vals <- cumsum( do.call("c", lapply(BD$B, function(x) x$n) ) )
#Limits <- cbind( c(1,(vals[1:6]+1)), vals[1:7])
#Lambda <- exp( predict(param.births$fitting.obj)+log(param.births$rho0))  

#pfc.inhom.list <- lapply(as.list(1:Time.steps), function(i) cbind( pcfinhom(BD$B[[i]], r = seq(0,200), lambda = Lambda[Limits[i,1]:Limits[i,2]])$iso ))
#pcf.mat <- do.call("cbind", pfc.inhom.list)

#pdf("inhom_pcf.pdf")
#matplot(seq(0,200), pcf.mat, type = "l", lty = 3, lwd = 2, ylab = expression( hat(g)(r) ), xlab = "r", xlim = c(3, 100) )
#legend("topright", legend = names(BD$B), lwd = 2, lty = 3, col = 1:length(names(BD$B)) )
#lines(seq(0,200),rowMeans(pcf.mat), col = "green", lwd = 3, lty = 1)
#abline(h = 1, lwd = 2)
#dev.off()


#pfc.inhom.list <- lapply(as.list(1:Time.steps), function(i) cbind( pcfinhom(BD$D[[i]], r = seq(0,200), lambda = Lambda.D[[i]])$iso ))
#pcf.mat <- do.call("cbind", pfc.inhom.list)
#matplot( pcf.mat, type = "l", lty = 1, lwd = 2, ylab = expression( hat(g)(r) ), xlab = "r", xlim = c(5, 100) )
#legend("topright", legend = names(BD$D), lwd = 2, lty = 1, col = 1:length(names(BD$D)) )

#i <- 6
#Low.lim.Deaths <- param.deaths$coefficients - sqrt(diag(Cova.estim.deaths[[i]]))*qnorm(0.975)
#Upp.lim.Deaths <- param.deaths$coefficients + sqrt(diag(Cova.estim.deaths[[i]]))*qnorm(0.975)

#Low.lim.Recruits <- param.births$coef - sqrt(diag(Cova.estim.births[[i]]))*qnorm(0.975)
#Upp.lim.Recruits <- param.births$coef + sqrt(diag(Cova.estim.births[[i]]))*qnorm(0.975)

#Z_score.B <- param.births$coef/sqrt(diag(Cova.estim.births[[i]]))
#Z_score.D <- param.deaths$coefficients/sqrt(diag(Cova.estim.deaths[[i]]))

#P_vals.B <- 2*(1 - pnorm( abs(Z_score.B) ))  
#P_vals.D <- 2*(1 - pnorm( abs(Z_score.D) )) 

#recruits <- cbind(param.births$coef, Low.lim.Recruits, Upp.lim.Recruits, Z_score.B, P_vals.B)
#rownames(recruits)[1:13] <- c(names(Z.births[[i]]) )
#round(recruits, 4)

#deaths <- cbind(param.deaths$coefficients, Low.lim.Deaths, Upp.lim.Deaths, Z_score.D, P_vals.D)
#rownames(deaths)[1:12] <- c(names(Z.deaths[[i]]))
#round(deaths, 4)[,c(1,4,5)]
 
############# INLA.fit #############
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
