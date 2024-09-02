set.seed(123456)
source("CE_Method.R")
M <- 101
N <- 201
N1 <- max(N,M)
par(mfrow = c(1,2))
param.B <- list(mean = 0, smooth = 1.75, scale = 4, sill  = 1, nugget = 0)
param.D <- list(mean = 0, smooth = 0.5, scale = 7, sill  = 1, nugget = 0)

param.ini <- param.B
param.ini$scale <- 4*param.ini$scale

Covariate.sim.B <- GRF.CE.Sim(nsim = 1, xrange = W.sim$xrange, yrange = W.sim$xrange, M = N1, N = N1, 
                              corrmodel = "Matern", param = param.ini, mean.val = 0, max.ext = 4) 

Covariate.B <- im(yrow = unique(Covariate.sim.B$grid.points[,1]), 
                  xcol = unique(Covariate.sim.B$grid.points[,2]), 
                  matrix(Covariate.sim.B$X[,1]/3 , ncol = N1, nrow = N1) )


expand.grid <- expand.grid(seq(W.sim$xrange[1], W.sim$xrange[2], l = N), seq(W.sim$yrange[1], W.sim$yrange[2], l = M))
Covariate.B.vec <- Covariate.B[expand.grid]
Covariate.B <- im( matrix(Covariate.B.vec , ncol = N, nrow = M, byrow = T), xrange = W.sim$xrange, yrange = W.sim$yrange)

param.ini <- param.D
param.ini$scale <- 4*param.ini$scale
Covariate.sim.D <- GRF.CE.Sim(nsim = 1, W.sim$xrange, W.sim$xrange, M = N1, N = N1, 
                              corrmodel = "Matern", param = param.ini, mean.val = 0, max.ext = 4) 
Covariate.D <- im(yrow = unique(Covariate.sim.B$grid.points[,1]), 
                  xcol = unique(Covariate.sim.B$grid.points[,2]),
                  matrix(Covariate.sim.D$X[,1]/3 , ncol = N1, nrow = N1))

Covariate.D.vec <- Covariate.D[expand.grid]
Covariate.D <- im( matrix(Covariate.D.vec , ncol = N, nrow = M, byrow = T), xrange = W.sim$xrange, yrange = W.sim$yrange)
plot(Covariate.D)

Covariate.B <- Covariate.B
Covariate.D <- Covariate.D

Z.deaths <- Z.births <- list()
Z.deaths <- lapply(1:Time.steps, function(x) Covariate.D ) #Using Cu as covariate
Z.births <- lapply(1:Time.steps, function(x) Covariate.B ) #Using Cu as covariate

rm(Covariate.sim.B); rm(Covariate.D)
rm(Covariate.sim.D); rm(Covariate.B)

x <- seq(W.sim$xrange[1], W.sim$xrange[2], l = N )
y <- seq(W.sim$yrange[1], W.sim$yrange[2], l = M )


#plot(Covariate.B, main = "")
#lines(W.sim)#
#plot(Covariate.B, main = "")

