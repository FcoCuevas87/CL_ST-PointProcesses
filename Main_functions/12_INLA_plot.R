rm(list = ls())
inla.result.list <- readRDS("inla_results.RDS")

Param1 <- do.call("rbind",  lapply(inla.result.list, function(x) x[1:5,1] ))
Param2 <- do.call("rbind",  lapply(inla.result.list, function(x) x[6:10,1] ))

N0 <- 500
prefix <- "W2"
experiment <- "Experiment5"
setting <- readRDS(paste(experiment,"/SimSetting_",N0,"_",prefix,"_all.RDS", sep = ""))
Births.fit <- readRDS(file = paste(experiment, "/Births_fit_",N0,"_",prefix,"_",setting$Kappa_fun,"_all_biv.RDS", sep = ""))
Deaths.fit <- readRDS(file = paste(experiment, "/Deaths_fit_",N0,"_",prefix,"_",setting$Kappa_fun,"_all_biv.RDS", sep = ""))

births.fit <- do.call("rbind", lapply(Births.fit, function(x) x$coef)); rm(Births.fit)
deaths.fit <- Deaths.fit

cex.val <- 2
pdf( paste("INLA_beta0_deaths_",prefix,".pdf", sep = "") )
#hist(deaths.fit[,1], col = "blue", border = "white", main = "", xlab = "", ylab = "", lwd = 3, breaks = 9, cex.lab = cex.val, cex.axis = cex.val, xlim = c(-0.45,-0.074)); #

plot( density(deaths.fit[1:100,1]), col = "blue", main = "", xlab = "", ylab = "", lwd = 3, cex.lab = cex.val, cex.axis = cex.val, xlim = c(-0.65,0.10)); #
abline(v = setting$etha.vec[1], col = "red", lwd = 3)

lines( density(Param2[,1]), col = "dark green", lty = 2, lwd = 3)
dev.off()


pdf( paste("INLA_beta1_deaths_",prefix,".pdf", sep = "") )
#hist(deaths.fit[,2], col = "blue", border = "white", main = "", xlab = "", ylab = "", lwd = 3, breaks = 9, cex.lab = cex.val, cex.axis = cex.val, xlim = c(-0.12,0.74)); #, breaks = 5
plot( density(deaths.fit[1:100,2]), col = "blue", main = "", xlab = "", ylab = "", lwd = 3, cex.lab = cex.val, cex.axis = cex.val, xlim = c(-0.30,0.7)); #, breaks = 5
abline(v = 0.25, col = "red", lwd = 3)

lines( density(Param2[,3]), col = "dark green", lty = 2, lwd = 3)
dev.off()


pdf( paste("INLA_beta2_deaths_",prefix,".pdf", sep = "") )
#hist(deaths.fit[,3], col = "blue", border = "white", main = "", xlab = "", ylab = "", lwd = 3, breaks = 9, cex.lab = cex.val, cex.axis = cex.val, xlim = c(-0.33,0.39)); #breaks = 5,
plot( density(deaths.fit[1:100,3]), col = "blue", main = "", xlab = "", ylab = "", lwd = 3, cex.lab = cex.val, cex.axis = cex.val, xlim = c(-0.61,0.50)); #breaks = 5,
abline(v = 0, col = "red", lwd = 3)

lines( density(Param2[,2]), col = "dark green", lty = 2, lwd = 3)
dev.off()


pdf( paste("INLA_beta3_deaths_",prefix,".pdf", sep = "") )
#hist(deaths.fit[,4], col = "blue", border = "white", main = "", xlab = "", ylab = "", lwd = 3, breaks = 9, cex.lab = cex.val, cex.axis = cex.val, xlim = c(-0.12,0.12)); #breaks = 5, 
plot( density(deaths.fit[1:100,4]), col = "blue", main = "", xlab = "", ylab = "", lwd = 3, cex.lab = cex.val, cex.axis = cex.val, xlim = c(-0.5,-0)); #breaks = 5, 
abline(v = -0.25, col = "red", lwd = 3)

lines( density(Param2[,4]), col = "dark green", lty = 2, lwd = 3)
dev.off()
#


pdf( paste("INLA_beta4_deaths_",prefix,".pdf", sep = "") )
#hist(deaths.fit[,4], col = "blue", border = "white", main = "", xlab = "", ylab = "", lwd = 3, breaks = 9, cex.lab = cex.val, cex.axis = cex.val, xlim = c(-0.12,0.12)); #breaks = 5, 
plot( density(deaths.fit[1:100,5]), col = "blue", main = "", xlab = "", ylab = "", lwd = 3, cex.lab = cex.val, cex.axis = cex.val, xlim = c(-0.1,0.1)); #breaks = 5, 
abline(v = 0, col = "red", lwd = 3)

lines( density(Param2[,5]), col = "dark green", lty = 2, lwd = 3)
dev.off()


pdf( paste("INLA_beta0_births_",prefix,".pdf", sep = "") )
#hist(births.fit[,1], col = "blue", border = "white", main = "", xlab = "", ylab = "", lwd = 3, breaks = 9, cex.lab = cex.val, cex.axis = cex.val, xlim = c(-6.4, -6.06)); #, breaks = 5
plot(  density(births.fit[1:100,1]), col = "blue", main = "", xlab = "", ylab = "", lwd = 3, cex.lab = cex.val, cex.axis = cex.val, xlim = c(-7.5, -6.68)); #, breaks = 5
abline(v = setting$beta.vec[1], col = "red", lwd = 3)

lines(  density(Births.coef[1:100,1]), col = "blue", main = "", xlab = "", ylab = "", lwd = 3, lty = 2, cex.lab = cex.val, cex.axis = cex.val, xlim = c(-7.5, -6.68)); #, breaks = 5

W <- owin(xrange = c(-2.5,1000), yrange = c(-2.5, 500))
lines( density(Param1[,1] - log(area(W)/(50*100))), col = "dark green", lty = 2, lwd = 3)
dev.off()


pdf( paste("INLA_beta1_births_",prefix,".pdf", sep = "") )
#hist(births.fit[,2], col = "blue", border = "white", main = "", xlab = "", ylab = "", lwd = 3, breaks = 9, cex.lab = cex.val, cex.axis = cex.val, xlim = c(-0.24, 0.23)); #, breaks = 5
plot( density(births.fit[1:100,2]), col = "blue", main = "", xlab = "", ylab = "", lwd = 3, cex.lab = cex.val, cex.axis = cex.val, xlim = c(-0.37, 0.33), ylim = c(0,7)); #, breaks = 5
abline(v = 0, col = "red", lwd = 3)

lines(  density(Births.coef[1:100,2]), col = "blue", main = "", lty = 2, xlab = "", ylab = "", lwd = 3, cex.lab = cex.val, cex.axis = cex.val, xlim = c(-7.5, -6.68)); #, breaks = 5

lines( density(Param1[,3]), col = "dark green", lty = 2, lwd = 3)
dev.off()


pdf( paste("INLA_beta2_births_",prefix,".pdf", sep = "") )
#hist(births.fit[,3], col = "blue", border = "white", main = "", xlab = "", ylab = "", lwd = 3, breaks = 9, cex.lab = cex.val, cex.axis = cex.val, xlim = c(-0.15,0.38)); #, breaks = 5
plot( density(births.fit[1:100,3]), col = "blue", main = "", xlab = "", ylab = "", lwd = 3, cex.lab = cex.val, cex.axis = cex.val, xlim = c(-0.30,0.47)); #, breaks = 5
abline(v = 0.1, col = "red", lwd = 3)

lines(  density(Births.coef[1:100,3]), col = "blue", main = "", xlab = "", ylab = "", lty = 2, lwd = 3, cex.lab = cex.val, cex.axis = cex.val, xlim = c(-7.5, -6.68)); #, breaks = 5

lines( density(Param1[,2]), col = "dark green", lty = 2, lwd = 3)
dev.off()


pdf( paste("INLA_beta3_births_",prefix,".pdf", sep = "") )
#hist(births.fit[,4], col = "blue", border = "white", main = "", xlab = "", ylab = "", lwd = 3, breaks = 9, cex.lab = cex.val, cex.axis = cex.val, xlim = c(-0.3, 0.49)); #, breaks = 5
plot( density(births.fit[1:100,4]), col = "blue", main = "", xlab = "", ylab = "", lwd = 3, cex.lab = cex.val, cex.axis = cex.val, xlim = c(-0.8, 0.72)); #, breaks = 5
abline(v = 0.1, col = "red", lwd = 3)

lines(  density(Births.coef[1:100,4]), col = "blue", lty = 2, main = "", xlab = "", ylab = "", lwd = 3, cex.lab = cex.val, cex.axis = cex.val, xlim = c(-7.5, -6.68)); #, breaks = 5

lines( density(Param1[,4]), col = "dark green", lty = 2, lwd = 3)
dev.off()


pdf( paste("INLA_beta4_births_",prefix,".pdf", sep = "") )
#hist(births.fit[,4], col = "blue", border = "white", main = "", xlab = "", ylab = "", lwd = 3, breaks = 9, cex.lab = cex.val, cex.axis = cex.val, xlim = c(-0.3, 0.49)); #, breaks = 5
plot( density(births.fit[1:100,5]), col = "blue", main = "", xlab = "", ylab = "", lwd = 3, cex.lab = cex.val, cex.axis = cex.val, xlim = c(-2.65, -1.0)); #, breaks = 5
abline(v = -2, col = "red", lwd = 3)

lines(  density(Births.coef[1:100,5]), col = "blue", main = "", xlab = "", ylab = "", lty = 2, lwd = 3, cex.lab = cex.val, cex.axis = cex.val, xlim = c(-7.5, -6.68)); #, breaks = 5

lines( density(Param1[,5]), col = "dark green", lty = 2, lwd = 3)
dev.off()

############### INLA PLOTS #########################


##########################################################################################
##########################################################################################
##########################################################################################

load("Setting_W.Rdat")

experiment <- "Experiment5"
N.ini <- 500
prefix <- "W2"
setting <- readRDS(paste(experiment,"/SimSetting_",N.ini,"_",prefix,"_all.RDS", sep = ""))

time.vec <- readRDS("time_vec.RDS")
boxplot(time.vec)

colMeans(time.vec)
apply(time.vec,2,sd)

Time.val.births.var <- readRDS(paste(experiment, "/Births_Cova_N0_",N0,"_",prefix,"_",setting$Kappa_fun,"_Time.RDS", sep = "") )
Time.val.deaths.var <- readRDS(paste(experiment, "/Deaths_Cova_N0_",N0,"_",prefix,"_",setting$Kappa_fun,"_Time.RDS", sep = "") )

Time.val.births.var <- readRDS(paste(experiment, "/Births_Cova_N0_",N0,"_",prefix,"_",setting$Kappa_fun,"_Time.RDS", sep = "") )
Time.val.deaths.var <- readRDS(paste(experiment, "/Deaths_Cova_N0_",N0,"_",prefix,"_",setting$Kappa_fun,"_Time.RDS", sep = "") )
time.vec <- readRDS("time_vec.RDS")

boxplot(time.vec)
range(time.vec)
boxplot(Time.val.deaths.var)

estim.time <- readRDS( paste(experiment,"/Both_fit_",N0,"_",prefix,"_",setting$Kappa_fun,"_time_estim.RDS", sep = ""))

boxplot( Time.val.deaths.var + do.call("c", Time.val.births.var) + rowSums(estim.time) )

our.time <- Time.val.deaths.var + do.call("c", Time.val.births.var) + rowSums(estim.time)

mean(time.vec);
sd(time.vec)

mean(our.time);
sd(our.time)


mean( Time.val.deaths.var + estim.time[,2])
sd( Time.val.deaths.var + estim.time[,2])

mean( do.call("c", Time.val.births.var) + estim.time[,1] )
sd( do.call("c", Time.val.births.var) + estim.time[,1] )
