rm(list = ls())
Coverage.results <- function(experiment, prefix, N.ini = 500){     
  setting <- readRDS(paste(experiment,"/SimSetting_",N.ini,"_",prefix,"_all.RDS", sep = ""))
  W.sim <- setting$window
  N0 <- setting$N0
  bw.seq <- setting$bw.seq
  Time.steps <- setting$Time.steps
  beta0b <- setting$beta.vec[1]
  beta1b <- setting$beta.vec[2]
  beta1k <- setting$beta.vec[3]
  
  results <- list()
  
  Births.fit <- readRDS(file = paste(experiment, "/Births_fit_",N0,"_",prefix,"_",setting$Kappa_fun,"_all_biv.RDS", sep = ""))
  Births.cova <- readRDS(file = paste(experiment, "/Births_Cova_N0_",N0,"_",prefix,"_",setting$Kappa_fun,"_all_biv.RDS", sep = ""))
  
  Births.coef <- do.call("rbind", lapply(Births.fit, function(x) x$coef ) )
  births.mean <- colMeans(Births.coef)
  
  ######################################
  ############### Births ###############
  ######################################
  
  Cova.B.11 <- do.call("rbind", lapply(Births.cova, function(x) do.call("c", lapply(x, function(z) z[1,1]) ) ) )
  Cova.B.22 <- do.call("rbind", lapply(Births.cova, function(x) do.call("c", lapply(x, function(z) z[2,2]) ) ) )
  Cova.B.33 <- do.call("rbind", lapply(Births.cova, function(x) do.call("c", lapply(x, function(z) z[3,3]) ) ) )
  Cova.B.44 <- do.call("rbind", lapply(Births.cova, function(x) do.call("c", lapply(x, function(z) z[4,4]) ) ) )
  Cova.B.55 <- do.call("rbind", lapply(Births.cova, function(x) do.call("c", lapply(x, function(z) z[5,5]) ) ) )
  
  Births.cova.val <- cov(Births.coef)
  
  results$Births <- list()
  results$Births$mean <- births.mean
  results$Births <- list(Cova1 = Cova.B.11, 
                         Cova2 = Cova.B.22, 
                         Cova3 = Cova.B.33,
                         Cova4 = Cova.B.44,
                         Cova5 = Cova.B.55,
                         Cova.val = Births.cova.val)
  
  
  CI.births <- lapply(as.list(1:length(bw.seq)), function(i){
    v1 <- sum( (Births.coef[,1] - sqrt(Cova.B.11[,i])*qnorm(0.975)) <= beta0b & beta0b <= (Births.coef[,1] + sqrt(Cova.B.11[,i])*qnorm(0.975)) )/1000
    v2 <- sum( (Births.coef[,2] - sqrt(Cova.B.22[,i])*qnorm(0.975)) <= 0 & 0 <= (Births.coef[,2] + sqrt(Cova.B.22[,i])*qnorm(0.975)) )/1000
    v3 <- sum( (Births.coef[,3] - sqrt(Cova.B.33[,i])*qnorm(0.975)) <= beta1b & beta1b <= (Births.coef[,3] + sqrt(Cova.B.33[,i])*qnorm(0.975)) )/1000
    v4 <- sum( (Births.coef[,4] - sqrt(Cova.B.44[,i])*qnorm(0.975)) <= beta1k & beta1k <= (Births.coef[,4] + sqrt(Cova.B.44[,i])*qnorm(0.975)) )/1000
    v5 <- sum( (Births.coef[,5] - sqrt(Cova.B.55[,i])*qnorm(0.975)) <= -2 & -2 <= (Births.coef[,5] + sqrt(Cova.B.55[,i])*qnorm(0.975)) )/1000
    return(c(v1, v2, v3, v4, v5))
  })
  
  CoverageB0 <- do.call("c",lapply(CI.births, function(x) x[1]))
  CoverageB1 <- do.call("c",lapply(CI.births, function(x) x[2]))
  CoverageB2 <- do.call("c",lapply(CI.births, function(x) x[3]))
  CoverageB3 <- do.call("c",lapply(CI.births, function(x) x[4]))
  CoverageB4 <- do.call("c",lapply(CI.births, function(x) x[5]))
  
  results$Births$coverageB0 <- CoverageB0
  results$Births$coverageB1 <- CoverageB1
  results$Births$coverageB2 <- CoverageB2
  results$Births$coverageB3 <- CoverageB3
  results$Births$coverageB4 <- CoverageB4
  
  rm(Births.fit)
  rm(Births.cova)
  
  ######################################
  ############### Deaths ###############
  ######################################
  etha0 <- setting$etha.vec[1]
  etha1 <- setting$etha.vec[2]
  
  #Loading parameters estimates
  Deaths.coef <- readRDS(file = paste(experiment, "/Deaths_fit_",N0,"_",prefix,"_",setting$Kappa_fun,"_all_biv.RDS", sep = ""))
  Deaths.cova <- readRDS(paste(experiment, "/Deaths_Cova_N0_",N0,"_",prefix,"_",setting$Kappa_fun,"_all_biv.RDS", sep = ""))
  Deaths.mean <- colMeans(Deaths.coef)
  
  #Loading variance of the parameters estimates
  
  Cova.D.11 <- do.call("rbind", lapply(Deaths.cova, function(x) do.call("c", lapply(x, function(z) z[1,1]) ) ) )
  Cova.D.22 <- do.call("rbind", lapply(Deaths.cova, function(x) do.call("c", lapply(x, function(z) z[2,2]) ) ) )
  Cova.D.33 <- do.call("rbind", lapply(Deaths.cova, function(x) do.call("c", lapply(x, function(z) z[3,3]) ) ) )
  Cova.D.44 <- do.call("rbind", lapply(Deaths.cova, function(x) do.call("c", lapply(x, function(z) z[4,4]) ) ) )
  Cova.D.55 <- do.call("rbind", lapply(Deaths.cova, function(x) do.call("c", lapply(x, function(z) z[5,5]) ) ) )
  
  Deaths.cova.val <- cov(Deaths.coef)
  
  results$Deaths <- list(Cova1 = Cova.D.11, 
                         Cova2 = Cova.D.22,
                         Cova3 = Cova.D.33,
                         Cova4 = Cova.D.44,
                         Cova5 = Cova.D.55,
                         Cova.val = Deaths.cova.val)
  
  results$Deaths$mean <- Deaths.mean
  
  CI.Deaths <- lapply(as.list(1:length(bw.seq)), function(i){
    v1 <- sum( (Deaths.coef[,1] - sqrt(Cova.D.11[,i])*qnorm(0.975)) <= etha0 & etha0 <= (Deaths.coef[,1] + sqrt(Cova.D.11[,i])*qnorm(0.975)) )/1000
    v2 <- sum( (Deaths.coef[,2] - sqrt(Cova.D.22[,i])*qnorm(0.975)) <= etha1 & etha1 <= (Deaths.coef[,2] + sqrt(Cova.D.22[,i])*qnorm(0.975)) )/1000
    v3 <- sum( (Deaths.coef[,3] - sqrt(Cova.D.33[,i])*qnorm(0.975)) <= 0 & 0 <= (Deaths.coef[,3] + sqrt(Cova.D.33[,i])*qnorm(0.975)) )/1000
    v4 <- sum( (Deaths.coef[,4] - sqrt(Cova.D.44[,i])*qnorm(0.975)) <= -0.25 & -0.25 <= (Deaths.coef[,4] + sqrt(Cova.D.44[,i])*qnorm(0.975)) )/1000
    v5 <- sum( (Deaths.coef[,5] - sqrt(Cova.D.55[,i])*qnorm(0.975)) <= 0 & 0 <= (Deaths.coef[,5] + sqrt(Cova.D.55[,i])*qnorm(0.975)) )/1000
    return(c(v1, v2, v3, v4, v5))
  })
  
  CoverageB0 <- do.call("c",lapply(CI.Deaths, function(x) x[1]))
  CoverageB1 <- do.call("c",lapply(CI.Deaths, function(x) x[2]))
  CoverageB2 <- do.call("c",lapply(CI.Deaths, function(x) x[3]))
  CoverageB3 <- do.call("c",lapply(CI.Deaths, function(x) x[4]))
  CoverageB4 <- do.call("c",lapply(CI.Deaths, function(x) x[5]))
  
  results$Deaths$coverageB0 <- CoverageB0
  results$Deaths$coverageB1 <- CoverageB1
  results$Deaths$coverageB2 <- CoverageB2
  results$Deaths$coverageB3 <- CoverageB3
  results$Deaths$coverageB4 <- CoverageB4
  results$bw.seq <- bw.seq
  results$N0 <- N0
  rm(Deaths.coef)
  rm(Deaths.cova)
  
  return(results)
}

experiment <- "Experiment5"
N.ini <- 500
prefix <- "W1"
Results.experiments <- Coverage.results(experiment, prefix, N.ini = N.ini)

setting <- readRDS(paste(experiment,"/SimSetting_",N.ini,"_",prefix,"_all.RDS", sep = ""))

cex.val <- 1.75
pdf( paste("Variance_estim_beta0_births_",prefix,".pdf", sep = "") )
boxplot(Results.experiments$Births$Cova1, type = "l", names = round(Results.experiments$bw.seq, 2), ylab = "", cex.lab = cex.val, cex.axis = cex.val);
points( rep(Results.experiments$Births$Cova.val[1,1], length(Results.experiments$bw.seq)), col = "red", pch = 16)
dev.off()

pdf( paste("Variance_estim_beta1_births_",prefix,".pdf", sep = "") )
boxplot(Results.experiments$Births$Cova2, type = "l", names = round(Results.experiments$bw.seq, 2), ylab = "", cex.lab = cex.val, cex.axis = cex.val);
points( rep(Results.experiments$Births$Cova.val[2,2], length(Results.experiments$bw.seq)), col = "red", pch = 16)
dev.off()

pdf( paste("Variance_estim_beta2_births_",prefix,".pdf", sep = "") )
boxplot(Results.experiments$Births$Cova3, type = "l", names = round(Results.experiments$bw.seq, 2), ylab = "", cex.lab = cex.val, cex.axis = cex.val);
points( rep(Results.experiments$Births$Cova.val[3,3], length(Results.experiments$bw.seq)), col = "red", pch = 16)
dev.off()

pdf( paste("Variance_estim_beta3_births_",prefix,".pdf", sep = "") )
boxplot(Results.experiments$Births$Cova4, type = "l", names = round(Results.experiments$bw.seq, 2), ylab = "", cex.lab = cex.val, cex.axis = cex.val);
points( rep(Results.experiments$Births$Cova.val[4,4], length(Results.experiments$bw.seq)), col = "red", pch = 16)
dev.off()

pdf( paste("Variance_estim_beta4_births_",prefix,".pdf", sep = "") )
boxplot(Results.experiments$Births$Cova5, type = "l", names = round(Results.experiments$bw.seq, 2), ylab = "", cex.lab = cex.val, cex.axis = cex.val);
points( rep(Results.experiments$Births$Cova.val[5,5], length(Results.experiments$bw.seq)), col = "red", pch = 16)
dev.off()

format( round( diag(Results.experiments$Births$Cova.val), 6), scientific = TRUE)

pdf( paste("Coverage_births_",prefix,".pdf", sep = "") )
matplot(Results.experiments$bw.seq, cbind(Results.experiments$Births$coverageB0, 
                                          Results.experiments$Births$coverageB1, 
                                          Results.experiments$Births$coverageB2,
                                          Results.experiments$Births$coverageB3,
                                          Results.experiments$Births$coverageB4), type = "l", lty = 1 ,lwd = 2, 
        col = c("red", "blue", "green", "purple", "orange"), main = "", ylim = c(0.85,1), ylab = "", xlab = "" , cex.lab = cex.val, cex.axis = cex.val)
abline(h = 0.95, lwd = 2); 

legend("bottom", c(expression(beta[0][b], beta[1][b], beta[2][b], gamma[1][b], gamma[2][b])), 
                    lty = 1, lwd = 2, 
                    col = c("red", "blue","green","purple", "orange", "white"), cex = 1.5, ncol=2)
dev.off()

pdf( paste("Variance_estim_beta0_deaths_",prefix,".pdf", sep = "") )
boxplot(Results.experiments$Deaths$Cova1, type = "l", names = round(Results.experiments$bw.seq, 2), ylab = "" , cex.lab = cex.val, cex.axis = cex.val);
points( rep(Results.experiments$Deaths$Cova.val[1,1], length(Results.experiments$bw.seq)), col = "red", pch = 16)
dev.off()

pdf( paste("Variance_estim_beta1_deaths_",prefix,".pdf", sep = "") )
boxplot(Results.experiments$Deaths$Cova2, type = "l", names = round(Results.experiments$bw.seq, 2), ylab = "" , cex.lab = cex.val, cex.axis = cex.val);
points( rep(Results.experiments$Deaths$Cova.val[2,2], length(Results.experiments$bw.seq)), col = "red", pch = 16)
dev.off()

pdf( paste("Variance_estim_beta2_deaths_",prefix,".pdf", sep = "") )
boxplot(Results.experiments$Deaths$Cova3, type = "l", names = round(Results.experiments$bw.seq, 2), ylab = "" , cex.lab = cex.val, cex.axis = cex.val);
points( rep(Results.experiments$Deaths$Cova.val[3,3], length(Results.experiments$bw.seq)), col = "red", pch = 16)
dev.off()

pdf( paste("Variance_estim_beta3_deaths_",prefix,".pdf", sep = "") )
boxplot(Results.experiments$Deaths$Cova4, type = "l", names = round(Results.experiments$bw.seq, 2), ylab = "" , cex.lab = cex.val, cex.axis = cex.val);
points( rep(Results.experiments$Deaths$Cova.val[4,4], length(Results.experiments$bw.seq)), col = "red", pch = 16)
dev.off()

pdf( paste("Variance_estim_beta4_deaths_",prefix,".pdf", sep = "") )
boxplot(Results.experiments$Deaths$Cova5, type = "l", names = round(Results.experiments$bw.seq, 2), ylab = "" , cex.lab = cex.val, cex.axis = cex.val);
points( rep(Results.experiments$Deaths$Cova.val[5,5], length(Results.experiments$bw.seq)), col = "red", pch = 16)
dev.off()

pdf( paste("Coverage_deaths_",prefix,".pdf", sep = "") )
matplot(Results.experiments$bw.seq, cbind(Results.experiments$Deaths$coverageB0, 
                                          Results.experiments$Deaths$coverageB1,
                                          Results.experiments$Deaths$coverageB2,
                                          Results.experiments$Deaths$coverageB3,
                                          Results.experiments$Deaths$coverageB4), type = "l", lty = 1 ,lwd = 2, 
        col = c("red", "blue","green","purple", "orange"), ylim = c(0.85,1), ylab = "", xlab = "" , cex.lab = cex.val, cex.axis = cex.val) 
abline(h = 0.95, lwd = 2); 
legend("bottom", c(expression(beta[0][d], beta[1][d], beta[2][d], gamma[1][d], gamma[2][d])), 
       lty = 1, lwd = 2, 
       col = c("red", "blue","green","purple", "orange", "white"), cex = 1.5, ncol=2)
dev.off()



N0 <- 500
prefix <- "W2"
setting <- readRDS(paste(experiment,"/SimSetting_",N.ini,"_",prefix,"_all.RDS", sep = ""))
Births.fit <- readRDS(file = paste(experiment, "/Births_fit_",N0,"_",prefix,"_",setting$Kappa_fun,"_all_biv.RDS", sep = ""))
Deaths.fit <- readRDS(file = paste(experiment, "/Deaths_fit_",N0,"_",prefix,"_",setting$Kappa_fun,"_all_biv.RDS", sep = ""))

births.fit <- do.call("rbind", lapply(Births.fit, function(x) x$coef)); rm(Births.fit)
deaths.fit <- Deaths.fit

round(apply(births.fit, 2, sd),3)
round(apply(deaths.fit, 2, sd),3)

cor.w2.b <- cor(births.fit)
cor.w2.d <- cor(deaths.fit)


pdf( paste("beta0_deaths_",prefix,".pdf", sep = "") )
#hist(deaths.fit[,1], col = "blue", border = "white", main = "", xlab = "", ylab = "", lwd = 3, breaks = 9, cex.lab = cex.val, cex.axis = cex.val, xlim = c(-0.45,-0.074)); #
plot( density(deaths.fit[,1]), col = "blue", main = "", xlab = "", ylab = "", lwd = 3, cex.lab = cex.val, cex.axis = cex.val, xlim = c(-1.0,0.5)); #
abline(v = setting$etha.vec[1], col = "red", lwd = 3)
abline(v = mean(deaths.fit[,1]), col = "green", lwd = 3)
dev.off()


pdf( paste("beta1_deaths_",prefix,".pdf", sep = "") )
#hist(deaths.fit[,2], col = "blue", border = "white", main = "", xlab = "", ylab = "", lwd = 3, breaks = 9, cex.lab = cex.val, cex.axis = cex.val, xlim = c(-0.12,0.74)); #, breaks = 5
plot( density(deaths.fit[,2]), col = "blue", main = "", xlab = "", ylab = "", lwd = 3, cex.lab = cex.val, cex.axis = cex.val, xlim = c(-0.40,0.8)); #, breaks = 5
abline(v = 0.25, col = "red", lwd = 3)
abline(v = mean(deaths.fit[,2]), col = "green", lwd = 3)
dev.off()

pdf( paste("beta2_deaths_",prefix,".pdf", sep = "") )
#hist(deaths.fit[,3], col = "blue", border = "white", main = "", xlab = "", ylab = "", lwd = 3, breaks = 9, cex.lab = cex.val, cex.axis = cex.val, xlim = c(-0.33,0.39)); #breaks = 5,
plot( density(deaths.fit[,3]), col = "blue", main = "", xlab = "", ylab = "", lwd = 3, cex.lab = cex.val, cex.axis = cex.val, xlim = c(-0.6,0.60)); #breaks = 5,
abline(v = 0, col = "red", lwd = 3)
abline(v = mean(deaths.fit[,3]), col = "green", lwd = 3)
dev.off()

pdf( paste("beta3_deaths_",prefix,".pdf", sep = "") )
#hist(deaths.fit[,4], col = "blue", border = "white", main = "", xlab = "", ylab = "", lwd = 3, breaks = 9, cex.lab = cex.val, cex.axis = cex.val, xlim = c(-0.12,0.12)); #breaks = 5, 
plot( density(deaths.fit[,4]), col = "blue", main = "", xlab = "", ylab = "", lwd = 3, cex.lab = cex.val, cex.axis = cex.val, xlim = c(-0.5,0.0)); #breaks = 5, 
abline(v = -0.25, col = "red", lwd = 3)
abline(v = mean(deaths.fit[,4]), col = "green", lwd = 3)
dev.off()
#

pdf( paste("beta4_deaths_",prefix,".pdf", sep = "") )
#hist(deaths.fit[,4], col = "blue", border = "white", main = "", xlab = "", ylab = "", lwd = 3, breaks = 9, cex.lab = cex.val, cex.axis = cex.val, xlim = c(-0.12,0.12)); #breaks = 5, 
plot( density(deaths.fit[,5]), col = "blue", main = "", xlab = "", ylab = "", lwd = 3, cex.lab = cex.val, cex.axis = cex.val, xlim = c(-0.25,0.25)); #breaks = 5, 
abline(v = 0, col = "red", lwd = 3)
abline(v = mean(deaths.fit[,5]), col = "green", lwd = 3)
dev.off()

pdf( paste("beta0_births_",prefix,".pdf", sep = "") )
#hist(births.fit[,1], col = "blue", border = "white", main = "", xlab = "", ylab = "", lwd = 3, breaks = 9, cex.lab = cex.val, cex.axis = cex.val, xlim = c(-6.4, -6.06)); #, breaks = 5
plot(  density(births.fit[,1]), col = "blue", main = "", xlab = "", ylab = "", lwd = 3, cex.lab = cex.val, cex.axis = cex.val, xlim = c(-7.2, -6.68)); #, breaks = 5
abline(v = setting$beta.vec[1], col = "red", lwd = 3)
abline(v = mean(births.fit[,1]), col = "green", lwd = 3)
dev.off()

pdf( paste("beta1_births_",prefix,".pdf", sep = "") )
#hist(births.fit[,2], col = "blue", border = "white", main = "", xlab = "", ylab = "", lwd = 3, breaks = 9, cex.lab = cex.val, cex.axis = cex.val, xlim = c(-0.24, 0.23)); #, breaks = 5
plot( density(births.fit[,2]), col = "blue", main = "", xlab = "", ylab = "", lwd = 3, cex.lab = cex.val, cex.axis = cex.val, xlim = c(-0.37, 0.33)); #, breaks = 5
abline(v = 0, col = "red", lwd = 3)
abline(v = mean(births.fit[,2]), col = "green", lwd = 3)
dev.off()

pdf( paste("beta2_births_",prefix,".pdf", sep = "") )
#hist(births.fit[,3], col = "blue", border = "white", main = "", xlab = "", ylab = "", lwd = 3, breaks = 9, cex.lab = cex.val, cex.axis = cex.val, xlim = c(-0.15,0.38)); #, breaks = 5
plot( density(births.fit[,3]), col = "blue", main = "", xlab = "", ylab = "", lwd = 3, cex.lab = cex.val, cex.axis = cex.val, xlim = c(-0.30,0.47)); #, breaks = 5
abline(v = 0.1, col = "red", lwd = 3)
abline(v = mean(births.fit[,3]), col = "green", lwd = 3)
dev.off()

pdf( paste("beta3_births_",prefix,".pdf", sep = "") )
#hist(births.fit[,4], col = "blue", border = "white", main = "", xlab = "", ylab = "", lwd = 3, breaks = 9, cex.lab = cex.val, cex.axis = cex.val, xlim = c(-0.3, 0.49)); #, breaks = 5
plot( density(births.fit[,4]), col = "blue", main = "", xlab = "", ylab = "", lwd = 3, cex.lab = cex.val, cex.axis = cex.val, xlim = c(-0.8, 0.72)); #, breaks = 5
abline(v = 0.1, col = "red", lwd = 3)
abline(v = mean(births.fit[,4]), col = "green", lwd = 3)
dev.off()

pdf( paste("beta4_births_",prefix,".pdf", sep = "") )
#hist(births.fit[,4], col = "blue", border = "white", main = "", xlab = "", ylab = "", lwd = 3, breaks = 9, cex.lab = cex.val, cex.axis = cex.val, xlim = c(-0.3, 0.49)); #, breaks = 5
plot( density(births.fit[,5]), col = "blue", main = "", xlab = "", ylab = "", lwd = 3, cex.lab = cex.val, cex.axis = cex.val, xlim = c(-3.0, -1.23)); #, breaks = 5
abline(v = -2, col = "red", lwd = 3)
abline(v = mean(births.fit[,5]), col = "green", lwd = 3)
dev.off()


# colMeans(Births.count.W2)/colMeans(Births.count.W1)
# 
# 
# # 
# experiment <- "Experiment5"
# prefix <- "W1"
# N.ini <- 750
# bw.seq <- seq(5,155, l = 16)
# Results.experiments.W1 <- Coverage.results(experiment, prefix, N.ini = N.ini)
# prefix <- "W2"
# Results.experiments.W2 <- Coverage.results(experiment, prefix, N.ini = N.ini)
# 
# sqrt(diag( Results.experiments.W1$Births$Cova.val )/diag( Results.experiments.W2$Births$Cova.val ))
# sqrt(diag( Results.experiments.W1$Deaths$Cova.val )/diag( Results.experiments.W2$Deaths$Cova.val ))


