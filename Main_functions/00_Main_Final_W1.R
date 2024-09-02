#Main
rm(list = ls())
setwd("/data/fcp-X79-03/Spatio_temporal_PP/stBCIcode/SimulationStudy/")

experiment <- "Experiment5"
my.win <- "Window2"
if(my.win == "Window1") prefix <- "W1"
if(my.win == "Window2") prefix <- "W2"

Kappa_select <- "Kappa_1"
N0 <- 500
Time.steps <- 10
bw.seq <- seq(5,105, l = 11)
save.image("Setting_W.Rdat")

source("00_Simulation_BDH_model_time_Biv.R")


rm(list = ls())
load("Setting_W.Rdat")
source("07_Fitting_parameters_all_biv.R")
#source("07a_Fitting_parameters_all.R")

rm(list = ls())
load("Setting_W.Rdat")
source("08_Deaths_Variance_biv.R")
#source("08a_Deaths_Variance.R")

rm(list = ls())
load("Setting_W.Rdat")
source("09_Births_Variance_biv.R")
#source("09a_Births_Variance.R")

rm(list = ls())
load("Setting_W.Rdat")

#source("10_Coverage_full.R")
#source("10a_Coverage_full.R")
