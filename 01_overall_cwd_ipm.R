###
### Alison C. Ketz 05/31/2023
###
###

###########################################################
### Preliminaries
###########################################################

rm(list = ls())

setwd("~/Documents/ipm/cwd_ipm_weekly_ew")

library(viridis)
library(RColorBrewer)
library(Hmisc)
library(lubridate)
library(readxl)
library(Matrix)
library(ggplot2)
library(gridExtra)
library(xtable)
library(nimble)
library(tidyverse)
library(dplyr)
library(lattice)
library(foreign)
library(rgdal)
library(rgeos)
library(zoo)
library(spdep)
library(parallel)
library(doParallel)
library(coda)
library(INLA)
library(sf)
library(terra)
library(splines)
library(MetBrewer)
library(ggforce)
library(tidyr)
library(ggridges)
library(ggh4x)

###########################################################
### Source summary function for posteriors
###########################################################

source("support_functions.R")

###########################################################
# Load Data
###########################################################

source("02_load_clean_data_foi.R")
source("03_load_clean_data_survival.R")

###############################################################
### Load/clean Spatially referenced data
###############################################################

source("04_homerange_data_survival.R")

###########################################################
### Format data for fitting age at harvest data
###########################################################

source("05_load_calibrate_data_aah.R")
source("06_load_agg_data_aah.R")

###########################################################
### Format data for fitting age-period survival models
###########################################################

source("07_format_data_survival.R")

##########################################################
### Setup collar data for FOI + Survival
##########################################################

source("08_format_data_combine_foi_surv.R")

##########################################################
### Setup collar data for FOI + Survival
##########################################################

source("09_format_data_cause.R")

###########################################################
### Setup consts etc for running the model
###########################################################

source("10_prelim_survival.R")
source("11_prelim_foi.R")
source("12_prelim_aah.R")

###########################################################
### Likelihoods
###########################################################

# source("13_distributions_check.R")
source("13_distributions.R")
# source("13_distributions_hunt.R")

###########################################################
### Functions for Efficient Calculations
###########################################################

source("14_calculations.R")

###########################################################
### Run model
###########################################################

source("15_modelcode.R")

###########################################################
### Run model
###########################################################

source("16_run_model.R")

# source("17_run_model_par.R")
# source("17_run_model_par_fun.R")

###########################################################
### Post processing
###########################################################

source("17_post_process.R")
