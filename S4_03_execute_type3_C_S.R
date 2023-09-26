######################################################################################################################3
###
### Load Dependent Libraries
###
######################################################################################################################3

library(nimble)
library(Matrix)
library(splines)
library(coda)
library(ggplot2)

########################################################
###
### Load the data generating function,
### Generate the data,
### Format the data for running in these models
###
###
########################################################

source("../S4_01_age_period_survival_generate_data_type3.R")

########################################################
### Generate the data and format to run in models
########################################################

source("../S4_02_format_data_type3.R")

########################################################
### Set constants needed to fit the model
########################################################

source("S4_04_prelim_consts_type3_C_S.R")

########################################################
### Run the model
########################################################

source("S4_05_run_model_type3_C_S.R")

########################################################
### Summary statistics, plot results
########################################################

source("S4_06_post_plots_type3_C_S.R")
