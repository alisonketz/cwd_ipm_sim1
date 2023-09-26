s
##########################################################################################################
###
### Survival Simulation Generate Data and Format
###
##########################################################################################################

#set preliminaries
# setwd("~/Documents/Survival/surv_app/")
set.seed(10000)

###simulating age period data
source("../S4_01_age_period_survival_generate_data_type3.R")
ageperiod_out <- ageperiod_surv_sim_data_type3()

### Setting constants
nT_age <- max(ageperiod_out$right_age) - 1  #ageperiod_out$nT_age-1
nT_period <- max(ageperiod_out$right_period) - 1

### Setting true values for log hazard
beta0_true <- ageperiod_out$beta0
age_effect_true <- ageperiod_out$age_effect
period_effect_true <- ageperiod_out$period_effect

#Initialize vectors for sub-calculations of Survival
llambda_age_true <- rep(NA, nT_age)
UCH0_age_true <- rep(NA, nT_age)
S0_age_true <- rep(NA, nT_age)

for (t in 1:nT_age) {
  llambda_age_true[t] <- beta0_true + age_effect_true[t]
  UCH0_age_true[t] <- exp(llambda_age_true[t])
  S0_age_true[t] <- exp(-sum(UCH0_age_true[1:t]))
}

#Initialize vectors for sub-calculations of Survival
llambda_period_true <- rep(NA, nT_period)
UCH0_period_true <- rep(NA, nT_period)
S0_period_true <- rep(NA, nT_period)

for (t in 1:nT_period) {
  llambda_period_true[t] <- beta0_true + period_effect_true[t] #female
  UCH0_period_true[t] <- exp(llambda_period_true[t])
  S0_period_true[t] <- exp(-sum(UCH0_period_true[1:t]))
}

###############################
###
### formatting data
###
###############################

### Set data generated from the data generating function
left_age <- ageperiod_out$left_age
right_age <- ageperiod_out$right_age
left_period <- ageperiod_out$left_period
right_period <- ageperiod_out$right_period
rt_censor <- ageperiod_out$rt_censor
n <- ageperiod_out$n


### Formatting data for entering into the model
base <- rep(0, 5)
for (i in 1:n) {
  if (rt_censor[i] == 0) {
    temp1 <- c(left_age[i],
               right_age[i] - 1,
               1,
               left_period[i],
               right_period[i] - 1)
    temp2 <- c(right_age[i] - 1,
               right_age[i],
               0,
               right_period[i] - 1,
               right_period[i])
    if (left_age[i] == (right_age[i] - 1)) {
      base <- rbind(base, temp2)
    } else {
      base <- rbind(base, temp1, temp2)
    }
  } else {
    base <- rbind(base, c(left_age[i],
                         right_age[i],
                         1,
                         left_period[i],
                         right_period[i]))
  }
}
base <- base[-1, ]
rownames(base) <- NULL
colnames(base) <- c("left_age",
                    "right_age",
                    "censored",
                    "left_period",
                    "right_period")
df_fit <- as.data.frame(base)
n_fit <- dim(df_fit)[1]

### Create age to date conversion vector
### This aligns the age/period intervals when
### looping over the age and period hazards

age2date <- df_fit$left_age - df_fit$left_period