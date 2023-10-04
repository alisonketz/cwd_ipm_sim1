##################################################################
###
### run the function to generate data
###
##################################################################

processnum <- 0
set.seed(10000 + processnum)

source("03_fun_generate_data.R")

dat <- ageperiod_surv_foi_sim_data(
            beta0_survival_sus = beta0_survival_sus,
            beta0_survival_inf = beta0_survival_inf,
            foi_age_effect = foi_age_effect,
            age_effect = age_effect_true,
            period_effect = period_effect_true,
            nT_age = nT_age,
            nT_period = nT_period,
            processnum = processnum
            )

# table(dat$pos1)
# table(dat$pos2)
# hist(dat$right_age,breaks = 100)
# hist(dat$right_age)
# sum(dat$rt_censor)
# n_ind

# hist(dat$right_period,breaks = 100)
# dat$right_age
# # prevalence
# sum(!is.na(dat$inf_age))/n_ind
# dat$inf_status[1,]
# which(dat$inf_age[!is.na(inf_age)] < dat$left_age[!is.na(inf_age)])


### Set data generated from the data generating function
# left_age <- dat$left_age
# right_age <- dat$right_age
# left_period <- dat$left_period
# right_period <- dat$right_period
# rt_censor <- dat$rt_censor
# n_ind <- dat$n_ind
# inf_age <- dat$inf_age

##############################################################
### there are too many deer infected at to capture
### because we aren't accounting for 
### disease-associated mortality in the generating function
##############################################################

n_remove <- sum(dat$pos1) - round(sum(dat$pos1)/5)
png("hist_icap_leftage_toweight.png")
hist(dat$left_age[dat$pos1==1],breaks = 100)
dev.off()


indx_remove <- sample(which(dat$pos1 == TRUE),
                      n_remove,
                      replace = FALSE,
                      prob = dat$left_age[dat$pos1==1]/sum(dat$left_age[dat$pos1==1])
                      ) 
n_ind <- dat$n_ind - n_remove


##############################################################
###
### setup overall data frame
###
##############################################################

df_fit <- data.frame(id = 1:n_ind,
                     left_age = dat$left_age[-indx_remove],
                     right_age = dat$right_age[-indx_remove],
                     left_period = dat$left_period[-indx_remove],
                     right_period = dat$right_period[-indx_remove],
                     rt_censor = dat$rt_censor[-indx_remove],
                     cwd_cap = dat$pos1[-indx_remove],
                     cwd_mort = dat$pos2[-indx_remove],
                     inf_age = dat$inf_age[-indx_remove]
                     )

df_fit$inf_period = df_fit$right_age - df_fit$inf_age + df_fit$left_period
n_fit <- nrow(df_fit)
# hist(df_fit$left_period, breaks = 100)
# hist(df_fit$inf_period[df_fit$inf_period>0],breaks = 100)
# sum(df_fitinf_cap)
# sum(1-df_fit$rt_censor)
# sum(df_fit$rt_censor)


### setting up in terms of e/r/s
df_fit$e_age <- df_fit$left_age
df_fit$r_age <- df_fit$right_age
df_fit$s_age <- df_fit$right_age
df_fit$r_age[df_fit$rt_censor == 0] <- df_fit$r_age[df_fit$rt_censor == 0] - 1
df_fit$s_age[df_fit$rt_censor == 1] <- NA


df_fit$e_period <- df_fit$left_period
df_fit$r_period <- df_fit$right_period
df_fit$s_period <- df_fit$right_period
df_fit$r_period[df_fit$rt_censor == 0] <- df_fit$r_period[df_fit$rt_censor == 0] - 1
df_fit$s_period[df_fit$rt_censor == 1] <- NA


########################################################################
### calibrating age of deer with the study time 
### for indexing in the likelihood for loops
### age2date = left_period - left_age
########################################################################

df_fit$age2date <- df_fit$e_period - df_fit$e_age