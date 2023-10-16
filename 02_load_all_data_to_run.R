
#######################################################################
###
### loading results from working data model
### to set true values for the mortality and infection hazards
###
#######################################################################

load("datafiles/fit_sum_finaldat.Rdata")

n_year  <- 5

n_ageclassm <- 6
n_ageclassf <- 7
n_ageclass <- n_ageclassm

n_agef <- 10
n_agem <- 7
n_age <- n_agem

nT_age_short_f <- 52 * (n_agef - 1) + 2
nT_age_surv_aah_f <- 52 * n_agef + 2
nT_age_short_m <- 52 * (n_agem - 1) + 1
nT_age_surv_aah_m <- 52 * n_agem + 1 

nT_age_short <- nT_age_short_m
nT_age_surv_aah <- nT_age_surv_aah_m

period_effect_true <- fit_sum[grep("period_effect_surv", rownames(fit_sum)), 1]
nT_period <- length(period_effect_true)
age_effect_true <- fit_sum[ grep("age_effect", rownames(fit_sum)), 1]
nT_age <- length(age_effect_true)
beta0_survival_sus_true <- fit_sum[grep("beta0_survival_sus",
                                   rownames(fit_sum)), 1] 

beta0_survival_inf_true <- fit_sum[grep("beta0_survival_inf",
                                   rownames(fit_sum)), 1]

beta_male_true <- fit_sum[grep("beta_male",
                          rownames(fit_sum)), 1]

f_age_foi <- fit_sum[grep("f_age_foi",
                         rownames(fit_sum)), 1][1:n_ageclassf]

m_age_foi <- fit_sum[grep("m_age_foi",
                         rownames(fit_sum)), 1][1:n_ageclassm]

f_age_foi_true <- c(rep(f_age_foi[1:4], each = 52),
                    rep(f_age_foi[5], each = 2 * 52),
                    rep(f_age_foi[6], each = 3 * 52))

f_age_foi_true <- c(f_age_foi_true,rep(f_age_foi[7],
                    nT_age - length(f_age_foi_true)))

m_age_foi_true <- c(rep(m_age_foi[1:4], each = 52),
                    rep(m_age_foi[5], each = 2 * 52),
                    rep(m_age_foi[6], each = 3 * 52))

m_age_foi_true <- c(m_age_foi_true,rep(m_age_foi[6],
                    nT_age - length(m_age_foi_true)))

age_foi <- m_age_foi

#scale these for the males, since using male foi
beta0_survival_sus_true <- beta0_survival_sus_true + beta_male_true
beta0_survival_inf_true <- beta0_survival_inf_true + beta_male_true


##################################################################
###
### setting params to run within
### the function, only run when developing simulation function
###
##################################################################

beta0_survival_sus <- beta0_survival_sus_true
beta0_survival_inf <- beta0_survival_inf_true
foi_age_effect <- m_age_foi_true
age_effect <- age_effect_true
period_effect <- period_effect_true