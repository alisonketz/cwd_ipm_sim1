####################################################################
###
### Data and paramters needed to calculate likelihoods
### Fake parameters or loaded from runs w/o integrated models
###
####################################################################

beta0_survival_sus <- -8.75
beta0_survival_inf <- -8
beta_male <- .5

age_effect_survival_test <- 2*exp(-.09*seq(1:nT_age_surv))
age_effect_survival_test[600:nT_age_surv] <- .00001*exp(.009*seq(600:nT_age_surv))
age_effect_survival_test <- age_effect_survival_test- mean(age_effect_survival_test)

period_effect_survival_test <- 1 * sin(2/52 * pi * (1:nT_period_overall_ext)) + rnorm(nT_period_overall_ext,0,.1)
period_effect_survival_test <- period_effect_survival_test - mean(period_effect_survival_test)

##########################################################################
### FOI Parameters
### loading age and period effects from Transmission v3 w/o Survival
###########################################################################

# load("~/Documents/Transmission/Transmission_v5/01_urw2_ew_timeRW1/fit_sum.Rdata")
load("datafiles/fit_sum_foi.Rdata")

m_age_foi <- fit_sum[grep("m_age",rownames(fit_sum))[1:7],1]
f_age_foi <- fit_sum[grep("f_age",rownames(fit_sum))[1:7],1]
m_period_foi <- fit_sum[grep("m_time",rownames(fit_sum)),1]
f_period_foi <- fit_sum[grep("f_time",rownames(fit_sum)),1]
x <- 2002:2022
lmfperiod <- lm(f_period_foi~x)
lmmperiod <- lm(m_period_foi~x)
pred_foiperiod_f <- predict(lmfperiod,newdata=data.frame(x=1994:2001))
pred_foiperiod_m <- predict(lmmperiod,newdata=data.frame(x=1994:2001))
f_period_foi <- c(pred_foiperiod_f,f_period_foi)
m_period_foi <- c(pred_foiperiod_m,m_period_foi)
f_period_foi <- f_period_foi[-length(f_period_foi)]
m_period_foi <- m_period_foi[-length(m_period_foi)]




##########################################################################
### testing based on initial values for mcmc
###########################################################################



beta_male = rnorm(1, -.5, .01)
beta0_survival_sus = rnorm(1, -8, 0.01)
beta0_survival_inf = rnorm(1, -7, 0.0001)
ln_b_age_survival = rnorm(nknots_age) * 10^-4
b_period_survival = rnorm(nknots_period) * 10^-4

age_effect_survival_temp  <- c()
b_age_survival <- c()
age_effect_survival_test  <- c()

for(k in 1:nknots_age){
    b_age_survival[k] <- exp(ln_b_age_survival[k])
}
for (t in 1:nT_age_surv) {
    age_effect_survival_temp[t] <- inprod(b_age_survival[1:nknots_age],
                                     Z_age[t, 1:nknots_age])
  }
mu_age_effect_survival_temp <- mean(age_effect_survival_temp[1:nT_age_surv])

for (t in 1:nT_age_surv) {
    age_effect_survival_test[t] <-  age_effect_survival_temp[t] -
                               mu_age_effect_survival_temp
}

tau_period_survival = runif(1, .1, 1)
tau_age_survival = runif(1, .1, .4)
tau_age_foi_male = runif(1, 1.5, 1.7)
tau_age_foi_female = runif(1, 2.7, 4.2)
tau_period_foi_male = runif(1, 4.2, 6.8)
tau_period_foi_female = runif(1, 2.27, 3.44)
m_period_foi = seq(-.25, .25, length = n_year)
f_period_foi = seq(-.5, .5, length = n_year)
m_age_foi = c(rnorm(1, -6, sd = .1),
                rnorm(1, -5.5, sd = .1),
                rnorm(1, -5, sd = .1),
                rnorm(1, -5.5, sd = .1),
                rnorm(1, -6, sd = .1),
                rnorm(1, -7.2, sd = .1)) - 1.5
f_age_foi = c(rnorm(1, -6, sd = .1),
                rnorm(1, -5.5, sd = .1),
                rnorm(1, -6, sd = .1),
                rnorm(1, -6.5, sd = .1),
                rnorm(1, -6.8, sd = .1),
                rnorm(1, -7.2, sd = .1),
                rnorm(1, -8, sd = .1)) - 1.5
tau_period_precollar = rgamma(1, 1, 1)
period_annual_survival = rnorm(n_year_precollar + 1, .1)
beta0_cause = rnorm(1, -2.8, .1)
beta_cause_male = rnorm(1, 0, .1)
beta_cause_gun = rnorm(1, 1.5, .1)
beta_cause_ng = rnorm(1, 3, .1)
tau_obs = runif(1, .01, .5)
tau_pop = rgamma(2, 10, 1)
report_overall = report_overall_init
report = report_init
fec_epsilon = fec_eps_init
mu_fec = rnorm(1, mu_fec_init, .01)
fec_prec_eps = runif(1, 5, 10)
space_temp = rnorm(1, -.55, .01)
space_mix = 1
space = c(0, space_temp * space_mix)

# period_effect_survival_test <- rep(0,nT_period_overall_ext)
# age_effect_survival_test <- rep(0, nT_age_surv)


#Period effects for survival from collar data
period_effect_surv <- c()
for (t in 1:nT_period_collar) {
period_effect_surv[t] <- inprod(b_period_survival[1:nknots_period],
                                Z_period[t, 1:nknots_period])
}

#Period effects from aah data
# period_effect_survival <- rep(NA,nT_period_overall_ext)

period_effect_survival_test <- set_period_effects_constant(
    n_year_precollar = n_year_precollar,
    n_year_precollar_ext = n_year_precollar_ext,
    n_year_prestudy_ext = n_year_prestudy_ext,
    nT_period_precollar_ext = nT_period_precollar_ext,
    nT_period_precollar = nT_period_precollar,
    nT_period_collar = nT_period_collar,
    nT_period_overall_ext = nT_period_overall_ext,
    nT_period_prestudy_ext = nT_period_prestudy_ext,
    yr_start = d_fit_season$yr_start[1:n_year],
    yr_end = d_fit_season$yr_end[1:n_year],
    period_effect_surv = period_effect_surv[1:nT_period_collar],
    period_annual_survival = period_annual_survival[1:(n_year_precollar + 1)]
)

#checking initial values of population size breakdown based on initial log vectors

#should this be different for pos/neg or m/f for study area?
# for(i in 1:2) {
# tau_pop[i] ~ dgamma(1, 1)
# }

llpop_sus  <- pop_sus <- array(NA,c(2,2,n_agef,n_year))
llpop_inf  <- pop_inf <- array(NA,c(2,2,n_agef,n_year))
for(k in 1:n_study_area) {
    for (a in 1:n_agef) {

    #Initial population structure pop[sex,age,year] for susceptible deer
    llpop_sus[k, 1, a, 1]  <- rnorm(1,f_logpop_sus[k, a], 1/sqrt(tau_pop[1]))#tau_pop[1]
    pop_sus[k, 1, a, 1] <- exp(llpop_sus[k, 1, a, 1])

    #Initial population structure pop[study_area=k,sex=i,year=t,age=a]
    llpop_inf[k, 1, a, 1] <- rnorm(1, f_logpop_inf[k, a], 1/sqrt(tau_pop[2]))#tau_pop[1]
    pop_inf[k, 1, a, 1] <- exp(llpop_inf[k, 1, a, 1])
    }
    for (a in 1:n_agem) {
        ### Initial population structure pop
        ### [study_area,sex,age,period(year)] for susceptible deer
        llpop_sus[k, 2, a, 1]  <- rnorm(1, m_logpop_sus[k, a],  1/sqrt(tau_pop[1]))#tau_pop[2]
        pop_sus[k, 2, a, 1] <- exp(llpop_sus[k, 2, a, 1])

        #Initial population structure pop for infected deer
        llpop_inf[k, 2, a, 1] <- rnorm(1,m_logpop_inf[k, a],  1/sqrt(tau_pop[2]))#tau_pop[2]
        pop_inf[k, 2, a, 1] <- exp(llpop_inf[k, 2, a, 1])
    }
}

pop_sus_rand <- pop_sus
pop_inf_rand <- pop_inf
round(pop_sus_rand[1,1,,1])
round(pop_inf_rand[1,1,,1])


llpop_sus  <- pop_sus <- array(NA,c(2,2,n_agef,n_year))
llpop_inf  <- pop_inf <- array(NA,c(2,2,n_agef,n_year))
for(k in 1:n_study_area) {
    for (a in 1:n_agef) {

    #Initial population structure pop[sex,age,year] for susceptible deer
    llpop_sus[k, 1, a, 1]  <- f_logpop_sus[k, a]
    pop_sus[k, 1, a, 1] <- exp(llpop_sus[k, 1, a, 1])

    #Initial population structure pop[study_area=k,sex=i,year=t,age=a]
    llpop_inf[k, 1, a, 1] <- f_logpop_inf[k, a]
    pop_inf[k, 1, a, 1] <- exp(llpop_inf[k, 1, a, 1])
    }
    for (a in 1:n_agem) {
        ### Initial population structure pop
        ### [study_area,sex,age,period(year)] for susceptible deer
        llpop_sus[k, 2, a, 1]  <- m_logpop_sus[k, a]
        pop_sus[k, 2, a, 1] <- exp(llpop_sus[k, 2, a, 1])

        #Initial population structure pop for infected deer
        llpop_inf[k, 2, a, 1] <- m_logpop_inf[k, a]
        pop_inf[k, 2, a, 1] <- exp(llpop_inf[k, 2, a, 1])
    }
}

round(pop_sus_rand[1,1,,1])
round(pop_sus[1,1,,1])
round(pop_inf_rand[1,1,,1])
round(pop_inf[1,1,,1])