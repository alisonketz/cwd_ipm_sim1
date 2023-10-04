############################################################################################
############################################################################################
############################################################################################
###
### Model Statement
###
############################################################################################
############################################################################################
############################################################################################

modelcode <- nimbleCode({

  ##############################
  ### Priors
  ##############################

  ##############################
  ### Force of infection model
  ##############################

  tau_age_foi  ~ dgamma(1, 1)
  tau1_age_foi <- .0000001 * tau_age_foi
  foi_age_effect[1] ~ dnorm(0, tau1_age_foi)
  foi_age_effect[2] ~ dnorm(0, tau1_age_foi)
  for (i in 3:n_ageclass) {
    foi_age_effect[i]~dnorm(2 * foi_age_effect[i-1] - foi_age_effect[i-2], tau_age_foi)
  }
  foi_age_mu <- mean(foi_age_effect[1:n_ageclass])

  ############################################################
  ############################################################
  ### Age/Period Survival Model
  ############################################################
  ############################################################

  ####################################
  ### Susceptibles survival intercept
  ####################################

  beta0_sus_temp ~ dnorm(0, .1)
  sus_mix ~ dunif(-1, 1)
  beta0_survival_sus <- beta0_sus_temp * sus_mix

  ##################################
  ### Infected survival intercept
  ##################################

  beta0_inf_temp ~ dnorm(0, .1)
  inf_mix ~ dunif(-1, 1)
  beta0_survival_inf <- beta0_inf_temp * inf_mix

  ########################################
  ### Priors for Age Effects Survival
  ########################################

  ### Age effects
  for (k in 1:nknots_age) {
    ln_b_age_survival[k] ~ dnorm(0, tau_age_survival)
    b_age_survival[k] <- exp(ln_b_age_survival[k])
  }
  tau_age_survival ~ dgamma(1, 1)

  for (t in 1:nT_age) {
    age_effect_survival_temp[t] <- inprod(b_age_survival[1:nknots_age],
                                     Z_age[t, 1:nknots_age])
  }
  mu_age_effect_survival_temp <- mean(age_effect_survival_temp[1:nT_age])
  
  for (t in 1:nT_age) {
    age_effect_survival[t] <-  age_effect_survival_temp[t] -
                               mu_age_effect_survival_temp
  }

  ########################################
  ### Priors for Period Effects Survival
  ########################################

  mix_survival ~ dunif(-1, 1)
  ln_sk_period ~ dnorm(0, sd = 1)
  sdk_period <- exp(mix_survival * ln_sk_period)
  tauk_period <- 1 / sdk_period^2
  stauk_period <- sqrt(tauk_period)
  sda_period ~ T(dnorm(0, sd = 1), 0, Inf)#<- 1/sqrt(taua_period)
  taua_period <- 1 / sda_period^2
  for (i in 1:(nknots_period)) {
    alpha_period[i] ~ dnorm(0, 1)
    alphau_period[i] <- sda_period * alpha_period[i]
  }
  ratioinf_period <- sdk_period / sda_period #ratio of variability

  period_effect_surv[1:nT_period] <- kernel_conv(
    nT = nT_period,
    Z = Z_period[1:nT_period, 1:nknots_period],
    stauk = stauk_period,
    nconst = nconst,
    tauk = tauk_period,
    nknots = nknots_period,
    alphau = alphau_period[1:nknots_period]
  )


  #######################################################################
  #######################################################################
  ## Likelihoods of Joint Model
  #######################################################################
  #######################################################################


  #######################################################################
  ###
  ###   User defined distribution for likelihood for
  ###   infected deer mortalities for radio marked deer that
  ###   enter the study as test positive at capture
  ###
  ###   d_fit_icap_mort
  ###
  ###   Overleaf Equation (12)
  ###
  #######################################################################

    y_icap_mort ~ dIcapMort(
      n_samples = nIcapMort,
      e = icap_mort_e_age[1:nIcapMort],
      r = icap_mort_r_age[1:nIcapMort],
      s = icap_mort_s_age[1:nIcapMort],
      age2date = icap_mort_age2date[1:nIcapMort],
      beta0_inf = beta0_survival_inf,
      age_effect_surv = age_effect_survival[1:nT_age],
      period_effect_surv = period_effect_surv[1:nT_period]
      )


    #######################################################################
    ###
    ###   User defined distribution for likelihood for
    ###   infected deer mortalities for radio marked deer that
    ###   enter the study as test positive at capture
    ###
    ###   d_fit_icap_cens
    ###
    ###   Overleaf Equation (10)
    ###
    #######################################################################

    y_icap_cens ~ dIcapCens(
      n_samples = nIcapCens,
      e = icap_cens_e_age[1:nIcapCens],
      r = icap_cens_r_age[1:nIcapCens],
      age2date = icap_cens_age2date[1:nIcapCens],
      beta0_inf = beta0_survival_inf,
      age_effect_surv = age_effect_survival[1:nT_age],
      period_effect_surv = period_effect_surv[1:nT_period]
      )


    #######################################################################
    ###
    ###   User defined distribution for likelihood for
    ###   infected deer mortalities for radio marked deer that
    ###   enter the study as test negative at capture
    ###
    ###   d_fit_idead
    ###
    ###   Overleaf Equation (24)
    ###
    #######################################################################

    y_idead ~ dNegCapPosMort(
        n_samples = nNegCapPosMort,
        e = idead_e_age[1:nNegCapPosMort],
        r = idead_r_age[1:nNegCapPosMort],
        s = idead_s_age[1:nNegCapPosMort],
        age2date = idead_age2date[1:nNegCapPosMort],
        beta0_sus = beta0_survival_sus,
        beta0_inf = beta0_survival_inf,
        age_effect_surv = age_effect_survival[1:nT_age],
        period_effect_surv = period_effect_surv[1:nT_period],
        foi_age_effect = foi_age_effect[1:n_ageclass],
        age_lookup = age_lookup[1:nT_age]
        )


    #######################################################################
    ###
    ###   User defined distribution for likelihood for
    ###   uninfected radio-marked deer mortalities:
    ###   test neg at cap and tested mort
    ###
    ###   d_fit_sus_mort_posttest
    ###
    ###   Overleaf Equation (11)
    ###
    #######################################################################

    y_sus_mort_posttest ~ dSusMortTest(
        n_samples = nSusMortTest,
        e = sus_mort_posttest_e_age[1:nSusMortTest],
        r = sus_mort_posttest_r_age[1:nSusMortTest],
        s = sus_mort_posttest_s_age[1:nSusMortTest],
        age2date = sus_mort_posttest_age2date[1:nSusMortTest],
        beta0_sus = beta0_survival_sus,
        age_effect_surv = age_effect_survival[1:nT_age],
        period_effect_surv = period_effect_surv[1:nT_period],
        foi_age_effect =  foi_age_effect[1:n_ageclass],
        age_lookup = age_lookup[1:nT_age]
        )


    #######################################################################
    ###
    ###   User defined distribution for likelihood for
    ###   Uninfected radio-marked deer right censored:
    ###   Test neg at cap and censoring
    ###
    ###   d_fit_sus_cens_postno
    ###
    ###   Overleaf Equation (9)
    ###
    #######################################################################

    y_sus_cens_postno ~ dSusCensNo(
          n_samples = nSusCensNo,
          e = sus_cens_postno_e_age[1:nSusCensNo],
          r = sus_cens_postno_r_age[1:nSusCensNo],
          age2date = sus_cens_postno_age2date[1:nSusCensNo],
          beta0_sus = beta0_survival_sus,
          beta0_inf = beta0_survival_inf,
          age_effect_surv = age_effect_survival[1:nT_age],
          period_effect_surv = period_effect_surv[1:nT_period],
          foi_age_effect = foi_age_effect[1:n_ageclass],
          age_lookup = age_lookup[1:nT_age]
          )

    #######################################################################
    ###
    ###   Derived parameters: 
    ###   Annual survival estimates "May 15 - May 14"
    ###   (not using hazards from Jan - May of the first year)
    ###
    ###
    #######################################################################

    sn_sus[1:n_age, 1:n_year]  <- calc_surv_aah(nT_age = nT_age,
        nT_period = nT_period,
        nT_age_short = nT_age_short,
        nT_age_surv_aah = nT_age_surv_aah,
        beta0 = beta0_survival_sus,
        age_effect = age_effect_survival[1:nT_age],
        period_effect = period_effect_surv[1:nT_period],
        yr_start_age = yr_start_age[1:n_yr_start_age],
        yr_start_pop = yr_start_pop[1:n_year],
        n_year = n_year,
        n_age = n_age)

    sn_inf[1:n_age, 1:n_year]  <- calc_surv_aah(nT_age = nT_age,
        nT_period = nT_period,
        nT_age_short = nT_age_short,
        nT_age_surv_aah = nT_age_surv_aah,
        beta0 = beta0_survival_inf,
        age_effect = age_effect_survival[1:nT_age],
        period_effect = period_effect_surv[1:nT_period], 
        yr_start_age = yr_start_age[1:n_yr_start_age],
        yr_start_pop = yr_start_pop[1:n_year],
        n_year = n_year,
        n_age = n_age)

    psi[1:n_age, 1:n_year] <- calc_infect_prob(age_lookup = age_lookup[1:nT_age],
                        n_age = n_age,
                        yr_start = yr_start_age[1:n_yr_start_age],
                        age_foi = foi_age_effect[1:n_ageclass],
                        nT_period = nT_period,
                        n_year = n_year,
                        nT_age_surv_aah = nT_age_surv_aah
                        )

})#end model statement
