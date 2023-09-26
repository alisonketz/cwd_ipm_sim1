

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
  
  beta_male ~ dnorm(0, .01)

  ##############################
  ### Force of infection model
  ##############################
  tau_age_foi_male  ~ dgamma(1, 1)
  tau_age_foi_female  ~ dgamma(1, 1)
  tau1_age_foi_male <- .0000001 * tau_age_foi_male
  tau1_age_foi_female <- .0000001 * tau_age_foi_female
  m_age_foi[1] ~ dnorm(0, tau1_age_foi_male)
  f_age_foi[1] ~ dnorm(0, tau1_age_foi_female)
  m_age_foi[2] ~ dnorm(0, tau1_age_foi_male)
  f_age_foi[2] ~ dnorm(0, tau1_age_foi_female)
  for (i in 3:n_ageclassm) {
    m_age_foi[i]~dnorm(2 * m_age_foi[i-1] - m_age_foi[i-2], tau_age_foi_male)
  }
  for (i in 3:n_ageclassf) {
    f_age_foi[i]~dnorm(2 * f_age_foi[i-1] - f_age_foi[i-2], tau_age_foi_female)
  }
  m_age_foi_mu <- mean(m_age_foi[1:n_ageclassm])
  f_age_foi_mu <- mean(f_age_foi[1:n_ageclassf])

  # Period effects
  tau_period_foi_male  ~ dgamma(1, 1)
  tau_period_foi_female  ~ dgamma(1, 1)

  f_period_foi[1:n_year] ~ dcar_normal(adj = adj_period[1:n_adj_period],
                                  weights = weights_period[1:n_adj_period],
                                  num = num_period[1:n_year],
                                  tau = tau_period_foi_female,
                                  zero_mean = 1)
  
  m_period_foi[1:n_year] ~ dcar_normal(adj = adj_period[1:n_adj_period],
                                  weights = weights_period[1:n_adj_period],
                                  num = num_period[1:n_year],
                                  tau = tau_period_foi_male,
                                  zero_mean = 1)

  ### random effect for East/West spatial model
  space[1] <- 0
  space[2]  <- space_temp * space_alpha

  space_temp ~ dnorm(0, 1)
  space_alpha ~ dunif(-1, 1)

  ############################################################
  ############################################################
  ### Age/Period Survival Model
  ############################################################
  ############################################################

  ##############################
  ### Susceptibles
  ##############################

  #Priors for intercept and covariate
  beta0_sus_temp ~ dnorm(0, .01)
  sus_mix ~ dunif(-1, 1)
  beta0_survival_sus <- beta0_sus_temp * sus_mix

  #Priors for Age and Period effects
  #Age effects
  for (k in 1:nknots_age) {
    b_age_survival[k] ~ dnorm(0, tau_age_survival)
  }
  tau_age_survival ~ dgamma(1, 1)

  for (t in 1:nT_age_surv) {
    age_effect_survival[t] <- inprod(b_age_survival[1:nknots_age],
                                     Z_age[t, 1:nknots_age])
  }

  #Period effects from collar data
  for (k in 1:nknots_period) {
    b_period_survival[k] ~ dnorm(0, tau_period_survival)
  }
  tau_period_survival ~ dgamma(1, 1)
  for (t in 1:nT_period_collar) {
    period_effect_surv[t] <- inprod(b_period_survival[1:nknots_period],
                                    Z_period[t, 1:nknots_period])
  }

  #Period effects from aah data
  tau_period_precollar ~ dgamma(1,1)
  for (k in 1:(n_year_precollar + 1)) {
    period_annual_survival[k] ~ dnorm(0, tau_period_precollar)
  }

  period_effect_survival[1:nT_period_overall] <- set_period_effects_constant(
        n_year_precollar = n_year_precollar,
        nT_period_precollar = nT_period_precollar,
        nT_period_collar = nT_period_collar,
        nT_period_overall = nT_period_overall,
        yr_start = yr_start[1:n_year],
        yr_end = yr_end[1:n_year],
        period_effect_surv = period_effect_surv[1:nT_period_collar],
        period_annual_survival = period_annual_survival[1:(n_year_precollar + 1)]
  )
  ##################################
  ## Infected survival intercept
  ##################################

  #Priors for intercept and covariate
  beta0_inf_temp ~ dnorm(0, .01)
  inf_mix ~ dunif(-1, 1)
  beta0_survival_inf <- beta0_inf_temp * inf_mix

  #######################################################################
  #######################################################################
  ## Likelihoods of Joint Model
  #######################################################################
  #######################################################################

  #######################################################################
  ###
  ###   User defined distribution for likelihood for
  ###   infected harvest deer
  ###
  ###   d_fit_hunt_pos
  ###   Overleaf Equation 3
  ###
  #######################################################################

  y_hunt_pos ~ dInfHarvest(
                  n_cases = hunt_pos_n_cases[1:nInfHarvest],
                  n_samples = nInfHarvest,
				          a = hunt_pos_ageweeks[1:nInfHarvest], #age (weeks) at harvest
                  sex = hunt_pos_sex[1:nInfHarvest],
                  age2date = hunt_pos_age2date[1:nInfHarvest],
                  beta_male = beta_male,
                  beta0_sus = beta0_survival_sus,
                  beta0_inf = beta0_survival_inf,
                  age_effect_surv = age_effect_survival[1:nT_age_surv],
                  period_effect_surv = period_effect_survival[1:nT_period_overall],
                  period_lookup_surv = lookup_pe_surv[1:nT_period_overall_ext],
                  f_age_foi = f_age_foi[1:n_ageclassf],
                  m_age_foi = m_age_foi[1:n_ageclassm],
                  age_lookup_f = age_lookup_f[1:n_age_lookup_f],
                  age_lookup_m = age_lookup_m[1:n_age_lookup_m],
                  period_lookup = period_lookup[1:n_period_lookup],
                  f_period_foi = f_period_foi[1:n_year],
                  m_period_foi = m_period_foi[1:n_year],
                  space = space[1:n_study_area],
                  sect = sect_hunt_pos[1:nInfHarvest]
                  )

#######################################################################
###
###   User defined distribution for likelihood for
###   uninfected harvest deer
###   d_fit_hunt_neg
###   Overleaf Equation (5)
###
#######################################################################

  y_hunt_neg ~ dSusHarvest(
                  n_cases = hunt_neg_n_cases[1:nInfHarvest],
                  n_samples = nSusHarvest,
				          a = hunt_neg_ageweeks[1:nSusHarvest], #age (weeks) at harvest
                  sex = hunt_neg_sex[1:nSusHarvest],
                  age2date = hunt_neg_age2date[1:nSusHarvest],
                  beta_male = beta_male,
                  beta0_sus = beta0_survival_sus,
                  age_effect_surv = age_effect_survival[1:nT_age_surv],
                  period_effect_surv = period_effect_survival[1:nT_period_overall],
                  period_lookup_surv = lookup_pe_surv[1:nT_period_overall_ext],
                  f_age_foi = f_age_foi[1:n_ageclassf],
                  m_age_foi = m_age_foi[1:n_ageclassm],
                  age_lookup_f = age_lookup_f[1:n_age_lookup_f],
                  age_lookup_m = age_lookup_m[1:n_age_lookup_m],
                  period_lookup = period_lookup[1:n_period_lookup],
                  f_period_foi = f_period_foi[1:n_year],
                  m_period_foi = m_period_foi[1:n_year],
                  space = space[1:n_study_area],
                  sect = sect_hunt_neg[1:nSusHarvest]
                  )

#######################################################################
###
###   User defined distribution for likelihood for
###   Uninfected radio-marked deer right censor:
###   Test neg at cap and censoring
###
###   d_fit_sus_cens_posttest
###   Overleaf Equation 7
###
#######################################################################


    y_sus_cens_posttest ~ dSusCensTest(
        n_samples = nSusCensTest,
        e = sus_cens_posttest_left_age_e[1:nSusCensTest],
        r = sus_cens_posttest_right_age_r[1:nSusCensTest],
        sex = sus_cens_posttest_sex[1:nSusCensTest],
        age2date = sus_cens_posttest_age2date[1:nSusCensTest],
        beta_male = beta_male,
        beta0_sus = beta0_survival_sus,
        age_effect_surv = age_effect_survival[1:nT_age_surv],
        period_effect_surv = period_effect_survival[1:nT_period_overall],
        f_age_foi = f_age_foi[1:n_ageclassf],
        m_age_foi = m_age_foi[1:n_ageclassm],
        age_lookup_f = age_lookup_f[1:n_age_lookup_f],
        age_lookup_m = age_lookup_m[1:n_age_lookup_m],
        period_lookup = period_lookup[1:n_period_lookup],
        f_period_foi = f_period_foi[1:n_year],
        m_period_foi = m_period_foi[1:n_year],
        space = space[1:n_study_area],
        sect_sus_cens_posttest[1:nSusCensTest]
        )



#######################################################################
###
###   User defined distribution for likelihood for
###   Uninfected radio-marked deer right censored:
###   Test neg at cap and censoring
###
###   d_fit_sus_cens_postno
###   d_fit_endlive
###
###   Overleaf Equation (9)
###
#######################################################################

  y_sus_cens_postno ~ dSusCensNo(
        n_samples = nSusCensNo,
        e = sus_cens_postno_left_age_e[1:nSusCensNo],
        r = sus_cens_postno_right_age_r[1:nSusCensNo],
        sex = sus_cens_postno_sex[1:nSusCensNo],
        age2date = sus_cens_postno_age2date[1:nSusCensNo],
        beta_male = beta_male,
        beta0_sus = beta0_survival_sus,
        beta0_inf = beta0_survival_inf,
        age_effect_surv = age_effect_survival[1:nT_age_surv],
        period_effect_surv = period_effect_survival[1:nT_period_overall],
        f_age_foi = f_age_foi[1:n_ageclassf],
        m_age_foi = m_age_foi[1:n_ageclassm],
        age_lookup_f = age_lookup_f[1:n_age_lookup_f],
        age_lookup_m = age_lookup_m[1:n_age_lookup_m],
        period_lookup = period_lookup[1:n_period_lookup],
        f_period_foi = f_period_foi[1:n_year],
        m_period_foi = m_period_foi[1:n_year],
        space = space[1:n_study_area],
        sect = sect_sus_cens_postno[1:nSusCensNo]
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
      e = sus_mort_posttest_left_age_e[1:nSusMortTest],
      r = sus_mort_posttest_right_age_r[1:nSusMortTest],
      s = sus_mort_posttest_right_age_s[1:nSusMortTest],
      sex = sus_mort_posttest_sex[1:nSusMortTest],
      age2date = sus_mort_posttest_age2date[1:nSusMortTest],
      beta_male = beta_male,
      beta0_sus = beta0_survival_sus,
      age_effect_surv = age_effect_survival[1:nT_age_surv],
      period_effect_surv = period_effect_survival[1:nT_period_overall],
      f_age_foi = f_age_foi[1:n_ageclassf],
      m_age_foi = m_age_foi[1:n_ageclassm],
      age_lookup_f = age_lookup_f[1:n_age_lookup_f],
      age_lookup_m = age_lookup_m[1:n_age_lookup_m],
      period_lookup = period_lookup[1:n_period_lookup],
      f_period_foi = f_period_foi[1:n_year],
      m_period_foi = m_period_foi[1:n_year],
      space = space[1:n_study_area],
      sect = sect_sus_mort_posttest[1:nSusMortTest]
      )

#######################################################################
###
###   User defined distribution for likelihood for
###   uninfected radio-marked deer mortalities:
###   test neg at cap and no test at mortality
###
###   d_fit_sus_mort_postno
###
###   Overleaf Equation (13)
###
#######################################################################


  y_sus_mort_postno ~ dSusMortNoTest(
      n_samples = nSusMortNoTest,
      e = sus_mort_postno_left_age_e[1:nSusMortNoTest],
      r = sus_mort_postno_right_age_r[1:nSusMortNoTest],
      s = sus_mort_postno_right_age_s[1:nSusMortNoTest],
      dn1 = sus_mort_postno_dn1[1:nSusMortNoTest],
      sex = sus_mort_postno_sex[1:nSusMortNoTest],
      age2date = sus_mort_postno_age2date[1:nSusMortNoTest],
      beta_male = beta_male,
      beta0_sus = beta0_survival_sus,
      beta0_inf = beta0_survival_inf,
      age_effect_surv = age_effect_survival[1:nT_age_surv],
      period_effect_surv = period_effect_survival[1:nT_period_overall],
      f_age_foi = f_age_foi[1:n_ageclassf],
      m_age_foi = m_age_foi[1:n_ageclassm],
      age_lookup_f = age_lookup_f[1:n_age_lookup_f],
      age_lookup_m = age_lookup_m[1:n_age_lookup_m],
      period_lookup = period_lookup[1:n_period_lookup],
      f_period_foi = f_period_foi[1:n_year],
      m_period_foi = m_period_foi[1:n_year],
      space = space[1:n_study_area],
      sect = sect_sus_mort_postno[1:nSusMortNoTest]
      )

#######################################################################
###
###   User defined distribution for likelihood for
###   infected deer mortalities for radio marked deer that
###   enter the study as test positive at capture
###
###   d_fit_icap_cens
###
###   Overleaf Equation (15)
###
#######################################################################

  y_icap_cens ~ dIcapCens(
    n_samples = nIcapCens,
    e = icap_cens_left_age_e[1:nIcapCens],
    r = icap_cens_right_age_r[1:nIcapCens],
    sex = icap_cens_sex[1:nIcapCens],
    age2date = icap_cens_age2date[1:nIcapCens],
    beta_male = beta_male,
    beta0_sus = beta0_survival_sus,
    beta0_inf = beta0_survival_inf,
    age_effect_surv = age_effect_survival[1:nT_age_surv],
    period_effect_surv = period_effect_survival[1:nT_period_overall],
    f_age_foi = f_age_foi[1:n_ageclassf],
    m_age_foi = m_age_foi[1:n_ageclassm],
    age_lookup_f = age_lookup_f[1:n_age_lookup_f],
    age_lookup_m = age_lookup_m[1:n_age_lookup_m],
    period_lookup = period_lookup[1:n_period_lookup],
    f_period_foi = f_period_foi[1:n_year],
    m_period_foi = m_period_foi[1:n_year],
    space = space[1:n_study_area],
    sect_icap_cens[1:nIcapCens]
    )

#######################################################################
###
###   User defined distribution for likelihood for
###   infected deer mortalities for radio marked deer that
###   enter the study as test positive at capture
###
###   d_fit_icap_mort
###
###   Overleaf Equation (17)
###
#######################################################################

  y_icap_mort ~ dIcapMort(
    n_samples = nIcapMort,
    e = icap_mort_left_age_e[1:nIcapMort],
    r = icap_mort_right_age_r[1:nIcapMort],
    s = icap_mort_right_age_s[1:nIcapMort],
    sex = icap_mort_sex[1:nIcapMort],
    age2date = icap_mort_age2date[1:nIcapMort],
    beta_male = beta_male,
    beta0_sus = beta0_survival_sus,
    beta0_inf = beta0_survival_inf,
    age_effect_surv = age_effect_survival[1:nT_age_surv],
    period_effect_surv = period_effect_survival[1:nT_period_overall],
    f_age_foi = f_age_foi[1:n_ageclassf],
    m_age_foi = m_age_foi[1:n_ageclassm],
    age_lookup_f = age_lookup_f[1:n_age_lookup_f],
    age_lookup_m = age_lookup_m[1:n_age_lookup_m],
    period_lookup = period_lookup[1:n_period_lookup],
    f_period_foi = f_period_foi[1:n_year],
    m_period_foi = m_period_foi[1:n_year],
    space = space[1:n_study_area],
    sect = sect_icap_mort[1:nIcapMort]
    )

#######################################################################
###
###   User defined distribution for likelihood for
###   uninfected deer that were test neg at capture,
###   then test negative at recap, that are right censored, 
###   and have been tested post censoring
###
###   d_fit_rec_neg_cens_posttest
###
###   Overleaf Equation (19)
###
#######################################################################

    y_rec_neg_cens_posttest ~ dRecNegCensTest(
      n_samples = nRecNegCensTest,
      e = rec_neg_cens_posttest_left_age_e[1:nRecNegCensTest],
      r = rec_neg_cens_posttest_right_age_r[1:nRecNegCensTest],
      sex = rec_neg_cens_posttest_sex[1:nRecNegCensTest],
      age2date = rec_neg_cens_posttest_age2date[1:nRecNegCensTest],
      beta_male = beta_male,
      beta0_sus = beta0_survival_sus,
      age_effect_surv = age_effect_survival[1:nT_age_surv],
      period_effect_surv = period_effect_survival[1:nT_period_overall],
      f_age_foi = f_age_foi[1:n_ageclassf],
      m_age_foi = m_age_foi[1:n_ageclassm],
      age_lookup_f = age_lookup_f[1:n_age_lookup_f],
      age_lookup_m = age_lookup_m[1:n_age_lookup_m],
      period_lookup = period_lookup[1:n_period_lookup],
      f_period_foi = f_period_foi[1:n_year],
      m_period_foi = m_period_foi[1:n_year],
      space = space[1:n_study_area],
      sect = sect_rec_neg_cens_posttest[1:nRecNegCensTest]
      )

#######################################################################
###
###   User defined distribution for likelihood for
###   uninfected deer that were test neg at capture,
###   then test negative at recap,
###   that die
###
###   d_fit_rec_neg_mort
###
###   Overleaf Equation (23)
###
#######################################################################


    y_rec_neg_mort ~ dRecNegMort(
      n_samples = nRecNegMort,
      e = rec_neg_mort_left_age_e[1:nRecNegMort],
      r = rec_neg_mort_right_age_r[1:nRecNegMort],
      s = rec_neg_mort_right_age_s[1:nRecNegMort],
      sex = rec_neg_mort_sex[1:nRecNegMort],
      age2date = rec_neg_mort_age2date[1:nRecNegMort],
      beta_male = beta_male,
      beta0_sus = beta0_survival_sus,
      age_effect_surv = age_effect_survival[1:nT_age_surv],
      period_effect_surv = period_effect_survival[1:nT_period_overall],
      f_age_foi = f_age_foi[1:n_ageclassf],
      m_age_foi = m_age_foi[1:n_ageclassm],
      age_lookup_f = age_lookup_f[1:n_age_lookup_f],
      age_lookup_m = age_lookup_m[1:n_age_lookup_m],
      period_lookup = period_lookup[1:n_period_lookup],
      f_period_foi = f_period_foi[1:n_year],
      m_period_foi = m_period_foi[1:n_year],
      space = space[1:n_study_area],
      sect = sect_rec_neg_mort[1:nRecNegMort]
      )

#######################################################################
###
###   User defined distribution for likelihood for
###   deer that were test neg at capture,
###   then test positive at recap,
###   than die
###
###   d_fit_rec_pos_mort
###
###   Overleaf Equation (25)
###
#######################################################################

  y_rec_pos_mort ~ dRecPosMort(
      n_samples = nRecPosMort,
      e = rec_pos_mort_left_age_e[1:nRecPosMort],
      r = rec_pos_mort_right_age_r[1:nRecPosMort],
      s = rec_pos_mort_right_age_s[1:nRecPosMort],
      dn1 = rec_pos_mort_dn1[1:nRecPosMort],
      dn = rec_pos_mort_dn[1:nRecPosMort],
      sex = rec_pos_mort_sex[1:nRecPosMort],
      age2date = rec_pos_mort_age2date[1:nRecPosMort],
      beta_male = beta_male,
      beta0_sus = beta0_survival_sus,
      beta0_inf = beta0_survival_inf,
      age_effect_surv = age_effect_survival[1:nT_age_surv],
      period_effect_surv = period_effect_survival[1:nT_period_overall],
      f_age_foi = f_age_foi[1:n_ageclassf],
      m_age_foi = m_age_foi[1:n_ageclassm],
      age_lookup_f = age_lookup_f[1:n_age_lookup_f],
      age_lookup_m = age_lookup_m[1:n_age_lookup_m],
      period_lookup = period_lookup[1:n_period_lookup],
      f_period_foi = f_period_foi[1:n_year],
      m_period_foi = m_period_foi[1:n_year],
      space = space[1:n_study_area],
      sect = sect_rec_pos_mort[1:nRecPosMort]
      )

#######################################################################
###
###   User defined distribution for likelihood for
###   infected deer that were test neg at capture,
###   then test positive at recap,
###   that are right censored
###
###   d_fit_rec_pos_cens
###
###   Overleaf Equation (27)
###
#######################################################################

  y_rec_pos_cens ~ dRecPosCens(
      e = rec_pos_cens_left_age_e,
      r = rec_pos_cens_right_age_r,
      dn1 = rec_pos_cens_dn1,
      dn = rec_pos_cens_dn,
      sex = rec_pos_cens_sex,
      age2date = rec_pos_cens_age2date,
      beta_male = beta_male,
      beta0_sus = beta0_survival_sus,
      beta0_inf = beta0_survival_inf,
      age_effect_surv = age_effect_survival[1:nT_age_surv],
      period_effect_surv = period_effect_survival[1:nT_period_overall],
      f_age_foi = f_age_foi[1:n_ageclassf],
      m_age_foi = m_age_foi[1:n_ageclassm],
      age_lookup_f = age_lookup_f[1:n_age_lookup_f],
      age_lookup_m = age_lookup_m[1:n_age_lookup_m],
      period_lookup = period_lookup[1:n_period_lookup],
      f_period_foi = f_period_foi[1:n_year],
      m_period_foi = m_period_foi[1:n_year],
      space = space[sect_rec_pos_cens]
  )

#######################################################################
###
###   User defined distribution for likelihood for
###   infected deer mortalities for radio marked deer that
###   enter the study as test negative at capture
###
###   d_fit_idead
###
###   Overleaf Equation (29)
###
#######################################################################


 y_idead ~ dNegCapPosMort(
      n_samples = nNegCapPosMort,
      e = idead_left_age_e[1:nNegCapPosMort],
      r = idead_right_age_r[1:nNegCapPosMort],
      s = idead_right_age_s[1:nNegCapPosMort],
      dn1 = idead_dn1[1:nNegCapPosMort],
      dn = idead_dn[1:nNegCapPosMort],
      sex = idead_sex[1:nNegCapPosMort],
      age2date = idead_age2date[1:nNegCapPosMort],
      beta_male = beta_male,
      beta0_sus = beta0_survival_sus,
      beta0_inf = beta0_survival_inf,
      age_effect_surv = age_effect_survival[1:nT_age_surv],
      period_effect_surv = period_effect_survival[1:nT_period_overall],
      f_age_foi = f_age_foi[1:n_ageclassf],
      m_age_foi = m_age_foi[1:n_ageclassm],
      age_lookup_f = age_lookup_f[1:n_age_lookup_f],
      age_lookup_m = age_lookup_m[1:n_age_lookup_m],
      period_lookup = period_lookup[1:n_period_lookup],
      f_period_foi = f_period_foi[1:n_year],
      m_period_foi = m_period_foi[1:n_year],
      space = space[1:n_study_area],
      sect = sect_idead[1:nNegCapPosMort]
      )

#######################################################################
###
###   User defined distribution for likelihood for
###   age-at-harvest deer used to estimate period effects
###
###   d_aah
###
###   Overleaf Equation (31)
###
#######################################################################


  y_aah ~ dAAH(
      n_samples = nAAH,
      a = aah_ageweeks[1:nAAH],
      sex = aah_sex[1:nAAH],
      age2date = aah_age2date[1:nAAH],
      n_cases = aah_n[1:nAAH],
      beta_male = beta_male,
      beta0_sus = beta0_survival_sus,
      beta0_inf = beta0_survival_inf,
      age_effect_surv = age_effect_survival[1:nT_age_surv],
      period_effect_surv = period_effect_survival[1:nT_period_overall],
      f_age_foi = f_age_foi[1:n_ageclassf],
      m_age_foi = m_age_foi[1:n_ageclassm],
      age_lookup_f = age_lookup_f[1:n_age_lookup_f],
      age_lookup_m = age_lookup_m[1:n_age_lookup_m],
      period_lookup = period_lookup[1:n_period_lookup],
      f_period_foi = f_period_foi[1:n_year],
      m_period_foi = m_period_foi[1:n_year],
      space = space[1:n_study_area],
      sect = sect_aah[1:nAAH]
      )

  #######################################################
  #######################################################
  #######################################################
  ### Cause-specific mortality model
  #######################################################
  #######################################################
  #######################################################
  #priors
  #sex-specific hunt probability
  beta0_cause ~ dnorm(0, .01)
  beta_cause_gun ~ dnorm(0, .01)
  beta_cause_ng ~ dnorm(0, .01)
  beta_cause_male ~ dnorm(0, .01)

  #######################################################
  ###
  ### Likelihood for cause of death
  ###
  #######################################################

  for (i in 1:records_cause) {
      p_cause[i]  <- ilogit(beta0_cause +
                            Z_cause_ng[interval[i]] * beta_cause_ng +
                            Z_cause_gun[interval[i]] * beta_cause_gun +
                            sex_cause[i] * beta_cause_male)
      mort_hh[i] ~ dbin(size = 1, prob = p_cause[i])
  }

  #######################################################
  ###
  ### Derived parameter for cause of death by hunter harvest
  ### given 
  ###
  #######################################################

  p_nogun_f <- ilogit(beta0_cause +
                      beta_cause_ng)
  p_gun_f <- ilogit(beta0_cause +
                    beta_cause_ng +
                    beta_cause_gun)
  p_nogun_m <- ilogit(beta0_cause +
                      beta_cause_ng +
                      beta_cause_male)
  p_gun_m <- ilogit(beta0_cause +
                    beta_cause_ng +
                    beta_cause_gun +
                    beta_cause_male)

  ##############################################################################
  ##############################################################################
  ### Age-at-harvest population model
  ##############################################################################
  ##############################################################################


  #######################################
  #### Initial population, currently 
  ### based on empirical bayes approach 
  ### of using the first year of data
  ########################################

  #should this be different for pos/neg or m/f for study area?
  for(k in 1:n_study_area){
    for(i in 1:2){
      tau_pop[k,i] ~ dgamma(1,1)
    }
  }
  
  for(k in 1:n_study_area){
    for (a in 1:n_agef) {

      #Initial population structure pop[sex,age,year] for susceptible deer
      llpop_sus[k, 1, a, 1] ~ dnorm(f_logpop_sus[k, a], tau_pop[k, 1])
      pop_sus[k, 1, a, 1] <- exp(llpop_sus[k, 1, a, 1])

      #Initial population structure pop[study_area=k,sex=i,year=t,age=a]
      llpop_inf[k, 1, a, 1] ~ dnorm(f_logpop_inf[k, a], tau_pop[k, 2])
      pop_inf[k, 1, a, 1] <- exp(llpop_inf[k, 1, a, 1])
    }

    for (a in 1:n_agem) {
        ### East
        #Initial population structure pop[year=1,sex=i,age=a] for susceptible deer
        llpop_sus[k, 2, a, 1] ~ dnorm(m_logpop_sus[k, a], tau_pop[k, 1])
        pop_sus[k, 2, a, 1] <- exp(llpop_sus[k, 2, a, 1])

        #Initial population structure pop for infected deer
        llpop_inf[k, 2, a, 1] ~ dnorm(m_logpop_inf[k, a], tau_pop[k, 2])
        pop_inf[k, 2, a, 1] <- exp(llpop_inf[k, 2, a, 1])
    }
  }

  ############################
  ####Reporting Rates
  ############################

  report_overall ~ dbeta(report_hyp_all[1], report_hyp_all[2])
  for(t in 1:22){#1992-2014
    report[t]  <- report_overall
  }
  for(t in 23:27){ #2015-2020
    report[t] ~ dbeta(report_hyp_y[t - 22, 1], report_hyp_y[t - 22, 2])
  }
  report[28]  <- report_overall #2021

  ############################
  #### Fecundity
  ############################

  mu_fec ~ dnorm(0, 1)
  fec_prec_eps ~ dgamma(1, 1)

  #Observations of fawns & does overall from the 3 counties
  for(t in 1:n_year_fec_early) {
    fec_epsilon[t] ~ dnorm(0, fec_prec_eps)
    fec[t] <- exp(mu_fec + fec_epsilon[t])
    Nfawn[t] ~ dpois(fec[t] * Ndoe[t])
  }

  #for 2017:2021
  for(t in (n_year_fec_early + 1):n_year){
    fec[t] ~ dgamma(obs_ct_fd_alpha[t], obs_ct_fd_beta[t])
  }

  ###################################################
  #### Overall Survival Susceptibles
  ###################################################

  sn_sus[1:n_sex, 1:n_agef, 1:n_year] <- calc_surv_aah(
          nT_age = nT_age_surv,
          nT_period = nT_period_overall,
          beta0 = beta0_survival_sus,
          beta_male = beta_male,
          age_effect = age_effect_survival[1:nT_age_surv_aah],
          period_effect = period_effect_survival[1:nT_period_overall],
          yr_start = yr_start[1:n_year],
          yr_end = yr_end[1:n_year],
          intvl_step_yr = intvl_step_yr,
          n_year = n_year,
          n_agef = n_agef,
          n_agem = n_agem)

  ###################################################
  #### Overall Survival CWD Infected
  ###################################################

  sn_inf[1:n_sex,1:n_agef,1:n_year] <- calc_surv_aah(
          nT_age = nT_age_surv,
          nT_period = nT_period_overall,
          beta0 = beta0_survival_inf,
          beta_male = beta_male,
          age_effect = age_effect_survival[1:nT_age_surv],
          period_effect = period_effect_survival[1:nT_period_overall],
          yr_start = yr_start[1:n_year],
          yr_end = yr_end[1:n_year],
          intvl_step_yr = intvl_step_yr,
          n_year = n_year,
          n_agef = n_agef,
          n_agem = n_agem)

  ###################################################
  #### Hunting Survival Susceptibles
  ###################################################

  sh_sus[1:n_sex,1:n_agef,1:n_year] <- calc_surv_harvest(
          nT_age = nT_age_surv, 
          nT_period = nT_period_overall,
          beta0 = beta0_survival_sus,
          beta_male = beta_male,
          age_effect = age_effect_survival[1:nT_age_surv],
          period_effect = period_effect_survival[1:nT_period_overall],
          intvl_step_yr = intvl_step_yr,
          n_year = n_year,
          n_agef = n_agef,
          n_agem = n_agem,
          ng_start = ng_start[1:n_year],
          gun_start = gun_start[1:n_year],
          gun_end = gun_end[1:n_year],
          ng_end = ng_end[1:n_year],
          yr_start = yr_start[1:n_year],
          yr_end = yr_end[1:n_year],
          p_nogun_f = p_nogun_f,
          p_nogun_m = p_nogun_m,
          p_gun_f = p_gun_f,
          p_gun_m = p_gun_m
          )

  ###################################################
  #### Hunting Survival Infected
  ###################################################

  sh_inf[1:n_sex,1:n_agef,1:n_year] <- calc_surv_harvest(
          nT_age = nT_age_surv,
          nT_period = nT_period_overall,
          beta0 = beta0_survival_inf,
          beta_male = beta_male,
          age_effect = age_effect_survival[1:nT_age_surv],
          period_effect = period_effect_survival[1:nT_period_overall],
          intvl_step_yr = intvl_step_yr,
          n_year = n_year,
          n_agef = n_agef,
          n_agem = n_agem,
          ng_start = ng_start[1:n_year],
          gun_start = gun_start[1:n_year],
          gun_end = gun_end[1:n_year],
          ng_end = ng_end[1:n_year],
          yr_start = yr_start[1:n_year],
          yr_end = yr_end[1:n_year],
          p_nogun_f = p_nogun_f,
          p_nogun_m = p_nogun_m,
          p_gun_f = p_gun_f,
          p_gun_m = p_gun_m
          )


  ###################################################
  #### Probability of Infection based on FOI hazards
  ###################################################

  psi[1:n_study_area, 1:n_sex, 1:n_agef, 1:n_year] <- 
      calc_infect_prob(
            age_lookup_f = age_lookup_f_conv[1:n_age_lookup_f_conv],
            age_lookup_m = age_lookup_m_conv[1:n_age_lookup_m_conv],
            Nage_lookup_f = n_age_lookup_f_conv,
            Nage_lookup_m = n_age_lookup_m_conv,
            n_agef = n_agef,
            n_agem = n_agem,
            yr_end = yr_end[1:n_year],
            f_age = f_age_foi[1:n_ageclassf],
            m_age = m_age_foi[1:n_ageclassm],
            f_period = f_period_foi[1:n_year],
            m_period = m_period_foi[1:n_year],
            n_year = n_year,
            n_sex = n_sex,
            n_study_area = n_study_area, 
            space = space[n_study_area]
            )


  ###################################################
  #### Earn-a-buck correction factor
  #### based on Van Deelen et al (2010)
  ###################################################

  # eab_antlerless_temp ~ dgamma(eab_antlerless_alpha,eab_antlerless_beta)
  # eab_antlered_temp ~ dgamma(eab_antlered_alpha,eab_antlered_beta)
  # for(t in 1:n_year) {
  #   eab_antlerless[t] <- eab_antlerless_temp^x_eab[t]
  #   eab_antlered[t]  <- eab_antlered_temp^x_eab[t]
  # }

  ######################################################################
  ###
  ### Population Process Model
  ### Population Projection
  ### pop_proj temporarily holds the projected age class
  ###
  ######################################################################

    for (k in 1:n_study_area){

      ##################
      ### Susceptible
      ##################

      for (t in 2:n_year) {
        ###########
        # Females
        ###########
        #Female: project forward anually
        for (a in 1:(n_agef - 1)) {
          pop_sus_proj[k, 1, a, t] <- pop_sus[k, 1, a, t - 1] * sn_sus[1, a, t - 1] * (1 - psi[k, 1, a, t - 1])
        }

        #female max age class
        pop_sus_proj[k, 1, n_agef, t] <- pop_sus_proj[k, 1,(n_agef - 1), t] +
                                      pop_sus[k, 1,  n_agef, t - 1] * sn_sus[1, n_agef,t - 1] * (1 - psi[k, 1, n_agef, t - 1])

        #Female: set projection into population model matrix
        for (a in 2:n_agef) {
          pop_sus[k, 1, a, t] <- pop_sus_proj[k, 1, (a - 1), t]
        }

        #Male: fawn class = total #females * unisex fawns per female/2
        #(should this be divided by 2?)
        pop_sus[k, 1, 1, t] <- (sum(pop_sus_proj[k, 1, 1:n_agef, t]) +
                            sum(pop_inf_proj[k, 1, 1:n_agef, t])) * fec[t] * (1 - psi[k, 2, 1, t]) / 2 
        ###########
        # Males
        ###########

        #Male: project forward anually
        for (a in 1:(n_agem - 1)) {
          pop_sus_proj[k, 2, a, t] <- pop_sus[k, 2, a, t - 1] * sn_sus[2, a, t - 1] * (1 - psi[k, 2, a, t - 1])
        }
        
        #Male: max age class
        pop_sus_proj[k, 2, n_agem, t] <- pop_sus_proj[k, 2, (n_agem - 1), t] +
                                      pop_sus[k, 2, n_agem, t - 1] * sn_sus[2, n_agem, t - 1] * (1 - psi[k, 2, n_agem, t - 1])


        #Male: set projection into population model matrix
        for (a in 2:n_agem) {
          pop_sus[k, 2, a, t] <- pop_sus_proj[k, 2, (a - 1), t]
        }

        # Male: fawn class = total #females * unisex fawns per female/2
        # (should this be divided by 2?)
        pop_sus[k, 2, 1, t] <- (sum(pop_sus_proj[k, 1, 1:n_agef, t]) +
                            sum(pop_inf_proj[k, 1, 1:n_agef, t])) * fec[t] * (1 - psi[k, 2, 1, t]) / 2 

        ###################################################
        ### Infected/Infectious
        ###################################################

        ###########
        # Females
        ###########

        #Female: project forward anually
        for (a in 1:(n_agef - 1)) {
          pop_inf_proj[k, 1, a, t] <- pop_inf[k, 1, a, t - 1] * sn_inf[1, a, t - 1] +
                                  pop_sus[k, 1, a, t - 1] * sn_sus[1, a, t - 1] * psi[k, 1, a, t - 1]
        }
        #Female: max age = 9.5+ years
        pop_inf_proj[k, 1, n_agef, t] <- pop_inf_proj[k, 1, (n_agef - 1), t] +
                                      pop_inf[k, 1, n_agef, t - 1] * sn_inf[1, n_agef, t - 1] +
                                      # pop_sus_proj[1, (n_agef - 1), t] * psi[1, (n_agef - 1), t] + #need to double check this
                                      pop_sus[k, 1, n_agef, t - 1] * sn_sus[1, n_agef, t - 1] * psi[k, 1, n_agef, t - 1]

        #Female: set projection into population model matrix
        for (a in 2:n_agef) {
          pop_inf[k, 1, a, t] <- pop_inf_proj[k, 1, (a - 1), t]
        }
        
        #Female: fawn class = total #females * unisex fawns per female/2
        #(should this be divided by 2?)
        pop_inf[k, 1, 1, t] <- (sum(pop_sus_proj[k, 1, 1:n_agef, t]) +
                            sum(pop_inf_proj[k, 1, 1:n_agef, t])) * fec[t] * psi[k, 1, 1, t] / 2

        ###########
        # Males
        ###########

        #Male: project forward anually
        for (a in 1:(n_agem - 1)) {
            pop_inf_proj[k, 2, a, t] <- pop_inf[k, 2, a, t - 1] * sn_inf[2, a, t - 1] +
                                    pop_sus[k, 2, a, t - 1] * sn_sus[2, a, t - 1] * psi[k, 2, a, t - 1]
        }

        #Male: max age class
        pop_inf_proj[k, 2, n_agem, t] <- pop_inf_proj[k, 2, (n_agem - 1), t] +
                                          pop_inf[k, 2, n_agem, t - 1] * sn_inf[2, n_agem, t - 1] +
                                          # pop_sus_proj[2, (n_agem - 1), t] * psi +#need to double check this
                                          pop_sus[k, 2, n_agem, t - 1] *  sn_sus[2, n_agem, t - 1] * psi[k, 2, n_agem, t - 1]

        #Male: set projection into population model matrix
        for (a in 2:n_agem) {
          pop_inf[k, 2, a, t] <- pop_inf_proj[k, 2, (a - 1), t]
        }

        #Male: fawn class = total #females * unisex fawns per female/2#(should this be divided by 2?)
        pop_inf[k, 2, 1, t] <- (sum(pop_sus_proj[k, 1, 1:n_agef, t]) +
                            sum(pop_inf_proj[k, 1, 1:n_agef, t])) * fec[t] * psi[k, 2, 1, t] / 2

    }#end t
    }#end study_area

  ######################################################################
  ### Observation Model
  ######################################################################

  for(k in 1:n_study_area){

    for (i in 1:n_sex) {
      tau_obs[k, i] ~ dgamma(1, 1)
    }#end i


  for (t in 1:n_year) {
    for (a in 1:n_agef) {
      harv_pop[k, 1, a, t] <- (pop_inf[k, 1, a, t] * (1 - sh_inf[1, a, t]) + pop_sus[k, 1, a, t] * (1 - sh_sus[1, a, t])) * report[t]
    }
    for (a in 1:n_agem) {
      harv_pop[k, 2, a, t] <- (pop_inf[k, 2, a, t] * (1 - sh_inf[2, a, t]) + pop_sus[k, 2, a, t] * (1 - sh_sus[2, a, t])) * report[t]
    }

    #Total Antlerless Harvest
    #adding in male fawns
    mu_obs[k, 1, t] <- (sum(harv_pop[k, 1, 1:n_agef, t]) + harv_pop[k, 2, 1, t]) # eab_antlerless[t] * 

    #Total Antlered Harvest
    mu_obs[k, 2, t] <- sum(harv_pop[k, 2, 2:n_agem, t])#excludes male fawns eab_antlered[t] * 

    ###################################
    #Likelihood for overall total
    ###################################
    
    for (j in 1:n_sex) {
      O[k, j, t] ~ dnorm(mu_obs[k, j, t], tau_obs[k, j])
    }#end i

    ###################################
    #Likelihood for overall total
    ###################################

    ###parameters for likelihood harvest data by antlerless group
    p_less[k, 1, t] <- harv_pop[k, 1, 1, t] / mu_obs[k, 1, t]#proportion female fawns
    p_less[k, 2, t] <- harv_pop[k, 1, 2, t] / mu_obs[k, 1, t]#1
    p_less[k, 3, t] <- harv_pop[k, 1, 3, t] / mu_obs[k, 1, t]#2
    p_less[k, 4, t] <- harv_pop[k, 1, 4, t] / mu_obs[k, 1, t]#3
    p_less[k, 5, t] <- sum(harv_pop[k, 1, 5:6, t]) / mu_obs[k, 1, t]#4-5
    p_less[k, 6, t] <- sum(harv_pop[k, 1, 7:9, t]) / mu_obs[k, 1, t]#6-8
    p_less[k, 7, t] <- harv_pop[k, 1, 10, t] / mu_obs[k, 1, t]#9+
    p_less[k, 8, t] <- 1 - sum(p_less[k, 1:7, t])#proportion male fawns, antlerless

    #harvest data bt antlered group
    p_ant[k, 1, t] <- harv_pop[k, 2, 2, t] / mu_obs[k, 2, t]#1
    p_ant[k, 2, t] <- harv_pop[k, 2, 3, t] / mu_obs[k, 2, t]#2
    p_ant[k, 3, t] <- harv_pop[k, 2, 4, t] / mu_obs[k, 2, t]#3
    p_ant[k, 4, t] <- sum(harv_pop[k, 2, 5:6, t]) / mu_obs[k, 2, t]#4-5
    p_ant[k, 5, t] <- 1 - sum(p_ant[k, 1:4, t]) #6+

  }# end t

  for(t in 1:n_year){
      #antlerless, male fawns
      Cage_less[k, 1:(n_ageclassf + 1),t] ~ dmulti(prob = p_less[k,1:(n_ageclassf + 1),t],
                            size = sizeCage_f[k, t])  
      #antlered
      Cage_ant[k, 1:(n_ageclassm - 1),t] ~ dmulti(prob = p_ant[k, 1:(n_ageclassm - 1),t],
                size = sizeCage_m[k, t])
    }
  }

})#end model statement


##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################

#######################################
### Data for Model Fitting
#######################################

nimData <- list(Z_period = Z_period,
                Z_age = Z_age,
                num_period = num_period,
                adj_period = adj_period,
                weights_period = weights_period,
                age_lookup_f = age_lookup_f,
                age_lookup_m = age_lookup_m,
                # period_effect_survival = rep(NA,nT_period_overall),
                y_hunt_pos = 1,
                hunt_pos_ageweeks = d_fit_hunt_pos$ageweeks,
                hunt_pos_sex = d_fit_hunt_pos$sex,
                hunt_pos_age2date = d_fit_hunt_pos$birthweek - 1,
                hunt_pos_n_cases = d_fit_hunt_pos$n_cases,
                y_hunt_neg = 1,
                hunt_neg_ageweeks = d_fit_hunt_neg$ageweeks,
                hunt_neg_sex = d_fit_hunt_neg$sex,
                hunt_neg_age2date = d_fit_hunt_neg$birthweek - 1,
                hunt_neg_n_cases = d_fit_hunt_neg$n_cases,
                y_sus_cens_posttest = 1,
                sus_cens_posttest_left_age_e = d_fit_sus_cens_posttest$left_age_e,
                sus_cens_posttest_right_age_r = d_fit_sus_cens_posttest$right_age_r,
                sus_cens_posttest_sex = d_fit_sus_cens_posttest$sex,
                sus_cens_posttest_age2date = sus_cens_posttest_age2date,
                y_sus_cens_postno = 1,
                sus_cens_postno_left_age_e = c(d_fit_sus_cens_postno$left_age_e, d_fit_endlive$left_age_e),
                sus_cens_postno_right_age_r = c(d_fit_sus_cens_postno$right_age_r, d_fit_endlive$right_age_r),
                sus_cens_postno_sex = c(d_fit_sus_cens_postno$sex, d_fit_endlive$sex),
                sus_cens_postno_age2date = c(sus_cens_postno_age2date, endlive_age2date),
                y_sus_mort_posttest = 1,
                sus_mort_posttest_left_age_e = d_fit_sus_mort_posttest$left_age_e,
                sus_mort_posttest_right_age_r = d_fit_sus_mort_posttest$right_age_r,
                sus_mort_posttest_right_age_s = d_fit_sus_mort_posttest$right_age_s,
                sus_mort_posttest_sex = d_fit_sus_mort_posttest$sex,
                sus_mort_posttest_age2date = sus_mort_posttest_age2date,
                y_sus_mort_postno = 1,
                sus_mort_postno_left_age_e = d_fit_sus_mort_postno$left_age_e,
                sus_mort_postno_right_age_r = d_fit_sus_mort_postno$right_age_r,
                sus_mort_postno_right_age_s = d_fit_sus_mort_postno$right_age_s,
                sus_mort_postno_dn1 = d_fit_sus_mort_postno$left_age_e,
                sus_mort_postno_sex = d_fit_sus_mort_postno$sex,
                sus_mort_postno_age2date = sus_mort_postno_age2date,
                y_icap_cens = 1,
                icap_cens_left_age_e = d_fit_icap_cens$left_age_e,
                icap_cens_right_age_r = d_fit_icap_cens$right_age_r,
                icap_cens_sex = d_fit_icap_cens$sex,
                icap_cens_age2date = icap_cens_age2date,
                y_icap_mort = 1,
                icap_mort_left_age_e = d_fit_icap_mort$left_age_e,
                icap_mort_right_age_r = d_fit_icap_mort$right_age_r,
                icap_mort_right_age_s = d_fit_icap_mort$right_age_s,
                icap_mort_sex = d_fit_icap_mort$sex,
                icap_mort_age2date = icap_mort_age2date,
                y_rec_neg_cens_posttest = 1,
                rec_neg_cens_posttest_left_age_e = d_fit_rec_neg_cens_posttest$left_age_e,
                rec_neg_cens_posttest_right_age_r = d_fit_rec_neg_cens_posttest$right_age_r,
                rec_neg_cens_posttest_sex = d_fit_rec_neg_cens_posttest$sex,
                rec_neg_cens_posttest_age2date = rec_neg_cens_posttest_age2date,
                y_rec_neg_mort = 1,
                rec_neg_mort_left_age_e = d_fit_rec_neg_mort$left_age_e,
                rec_neg_mort_right_age_r = d_fit_rec_neg_mort$right_age_r,
                rec_neg_mort_right_age_s = d_fit_rec_neg_mort$right_age_s,
                rec_neg_mort_sex = d_fit_rec_neg_mort$sex,
                rec_neg_mort_age2date = rec_neg_mort_age2date,
                y_rec_pos_mort = 1,
                rec_pos_mort_left_age_e = d_fit_rec_pos_mort$left_age_e,
                rec_pos_mort_right_age_r = d_fit_rec_pos_mort$right_age_r,
                rec_pos_mort_right_age_s = d_fit_rec_pos_mort$right_age_s,
                rec_pos_mort_dn1 = d_fit_rec_pos_mort$left_age_e,
                rec_pos_mort_dn = d_fit_rec_pos_mort$ageweek_recap,
                rec_pos_mort_sex = d_fit_rec_pos_mort$sex,
                rec_pos_mort_age2date = rec_pos_mort_age2date,
                y_rec_pos_cens = 1,
                rec_pos_cens_left_age_e = d_fit_rec_pos_cens$left_age_e,
                rec_pos_cens_right_age_r = d_fit_rec_pos_cens$right_age_r,
                rec_pos_cens_dn1 = d_fit_rec_pos_cens$ageweek_recap,
                rec_pos_cens_dn = d_fit_rec_pos_cens$right_age_r,
                rec_pos_cens_sex = d_fit_rec_pos_cens$sex,
                rec_pos_cens_age2date = rec_pos_cens_age2date,
                y_idead = 1,
                idead_left_age_e = d_fit_idead$left_age_e,
                idead_right_age_r = d_fit_idead$right_age_r,
                idead_right_age_s = d_fit_idead$right_age_s,
                idead_dn1 = d_fit_idead$left_age_e,
                idead_dn = d_fit_idead$right_age_s,
                idead_sex = d_fit_idead$sex,
                idead_age2date = idead_age2date,
                y_aah = nrow(d_fit_aah),
                aah_ageweeks = d_fit_aah$ageweeks,
                aah_sex = d_fit_aah$sexnum,
                aah_age2date = d_fit_aah$age2date_weeks,
                aah_n = d_fit_aah$n,
                mort_hh = d_fit_hh$mort_h,
                sex_cause = d_fit_hh$sex,
                Z_cause_gun = Z_cause_gun,
                Z_cause_ng = Z_cause_ng,
                Cage_less = Cage_less,
                Cage_ant = Cage_ant,
                O = Ototal,
                f_logpop_sus = f_logpop_sus,
                f_logpop_inf = f_logpop_inf,
                m_logpop_sus = m_logpop_sus,
                m_logpop_inf = m_logpop_inf,
                obs_ct_fd_alpha = obs_ct_fd_alpha,
                obs_ct_fd_beta = obs_ct_fd_beta,
                Nfawn = fawndoe_df$overall_fawn,
                Ndoe = fawndoe_df$overall_doe
                )

#######################################
### Constants for MCMC
#######################################

nimConsts <- list(n_year = n_year,
    n_year_precollar = n_year_precollar,
    n_study_area = n_study_area,
    n_sex = n_sex,
    n_agef = n_agef,
    n_agem = n_agem,
    n_ageclassf = n_ageclassf,
    n_ageclassm = n_ageclassm,
    n_sex = n_sex,
    sizeCage_f = sizeCage_f,
    sizeCage_m = sizeCage_m,
    report_hyp_all = report_hyp_all,
    report_hyp_y = report_hyp_y,
    nT_period_overall = nT_period_overall,
    nT_period_overall_ext = nT_period_overall_ext,
    nT_period_precollar = nT_period_precollar,
    nT_period_collar = nT_period_collar,
    nT_age_surv = nT_age_surv,
    nT_age_surv_aah = nT_age_surv_aah,
    n_year_fec_early = n_year_fec_early,
    nknots_age = nknots_age,
    nknots_period = nknots_period,
    n_period_lookup = n_period_lookup,
    n_adj_period = n_adj_period,
    n_age_lookup_f = length(age_lookup_f),
    n_age_lookup_m = length(age_lookup_m),
    n_age_lookup_f_conv = n_age_lookup_f_conv,
    n_age_lookup_m_conv = n_age_lookup_m_conv,
    period_lookup = period_lookup,
    age_lookup_f_conv = age_lookup_f_conv,
    age_lookup_m_conv = age_lookup_m_conv,
    ng_start = d_fit_season$ng_start,
    gun_start = d_fit_season$gun_start,
    gun_end = d_fit_season$gun_end,
    ng_end = d_fit_season$ng_end,
    yr_start = d_fit_season$yr_start,
    yr_end = d_fit_season$yr_end,
    nInfHarvest = nrow(d_fit_hunt_pos),
    nSusHarvest = nrow(d_fit_hunt_neg),
    nSusCensTest = nrow(d_fit_sus_cens_posttest),
    nSusCensNo = nrow(d_fit_sus_cens_postno) +
                  nrow(d_fit_endlive),
    nSusMortTest = nrow(d_fit_sus_mort_posttest),
    nSusMortNoTest = nrow(d_fit_sus_mort_postno),
    nIcapCens = nrow(d_fit_icap_cens),
    nIcapMort = nrow(d_fit_icap_mort),
    nRecNegCensTest = nrow(d_fit_rec_neg_cens_posttest),
    nRecNegMort = nrow(d_fit_rec_neg_mort),
    nRecPosMort = nrow(d_fit_rec_pos_mort),
    nNegCapPosMort = nrow(d_fit_idead),
    nAAH = nrow(d_fit_aah),
    sect_hunt_pos = d_fit_hunt_pos$ew,
    sect_hunt_neg = d_fit_hunt_neg$ew,
    sect_sus_cens_posttest = d_fit_sus_cens_posttest$study_area,
    sect_sus_cens_postno = c(d_fit_sus_cens_postno$study_area,
                            d_fit_endlive$study_area),
    sect_sus_mort_posttest = d_fit_sus_mort_posttest$study_area,
    sect_sus_mort_postno = d_fit_sus_mort_postno$study_area,
    sect_icap_cens = d_fit_icap_cens$study_area,
    sect_icap_mort = d_fit_icap_mort$study_area,
    sect_rec_neg_cens_posttest = d_fit_rec_neg_cens_posttest$study_area,
    sect_rec_neg_mort = d_fit_rec_neg_mort$study_area,
    sect_rec_pos_mort = d_fit_rec_pos_mort$study_area,
    sect_rec_pos_cens = d_fit_rec_pos_cens$study_area,
    sect_idead = d_fit_idead$study_area,
    sect_aah = d_fit_aah$study_area,
    records_cause = records_cause,
    interval = d_fit_hh$right_period_s - 1,
    intvl_step_yr = intvl_step_yr_weekly,
    lookup_pe_surv = lookup_pe_surv
    )


#######################################
### Initial Values for MCMC
#######################################

initsFun <- function()list(
                          beta_male = rnorm(1, -.5, .01),
                          beta0_sus_temp = rnorm(1, -5.5, 0.0001),
                          sus_mix = 1,
                          beta0_inf_temp = rnorm(1, -4, 0.0001),
                          inf_mix = 1,
                          b_age_survival = rnorm(nknots_age) * 10^-4,
                          b_period_survival = rnorm(nknots_period) * 10^-4,
                          tau_period_survival = runif(1, .1, 1),
                          tau_age_survival = runif(1, .1, .4),
                          tau_age_foi_male = runif(1, 1.5, 1.7),
                          tau_age_foi_female = runif(1, 2.7, 4.2),
                          tau_period_foi_male = runif(1, 4.2, 6.8),
                          tau_period_foi_female = runif(1, 2.27, 3.44),
                          m_period_foi = seq(-2, 2, length = n_year),
                          f_period_foi = seq(-2, 2, length = n_year),
                          m_age_foi = seq(-6, -4, length = n_agem),
                          f_age_foi = seq(-7, -5, length = n_agef),
                          tau_period_precollar = rgamma(1,1,1),
                          period_annual_survival = rnorm(n_year_precollar + 1),
                          beta0_cause = rnorm(1,0,1),
                          beta_cause_male = rnorm(1,0,1),
                          beta_cause_gun = rnorm(1,0,1),
                          beta_cause_ng = rnorm(1,0,1),
                          tau_obs = matrix(runif(4, 1, 3),2,2),
                          tau_pop = matrix(runif(4, .5, 1),2,2),
                          report_overall = report_overall_init,
                          report = report_init,
                          fec_epsilon = fec_eps_init,
                          mu_fec = rnorm(1, mu_fec_init, .01),
                          fec_prec_eps = runif(1, 5, 10),
                          space_temp = rnorm(1,-.55,.01)
                          )
nimInits <- initsFun()

# start_Rmodel <- Sys.time()
Rmodel <- nimbleModel(code = modelcode,
                      constants = nimConsts,
                      data = nimData,
                      inits = initsFun(),
                      calculate = FALSE,
                      check = FALSE
                      )
# end_Rmodel <- Sys.time() - start_Rmodel
Rmodel$initializeInfo()

Cnim <- compileNimble(Rmodel)

for(i in 1:10){beepr::beep(1)}

#######################################
### Parameters to trace in MCMC
#######################################

parameters <- c(
              "beta_male",
              "tau_age_foi_male",
              "tau_age_foi_female",
              "tau1_age_foi_male",
              "tau1_age_foi_female",
              "m_age_foi",
              "f_age_foi",
              "m_age_foi_mu",
              "f_age_foi_mu",
              "tau_period_foi_male",
              "tau_period_foi_female",
              "f_period_foi",
              "m_period_foi",
              "space",
              "space_temp",
              "space_alpha",
              "beta0_sus_temp",
              "sus_mix",
              "beta0_survival_sus",
              "tau_age_survival",
              "age_effect_survival",
              "tau_period_survival",
              "tau_period_precollar",
              "period_effect_survival",
              "beta0_inf_temp",
              "inf_mix",
              "beta0_survival_inf",
              "beta0_cause",
              "beta_cause_gun",
              "beta_cause_ng",
              "beta_cause_male",
              "p_nogun_f",
              "p_gun_f",
              "p_nogun_m",
              "p_gun_m",
              "report",
              "fec",
              "mu_fec",
              "fec_prec_eps",
              "fec_epsilon",
              "sn_inf",
              "sn_sus",
              "sh_inf",
              "sh_sus",
              "mu_obs",
              "tau_obs",
              "tau_pop"
               )

confMCMC <- configureMCMC(Rmodel,
                         monitors = parameters,
                         thin = 1,
                         # enableWAIC = TRUE,
                         useConjugacy = FALSE)
nimMCMC <- buildMCMC(confMCMC)
CnimMCMC <- compileNimble(nimMCMC,
                         project = Rmodel)
for(i in 1:10){beepr::beep(1)}

set.seed(7654321)
starttime <- Sys.time()
mcmcout <- runMCMC(CnimMCMC,
                  niter = 1000,
                  nburnin = 0,
                  nchains = 1,
                  inits = initsFun,
                  samplesAsCodaMCMC = TRUE,
                  summary = TRUE
                  )
runtime <- difftime(Sys.time(),
                    starttime,
                    units = "min")
for (i in 1:10) {beepr::beep(1)}


# end_Rmodel
# endtime_rmodel_compile
# endtime_mcmc
# runtime

sink("runtime_allsteps.txt")
# cat("Rmodel:\n")
# end_Rmodel
# cat("\nCompile Rmodel:\n")
# endtime_rmodel_compile
# cat("\nCompile MCMC:\n")
# endtime_mcmc
cat("\nRun MCMC 1000 iter: ",runtime)
sink()
#############################################################
###
### Running in parallel, restartable
###
#############################################################

reps  <- 1000
bin <- reps * .5
n_thin <- 1
n_chains <- 3
starttime <- Sys.time()
cl <- makeCluster(n_chains, timeout = 5184000)

clusterExport(cl, c("modelcode",
                    "initsFun",
                    "nimData",
                    "nimConsts",
                    "parameters",
                    "reps",
                    "bin",
                    "n_thin",
                    "set_period_effects_constant",
                    "dInfHarvest",
                    "dSusHarvest",
                    "dSusCensTest",
                    "dSusCensNo",
                    "dSusMortTest",
                    "dSusMortNoTest",
                    "dIcapCens",
                    "dIcapMort",
                    "dRecNegCensTest",
                    "dRecNegMort",
                    "dRecPosMort",
                    "dRecPosCens",
                    "dNegCapPosMort",
                    "dAAH",
                    "calc_surv_aah",
                    "calc_surv_harvest",
                    "calc_infect_prob"
                    ))
for (j in seq_along(cl)) {
  set.seed(j + 1000)
  init <- initsFun()
  clusterExport(cl[j], "init")
}
for (i in 1:10) {beepr::beep(1)}

starttime <- Sys.time()
mcmcout1 <-  mcmc.list(clusterEvalQ(cl, {
  library(nimble)
  library(coda)

  assign("set_period_effects_constant", set_period_effects_constant, envir = .GlobalEnv)
  assign("dInfHarvest", dInfHarvest, envir = .GlobalEnv)
  assign("dSusHarvest", dSusHarvest, envir = .GlobalEnv)
  assign("dSusCensTest", dSusCensTest, envir = .GlobalEnv)
  assign("dSusCensNo", dSusCensNo, envir = .GlobalEnv)
  assign("dSusMortTest", dSusMortTest, envir = .GlobalEnv)
  assign("dSusMortNoTest", dSusMortNoTest, envir = .GlobalEnv)
  assign("dIcapCens", dIcapCens, envir = .GlobalEnv)
  assign("dIcapMort", dIcapMort, envir = .GlobalEnv)
  assign("dRecNegCensTest", dRecNegCensTest, envir = .GlobalEnv)
  assign("dRecNegMort", dRecNegMort, envir = .GlobalEnv)
  assign("dRecPosMort", dRecPosMort, envir = .GlobalEnv)
  assign("dRecPosCens", dRecPosCens, envir = .GlobalEnv)
  assign("dNegCapPosMort", dNegCapPosMort, envir = .GlobalEnv)
  assign("dAAH", dAAH, envir = .GlobalEnv)
  assign("calc_surv_aah", calc_surv_aah, envir = .GlobalEnv)
  assign("calc_surv_harvest", calc_surv_harvest, envir = .GlobalEnv)
  assign("calc_infect_prob", calc_infect_prob, envir = .GlobalEnv)

  ##############################################################
  ###
  ### Execute MCMC
  ###
  ##############################################################

  Rmodel <- nimbleModel(code = modelcode,
                        name = "modelcode",
                        constants = nimConsts,
                        data = nimData,
                        inits = initsFun)
  Cnim <- compileNimble(Rmodel)
  confMCMC <- configureMCMC(Rmodel,
                            monitors = parameters,
                            thin = n_thin,
                            useConjugacy = FALSE)
  nimMCMC <- buildMCMC(confMCMC)
  CnimMCMC <- compileNimble(nimMCMC,
                            project = Rmodel)

  CnimMCMC$run(reps, reset = FALSE)

  return(as.mcmc(as.matrix(CnimMCMC$mvSamples)))
}))
runtime1 <- difftime(Sys.time(), starttime, units = "min")

# for(chn in 1:nc) { # nc must be > 1
#   ind.keep <- c()
#   for(p in 1:length(parameters)) ind.keep <-
#       c(ind.keep, which(str_detect(dimnames(out1[[chn]])[[2]], parameters[p]))) %>% unique()
#   out1[[chn]] <- out1[[chn]][,ind.keep]
# }

# ## Check convergence ##
# out2 <- out1
# ni.saved <- nrow(out2[[1]])
# for(chn in 1:nc) { # nc must be > 1
  
#   if(nb < 1) {
#     nb.real <- (round(ni.saved * nb)+1)
#   } else {
#     nb.real <- (round(nb/nt)+1)
#   }
#   out2[[chn]] <- out2[[chn]][nb.real:ni.saved,]
# }
# out.mcmc <- coda::as.mcmc.list(lapply(out2, coda::as.mcmc))
stopCluster(cl)

save(mcmcout1, file = "mcmcout1.Rdata")
save(runtime1, file = "runtime1.Rdata")
# save(endtime_rmodel_compile, file = "endtime_rmodel_compile.Rdata")
# save(endtime_mcmc, file = "endtime_mcmc.Rdata")

#not calculating waic, because too many params would need to be traced
# posteriorSamplesMatrix <- rbind(mcmcout[[1]], mcmcout[[2]], mcmcout[[3]])
# CnimMCMC$run(5)   ## non-zero number of iterations
# nimble:::matrix2mv(posteriorSamplesMatrix, CnimMCMC$mvSamples)
# # CnimMCMC$enableWAIC <- TRUE
# waic_spline <- calculateWAIC(posteriorSamplesMatrix,Rmodel)


# waic_spline_covs <- mcmcout$WAIC
# save(waic_spline, file = "waic_spline.Rdata")



###
### save model run
###

# save(runtime,file="results/runtime.Rdata")
# save(mcmcout,file="results/mcmcout.Rdata")