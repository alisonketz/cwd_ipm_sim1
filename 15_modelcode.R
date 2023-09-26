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

  beta_male ~ dnorm(0, .1)

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

  ### Period effects 
  # tau_period_foi_male  ~ dgamma(1, 1)
  # tau_period_foi_female  ~ dgamma(1, 1)

  ###ICAR specification
  # f_period_foi[1:n_year] ~ dcar_normal(adj = adj_period[1:n_adj_period],
  #                                 weights = weights_period[1:n_adj_period],
  #                                 num = num_period[1:n_year],
  #                                 tau = tau_period_foi_female,
  #                                 zero_mean = 1)
  
  # m_period_foi[1:n_year] ~ dcar_normal(adj = adj_period[1:n_adj_period],
  #                                 weights = weights_period[1:n_adj_period],
  #                                 num = num_period[1:n_year],
  #                                 tau = tau_period_foi_male,
  #                                 zero_mean = 1)

  ### RW1 Specification
  # tau1_period_foi_f <- .0000001 * tau_period_foi_female
  # tau1_period_foi_m <- .0000001 * tau_period_foi_male
  # f_period_foi_temp[1] ~ dnorm(0, tau1_period_foi_f)
  # m_period_foi_temp[1] ~ dnorm(0, tau1_period_foi_m)
  # for (t in 2:n_year) {
  #   f_period_foi_temp[t] ~ dnorm(f_period_foi_temp[t - 1], tau_period_foi_female)
  #   m_period_foi_temp[t] ~ dnorm(m_period_foi_temp[t - 1], tau_period_foi_male)
  # }
  # f_period_foi_mu <- mean(f_period_foi_temp[1:n_year])
  # m_period_foi_mu <- mean(m_period_foi_temp[1:n_year])
  # for (t in 1:n_year) {
  #   f_period_foi[t] <- f_period_foi_temp[t] - f_period_foi_mu
  #   m_period_foi[t] <- m_period_foi_temp[t] - m_period_foi_mu
  # }
  
  ### removing first years prior to surveillance data
  # tau1_period_foi_f <- .0000001 * tau_period_foi_female
  # tau1_period_foi_m <- .0000001 * tau_period_foi_male
  # f_period_foi_temp[9] ~ dnorm(0, tau1_period_foi_f)
  # m_period_foi_temp[9] ~ dnorm(0, tau1_period_foi_m)
  # for (t in 10:n_year) {
  #   f_period_foi_temp[t] ~ dnorm(f_period_foi_temp[t - 1], tau_period_foi_female)
  #   m_period_foi_temp[t] ~ dnorm(m_period_foi_temp[t - 1], tau_period_foi_male)
  # }
  # f_period_foi_mu <- mean(f_period_foi_temp[9:n_year])
  # m_period_foi_mu <- mean(m_period_foi_temp[9:n_year])
  # for (t in 9:n_year) {
  #   f_period_foi[t] <- f_period_foi_temp[t] - f_period_foi_mu
  #   m_period_foi[t] <- m_period_foi_temp[t] - m_period_foi_mu
  # }

  ### setting early years prior to surveillance data
  # tau1_period_foi_f <- .0000001 * tau_period_foi_female
  # tau1_period_foi_m <- .0000001 * tau_period_foi_male
  # f_period_foi_temp[9] ~ dnorm(0, tau1_period_foi_f)
  # m_period_foi_temp[9] ~ dnorm(0, tau1_period_foi_m)
  # for (t in 10:n_year) {
  #   f_period_foi_temp[t] ~ dnorm(f_period_foi_temp[t - 1], tau_period_foi_female)
  #   m_period_foi_temp[t] ~ dnorm(m_period_foi_temp[t - 1], tau_period_foi_male)
  # }
  # for(t in 1:8) {
  #   f_period_foi_temp[t] <- f_period_foi_temp[9]
  #   m_period_foi_temp[t] <- m_period_foi_temp[9]
  # }
  # f_period_foi_mu <- mean(f_period_foi_temp[1:n_year])
  # m_period_foi_mu <- mean(m_period_foi_temp[1:n_year])
  # for (t in 1:n_year) {
  #   f_period_foi[t] <- f_period_foi_temp[t] - f_period_foi_mu
  #   m_period_foi[t] <- m_period_foi_temp[t] - m_period_foi_mu
  # }

  ##########################################################
  ### FOI period effect cgam, sex-specific
  ##########################################################

  # tau_period_foi_male ~ dgamma(1, 1)
  # tau_period_foi_female ~ dgamma(1, 1)

  # for (k in 1:nknots_foi_cgam) {
  #   ln_b_foi_cgam_f[k] ~ dnorm(0, tau_period_foi_female)
  #   b_foi_cgam_f[k] <- exp(ln_b_foi_cgam_f[k])

  #   ln_b_foi_cgam_m[k] ~ dnorm(0, tau_period_foi_male)
  #   b_foi_cgam_m[k] <- exp(ln_b_foi_cgam_m[k])
  # }
  # for (t in 1:n_year) {
  #   m_period_foi_cgam_temp[t] <- inprod(b_foi_cgam_m[1:nknots_foi_cgam],
  #                                Z_foi_cgam[t, 1:nknots_foi_cgam])
  #   m_period_foi_cgam[t] <- m_period_foi_cgam_temp[t] - m_mu_foi_cgam

  #   f_period_foi_cgam_temp[t] <- inprod(b_foi_cgam_f[1:nknots_foi_cgam],
  #                                Z_foi_cgam[t, 1:nknots_foi_cgam])
  #   f_period_foi_cgam[t] <- f_period_foi_cgam_temp[t] - f_mu_foi_cgam

  # }
  # m_mu_foi_cgam <- mean(m_period_foi_cgam_temp[1:n_year])
  # f_mu_foi_cgam <- mean(f_period_foi_cgam_temp[1:n_year])

  # ##########################################################
  # ### FOI period effect additive spline
  # ##########################################################

  # tau_foi_spline_male ~ dgamma(.01, .01)
  # tau_foi_spline_female ~ dgamma(.01, .01)
  # for (k in 1:nknots_foi_spline) {
  #   f_b_foi_spline[k] ~ ddexp(0, tau_foi_spline_female)
  #   m_b_foi_spline[k] ~ ddexp(0, tau_foi_spline_male)
  # }
  # for (t in 1:(n_year/2)) {
  #   m_period_foi_spline[t] <- inprod(m_b_foi_spline[1:nknots_foi_spline],
  #                                  Z_foi_spline[t, 1:nknots_foi_spline])
  #   f_period_foi_spline[t] <- inprod(f_b_foi_spline[1:nknots_foi_spline],
  #                                  Z_foi_spline[t, 1:nknots_foi_spline])
  # }

  # f_period_foi[1:(n_year/2)] <- f_period_foi_cgam[1:(n_year/2)]
  # m_period_foi[1:(n_year/2)] <- m_period_foi_cgam[1:(n_year/2)]
  # f_period_foi[(n_year * .5 + 1):n_year] <- f_period_foi_cgam[(n_year * .5 + 1):n_year] + f_period_foi_spline[1:(n_year * .5)]
  # m_period_foi[(n_year * .5 + 1):n_year] <- m_period_foi_cgam[(n_year * .5 + 1):n_year] + m_period_foi_spline[1:(n_year * .5)]

  ##########################################################
  ### Random effect for East/West spatial model
  ##########################################################

  space[1] <- 0
  space[2]  <- space_temp * space_mix

  space_temp ~ dnorm(0, 1)
  space_mix ~ dunif(-1, 1)

  ############################################################
  ############################################################
  ### Age/Period Survival Model
  ############################################################
  ############################################################

  ####################################
  ### Susceptibles survival intercept
  ####################################

  # beta0_survival_sus ~ dnorm(0, .1)
  # beta0_survival_sus ~ T(dnorm(-6, .1),,0)
  beta0_sus_temp ~ dnorm(0, .1)
  sus_mix ~ dunif(-1, 1)
  beta0_survival_sus <- beta0_sus_temp * sus_mix

  ##################################
  ### Infected survival intercept
  ##################################

  # beta0_survival_inf ~ dnorm(0, .1)
  # beta0_survival_inf ~ T(dnorm(-6, .1),,0)
  beta0_inf_temp ~ dnorm(0, .1)
  inf_mix ~ dunif(-1, 1)
  beta0_survival_inf <- beta0_inf_temp * inf_mix

  ########################################
  ### Priors for Age Effects Survival
  ########################################

  ### Age effects
  for (k in 1:nknots_age) {
    # ln_b_age_survival[k] ~ dnorm(0, tau_age_survival)
    ln_b_age_survival[k] ~ dnorm(0, 1)
    b_age_survival[k] <- exp(ln_b_age_survival[k])
  }
  # tau_age_survival ~ dgamma(1, 1)

  for (t in 1:nT_age_surv) {
    age_effect_survival_temp[t] <- inprod(b_age_survival[1:nknots_age],
                                     Z_age[t, 1:nknots_age])
  }
  mu_age_effect_survival_temp <- mean(age_effect_survival_temp[1:nT_age_surv])
  
  for (t in 1:nT_age_surv) {
    age_effect_survival[t] <-  age_effect_survival_temp[t] -
                               mu_age_effect_survival_temp
  }

  ########################################
  ### Priors for Period Effects Survival
  ########################################

  # ### Period effects from collar data
  # for (k in 1:nknots_period) {
  #   b_period_survival[k] ~ dnorm(0, tau_period_survival)
  # }
  # tau_period_survival ~ dgamma(1, 1)
  # for (t in 1:nT_period_collar) {
  #   period_effect_surv[t] <- inprod(b_period_survival[1:nknots_period],
  #                                   Z_period[t, 1:nknots_period]) +
  #                             Z_collar_gun[t] * beta_harvest_gun +
  #                             Z_collar_ng[t] * beta_harvest_ng
  # }

  # #additive effect of harvest on total mortality hazard
  # beta_harvest_gun ~ dnorm(0, .01)
  # beta_harvest_ng ~ dnorm(0, .01)

  # ### Period effects from aah data
  # tau_period_precollar ~ dgamma(1,1)
  # for (k in 1:(n_year_precollar + 1)) {
  #   period_annual_survival[k] ~ dnorm(0, tau_period_precollar)
  # }

  #Period effects from aah data - multiple intercepts version
  # tau_period_precollar ~ dgamma(1, 1)
  # period_int_survival ~ dnorm(0, tau_period_precollar)
  # period_annual_survival[1:18] <- period_int_survival
  # period_annual_survival[19:(n_year_precollar+1)] <- 0
                                                                                                                                                                                                                                                                       
  # period_effect_survival[1:nT_period_overall_ext] <- set_period_effects_constant(
  #       n_year_precollar = n_year_precollar,
  #       n_year_precollar_ext = n_year_precollar_ext,
  #       n_year_prestudy_ext = n_year_prestudy_ext,
  #       nT_period_precollar_ext = nT_period_precollar_ext,
  #       nT_period_precollar = nT_period_precollar,
  #       nT_period_collar = nT_period_collar,
  #       nT_period_overall_ext = nT_period_overall_ext,
  #       nT_period_prestudy_ext = nT_period_prestudy_ext,
  #       yr_start = yr_start[1:n_year],
  #       yr_end = yr_end[1:n_year],
  #       period_effect_surv = period_effect_surv[1:nT_period_collar],
  #       period_annual_survival = period_annual_survival[1:(n_year_precollar + 1)]
  # )

#   period_effect_survival[1:nT_period_overall_ext] <- set_period_effects_ave(
#         n_year_precollar = n_year_precollar,
#         nT_period_precollar = nT_period_precollar,
#         nT_period_collar = nT_period_collar,
#         nT_period_overall = nT_period_overall,
#         nT_period_overall_ext = nT_period_overall_ext,
#         nT_period_prestudy_ext = nT_period_prestudy_ext,
#         yr_start = yr_start[1:n_year],
#         yr_end = yr_end[1:n_year],
#         period_effect_surv = period_effect_surv[1:nT_period_collar],
#         period_annual_survival = period_annual_survival[1:(n_year_precollar + 1)],
#         indx_mat_pe_surv = indx_mat_pe_surv[1:6,1:intvl_step_yr],
#         intvl_step_yr = intvl_step_yr
#   )



  #######################################################################
  ###
  ###   User defined distribution for likelihood for
  ###   all harvested deer without joint estimation w/ survival
  ###
  ###   d_fit_hunt
  ###
  #######################################################################

  # y_hunt ~ dFOIhunt(
  #                 test_status = hunt_test_status[1:nHarvest],
  #                 n_cases = hunt_n_cases[1:nHarvest],
  #                 n_samples = nHarvest,
  #                 a = hunt_ageweeks[1:nHarvest], #age (weeks) at harvest
  #                 sex = hunt_sex[1:nHarvest],
  #                 age2date = hunt_age2date[1:nHarvest],
  #                 f_age_foi = f_age_foi[1:n_ageclassf],
  #                 m_age_foi = m_age_foi[1:n_ageclassm],
  #                 age_lookup_f = age_lookup_f[1:nT_age_surv],
  #                 age_lookup_m = age_lookup_m[1:nT_age_surv],
  #                 period_lookup_foi = period_lookup_foi[1:nT_period_overall_ext],
  #                 f_period_foi = f_period_foi[1:n_year],
  #                 m_period_foi = m_period_foi[1:n_year],
  #                 space = space[1:n_study_area],
  #                 sect = sect_hunt[1:nHarvest]
  #                 )


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
                  period_effect_surv = period_effect_survival[1:nT_period_overall_ext],
                  f_age_foi = f_age_foi[1:n_ageclassf],
                  m_age_foi = m_age_foi[1:n_ageclassm],
                  age_lookup_f = age_lookup_f[1:nT_age_surv],
                  age_lookup_m = age_lookup_m[1:nT_age_surv],
                  period_lookup_foi = period_lookup_foi[1:nT_period_overall_ext],
                  f_period_foi = f_period_foi[1:n_year],
                  m_period_foi = m_period_foi[1:n_year],
                  space = space[1:n_study_area],
                  sect = sect_hunt_pos[1:nInfHarvest]
                  )

######################################################################
##
##   User defined distribution for likelihood for
##   uninfected harvest deer
##   d_fit_hunt_neg
##   Overleaf Equation (5)
##
######################################################################

  y_hunt_neg ~ dSusHarvest(
                  n_cases = hunt_neg_n_cases[1:nInfHarvest],
                  n_samples = nSusHarvest,
                  a = hunt_neg_ageweeks[1:nSusHarvest], #age (weeks) at harvest
                  sex = hunt_neg_sex[1:nSusHarvest],
                  age2date = hunt_neg_age2date[1:nSusHarvest],
                  beta_male = beta_male,
                  beta0_sus = beta0_survival_sus,
                  age_effect_surv = age_effect_survival[1:nT_age_surv],
                  period_effect_surv = period_effect_survival[1:nT_period_overall_ext],
                  f_age_foi = f_age_foi[1:n_ageclassf],
                  m_age_foi = m_age_foi[1:n_ageclassm],
                  age_lookup_f = age_lookup_f[1:nT_age_surv],
                  age_lookup_m = age_lookup_m[1:nT_age_surv],
                  period_lookup_foi = period_lookup_foi[1:nT_period_overall_ext],
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
        period_effect_surv = period_effect_survival[1:nT_period_overall_ext],
        f_age_foi = f_age_foi[1:n_ageclassf],
        m_age_foi = m_age_foi[1:n_ageclassm],
        age_lookup_f = age_lookup_f[1:nT_age_surv],
        age_lookup_m = age_lookup_m[1:nT_age_surv],
        period_lookup_foi = period_lookup_foi[1:nT_period_overall_ext],
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
        period_effect_surv = period_effect_survival[1:nT_period_overall_ext],
        f_age_foi = f_age_foi[1:n_ageclassf],
        m_age_foi = m_age_foi[1:n_ageclassm],
        age_lookup_f = age_lookup_f[1:nT_age_surv],
        age_lookup_m = age_lookup_m[1:nT_age_surv],
        period_lookup_foi = period_lookup_foi[1:nT_period_overall_ext],
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
      period_effect_surv = period_effect_survival[1:nT_period_overall_ext],
      f_age_foi = f_age_foi[1:n_ageclassf],
      m_age_foi = m_age_foi[1:n_ageclassm],
      age_lookup_f = age_lookup_f[1:nT_age_surv],
      age_lookup_m = age_lookup_m[1:nT_age_surv],
      period_lookup_foi = period_lookup_foi[1:nT_period_overall_ext],
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
      period_effect_surv = period_effect_survival[1:nT_period_overall_ext],
      f_age_foi = f_age_foi[1:n_ageclassf],
      m_age_foi = m_age_foi[1:n_ageclassm],
      age_lookup_f = age_lookup_f[1:nT_age_surv],
      age_lookup_m = age_lookup_m[1:nT_age_surv],
      period_lookup_foi = period_lookup_foi[1:nT_period_overall_ext],
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
    period_effect_surv = period_effect_survival[1:nT_period_overall_ext],
    f_age_foi = f_age_foi[1:n_ageclassf],
    m_age_foi = m_age_foi[1:n_ageclassm],
    age_lookup_f = age_lookup_f[1:nT_age_surv],
    age_lookup_m = age_lookup_m[1:nT_age_surv],
    period_lookup_foi = period_lookup_foi[1:nT_period_overall_ext],
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
    period_effect_surv = period_effect_survival[1:nT_period_overall_ext],
    f_age_foi = f_age_foi[1:n_ageclassf],
    m_age_foi = m_age_foi[1:n_ageclassm],
    age_lookup_f = age_lookup_f[1:nT_age_surv],
    age_lookup_m = age_lookup_m[1:nT_age_surv],
    period_lookup_foi = period_lookup_foi[1:nT_period_overall_ext],
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
      period_effect_surv = period_effect_survival[1:nT_period_overall_ext],
      f_age_foi = f_age_foi[1:n_ageclassf],
      m_age_foi = m_age_foi[1:n_ageclassm],
      age_lookup_f = age_lookup_f[1:nT_age_surv],
      age_lookup_m = age_lookup_m[1:nT_age_surv],
      period_lookup_foi = period_lookup_foi[1:nT_period_overall_ext],
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
      period_effect_surv = period_effect_survival[1:nT_period_overall_ext],
      f_age_foi = f_age_foi[1:n_ageclassf],
      m_age_foi = m_age_foi[1:n_ageclassm],
      age_lookup_f = age_lookup_f[1:nT_age_surv],
      age_lookup_m = age_lookup_m[1:nT_age_surv],
      period_lookup_foi = period_lookup_foi[1:nT_period_overall_ext],
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
      period_effect_surv = period_effect_survival[1:nT_period_overall_ext],
      f_age_foi = f_age_foi[1:n_ageclassf],
      m_age_foi = m_age_foi[1:n_ageclassm],
      age_lookup_f = age_lookup_f[1:nT_age_surv],
      age_lookup_m = age_lookup_m[1:nT_age_surv],
      period_lookup_foi = period_lookup_foi[1:nT_period_overall_ext],
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
      period_effect_surv = period_effect_survival[1:nT_period_overall_ext],
      f_age_foi = f_age_foi[1:n_ageclassf],
      m_age_foi = m_age_foi[1:n_ageclassm],
      age_lookup_f = age_lookup_f[1:nT_age_surv],
      age_lookup_m = age_lookup_m[1:nT_age_surv],
      period_lookup_foi = period_lookup_foi[1:nT_period_overall_ext],
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
      period_effect_surv = period_effect_survival[1:nT_period_overall_ext],
      f_age_foi = f_age_foi[1:n_ageclassf],
      m_age_foi = m_age_foi[1:n_ageclassm],
      age_lookup_f = age_lookup_f[1:nT_age_surv],
      age_lookup_m = age_lookup_m[1:nT_age_surv],
      period_lookup_foi = period_lookup_foi[1:nT_period_overall_ext],
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
      period_effect_surv = period_effect_survival[1:nT_period_overall_ext],
      f_age_foi = f_age_foi[1:n_ageclassf],
      m_age_foi = m_age_foi[1:n_ageclassm],
      age_lookup_f = age_lookup_f[1:nT_age_surv],
      age_lookup_m = age_lookup_m[1:nT_age_surv],
      period_lookup_foi = period_lookup_foi[1:nT_period_overall_ext],
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

  ### priors
  ### sex-specific hunt probability given mortalities
  beta0_cause ~ dnorm(0, .1)
  beta_cause_gun ~ dnorm(0, .1)
  beta_cause_ng ~ dnorm(0, .1)
  beta_cause_male ~ dnorm(0, .1)

  #######################################################
  ###
  ### Likelihood for cause of death
  ###
  #######################################################

  for (i in 1:records_cause) {
      p_cause[i]  <- ilogit(beta0_cause +
                            Z_cause_ng[interval_cause[i]] * beta_cause_ng +
                            Z_cause_gun[interval_cause[i]] * beta_cause_gun +
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


  # #######################################
  # ### Initial population
  # ########################################
  # # tau_pop ~ dgamma(1, 1)

  # #should this be different for pos/neg or m/f for study area?
  # # for(i in 1:2) {
  # #   tau_pop[i] ~ dgamma(10, 1)
  # # }
  # for(k in 1:n_study_area) {
  #     for (a in 1:n_agef) {

  #     #Initial population structure pop[sex,age,year] for susceptible deer
  #     llpop_sus[k, 1, a, 1]  <- f_logpop_sus[k, a]
  #     pop_sus[k, 1, a, 1] <- exp(llpop_sus[k, 1, a, 1])

  #     #Initial population structure pop[study_area=k,sex=i,year=t,age=a]
  #     llpop_inf[k, 1, a, 1] <- f_logpop_inf[k, a]
  #     pop_inf[k, 1, a, 1] <- exp(llpop_inf[k, 1, a, 1])
  #     }
  #     for (a in 1:n_agem) {
  #         ### Initial population structure pop
  #         ### [study_area,sex,age,period(year)] for susceptible deer
  #         llpop_sus[k, 2, a, 1]  <- m_logpop_sus[k, a]
  #         pop_sus[k, 2, a, 1] <- exp(llpop_sus[k, 2, a, 1])

  #         #Initial population structure pop for infected deer
  #         llpop_inf[k, 2, a, 1] <- m_logpop_inf[k, a]
  #         pop_inf[k, 2, a, 1] <- exp(llpop_inf[k, 2, a, 1])
  #     }
  # }


  # # for(k in 1:n_study_area) {
  # #   for (a in 1:n_agef) {

  # #     #Initial population structure pop[sex,age,year] for susceptible deer
  # #     llpop_sus[k, 1, a, 1] ~ dnorm(f_logpop_sus[k, a], tau_pop[1])#tau_pop[1]
  # #     pop_sus[k, 1, a, 1] <- exp(llpop_sus[k, 1, a, 1])

  # #     #Initial population structure pop[study_area=k,sex=i,year=t,age=a]
  # #     llpop_inf[k, 1, a, 1] ~ dnorm(f_logpop_inf[k, a], tau_pop[2])#tau_pop[1]
  # #     pop_inf[k, 1, a, 1] <- exp(llpop_inf[k, 1, a, 1])
  # #   }
  # #   for (a in 1:n_agem) {
  # #       ### Initial population structure pop
  # #       ### [study_area,sex,age,period(year)] for susceptible deer
  # #       llpop_sus[k, 2, a, 1] ~ dnorm(m_logpop_sus[k, a], tau_pop[1])#tau_pop[2]
  # #       pop_sus[k, 2, a, 1] <- exp(llpop_sus[k, 2, a, 1])

  # #       #Initial population structure pop for infected deer
  # #       llpop_inf[k, 2, a, 1] ~ dnorm(m_logpop_inf[k, a], tau_pop[2])#tau_pop[2]
  # #       pop_inf[k, 2, a, 1] <- exp(llpop_inf[k, 2, a, 1])
  # #   }
  # # }

  # ############################
  # ####Reporting Rates
  # ############################

  # report_overall ~ dbeta(report_hyp_all[1], report_hyp_all[2])
  # for(t in 1:22){#1992-2014
  #   report[t]  <- report_overall
  # }
  # for(t in 23:27){ #2015-2020
  #   report[t] ~ dbeta(report_hyp_y[t - 22, 1], report_hyp_y[t - 22, 2])
  # }
  # report[28]  <- report_overall #2021

  # ############################
  # #### Fecundity
  # ############################

  # mu_fec ~ dnorm(0, 1)
  # fec_prec_eps ~ dgamma(1, 1)

  # #Observations of fawns & does overall from the 3 counties
  # for(t in 1:n_year_fec_early) {
  #   fec_epsilon[t] ~ dnorm(0, fec_prec_eps)
  #   fec[t] <- exp(mu_fec + fec_epsilon[t])
  #   Nfawn[t] ~ dpois(fec[t] * Ndoe[t])
  # }

  # #for 2017:2021
  # for(t in (n_year_fec_early + 1):n_year){
  #   fec[t] ~ dgamma(obs_ct_fd_alpha[t], obs_ct_fd_beta[t])
  # }

  ###################################################
  #### Overall Survival Susceptibles
  ###################################################

  sn_sus[1:n_sex, 1:n_agef, 1:n_year] <- calc_surv_aah(
      nT_age = nT_age_surv,
      nT_period_overall = nT_period_overall,
      nT_age_short_f = nT_age_short_f,
      nT_age_short_m = nT_age_short_m,
      nT_age_surv_aah_f = nT_age_surv_aah_f,
      nT_age_surv_aah_m = nT_age_surv_aah_m,
      beta0 = beta0_survival_sus,
      beta_male = beta_male,
      age_effect = age_effect_survival[1:nT_age_surv],
      period_effect = period_effect_survival[(nT_period_prestudy_ext + 1):
                                              nT_period_overall_ext],
      yr_start = yr_start[1:n_year],
      yr_end = yr_end[1:n_year],
      n_year = n_year,
      n_agef = n_agef,
      n_agem = n_agem)

  ###################################################
  #### Overall Survival CWD Infected
  ###################################################

  sn_inf[1:n_sex, 1:n_agef, 1:n_year] <- calc_surv_aah(
      nT_age = nT_age_surv,
      nT_period_overall = nT_period_overall,
      nT_age_short_f = nT_age_short_f,
      nT_age_short_m = nT_age_short_m,
      nT_age_surv_aah_f = nT_age_surv_aah_f,
      nT_age_surv_aah_m = nT_age_surv_aah_m,
      beta0 = beta0_survival_inf,
      beta_male = beta_male,
      age_effect = age_effect_survival[1:nT_age_surv],
      period_effect = period_effect_survival[(nT_period_prestudy_ext + 1):
                                              nT_period_overall_ext],
      yr_start = yr_start[1:n_year],
      yr_end = yr_end[1:n_year],
      n_year = n_year,
      n_agef = n_agef,
      n_agem = n_agem)

  ###################################################
  #### Hunting Survival Susceptibles
  ###################################################

  sh_sus[1:n_sex, 1:n_agef, 1:n_year] <- calc_surv_harvest(
      nT_age = nT_age_surv,
      nT_period_overall = nT_period_overall,
      nT_age_short_f = nT_age_short_f,
      nT_age_short_m = nT_age_short_m,
      nT_age_surv_aah_f = nT_age_surv_aah_f,
      nT_age_surv_aah_m = nT_age_surv_aah_m,
      beta0 = beta0_survival_sus,
      beta_male = beta_male,
      age_effect = age_effect_survival[1:nT_age_surv],
      period_effect = period_effect_survival[(nT_period_prestudy_ext + 1):
                                             nT_period_overall_ext],
      n_sex = n_sex,
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

  sh_inf[1:n_sex, 1:n_agef, 1:n_year] <- calc_surv_harvest(
      nT_age = nT_age_surv,
      nT_period_overall = nT_period_overall,
      nT_age_short_f = nT_age_short_f,
      nT_age_short_m = nT_age_short_m,
      nT_age_surv_aah_f = nT_age_surv_aah_f,
      nT_age_surv_aah_m = nT_age_surv_aah_m,
      beta0 = beta0_survival_inf,
      beta_male = beta_male,
      age_effect = age_effect_survival[1:nT_age_surv],
      period_effect = period_effect_survival[(nT_period_prestudy_ext + 1):
                                              nT_period_overall_ext],
      n_sex = n_sex,
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


  ###########################################################
  #### Annual Probability of Infection based on FOI hazards
  ###########################################################

  psi[1:n_study_area, 1:n_sex, 1:n_agef, 1:n_year] <-
      calc_infect_prob(age_lookup_f = age_lookup_f[1:nT_age_surv],
            age_lookup_m = age_lookup_m[1:nT_age_surv],
            n_agef = n_agef,
            n_agem = n_agem,
            yr_start = yr_start[1:n_year],
            yr_end = yr_end[1:n_year],
            f_age = f_age_foi[1:n_ageclassf],
            m_age = m_age_foi[1:n_ageclassm],
            f_period = f_period_foi[1:n_year],
            m_period = m_period_foi[1:n_year],
            nT_period_overall = nT_period_overall,
            period_lookup_foi_study = period_lookup_foi_study[1:nT_period_overall],
            n_year = n_year,
            n_sex = n_sex,
            n_study_area = n_study_area,
            space = space[n_study_area],
            nT_age_surv_aah_f = nT_age_surv_aah_f,
            nT_age_surv_aah_m = nT_age_surv_aah_m
            )

  ##################################################################
  ### Probability of Infection from birth pulse to end of harvest
  ### based on FOI hazards
  ##################################################################

  psi_hat[1:n_study_area, 1:n_sex, 1:n_agef, 1:n_year] <-
      calc_infect_prob_hunt(age_lookup_f = age_lookup_f[1:nT_age_surv],
          age_lookup_m = age_lookup_m[1:nT_age_surv],
          n_agef = n_agef,
          n_agem = n_agem,
          yr_start = yr_start[1:n_year],
          yr_end = yr_end[1:n_year],
          ng_end = ng_end[1:n_year],
          f_age = f_age_foi[1:n_ageclassf],
          m_age = m_age_foi[1:n_ageclassm],
          f_period = f_period_foi[1:n_year],
          m_period = m_period_foi[1:n_year],
          nT_period_overall = nT_period_overall,
          period_lookup_foi_study = period_lookup_foi_study[1:nT_period_overall],
          n_year = n_year,
          n_sex = n_sex,
          n_study_area = n_study_area,
          space = space[n_study_area],
          nT_age_surv_aah_f = nT_age_surv_aah_f,
          nT_age_surv_aah_m = nT_age_surv_aah_m,
          fudge_factor = .5
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

  #   #####################################
  #   ### version without projection matrix
  #   #####################################

  # for (k in 1:n_study_area) {
  #     for (t in 2:n_year) {

  #       #Female: set projection into population model matrix
  #       for (a in 2:(n_agef - 1)) {
  #         pop_sus[k, 1, a, t] <-  pop_sus[k, 1, a - 1, t - 1] *
  #                                     sn_sus[1, a - 1, t - 1] *
  #                                     (1 - psi[k, 1, a - 1, t - 1])
  #       }
  #       pop_sus[k, 1, n_agef, t] <- pop_sus[k, 1, n_agef - 1, t - 1] *
  #                                           sn_sus[1, n_agef - 1, t - 1] *
  #                                           (1 - psi[k, 1, n_agef - 1, t - 1]) +
  #                                           pop_sus[k, 1,  n_agef, t - 1] *
  #                                           sn_sus[1, n_agef, t - 1] *
  #                                           (1 - psi[k, 1, n_agef, t - 1])

  #       #Female: fawn class = total #females * unisex fawns per female/2
  #       pop_sus[k, 1, 1, t] <- (sum(pop_sus[k, 1, 2:n_agef, t]) +
  #                               sum(pop_inf[k, 1, 2:n_agef, t])) *
  #                               fec[t] * .5

  #       ##########
  #       ### Males
  #       ##########

  #       #Male: set projection into population model matrix
  #       for (a in 2:(n_agem - 1)) {
  #         pop_sus[k, 2, a, t] <- pop_sus[k, 2, a - 1, t - 1] *
  #                                sn_sus[2, a - 1, t - 1] *
  #                                (1 - psi[k, 2, a - 1, t - 1])
  #       }
  #       #Male: accumulating age class
  #       pop_sus_[k, 2, n_agem, t] <- pop_sus[k, 2, n_agem - 1, t - 1] *
  #                                    sn_sus[2, n_agem - 1, t - 1] *
  #                                    (1 - psi[k, 2, n_agem - 1, t - 1]) +
  #                                    pop_sus[k, 2,  n_agem, t - 1] *
  #                                    sn_sus[2, n_agem, t - 1] *
  #                                    (1 - psi[k, 2, n_agem, t - 1])

  #       # Male: fawn class = # female fawns
  #       pop_sus[k, 2, 1, t] <- pop_sus[k, 1, 1, t]


  #       ###################################################
  #       ### Infected/Infectious
  #       ###################################################

  #       ###########
  #       ### Females
  #       ###########

  #       #Female: set projection into population model matrix
  #       for (a in 2:(n_agef - 1)) {
  #         pop_inf[k, 1, a, t] <- pop_inf[k, 1, a - 1, t - 1] *
  #                                     sn_inf[1, a - 1, t - 1] +
  #                                     pop_sus[k, 1, a - 1, t - 1] *
  #                                     sn_sus[1, a - 1, t - 1] *
  #                                     psi[k, 1, a - 1, t - 1]
  #       }
  #       ##Female: accumulating age = 9.5+ years
  #       pop_inf[k, 1, n_agef, t] <- pop_inf[k, 1, n_agef - 1, t - 1] *
  #                                   sn_inf[1, n_agef - 1, t - 1] +
  #                                   pop_inf[k, 1, n_agef, t - 1] *
  #                                   sn_inf[1, n_agef, t - 1] +
  #                                   pop_sus[k, 1, n_agef - 1, t - 1] *
  #                                   sn_sus[1, n_agef - 1, t - 1] *
  #                                   psi[k, 1, n_agef - 1, t - 1] +
  #                                   pop_sus[k, 1,  n_agef, t - 1] *
  #                                   sn_sus[1, n_agef, t - 1] *
  #                                   psi[k, 1, n_agef, t - 1]

  #       ##Female: fawn class
  #       ##there are no infected fawns at birth
  #       pop_inf[k, 1, 1, t] <- 0

  #       ##########
  #       ### Males
  #       ##########

  #       ### Male: project forward anually
  #       for (a in 2:(n_agem - 1)) {
  #           pop_inf[k, 2, a, t] <- pop_inf[k, 2, a - 1, t - 1] *
  #                                         sn_inf[2, a - 1, t - 1] +
  #                                         pop_sus[k, 2, a - 1, t - 1] *
  #                                         sn_sus[2, a - 1, t - 1] *
  #                                         psi[k, 2, a - 1, t - 1]
  #       }

  #       ### Male: accumulating age class = 6.5+
  #       pop_inf[k, 2, n_agem, t] <- pop_inf[k, 2, n_agem - 1, t - 1] *
  #                                            sn_inf[2, n_agem - 1, t - 1] +
  #                                            pop_inf[k, 2, n_agem, t - 1] *
  #                                            sn_inf[2, n_agem, t - 1] +
  #                                            pop_sus[k, 2, n_agem - 1, t - 1] *
  #                                            sn_sus[2, n_agem - 1, t - 1] *
  #                                            psi[k, 2, n_agem - 1, t - 1] +
  #                                            pop_sus[k, 2,  n_agem, t - 1] *
  #                                            sn_sus[2, n_agem, t - 1] *
  #                                            psi[k, 2, n_agem, t - 1] 

  #       #Male: fawn class
  #       #there are no infected fawns at birth
  #       pop_inf[k, 2, 1, t] <- 0

  #       }
  #       }

    #####################################
    ### version with projection matrix
    #####################################

  #   for (k in 1:n_study_area) {

  #     ##################
  #     ### Susceptible
  #     ##################

  #     for (t in 2:n_year) {
        

  #       ##############################################################
  #       ##############################################################
  #       ### AAH population model version with a projection matrix
  #       ##############################################################
  #       ##############################################################

  #       #Female: project forward annually
  #       for (a in 1:(n_agef - 2)) {
  #         pop_sus_proj[k, 1, a, t] <- pop_sus[k, 1, a, t - 1] *
  #                                     sn_sus[1, a, t - 1] *
  #                                     (1 - psi[k, 1, a, t - 1])
  #       }

  #       #female accumulating age class
  #       pop_sus_proj[k, 1, n_agef - 1, t] <- pop_sus[k, 1, n_agef - 1, t - 1] *
  #                                           sn_sus[1, n_agef - 1, t - 1] *
  #                                           (1 - psi[k, 1, n_agef - 1, t - 1]) +
  #                                           pop_sus[k, 1,  n_agef, t - 1] *
  #                                           sn_sus[1, n_agef, t - 1] *
  #                                           (1 - psi[k, 1, n_agef, t - 1])

  #       #Female: set projection into population model matrix
  #       for (a in 2:n_agef) {
  #         pop_sus[k, 1, a, t] <- pop_sus_proj[k, 1, (a - 1), t]
  #       }

  #       #Female: fawn class = total #females * unisex fawns per female/2
  #       pop_sus[k, 1, 1, t] <- (sum(pop_sus_proj[k, 1, 1:(n_agef - 1), t]) +
  #                               sum(pop_inf_proj[k, 1, 1:(n_agef - 1), t])) *
  #                               fec[t] * .5
  #       ##########
  #       ### Males
  #       ##########

  #       #Male: project forward anually
  #       for (a in 1:(n_agem - 2)) {
  #         pop_sus_proj[k, 2, a, t] <- pop_sus[k, 2, a, t - 1] *
  #                                     sn_sus[2, a, t - 1] *
  #                                     (1 - psi[k, 2, a, t - 1])
  #       }

  #       #Male: accumulating age class
  #       pop_sus_proj[k, 2, n_agem - 1, t] <- pop_sus[k, 2, n_agem - 1, t - 1] *
  #                                           sn_sus[2, n_agem - 1, t - 1] *
  #                                           (1 - psi[k, 2, n_agem - 1, t - 1]) +
  #                                           pop_sus[k, 2,  n_agem, t - 1] *
  #                                           sn_sus[2, n_agem, t - 1] *
  #                                           (1 - psi[k, 2, n_agem, t - 1])

  #       #Male: set projection into population model matrix
  #       for (a in 2:n_agem) {
  #         pop_sus[k, 2, a, t] <- pop_sus_proj[k, 2, (a - 1), t]
  #       }

  #       # Male: fawn class = total #females * unisex fawns per female/2
  #       pop_sus[k, 2, 1, t] <- pop_sus[k, 1, 1, t]

  #       ###################################################
  #       ### Infected/Infectious
  #       ###################################################

  #       ###########
  #       ### Females
  #       ###########

  #       #Female: project forward anually
  #       for (a in 1:(n_agef - 2)) {
  #         pop_inf_proj[k, 1, a, t] <- pop_inf[k, 1, a, t - 1] *
  #                                     sn_inf[1, a, t - 1] +
  #                                     pop_sus[k, 1, a, t - 1] *
  #                                     sn_sus[1, a, t - 1] *
  #                                     psi[k, 1, a, t - 1]
  #       }
  #       ##Female: accumulating age = 9.5+ years
  #       pop_inf_proj[k, 1, n_agef - 1, t] <- pop_inf[k, 1, n_agef - 1, t - 1] *
  #                                            sn_inf[1, n_agef - 1, t - 1] +
  #                                            pop_inf[k, 1, n_agef, t - 1] *
  #                                            sn_inf[1, n_agef, t - 1] +
  #                                            pop_sus[k, 1, n_agef - 1, t - 1] *
  #                                            sn_sus[1, n_agef - 1, t - 1] *
  #                                            psi[k, 1, n_agef - 1, t - 1] +
  #                                            pop_sus[k, 1,  n_agef, t - 1] *
  #                                            sn_sus[1, n_agef, t - 1] *
  #                                            psi[k, 1, n_agef, t - 1]


  #       #Female: set projection into population model matrix
  #       for (a in 2:n_agef) {
  #         pop_inf[k, 1, a, t] <- pop_inf_proj[k, 1, (a - 1), t]
  #       }
  #       ##Female: fawn class
  #       ##there are no infected fawns at birth
  #       pop_inf[k, 1, 1, t] <- 0

  #       ##########
  #       ### Males
  #       ##########

  #       ### Male: project forward anually
  #       for (a in 1:(n_agem - 2)) {
  #           pop_inf_proj[k, 2, a, t] <- pop_inf[k, 2, a, t - 1] *
  #                                         sn_inf[2, a, t - 1] +
  #                                         pop_sus[k, 2, a, t - 1] *
  #                                         sn_sus[2, a, t - 1] *
  #                                         psi[k, 2, a, t - 1]
  #       }

  #       ### Male: accumulating age class
  #       pop_inf_proj[k, 2, n_agem - 1, t] <- pop_inf[k, 2, n_agem - 1, t - 1] *
  #                                            sn_inf[2, n_agem - 1, t - 1] +
  #                                            pop_inf[k, 2, n_agem, t - 1] *
  #                                            sn_inf[2, n_agem, t - 1] +
  #                                            pop_sus[k, 2, n_agem - 1, t - 1] *
  #                                            sn_sus[2, n_agem - 1, t - 1] *
  #                                            psi[k, 2, n_agem - 1, t - 1] +
  #                                            pop_sus[k, 2,  n_agem, t - 1] *
  #                                            sn_sus[2, n_agem, t - 1] *
  #                                            psi[k, 2, n_agem, t - 1] 

  #       ### Male: set projection into population model matrix
  #       for (a in 2:n_agem) {
  #         pop_inf[k, 2, a, t] <- pop_inf_proj[k, 2, (a - 1), t]
  #       }

  #       #Male: fawn class
  #       #there are no infected fawns at birth
  #       pop_inf[k, 2, 1, t] <- 0

  #     }#end t
  #   }#end study_area

  #   ######################################################################
  #   ### Observation Model
  #   ######################################################################

  #   tau_obs ~ dgamma(1, 1)

  #   # for (i in 1:n_sex) {
  #   #   tau_obs[i] ~ dgamma(1, 1)
  #   # }#end i

  #   for(k in 1:n_study_area) {

  #     # for (i in 1:n_sex) {
  #     #   tau_obs[k, i] ~ dgamma(1, 1)
  #     # }#end i

  #     for (t in 1:n_year) {
  #       for (a in 1:n_agef) {
  #         harv_pop[k, 1, a, t] <- (pop_inf[k, 1, a, t] *
  #                                 (1 - sh_inf[1, a, t]) +
  #                                 pop_sus[k, 1, a, t] *
  #                                 (1 - sh_sus[1, a, t]) *
  #                                 (1 - psi_hat[k, 1, a, t]) +
  #                                 pop_sus[k, 1, a, t] *
  #                                 (1 - sh_inf[1, a, t]) *
  #                                 psi_hat[k, 1, a, t]) *
  #                                 report[t]
  #       }
  #       for (a in 1:n_agem) {
  #         harv_pop[k, 2, a, t] <- (pop_inf[k, 2, a, t] *
  #                                 (1 - sh_inf[2, a, t]) +
  #                                 pop_sus[k, 2, a, t] *
  #                                 (1 - sh_sus[2, a, t]) *
  #                                 (1 - psi_hat[k, 2, a, t]) +
  #                                 pop_sus[k, 2, a, t] *
  #                                 (1 - sh_inf[2, a, t]) *
  #                                 psi_hat[k, 2, a, t]) *
  #                                 report[t]
  #       }

  #       #Total Antlerless Harvest
  #       #adding in male fawns
  #       mu_obs[k, 1, t] <- (sum(harv_pop[k, 1, 1:n_agef, t]) +
  #                           harv_pop[k, 2, 1, t]) # eab_antlerless[t] *

  #       #Total Antlered Harvest
  #       #excludes male fawns 
  #       mu_obs[k, 2, t] <- sum(harv_pop[k, 2, 2:n_agem, t]) #eab_antlered[t] *

  #       ###################################
  #       #Likelihood for overall total
  #       ###################################

  #       for (j in 1:n_sex) {
  #         lobs[k, j, t] ~ dnorm(log(mu_obs[k, j, t]), tau_obs)# tau_obs[j]
  #       }#end i

  #     #   ###################################
  #     #   #Likelihood for overall total
  #     #   ###################################

  #       # ### parameters for likelihood harvest data by antlerless group
  #       # ### proportion of each age class
  #       # ### Antlerless deer
  #       # p_less[k, 1, t] <- harv_pop[k, 1, 1, t] / mu_obs[k, 1, t]#f
  #       # p_less[k, 2, t] <- harv_pop[k, 1, 2, t] / mu_obs[k, 1, t]#1
  #       # p_less[k, 3, t] <- harv_pop[k, 1, 3, t] / mu_obs[k, 1, t]#2
  #       # p_less[k, 4, t] <- harv_pop[k, 1, 4, t] / mu_obs[k, 1, t]#3
  #       # p_less[k, 5, t] <- sum(harv_pop[k, 1, 5:6, t]) / mu_obs[k, 1, t]#4-5
  #       # p_less[k, 6, t] <- sum(harv_pop[k, 1, 7:9, t]) / mu_obs[k, 1, t]#6-8
  #       # p_less[k, 7, t] <- harv_pop[k, 1, 10, t] / mu_obs[k, 1, t]#9+
  #       # p_less[k, 8, t] <- 1 - sum(p_less[k, 1:7, t])#male f, antlerless

  #       # ### harvest data bt antlered group
  #       # p_ant[k, 1, t] <- harv_pop[k, 2, 2, t] / mu_obs[k, 2, t]#1
  #       # p_ant[k, 2, t] <- harv_pop[k, 2, 3, t] / mu_obs[k, 2, t]#2
  #       # p_ant[k, 3, t] <- harv_pop[k, 2, 4, t] / mu_obs[k, 2, t]#3
  #       # p_ant[k, 4, t] <- sum(harv_pop[k, 2, 5:6, t]) / mu_obs[k, 2, t]#4-5
  #       # p_ant[k, 5, t] <- 1 - sum(p_ant[k, 1:4, t]) #6+

  #     # }# end t

  #     # for(t in 1:n_year){
  #     #     #antlerless, male fawns
  #     #     Cage_less[k, 1:(n_ageclassf + 1),t] ~ dmulti(prob = p_less[k,1:(n_ageclassf + 1),t],
  #     #                           size = sizeCage_f[k, t])
  #     #     #antlered
  #     #     Cage_ant[k, 1:(n_ageclassm - 1),t] ~ dmulti(prob = p_ant[k, 1:(n_ageclassm - 1),t],
  #     #               size = sizeCage_m[k, t])
  #     #   }
  #     }# end t
  # }

})#end model statement
