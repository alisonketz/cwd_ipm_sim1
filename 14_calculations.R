#######################################################################
###
### Function to calculate Annual survival probability based on 
### age effects and period effects
###
#######################################################################

calc_surv_aah <- nimble::nimbleFunction(
    run = function(
        ### argument type declarations
        nT_age = double(0),
        nT_period_overall = double(0),
        nT_age_short_f = double(0),
        nT_age_short_m = double(0),
        nT_age_surv_aah_f = double(0),
        nT_age_surv_aah_m = double(0),
        beta0 = double(0),
        beta_male = double(0),
        age_effect = double(1),
        period_effect = double(1),
        yr_start = double(1),
        yr_end = double(1),
        n_year = double(0),
        n_agef = double(0),
        n_agem = double(0)
        ) {

    ###################################################################
    ###
    ### General Survival Surfaces for Susceptible/Infected Individuals
    ###
    ###################################################################

    mu_old_age_effect_f <- mean(age_effect[(nT_age_short_f + 1):nT_age])
    mu_old_age_effect_m <- mean(age_effect[(nT_age_short_m + 1):nT_age])

	############################################
	# Calculate hazards
	############################################

    UCH <- nimArray(NA, c(2, nT_age_surv_aah_f, nT_period_overall))
    s_aah <- nimArray(NA, c(2, n_agef, n_year))

    ### Females
    for(j in 1:nT_period_overall) {
        for(i in 1:nT_age_short_f) {
            UCH[1, i, j] <- exp(beta0 +
                                age_effect[i] +
                                period_effect[j])
        }
        for(i in (nT_age_short_f + 1):(nT_age_surv_aah_f)) {
            UCH[1, i, j] <- exp(beta0 +
                                mu_old_age_effect_f +
                                period_effect[j])
        }
    }

    ### Males
    for(j in 1:nT_period_overall) {
        for(i in 1:nT_age_short_m) {
            UCH[2, i, j] <- exp(beta0 +
                                beta_male +
                                age_effect[i] +
                                period_effect[j])
        }
        for(i in (nT_age_short_m + 1):(nT_age_surv_aah_m)) {
            UCH[2, i, j] <- exp(beta0 +
                                beta_male +
                                mu_old_age_effect_m +
                                period_effect[j])
        }
    }
    ############################################
    # calculate survival from cummulative haz
    ############################################
    for (t in 1:n_year) {
        for (a in 1:n_agef) {
            s_aah[1, a, t] <- exp(-sum(diag(UCH[1,
                               yr_start[a]:(yr_start[a] + 52),
                               yr_start[t]:(yr_start[t] + 52)])))
        }
        for(a in 1:n_agem) {
            s_aah[2, a, t] <- exp(-sum(diag(UCH[2,
                               yr_start[a]:(yr_start[a] + 52),
                               yr_start[t]:(yr_start[t] + 52)])))
        }
    }
    # for (t in 1:n_year) {
    #     for (a in 1:n_agef) {
    #         s_aah[1, a, t] <- exp(-sum(diag(UCH[1,
    #                            yr_start[a]:(yr_start[a] + length(yr_start[t]:yr_end[t]) - 1),
    #                            yr_start[t]:yr_end[t]])))
    #     }
    #     for(a in 1:n_agem) {
    #         s_aah[2, a, t] <- exp(-sum(diag(UCH[2,
    #                            yr_start[a]:(yr_start[a] + length(yr_start[t]:yr_end[t]) - 1),
    #                            yr_start[t]:(yr_end[t])])))
    #     }
    # }

  returnType(double(3))
  return(s_aah[1:2, 1:n_agef, 1:n_year])
})

Ccalc_surv_aah <- compileNimble(calc_surv_aah)

assign("calc_surv_aah", calc_surv_aah, envir = .GlobalEnv)

# starttime <- Sys.time()
# sn_sus <- Ccalc_surv_aah(
# 	nT_age = nT_age_surv,
#     nT_period_overall = nT_period_overall,
#     nT_age_short_f = nT_age_short_f,
#     nT_age_short_m = nT_age_short_m,
#     nT_age_surv_aah_f = nT_age_surv_aah_f,
#     nT_age_surv_aah_m = nT_age_surv_aah_m,
#     beta0 = beta0_survival_sus+2.3,
#     beta_male = beta_male,
#     age_effect = age_effect_survival_test,        # length = 962
#     period_effect = period_effect_survival_test[(nT_period_prestudy_ext + 1):(nT_period_overall_ext)], 
# 	yr_start = d_fit_season$yr_start,
#     yr_end = d_fit_season$yr_end,
#     n_year = n_year,
#     n_agef = n_agef,
#     n_agem = n_agem)
# (endtime1 <- Sys.time() - starttime)
# sn_sus[1,3,]
# sn_sus[2,,]

# starttime <- Sys.time()
# sn_inf <- Ccalc_surv_aah(
#         nT_age = nT_age_surv,
#         nT_period_overall = nT_period_overall,
#         nT_age_short_f = nT_age_short_f,
#         nT_age_short_m = nT_age_short_m,
#         nT_age_surv_aah_f = nT_age_surv_aah_f,
#         nT_age_surv_aah_m = nT_age_surv_aah_m,
#         beta0 = beta0_survival_inf + 3,
#         beta_male = beta_male,
#         age_effect = age_effect_survival_test,
#         period_effect = period_effect_survival_test[(nT_period_prestudy_ext + 1):(nT_period_overall_ext)],
#         yr_end = d_fit_season$yr_end,
#         yr_start = d_fit_season$yr_start,
#         n_year = n_year,
#         n_agef = n_agef,
#         n_agem = n_agem)
# (endtime2 <- Sys.time() - starttime)
# round(sn_inf[1,3,],2)
# round(sn_inf[2,3,],2)

# plot(1:28,sn_inf[1,3,])

# plot(1:28,sn_inf[2,3,])
# sn_inf[2,3,]

# starttime <- Sys.time()
# sn_sus <- calc_surv_aah(
# 	nT_age = nT_age_surv,
#     nT_period_overall = nT_period_overall,
#     nT_age_short_f = nT_age_short_f,
#     nT_age_short_m = nT_age_short_m,
#     nT_age_surv_aah_f = nT_age_surv_aah_f,
#     nT_age_surv_aah_m = nT_age_surv_aah_m,
#     beta0 = -10.5,
#     beta_male = .5,
#     age_effect = rep(0,nT_age_surv),        # length = 962
#     period_effect = rep(0,length((nT_period_prestudy_ext + 1):(nT_period_overall_ext))),
# 	yr_start = d_fit_season$yr_start,
#     yr_end = d_fit_season$yr_end,
#     n_year = n_year,
#     n_agef = n_agef,
#     n_agem = n_agem)
# sn_sus[1,2,]
# (endtime1 <- Sys.time() - starttime)

#######################################################################
###
### Function to calculate annual survival based on age effects and 
### period effects
###
#######################################################################

calc_surv_harvest <- nimble::nimbleFunction(
    run = function(
        ### argument type declarations
        nT_age = double(0),
        nT_age_short_f = double(0),
        nT_age_short_m = double(0),
        nT_age_surv_aah_f = double(0),
        nT_age_surv_aah_m = double(0),
        nT_period_overall = double(0),
        beta0 = double(0),
        beta_male = double(0),
        age_effect = double(1),
        period_effect = double(1),
        n_sex = double(0),
        n_year = double(0),
        n_agef = double(0),
        n_agem = double(0),
        ng_start = double(1),
        gun_start = double(1),
        gun_end = double(1),
        ng_end = double(1),
        yr_start = double(1),
        yr_end = double(1),
        p_nogun_f = double(0),
        p_nogun_m = double(0),
        p_gun_f = double(0),
        p_gun_m = double(0)
        ) {

    ###################################################################
    ###
    ### General Survival Surfaces for Susceptible/Infected Individuals
    ###
    ###################################################################

    mu_old_age_effect_f <- mean(age_effect[(nT_age_short_f + 1):nT_age])
    mu_old_age_effect_m <- mean(age_effect[(nT_age_short_m + 1):nT_age])

    ############################################
    ### Initialize hazard array
    ############################################

    UCH <- nimArray(NA, c(n_sex, nT_age_surv_aah_f, nT_period_overall))
    UCH_hunt <- nimArray(NA, c(n_sex, nT_age_surv_aah_f, nT_period_overall))
    s_hunt <- nimArray(NA, c(n_sex, n_agef, n_year))

    ############################################
    ### Calculate hazards
    ############################################

    # for(j in 1:nT_period_overall) {
    #     for(i in 1:nT_age_short_m) {
    #         UCH[2, i, j] <- exp(beta0 +
    #                             beta_male +
    #                             age_effect[i] +
    #                             period_effect[j])
    #     }
    #     for(i in (nT_age_short_m + 1):(nT_age_surv_aah_m)) {
    #         UCH[2, i, j] <- exp(beta0 +
    #                             beta_male +
    #                             mu_old_age_effect_m +
    #                             period_effect[j])
    #     }
    # }
    ### Females
    for(j in 1:nT_period_overall) {
        for(i in 1:nT_age_short_f) {
            UCH[1, i, j] <- exp(beta0 +
                                age_effect[i] +
                                period_effect[j])
        }
        for(i in (nT_age_short_f + 1):(nT_age_surv_aah_f)) {
            UCH[1, i, j] <- exp(beta0 +
                                mu_old_age_effect_f +
                                period_effect[j])
        }
    }
    ### Males
    for(j in 1:nT_period_overall) {
        for(i in 1:nT_age_short_m) {
            UCH[2, i, j] <- exp(beta0 +
                                beta_male +
                                age_effect[i] +
                                period_effect[j])
        }
        for(i in (nT_age_short_m + 1):(nT_age_surv_aah_m)) {
            UCH[2, i, j] <- exp(beta0 +
                                beta_male +
                                mu_old_age_effect_m +
                                period_effect[j])
        }
    }
	############################################
	# adjust hazards to remove harvest hazards
	############################################
	UCH_hunt[1, 1:nT_age_surv_aah_f, 1:nT_period_overall] <- UCH[1, 1:nT_age_surv_aah_f, 1:nT_period_overall]
	UCH_hunt[2, 1:nT_age_surv_aah_f, 1:nT_period_overall] <- UCH[2, 1:nT_age_surv_aah_f, 1:nT_period_overall]

    for(i in 1:n_year){
    #    for(j in yr_start[i]:(ng_start[i] - 1)) {
    #        UCH_hunt[1, 1:nT_age_surv_aah_f, j] <- UCH[1, 1:nT_age_surv_aah_f, j]
    #        UCH_hunt[2, 1:nT_age_surv_aah_m, j] <- UCH[2, 1:nT_age_surv_aah_m, j]
    #    }
    #    for(j in (ng_end[i] + 1):(yr_end[i])){
    #        UCH_hunt[1, 1:nT_age_surv_aah_f, j] <- UCH[1, 1:nT_age_surv_aah_f, j]
    #        UCH_hunt[2, 1:nT_age_surv_aah_m, j] <- UCH[2, 1:nT_age_surv_aah_m, j]
    #    }
        for(j in ng_start[i]:(gun_start[i] - 1)){
            UCH_hunt[1, 1:nT_age_surv_aah_f, j] <- UCH[1, 1:nT_age_surv_aah_f, j] * p_nogun_f
            UCH_hunt[2, 1:nT_age_surv_aah_m, j] <- UCH[2, 1:nT_age_surv_aah_m, j] * p_nogun_m
        }
        for(j in gun_start[i]:(gun_end[i])){
            UCH_hunt[1, 1:nT_age_surv_aah_f, j] <- UCH[1, 1:nT_age_surv_aah_f, j] * p_gun_f
            UCH_hunt[2, 1:nT_age_surv_aah_m, j] <- UCH[2, 1:nT_age_surv_aah_m, j] * p_gun_m
        }
		if(gun_end[i] < ng_end[i]){
        for(j in (gun_end[i] + 1):(ng_end[i])){
            UCH_hunt[1, 1:nT_age_surv_aah_f, j] <- UCH[1, 1:nT_age_surv_aah_f, j] * p_nogun_f
            UCH_hunt[2, 1:nT_age_surv_aah_m, j] <- UCH[2, 1:nT_age_surv_aah_m, j] * p_nogun_m
       }}
    }
	
	############################################
	# calculate survival from cum haz
	############################################
    for (t in 1:n_year) {
        for (a in 1:n_agef) {
            s_hunt[1, a, t] <- exp(-sum(diag(UCH_hunt[1,
                                yr_start[a]:(yr_start[a] + length(yr_start[t]:ng_end[t]) - 1),
                                yr_start[t]:ng_end[t]])))
        }
        for(a in 1:n_agem) {
            s_hunt[2, a, t] <- exp(-sum(diag(UCH_hunt[2,
                                yr_start[a]:(yr_start[a] + length(yr_start[t]:ng_end[t]) - 1),
                                yr_start[t]:ng_end[t]])))
        }
    }
  returnType(double(3))
  return(s_hunt[1:n_sex,1:n_agef,1:n_year])
})

Ccalc_surv_harvest <- compileNimble(calc_surv_harvest)

assign("calc_surv_harvest", calc_surv_harvest, envir = .GlobalEnv)

# starttime <- Sys.time()
# sh_sus <- calc_surv_harvest(
#         nT_age = nT_age_surv,
#         nT_period_overall = nT_period_overall,
#         nT_age_short_f = nT_age_short_f,
#         nT_age_short_m = nT_age_short_m,
#         nT_age_surv_aah_f = nT_age_surv_aah_f,
#         nT_age_surv_aah_m = nT_age_surv_aah_m,
#         beta0 = beta0_survival_sus,
#         beta_male = beta_male,
#         age_effect = age_effect_survival_test,
#         period_effect = period_effect_survival_test[(nT_period_prestudy_ext + 1):(nT_period_overall_ext)], 
#         n_year = n_year,
#         n_agef = n_agef,
#         n_agem = n_agem,
#         ng_start = d_fit_season$ng_start,
#         gun_start = d_fit_season$gun_start,
#         gun_end = d_fit_season$gun_end,
#         ng_end = d_fit_season$ng_end,
#         yr_start = d_fit_season$yr_start,
#         yr_end = d_fit_season$yr_end,
#         p_nogun_f = .55,
#         p_nogun_m = .65,
#         p_gun_f = .8,
#         p_gun_m = .95
#         # p_nogun_f = p_ng_f,
#         # p_nogun_m = p_ng_m,
#         # p_gun_f = p_gun_f,
#         # p_gun_m = p_gun_m
#         )
# (endtime3 <- Sys.time() - starttime)
# sh_sus[1,4,]

#######################################################################
###
### Function to calculate probability of infection
### based on FOI age and period effects
### Weekly Version
###
#######################################################################

calc_infect_prob <- nimbleFunction(
  run = function(age_lookup_f = double(1),
                 age_lookup_m = double(1),
                 n_agef = double(0),
                 n_agem = double(0),
                 yr_start = double(1),
                 yr_end = double(1),
                 f_age = double(1),
                 m_age = double(1),
                 f_period = double(1),
                 m_period = double(1),
                 nT_period_overall = double(0),
                 period_lookup_foi_study = double(1),
                 n_year = double(0),
                 n_sex = double(0),
                 n_study_area = double(0),
                 space = double(0),
                 nT_age_surv_aah_f = double(0),
                 nT_age_surv_aah_m = double(0)) {

    gam <- nimArray(value = 0, c(n_study_area,
                                 n_sex,
                                 nT_age_surv_aah_f,
                                 nT_period_overall))
    p_inf <- nimArray(value = 0, c(n_study_area,
                                   n_sex,
                                   n_agef,
                                   n_year))

    for (t in 1:nT_period_overall) {
        for (i in 1:nT_age_surv_aah_f) {
            ### Female
            ### East
            gam[1, 1, i, t] <- exp(f_age[age_lookup_f[i]] +
                               f_period[period_lookup_foi_study[t]])
            ### West
            gam[2, 1, i, t] <- exp(f_age[age_lookup_f[i]] +
                               f_period[period_lookup_foi_study[t]] +
                               space)

        }
        for (i in 1:nT_age_surv_aah_m) {
            ### Male
            ### East
            gam[1, 2, i, t] <- exp(m_age[age_lookup_m[i]] +
                               m_period[period_lookup_foi_study[t]])
            ### West
            gam[2, 2, i, t] <- exp(m_age[age_lookup_m[i]] +
                               m_period[period_lookup_foi_study[t]] +
                               space)
        }
    }
   # infection probability all ages all years 
    for(k in 1:n_study_area) {
        for (t in 1:n_year) {
            for (a in 1:n_agef) {
                p_inf[k, 1, a, t] <- 
                    1 - exp(-sum(diag(gam[k,
                                          1,
                                          yr_start[a]:(yr_start[a] + 52),
                                          yr_start[t]:(yr_start[t] + 52)])))
            }
            for (a in 1:n_agem) {
                p_inf[k, 2, a, t] <- 
                    1 - exp(-sum(diag(gam[k,
                                          2,
                                          yr_start[a]:(yr_start[a] + 52),
                                          yr_start[t]:(yr_start[t] + 52)])))
            }
        }
    }

    returnType(double(4))
    return(p_inf[1:n_study_area, 1:n_sex, 1:n_agef, 1:n_year])
  })

Ccalc_infect_prob <- compileNimble(calc_infect_prob)

assign("calc_infect_prob", calc_infect_prob, envir = .GlobalEnv)

# yr_start = d_fit_season$yr_start
# yr_end = d_fit_season$yr_end
# f_age = f_age_foi-2
# m_age = m_age_foi-2
# f_period = f_period_foi-2
# m_period = m_period_foi-2
# space <- space[2]
# ##testing state.transition function as R function
# starttime <- Sys.time()
# psi <- Ccalc_infect_prob(age_lookup_f = age_lookup_f,
#                         age_lookup_m = age_lookup_m,
#                         n_agef = n_agef,
#                         n_agem = n_agem,
#                         yr_start = d_fit_season$yr_start,
#                         yr_end = d_fit_season$yr_end,
#                         f_age = f_age_foi,
#                         m_age = m_age_foi,
#                         f_period = f_period_foi,
#                         m_period = m_period_foi,
#                         nT_period_overall = nT_period_overall,
#                         period_lookup_foi_study = period_lookup_foi_study[1:nT_period_overall],
#                         n_year = n_year,
#                         n_sex = n_sex,
#                         n_study_area = n_study_area, 
#                         space = -.55,
#                         nT_age_surv_aah_f = nT_age_surv_aah_f,
#                         nT_age_surv_aah_m = nT_age_surv_aah_m
#                         )
# (endtime5 <- Sys.time() - starttime)
# psi[1,1,2,]
# psi[1,2,2,]

#######################################################################
###
### Function to calculate probability of infection 
### from birth to end of hunter harvest season
### based on FOI age and period effects
### Weekly Version
###
#######################################################################

calc_infect_prob_hunt <- nimbleFunction(
  run = function(age_lookup_f = double(1),
                 age_lookup_m = double(1),
                 n_agef = double(0),
                 n_agem = double(0),
                 yr_start = double(1),
                 yr_end = double(1),
                 ng_end = double(1),
                 f_age = double(1),
                 m_age = double(1),
                 f_period = double(1),
                 m_period = double(1),
                 nT_period_overall = double(0),
                 period_lookup_foi_study = double(1),
                 n_year = double(0),
                 n_sex = double(0),
                 n_study_area = double(0),
                 space = double(0),
                 nT_age_surv_aah_f = double(0),
                 nT_age_surv_aah_m = double(0),
                 fudge_factor = double(0)) {

    gam <- nimArray(value = 0, c(n_study_area,
                                 n_sex,
                                 nT_age_surv_aah_f,
                                 nT_period_overall))
    p_inf <- nimArray(value = 0, c(n_study_area,
                                   n_sex,
                                   n_agef,
                                   n_year))

    for (t in 1:nT_period_overall) {
        for (i in 1:nT_age_surv_aah_f) {
            ### Female
            ### East
            gam[1, 1, i, t] <- exp(f_age[age_lookup_f[i]] +
                               f_period[period_lookup_foi_study[t]])
            ### West
            gam[2, 1, i, t] <- exp(f_age[age_lookup_f[i]] +
                               f_period[period_lookup_foi_study[t]] +
                               space)

        }
        for (i in 1:nT_age_surv_aah_m) {
            ### Male
            ### East
            gam[1, 2, i, t] <- exp(m_age[age_lookup_m[i]] +
                               m_period[period_lookup_foi_study[t]])
            ### West
            gam[2, 2, i, t] <- exp(m_age[age_lookup_m[i]] +
                               m_period[period_lookup_foi_study[t]] +
                               space)
        }

    }

    #Probability of getting infected with CWD
    for(k in 1:n_study_area) {
        # infection probability all ages all years 
        for (t in 1:n_year) {
            for (a in 1:n_agef) {
                p_inf[k, 1, a, t] <- 
                    (1 - exp(-sum(diag(gam[k,
                                           1,
                                           yr_start[a]:(yr_start[a] + length(yr_start[t]:ng_end[t]) - 1),
                                           yr_start[t]:ng_end[t]])))) * fudge_factor
            }
            for (a in 1:n_agem) {
                p_inf[k, 2, a, t] <- 
                    (1 - exp(-sum(diag(gam[k,
                                           2,
                                           yr_start[a]:(yr_start[a] + length(yr_start[t]:ng_end[t]) - 1),
                                           yr_start[t]:ng_end[t]])))) * fudge_factor
            }
        }
    }
    returnType(double(4))
    return(p_inf[1:n_study_area, 1:n_sex, 1:n_agef, 1:n_year])
  })

Ccalc_infect_prob_hunt <- compileNimble(calc_infect_prob_hunt)

assign("calc_infect_prob_hunt", calc_infect_prob_hunt, envir = .GlobalEnv)

##testing state.transition function as R function
# starttime <- Sys.time()
# psi_hat <- calc_infect_prob_hunt(age_lookup_f = age_lookup_f,
#                         age_lookup_m = age_lookup_m,
#                         n_agef = n_agef,
#                         n_agem = n_agem,
#                         yr_start = d_fit_season$yr_start,
#                         yr_end = d_fit_season$yr_end,
#                         ng_end = d_fit_season$ng_end,
#                         f_age = f_age_foi-2,
#                         m_age = m_age_foi-2,
#                         f_period = f_period_foi-1.5,
#                         m_period = m_period_foi-1.5,
#                         nT_period_overall = nT_period_overall,
#                         period_lookup_foi_study = period_lookup_foi_study[1:nT_period_overall],
#                         n_year = n_year,
#                         n_sex = n_sex,
#                         n_study_area = n_study_area, 
#                         space = -.55,
#                         nT_age_surv_aah_f = nT_age_surv_aah_f,
#                         nT_age_surv_aah_m = nT_age_surv_aah_m,
#                         fudge_factor = .5
#                         )
# (endtime6 <- Sys.time() - starttime)
# psi_hat[1,1,,]
            

#######################################################################
### INCORRECT
### Function to calculate probability of infection
### based on FOI age and period effects
### Weekly Version
###
#######################################################################

# calc_infect_prob <- nimbleFunction(
#   run = function(age_lookup_f = double(1),
#                  age_lookup_m = double(1),
#                  n_agef = double(0),
#                  n_agem = double(0),
#                  yr_end = double(1),
#                  f_age = double(1),
#                  m_age = double(1),
#                  f_period = double(1),
#                  m_period = double(1),
#                  n_year = double(0),
#                  n_sex = double(0),
#                  n_study_area = double(0),
#                  space = double(0),
#                  intvl_step_yr = double(0)) {

#     np <- intvl_step_yr * n_agef
#     p <- nimArray(value = 0, c(n_study_area, 
#                                n_sex,
#                                np,
#                                n_year))
#     gam <- nimArray(value = 0, c(n_study_area,
#                                  n_sex,
#                                  np,
#                                  n_year))
#     p_inf <- nimArray(value = 0, c(n_study_area,
#                                    n_sex,
#                                    n_agef,
#                                    n_year))

#     for (t in 1:n_year) {
#         for (i in 1:np) {
#             ### Female
#             ### East
#             gam[1, 1, i, t] <- f_age[age_lookup_f[i]] + f_period[t]
#             p[1, 1, i, t] <- exp(-sum(exp(gam[1, 1, 1:i, t])))
#             ### West
#             gam[2, 1, i, t] <- f_age[age_lookup_f[i]] + f_period[t] + space
#             p[2, 1, i, t] <- exp(-sum(exp(gam[1, 1, 1:i, t])))

#             ### Male
#             ### East
#             gam[1, 2, i, t] <- m_age[age_lookup_m[i]] + m_period[t]
#             p[1, 2, i, t] <- exp(-sum(exp(gam[1, 2, 1:i, t])))
#             ### West
#             gam[2, 2, i, t] <- m_age[age_lookup_m[i]] + m_period[t] + space
#             p[2, 2, i, t] <- exp(-sum(exp(gam[2, 2, 1:i, t])))
#         }
#     }

#     for(k in 1:n_study_area) {
#         #fawn probability of infection all years
#         #both sexes
#         for(t in 1:n_year) {
#             p_inf[k, 1, 1, t] <- 1 - p[k, 1, yr_end[1], t]
#             p_inf[k, 2, 1, t] <- 1 - p[k, 2, yr_end[1], t]
#         }
#         #all non-fawn prob of infection first year
#         for (a in 2:n_agef) {
#             p_inf[k, 1, a, 1] <- 1 - p[k, 1, yr_end[a], 1]
#         }
#         for (a in 2:n_agem) {
#             p_inf[k, 2, a, 1] <- 1 - p[k, 2, yr_end[a], 1]
#         }

#         #non-fawn infection probability all years except first year
#         for (t in 2:n_year) {
#             for (a in 2:n_agef) {
#                 p_inf[k, 1, a, t] <- 
#                     1 - (p[k, 1, yr_end[a], t] /
#                          p[k, 1, yr_end[a - 1], t])
#             }
#             for (a in 2:n_agem) {
#                 p_inf[k, 2, a, t] <-
#                     1 - (p[k, 2, yr_end[a], t] /
#                          p[k, 1, yr_end[a - 1], t])
#             }
#         }
#     }

#     returnType(double(4))
#     return(p_inf[1:n_study_area, 1:n_sex, 1:n_agef, 1:n_year])
#   })

# Ccalc_infect_prob <- compileNimble(calc_infect_prob)

# assign("calc_infect_prob", calc_infect_prob, envir = .GlobalEnv)

# ##testing state.transition function as R function
# starttime <- Sys.time()
# psi <- calc_infect_prob(age_lookup_f = age_lookup_f,
#                         age_lookup_m = age_lookup_m,
#                         n_agef = n_agef,
#                         n_agem = n_agem,
#                         yr_end = d_fit_season$yr_end,
#                         f_age = f_age_foi-.5,
#                         m_age = m_age_foi-.5,
#                         f_period = f_period_foi-.5,
#                         m_period = m_period_foi-.5,
#                         n_year = n_year,
#                         n_sex = n_sex,
#                         n_study_area = n_study_area, 
#                         space = -.55,
#                         intvl_step_yr = intvl_step_yr_weekly
#                         )
# (endtime5 <- Sys.time() - starttime)
# psi[1,,,]


#######################################################################
###
### Function to set and center period effects across full study time
### where period effects for non-harvest period effects
### based on aah data were different but constant each year
###
#######################################################################

set_period_effects_constant <- nimble::nimbleFunction(
    run = function(
        ### argument type declarations
        n_year_precollar = double(0),
        n_year_precollar_ext = double(0),
        n_year_prestudy_ext = double(0),
        nT_period_precollar_ext = double(0),
        nT_period_precollar = double(0),
        nT_period_collar = double(0),
        nT_period_overall_ext = double(0),
        nT_period_prestudy_ext = double(0),
        yr_start = double(1),
        yr_end = double(1),
        period_effect_surv = double(1),
        period_annual_survival = double(1)) {

  period_effect_survival_temp <- nimNumeric(nT_period_overall_ext)
  period_effect_survival <- nimNumeric(nT_period_overall_ext)

  for(i in 1:nT_period_prestudy_ext) {
    period_effect_survival_temp[i] <- period_annual_survival[1]
  }

  for(i in 1:(n_year_precollar + 1)) {
    period_effect_survival_temp[(yr_start[i] + nT_period_prestudy_ext):(yr_end[i] + nT_period_prestudy_ext)] <- period_annual_survival[i]
  }

  ### for the year of the study when switching from
  ### estimating period effects from collar data
  ### versus estimating from aah data

  period_effect_survival_temp[
          (yr_start[n_year_precollar] + nT_period_prestudy_ext):nT_period_precollar_ext] <-
      period_annual_survival[n_year_precollar + 1]

  ############################################################
  ## incorporating period effects from collar data
  ############################################################

  period_effect_survival_temp[(nT_period_precollar_ext + 1):
                               nT_period_overall_ext] <-
                               period_effect_surv[1:nT_period_collar]

  #making the period effects sum to zero using centering
#   mu_period <- mean(period_effect_survival_temp[1:nT_period_overall_ext])
#   period_effect_survival[1:nT_period_overall_ext] <-
#       period_effect_survival_temp -
#       mu_period

  returnType(double(1))
  return(period_effect_survival_temp[1:nT_period_overall_ext])
})

Cset_period_effects_constant <- compileNimble(set_period_effects_constant)

assign("set_period_effects_constant", set_period_effects_constant, envir = .GlobalEnv)


#test <- set_period_effects_constant(
#         n_year_precollar = n_year_precollar,
#         n_year_precollar_ext = n_year_precollar_ext,
#         n_year_prestudy_ext = n_year_prestudy_ext,
#         nT_period_precollar_ext = nT_period_precollar_ext,
#         nT_period_precollar = nT_period_precollar,
#         nT_period_collar = nT_period_collar,
#         nT_period_overall = nT_period_overall,
#         nT_period_overall_ext = nT_period_overall_ext,
#         nT_period_prestudy_ext = nT_period_prestudy_ext,
#         yr_start = d_fit_season$yr_start,
#         yr_end = d_fit_season$yr_end,
#         period_effect_surv = period_effect_survival_test[(nT_period_precollar_ext + 1):nT_period_overall_ext],
#         period_annual_survival = period_annual_survival)




#######################################################################
###
### Function to set and center period effects across full study time
### where period effects for non-harvest period effects
### based on aah data were different but constant each year
###
#######################################################################

set_period_effects_ave <- nimble::nimbleFunction(
    run = function(
        ### argument type declarations
        n_year_precollar = double(0),
        nT_period_precollar = double(0),
        nT_period_precollar_ext = double(0),
        nT_period_collar = double(0),
        nT_period_overall = double(0),
        nT_period_overall_ext = double(0),
        nT_period_prestudy_ext = double(0),
        yr_start = double(1),
        yr_end = double(1),
        period_effect_surv = double(1),
        period_annual_survival = double(1),
        indx_mat_pe_surv = double(2),
        intvl_step_yr = double(0)
        ) {

    period_effect_survival_temp <- nimNumeric(nT_period_overall_ext)
    period_effect_survival <- nimNumeric(nT_period_overall_ext)
    period_effect_within_yr <- nimNumeric(intvl_step_yr)


    ### Calculating the period effects by week, 
    ### averaged across all years with collar data

    for(i in 1:34){
        period_effect_within_yr[i] <- mean(period_effect_surv[indx_mat_pe_surv[2:6,i]])
    }
    for(i in 35:intvl_step_yr){
        period_effect_within_yr[i] <- mean(period_effect_surv[indx_mat_pe_surv[1:6,i]])
    }

    for(i in 1:n_year_precollar) {
        period_effect_survival_temp[(yr_start[i]+nT_period_prestudy_ext):
                                (yr_end[i]+nT_period_prestudy_ext)] <- 
            period_annual_survival[i] + period_effect_within_yr
    }

    ### for the year of the study when switching from
    ### estimating period effects from collar data
    ### versus estimating from aah data

    period_effect_survival_temp[
            (yr_start[n_year_precollar] + nT_period_prestudy_ext):
            (nT_period_precollar + nT_period_prestudy_ext)] <-
                period_annual_survival[n_year_precollar + 1] +
                period_effect_within_yr[1:length((yr_start[n_year_precollar]:nT_period_precollar))]

    ############################################################
    ## incorporating period effects from collar data
    ############################################################
    period_effect_survival_temp[(nT_period_precollar_ext+1): nT_period_overall_ext] <- period_effect_surv[1:nT_period_collar]

    #making the period effects sum to zero using centering
    ### only during the study though.... ignoring centering outside the study?
    # mu_period <- mean(period_effect_survival_temp[(1 + nT_period_prestudy_ext):
    #                                                (nT_period_overall + nT_period_prestudy_ext)])
    
    # period_effect_survival[1:nT_period_prestudy_ext] <-
    #     period_effect_survival_temp[1:nT_period_prestudy_ext] 

    # period_effect_survival[(1+nT_period_prestudy_ext):(nT_period_overall+nT_period_prestudy_ext)] <-
    #     period_effect_survival_temp - mu_period

    #### centering across all period effects, including years prior to study
    # mu_period <- mean(period_effect_survival_temp[1:nT_period_overall_ext])
    
    # period_effect_survival[1:nT_period_overall_ext] <- 
    #             period_effect_survival_temp[1:nT_period_overall_ext] - 
    #             mu_period

    returnType(double(1))
    return(period_effect_survival_temp[1:(nT_period_overall_ext)])
})

Cset_period_effects_ave <- compileNimble(set_period_effects_ave)

assign("set_period_effects_ave", set_period_effects_ave, envir = .GlobalEnv)


# set_period_effects_ave(
#         n_year_precollar = n_year_precollar,
#         nT_period_precollar = nT_period_precollar_ext,
#         nT_period_collar = nT_period_collar,
#         nT_period_overall = nT_period_overall,
#         nT_period_overall_ext = nT_period_overall_ext,
#         nT_period_prestudy_ext = nT_period_prestudy_ext,
#         yr_start = yr_start[1:n_year],
#         yr_end = yr_end[1:n_year],
#         period_effect_surv = period_effect_survival_test[nT_period_precollar_ext:nT_period_overall_ext],
#         period_annual_survival = period_annual_survival[1:(n_year_precollar + 1)],
#         indx_mat_pe_surv = indx_mat_pe_surv[1:6,1:intvl_step_yr_weekly],
#         intvl_step_yr = intvl_step_yr_weekly
#   )

# length(nT_period_prestudy_ext:nT_period_overall_ext)
# nT_period_precollar_ext:nT_period_overall_ext