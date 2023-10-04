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
        nT_period = double(0),
        nT_age_short = double(0),
        nT_age_surv_aah = double(0),
        beta0 = double(0),
        age_effect = double(1),
        period_effect = double(1),
        yr_start_age = double(1),
        yr_start_pop = double(1),
        n_year = double(0),
        n_age = double(0)
        ) {

    ###################################################################
    ###
    ### General Survival Surfaces for Susceptible/Infected Individuals
    ###
    ###################################################################

    mu_old_age_effect <- mean(age_effect[(nT_age_short + 1):nT_age])

    ############################################
    # Calculate hazards
    ############################################

    UCH <- nimMatrix(NA, nrow = nT_age_surv_aah, ncol = nT_period)
    s_aah <- nimMatrix(NA, nrow = n_age, ncol = n_year)

    ### Females
    for(j in 1:nT_period) {
        for(i in 1:nT_age_short) {
            UCH[i, j] <- exp(beta0 +
                             age_effect[i] +
                             period_effect[j])
        }
        for(i in (nT_age_short + 1):(nT_age_surv_aah)) {
            UCH[i, j] <- exp(beta0 +
                                mu_old_age_effect +
                                period_effect[j])
        }
    }
    ############################################
    # calculate survival from cummulative haz
    ############################################

    for (t in 1:n_year) {
        for (a in 1:n_age) {
            s_aah[a, t] <- exp(-sum(diag(UCH[
                               yr_start_age[a]:(yr_start_age[a] + 51),
                               yr_start_pop[t]:(yr_start_pop[t] + 51)])))
        }
    }

  returnType(double(2))
  return(s_aah[1:n_age, 1:n_year])
})

Ccalc_surv_aah <- compileNimble(calc_surv_aah)


# starttime <- Sys.time()
# sn_sus <- calc_surv_aah(
#     nT_age = nT_age,
#     nT_period = nT_period,
#     nT_age_short = nT_age_short,
#     nT_age_surv_aah = nT_age_surv_aah,
#     beta0 = beta0_survival_sus_true,
#     age_effect = age_effect_true,
#     period_effect = period_effect_true,
#     yr_start_age = yr_start_age,
#     yr_start_pop = yr_start_pop,
#     n_year = n_year,
#     n_age = n_age)
# (endtime1 <- Sys.time() - starttime)
# sn_sus

#######################################################################
###
### Function to calculate probability of infection
### based on FOI age and period effects
### Weekly Version
###
#######################################################################

calc_infect_prob <- nimbleFunction(
  run = function(age_lookup = double(1),
                 n_age = double(0),
                 yr_start = double(1),
                 age_foi = double(1),
                 nT_period = double(0),
                 n_year = double(0),
                 nT_age_surv_aah = double(0)) {

    gam <- nimMatrix(value = 0, nrow = nT_age_surv_aah, ncol = nT_period)
    p_inf <- nimMatrix(value = 0, nrow = n_age, ncol = n_year)

    for (t in 1:nT_period) {
        for (i in 1:nT_age_surv_aah) {
            gam[i, t] <- exp(age_foi[age_lookup[i]])

        }
    }
    # infection probability all ages all years 
    for (t in 1:n_year) {
        for (a in 1:n_age) {
            p_inf[a, t] <- 1 - exp(-sum(diag(gam[yr_start[a]:(yr_start[a] + 51),
                                        yr_start[t]:(yr_start[t] + 51)])))
        }
    }
    returnType(double(2))
    return(p_inf[1:n_age, 1:n_year])
  })

Ccalc_infect_prob <- compileNimble(calc_infect_prob)

assign("calc_infect_prob", calc_infect_prob, envir = .GlobalEnv)

# #testing state.transition function as R function
# starttime <- Sys.time()
# psi <- Ccalc_infect_prob(age_lookup = age_lookup,
#                         n_age = n_age,
#                         yr_start = yr_start_age,
#                         age_foi = age_foi,
#                         nT_period = nT_period,
#                         n_year = n_year,
#                         nT_age_surv_aah = nT_age_surv_aah
#                         )

# (endtime5 <- Sys.time() - starttime)
# psi


##################################################################
###
### Calculating true values of derived parameters
###
##################################################################


sn_sus_true <- calc_surv_aah(
    nT_age = nT_age,
    nT_period = nT_period,
    nT_age_short = nT_age_short,
    nT_age_surv_aah = nT_age_surv_aah,
    beta0 = beta0_survival_sus_true,
    age_effect = age_effect_true,
    period_effect = period_effect_true,
    yr_start_age = yr_start_age,
    yr_start_pop = yr_start_pop,
    n_year = n_year,
    n_age = n_age)


sn_inf_true <- calc_surv_aah(
    nT_age = nT_age,
    nT_period = nT_period,
    nT_age_short = nT_age_short,
    nT_age_surv_aah = nT_age_surv_aah,
    beta0 = beta0_survival_inf_true,
    age_effect = age_effect_true,
    period_effect = period_effect_true,
    yr_start_age = yr_start_age,
    yr_start_pop = yr_start_pop,
    n_year = n_year,
    n_age = n_age)

psi_true <- Ccalc_infect_prob(age_lookup = age_lookup,
                        n_age = n_age,
                        yr_start = yr_start_age,
                        age_foi = age_foi,
                        nT_period = nT_period,
                        n_year = n_year,
                        nT_age_surv_aah = nT_age_surv_aah
                        )
psi_true <- psi_true[,1]
psi_true <- psi_true[-(n_age-1)]
