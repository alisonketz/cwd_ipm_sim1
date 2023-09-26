#######################################################################
###
###   Likelihoods for each kind of data listed below
###
#######################################################################

# d_fit_hunt_neg
# d_fit_hunt_pos
# d_fit_sus_cens_posttest
# d_fit_sus_cens_postno
# d_fit_sus_mort_posttest
# d_fit_sus_mort_postno
# d_fit_icap_cens
# d_fit_icap_mort
# d_fit_rec_neg_mort
# d_fit_rec_neg_cens
# d_fit_rec_pos_mort
# d_fit_rec_pos_cens
# d_fit_idead
# d_fit_endlive
# d_aah

#######################################################################
###
###   User defined distribution for likelihood for
###   infected harvest deer, revised for multiple deer
###
###   d_fit_hunt_pos
###   Overleaf Equation 3
###
#######################################################################

# dInfHarvest <- nimble::nimbleFunction( # nolint
#     run = function( # nolint
#                    ### argument type declarations
#                    x = integer(0),
#                    n_cases = double(1),
#                    n_samples = integer(0), # number of samples in dataset
#                    a = double(1), # age (weeks) at harvest
#                    sex = double(1),
#                    age2date = double(1),
#                    beta_male = double(0),
#                    beta0_sus = double(0),
#                    beta0_inf = double(0),
#                    age_effect_surv = double(1),
#                    period_effect_surv = double(1),
#                    f_age_foi = double(1),
#                    m_age_foi = double(1),
#                    age_lookup_f = double(1),
#                    age_lookup_m = double(1),
#                    period_lookup_foi = double(1),
#                    f_period_foi = double(1),
#                    m_period_foi = double(1),
#                    space = double(1),
#                    sect = double(1),
#                    log = double(0)) {

#         # start the loop through individuals
#         sumllik <- 0
#         for (i in 1:n_samples) {

#             # intitialize scalars
#             lam_foi <- 0
#             lam_sus <- 0
#             lam_inf <- 0
#             lik_temp <- 0

#             #################################################
#             ### loop over ages and accumulate sums over 1:a-1
#             ### have to loop separately for lam_inf
#             #################################################

#             if (sex[i] == 0) { # age loops for females
#                 for (j in 1:(a[i] - 1)) {
#                     # sum up infected hazard from 1 to a-1
#                     lam_inf <- lam_inf +
#                         exp(beta0_inf +
#                             age_effect_surv[j] +
#                             period_effect_surv[age2date[i] + j])
#                 }
#                 for (j in 1:a[i]) {
#                     # sum up foi hazard from 1  to j
#                     lam_foij <- exp(space[sect[i]] +
#                         f_age_foi[age_lookup_f[j]] +
#                         f_period_foi[period_lookup_foi[age2date[i] + j]])
#                     # sum up like_temp (no sus hazard when j=1)
#                     lik_temp <- lik_temp +
#                         lam_foij *
#                             exp(-(lam_sus + lam_foi)) *
#                             exp(-lam_inf)
#                     # now increment the hazards (subtract for lam_inf)
#                     lam_foi <- lam_foi + lam_foij
#                     if (j < a[i]) {
#                         lam_sus <- lam_sus +
#                             exp(beta0_sus +
#                                 age_effect_surv[j] +
#                                 period_effect_surv[age2date[i] + j])
#                         lam_inf <- lam_inf -
#                             exp(beta0_inf +
#                                 age_effect_surv[j] +
#                                 period_effect_surv[age2date[i] + j])
#                     } else { # now calculate infected hazard for final age
#                         lam_inf <- exp(beta0_inf +
#                             age_effect_surv[j] +
#                             period_effect_surv[age2date[i] + j])
#                     }
#                 }
#             } else { # age loops for males
#                 for (j in 1:(a[i] - 1)) {
#                     # sum up infected hazard from 1 to a-1
#                     lam_inf <- lam_inf +
#                         exp(beta0_inf +
#                             age_effect_surv[j] +
#                             period_effect_surv[age2date[i] + j] +
#                             beta_male)
#                 }
#                 for (j in 1:(a[i])) {
#                     # sum up foi hazard from 1  to j
#                     lam_foij <- exp(space[sect[i]] +
#                         m_age_foi[age_lookup_f[j]] +
#                         m_period_foi[period_lookup_foi[age2date[i] + j]])
#                     # sum up like_temp (no sus hazard when j=1)
#                     lik_temp <- lik_temp +
#                         lam_foij * exp(-(lam_sus + lam_foi)) *
#                             exp(-lam_inf)
#                     # now increment the hazards (subtract for lam_inf)
#                     lam_foi <- lam_foi + lam_foij
#                     if (j < a[i]) {
#                         lam_sus <- lam_sus +
#                             exp(beta0_sus +
#                                 age_effect_surv[j] +
#                                 period_effect_surv[age2date[i] + j] +
#                                 beta_male)
#                         lam_inf <- lam_inf -
#                             exp(beta0_inf +
#                                 age_effect_surv[j] +
#                                 period_effect_surv[age2date[i] + j] +
#                                 beta_male)
#                     } else {
#                         # now calculate infected hazard for final age
#                         lam_inf <- exp(beta0_inf +
#                             age_effect_surv[j] +
#                             period_effect_surv[age2date[i] + j] +
#                             beta_male)
#                     }
#                 }
#             } # end if(sex)
#             #######################################
#             ### accumulate the joint likelihood
#             #######################################
#             sumllik <- sumllik +
#                 (log(lam_inf) + log(lik_temp)) * n_cases[i]
#         } # end the loop for individual i
#         returnType(double(0))
#         if (log) {
#             return(sumllik)
#         } else {
#             return(exp(sumllik))
#         }
#     }
# )


dInfHarvest <- nimble::nimbleFunction( # nolint
    run = function( # nolint
                   ### argument type declarations
                   x = integer(0),
                   n_cases = double(1),
                   n_samples = integer(0), # number of samples in dataset
                   a = double(1), # age (weeks) at harvest
                   sex = double(1),
                   age2date = double(1),
                   beta_male = double(0),
                   beta0_sus = double(0),
                   beta0_inf = double(0),
                   age_effect_surv = double(1),
                   period_effect_surv = double(1),
                   f_age_foi = double(1),
                   m_age_foi = double(1),
                   age_lookup_f = double(1),
                   age_lookup_m = double(1),
                   period_lookup_foi = double(1),
                   f_period_foi = double(1),
                   m_period_foi = double(1),
                   space = double(1),
                   sect = double(1),
                   log = double(0)) {

        # start the loop through individuals
        sumllik <- 0
        for (i in 1:n_samples) {

            # intitialize scalars
            lam_foi <- 0
            lam_sus <- 0
            lam_inf <- 0
            lik_temp <- 0

            #################################################
            ### loop over ages and accumulate sums over 1:a-1
            ### have to loop separately for lam_inf
            #################################################

            if (sex[i] == 0) { # age loops for females
                for (j in 1:(a[i] - 2)) {
                    # sum up infected hazard from 1 to a-1
                    lam_inf <- lam_inf +
                        exp(beta0_inf +
                            age_effect_surv[j] +
                            period_effect_surv[age2date[i] + j])
                }
                for (j in 1:(a[i] - 1)) {
                    # sum up foi hazard from 1  to j
                    lam_foij <- exp(space[sect[i]] +
                        f_age_foi[age_lookup_f[j]] +
                        f_period_foi[period_lookup_foi[age2date[i] + j]])
                    # sum up like_temp (no sus hazard when j=1)
                    lik_temp <- lik_temp +
                        lam_foij *
                            exp(-(lam_sus + lam_foi)) *
                            exp(-lam_inf)
                    # now increment the hazards (subtract for lam_inf)
                    lam_foi <- lam_foi + lam_foij
                    if (j < (a[i] - 1)) {
                        lam_sus <- lam_sus +
                            exp(beta0_sus +
                                age_effect_surv[j] +
                                period_effect_surv[age2date[i] + j])
                        lam_inf <- lam_inf -
                            exp(beta0_inf +
                                age_effect_surv[j] +
                                period_effect_surv[age2date[i] + j])
                    } else { # now calculate infected hazard for final age
                        lam_inf <- exp(beta0_inf +
                            age_effect_surv[j] +
                            period_effect_surv[age2date[i] + j])
                    }
                }
            } else { # age loops for males
                for (j in 1:(a[i] - 2)) {
                    # sum up infected hazard from 1 to a-1
                    lam_inf <- lam_inf +
                        exp(beta0_inf +
                            age_effect_surv[j] +
                            period_effect_surv[age2date[i] + j] +
                            beta_male)
                }
                for (j in 1:(a[i] - 1)) {
                    # sum up foi hazard from 1  to j
                    lam_foij <- exp(space[sect[i]] +
                        m_age_foi[age_lookup_f[j]] +
                        m_period_foi[period_lookup_foi[age2date[i] + j]])
                    # sum up like_temp (no sus hazard when j=1)
                    lik_temp <- lik_temp +
                        lam_foij * exp(-(lam_sus + lam_foi)) *
                            exp(-lam_inf)
                    # now increment the hazards (subtract for lam_inf)
                    lam_foi <- lam_foi + lam_foij
                    if (j < (a[i] - 1)) {
                        lam_sus <- lam_sus +
                            exp(beta0_sus +
                                age_effect_surv[j] +
                                period_effect_surv[age2date[i] + j] +
                                beta_male)
                        lam_inf <- lam_inf -
                            exp(beta0_inf +
                                age_effect_surv[j] +
                                period_effect_surv[age2date[i] + j] +
                                beta_male)
                    } else {
                        # now calculate infected hazard for final age
                        lam_inf <- exp(beta0_inf +
                            age_effect_surv[j] +
                            period_effect_surv[age2date[i] + j] +
                            beta_male)
                    }
                }
            } # end if(sex)
            #######################################
            ### accumulate the joint likelihood
            #######################################
            sumllik <- sumllik +
                (log(lam_inf) + log(lik_temp)) * n_cases[i]
        } # end the loop for individual i
        returnType(double(0))
        if (log) {
            return(sumllik)
        } else {
            return(exp(sumllik))
        }
    }
)



nimble::registerDistributions(list(
    dInfHarvest = list(
        BUGSdist = "dInfHarvest(n_cases,n_samples,a,sex,age2date,beta_male,beta0_sus,beta0_inf,age_effect_surv,period_effect_surv,f_age_foi,m_age_foi,age_lookup_f,age_lookup_m,f_period_foi,m_period_foi,period_lookup_foi,space,sect)",
        types = c(
            "value = integer(0)",
            "n_cases = double(1)",
            "n_samples = integer(0)",
            "a = double(1)",
            "sex = double(1)",
            "age2date = double(1)",
            "beta_male = double(0)",
            "beta0_sus = double(0)",
            "beta0_inf = double(0)",
            "age_effect_surv = double(1)",
            "period_effect_surv = double(1)",
            "f_age_foi = double(1)",
            "m_age_foi = double(1)",
            "age_lookup_f = double(1)",
            "age_lookup_m = double(1)",
            "f_period_foi = double(1)",
            "m_period_foi = double(1)",
            "period_lookup_foi = double(1)",
            "space = double(1)",
            "sect = double(1)",
            "returnType = double(0)",
            "log = double(0)"
        ),
        discrete = TRUE
    )
))

# for a user-defined distribution
assign("dInfHarvest",
    dInfHarvest,
    envir = .GlobalEnv
)

# starttime <- Sys.time()
# test <- dInfHarvest(
#     x = 1,
#     n_cases = d_fit_hunt_pos$n_cases,
#     n_samples = nrow(d_fit_hunt_pos),
#     a = d_fit_hunt_pos$ageweeks,
#     sex = d_fit_hunt_pos$sex,
#     age2date = d_fit_hunt_pos$birthweek - 1,
#     beta_male = beta_male,
#     beta0_sus = beta0_survival_sus,
#     beta0_inf = beta0_survival_inf,
#     age_effect_surv = age_effect_survival_test,
#     period_effect_surv = period_effect_survival_test,
#     f_age_foi = f_age_foi,
#     m_age_foi = m_age_foi,
#     age_lookup_f = age_lookup_f,
#     age_lookup_m = age_lookup_m,
#     period_lookup_foi = period_lookup_foi,
#     f_period_foi = f_period_foi,
#     m_period_foi = m_period_foi,
#     space = c(0, -.55),
#     sect = d_fit_hunt_pos$ew,
#     log = TRUE
# )
# (end <- Sys.time() - starttime)
# test

#######################################################################
###
###   User defined distribution for likelihood for
###   uninfected harvest deer
###
###   d_fit_hunt_neg
###
###   Overleaf Equation (5)
###
#######################################################################

# dSusHarvest <- nimble::nimbleFunction(
#     run = function( ### argument type declarations
#                    x = double(0),
#                    n_cases = double(1),
#                    n_samples = integer(0), # number of samples in dataset
#                    a = double(1), # age (weeks) at harvest
#                    sex = double(1),
#                    age2date = double(1),
#                    beta_male = double(0),
#                    beta0_sus = double(0),
#                    age_effect_surv = double(1),
#                    period_effect_surv = double(1),
#                    f_age_foi = double(1),
#                    m_age_foi = double(1),
#                    age_lookup_f = double(1),
#                    age_lookup_m = double(1),
#                    period_lookup_foi = double(1),
#                    f_period_foi = double(1),
#                    m_period_foi = double(1),
#                    space = double(1),
#                    sect = double(1),
#                    log = double(0)) {

#         # start the loop through individuals
#         sumllik <- 0
#         for (i in 1:n_samples) {
#             # intitialize vectors
#             lam_foi <- 0
#             lam_sus <- 0
#             sum_haz <- 0

#             #############################################
#             # preliminary hazards for the likelihood
#             #############################################
#             if (sex[i] == 0) { # age loops for females
#                 for (j in 1:(a[i] - 1)) {
#                     # sum up foi and sus hazards from 1  to a-1
#                     sum_haz <- sum_haz +
#                         exp(space[sect[i]] +
#                             f_age_foi[age_lookup_f[j]] +
#                             f_period_foi[period_lookup_foi[age2date[i] + j]]) +
#                         exp(beta0_sus +
#                             age_effect_surv[j] +
#                             period_effect_surv[age2date[i] + j])
#                 }
#                 # now hazards for a
#                 lam_sus <- exp(beta0_sus +
#                     age_effect_surv[a[i]] +
#                     period_effect_surv[age2date[i] + a[i]])
#                 lam_foi <- exp(space[sect[i]] +
#                     f_age_foi[age_lookup_f[a[i]]] +
#                     f_period_foi[period_lookup_foi[age2date[i] + a[i]]])
#             } else { # age loops for males
#                 for (j in 1:(a[i] - 1)) {
#                     # sum up foi and sus hazards from 1  to a-1
#                     sum_haz <- sum_haz +
#                         exp(space[sect[i]] +
#                             m_age_foi[age_lookup_f[j]] +
#                             m_period_foi[period_lookup_foi[age2date[i] + j]]) +
#                         exp(beta0_sus +
#                             age_effect_surv[j] +
#                             period_effect_surv[age2date[i] + j]+
#                             beta_male)
#                 }
#                 # now hazards for a
#                 lam_sus <- exp(beta0_sus +
#                     age_effect_surv[a[i]] +
#                     period_effect_surv[age2date[i] + a[i]] +
#                     beta_male)
#                 lam_foi <- exp(space[sect[i]] +
#                     m_age_foi[age_lookup_f[a[i]]] +
#                     m_period_foi[period_lookup_foi[age2date[i] + a[i]]])
#             }

#             #######################################
#             ###
#             ### calculating the joint likelihood
#             ###
#             #######################################

#             sumllik <- sumllik +
#                 (-sum_haz - lam_foi + log(lam_sus)) * n_cases[i]
#         }
#         returnType(double(0))
#         if (log) {
#             return(sumllik)
#         } else {
#             return(exp(sumllik))
#         } ## return log-likelihood
#     }
# )

dSusHarvest <- nimble::nimbleFunction(
    run = function( ### argument type declarations
                   x = double(0),
                   n_cases = double(1),
                   n_samples = integer(0), # number of samples in dataset
                   a = double(1), # age (weeks) at harvest
                   sex = double(1),
                   age2date = double(1),
                   beta_male = double(0),
                   beta0_sus = double(0),
                   age_effect_surv = double(1),
                   period_effect_surv = double(1),
                   f_age_foi = double(1),
                   m_age_foi = double(1),
                   age_lookup_f = double(1),
                   age_lookup_m = double(1),
                   period_lookup_foi = double(1),
                   f_period_foi = double(1),
                   m_period_foi = double(1),
                   space = double(1),
                   sect = double(1),
                   log = double(0)) {

        # start the loop through individuals
        sumllik <- 0
        for (i in 1:n_samples) {
            # intitialize vectors
            lam_foi <- 0
            lam_sus <- 0
            sum_haz <- 0

            #############################################
            # preliminary hazards for the likelihood
            #############################################
            if (sex[i] == 0) { # age loops for females
                for (j in 1:(a[i] - 2)) {
                    # sum up foi and sus hazards from 1  to a-1
                    sum_haz <- sum_haz +
                        exp(space[sect[i]] +
                            f_age_foi[age_lookup_f[j]] +
                            f_period_foi[period_lookup_foi[age2date[i] + j]]) +
                        exp(beta0_sus +
                            age_effect_surv[j] +
                            period_effect_surv[age2date[i] + j])
                }
                # now hazards for a
                lam_sus <- exp(beta0_sus +
                    age_effect_surv[a[i] - 1] +
                    period_effect_surv[age2date[i] + a[i] - 1])
            } else { # age loops for males
                for (j in 1:(a[i] - 2)) {
                    # sum up foi and sus hazards from 1  to a-1
                    sum_haz <- sum_haz +
                        exp(space[sect[i]] +
                            m_age_foi[age_lookup_f[j]] +
                            m_period_foi[period_lookup_foi[age2date[i] + j]]) +
                        exp(beta0_sus +
                            age_effect_surv[j] +
                            period_effect_surv[age2date[i] + j]+
                            beta_male)
                }
                # now hazards for a
                lam_sus <- exp(beta0_sus +
                    age_effect_surv[a[i] - 1] +
                    period_effect_surv[age2date[i] + a[i] - 1] +
                    beta_male)
            }

            #######################################
            ###
            ### calculating the joint likelihood
            ###
            #######################################

            sumllik <- sumllik +
                (-sum_haz + log(lam_sus)) * n_cases[i]
        }
        returnType(double(0))
        if (log) {
            return(sumllik)
        } else {
            return(exp(sumllik))
        } ## return log-likelihood
    }
)

nimble::registerDistributions(list(
    dSusHarvest = list(
        BUGSdist = "dSusHarvest(n_cases,n_samples,a,sex,age2date,beta_male,beta0_sus,age_effect_surv,period_effect_surv,f_age_foi,m_age_foi,age_lookup_f,age_lookup_m,f_period_foi,m_period_foi,period_lookup_foi,space,sect)",
        types = c(
            "value = double(0)",
            "n_cases = double(1)",
            "a = double(1)",
            "n_samples = integer(0)",
            "sex = double(1)",
            "age2date = double(1)",
            "beta_male = double(0)",
            "beta0_sus = double(0)",
            "age_effect_surv = double(1)",
            "period_effect_surv = double(1)",
            "f_age_foi = double(1)",
            "m_age_foi = double(1)",
            "age_lookup_f = double(1)",
            "age_lookup_m = double(1)",
            "period_lookup_foi=double(1)",
            "f_period_foi=double(1)",
            "m_period_foi=double(1)",
            "space = double(1)",
            "sect = double(1)",
            "log = double(0)"
        ),
        discrete = TRUE
    )
))

# for a user-defined distribution
assign("dSusHarvest", dSusHarvest, envir = .GlobalEnv)

# starttime <- Sys.time()
#  test  <- dSusHarvest(
#     x = 1,
#     n_cases = d_fit_hunt_neg$n_cases,
#     n_samples = nrow(d_fit_hunt_neg),
#     a = d_fit_hunt_neg$ageweeks,
#     sex = d_fit_hunt_neg$sex,
#     age2date = d_fit_hunt_neg$birthweek - 1,
#     beta_male = beta_male,
#     beta0_sus = beta0_survival_sus,
#     age_effect_surv = age_effect_survival_test,
#     period_effect_surv = period_effect_survival_test,
#     f_age_foi = f_age_foi,
#     m_age_foi = m_age_foi,
#     age_lookup_f = age_lookup_f,
#     age_lookup_m = age_lookup_m,
#     period_lookup_foi = period_lookup_foi,
#     f_period_foi = f_period_foi,
#     m_period_foi = m_period_foi,
#     space = c(0, -.55),
#     sect = d_fit_hunt_neg$ew,
#     log = TRUE
# )
# (end <- Sys.time() - starttime)
# test


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

dSusCensTest <- nimble::nimbleFunction(
    run = function( ### argument type declarations
                   x = integer(0),
                   n_samples = integer(0),
                   e = double(1), # e, age of entry
                   r = double(1), # r, age of last known alive
                   sex = double(1),
                   age2date = double(1),
                   beta_male = double(0),
                   beta0_sus = double(0),
                   age_effect_surv = double(1),
                   period_effect_surv = double(1),
                   f_age_foi = double(1),
                   m_age_foi = double(1),
                   age_lookup_f = double(1),
                   age_lookup_m = double(1),
                   period_lookup_foi = double(1),
                   f_period_foi = double(1),
                   m_period_foi = double(1),
                   sect = double(1),
                   space = double(1),
                   log = double()) {

        # starttime the loop through individuals
        sumllik <- 0
        for (i in 1:n_samples) {
            # intitialize vectors
            lam_foi <- 0
            lam_sus <- 0

            #############################################
            # preliminary hazards for the likelihood
            #############################################
            if (sex[i] == 0) { # age loops for females
                for (j in 1:(r[i] - 1)) {
                    # sum up foi
                    lam_foi <- lam_foi +
                        exp(space[sect[i]] +
                            f_age_foi[age_lookup_f[j]] +
                            f_period_foi[period_lookup_foi[age2date[i] + j]])
                    if (j > (e[i] - 1)) {
                        lam_sus <- lam_sus +
                            exp(beta0_sus +
                                age_effect_surv[j] +
                                period_effect_surv[age2date[i] + j])
                    }
                }
            } else { # age loops for males
                for (j in 1:(r[i] - 1)) {
                    # sum up foi
                    lam_foi <- lam_foi +
                        exp(space[sect[i]] +
                            m_age_foi[age_lookup_m[j]] +
                            m_period_foi[period_lookup_foi[age2date[i] + j]])
                    if (j > (e[i] - 1)) {
                        lam_sus <- lam_sus +
                            exp(beta0_sus +
                                age_effect_surv[j] +
                                period_effect_surv[age2date[i] + j] +
                                beta_male)
                    }
                }
            }
            #######################################
            ### calculating the joint likelihood
            #######################################

            sumllik <- sumllik - lam_sus - lam_foi
        } # end loop over individuals

        returnType(double(0))
        if (log) {
            return(sumllik)
        } else {
            return(exp(sumllik))
        } ## return log-likelihood
    }
)


nimble::registerDistributions(list(
    dSusCensTest = list(
        BUGSdist = "dSusCensTest(n_samples,e,r,sex,age2date,beta_male,beta0_sus,age_effect_surv,period_effect_surv,f_age_foi,m_age_foi,age_lookup_f,age_lookup_m,f_period_foi,m_period_foi,period_lookup_foi,sect,space)",
        types = c(
            "value = integer(0)",
            "n_samples = integer(0)",
            "e = double(1)",
            "r = double(1)",
            "sex = double(1)",
            "age2date = double(1)",
            "beta_male = double(0)",
            "beta0_sus = double(0)",
            "age_effect_surv = double(1)",
            "period_effect_surv = double(1)",
            "f_age_foi = double(1)",
            "m_age_foi = double(1)",
            "age_lookup_f = double(1)",
            "age_lookup_m = double(1)",
            "period_lookup_foi = double(1)",
            "f_period_foi = double(1)",
            "m_period_foi = double(1)",
            "sect = double(1)",
            "space = double(1)",
            "log = double(0)"
        ),
        discrete = TRUE
    )
))


### for a user-defined distribution
assign("dSusCensTest", dSusCensTest, envir = .GlobalEnv)

# 		n_samples = nrow(d_fit_sus_cens_posttest)
#         e = d_fit_sus_cens_posttest$left_age_e
#         r = d_fit_sus_cens_posttest$right_age_r
#         sex = d_fit_sus_cens_posttest$sex
#         age2date = sus_cens_posttest_age2date
#         beta_male = beta_male
#         beta0_sus = beta0_survival_sus
#         age_effect_surv = age_effect_survival_test
#         period_effect_surv = period_effect_survival_test
#         f_age_foi = f_age_foi
#         m_age_foi = m_age_foi
#         age_lookup_f = age_lookup_f
#         age_lookup_m = age_lookup_m
#         period_lookup_foi = period_lookup_foi
#         f_period_foi = f_period_foi
#         m_period_foi = m_period_foi
#         sect = d_fit_sus_cens_posttest$study_area
#         space = c(0,-.55)

# starttime <- Sys.time()
# test <- dSusCensTest(
#         x = 1,
# 		n_samples = nrow(d_fit_sus_cens_posttest),
#         e = d_fit_sus_cens_posttest$left_age_e,
#         r = d_fit_sus_cens_posttest$right_age_r,
#         sex = d_fit_sus_cens_posttest$sex,
#         age2date = sus_cens_posttest_age2date,
#         beta_male = beta_male,
#         beta0_sus = beta0_survival_sus,
#         age_effect_surv = age_effect_survival_test,
#         period_effect_surv = period_effect_survival_test,
#         f_age_foi = f_age_foi,
#         m_age_foi = m_age_foi,
#         age_lookup_f = age_lookup_f,
#         age_lookup_m = age_lookup_m,
#         period_lookup_foi = period_lookup_foi,
#         f_period_foi = f_period_foi,
#         m_period_foi = m_period_foi,
#         sect = d_fit_sus_cens_posttest$study_area,
#         space = c(0,-.55),
#         log = TRUE
#         )
# (end<- Sys.time()-starttime)
# test


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

dSusCensNo <- nimble::nimbleFunction(
    run = function( ### argument type declarations
                   x = integer(0),
                   n_samples = integer(0), # number of samples in dataset
                   e = double(1), # e, age of entry
                   r = double(1), # r, age of last known alive
                   sex = double(1),
                   age2date = double(1),
                   beta_male = double(0),
                   beta0_sus = double(0),
                   beta0_inf = double(0),
                   age_effect_surv = double(1),
                   period_effect_surv = double(1),
                   f_age_foi = double(1),
                   m_age_foi = double(1),
                   age_lookup_f = double(1),
                   age_lookup_m = double(1),
                   period_lookup_foi = double(1),
                   f_period_foi = double(1),
                   m_period_foi = double(1),
                   space = double(1),
                   sect = double(1),
                   log = double()) {
        sumllik <- 0 # intialize log-likelihood
        for (i in 1:n_samples) {
            lam_foi <- 0
            lam_sus <- 0
            lam_inf <- 0
            lik_temp <- 0

            #############################################
            # preliminary hazards for the likelihood
            #############################################
            if (sex[i] == 0) { # age loops for females
                for (j in (e[i] + 1):(r[i] - 1)) {
                    # sum up infected hazard from e to r
                    lam_inf <- lam_inf +
                        exp(beta0_inf +
                            age_effect_surv[j] +
                            period_effect_surv[age2date[i] + j])
                }
                for (j in 1:e[i]) {
                    # sum up infected hazard from 1 to e
                    lam_foi <- lam_foi +
                        exp(space[sect[i]] +
                            f_age_foi[age_lookup_f[j]] +
                            f_period_foi[period_lookup_foi[age2date[i] + j]])
                    if (j == e[i]) {
                        lam_sus_e <- exp(beta0_sus +
                            age_effect_surv[j] +
                            period_effect_surv[age2date[i] + j])
                        lik_temp_e <- exp(-lam_sus_e) * exp(-lam_foi)
                    }
                }
                # for(j = (e[i]+1)) {
                # sum up foi hazard from e to r
                lam_foij <- exp(space[sect[i]] +
                    f_age_foi[age_lookup_f[e[i] + 1]] +
                    f_period_foi[period_lookup_foi[age2date[i] + e[i] + 1]])
                # sum up like_temp
                lik_temp <- lik_temp + lam_foij *
                    exp(-lam_inf)
                # now redefine hazards (subtract for lam_inf)
                lam_sus <- exp(beta0_sus +
                    age_effect_surv[e[i] + 1] +
                    period_effect_surv[age2date[i] + e[i] + 1])
                lam_foie <- lam_foij
                lam_foi <- lam_foi + lam_foij
                lam_inf <- lam_inf -
                    exp(beta0_inf +
                        age_effect_surv[e[i] + 1] +
                        period_effect_surv[age2date[i] + e[i] + 1])

                if ((r[i] - e[i]) > 2) {
                    for (j in (e[i] + 2):(r[i] - 1)) {
                        # sum up foi hazard from e to r
                        lam_foij <- exp(space[sect[i]] +
                            f_age_foi[age_lookup_f[j]] +
                            f_period_foi[period_lookup_foi[age2date[i] + j]])
                        # sum up like_temp
                        lik_temp <- lik_temp + lam_foij *
                            exp(-(lam_sus + lam_foie)) *
                            exp(-lam_inf)
                        # now increment the hazards (subtract for lam_inf)
                        lam_foie <- lam_foie + lam_foij
                        lam_foi <- lam_foi + lam_foij
                        lam_sus <- lam_sus +
                            exp(beta0_sus +
                                age_effect_surv[j] +
                                period_effect_surv[age2date[i] + j])
                        lam_inf <- lam_inf -
                            exp(beta0_inf +
                                age_effect_surv[j] +
                                period_effect_surv[age2date[i] + j])
                    }
                }
            } else { # age loops for males
                for (j in (e[i] + 1):(r[i] - 1)) {
                    # sum up infected hazard from e to r
                    lam_inf <- lam_inf +
                               exp(beta0_inf +
                                   age_effect_surv[j] +
                                   period_effect_surv[age2date[i] + j] +
                                   beta_male)
                }
                for (j in 1:e[i]) {
                    # sum up infected hazard from 1 to e
                    lam_foi <- lam_foi +
                               exp(space[sect[i]] +
                               m_age_foi[age_lookup_m[j]] +
                               m_period_foi[period_lookup_foi[age2date[i] + j]])
                    if (j == e[i]) {
                        lam_sus_e <- exp(beta0_sus +
                                         age_effect_surv[j] +
                                         period_effect_surv[age2date[i] + j] +
                                         beta_male)
                        lik_temp_e <- exp(-lam_sus_e) * exp(-lam_foi)
                    }
                }
                # for(j = (e[i]+1)) {
                # sum up foi hazard from e to r
                lam_foij <- exp(space[sect[i]] +
                    m_age_foi[age_lookup_m[e[i] + 1]] +
                    m_period_foi[period_lookup_foi[age2date[i] + e[i] + 1]])
                # sum up like_temp
                lik_temp <- lik_temp +
                            lam_foij * exp(-lam_inf)
                # now redefine the hazards (subtract for lam_inf)
                lam_sus <- exp(beta0_sus +
                               age_effect_surv[e[i] + 1] +
                               period_effect_surv[age2date[i] + e[i] + 1] +
                               beta_male)
                lam_foie <- lam_foij
                lam_foi <- lam_foi + lam_foij
                lam_inf <- lam_inf -
                           exp(beta0_inf +
                               age_effect_surv[e[i] + 1] +
                               period_effect_surv[age2date[i] + e[i] + 1] +
                               beta_male)

                if ((r[i] - e[i]) > 2) {
                    for (j in (e[i] + 2):(r[i] - 1)) {
                        # sum up foi hazard from e to r
                        lam_foij <- exp(space[sect[i]] +
                                        m_age_foi[age_lookup_m[j]] +
                                        m_period_foi[period_lookup_foi[age2date[i] + j]])
                        # sum up like_temp
                        lik_temp <- lik_temp +
                        lam_foij * exp(-(lam_sus + lam_foie)) * exp(-lam_inf)
                        # now increment the hazards (subtract for lam_inf)
                        lam_foie <- lam_foie + lam_foij
                        lam_foi <- lam_foi + lam_foij
                        lam_sus <- lam_sus +
                                   exp(beta0_sus +
                                       age_effect_surv[j] +
                                       period_effect_surv[age2date[i] + j] +
                                       beta_male)
                        lam_inf <- lam_inf -
                                   exp(beta0_inf +
                                       age_effect_surv[j] +
                                       period_effect_surv[age2date[i] + j] +
                                       beta_male)
                    }
                }
            }
            #######################################
            ### calculating the joint likelihood
            #######################################
            sumllik <- sumllik +
                       log(lik_temp * lik_temp_e +
                       exp(-(lam_sus + lam_sus_e)) * exp(-lam_foi))
        }
        returnType(double(0))
        if (log) {
            return(sumllik)
        } else {
            return(exp(sumllik))
        } ## return log-likelihood
    }
)

nimble::registerDistributions(list(
    dSusCensNo = list(
        BUGSdist = "dSusCensNo(n_samples,e,r,sex,age2date,beta_male,beta0_sus,beta0_inf,age_effect_surv,period_effect_surv,f_age_foi,m_age_foi,age_lookup_f,age_lookup_m,f_period_foi,m_period_foi,period_lookup_foi,space,sect)",
        types = c(
            "value=integer(0)",
            "n_samples = integer(0)",
            "e = double(1)",
            "r = double(1)",
            "sex = double(1)",
            "age2date = double(1)",
            "beta_male = double(0)",
            "beta0_sus = double(0)",
            "beta0_inf = double(0)",
            "age_effect_surv = double(1)",
            "period_effect_surv = double(1)",
            "f_age_foi = double(1)",
            "m_age_foi = double(1)",
            "age_lookup_f = double(1)",
            "age_lookup_m = double(1)",
            "period_lookup_foi=double(1)",
            "f_period_foi=double(1)",
            "m_period_foi=double(1)",
            "space = double(1)",
            "sect = double(1)",
            "log = double()"
        ),
        discrete = TRUE
    )
))

### for a user-defined distribution
assign("dSusCensNo", dSusCensNo, envir = .GlobalEnv)

#         n_samples = nrow(d_fit_sus_cens_postno)
#         e = d_fit_sus_cens_postno$left_age_e
#         r = d_fit_sus_cens_postno$right_age_r
#         sex = d_fit_sus_cens_postno$sex
#         age2date = sus_cens_postno_age2date
#         beta_male = beta_male
#         beta0_sus = beta0_survival_sus
#         beta0_inf = beta0_survival_inf
#         age_effect_surv = age_effect_survival_test
#         period_effect_surv = period_effect_survival_test
#         f_age_foi = f_age_foi
#         m_age_foi = m_age_foi
#         age_lookup_f = age_lookup_f
#         age_lookup_m = age_lookup_m
#         period_lookup_foi = period_lookup_foi
#         f_period_foi = f_period_foi
#         m_period_foi = m_period_foi
#         space = c(0,-.55)
#         sect = d_fit_sus_cens_postno$study_area


# starttime <- Sys.time()
# test <- dSusCensNo(
#         x = 1,
#         n_samples = nrow(d_fit_sus_cens_postno),
#         e = d_fit_sus_cens_postno$left_age_e,
#         r = d_fit_sus_cens_postno$right_age_r,
#         sex = d_fit_sus_cens_postno$sex,
#         age2date = sus_cens_postno_age2date,
#         beta_male = beta_male,
#         beta0_sus = beta0_survival_sus,
#         beta0_inf = beta0_survival_inf,
#         age_effect_surv = age_effect_survival_test,
#         period_effect_surv = period_effect_survival_test,
#         f_age_foi = f_age_foi,
#         m_age_foi = m_age_foi,
#         age_lookup_f = age_lookup_f,
#         age_lookup_m = age_lookup_m,
#         period_lookup_foi = period_lookup_foi,
#         f_period_foi = f_period_foi,
#         m_period_foi = m_period_foi,
#         space = c(0,-.55),
#         sect = d_fit_sus_cens_postno$study_area,
#         log = TRUE
#         )
# (end<- Sys.time()-starttime)
# test

# starttime <- Sys.time()
# test2 <- dSusCensNo(
#         x = 1,
# 		n_samples = nrow(d_fit_endlive),
#         e = d_fit_endlive$left_age_e,
#         r = d_fit_endlive$right_age_r,
#         sex = d_fit_endlive$sex,
#         age2date = endlive_age2date,
#         beta_male = beta_male,
#         beta0_sus = beta0_survival_sus,
#         beta0_inf = beta0_survival_inf,
#         age_effect_surv = age_effect_survival_test,
#         period_effect_surv = period_effect_survival_test,
#         f_age_foi = f_age_foi,
#         m_age_foi = m_age_foi,
#         age_lookup_f = age_lookup_f,
#         age_lookup_m = age_lookup_m,
#         period_lookup_foi = period_lookup_foi,
#         f_period_foi = f_period_foi,
#         m_period_foi = m_period_foi,
#         space = c(0,-.55),
# 		sect = d_fit_endlive$study_area,
#         log = TRUE
#         )
# (end<- Sys.time()-starttime)
# test2

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

dSusMortTest <- nimble::nimbleFunction(
    run = function( ### argument type declarations
                   x = double(0),
                   n_samples = integer(0), # number of samples in dataset
                   e = double(1), # e, age of entry
                   r = double(1), # r, age of last known alive
                   s = double(1), # s, age of known mortality
                   sex = double(1),
                   age2date = double(1),
                   beta_male = double(0),
                   beta0_sus = double(0),
                   age_effect_surv = double(1),
                   period_effect_surv = double(1),
                   f_age_foi = double(1),
                   m_age_foi = double(1),
                   age_lookup_f = double(1),
                   age_lookup_m = double(1),
                   period_lookup_foi = double(1),
                   f_period_foi = double(1),
                   m_period_foi = double(1),
                   space = double(1),
                   sect = double(1),
                   log = double()) {
        sumllik <- 0 # intialize log-likelihood
        for (i in 1:n_samples) {
            lam_foi <- 0
            lam_sus <- 0
            lam_susD <- 0

            #############################################
            # preliminary hazards for the likelihood
            #############################################
            if (sex[i] == 0) {
                for (j in r[i]:(s[i] - 1)) {
                    # sum up surv hazard from r to s-1
                    lam_susD <- lam_susD +
					            exp(beta0_sus +
								    age_effect_surv[j] +
									period_effect_surv[age2date[i] + j])
                }

                for (j in 1:(s[i] - 1)) {
                    #  up foi hazard from 1  to j
                    lam_foi <- lam_foi +
					           exp(space[sect[i]] +
							       f_age_foi[age_lookup_f[j]] +
								   f_period_foi[period_lookup_foi[age2date[i] + j]])
                    if (j > (e[i] - 1) & j < r[i]) {
                        lam_sus <- lam_sus +
						           exp(beta0_sus +
								       age_effect_surv[j] +
									   period_effect_surv[age2date[i] + j])
                    }
                }
            } else {
                for (j in r[i]:(s[i] - 1)) {
                    # sum up surv hazard from r to s-1
                    lam_susD <- lam_susD +
					            exp(beta0_sus +
								    age_effect_surv[j] +
									period_effect_surv[age2date[i] + j] +
									beta_male)
                }

                for (j in 1:(s[i] - 1)) {
                    #  up foi hazard from 1  to j
                    lam_foi <- lam_foi +
					           exp(space[sect[i]] +
							       m_age_foi[age_lookup_f[j]] +
								   m_period_foi[period_lookup_foi[age2date[i] + j]])
                    if (j > (e[i] - 1) & j < r[i]) {
                        lam_sus <- lam_sus +
						           exp(beta0_sus +
								       age_effect_surv[j] +
									   period_effect_surv[age2date[i] + j] +
									   beta_male)
                    }
                }
            }
            #######################################
            ### calculating the joint likelihood
            #######################################

            sumllik <- sumllik +
			           log(1 - exp(-lam_susD)) -
					   lam_sus -
					   lam_foi
        }
        returnType(double(0))
        if (log) {
            return(sumllik)
        } else {
            return(exp(sumllik))
        } ## return log-likelihood
    }
)

nimble::registerDistributions(list(
    dSusMortTest = list(
        BUGSdist = "dSusMortTest(n_samples,e,r,s,sex,age2date,beta_male,beta0_sus,age_effect_surv,period_effect_surv,f_age_foi,m_age_foi,age_lookup_f,age_lookup_m,f_period_foi,m_period_foi,period_lookup_foi,space,sect)",
        types = c(
            "value = double(0)",
            "n_samples = integer(0)",
            "e = double(1)",
            "r = double(1)",
            "s = double(1)",
            "sex = double(1)",
            "age2date = double(1)",
            "beta_male = double(0)",
            "beta0_sus = double(0)",
            "age_effect_surv = double(1)",
            "period_effect_surv = double(1)",
            "f_age_foi = double(1)",
            "m_age_foi = double(1)",
            "age_lookup_f = double(1)",
            "age_lookup_m = double(1)",
            "period_lookup_foi=double(1)",
            "f_period_foi=double(1)",
            "m_period_foi=double(1)",
            "space = double(1)",
            "sect = double(1)",
            "log = double(0)"
        ),
        discrete = TRUE
    )
))


# ###for a user-defined distribution
assign("dSusMortTest", dSusMortTest, envir = .GlobalEnv)


# starttime <- Sys.time()
# test <- dSusMortTest(
#         x = 1,
# 		n_samples = nrow(d_fit_sus_mort_posttest),
#         e = d_fit_sus_mort_posttest$left_age_e,
#         r = d_fit_sus_mort_posttest$right_age_r,
#         s = d_fit_sus_mort_posttest$right_age_s,
#         sex = d_fit_sus_mort_posttest$sex,
#         age2date = sus_mort_posttest_age2date,
#         beta_male = beta_male,
#         beta0_sus = beta0_survival_sus,
#         age_effect_surv = age_effect_survival_test,
#         period_effect_surv = period_effect_survival_test,
#         f_age_foi = f_age_foi,
#         m_age_foi = m_age_foi,
#         age_lookup_f = age_lookup_f,
#         age_lookup_m = age_lookup_m,
#         period_lookup_foi = period_lookup_foi,
#         f_period_foi = f_period_foi,
#         m_period_foi = m_period_foi,
#         space = c(0,-.55),
# 		sect = d_fit_sus_mort_posttest$study_area,
#         log = TRUE
#         )
# (end<- Sys.time()-starttime)
# test


#######################################################################
###
###   User defined distribution for likelihood for
###   uninfected radio-marked deer mortalities:
###   test neg at cap and no test at mortality, no recap
###
###   d_fit_sus_mort_postno
###
###   Overleaf Equation (13)
###
#######################################################################

dSusMortNoTest <- nimble::nimbleFunction(
    run = function( ### argument type declarations
                   x = integer(0),
                   n_samples = integer(0), # number of samples in dataset
                   e = double(1), # e, age of entry
                   r = double(1), # r, age of last known alive
                   s = double(1), # s, age of known mortality
                   dn1 = double(1), # right of last test negative
                   sex = double(1),
                   age2date = double(1),
                   beta_male = double(0),
                   beta0_sus = double(0),
                   beta0_inf = double(0),
                   age_effect_surv = double(1),
                   period_effect_surv = double(1),
                   f_age_foi = double(1),
                   m_age_foi = double(1),
                   age_lookup_f = double(1),
                   age_lookup_m = double(1),
                   period_lookup_foi = double(1),
                   f_period_foi = double(1),
                   m_period_foi = double(1),
                   space = double(1),
                   sect = double(1),
                   log = double(0)) {
        sumllik <- 0 # intialize log-likelihood
        for (i in 1:n_samples) {
            lam_foi <- nimNumeric(s[i] - 1, init = FALSE)
            lam_sus <- nimNumeric(s[i] - 1, init = FALSE)
            lam_inf <- nimNumeric(s[i] - 1, init = FALSE)
            lik_temp <- 0

            #############################################
            # preliminary hazards for the likelihood
            #############################################
            if (sex[i] == 0) { # age loops for females
                # survival hazard for susceptible deer
                lam_sus[e[i]:(s[i] - 1)] <- exp(beta0_sus +
                    age_effect_surv[e[i]:(s[i] - 1)] +
                    period_effect_surv[(e[i] + age2date[i]):(s[i] - 1 + age2date[i])])

                # survival hazard while infected
                lam_inf[e[i]:(s[i] - 1)] <- exp(beta0_inf +
                    age_effect_surv[e[i]:(s[i] - 1)] +
                    period_effect_surv[(e[i] + age2date[i]):(s[i] - 1 + age2date[i])])

                # force of infection infection hazard
                lam_foi[1:(s[i] - 1)] <- exp(space[sect[i]] +
                    f_age_foi[age_lookup_f[1:(s[i] - 1)]] +
                    f_period_foi[period_lookup_foi[(1 + age2date[i]):(s[i] - 1 + age2date[i])]])
            } else { # males
                # survival hazard for susceptible deer
                lam_sus[e[i]:(s[i] - 1)] <- exp(beta0_sus +
                    age_effect_surv[e[i]:(s[i] - 1)] +
                    period_effect_surv[(e[i] + age2date[i]):(s[i] - 1 + age2date[i])] +
                    beta_male)

                # survival hazard while infected
                lam_inf[e[i]:(s[i] - 1)] <- exp(beta0_inf +
                    age_effect_surv[e[i]:(s[i] - 1)] +
                    period_effect_surv[(e[i] + age2date[i]):(s[i] - 1 + age2date[i])] +
                    beta_male)

                # force of infection infection hazard
                lam_foi[1:(s[i] - 1)] <- exp(space[sect[i]] +
                    m_age_foi[age_lookup_m[1:(s[i] - 1)]] +
                    m_period_foi[period_lookup_foi[(1 + age2date[i]):(s[i] - 1 + age2date[i])]])
            }
            #######################################
            ### calculating the joint likelihood
            #######################################

            # total probability of getting infected and dying before end of the study
            lik_temp <- lik_temp + lam_foi[dn1[i] + 1] *
                exp(-sum(lam_inf[(dn1[i] + 1):(r[i] - 1)])) *
                (1 - exp(-sum(lam_inf[r[i]:(s[i] - 1)])))
            if ((r[i] - dn1[i]) > 3) {
                for (k in (dn1[i] + 2):(r[i] - 2)) {
                    lik_temp <- lik_temp + lam_foi[k] *
                        exp(-sum(lam_foi[(dn1[i] + 1):(k - 1)])) *
                        exp(-sum(lam_sus[(dn1[i] + 1):(k - 1)])) *
                        exp(-sum(lam_inf[k:(r[i] - 1)])) *
                        (1 - exp(-sum(lam_inf[r[i]:(s[i] - 1)])))
                }
            }
            if ((r[i] - dn1[i]) > 2) {
                for (k in (r[i] - 1):(s[i] - 2)) {
                    lik_temp <- lik_temp + lam_foi[k] *
                        exp(-sum(lam_foi[(dn1[i] + 1):(k - 1)])) *
                        exp(-sum(lam_sus[(dn1[i] + 1):(k - 1)])) *
                        (1 - exp(-sum(lam_inf[k:(s[i] - 1)])))
                }
            }
            lik_temp <- lik_temp + lam_foi[(s[i] - 1)] *
                exp(-sum(lam_foi[(dn1[i] + 1):((s[i] - 1) - 1)])) *
                exp(-sum(lam_sus[(dn1[i] + 1):((s[i] - 1) - 1)])) *
                (lam_inf[(s[i] - 1)])

            sumllik <- sumllik + log(exp(-sum(lam_sus[e[i]:dn1[i]])) *
                exp(-sum(lam_foi[1:dn1[i]])) *
                lik_temp +
                exp(-sum(lam_foi[1:(s[i] - 1)])) *
                    exp(-sum(lam_sus[e[i]:(r[i] - 1)])) *
                    (1 - exp(-sum(lam_sus[r[i]:(s[i] - 1)]))))
        }
        returnType(double(0))
        if (log) {
            return(sumllik)
        } else {
            return(exp(sumllik))
        }
    }
)

nimble::registerDistributions(list(
    dSusMortNoTest = list(
        BUGSdist = "dSusMortNoTest(n_samples,e,r,s,dn1,sex,age2date,beta_male,beta0_sus,beta0_inf,age_effect_surv,period_effect_surv,f_age_foi,m_age_foi,age_lookup_f,age_lookup_m,f_period_foi,m_period_foi,period_lookup_foi,space,sect)",
        types = c(
            "value = integer(0)",
            "n_samples = integer(0)",
            "e = double(1)",
            "r = double(1)",
            "s = double(1)",
            "dn1 = double(1)",
            "sex = double(1)",
            "age2date = double(1)",
            "beta_male = double(0)",
            "beta0_inf = double(0)",
            "beta0_sus = double(0)",
            "age_effect_surv = double(1)",
            "period_effect_surv = double(1)",
            "f_age_foi = double(1)",
            "m_age_foi = double(1)",
            "age_lookup_f = double(1)",
            "age_lookup_m = double(1)",
            "period_lookup_foi=double(1)",
            "f_period_foi=double(1)",
            "m_period_foi=double(1)",
            "space = double(1)",
            "sect = double(1)",
            "log = double(0)"
        ),
        discrete = TRUE
    )
))

### for a user-defined distribution
assign("dSusMortNoTest", dSusMortNoTest, envir = .GlobalEnv)

# starttime <- Sys.time()
# test <-  dSusMortNoTest(
#         x = 1,
#         n_samples = nrow(d_fit_sus_mort_postno),
#         e = d_fit_sus_mort_postno$left_age_e,
#         r = d_fit_sus_mort_postno$right_age_r,
#         s = d_fit_sus_mort_postno$right_age_s,
#         dn1 = d_fit_sus_mort_postno$left_age_e,
#         sex = d_fit_sus_mort_postno$sex,
#         age2date = sus_mort_postno_age2date,
#         beta_male = beta_male,
#         beta0_sus = beta0_survival_sus,
#         beta0_inf = beta0_survival_inf,
#         age_effect_surv = age_effect_survival_test,
#         period_effect_surv = period_effect_survival_test,
#         f_age_foi = f_age_foi,
#         m_age_foi = m_age_foi,
#         age_lookup_f = age_lookup_f,
#         age_lookup_m = age_lookup_m,
#         period_lookup_foi = period_lookup_foi,
#         f_period_foi = f_period_foi,
#         m_period_foi = m_period_foi,
#         sect = d_fit_sus_mort_postno$study_area,
#         space = c(0,-.55),
#         log = TRUE
#         )
# (end<- Sys.time()-starttime)
# test


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

dIcapCens <- nimble::nimbleFunction(
    run = function( ### argument type declarations
                   x = integer(0),
                   n_samples = integer(0), # number of samples in dataset
                   e = double(1), # e, age of entry
                   r = double(1), # r, age of last known alive
                   sex = double(1),
                   age2date = double(1),
                   beta_male = double(0),
                   beta0_sus = double(0),
                   beta0_inf = double(0),
                   age_effect_surv = double(1),
                   period_effect_surv = double(1),
                   f_age_foi = double(1),
                   m_age_foi = double(1),
                   age_lookup_f = double(1),
                   age_lookup_m = double(1),
                   period_lookup_foi = double(1),
                   f_period_foi = double(1),
                   m_period_foi = double(1),
                   space = double(1),
                   sect = double(1),
                   log = double()) {
        sumllik <- 0 # intialize log-likelihood
        for (i in 1:n_samples) {
            lam_foi <- 0
            lam_sus <- 0
            lam_inf <- 0
            lik_temp <- 0

            #############################################
            # preliminary hazards for the likelihood
            #############################################
            if (sex[i] == 0) { # age loops for females
                # lam_inf from 1 to e-1
                for (j in 1:(e[i] - 1)) {
                    lam_inf <- lam_inf +
					           exp(beta0_inf +
							       age_effect_surv[j] +
								   period_effect_surv[age2date[i] + j])
                }
                for (j in 1:(e[i] - 1)) {
                    lam_foij <- exp(space[sect[i]] +
					                f_age_foi[age_lookup_f[j]] +
									f_period_foi[period_lookup_foi[age2date[i] + j]])
                    lik_temp <- lik_temp +
					            lam_foij *
								  exp(-lam_sus) *
								  exp(-lam_inf) *
								  exp(-lam_foi)
                    lam_foi <- lam_foi + lam_foij
                    lam_sus <- lam_sus +
					           exp(beta0_sus +
							       age_effect_surv[j] +
								   period_effect_surv[age2date[i] + j])
                    lam_inf <- lam_inf -
					           exp(beta0_inf +
							       age_effect_surv[j] +
								   period_effect_surv[age2date[i] + j])
                }
                # now recalculate lam_inf to sum over e to r-1
                lam_inf <- 0
                for (j in e[i]:(r[i] - 1)) {
                    lam_inf <- lam_inf +
					           exp(beta0_inf +
							       age_effect_surv[j] +
								   period_effect_surv[age2date[i] + j])
                }
            } else { # age loops for males
                # lam_inf from 1 to e-1
                for (j in 1:(e[i] - 1)) {
                    lam_inf <- lam_inf +
					           exp(beta0_inf +
							       age_effect_surv[j] +
								   period_effect_surv[age2date[i] + j] +
								   beta_male)
                }
                for (j in 1:(e[i] - 1)) {
                    lam_foij <- exp(space[sect[i]] +
					                m_age_foi[age_lookup_m[j]] +
									m_period_foi[period_lookup_foi[age2date[i] + j]])
                    lik_temp <- lik_temp +
					            lam_foij *
								    exp(-lam_sus) *
								    exp(-lam_inf) *
								    exp(-lam_foi)
                    lam_foi <- lam_foi + lam_foij
                    lam_sus <- lam_sus +
					           exp(beta0_sus +
							       age_effect_surv[j] +
								   period_effect_surv[age2date[i] + j] +
								   beta_male)
                    lam_inf <- lam_inf -
					           exp(beta0_inf +
							       age_effect_surv[j] +
								   period_effect_surv[age2date[i] + j] +
								   beta_male)
                }
                # now recalculate lam_inf to sum over e to r - 1
                lam_inf <- 0
                for (j in e[i]:(r[i] - 1)) {
                    lam_inf <- lam_inf +
					           exp(beta0_inf +
							       age_effect_surv[j] +
								   period_effect_surv[age2date[i] + j] +
								   beta_male)
                }
            }
            #######################################
            ### calculating the joint likelihood
            #######################################
            sumllik <- sumllik + log(lik_temp) - lam_inf
        }
        returnType(double(0))
        if (log) {
            return(sumllik)
        } else {
            return(exp(sumllik))
        } ## return log-likelihood
    }
)


nimble::registerDistributions(list(
    dIcapCens = list(
        BUGSdist = "dIcapCens(n_samples,e,r,sex,age2date,beta_male,beta0_inf,beta0_sus,age_effect_surv,period_effect_surv,f_age_foi,m_age_foi,age_lookup_f,age_lookup_m,f_period_foi,m_period_foi,period_lookup_foi,space,sect)",
        types = c(
            "value=integer(0)",
            "n_samples = integer(0)",
            "e = double(1)",
            "r = double(1)",
            "sex = double(1)",
            "age2date = double(1)",
            "beta_male = double(0)",
            "beta0_inf = double(0)",
            "beta0_sus = double(0)",
            "age_effect_surv = double(1)",
            "period_effect_surv = double(1)",
            "f_age_foi = double(1)",
            "m_age_foi = double(1)",
            "age_lookup_f = double(1)",
            "age_lookup_m = double(1)",
            "period_lookup_foi=double(1)",
            "f_period_foi=double(1)",
            "m_period_foi=double(1)",
            "space = double(1)",
            "sect = double(1)",
            "log = double()"
        ),
        discrete = TRUE
    )
))

# for a user-defined distribution
assign("dIcapCens", dIcapCens, envir = .GlobalEnv)

# starttime <- Sys.time()
# test <-  dIcapCens(
#         x = 1,
#         n_samples = nrow(d_fit_icap_cens),
#         e = d_fit_icap_cens$left_age_e,
#         r = d_fit_icap_cens$right_age_r,
#         sex = d_fit_icap_cens$sex,
#         age2date = icap_cens_age2date,
#         beta_male = beta_male,
#         beta0_sus = beta0_survival_sus,
#         beta0_inf = beta0_survival_inf,
#         age_effect_surv = age_effect_survival_test,
#         period_effect_surv = period_effect_survival_test,
#         f_age_foi = f_age_foi,
#         m_age_foi = m_age_foi,
#         age_lookup_f = age_lookup_f,
#         age_lookup_m = age_lookup_m,
#         period_lookup_foi = period_lookup_foi,
#         f_period_foi = f_period_foi,
#         m_period_foi = m_period_foi,
#         space = c(0,-.55),
#         sect = d_fit_icap_cens$study_area,
#         log = TRUE
#         )
# (endtime <- Sys.time() - starttime)
# test


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

dIcapMort <- nimble::nimbleFunction(
    run = function( ### argument type declarations
                   x = integer(0),
                   n_samples = integer(0), # number of samples in dataset
                   e = double(1), # e, age of entry
                   r = double(1), # r, age of last known alive
                   s = double(1), # s, age of known mortality
                   sex = double(1),
                   age2date = double(1),
                   beta_male = double(0),
                   beta0_sus = double(0),
                   beta0_inf = double(0),
                   age_effect_surv = double(1),
                   period_effect_surv = double(1),
                   f_age_foi = double(1),
                   m_age_foi = double(1),
                   age_lookup_f = double(1),
                   age_lookup_m = double(1),
                   period_lookup_foi = double(1),
                   f_period_foi = double(1),
                   m_period_foi = double(1),
                   space = double(1),
                   sect = double(1),
                   log = double()) {

        # start the loop through individuals
        sumllik <- 0
        for (i in 1:n_samples) {
            lam_inf <- 0
            lam_foi <- 0
            lam_sus <- 0
            lik_temp <- 0

            #############################################
            # preliminary hazards for the likelihood
            #############################################
            if (sex[i] == 0) { # females
                # sum up infected hazard from 1 to e-1
                for (j in 1:(e[i] - 1)) {
                    lam_inf <- lam_inf +
                               exp(beta0_inf +
                                   age_effect_surv[j] +
                                   period_effect_surv[age2date[i] + j])
                }
                # now accumulate the other sums up to e-1
                for (j in 1:(e[i] - 1)) {
                    # sum up foi hazard from 1  to j
                    lam_foij <- exp(space[sect[i]] +
                        f_age_foi[age_lookup_f[j]] +
                        f_period_foi[period_lookup_foi[age2date[i] + j]])
                    # sum up like_temp (no sus hazard when j=1)
                    lik_temp <- lik_temp +
                                lam_foij *
                                    exp(-lam_sus) *
                                    exp(-lam_foi) *
                                    exp(-lam_inf)
                    # now increment the hazards (subtract for lam_inf)
                    if (j < (e[i] - 1)) {
                        lam_foi <- lam_foi + lam_foij
                        lam_sus <- lam_sus +
                            exp(beta0_sus +
                                age_effect_surv[j] +
                                period_effect_surv[age2date[i] + j])
                        lam_inf <- lam_inf -
                            exp(beta0_inf +
                                age_effect_surv[j] +
                                period_effect_surv[age2date[i] + j])
                    } else { # now calculate infected hazard for later
                        lam_inf <- sum(exp(beta0_inf +
                            age_effect_surv[e[i]:(r[i] - 1)] +
                            period_effect_surv[age2date[i] + e[i]:(r[i] - 1)]))
                        lam_inf_s <- sum(exp(beta0_inf +
                            age_effect_surv[(r[i]):(s[i] - 1)] +
                            period_effect_surv[age2date[i] + (r[i]):(s[i] - 1)]))
                    }
                }
            } else {
                # sum up infected hazard from 1 to e-1
                for (j in 1:(e[i] - 1)) {
                    lam_inf <- lam_inf +
                               exp(beta0_inf +
                               age_effect_surv[j] +
                               period_effect_surv[age2date[i] + j] +
                               beta_male)
                }
                # now accumulate the other sums up to e-1
                for (j in 1:(e[i] - 1)) {
                    # sum up foi hazard from 1  to j
                    lam_foij <- exp(space[sect[i]] +
                                    m_age_foi[age_lookup_m[j]] +
                                    m_period_foi[period_lookup_foi[age2date[i] + j]])
                    # sum up like_temp (no sus hazard when j=1)
                    lik_temp <- lik_temp + lam_foij *
                        exp(-lam_sus) *
                        exp(-lam_foi) *
                        exp(-lam_inf)
                    # now increment the hazards (subtract for lam_inf)
                    if (j < (e[i] - 1)) {
                        lam_foi <- lam_foi + lam_foij
                        lam_sus <- lam_sus +
                                   exp(beta0_sus +
                                       age_effect_surv[j] +
                                       period_effect_surv[age2date[i] + j] +
                                       beta_male)
                        lam_inf <- lam_inf -
                                    exp(beta0_inf +
                                       age_effect_surv[j] +
                                       period_effect_surv[age2date[i] + j] +
                                       beta_male)
                    } else { # now calculate infected hazard for later
                        lam_inf <- sum(exp(beta0_inf +
                                           age_effect_surv[e[i]:(r[i] - 1)] +
                                           period_effect_surv[age2date[i] + e[i]:(r[i] - 1)] +
                                           beta_male))
                        lam_inf_s <- sum(exp(beta0_inf +
                                             age_effect_surv[r[i]:(s[i] - 1)] +
                                             period_effect_surv[age2date[i] + r[i]:(s[i] - 1)] +
                                             beta_male))
                    }
                }
            }
            #######################################
            ### calculating the joint likelihood
            #######################################

            sumllik <- sumllik -
                       lam_inf +
                       log(1 - exp(-lam_inf_s)) +
                       log(lik_temp)
        }
        returnType(double(0))
        if (log) {
            return(sumllik)
        } else {
            return(exp(sumllik))
        } ## return log-likelihood
    }
)


nimble::registerDistributions(list(
    dIcapMort = list(
        BUGSdist = "dIcapMort(n_samples,e,r,s,sex,age2date,beta_male,beta0_sus,beta0_inf,age_effect_surv,period_effect_surv,f_age_foi,m_age_foi,age_lookup_f,age_lookup_m,f_period_foi,m_period_foi,period_lookup_foi,space,sect)",
        types = c(
            "value = integer(0)",
            "n_samples = integer(0)",
            "e = double(1)",
            "r = double(1)",
            "s = double(1)",
            "sex = double(1)",
            "age2date = double(1)",
            "beta_male = double(0)",
            "beta0_sus = double(0)",
            "beta0_inf = double(0)",
            "age_effect_surv = double(1)",
            "period_effect_surv = double(1)",
            "f_age_foi = double(1)",
            "m_age_foi = double(1)",
            "age_lookup_f = double(1)",
            "age_lookup_m = double(1)",
            "f_period_foi=double(1)",
            "m_period_foi=double(1)",
            "period_lookup_foi=double(1)",
            "space = double(1)",
            "sect = double(1)",
            "log = double()"
        ),
        discrete = TRUE
    )
))

### for a user-defined distribution
assign("dIcapMort", dIcapMort, envir = .GlobalEnv)

# starttime <- Sys.time()
# test <- dIcapMort(
#         x = 1,
#         n_samples = nrow(d_fit_icap_mort),
#         e = d_fit_icap_mort$left_age_e,
#         r = d_fit_icap_mort$right_age_r,
#         s = d_fit_icap_mort$right_age_s,
#         sex = d_fit_icap_mort$sex,
#         age2date = icap_mort_age2date,
#         beta_male = beta_male,
#         beta0_sus = beta0_survival_sus,
#         beta0_inf = beta0_survival_inf,
#         age_effect_surv = age_effect_survival_test,
#         period_effect_surv = period_effect_survival_test,
#         f_age_foi = f_age_foi,
#         m_age_foi = m_age_foi,
#         age_lookup_f = age_lookup_f,
#         age_lookup_m = age_lookup_m,
#         period_lookup_foi = period_lookup_foi,
#         f_period_foi = f_period_foi,
#         m_period_foi = m_period_foi,
#         space = c(0, -.55),
#         sect = d_fit_icap_mort$study_area,
#         log = TRUE
#         )
# (end <- Sys.time() - starttime)
# test

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

dRecNegCensTest <- nimble::nimbleFunction(
    run = function( ### argument type declarations
                   x = integer(0),
                   n_samples = integer(0), # number of samples in dataset
                   e = double(1), # e, age of entry
                   r = double(1), # r, age of last known alive
                   sex = double(1),
                   age2date = double(1),
                   beta_male = double(0),
                   beta0_sus = double(0),
                   age_effect_surv = double(1),
                   period_effect_surv = double(1),
                   f_age_foi = double(1),
                   m_age_foi = double(1),
                   age_lookup_f = double(1),
                   age_lookup_m = double(1),
                   period_lookup_foi = double(1),
                   f_period_foi = double(1),
                   m_period_foi = double(1),
                   sect = double(1),
                   space = double(1),
                   log = double()) {
        sumllik <- 0 # intialize log-likelihood
        for (i in 1:n_samples) {
            lam_foi <- 0
            lam_sus <- 0
            # lam_inf <- 0
            lik_temp <- 0

            #############################################
            # preliminary hazards for the likelihood
            #############################################
            if (sex[i] == 0) { # females
                for (j in e[i]:(r[i] - 1)) {
                    # survival hazard for susceptible deer
                    lam_sus <- lam_sus + exp(beta0_sus +
                        age_effect_surv[j] +
                        period_effect_surv[j + age2date[i]])
                }
                # force of infection infection hazard
                for (j in 1:(r[i] - 1)) {
                    lam_foi <- lam_foi + exp(space[sect[i]] +
                        f_age_foi[age_lookup_f[j]] +
                        f_period_foi[period_lookup_foi[j + age2date[i]]])
                }
            } else { # males
                for (j in e[i]:(r[i] - 1)) {
                    # survival hazard for susceptible deer
                    lam_sus <- lam_sus + exp(beta0_sus +
                        age_effect_surv[j] +
                        period_effect_surv[j + age2date[i]] +
                        beta_male)
                }
                # force of infection infection hazard
                for (j in 1:(r[i] - 1)) {
                    lam_foi <- lam_foi + exp(space[sect[i]] +
                        m_age_foi[age_lookup_m[j]] +
                        m_period_foi[period_lookup_foi[j + age2date[i]]])
                }
            }
            #######################################
            ### calculating the joint likelihood
            #######################################
            sumllik <- sumllik - lam_sus - lam_foi
        }
        returnType(double(0))
        if (log) {
            return(sumllik)
        } else {
            return(exp(sumllik))
        } ## return log-likelihood
    }
)

nimble::registerDistributions(list(
    dRecNegCensTest = list(
        BUGSdist = "dRecNegCensTest(n_samples,e,r,sex,age2date,beta_male,beta0_sus,age_effect_surv,period_effect_surv,f_age_foi,m_age_foi,age_lookup_f,age_lookup_m,f_period_foi,m_period_foi,period_lookup_foi,space,sect)",
        types = c(
            "value = integer(0)",
            "n_samples = integer(0)",
            "e = double(1)",
            "r = double(1)",
            "sex = double(1)",
            "age2date = double(1)",
            "beta_male = double(0)",
            "beta0_sus = double(0)",
            "age_effect_surv = double(1)",
            "period_effect_surv = double(1)",
            "f_age_foi = double(1)",
            "m_age_foi = double(1)",
            "age_lookup_f = double(1)",
            "age_lookup_m = double(1)",
            "f_period_foi=double(1)",
            "m_period_foi=double(1)",
            "period_lookup_foi=double(1)",
            "space = double(1)",
            "sect = double(1)",
            "log = double(0)"
        ),
        discrete = TRUE
    )
))

## for a user-defined distribution
assign("dRecNegCensTest", dRecNegCensTest, envir = .GlobalEnv)

# starttime <- Sys.time()
# test <- dRecNegCensTest(
#         x = 1,
#         n_samples = nrow(d_fit_rec_neg_cens_posttest),
#         e = d_fit_rec_neg_cens_posttest$left_age_e,
#         r = d_fit_rec_neg_cens_posttest$right_age_r,
#         sex = d_fit_rec_neg_cens_posttest$sex,
#         age2date = rec_neg_cens_posttest_age2date,
#         beta_male = beta_male,
#         beta0_sus = beta0_survival_sus,
#         age_effect_surv = age_effect_survival_test,
#         period_effect_surv = period_effect_survival_test,
#         f_age_foi = f_age_foi,
#         m_age_foi = m_age_foi,
#         age_lookup_f = age_lookup_f,
#         age_lookup_m = age_lookup_m,
#         period_lookup_foi = period_lookup_foi,
#         f_period_foi = f_period_foi,
#         m_period_foi = m_period_foi,
#         sect = d_fit_rec_neg_cens_posttest$study_area,
#         space = c(0,-.55),
#         log = TRUE
#         )
# (end<- Sys.time()-starttime)
# test


#######################################################################
###
###   User defined distribution for likelihood for
###   uninfected deer that were test neg at capture,
###   then test negative at recap, that are right censored,
###   and have not been tested post censoring
###
###   d_fit_rec_neg_cens_postno
###
###   Overleaf Equation (21)
###
#######################################################################

dRecNegCensPostNo <- nimble::nimbleFunction(
    run = function( ### argument type declarations
                   x = integer(0),
                   n_samples = integer(0), # number of samples in dataset
                   e = double(1), # e, age of entry
                   r = double(1), # r, age of last known alive
                   dn1 = double(1), # right of last test negative
                   sex = double(1),
                   age2date = double(1),
                   beta_male = double(0),
                   beta0_sus = double(0),
                   beta0_inf = double(0),
                   age_effect_surv = double(1),
                   period_effect_surv = double(1),
                   f_age_foi = double(1),
                   m_age_foi = double(1),
                   age_lookup_f = double(1),
                   age_lookup_m = double(1),
                   period_lookup_foi = double(1),
                   f_period_foi = double(1),
                   m_period_foi = double(1),
                   sect = double(1),
                   space = double(1),
                   log = double()) {
        sumllik <- 0 # intialize log-likelihood
        for (i in 1:n_samples) {
            lam_foi <- nimNumeric(r[i], init = FALSE)
            lam_sus <- nimNumeric(r[i], init = FALSE)
            lam_inf <- nimNumeric(r[i], init = FALSE)
            # lik_temp <- 0

            #############################################
            # preliminary hazards for the likelihood
            #############################################
            if (sex[i] == 0) { # females
                # survival hazard for susceptible deer
                lam_sus[e[i]:(r[i] - 1)] <- exp(beta0_sus +
                    age_effect_surv[e[i]:(r[i] - 1)] +
                    period_effect_surv[(e[i] + age2date[i]):(r[i] - 1 + age2date[i])])
                # survival hazard while infected
                lam_inf[(dn1[i] + 1):(r[i] - 1)] <- exp(beta0_inf +
                    age_effect_surv[(dn1[i] + 1):(r[i] - 1)] +
                    period_effect_surv[(dn1[i] + 1 + age2date[i]):(r[i] - 1 + age2date[i])])
                # force of infection infection hazard
                lam_foi[1:(r[i] - 1)] <- exp(space[sect[i]] +
                    f_age_foi[age_lookup_f[1:(r[i] - 1)]] +
                    f_period_foi[period_lookup_foi[(1 + age2date[i]):(r[i] - 1 + age2date[i])]])
            } else {
                # survival hazard for susceptible deer
                lam_sus[e[i]:(r[i] - 1)] <- exp(beta0_sus +
                    age_effect_surv[e[i]:(r[i] - 1)] +
                    period_effect_surv[(e[i] + age2date[i]):(r[i] - 1 + age2date[i])] +
                    beta_male)
                # survival hazard while infected
                lam_inf[(dn1[i] + 1):(r[i] - 1)] <- exp(beta0_inf +
                    age_effect_surv[(dn1[i] + 1):(r[i] - 1)] +
                    period_effect_surv[(dn1[i] + 1 + age2date[i]):(r[i] - 1 + age2date[i])] +
                    beta_male)
                # force of infection infection hazard
                lam_foi[1:(r[i] - 1)] <- exp(space[sect[i]] +
                    m_age_foi[age_lookup_m[1:(r[i] - 1)]] +
                    m_period_foi[period_lookup_foi[(1 + age2date[i]):(r[i] - 1 + age2date[i])]])
            }

            #######################################
            ### calculating the joint likelihood
            #######################################

            lik_temp <- lam_foi[(dn1[i] + 1)] *
						exp(-sum(lam_inf[(dn1[i] + 1):(r[i] - 1)]))

            for (k in (dn1[i] + 2):(r[i] - 1)) {
                lik_temp <- lik_temp +
				            lam_foi[k] *
                            exp(-sum(lam_sus[dn1[i]:(k - 1)])) *
                            exp(-sum(lam_inf[k:(r[i] - 1)]))
            }

            sumllik <- sumllik +
			           log(exp(-sum(lam_sus[e[i]:dn1[i]])) *
                           exp(-sum(lam_foi[1:dn1[i]])) *
                           lik_temp +
                       exp(-sum(lam_foi[1:(r[i] - 1)])) *
                           exp(-sum(lam_sus[e[i]:(r[i] - 1)])))
        }
        returnType(double(0))
        if (log) {
            return(sumllik)
        } else {
            return(exp(sumllik))
        } ## return log-likelihood
    }
)

nimble::registerDistributions(list(
    dRecNegCensPostNo = list(
        BUGSdist = "dRecNegCensPostNo(n_samples,e,r,dn1,sex,age2date,beta_male,beta0_inf,beta0_sus,age_effect_surv,period_effect_surv,f_age_foi,m_age_foi,age_lookup_f,age_lookup_m,f_period_foi,m_period_foi,period_lookup_foi,sect,space)",
        types = c(
            "value = integer(0)",
            "n_samples = integer(0)",
            "e = double(1)",
            "r = double(1)",
            "dn1 = double(1)",
            "sex = double(1)",
            "age2date = double(1)",
            "beta_male = double(0)",
            "beta0_inf = double(0)",
            "beta0_sus = double(0)",
            "age_effect_surv = double(1)",
            "period_effect_surv = double(1)",
            "f_age_foi = double(1)",
            "m_age_foi = double(1)",
            "age_lookup_f = double(1)",
            "age_lookup_m = double(1)",
            "period_lookup_foi=double(1)",
            "f_period_foi=double(1)",
            "m_period_foi=double(1)",
            "space = double(0)",
            "log = double(0)"
        ),
        discrete = TRUE
    )
))

### for a user-defined distribution
assign("dRecNegCensPostNo", dRecNegCensPostNo, envir = .GlobalEnv)

# starttime <- Sys.time()
# test <-  dRecNegCensPostNo(
#         x = 1,
#         n_samples = nrow(d_fit_rec_neg_cens_postno),
#         e = d_fit_rec_neg_cens_postno$left_age_e,
#         r = d_fit_rec_neg_cens_postno$right_age_r,
#         dn1 = d_fit_rec_neg_cens_postno$ageweek_recap,
#         sex = d_fit_rec_neg_cens_postno$sex,
#         age2date = rec_neg_cens_postno_age2date,
#         beta_male = beta_male,
#         beta0_sus = beta0_survival_sus,
#         beta0_inf = beta0_survival_inf,
#         age_effect_surv = age_effect_survival_test,
#         period_effect_surv = period_effect_survival_test,
#         f_age_foi = f_age_foi,
#         m_age_foi = m_age_foi,
#         age_lookup_f = age_lookup_f,
#         age_lookup_m = age_lookup_m,
#         period_lookup_foi = period_lookup_foi,
#         f_period_foi = f_period_foi,
#         m_period_foi = m_period_foi,
#         sect = d_fit_rec_neg_cens_postno$study_area,
#         space = c(0,-.55),
#         log = TRUE
#         )
# (end<- Sys.time()-starttime)
# test


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

dRecNegMort <- nimble::nimbleFunction(
    run = function( ### argument type declarations
                   x = integer(0),
                   n_samples = integer(0), # number of samples in dataset
                   e = double(1), # e, age of entry
                   r = double(1), # r, age of last known alive
                   s = double(1), # s, age of known mortality
                   sex = double(1),
                   age2date = double(1),
                   beta_male = double(0),
                   beta0_sus = double(0),
                   age_effect_surv = double(1),
                   period_effect_surv = double(1),
                   f_age_foi = double(1),
                   m_age_foi = double(1),
                   age_lookup_f = double(1),
                   age_lookup_m = double(1),
                   period_lookup_foi = double(1),
                   f_period_foi = double(1),
                   m_period_foi = double(1),
                   sect = double(1),
                   space = double(1),
                   log = double(0)) {
        sumllik <- 0 # intialize log-likelihood
        for (i in 1:n_samples) {
            lam_foi <- 0
            lam_sus <- 0

            #############################################
            # preliminary hazards for the likelihood
            #############################################
            if (sex[i] == 0) { # age loops for females
                # survival hazard for susceptible deer
                for (j in e[i]:(r[i] - 1)) {
                    lam_sus <- lam_sus + exp(beta0_sus +
                        age_effect_surv[j] +
                        period_effect_surv[(j + age2date[i])])
                }
                lam_sus_s <- sum(exp(beta0_sus +
                    age_effect_surv[r[i]:(s[i] - 1)] +
                    period_effect_surv[(r[i] + age2date[i]):(s[i] - 1 + age2date[i])]))
                # force of infection infection hazard
                for (j in 1:(s[i] - 1)) {
                    lam_foi <- lam_foi + exp(space[sect[i]] +
                        f_age_foi[age_lookup_f[j]] +
                        f_period_foi[period_lookup_foi[(j + age2date[i])]])
                }
            } else { # males
                # survival hazard for susceptible deer
                for (j in e[i]:(r[i] - 1)) {
                    lam_sus <- lam_sus + exp(beta0_sus +
                        age_effect_surv[j] +
                        period_effect_surv[(j + age2date[i])] +
                        beta_male)
                }
                lam_sus_s <- sum(exp(beta0_sus +
                    age_effect_surv[r[i]:(s[i] - 1)] +
                    period_effect_surv[(r[i] + age2date[i]):(s[i] - 1 + age2date[i])] +
                        beta_male))

                # force of infection infection hazard
                for (j in 1:(s[i] - 1)) {
                    lam_foi <- lam_foi + exp(space[sect[i]] +
                        m_age_foi[age_lookup_m[j]] +
                        m_period_foi[period_lookup_foi[(j + age2date[i])]])
                }
            }

            #######################################
            ### calculating the joint likelihood
            #######################################

            sumllik <- sumllik - lam_sus + log(1 - exp(-lam_sus_s)) - lam_foi
        }
        returnType(double(0))
        if (log) {
            return(sumllik)
        } else {
            return(exp(sumllik))
        } ## return log-likelihood
    }
)

nimble::registerDistributions(list(
    dRecNegMort = list(
        BUGSdist = "dRecNegMort(n_samples,e,r,s,sex,age2date,beta_male,beta0_sus,age_effect_surv,period_effect_surv,f_age_foi,m_age_foi,age_lookup_f,age_lookup_m,f_period_foi,m_period_foi,period_lookup_foi,space,sect)",
        types = c(
            "value = integer(0)",
            "n_samples = integer(0)",
            "e = double(1)",
            "r = double(1)",
            "s = double(1)",
            "sex =  double(1)",
            "age2date =  double(1)",
            "beta_male = double(0)",
            "beta0_sus = double(0)",
            "age_effect_surv = double(1)",
            "period_effect_surv = double(1)",
            "f_age_foi = double(1)",
            "m_age_foi = double(1)",
            "age_lookup_f = double(1)",
            "age_lookup_m = double(1)",
            "period_lookup_foi = double(1)",
            "f_period_foi = double(1)",
            "m_period_foi = double(1)",
            "space = double(1)",
            "sect = double(1)",
            "log = double(0)"
        ),
        discrete = TRUE
    )
))

# ###for a user-defined distribution
assign("dRecNegMort", dRecNegMort, envir = .GlobalEnv)

# starttime <- Sys.time()
# test <-  dRecNegMort(
#         x = 1,
#         n_samples = nrow(d_fit_rec_neg_mort),
#         e = d_fit_rec_neg_mort$left_age_e,
#         r = d_fit_rec_neg_mort$right_age_r,
#         s = d_fit_rec_neg_mort$right_age_s,
#         sex = d_fit_rec_neg_mort$sex,
#         age2date = rec_neg_mort_age2date,
#         beta_male = beta_male,
#         beta0_sus = beta0_survival_sus,
#         age_effect_surv = age_effect_survival_test,
#         period_effect_surv = period_effect_survival_test,
#         f_age_foi = f_age_foi,
#         m_age_foi = m_age_foi,
#         age_lookup_f = age_lookup_f,
#         age_lookup_m = age_lookup_m,
#         period_lookup_foi = period_lookup_foi,
#         f_period_foi = f_period_foi,
#         m_period_foi = m_period_foi,
#         sect = d_fit_rec_neg_mort$study_area,
#         space = c(0,-.55),
#         log = TRUE
#         )
# (end<- Sys.time()-starttime)
# test


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

dRecPosMort <- nimble::nimbleFunction(
    run = function( ### argument type declarations
                   x = integer(0),
                   n_samples = integer(0), # number of samples in dataset
                   e = double(1), # e, age of entry
                   r = double(1), # r, age of last known alive
                   s = double(1), # s, age of known mortality
                   dn1 = double(1), # last interval of test negative
                   dn = double(1), # interval of test positive (could be interval of mortality)
                   sex = double(1),
                   age2date = double(1),
                   beta_male = double(0),
                   beta0_sus = double(0),
                   beta0_inf = double(0),
                   age_effect_surv = double(1),
                   period_effect_surv = double(1),
                   f_age_foi = double(1),
                   m_age_foi = double(1),
                   age_lookup_f = double(1),
                   age_lookup_m = double(1),
                   period_lookup_foi = double(1),
                   f_period_foi = double(1),
                   m_period_foi = double(1),
                   space = double(1),
                   sect = double(1),
                   log = double()) {
        sumllik <- 0 # intialize log-likelihood
        for (i in 1:n_samples) {
            lam_foi <- nimNumeric(s[i] - 1, init = FALSE)
            lam_sus <- nimNumeric(s[i] - 1, init = FALSE)
            lam_inf <- nimNumeric(s[i] - 1, init = FALSE)
            lik_temp <- 0

            #############################################
            # preliminary hazards for the likelihood
            #############################################
            if (sex[i] == 0) { # females
                # survival hazard for susceptible deer
                lam_sus[e[i]:(dn[i] - 1)] <- exp(beta0_sus +
                    age_effect_surv[e[i]:(dn[i] - 1)] +
                    period_effect_surv[(e[i] + age2date[i]):
									   (dn[i] - 1 + age2date[i])])
                # survival hazard while infected
                lam_inf[dn1[i]:(s[i] - 1)] <- exp(beta0_inf +
                    age_effect_surv[dn1[i]:(s[i] - 1)] +
                    period_effect_surv[(dn1[i] + age2date[i]):
                    (s[i] - 1 + age2date[i])] +
                    beta_male)
                # force of infection infection hazard
                lam_foi[1:(dn[i] - 1)] <- exp(space[sect[i]] +
                    f_age_foi[age_lookup_f[1:(dn[i] - 1)]] +
                    f_period_foi[period_lookup_foi[(1 + age2date[i]):
					                           ((dn[i] - 1) + age2date[i])]])
            } else { # males
                # survival hazard for susceptible deer
                lam_sus[e[i]:(dn[i] - 1)] <- exp(beta0_sus +
                    age_effect_surv[e[i]:(dn[i] - 1)] +
                    period_effect_surv[(e[i] + age2date[i]):
					                   (dn[i] - 1 + age2date[i])] +
                    beta_male)
                # survival hazard while infected
                lam_inf[dn1[i]:(s[i] - 1)] <- exp(beta0_inf +
                    age_effect_surv[dn1[i]:(s[i] - 1)] +
                    period_effect_surv[(dn1[i] + age2date[i]):
                    				   (s[i] - 1 + age2date[i])] +
                    beta_male)
                # force of infection infection hazard
                lam_foi[1:(dn[i] - 1)] <- exp(space[sect[i]] +
                    m_age_foi[age_lookup_m[1:(dn[i] - 1)]] +
                    m_period_foi[period_lookup_foi[(1 + age2date[i]):
					                           ((dn[i] - 1) + age2date[i])]])
            }

            #######################################
            ### calculating the joint likelihood
            #######################################
            lik_temp <- lik_temp + lam_foi[(dn1[i] + 1)] *
                exp(-sum(lam_foi[1:dn1[i]])) *
                exp(-sum(lam_inf[(dn1[i] + 1):(r[i] - 1)]))

            for (k in (dn1[i] + 2):(dn[i] - 1)) {
                lik_temp <- lik_temp + lam_foi[k] *
                    exp(-sum(lam_sus[e[i]:(k - 1)])) *
                    exp(-sum(lam_foi[1:(k - 1)])) *
                    exp(-sum(lam_inf[k:(r[i] - 1)]))
            }

            sumllik <- sumllik + log(1 - exp(-sum(lam_inf[r[i]:(s[i] - 1)]))) +
                log(lik_temp)
        }
        returnType(double(0))
        if (log) {
            return(sumllik)
        } else {
            return(exp(sumllik))
        } ## return log-likelihood
    }
)

nimble::registerDistributions(list(
    dRecPosMort = list(
        BUGSdist = "dRecPosMort(n_samples,e,r,s,dn1,dn,sex,age2date,beta_male,beta0_sus,beta0_inf,age_effect_surv,period_effect_surv,f_age_foi,m_age_foi,age_lookup_f,age_lookup_m,f_period_foi,m_period_foi,period_lookup_foi,space,sect)",
        types = c(
            "value = integer(0)",
            "n_samples = integer(0)",
            "e = double(1)",
            "r = double(1)",
            "s = double(1)",
            "dn1 = double(1)",
            "dn = double(1)",
            "sex = double(1)",
            "age2date = double(1)",
            "beta_male = double(0)",
            "beta0_sus = double(0)",
            "beta0_inf = double(0)",
            "age_effect_surv = double(1)",
            "period_effect_surv = double(1)",
            "f_age_foi = double(1)",
            "m_age_foi = double(1)",
            "age_lookup_f = double(1)",
            "age_lookup_m = double(1)",
            "period_lookup_foi=double(1)",
            "f_period_foi=double(1)",
            "m_period_foi=double(1)",
            "space = double(1)",
            "sect = double(1)",
            "log = double(0)"
        ),
        discrete = TRUE
    )
))

### Global Declaration so Nimble can access
assign("dRecPosMort", dRecPosMort, envir = .GlobalEnv)

# starttime <- Sys.time()
# test <-  dRecPosMort(
#         x = 1,
#         n_samples = nrow(d_fit_rec_pos_mort),
#         e = d_fit_rec_pos_mort$left_age_e,
#         r = d_fit_rec_pos_mort$right_age_r,
#         s = d_fit_rec_pos_mort$right_age_s,
#         dn1 = d_fit_rec_pos_mort$left_age_e,
#         dn = d_fit_rec_pos_mort$ageweek_recap,
#         sex = d_fit_rec_pos_mort$sex,
#         age2date = rec_pos_mort_age2date,
#         beta_male = beta_male,
#         beta0_sus = beta0_survival_sus,
#         beta0_inf = beta0_survival_inf,
#         age_effect_surv = age_effect_survival_test,
#         period_effect_surv = period_effect_survival_test,
#         f_age_foi = f_age_foi,
#         m_age_foi = m_age_foi,
#         age_lookup_f = age_lookup_f,
#         age_lookup_m = age_lookup_m,
#         period_lookup_foi = period_lookup_foi,
#         f_period_foi = f_period_foi,
#         m_period_foi = m_period_foi,
#         sect = d_fit_rec_pos_mort$study_area,
#         space = c(0,-.55),
#         log = TRUE
#         )
# (end<- Sys.time()-starttime)
# test


#######################################################################
###
###   User defined distribution for likelihood for
###   infected deer that were test neg at capture,
###   then test positive at recap,
###   than these were right censored
###
###   d_fit_rec_pos_cens
###
###   Overleaf Equation (27)
###
#######################################################################

dRecPosCens <- nimble::nimbleFunction(
    run = function( ### argument type declarations
                   x = integer(0),
                   # n_samples = integer(0), # number of samples in dataset
                   e = double(0), # e, age of entry
                   r = double(0), # r, age of last known alive
                   dn1 = double(0), # date of last test negative
                   dn = double(0), # date of positive test
                   sex = double(0),
                   age2date = double(0),
                   beta_male = double(0),
                   beta0_sus = double(0),
                   beta0_inf = double(0),
                   age_effect_surv = double(1),
                   period_effect_surv = double(1),
                   f_age_foi = double(1),
                   m_age_foi = double(1),
                   age_lookup_f = double(1),
                   age_lookup_m = double(1),
                   period_lookup_foi = double(1),
                   f_period_foi = double(1),
                   m_period_foi = double(1),
                   space = double(0),
                   # sect = double(0),
                   log = double()) {
        sumllik <- 0 # intialize log-likelihood
        lam_foi <- nimNumeric(r - 1, init = FALSE)
        lam_sus <- nimNumeric(r - 1, init = FALSE)
        lam_inf <- nimNumeric(r - 1, init = FALSE)
        lik_temp <- 0

        #############################################
        # preliminary hazards for the likelihood
        #############################################
        if (sex == 1) { # females
            # survival hazard for susceptible deer
            lam_sus[e:(dn - 1)] <- exp(beta0_sus +
                age_effect_surv[e:(dn - 1)] +
                period_effect_surv[(e + age2date):(dn - 1 + age2date)])

            # survival hazard while infected
            lam_inf[dn1:(r - 1)] <- exp(beta0_inf +
                age_effect_surv[dn1:(r - 1)] +
                period_effect_surv[(dn1 + age2date):
                				   (r - 1 + age2date)])
            # force of infection infection hazard
            lam_foi[1:(dn - 1)] <- exp(space +
                f_age_foi[age_lookup_f[1:(dn - 1)]] +
                f_period_foi[period_lookup_foi[(1 + age2date):
				                           ((dn - 1) + age2date)]])
        } else { # males
            # survival hazard for susceptible deer
            lam_sus[e:(dn - 1)] <- exp(beta0_sus +
                age_effect_surv[e:(dn - 1)] +
                period_effect_surv[(e + age2date):(dn - 1 + age2date)] +
                beta_male)
            # survival hazard while infected
            lam_inf[dn1:(r - 1)] <- exp(beta0_inf +
                age_effect_surv[dn1:(r - 1)] +
                period_effect_surv[(dn1 + age2date):
                (r - 1 + age2date)] +
                beta_male)
            # force of infection infection hazard
            lam_foi[1:(dn - 1)] <- exp(space +
                m_age_foi[age_lookup_m[1:(dn - 1)]] +
                m_period_foi[period_lookup_foi[(1 + age2date):
				                           ((dn - 1) + age2date)]])
        }

        #######################################
        ### calculating the joint likelihood
        #######################################

        lik_temp <- lik_temp + lam_foi[(dn1 + 1)] *
            exp(-sum(lam_foi[1:dn1])) *
            exp(-sum(lam_inf[(dn1 + 1):(r - 1)]))

        for (k in (dn1 + 2):(dn - 1)) {
            lik_temp <- lik_temp + lam_foi[k] *
                exp(-sum(lam_sus[(dn1 + 1):(k - 1)])) *
                exp(-sum(lam_foi[1:(k - 1)])) *
                exp(-sum(lam_inf[k:(r - 1)]))
        }
        sumllik <- sumllik + log(lik_temp) -
            sum(lam_sus[e:dn1]) -
            sum(lam_foi[1:dn1])

        returnType(double(0))
        if (log) {
            return(sumllik)
        } else {
            return(exp(sumllik))
        } ## return log-likelihood
    }
)

nimble::registerDistributions(list(
    dRecPosCens = list(
        BUGSdist = "dRecPosCens(e,r,dn1,dn,sex,age2date,beta_male,beta0_sus,beta0_inf,age_effect_surv,period_effect_surv,f_age_foi,m_age_foi,age_lookup_f,age_lookup_m,f_period_foi,m_period_foi,period_lookup_foi,space)",
        types = c(
            "value = integer(0)",
            "e = double(0)",
            "r = double(0)",
            "dn1 = double(0)",
            "dn = double(0)",
            "sex = double(0)",
            "age2date = double(0)",
            "beta_male = double(0)",
            "beta0_inf = double(0)",
            "beta0_sus = double(0)",
            "age_effect_surv = double(1)",
            "period_effect_surv = double(1)",
            "f_age_foi = double(1)",
            "m_age_foi = double(1)",
            "age_lookup_f = double(1)",
            "age_lookup_m = double(1)",
            "period_lookup_foi=double(1)",
            "f_period_foi=double(1)",
            "m_period_foi=double(1)",
            "space = double(0)",
            #   "sect = double(0)",
            "log = double(0)"
        ),
        discrete = TRUE
    )
))

### for a user-defined distribution
assign("dRecPosCens", dRecPosCens, envir = .GlobalEnv)

# starttime <- Sys.time()
# test <-  dRecPosCens(
#         x = 1,
#         #n_samples = nrow(d_fit_rec_pos_cens),
#         e = d_fit_rec_pos_cens$left_age_e,
#         r = d_fit_rec_pos_cens$right_age_r,
#         dn1 = d_fit_rec_pos_cens$left_age_e,
#         dn = d_fit_rec_pos_cens$ageweek_recap,
#         sex = d_fit_rec_pos_cens$sex,
#         age2date = rec_pos_cens_age2date,
#         beta_male = beta_male,
#         beta0_sus = beta0_survival_sus,
#         beta0_inf = beta0_survival_inf,
#         age_effect_surv = age_effect_survival_test,
#         period_effect_surv = period_effect_survival_test,
#         f_age_foi = f_age_foi,
#         m_age_foi = m_age_foi,
#         age_lookup_f = age_lookup_f,
#         age_lookup_m = age_lookup_m,
#         period_lookup_foi = period_lookup_foi,
#         f_period_foi = f_period_foi,
#         m_period_foi = m_period_foi,
#         space = 0,#c(0,-.55),
#         # sect = d_fit_rec_pos_cens$study_area,
#         log = TRUE
#         )
# (end<- Sys.time()-starttime)
# test

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

dNegCapPosMort <- nimble::nimbleFunction(
    run = function( ### argument type declarations
                   x = integer(0),
                   n_samples = integer(0), # number of samples in dataset
                   e = double(1), # e, age of entry
                   r = double(1), # r, age of last known alive
                   s = double(1), # s, age of known mortality
                   dn1 = double(1), # last interval of test negative
                   dn = double(1), # interval of test positive (could be interval of mortality)
                   sex = double(1),
                   age2date = double(1),
                   beta_male = double(0),
                   beta0_sus = double(0),
                   beta0_inf = double(0),
                   age_effect_surv = double(1),
                   period_effect_surv = double(1),
                   f_age_foi = double(1),
                   m_age_foi = double(1),
                   age_lookup_f = double(1),
                   age_lookup_m = double(1),
                   period_lookup_foi = double(1),
                   f_period_foi = double(1),
                   m_period_foi = double(1),
                   space = double(1),
                   sect = double(1),
                   log = double()) {
        sumllik <- 0 # intialize log-likelihood
        for (i in 1:n_samples) {
            lam_foi <- nimNumeric(s[i] - 1, init = FALSE)
            lam_sus <- nimNumeric(s[i] - 1, init = FALSE)
            lam_inf <- nimNumeric(s[i] - 1, init = FALSE)
            lik_temp <- 0

            #############################################
            # preliminary hazards for the likelihood
            #############################################
            if (sex[i] == 0) { # females
                # survival hazard for susceptible deer
                lam_sus[e[i]:(dn[i] - 1)] <- exp(beta0_sus +
                    age_effect_surv[e[i]:(dn[i] - 1)] +
                    period_effect_surv[(e[i] + age2date[i]):
					                   (dn[i] - 1 + age2date[i])])
                # survival hazard while infected
                lam_inf[dn1[i]:(s[i] - 1)] <- exp(beta0_inf +
                    age_effect_surv[dn1[i]:(s[i] - 1)] +
                    period_effect_surv[(dn1[i] + age2date[i]):
					                   (s[i] - 1 + age2date[i])])
                # force of infection infection hazard
                lam_foi[1:(dn[i] - 1)] <- exp(space[sect[i]] +
                    f_age_foi[age_lookup_f[1:(dn[i] - 1)]] +
                    f_period_foi[period_lookup_foi[(1 + age2date[i]):
					                           (dn[i] - 1 + age2date[i])]])
            } else { # males
                # survival hazard for susceptible deer
                lam_sus[e[i]:(dn[i] - 1)] <- exp(beta0_sus +
                    age_effect_surv[e[i]:(dn[i] - 1)] +
                    period_effect_surv[(e[i] + age2date[i]):
					                   (dn[i] - 1 + age2date[i])] +
                    beta_male)
                # survival hazard while infected
                lam_inf[dn1[i]:(s[i] - 1)] <- exp(beta0_inf +
                    age_effect_surv[dn1[i]:(s[i] - 1)] +
                    period_effect_surv[(dn1[i] + age2date[i]):
					                   (s[i] - 1 + age2date[i])] +
                    beta_male)
                # force of infection infection hazard
                lam_foi[1:(dn[i] - 1)] <- exp(space[sect[i]] +
                    m_age_foi[age_lookup_m[1:(dn[i] - 1)]] +
                    m_period_foi[period_lookup_foi[(1 + age2date[i]):
					                           (dn[i] - 1 + age2date[i])]])
            }

            #######################################
            ### calculating the joint likelihood
            #######################################

            if (r[i] - dn1[i] > 1) {
                lik_temp <- lik_temp + lam_foi[dn1[i] + 1] *
                    exp(-sum(lam_inf[(dn1[i] + 1):(r[i] - 1)])) *
                    (1 - exp(-sum(lam_inf[r[i]:(s[i] - 1)])))
                lik_temp <- lik_temp + lam_foi[(s[i] - 1)] *
                    exp(-sum(lam_foi[(dn1[i] + 1):((s[i] - 1) - 1)])) *
                    exp(-sum(lam_sus[(dn1[i] + 1):((s[i] - 1) - 1)])) *
                    lam_inf[(s[i] - 1)]
            }

            if (r[i] - dn1[i] == 1) {
                lik_temp <- lik_temp + lam_foi[(s[i] - 1)] *
                    lam_inf[(s[i] - 1)]
            }

            if (r[i] - dn1[i] > 3) {
                for (k in (dn1[i] + 2):(r[i] - 2)) {
                    lik_temp <- lik_temp + lam_foi[k] *
                        exp(-sum(lam_foi[(dn1[i] + 1):(k - 1)])) *
                        exp(-sum(lam_sus[(dn1[i] + 1):(k - 1)])) *
                        exp(-sum(lam_inf[k:(r[i] - 1)])) *
                        (1 - exp(-sum(lam_inf[r[i]:(s[i] - 1)])))
                }
            }
            if (r[i] - dn1[i] > 2) {
                for (k in (r[i] - 1):(s[i] - 2)) {
                    lik_temp <- lik_temp + lam_foi[k] *
                        exp(-sum(lam_foi[(dn1[i] + 1):(k - 1)])) *
                        exp(-sum(lam_sus[(dn1[i] + 1):(k - 1)])) *
                        (1 - exp(-sum(lam_inf[k:(s[i] - 1)])))
                }
            }

            sumllik <- sumllik - sum(lam_sus[e[i]:dn1[i]]) -
                sum(lam_foi[1:dn1[i]]) +
                log(lik_temp)
        }
        returnType(double(0))
        if (log) {
            return(sumllik)
        } else {
            return(exp(sumllik))
        } ## return log-likelihood
    }
)


nimble::registerDistributions(list(
    dNegCapPosMort = list(
        BUGSdist = "dNegCapPosMort(n_samples,e,r,s,dn1,dn,sex,age2date,beta_male,beta0_sus,beta0_inf,age_effect_surv,period_effect_surv,f_age_foi,m_age_foi,age_lookup_f,age_lookup_m,f_period_foi,m_period_foi,period_lookup_foi,space,sect)",
        types = c(
            "value = integer(0)",
            "n_samples = integer(0)",
            "e = double(1)",
            "r = double(1)",
            "s = double(1)",
            "dn1 = double(1)",
            "dn = double(1)",
            "sex = double(1)",
            "age2date = double(1)",
            "beta_male = double(0)",
            "beta0_sus = double(0)",
            "beta0_inf = double(0)",
            "age_effect_surv = double(1)",
            "period_effect_surv = double(1)",
            "f_age_foi = double(1)",
            "m_age_foi = double(1)",
            "age_lookup_f = double(1)",
            "age_lookup_m = double(1)",
            "period_lookup_foi = double(1)",
            "f_period_foi = double(1)",
            "m_period_foi = double(1)",
            "space = double(1)",
            "sect = double(1)",
            "log = double()"
        ),
        discrete = TRUE
    )
))

### for a user-defined distribution
assign("dNegCapPosMort", dNegCapPosMort, envir = .GlobalEnv)

# starttime <- Sys.time()
# test <-  dNegCapPosMort(
#         x = 1,
#         n_samples = nrow(d_fit_idead),
#         e = d_fit_idead$left_age_e,
#         r = d_fit_idead$right_age_r,
#         s = d_fit_idead$right_age_s,
#         dn1 = d_fit_idead$left_age_e,
#         dn = d_fit_idead$right_age_s,
#         sex = d_fit_idead$sex,
#         age2date = idead_age2date,
#         beta_male = beta_male,
#         beta0_sus = beta0_survival_sus,
#         beta0_inf = beta0_survival_inf,
#         age_effect_surv = age_effect_survival_test,
#         period_effect_surv = period_effect_survival_test,
#         f_age_foi = f_age_foi,
#         m_age_foi = m_age_foi,
#         age_lookup_f = age_lookup_f,
#         age_lookup_m = age_lookup_m,
#         period_lookup_foi = period_lookup_foi,
#         f_period_foi = f_period_foi,
#         m_period_foi = m_period_foi,
#         space = c(0,-.55),
#         sect = d_fit_idead$study_area,
#         log = TRUE
#         )
# (end<- Sys.time()-starttime)
# test

#######################################################################
###
###   User defined distribution for likelihood for
###   Age-At-Harvest based period effects
###
###   Overleaf Equation
###
#######################################################################

dAAH <- nimble::nimbleFunction(
    run = function( ### argument type declarations
                   x = integer(),
                   n_samples = double(0),
                   a = double(1),
                   sex = double(1),
                   age2date = double(1),
                   n_cases = double(1),
                   beta_male = double(0),
                   beta0_sus = double(0),
                   beta0_inf = double(0),
                   age_effect_surv = double(1),
                   period_effect_surv = double(1),
                   f_age_foi = double(1),
                   m_age_foi = double(1),
                   age_lookup_f = double(1),
                   age_lookup_m = double(1),
                   period_lookup_foi = double(1),
                   f_period_foi = double(1),
                   m_period_foi = double(1),
                   space = double(1),
                   sect = double(1),
                   log = double()) {

    sumllik <- 0 # intialize log-likelihood
    for (i in 1:n_samples) {
        lam_inf <- nimNumeric(a[i])
        lam_foi <- nimNumeric(a[i])
        lam_sus <- nimNumeric(a[i])
        lik_temp <- nimNumeric(a[i])

        if (sex[i] == 0) { # females
            #############################################
            # preliminary hazards for the likelihood
            #############################################

            # survival hazard for susceptible deer
            lam_sus[1:a[i]] <- exp(beta0_sus +
                age_effect_surv[1:a[i]] +
                period_effect_surv[(1 + age2date[i]):(a[i] + age2date[i])])

            # survival hazard while infected
            lam_inf[1:a[i]] <- exp(beta0_inf +
                age_effect_surv[1:a[i]] +
                period_effect_surv[(1 + age2date[i]):(a[i] + age2date[i])])

            # force of infection infection hazard
            lam_foi[1:a[i]] <- exp(space[sect[i]] +
                f_age_foi[age_lookup_f[1:a[i]]] +
                f_period_foi[period_lookup_foi[(1 + age2date[i]):
                                           (a[i] + age2date[i])]])
            } else { # males
            # survival hazard for susceptible deer
            lam_sus[1:a[i]] <- exp(beta0_sus +
                age_effect_surv[1:a[i]] +
                period_effect_surv[(1 + age2date[i]):(a[i] + age2date[i])] +
                beta_male)

            # survival hazard while infected
            lam_inf[1:a[i]] <- exp(beta0_inf +
                age_effect_surv[1:a[i]] +
                period_effect_surv[(1 + age2date[i]):(a[i] + age2date[i])] +
                beta_male)

            # force of infection infection hazard
            lam_foi[1:a[i]] <- exp(space[sect[i]] + 
                m_age_foi[age_lookup_m[1:a[i]]] +
                m_period_foi[period_lookup_foi[(1 + age2date[i]):
                                          (a[i] + age2date[i])]])

            }
            #######################################
            ### calculating the joint likelihood
            #######################################

            lik_temp[1] <- lam_foi[1] *
                exp(-sum(lam_inf[1:(a[i] - 1)])) *
                lam_inf[a[i]]

            for (k in (2):(a[i] - 1)) {
                lik_temp[k] <- lam_foi[k] *
                    exp(-sum(lam_foi[1:(k - 1)])) *
                    exp(-sum(lam_sus[1:(k - 1)])) *
                    exp(-sum(lam_inf[k:(a[i] - 1)])) *
                    lam_inf[a[i]]
            }
            lik_temp[(a[i])] <- lam_foi[a[i]] *
                exp(-sum(lam_foi[1:(a[i] - 1)])) *
                exp(-sum(lam_sus[1:(a[i] - 1)])) *
                lam_inf[a[i]]

            sumllik <- sumllik + n_cases[i] * log(lam_sus[a[i]] *
                                        exp(-sum(lam_sus[1:(a[i] - 1)])) *
                                        exp(-sum(lam_foi[1:(a[i] - 1)])) +
                                        sum(lik_temp[1:a[i]]))
        }
        returnType(double(0))
        if (log) {
            return(sumllik)
        } else {
            return(exp(sumllik))
        } ## return log-likelihood
    }
)

nimble::registerDistributions(list(
    dAAH = list(
        BUGSdist = "dAAH(n_samples,a,sex,age2date,n_cases,beta_male,beta0_sus,beta0_inf,age_effect_surv,period_effect_surv,f_age_foi,m_age_foi,age_lookup_f,age_lookup_m,f_period_foi,m_period_foi,period_lookup_foi,space,sect)",
        types = c(
            "n_samples = double(0)",
            "a = double(1)",
            "sex = double(1)",
            "age2date = double(1)",
            "n_cases = double(1)",
            "beta_male = double(0)",
            "beta0_sus = double(0)",
            "beta0_inf = double(0)",
            "age_effect_surv = double(1)",
            "period_effect_surv = double(1)",
            "f_age_foi = double(1)",
            "m_age_foi = double(1)",
            "age_lookup_f = double(1)",
            "age_lookup_m = double(1)",
            "period_lookup_foi = double(1)",
            "f_period_foi = double(1)",
            "m_period_foi = double(1)",
            "space = double(1)",
            "sect = double(1)",
            "log = double()"
        ),
        discrete = TRUE
    )
))

### for a user-defined distribution
assign("dAAH", dAAH, envir = .GlobalEnv)

# starttime <- Sys.time()
# test <- dAAH(x = 1,
#     n_samples = nrow(d_fit_notest),
#     a = d_fit_notest$ageweeks,
#     sex = d_fit_notest$sexnum,
#     age2date = d_fit_notest$age2date_weeks,
#     n_cases = d_fit_notest$n,
#     beta_male = beta_male,
#     beta0_sus = beta0_survival_sus,
#     beta0_inf = beta0_survival_inf,
#     age_effect_surv = age_effect_survival_test,
#     period_effect_surv = period_effect_survival_test,
#     f_age_foi = f_age_foi,
#     m_age_foi = m_age_foi,
#     age_lookup_f = age_lookup_f,
#     age_lookup_m = age_lookup_m,
#     period_lookup_foi = period_lookup_foi,
#     f_period_foi = f_period_foi,
#     m_period_foi = m_period_foi,
#     space = c(0,-.55),
#     sect =  d_fit_notest$study_area,
#     log = TRUE
#     )
# (endtime <- Sys.time()-starttime)
# test