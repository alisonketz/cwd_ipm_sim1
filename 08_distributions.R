#######################################################################
###
###   Likelihoods for each kind of data listed below
###
#######################################################################

# d_fit_icap_mort
# d_fit_icap_cens
# d_fit_idead
# d_fit_sus_mort_posttest
# d_fit_sus_cens_postno

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

dIcapMort <- nimble::nimbleFunction(
    run = function( ### argument type declarations
                   x = integer(0),
                   n_samples = integer(0), # number of samples in dataset
                   e = double(1), # e, age of entry
                   r = double(1), # r, age of last known alive
                   s = double(1), # s, age of known mortality
                   age2date = double(1),
                   beta0_inf = double(0),
                   age_effect_surv = double(1),
                   period_effect_surv = double(1),
                   log = double()) {

        # start the loop through individuals
        sumllik <- 0
        for (i in 1:n_samples) {
            lam_inf <- 0
            #############################################
            # preliminary hazards for the likelihood
            #############################################

            if(r[i]>e[i]) {
                    for (j in e[i]:(r[i] - 1)) {
                    lam_inf <- lam_inf +
                               exp(beta0_inf +
                                   age_effect_surv[j] +
                                   period_effect_surv[age2date[i] + j]
                                   )
                    }
                }
            lam_inf_s <- sum(exp(beta0_inf +
                age_effect_surv[r[i]:(s[i] - 1)] +
                period_effect_surv[age2date[i] + r[i]:(s[i] - 1)]))
                
            #######################################
            ### calculating the joint likelihood
            #######################################

            sumllik <- sumllik - lam_inf +
                       log(1 - exp(-lam_inf_s))
        }#end n_samples

        returnType(double(0))
        if (log) {
            return(sumllik)
        } else {
            return(exp(sumllik))
        } ## return log-likelihood
})


nimble::registerDistributions(list(
    dIcapMort = list(
        BUGSdist = "dIcapMort(n_samples,e,r,s,age2date,beta0_inf,age_effect_surv,period_effect_surv)",
        types = c(
            "value = integer(0)",
            "n_samples = integer(0)",
            "e = double(1)",
            "r = double(1)",
            "s = double(1)",
            "age2date = double(1)",
            "beta0_inf = double(0)",
            "age_effect_surv = double(1)",
            "period_effect_surv = double(1)",
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
#         e = d_fit_icap_mort$e_age,
#         r = d_fit_icap_mort$r_age,
#         s = d_fit_icap_mort$s_age,
#         age2date = d_fit_icap_mort$age2date,
#         beta0_inf = beta0_survival_inf,
#         age_effect_surv = age_effect_survival_test,
#         period_effect_surv = period_effect_surv,
#         log = TRUE
#         )
# (end <- Sys.time() - starttime)
# test


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

dIcapCens <- nimble::nimbleFunction(
    run = function( ### argument type declarations
                   x = integer(0),
                   n_samples = integer(0), # number of samples in dataset
                   e = double(1), # e, age of entry
                   r = double(1), # r, age of last known alive
                   age2date = double(1),
                   beta0_inf = double(0),
                   age_effect_surv = double(1),
                   period_effect_surv = double(1),
                   log = double()) {
        sumllik <- 0 # intialize log-likelihood
        for (i in 1:n_samples) {
            lam_inf <- 0

            #############################################
            # preliminary hazards for the likelihood
            #############################################

            # recalculate lam_inf to sum over e to r-1
            for (j in e[i]:(r[i] - 1)) {
                lam_inf <- lam_inf +
                            exp(beta0_inf +
                                age_effect_surv[j] +
                                period_effect_surv[age2date[i] + j])
            }

            #######################################
            ### calculating the joint likelihood
            #######################################
            sumllik <- sumllik - lam_inf
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
        BUGSdist = "dIcapCens(n_samples,e,r,age2date,beta0_inf,age_effect_surv,period_effect_surv)",
        types = c(
            "value=integer(0)",
            "n_samples = integer(0)",
            "e = double(1)",
            "r = double(1)",
            "age2date = double(1)",
            "beta0_inf = double(0)",
            "age_effect_surv = double(1)",
            "period_effect_surv = double(1)",
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
#         e = d_fit_icap_cens$e_age,
#         r = d_fit_icap_cens$r_age,
#         age2date = icap_cens_age2date,
#         beta0_inf = beta0_survival_inf,
#         age_effect_surv = age_effect_survival_test,
#         period_effect_surv = period_effect_surv,
#         log = TRUE
#         )
# (endtime <- Sys.time() - starttime)
# test


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

dNegCapPosMort <- nimble::nimbleFunction(
    run = function( ### argument type declarations
                   x = integer(0),
                   n_samples = integer(0), # number of samples in dataset
                   e = double(1), # e, age of entry
                   r = double(1), # r, age of last known alive
                   s = double(1), # s, age of known mortality
                   age2date = double(1),
                   beta0_sus = double(0),
                   beta0_inf = double(0),
                   age_effect_surv = double(1),
                   period_effect_surv = double(1),
                   foi_age_effect = double(1),
                   age_lookup = double(1),
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
            # survival hazard for susceptible deer
            lam_sus[e[i]:(s[i] - 1)] <- exp(beta0_sus +
                age_effect_surv[e[i]:(s[i] - 1)] +
                period_effect_surv[(e[i] + age2date[i]):
                                    (s[i] - 1 + age2date[i])])
            # survival hazard while infected
            lam_inf[e[i]:(s[i] - 1)] <- exp(beta0_inf +
                age_effect_surv[e[i]:(s[i] - 1)] +
                period_effect_surv[(e[i] + age2date[i]):
                                    (s[i] - 1 + age2date[i])])
            # force of infection infection hazard
            lam_foi[e[i]:(s[i] - 1)] <- exp(foi_age_effect[age_lookup[e[i]:(s[i] - 1)]])
        

            #######################################
            ### calculating the joint likelihood
            #######################################
            lik_temp <- lik_temp +
                        (1 - exp(-lam_foi[e[i]])) *
                        exp(-sum(lam_inf[e[i]:(r[i] - 1)]))
            if((r[i] - e[i])>2) {
                for (k in (e[i] + 1):(r[i] - 1)) {
                        lik_temp <- lik_temp +
                                    (1 - exp(-lam_foi[k])) *
                                    exp(-sum(lam_inf[(k):(r[i] - 1)])) *
                                    exp(-sum(lam_sus[e[i]:(k - 1)])) *
                                    exp(-sum(lam_foi[e[i]:(k - 1)]))
                }
            }
            lik_temp <- lik_temp +
                        (1 - exp(-lam_foi[s[i] - 1])) *
                        exp(-sum(lam_sus[e[i]:((s[i] - 1) - 1)])) *
                        exp(-sum(lam_foi[e[i]:((s[i] - 1) - 1)]))

            sumllik <- sumllik +
                       log(1 - exp(-sum(lam_inf[r[i]:(s[i] - 1)]))) +
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
        BUGSdist = "dNegCapPosMort(n_samples,e,r,s,age2date,beta0_sus,beta0_inf,age_effect_surv,period_effect_surv,foi_age_effect,age_lookup)",
        types = c(
            "value = integer(0)",
            "n_samples = integer(0)",
            "e = double(1)",
            "r = double(1)",
            "s = double(1)",
            "age2date = double(1)",
            "beta0_sus = double(0)",
            "beta0_inf = double(0)",
            "age_effect_surv = double(1)",
            "period_effect_surv = double(1)",
            "foi_age_effect = double(1)",
            "age_lookup = double(1)",
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
#         n_samples = nrow(df_fit_idead),
#         e = df_fit_idead$e_age,
#         r = df_fit_idead$r_age,
#         s = df_fit_idead$s_age,
#         age2date = df_fit_idead$age2date,
#         beta0_sus = beta0_survival_sus,
#         beta0_inf = beta0_survival_inf,
#         age_effect_surv = age_effect_true,
#         period_effect_surv = period_effect_true,
#         foi_age_effect = age_foi,
#         age_lookup = age_lookup,
#         log = TRUE
#         )
# (end<- Sys.time()-starttime)
# test


#######################################################################
###
###   User defined distribution for likelihood for
###   uninfected radio-marked deer mortalities:
###   test neg at cap and tested mort
###
###   d_fit_sus_mort_posttest
###
###   Overleaf Equation (6)
###
#######################################################################

dSusMortTest <- nimble::nimbleFunction(
    run = function( ### argument type declarations
                   x = double(0),
                   n_samples = integer(0), # number of samples in dataset
                   e = double(1), # e, age of entry
                   r = double(1), # r, age of last known alive
                   s = double(1), # s, age of known mortality
                   age2date = double(1),
                   beta0_sus = double(0),
                   age_effect_surv = double(1),
                   period_effect_surv = double(1),
                   foi_age_effect = double(1),
                   age_lookup = double(1),
                   log = double()) {
        sumllik <- 0 # intialize log-likelihood
        for (i in 1:n_samples) {
            lam_foi <- 0
            lam_sus <- 0
            lam_sus_s <- 0

            #############################################
            # preliminary hazards for the likelihood
            #############################################
            for (j in r[i]:(s[i] - 1)) {
                # sum up surv hazard from r to s-1
                lam_sus_s <- lam_sus_s +
                                exp(beta0_sus +
                                age_effect_surv[j] +
                                period_effect_surv[age2date[i] + j])
            }

            for (j in e[i]:(s[i] - 1)) {
                #  up foi hazard from 1  to j
                lam_foi <- lam_foi +
                            exp(foi_age_effect[age_lookup[j]])
                if (j < r[i]) {
                lam_sus <- lam_sus +
                            exp(beta0_sus +
                                age_effect_surv[j] +
                                period_effect_surv[age2date[i] + j])
                }
            }
            #######################################
            ### calculating the joint likelihood
            #######################################

            sumllik <- sumllik +
                       log(1 - exp(-lam_sus_s)) -
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
        BUGSdist = "dSusMortTest(n_samples,e,r,s,age2date,beta0_sus,age_effect_surv,period_effect_surv,foi_age_effect,age_lookup)",
        types = c(
            "value = double(0)",
            "n_samples = integer(0)",
            "e = double(1)",
            "r = double(1)",
            "s = double(1)",
            "age2date = double(1)",
            "beta0_sus = double(0)",
            "age_effect_surv = double(1)",
            "period_effect_surv = double(1)",
            "foi_age_effect = double(1)",
            "age_lookup = double(1)",
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
#         n_samples = nrow(d_fit_sus_mort_posttest),
#         e = d_fit_sus_mort_posttest$e_age,
#         r = d_fit_sus_mort_posttest$r_age,
#         s = d_fit_sus_mort_posttest$s_age,
#         age2date = sus_mort_posttest_age2date,
#         beta0_sus = beta0_survival_sus,
#         age_effect_surv = age_effect_survival_test,
#         period_effect_surv = period_effect_surv,
#         foi_age_effect = foi_age_effect,
#         age_lookup = age_lookup
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
###
###   Overleaf Equation (4)
###
#######################################################################

dSusCensNo <- nimble::nimbleFunction(
    run = function( ### argument type declarations
                   x = integer(0),
                   n_samples = integer(0), # number of samples in dataset
                   e = double(1), # e, age of entry
                   r = double(1), # r, age of last known alive
                   age2date = double(1),
                   beta0_sus = double(0),
                   beta0_inf = double(0),
                   age_effect_surv = double(1),
                   period_effect_surv = double(1),
                   foi_age_effect = double(1),
                   age_lookup = double(1),
                   log = double()) {
        sumllik <- 0 # intialize log-likelihood
        for (i in 1:n_samples) {
            lam_foi <- nimNumeric(r[i] - 1)
            lam_sus <- nimNumeric(r[i] - 1)
            lam_inf <- nimNumeric(r[i]- 1)
            lik_temp <- 0

            #############################################
            # preliminary hazards for the likelihood
            #############################################
            lam_foi[e[i]:(r[i] - 1)] <- exp(foi_age_effect[age_lookup[e[i]:(r[i]-1)]])

            lam_sus[e[i]:(r[i] - 1)] <- exp(beta0_sus +
                                        age_effect_surv[e[i]:(r[i] - 1)] +
                period_effect_surv[(e[i] + age2date[i]):
                                    (r[i] - 1 + age2date[i])])
            lam_inf[e[i]:(r[i] - 1)] <- exp(beta0_inf +
                    age_effect_surv[e[i]:(r[i] - 1)] +
                    period_effect_surv[(age2date[i] + e[i]):
                                        (age2date[i] + r[i] - 1)])

            lik_temp <- (1 - exp(-lam_foi[e[i]])) *
                        exp(-sum(lam_inf[e[i]:(r[i] - 1)]))

            if ((r[i] - e[i]) > 1) {
                for(k in (e[i] + 1):(r[i] - 1)) {
                    lik_temp <- lik_temp +
                                (1 - exp(-lam_foi[k])) *
                                exp(-sum(lam_inf[k:(r[i] - 1)])) *
                                exp(-sum(lam_foi[e[i]:(k - 1)])) *
                                exp(-sum(lam_sus[e[i]:(k - 1)])) 
                }
            }
            #######################################
            ### calculating the joint likelihood
            #######################################
            sumllik <- sumllik +
                       log(
                        exp(-sum(lam_foi[e[i]:(r[i] - 1)])) *
                        exp(-sum(lam_sus[e[i]:(r[i] - 1)])) +
                        lik_temp
                       )
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
        BUGSdist = "dSusCensNo(n_samples,e,r,age2date,beta0_sus,beta0_inf,age_effect_surv,period_effect_surv,foi_age_effect,age_lookup)",
        types = c(
            "value=integer(0)",
            "n_samples = integer(0)",
            "e = double(1)",
            "r = double(1)",
            "age2date = double(1)",
            "beta0_sus = double(0)",
            "beta0_inf = double(0)",
            "age_effect_surv = double(1)",
            "period_effect_surv = double(1)",
            "foi_age_effect = double(1)",
            "age_lookup = double(1)",
            "log = double()"
        ),
        discrete = TRUE
    )
))

### for a user-defined distribution
assign("dSusCensNo", dSusCensNo, envir = .GlobalEnv)

# starttime <- Sys.time()
# test <- dSusCensNo(
#         x = 1,
        # n_samples = nrow(d_fit_sus_cens_postno),
        # e = d_fit_sus_cens_postno$e_age,    
        # r = d_fit_sus_cens_postno$r_age,   
        # age2date = sus_cens_postno_age2date, 
#         beta0_sus = beta0_survival_sus,
#         beta0_inf = beta0_survival_inf,
#         age_effect_surv = age_effect_survival_test,
#         period_effect_surv = period_effect_surv,
#         foi_age_effect = foi_age_effect,
#         age_lookup = age_lookup
#         log = TRUE
#         )
# (end<- Sys.time()-starttime)
# test
