##################################################
###
###  Function state transition probability
###  Age-Period Model
###
###################################################

state_transition <- nimbleFunction(
  run = function(records = double(0),
                 left = double(1),
                 right = double(1),
                 age_effect = double(1),
                 period_effect = double(1),
                 age2date = double(1),
                 nT_age = double(0),
                 beta0 = double(0)
  ){

    SLR <- nimNumeric(records)
    UCH <- nimMatrix(value = 0, nrow = records, ncol = nT_age)
    for (j in 1:records) {
      for (k in left[j]:(right[j] - 1)) {
        UCH[j, k] <- exp(beta0 +
                         age_effect[k] +
                         period_effect[k - age2date[j]])
      }
      SLR[j] <- exp(-sum(UCH[j, left[j]:(right[j] - 1)]))
    }
    returnType(double(1))
    return(SLR[1:records])
  })

cstate_transition <- compileNimble(state_transition)

###########################
###
### Model Specification
###
###########################

modelcode <- nimbleCode({
  
  ### Prior for intercept
  ### using parameter expansion for convergence

  beta0_temp ~ dnorm(0,0.01)
  mix ~ dunif(-1, 1)
  beta0 <- beta0_temp * mix

  ###
  ### Age Effects & Period Effects
  ###

  for (k in 1:nknots_age) {
    ln_b_age[k] ~ dnorm(0, tau_age)
    b_age[k] <- exp(ln_b_age[k])
  }
  tau_age ~ dgamma(1,1)

  for (k in 1:nknots_period) {
    b_period[k] ~ dnorm(0, tau_period)
  }
  tau_period ~ dgamma(.1,.1)
  
  for (t in 1:nT_age) {
    age_effect_temp[t] <- inprod(b_age[1:nknots_age], Z_age[t, 1:nknots_age])
    age_effect[t] <- age_effect_temp[t] - mu_age
  }
  mu_age <- mean(age_effect_temp[1:nT_age])
  for (t in 1:nT_period) {
    period_effect[t] <- inprod(b_period[1:nknots_period],
                               Z_period[t, 1:nknots_period])
  }

  ###
  ### Computing state transisiton probability
  ###

  SLR[1:records] <- state_transition(records = records,
                                   left = left_age[1:records],
                                   right = right_age[1:records],
                                   nT_age = nT_age,
                                   age_effect = age_effect[1:nT_age],
                                   period_effect = period_effect[1:nT_period],
                                   age2date = age2date[1:records],
                                   beta0 = beta0)

  for (j in 1:records) {
    censor[j] ~ dbern(SLR[j])
  }

  ##########################
  ### Derived parameters
  ##########################

  for (t in 1:nT_age) {
    llambda_age[t] <- beta0 + age_effect[t]
    UCH0_age[t] <- exp(llambda_age[t])
    S0_age[t] <- exp(-sum(UCH0_age[1:t]))
  }
  for (t in 1:nT_period) {
    llambda_period[t] <- beta0 + period_effect[t]
    UCH0_period[t] <- exp(llambda_period[t])
    S0_period[t] <- exp(-sum(UCH0_period[1:t]))
  }

})#end model statement

#Data
nimData <- list(censor = df_fit[, 3],
                Z_period = Z_period,
                Z_age = Z_age,
                left_age = df_fit[, 1],
                right_age = df_fit[, 2],
                age2date = age2date
                )
#Constants
nimConsts <- list(records = n_fit,
                 nT_age = nT_age,
                 nT_period = nT_period,
                 nknots_age = nknots_age,
                 nknots_period = nknots_period)

#Specify initial values

initsFun <- function()list(
                          tau_period = runif(1),
                          beta0_temp = rnorm(1, beta0_true, .001),
                          mix = 1,
                          b_period = rnorm(nknots_period)*10^-4,
                          ln_b_age = runif(nknots_age, -10, -5),
                          tau_age = runif(1, .1, 1)
                          )
nimInits <- initsFun()


Rmodel <- nimbleModel(code = modelcode,
                      constants = nimConsts,
                      data = nimData,
                      inits = initsFun()
                      )

#identify params to monitor
parameters <- c(
              "beta0",
              "b_age",
              "S0_age",
              "age_effect",
              "period_effect",
              "S0_period",
              "b_age",
              "b_period",
              "tau_period",
              "tau_age",
              "mu_age"
)


starttime <- Sys.time()
confMCMC <- configureMCMC(Rmodel,
                          monitors = parameters,
                          thin = n_thin,
                          useConjugacy = FALSE,
                          enableWAIC = TRUE)
confMCMC$addSampler(target = "ln_b_age", type = "RW_block")
nimMCMC <- buildMCMC(confMCMC)
Cnim <- compileNimble(Rmodel)
CnimMCMC <- compileNimble(nimMCMC, project = Rmodel)
mcmcout <- runMCMC(CnimMCMC,
                   niter = reps,
                   nburnin = bin,
                   nchains = n_chains,
                   inits = initsFun,
                   samplesAsCodaMCMC = TRUE,
                   summary = TRUE,
                   WAIC = TRUE)

runtime <- difftime(Sys.time(), starttime, units = "min")

save(runtime, file = "runtime.Rdata")
save(mcmcout, file = "mcmcout.Rdata")