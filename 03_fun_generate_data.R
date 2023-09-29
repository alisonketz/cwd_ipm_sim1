ageperiod_surv_foi_sim_data <- function(
    beta0_survival_sus,
    beta0_survival_inf,
    foi_age_effect,
    age_effect,
    period_effect,
    nT_age,
    nT_period,
    prop_recap = 0.03
    ) {
  
  set.seed(12345)
  
  # sample size of individuals from each year
  n_years <- 5
  n_per_yr <- 250
  n_yr <- rep(n_per_yr, n_years)
  n_ind <- sum(n_yr)

  ### Intercepts for mortality
  beta0_sus <- beta0_survival_sus
  beta0_inf <- beta0_survival_inf

  ### maximum age interval
  nT_age <- nT_age

  ### maximum period interval
  nT_period <- nT_period

  ########################################
  # Mortality hazards
  ########################################

  age_effect_surv <- age_effect
  period_effect_surv <- period_effect

  ########################################
  # Age effects for the mortality hazard
  ########################################
  # age <- seq(1, nT_age - 1, by = 1)
  # age_effect <- age_effect
  # age_effect_surv <- age_effect
  
  ########################################
  # Age effects for the foi log hazard
  ########################################

  # logistic
  # age_effect_foi <- 1/(1+(1/.01-1)*exp(-0.7/52*age))
  # age_effect_foi <- age_effect_foi - mean(age_effect_foi)

  # png("figures/age_effects_foi_generating.png")
  # plot(age_effect_foi)
  # dev.off()


  #####################################################
  ### Age at entry with staggered entry,
  ### drawn from a stable age distribution
  ### based on the age effect hazard (but considering only susceptibles)
  ### for n_years 'years', and for the age that
  ### individuals that go from birth to censoring
  #####################################################

  hazard_se <- -exp((beta0_sus + age_effect_surv)) # not quite right because it only accounts for sus survival
  stat_se <- rep(NA, nT_age - 1)
  stat_se[1] <- (1 - exp(hazard_se[1]))
  for (j in 2:(nT_age - 1)) {
    stat_se[j] <- (1 - exp(hazard_se[j])) * exp(sum(hazard_se[1:(j - 1)]))
  }
  left_age <- nimble::rcat(n_ind, stat_se)
  maxages <- nT_age - left_age

  ###########################################
  # Staggered entry period effects
  ###########################################

  weeks_entry <- 12
  left_yr <- matrix(NA,n_years,n_per_yr)
  for(i in 1:n_years){
	left_yr[i,] <- nimble::rcat(n_per_yr,
          prob = rep(1 / weeks_entry, weeks_entry)) + (i - 1) * 52
  }
  left_period <- c(t(left_yr))

  ### dealing with maximum possible period individuals can
  ### be alive and in the study
  maxtimes <- nT_period - left_period
  maxtimes <- ifelse(maxtimes <= maxages, maxtimes, maxages)
 
  ########################################
  ### Period effects for the survival log hazard 
  ########################################

  # period <- seq(1, nT_period - 1, by = 1)
  # period_effect_surv <- 1 * sin(2/52 * pi * (period) + 1)
  # period_effect_surv <- period_effect_surv - mean(period_effect_surv)

  # png("figures/period_effects_survival_generating.png")
  # plot(period_effect_surv)
  # dev.off()

  ########################################
  ### Force of infection, log hazards 
  ########################################


  ########################################
  ### Period effects for the foi log hazard 
  ########################################

  # logistic
  # period_effect_foi <- 1/(1+(1/.001-1)*exp(-0.5/52*period))
  # period_effect_foi <- period_effect_foi - mean(period_effect_foi)

  # png("figures/period_effects_foi_generating.png")
  # plot(period_effect_foi)
  # dev.off()
  ########################################
  ### Log hazards 
  ########################################
  hazard_sus <- hazard_inf <- hazard_foi <- matrix(0, n_ind, nT_age)
  for (i in 1:n_ind) {
    hazard_sus[i,left_age[i]:(left_age[i] + maxtimes[i] - 1)] <-  exp(beta0_sus+
                      age_effect_surv[left_age[i]:(left_age[i] + maxtimes[i] - 1)] +
                      period_effect_surv[left_period[i]:(left_period[i] + maxtimes[i] - 1)])

    hazard_inf[i,left_age[i]:(left_age[i] + maxtimes[i] - 1)] <- exp(beta0_inf+
                      age_effect_surv[left_age[i]:(left_age[i] + maxtimes[i] - 1)] +
                      period_effect_surv[left_period[i]:(left_period[i] + maxtimes[i] - 1)])

    hazard_foi[i,1:(left_age[i] + maxtimes[i] - 1)] <- exp(
                      foi_age_effect[1:(left_age[i] + maxtimes[i] - 1)])
  }


  ########################################
  ### Probability of infection
  ########################################

 # incidence <- matrix(0, n_ind, nT_age)
 # for (i in 1:n_ind) {
 #   for (j in left_age[i]:(left_age[i] + maxages[i] - 1)) {
 #     if (j==left_age[i]) {
 #       incidence[i, j] <- (1 - exp(-exp(hazard_foi[i, j])))
 #     } else {
 #       incidence[i, j] <-
 #         (1 - exp(-exp(hazard_foi[i, j])))*exp(-sum(exp(hazard_foi[i, left_age[i]:(j - 1)])))
 #     }
 #   }
 #   incidence[i, left_age[i] + maxages[i]] <-
 #         exp(-sum(exp(hazard_foi[i, left_age[i]:(left_age[i] + maxages[i] - 1)])))
 # }

  incidence <- matrix(0, n_ind, nT_age)
  for (i in 1:n_ind) {
    for (j in 1:(left_age[i] + maxtimes[i] - 1)) {
        incidence[i, j] <- (1 - exp(-hazard_foi[i, j]))
      } 
  }


  ########################################
  ### Calculating infection status
  ########################################
  inf_age <- rep(NA, n_ind)
  inf_status <- matrix(0, n_ind, nT_age)
  for(i in 1:n_ind) {
    inf_status[i,] <- rbinom(nT_age,1,incidence[i,1:nT_age])
    if(max(inf_status[i,]) == 0) next
      inf_age[i] <- min(which(inf_status[i,] == 1))
      inf_status[i,inf_age[i]:nT_age] <- 1
  }

  ########################################
  ### Probability of mortality
  ########################################

  survival_prob <- matrix(0, n_ind, nT_age)
  for (i in 1:n_ind) {
    for (j in left_age[i]:(left_age[i] + maxtimes[i] - 1)) {
      if (j==left_age[i]) {
        survival_prob[i, j] <- 
          1 - ((1 - inf_status[i, j]) * exp(-hazard_sus[i, j]) + 
				  inf_status[i, j] * exp(-hazard_inf[i, j]))
      } else {
	      survival_prob[i, j] <- 
          (1 - ((1 - inf_status[i, j])*exp(-hazard_sus[i, j]) + 
				  inf_status[i, j] * exp(-hazard_inf[i, j]))) *
			    exp(-(sum((1 - inf_status[i, left_age[i]:(j - 1)]) * hazard_sus[i, left_age[i]:(j - 1)]) +
				      sum(inf_status[i, left_age[i]:(j - 1)] * hazard_inf[i, left_age[i]:(j - 1)])))
	    }
    }
    survival_prob[i, left_age[i] + maxtimes[i]] <-
		 (1 - inf_status[i, left_age[i] + maxtimes[i]])*exp(-sum(hazard_sus[i, left_age[i]:(left_age[i] + maxtimes[i] - 1)])) + 
			  inf_status[i, left_age[i] + maxtimes[i]]*exp(-sum(hazard_inf[i, left_age[i]:(left_age[i] + maxtimes[i] - 1)]))
  }


  ########################################
  ### Calculating right censoring
  ########################################
  fail_int <- rep(0, n_ind)
  right_age <- rep(0, n_ind)
  rt_censor <- rep(0, n_ind)
  test1 <- left_age
  pos1 <- rep(0, n_ind)
  test2 <- rep(NA, ) # after death or censoring (only if 1st test was neg)
  p_test2 <- 0.65
  pos2 <- rep(0, n_ind)
 
  for(i in 1:n_ind) {
    fail_int[i] <- which(rmultinom(1, 1, survival_prob[i,1:nT_age]) == 1)
    right_age[i] <- ifelse(fail_int[i] >= left_age[i] + maxtimes[i],
                           left_age[i] + maxtimes[i],
                           fail_int[i] + 1)
    rt_censor[i] <- ifelse(fail_int[i] < (left_age[i] + maxtimes[i]), 0, 1)
    if(!is.na(inf_age[i])){
      if(inf_age[i] >= right_age[i]){
        inf_age[i] <- NA
      } else {
        if(inf_age[i] < left_age[i]) {
          pos1[i] <- 1
          }
      }
    }
    # test2yes <- ifelse(pos1[i] == 1, 0, rbinom(1, 1, p_test2))
    # if(test2yes == 1) {
    #   test2[i] <- right_age[i] # for simplicity, assume test2 date is immediately after mortality or censoring
    #   pos2[i] <- ifelse(!is.na(inf_age[i]), 1, 0)
    # }
    if(pos1[i]==1){
      pos2[i] <- 1
    } else {
      if(!is.na(inf_age[i])) {
        if(right_age[i] > inf_age[i]) {
          pos2[i] <- 1
        }
      }
    }

  }

  right_period <- right_age - left_age + left_period


  ########################################
  ### setup data into overall data frame
  ########################################

  df_fit <- data.frame(id = 1:n_ind,
                      left_age = left_age,
                      right_age  = right_age,
                      left_period = left_period,
                      right_period = right_period,
                      rt_censor = rt_censor,
                      cwd_cap = pos1,
                      cwd_mort = pos2,
                      inf_age = inf_age
                      )

  df_fit$inf_period = df_fit$right_age - df_fit$inf_age + df_fit$left_period

  ### setting up in terms of e/r/s
  df_fit$e_age <- df_fit$left_age
  df_fit$r_age <- df_fit$right_age
  df_fit$s_age <- df_fit$right_age
  df_fit$r_age[df_fit$rt_censor == 0] <- df_fit$r_age[df_fit$rt_censor == 0] - 1
  df_fit$s_age[df_fit$rt_censor == 1] <- NA


  df_fit$e_period <- df_fit$left_period
  df_fit$r_period <- df_fit$right_period
  df_fit$s_period <- df_fit$right_period
  df_fit$r_period[df_fit$rt_censor == 0] <- df_fit$r_period[df_fit$rt_censor == 0] - 1
  df_fit$s_period[df_fit$rt_censor == 1] <- NA

  ########################################
  ### setup all data types
  ########################################



  ################################################
  ###
  ### remove the over-abundance infected
  ### at capture, build in 4 options
  ###   no removal, 
  ###   random sample,
  ###   weighted sample by left_age/sum(left_age)
  ###   weighted sample by (1-1/left_age)/sum(1-1/left_age))
  ###
  #################################################


  ########################################
  ### Sampling structuring recap values
  ########################################

  ### this is just a sketch and not functional yet

  # recap_sample_indx <- which((df_typr_period - e_period) > 52)
  # recap_sample <- sample(recap_sample_indx,prop_recap)

  # #catch if there's no recap_sample_indx
  # if(length(recap_sample == 0)) recap_sample <- recap_sample_indx[1]

  # for(i in recap_sample_indx){
  #   df$recap[i] <- sample(e_period[i]:(r_period[i]-2), size = 1, replace = FALSE)
  # }


  ########################################
  ### Return values
  ########################################

  return(list(n_ind = n_ind,
              nT_age = nT_age,
              nT_period = nT_period,
              beta0_sus = beta0_sus,
              beta0_inf = beta0_inf,
              period_effect_surv = period_effect_surv,
              age_effect_surv = age_effect_surv,
              foi_age_effect = foi_age_effect,
              hazard_sus = hazard_sus,
              hazard_inf = hazard_inf,
              hazard_foi = hazard_foi,
              survival_prob = survival_prob,
              incidence = incidence,
              right_age = right_age,
              left_age = left_age,
              right_period = right_period,
              left_period = left_period,
              rt_censor = rt_censor,
              inf_age = inf_age,
              inf_status = inf_status,
              test1 = test1,
              test2 = test2,
              pos1 = pos1,
              pos2 = pos2
              )
         )
}

