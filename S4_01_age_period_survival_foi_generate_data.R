ageperiod_surv_foi_sim_data <- function() {
  
  set.seed(12345)
  
  # sample size of individuals from each year
  n_years <- 5
  n_per_yr <- 500
  n_yr <- rep(n_per_yr, n_years)
  n <- sum(n_yr)
  
  # intercept
  beta0_sus <- -5
  beta0_inf <- -4
  beta0_foi <- -6

  ### maximum age interval
  nT_age <- 963

  ### maximum period interval
  nT_period <- n_years * 52

  ########################################
  # Age effects for the survival log hazard
  ########################################
  age <- seq(1, nT_age - 1, by = 1)
  
  #Inducing wiggles during the early age intervals
  # anthro <- .5 * cos(1 / 26 * pi * age)
  # anthro[1:19] <- 0
  # anthro[20:38] <- anthro[20:38] * seq(0, 1, length = length(20:38))
  # anthro[50:170] <- anthro[50:170] * seq(1, 0, length = length(50:170))
  # anthro[171:(nT_age - 1)] <- 0 

  #Age effects on log hazard
  #baseline is a quadratic function
  lam <- .5
  age_effect_surv <- 10 * lam * age^(-(1 - lam))
  ###quadratic
  # age_effect_surv <- 5e-05 * age^2 - 0.02 * age + 1
  # age_effect_surv <- age_effect_surv #+ anthro
  age_effect_surv <- age_effect_surv - mean(age_effect_surv)

  png("figures/age_effects_survival_generating.png")
  plot(age_effect_surv)
  dev.off()

  ########################################
  # Age effects for the foi log hazard
  ########################################

  # logistic
  age_effect_foi <- 1/(1+(1/.01-1)*exp(-0.7/52*age))

  # age_effect_foi <- age_effect_foi - mean(age_effect_foi)

  png("figures/age_effects_foi_generating.png")
  plot(age_effect_foi)
  dev.off()


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
  left_age <- nimble::rcat(n, stat_se)
  maxages <- nT_age - left_age

  ###########################################
  # Staggered entry period effects
  ###########################################

  weeks_entry <- 12
  left_yr <- matrix(NA,n_years,n_per_yr)
  for(i in 1:n_years){
	left_yr[i,] <- nimble::rcat(n_per_yr, prob = rep(1 / weeks_entry, weeks_entry))+(i-1)*52
  }
  left_period <- c(t(left_yr))

  # no staggered entry
  maxtimes <- nT_period - left_period
  maxtimes <- ifelse(maxtimes <= maxages,maxtimes,maxages)

 
  ########################################
  ### Period effects for the survival log hazard 
  ########################################

  period <- seq(1, nT_period - 1, by = 1)
  period_effect_surv <- 1 * sin(2/52 * pi * (period) + 1)
  period_effect_surv <- period_effect_surv - mean(period_effect_surv)

  png("figures/period_effects_survival_generating.png")
  plot(period_effect_surv)
  dev.off()


  ########################################
  ### Period effects for the foi log hazard 
  ########################################

  # logistic
  period_effect_foi <- 1/(1+(1/.001-1)*exp(-0.5/52*period))
  period_effect_foi <- period_effect_foi - mean(period_effect_foi)

  png("figures/period_effects_foi_generating.png")
  plot(period_effect_foi)
  dev.off()
  ########################################
  ### Log hazards
  ########################################
  hazard_sus <- hazard_inf <- hazard_foi <- matrix(0, n, nT_age)
  for (i in 1:n) {
    hazard_sus[i,left_age[i]:(left_age[i] + maxtimes[i] - 1)] <-  exp(beta0_sus+
                      age_effect_surv[left_age[i]:(left_age[i] + maxtimes[i] - 1)] +
                      period_effect_surv[left_period[i]:(left_period[i] + maxtimes[i] - 1)])
    hazard_inf[i,left_age[i]:(left_age[i] + maxtimes[i] - 1)] <- exp(beta0_inf+
                      age_effect_surv[left_age[i]:(left_age[i] + maxtimes[i] - 1)] +
                      period_effect_surv[left_period[i]:(left_period[i] + maxtimes[i] - 1)])
					  period_ind <- ifelse(left_period[i]-left_age[i] >=0,
										   left_period[i]-left_age[i]+1,
										   c(rep(1,left_age[i]-left_period[i]),1:(left_period[i] + maxtimes[i] - 1)))
    hazard_foi[i,1:(left_age[i] + maxtimes[i] - 1)] <- exp(beta0_foi +
                      age_effect_foi[1:(left_age[i] + maxtimes[i] - 1)] +
                      period_effect_foi[period_ind])
  }


  ########################################
  ### Probability of infection
  ########################################

 # test_stat_foi <- matrix(0, n, nT_age)
 # for (i in 1:n) {
 #   for (j in left_age[i]:(left_age[i] + maxages[i] - 1)) {
 #     if (j==left_age[i]) {
 #       test_stat_foi[i, j] <- (1 - exp(-exp(hazard_foi[i, j])))
 #     } else {
 #       test_stat_foi[i, j] <-
 #         (1 - exp(-exp(hazard_foi[i, j])))*exp(-sum(exp(hazard_foi[i, left_age[i]:(j - 1)])))
 #     }
 #   }
 #   test_stat_foi[i, left_age[i] + maxages[i]] <-
 #         exp(-sum(exp(hazard_foi[i, left_age[i]:(left_age[i] + maxages[i] - 1)])))
 # }

  test_stat_foi <- matrix(0, n, nT_age)
  for (i in 1:n) {
    for (j in 1:(left_age[i] + maxtimes[i] - 1)) {
        test_stat_foi[i, j] <- (1 - exp(-hazard_foi[i, j]))
      } 
  }


  ########################################
  ### Calculating infection status
  ########################################
  inf_age <- rep(NA, n)
  inf_status <- matrix(0, n, nT_age)
  for(i in 1:n) {
    inf_status[i,] <- rbinom(nT_age,1,test_stat_foi[i,1:nT_age])
	if(max(inf_status[i,])==0) next
	inf_age[i] <- min(which(inf_status[i,] == 1))
	inf_status[i,inf_age[i]:nT_age] <- 1
  }


  ########################################
  ### Probability of mortality
  ########################################

  test_stat_surv <- matrix(0, n, nT_age)
  for (i in 1:n) {
    for (j in left_age[i]:(left_age[i] + maxtimes[i] - 1)) {
      if (j==left_age[i]) {
        test_stat_surv[i, j] <- 1 - (
				(1 - inf_status[i,j])*exp(-hazard_sus[i, j]) + 
				inf_status[i,j]*exp(-hazard_inf[i, j]))
      } else {
	    test_stat_surv[i, j] <- (1 - (
				(1 - inf_status[i,j])*exp(-hazard_sus[i, j]) + 
				inf_status[i,j]*exp(-hazard_inf[i, j]))
				)*
			    exp(-(sum((1 - inf_status[i,left_age[i]:(j - 1)]) * hazard_sus[i, left_age[i]:(j - 1)]) +
				      sum(inf_status[i,left_age[i]:(j - 1)] * hazard_inf[i, left_age[i]:(j - 1)])))
	  }
    }
    test_stat_surv[i, left_age[i] + maxtimes[i]] <-
		 (1 - inf_status[i,left_age[i] + maxtimes[i]])*exp(-sum(hazard_sus[i, left_age[i]:(left_age[i] + maxtimes[i] - 1)])) + 
			  inf_status[i,left_age[i] + maxtimes[i]]*exp(-sum(hazard_inf[i, left_age[i]:(left_age[i] + maxtimes[i] - 1)]))
  }


  ########################################
  ### Calculating right censoring
  ########################################
  fail_int <- rep(0, n)
  right_age <- rep(0, n)
  rt_censor <- rep(0, n)
  test1 <- left_age
  pos1 <- rep(0, n)
  test2 <- rep(NA, ) # after death or censoring (only if 1st test was neg)
  p_test2 <- 0.65
  pos2 <- rep(NA, n)
 
  for(i in 1:n) {
    fail_int[i] <- which(rmultinom(1, 1, test_stat_surv[i,1:nT_age]) == 1)
    right_age[i] <- ifelse(fail_int[i] >= left_age[i] + maxtimes[i],
                           left_age[i] + maxtimes[i], fail_int[i] + 1)
    rt_censor[i] <- ifelse(fail_int[i] < (left_age[i] + maxtimes[i]), 0, 1)
	if(!is.na(inf_age[i])){
		if(inf_age[i] >= right_age[i]){
			inf_age[i] <- NA
		}else{
			if(inf_age[i] < left_age[i]){pos1[i] <- 1}
		}
	}
	test2yes <- ifelse(pos1[i] == 1, 0, rbinom(1,1,p_test2))
	if(test2yes == 1){
		test2[i] <- right_age[i] # for simplicity, assume test2 date is immediately after mortality or censoring
		pos2[i] <- ifelse(!is.na(inf_age[i]), 1, 0)
	}
  }

  right_period <- right_age - left_age + left_period

  ########################################
  ### Return values
  ########################################

  return(list(n = n,
              nT_age = nT_age,
              nT_period = nT_period,
              beta0_sus = beta0_sus,
              beta0_inf = beta0_inf,
              beta0_foi = beta0_foi,
              period_effect_surv = period_effect_surv,
			  age_effect_surv = age_effect_surv,
              period_effect_foi = period_effect_foi,
			  age_effect_foi = age_effect_foi,
              hazard_sus = hazard_sus,
              hazard_inf = hazard_inf,
              hazard_foi = hazard_foi,
              test_stat_surv = test_stat_surv,
              test_stat_foi = test_stat_foi,
              right_age = right_age,
              left_age = left_age,
              right_period = right_period,
              left_period = left_period,
              rt_censor = rt_censor,
			  test1 = test1,
			  test2 = test2,
			  pos1 = pos1,
			  pos2 = pos2,
              prop_right_cens = sum(right_age == nT_age) / length(right_age)
              )
         )
}

table(dat$pos1)
table(dat$pos2)