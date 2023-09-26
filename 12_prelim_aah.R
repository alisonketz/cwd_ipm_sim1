###########################################################################
###
### Setup for Age-At-Harvest (AAH population model)
###
###########################################################################

####################################
###
### Calculating sex ratio prior to 
### study starting
###
####################################

df_age_before_female <- df_age_before %>% 
                          group_by(age, sex) %>%
                          summarise(n = sum(n)) %>%
                          filter(sex == "Female")
df_age_before_male <- df_age_before %>%
                        group_by(age, sex) %>%
                        summarise(n = sum(n)) %>%
                        filter(sex == "Male")


df_age_before_antlerless <- df_age_before_female
df_age_before_antlered <- df_age_before_male

#adding male fawns to anterless females
df_age_before_antlerless$n[df_age_before_antlerless$age == 0] <-
              df_age_before_female$n[df_age_before_female$age == 0] +
              df_age_before_male$n[df_age_before_male$age == 0]

df_age_before_antlered <- df_age_before_male %>% filter(age != 0)
df_age_before_female_nofawn <- df_age_before_female %>% filter(age != 0)

#ratio to calculate sex-ratio for initial population
#  sum(df_age_before_female_nofawn$n)
#  sum(df_age_before_antlered$n)

doe_per <- df_age_before_female_nofawn$n[
                df_age_before_female_nofawn$age == 1] /
                sum(df_age_before_female_nofawn$n)
buck_per <- df_age_before_antlered$n[
                df_age_before_antlered$age == 1] /
                sum(df_age_before_antlered$n)
sex_ratio_early <- doe_per / buck_per
female_proportion_early <- 1 - (sex_ratio_early/(sex_ratio_early + 1))


#proportion of females in the total population size


df_age_before_antlerless$proportion <- 
    df_age_before_antlerless$n / sum(df_age_before_antlerless$n)
df_age_before_antlered$proportion <-
    df_age_before_antlered$n / sum(df_age_before_antlered$n)
df_age_before_female$proportion <-
    df_age_before_female$n / sum(df_age_before_female$n)
df_age_before_male$proportion <-
    df_age_before_male$n / sum(df_age_before_male$n)

########################################
### Hyper prior for  study-area- sex-age
### structured initial population size
### accounting for CWD infection
########################################

###################
### East
###################

prevalence_f_east <- apply(Cage_inf[1, 1, , ], 2, sum) /
                      apply(Cage[1, 1, , 9:28], 2, sum)
prevalence_m_east <- apply(Cage_inf[1, 2, , ], 2, sum) /
                      apply(Cage[1, 2, , 9:28],2,sum)

df_prev_f_east <- data.frame(x = 2002:2021,
                    prevalence_f_east = log(prevalence_f_east))
lm_f_east <- lm(prevalence_f_east ~ x,
                data = df_prev_f_east)
# summary(lm_f_east)
df_prev_m_east <- data.frame(x = 2002:2021, 
                    prevalence_m_east = log(prevalence_m_east))
lm_m_east <- lm(prevalence_m_east ~ x,
                data = df_prev_m_east)
# summary(lm_m_east)

pred_prev_f_east <- exp(predict(lm_f_east,
                      newdata = data.frame(x = 1994:2001)))
pred_prev_m_east <- exp(predict(lm_m_east,
                      newdata = data.frame(x = 1994:2001)))

###################
### West
###################

prevalence_f_west <- apply(Cage_inf[2, 1, , ], 2, sum) /
                     apply(Cage[2, 1, , 9:28], 2, sum)
prevalence_m_west <- apply(Cage_inf[2, 2, , ], 2, sum) /
                     apply(Cage[2, 2, , 9:28], 2, sum)
#tricking prevalence to be super small rather than 0 to fit the log linear model. 
prevalence_f_west[prevalence_f_west == 0] <- .0000001

df_prev_f_west <- data.frame(x = 2002:2021,
                    prevalence_f_west = log(prevalence_f_west))
lm_f_west <- lm(prevalence_f_west ~ x,
                    data = df_prev_f_west)
# summary(lm_f_west)
df_prev_m_west <- data.frame(x = 2002:2021,
                    prevalence_m_west = log(prevalence_m_west))
lm_m_west <- lm(prevalence_m_west ~ x,
                  data = df_prev_m_west)
# summary(lm_m_west)
# 
pred_prev_f_west <- exp(predict(lm_f_west,
                        newdata = data.frame(x = 1994:2001)))
pred_prev_m_west <- exp(predict(lm_m_west,
                        newdata = data.frame(x = 1994:2001)))
# pred_prev_f_east
# pred_prev_m_east
# pred_prev_f_west
# pred_prev_m_west

#############################################################
###
### Hyper prior for  sex-age
### structured initial population size
###
#############################################################

initN_sus <- array(0, dim = c(n_study_area, n_sex, n_agef))
initN_inf <- array(0, dim = c(n_study_area, n_sex, n_agef))

#######################
### East
#######################

### Total population initialized from SAK estimate from 1994
init0_east <- c()
init0_east$female <- female_proportion_early * pop_estimate_east
init0_east$male <- (1 - female_proportion_early) * pop_estimate_east


###susceptible initial population
###females
initN_sus[1, 1, 1] <- init0_east$female *
                   .5 *
                   df_age_before_antlerless$proportion[1] *
                   (1 - pred_prev_f_east[1])   #F
for(a in 2:4){
      initN_sus[1, 1, a] <- init0_east$female *
                         df_age_before_antlerless$proportion[a] *
                         (1 - pred_prev_f_east[1])   #1,2,3
}
initN_sus[1, 1, 5] <- init0_east$female * 
                   df_age_before_antlerless$proportion[5] *
                   (2/3) *
                   (1 - pred_prev_f_east[[1]])  #4
initN_sus[1, 1, 6] <- init0_east$female *
                   df_age_before_antlerless$proportion[5] *
                   (1/3) * (1 - pred_prev_f_east[[1]])  #5
initN_sus[1, 1, 7] <- init0_east$female *
                   df_age_before_antlerless$proportion[6] * 
                   (1/2) * 
                   (1 - pred_prev_f_east[[1]])  #6
initN_sus[1, 1, 8] <- init0_east$female *
                   df_age_before_antlerless$proportion[6] *
                   (1/4) *
                   (1 - pred_prev_f_east[[1]]) #7
initN_sus[1, 1, 9] <- init0_east$female * 
                   df_age_before_antlerless$proportion[6]*
                   (1/6) *
                   (1 - pred_prev_f_east[[1]]) #8
initN_sus[1, 1, 10] <- init0_east$female *
                    df_age_before_antlerless$proportion[6] *
                    (1/12) *
                    (1 - pred_prev_f_east[[1]])  #9+,same as 6-8, decaying from 6-8

###antlered
initN_sus[1, 2, 1] <- init0_east$female *
                   .5 *
                   df_age_before_antlerless$proportion[1] *
                   (1 - pred_prev_m_east[1]) 
for(a in 2:4){
      initN_sus[1, 2, a] <- init0_east$male *
                         df_age_before_antlered$proportion[a - 1] *
                         (1 - pred_prev_m_east[1])  #1,2,3,
}
initN_sus[1, 2, 5] <- 1 * (1 - pred_prev_m_east[1]) # 4
initN_sus[1, 2, 6] <- 1 * (1 - pred_prev_m_east[1]) # 5
initN_sus[1, 2, 7] <- 1 * (1 - pred_prev_m_east[1]) # 6+

####################################
###infected initial population
####################################

###antlerless
initN_inf[1, 1, 1] <- init0_east$female *
                   .5 *
                   df_age_before_antlerless$proportion[1] *
                   (pred_prev_f_east[1])   #F
for(a in 2:4){
      initN_inf[1, 1, a] <- init0_east$female *
                         df_age_before_antlerless$proportion[a] *
                         (pred_prev_f_east[1])   #1,2,3
}
initN_inf[1, 1, 5] <- init0_east$female *
                   df_age_before_antlerless$proportion[5] *
                   (2 / 3) *
                   (pred_prev_f_east[[1]])  #4
initN_inf[1, 1, 6] <- init0_east$female * 
                   df_age_before_antlerless$proportion[5] *
                   (1 / 3) *
                   (pred_prev_f_east[[1]])  #5
initN_inf[1, 1, 7] <- init0_east$female *
                   df_age_before_antlerless$proportion[6] *
                   (1 / 2) *
                   (pred_prev_f_east[[1]])  #6
initN_inf[1, 1, 8] <- init0_east$female *
                   df_age_before_antlerless$proportion[6] *
                   (1 / 4) *
                   (pred_prev_f_east[[1]]) #7
initN_inf[1, 1, 9] <- init0_east$female *
                   df_age_before_antlerless$proportion[6] *
                   (1 / 6) *
                   (pred_prev_f_east[[1]]) #8
initN_inf[1, 1, 10] <- init0_east$female *
                    df_age_before_antlerless$proportion[6] *
                    (1 / 12) *
                    (pred_prev_f_east[[1]])  # 9+, set as same proportion from 6-8
 
###antlered
initN_inf[1, 2, 1] <- init0_east$female *
                   .5 *
                   df_age_before_antlerless$proportion[1] *
                   (pred_prev_m_east[1])  #F
for(a in 2:4) {
      initN_inf[1, 2, a] <- init0_east$male *
                         df_age_before_antlered$proportion[a - 1] *
                         (pred_prev_m_east[1])  #1,2,3
}
initN_inf[1, 2, 5] <- 1  * (pred_prev_m_east[1])# 4
initN_inf[1, 2, 6] <- 1  * (pred_prev_m_east[1])# 5
initN_inf[1, 2, 7] <- 1  * (pred_prev_m_east[1])# 6+

#######################
### West Study Area
#######################

### Total population initialized from SAK estimate from 1994
init0_west <- c()
init0_west$female <- sex_ratio_early * pop_estimate_west
init0_west$male <- (1 - sex_ratio_early) * pop_estimate_west

###susceptible initial population
###females
initN_sus[2, 1, 1] <- init0_west$female *
                   .5 *
                   df_age_before_antlerless$proportion[1] *
                   (1 - pred_prev_f_west[1])   #F
for(a in 2:4){
      initN_sus[2, 1, a] <- init0_west$female *
                         df_age_before_antlerless$proportion[a] *
                         (1 - pred_prev_f_west[1])   #1,2,3
}
initN_sus[2, 1, 5] <- init0_west$female * 
                   df_age_before_antlerless$proportion[5] *
                   (2/3) *
                   (1 - pred_prev_f_west[[1]])  #4
initN_sus[2, 1, 6] <- init0_west$female *
                   df_age_before_antlerless$proportion[5] *
                   (1/3) * (1 - pred_prev_f_west[[1]])  #5
initN_sus[2, 1, 7] <- init0_west$female *
                   df_age_before_antlerless$proportion[6] * 
                   (1/2) * 
                   (1 - pred_prev_f_west[[1]])  #6
initN_sus[2, 1, 8] <- init0_west$female *
                   df_age_before_antlerless$proportion[6] *
                   (1/4) *
                   (1 - pred_prev_f_west[[1]]) #7
initN_sus[2, 1, 9] <- init0_west$female * 
                   df_age_before_antlerless$proportion[6]*
                   (1/6) *
                   (1 - pred_prev_f_west[[1]]) #8
initN_sus[2, 1, 10] <- init0_west$female *
                    df_age_before_antlerless$proportion[6] *
                    (1/12) *
                    (1 - pred_prev_f_west[[1]])  #9+,same as 6-8, decaying from 6-8

###antlered
initN_sus[2, 2, 1] <- init0_west$female *
                   .5 *
                   df_age_before_antlerless$proportion[1] *
                   (1 - pred_prev_m_west[1]) 
for(a in 2:4){
      initN_sus[2, 2, a] <- init0_west$male *
                         df_age_before_antlered$proportion[a - 1] *
                         (1 - pred_prev_m_west[1])  #1,2,3,
}
initN_sus[2, 2, 5] <- 1 * (1 - pred_prev_m_west[1]) # 4
initN_sus[2, 2, 6] <- 1 * (1 - pred_prev_m_west[1]) # 5
initN_sus[2, 2, 7] <- 1 * (1 - pred_prev_m_west[1]) # 6+

####################################
###infected initial population
####################################

###antlerless
initN_inf[2, 1, 1] <- init0_west$female *
                   .5 *
                   df_age_before_antlerless$proportion[1] *
                   (pred_prev_f_west[1])   #F
for(a in 2:4){
      initN_inf[2, 1, a] <- init0_west$female *
                         df_age_before_antlerless$proportion[a] *
                         (pred_prev_f_west[1])   #1,2,3
}
initN_inf[2, 1, 5] <- init0_west$female *
                   df_age_before_antlerless$proportion[5] *
                   (2 / 3) *
                   (pred_prev_f_west[[1]])  #4
initN_inf[2, 1, 6] <- init0_west$female * 
                   df_age_before_antlerless$proportion[5] *
                   (1 / 3) *
                   (pred_prev_f_west[[1]])  #5
initN_inf[2, 1, 7] <- init0_west$female *
                   df_age_before_antlerless$proportion[6] *
                   (1 / 2) *
                   (pred_prev_f_west[[1]])  #6
initN_inf[2, 1, 8] <- init0_west$female *
                   df_age_before_antlerless$proportion[6] *
                   (1 / 4) *
                   (pred_prev_f_west[[1]]) #7
initN_inf[2, 1, 9] <- init0_west$female *
                   df_age_before_antlerless$proportion[6] *
                   (1 / 6) *
                   (pred_prev_f_west[[1]]) #8
initN_inf[2, 1, 10] <- init0_west$female *
                    df_age_before_antlerless$proportion[6] *
                    (1 / 12) *
                    (pred_prev_f_west[[1]])  # 9+, set as same proportion from 6-8

###antlered
initN_inf[2, 2, 1] <- init0_west$female *
                   .5 *
                   df_age_before_antlerless$proportion[1] *
                   (pred_prev_m_west[1])  #F
for(a in 2:4) {
      initN_inf[2, 2, a] <- init0_west$male *
                         df_age_before_antlered$proportion[a - 1] *
                         (pred_prev_m_west[1])  #1,2,3
}
initN_inf[2, 2, 5] <- 1  * (pred_prev_m_west[1])# 4
initN_inf[2, 2, 6] <- 1  * (pred_prev_m_west[1])# 5
initN_inf[2, 2, 7] <- 1  * (pred_prev_m_west[1])# 6+


################################################
f_logpop_sus_east <- log(initN_sus[1, 1,])
f_logpop_inf_east <- log(initN_inf[1, 1,])
m_logpop_sus_east <- log(initN_sus[1, 2,1:n_agem])
m_logpop_inf_east <- log(initN_inf[1, 2,1:n_agem])


f_logpop_sus_west <- log(initN_sus[2, 1,])
f_logpop_inf_west <- log(initN_inf[2, 1,])
m_logpop_sus_west <- log(initN_sus[2, 2,1:n_agem])
m_logpop_inf_west <- log(initN_inf[2, 2,1:n_agem])



f_logpop_sus <- rbind(f_logpop_sus_east,
                      f_logpop_sus_west)

f_logpop_inf <- rbind(f_logpop_inf_east,
                      f_logpop_inf_west)

m_logpop_sus <- rbind(m_logpop_sus_east,
                      m_logpop_sus_west)
                      
m_logpop_inf <- rbind(m_logpop_inf_east,
                      m_logpop_inf_west)

################################################
###
### Setting up the Cage data
###
################################################

#adding antlerless male fawns to the antlerless aging data
Cage_less <- array(NA, c(n_study_area, n_ageclassf + 1, n_year))
Cage_less[1,,] <- rbind(Cage[1,1,,], Cage[1,2,1,])
Cage_less[2,,] <- rbind(Cage[2,1,,], Cage[2,2,1,])

#removing antlered male fawns from the male aging data
Cage_ant <- Cage[,2,2:n_ageclassm,]

sizeCage_f <- apply(Cage_less, c(1,3), sum)
sizeCage_m <- apply(Cage_ant, c(1,3), sum)


#############################################
###
### initial values for reporting rates
###
#############################################

report_overall_init <- rbeta(1, 
                        report_hyp_all[1],
                        report_hyp_all[2])

report_init <- rep(report_overall_init, n_year)
for(y in 23:n_year) {
    report_init[y] <- rbeta(1, report_hyp_y$alpha[y - 22],
                               report_hyp_y$beta[y - 22])
}

#########################################################################
###
### Fecundity/FDR Preliminaries
###
#########################################################################

###################################################################
### Gamma prior for camera trap, poisson for earlier data
###################################################################

fdr_ct_gam_moments_1992_2016 <- gamma.moments(fawndoe_df$overall_fd,
                                  mean(df_camtrap_fd$fdr_sd) * 2)
fdr_ct_gam_moments_2017_2021 <- gamma.moments(df_camtrap_fd$fdr_mean,
                                  df_camtrap_fd$fdr_sd)
obs_ct_fd_alpha  <- c(fdr_ct_gam_moments_1992_2016$alpha,
                      fdr_ct_gam_moments_2017_2021$alpha)
obs_ct_fd_beta <- c(fdr_ct_gam_moments_1992_2016$beta,
                    fdr_ct_gam_moments_2017_2021$beta)
fec_init <- c(fawndoe_df$overall_fd,
              df_camtrap_fd$fdr_mean)
mu_fec_init <-  mean(log(fec_init))
fec_eps_init <- log(fec_init) - mu_fec_init
n_year_fec_early <- nrow(fawndoe_df)

############################################################
###
### EAB prior allocation
###
############################################################

eab_antlerless_alpha <- gamma.moments(1.49, .1225^2)$alpha
eab_antlerless_beta <- gamma.moments(1.49, .1225^2)$beta
eab_antlered_alpha <- gamma.moments(.71, .1015^2)$alpha
eab_antlered_beta <- gamma.moments(.71, .1015^2)$beta

# pdf("figures/eab_prior_plot.pdf")
# hist(rgamma(10000,eab_antlerless_alpha,eab_antlerless_beta))
# hist(rgamma(10000,eab_antlered_alpha,eab_antlered_beta))
# dev.off()

# df_temp <- data.frame(year=1994:2021,
#              tot = Ototal$antlerless,
#              totant = Ototal$antlered,
#              eab = df_eab$EAB)
# df_temp$eab <- as.factor(df_temp$eab)
# eab_indicator_plot <- ggplot(data = df_temp) + 
#       geom_point(aes(x = year,
#                         y = tot,
#                         color = eab),
#                   size = 8) +
#       geom_point(aes(x = year,
#                      y = totant,
#                      color = eab),
#                   size = 8) +
#       geom_line(aes(x = year,y = tot), size = 1) + 
#       geom_line(aes(x = year, y = totant),
#                     size = 1,color = "green4") +
#       theme_bw()

# ggsave("eab_indicator_plot.png", 
#        eab_indicator_plot, height = 6, width = 6)

# eab_indicator_plot



####################################################
###
### setup age at harvest to estimate period effects
###
###################################################

d_fit_aah$birthweek <- ceiling(interval(study_origin,
                                d_fit_aah$birth_date) / weeks(1)) -
                        nT_period_prestudy_ext
d_fit_aah$birthmonth <- ceiling(interval(study_origin,
                                d_fit_aah$birth_date) / months(1)) -
                        nT_period_prestudy_ext_monthly
d_fit_aah$age2date_weeks <- ceiling(interval(study_origin,
                                d_fit_aah$birth_date) / weeks(1))
d_fit_aah$age2date_months <- ceiling(interval(study_origin,
                                d_fit_aah$birth_date) / months(1))

d_fit_aah$study_area <- ifelse(d_fit_aah$study_area == "east", 1, 2)

####################################################
###
### setup age at harvest to estimate period effects
###
###################################################

d_fit_notest$birthweek <- ceiling(interval(study_origin,
                                d_fit_notest$birth_date) / weeks(1)) -
                        nT_period_prestudy_ext
d_fit_notest$birthmonth <- ceiling(interval(study_origin,
                                d_fit_notest$birth_date) / months(1)) -
                        nT_period_prestudy_ext_monthly
d_fit_notest$age2date_weeks <- ceiling(interval(study_origin,
                                d_fit_notest$birth_date) / weeks(1))
d_fit_notest$age2date_months <- ceiling(interval(study_origin,
                                d_fit_notest$birth_date) / months(1))
d_fit_notest$study_area <- ifelse(d_fit_notest$study_area == "east", 1, 2)
d_fit_notest <- d_fit_notest[d_fit_notest$n > 0,]

