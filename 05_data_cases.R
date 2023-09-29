##############################################################
###
### Partition by data types by mortality/censor and infection
###
##############################################################

#infected at capture, mortality
df_fit_icap_mort <- df_fit[df_fit$cwd_cap==1 & df_fit$rt_censor==0,]
n_icap_mort <- nrow(df_fit_icap_mort)

#infected at capture, right censored
df_fit_icap_cens <- df_fit[df_fit$cwd_cap==1 & df_fit$rt_censor==1,]
n_icap_cens <- nrow(df_fit_icap_cens)

#neg at capture, infected at mort
df_fit_idead <- df_fit[which(df_fit$cwd_cap == 0 &
                  df_fit$cwd_mort == 1 &
                  df_fit$rt_censor == 0), ]
n_idead <- nrow(df_fit_idead)

#neg at capture, neg at mort
df_fit_sus_mort_posttest <- df_fit[which(df_fit$cwd_cap == 0 &
                              df_fit$cwd_mort == 0 &
                              df_fit$rt_censor == 0), ]
n_sus_mort_posttest <- nrow(df_fit_sus_mort_posttest)


# neg at capture
# all of these deer are right censored at the end of the study
# so we don't actually know their post-censoring test status

df_fit_sus_cens_postno <- df_fit[which(df_fit$cwd_cap == 0 &
                              df_fit$rt_censor == 1), ]

n_sus_cens_postno <- nrow(df_fit_sus_cens_postno)




################################################
### checking that these cases add up to the 
### total number of generated individuals
################################################

n_sus_mort_posttest +
n_idead +
n_sus_cens_postno +
n_icap_cens +
n_icap_mort
