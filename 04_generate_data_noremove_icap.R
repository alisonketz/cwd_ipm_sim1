##################################################################
###
### run the function to generate data
###
##################################################################

source("03_fun_generate_data.R")

dat <- ageperiod_surv_foi_sim_data(
            beta0_survival_sus = beta0_survival_sus,
            beta0_survival_inf = beta0_survival_inf,
            foi_age_effect = foi_age_effect,
            age_effect = age_effect_true,
            period_effect = period_effect_true,
            nT_age = nT_age,
            nT_period = nT_period,
            processnum = processnum
            )

# table(dat$pos1)
# table(dat$pos2)
# hist(dat$right_age,breaks = 100)
# hist(dat$right_age)
# sum(dat$rt_censor)
# n_ind

# hist(dat$right_period,breaks = 100)
# dat$right_age
# # prevalence
# sum(!is.na(dat$inf_age))/n_ind
# dat$inf_status[1,]
# which(dat$inf_age[!is.na(inf_age)] < dat$left_age[!is.na(inf_age)])


##############################################################
### there are too many deer infected at to capture
### because we aren't accounting for 
### disease-associated mortality in the generating function
##############################################################
# icap_eval_df <- data.frame(left_age = dat$left_age[dat$pos1==TRUE])
# icap_eval_plot <- ggplot(data = icap_eval_df) +
#       geom_histogram(aes(x = left_age),bins = 50) +
#       theme_bw() +
#       theme(axis.text = element_text(size = 12),
#         axis.title = element_text(size = 16),
#         strip.text = element_text(size = 16),
#         legend.title = element_text(size = 16),
#         legend.text = element_text(size = 14),
#         plot.title = element_text(size = 16)
#   )


# ggsave(paste0("figures/icap_eval_left_age.png"),
#       icap_eval_plot,
#       height = 6,
#       width = 8)

# icap_eval_df$weights <- icap_eval_df$left_age/sum(icap_eval_df$left_age)
# icap_eval_df$weights2 <- (1 - 1/icap_eval_df$left_age)/sum((1 - 1/icap_eval_df$left_age))

# icap_weights_plot <- ggplot(data = icap_eval_df) +
#       geom_point(aes(x = left_age,y = weights), size = 1.3, alpha = .5) +
#       theme_bw() +
#       xlab("Age at Entry") + 
#       theme(axis.text = element_text(size = 12),
#         axis.title = element_text(size = 16),
#         strip.text = element_text(size = 16),
#         legend.title = element_text(size = 16),
#         legend.text = element_text(size = 14),
#         plot.title = element_text(size = 16)
#   )

# icap_weights2_plot <- ggplot(data = icap_eval_df) +
#       geom_point(aes(x = left_age,y = weights2), size = 1.3, alpha = .5) +
#       theme_bw() +
#       xlab("Age at Entry") + 
#       theme(axis.text = element_text(size = 12),
#         axis.title = element_text(size = 16),
#         strip.text = element_text(size = 16),
#         legend.title = element_text(size = 16),
#         legend.text = element_text(size = 14),
#         plot.title = element_text(size = 16)
#   )

# ggsave(paste0("figures/icap_eval_weights.png"),
#       icap_weights_plot,
#       height = 6,
#       width = 8)



# ggsave(paste0("figures/icap_eval_weights2.png"),
#       icap_weights_plot,
#       height = 6,
#       width = 8)


# n_remove <- sum(dat$pos1) - round(sum(dat$pos1)/5)
# indx_remove <- sample(which(dat$pos1==TRUE),n_remove,replace=FALSE) 
# n_ind <- dat$n_ind - n_remove
n_ind <- dat$n_ind


##############################################################
###
### setup overall data frame
###
##############################################################

df_fit <- data.frame(id = 1:n_ind,
                     left_age = dat$left_age,
                     right_age = dat$right_age,
                     left_period = dat$left_period,
                     right_period = dat$right_period,
                     rt_censor = dat$rt_censor,
                     cwd_cap = dat$pos1,
                     cwd_mort = dat$pos2,
                     inf_age = dat$inf_age
                     )

df_fit$inf_period = df_fit$right_age - df_fit$inf_age + df_fit$left_period
n_fit <- nrow(df_fit)
# hist(df_fit$left_period, breaks = 100)
# hist(df_fit$inf_period[df_fit$inf_period>0],breaks = 100)
# sum(df_fitinf_cap)
# sum(1-df_fit$rt_censor)
# sum(df_fit$rt_censor)


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


########################################################################
### calibrating age of deer with the study time 
### for indexing in the likelihood for loops
### age2date = left_period - left_age
########################################################################

df_fit$age2date <- df_fit$e_period - df_fit$e_age
