###############################################################################################
###
### Post processing
###
###############################################################################################

# source("summarize.R")
# load("mcmcout.Rdata")
# load("runtime.Rdata")
# out <- mcmc.list(mcmcout1)
# fit_sum <- summarize(out)

# fit_sum <- mcmcout$summary
fit_sum <- mcmcout$summary$all.chains
out <- mcmcout$samples

modelid <- "A"

#############################
### Saving Model Description
#############################
sink(paste0("figures/", modelid, "/model_description_", modelid, ".txt"))
cat("Model description specifics:\n")
cat("niter:  ",reps,"\n")
cat("burnin:  ",bin,"\n")
cat("n_chains:  ",n_chains,"\n")
cat("Model Variation:\n")
cat("no FOI period effects\n")
cat("includes FOI age effects RW2 with implicit interept \n\n")
cat("includes survival period effects kernel convolution\n")
cat("includes survival age effects cgam convex \n\n")
cat("runtime:  ", runtime, "\n")
cat("Gelman diag:")
print(gelman.diag(out[,grep("beta",rownames(fit_sum))],multivariate = FALSE))
print(gelman.diag(out[,grep("tau",rownames(fit_sum))],multivariate = FALSE))
print(gelman.diag(out[,grep("sd",rownames(fit_sum))],multivariate = FALSE))
print(gelman.diag(out[,grep("foi_age_effect",rownames(fit_sum))],multivariate = FALSE))
cat("\nSummary Stats:  \n")
print(fit_sum)
sink()

#############################
### from single run
#############################

pdf(paste0("figures/",modelid,"/traceplots_",format(Sys.time(),"%y%m%d%m%s"),"_",modelid,".pdf"))
traceplot(out[, "tau_age_foi"], ylab = "tau_age_foi")
traceplot(out[, "beta0_survival_sus"], ylab = "beta0_survival_sus")
traceplot(out[, "beta0_survival_inf"], ylab = "beta0_survival_inf")
traceplot(out[, "tau_age_survival"], ylab = "tau_age_survival")
traceplot(out[, "sdk_period"], ylab = "sdk_period")
traceplot(out[, "tauk_period"], ylab = "tauk_period")
traceplot(out[, "stauk_period"], ylab = "stauk_period")
traceplot(out[, "sda_period"], ylab = "sda_period")
traceplot(out[, "taua_period"], ylab = "taua_period")
dev.off()

png(paste0("figures/",modelid,"/beta0_survival_sus_traceplot_",modelid,".png"))
traceplot(out[, "beta0_survival_sus"], ylab = "beta0_survival_sus")
dev.off()

png(paste0("figures/",modelid,"/beta0_survival_inf_traceplot_",modelid,".png"))
traceplot(out[, "beta0_survival_inf"], ylab = "beta0_survival_inf")
dev.off()

png(paste0("figures/",modelid,"/beta0_survival_sus_densityplot_",modelid,".png"))
densityplot(out[, "beta0_survival_sus"], ylab = "beta0_survival_sus")
dev.off()

png(paste0("figures/",modelid,"/beta0_survival_inf_densityplot_",modelid,".png"))
densityplot(out[, "beta0_survival_inf"], ylab = "beta0_survival_inf")
dev.off()

pdf(paste0("figures/",modelid,"/traceplot_foi_age_",modelid,".pdf"))
for(i in 1:n_ageclass){
    traceplot(out[, paste0("foi_age_effect[",i,"]")], ylab = paste0("foi_age_effect[",i,"]"))
}
dev.off()

pdf(paste0("figures/",modelid,"/traceplot_ln_b_age_survival_",modelid,".pdf"))
for(i in 1:nknots_age){
    traceplot(out[, paste0("ln_b_age_survival[",i,"]")], ylab = paste0("ln_b_age_survival[",i,"]"))
}
dev.off()

# pdf(paste0("figures/",modelid,"/traceplot_b_age_survival_",modelid,".pdf"))
# for(i in 1:nknots_age){
#     traceplot(out[, paste0("b_age_survival[",i,"]")], ylab = paste0("b_age_survival[",i,"]"))
# }
# dev.off()

#########################################################
### Color pallets
#########################################################

renoir_pal <- met.brewer(name="Renoir", n=10, type="discrete")
pillement_pal <- met.brewer(name="Pillement", n=6, type="discrete")
troy_pal <- met.brewer(name="Troy", n=8, type="discrete")


###############################################
###
### Plots of age effects for mortality hazard
### lineplot
###
###############################################

foi_age_indx <- grep("foi_age_effect",rownames(fit_sum))

age_effect_mean <- fit_sum[foi_age_indx,2]
age_effect_lower <- fit_sum[foi_age_indx,4]
age_effect_upper <- fit_sum[foi_age_indx,5]

ageclass <- 1:n_ageclass

foi_age_effect_out <- data.frame(ageclass,age_effect_mean,age_effect_lower,age_effect_upper)

foi_age_effect_plot <- ggplot(data =foi_age_effect_out, aes(x = ageclass))+
  geom_line(aes(x = ageclass,y=age_effect_mean), size = 1)+
  geom_ribbon(aes(ymin = age_effect_lower,
                  ymax = age_effect_upper),
              alpha = .2,
              linetype = 0)+
  ggtitle("Age Effect Posterior")+xlab("Age (Years)")+ylab("Effect Size")+
  theme_bw()#+
  #scale_x_continuous(breaks = seq(0,nT_age,by=104),labels=seq(0,n_year,by=2))+
  #scale_color_manual("Year",values = met.brewer("Kandinsky", 2)) +
  #scale_fill_manual("Year",values = met.brewer("Kandinsky", 2))
  # theme(axis.text.x = element_text(angle = 90, hjust = 1))

foi_age_effect_plot

# ggsave("figures/foi_age_effect.pdf",foi_age_effect_plot)
ggsave(paste0("figures/",modelid,"/foi_age_effect_",modelid,".png"),foi_age_effect_plot)


###################################################
###
### Force of infection age plots - step
###
##################################################


foi_age_indx <- grep("foi_age_effect",rownames(fit_sum))

foi_age_effect_out <- data.frame(agegroups = 1:n_ageclass,
                   age_effect_mean = fit_sum[foi_age_indx,2],
                   age_effect_lower = fit_sum[foi_age_indx,4],
                   age_effect_upper = fit_sum[foi_age_indx,5])

df_age_foi <- rbind(foi_age_effect_out,foi_age_effect_out[nrow(foi_age_effect_out),])
df_age_foi$agegroups[nrow(df_age_foi)] <- df_age_foi$agegroups[nrow(df_age_foi)] + 1

df_age_foi$agegroups <- df_age_foi$agegroups - 1

df_age_foi$truth <- c(age_foi,age_foi[length(age_foi)])


foi_age_step_plot <- ggplot(data=df_age_foi)+
  geom_rect(aes(xmin=agegroups,
                xmax=lead(agegroups),
                ymin=age_effect_lower,
                ymax=age_effect_upper),
            fill=pillement_pal[2],alpha=.4)+
  geom_step(aes(x=agegroups,y=age_effect_mean),size = 1.5, color = pillement_pal[2])+
  geom_step(aes(x=agegroups,y=truth),linetype = "dotted",size = 1)+
  ggtitle("Female Force of Infection")+
  theme_bw()+xlab("Age Class")+ylab("Age Effects Weekly Conversion Hazard (Log)")+
  scale_x_continuous(breaks = (1:n_ageclassf)-.5,labels =c("Fawn","1.5","2.5","3.5","4.5-5.5","6.5-8.5","9.5+"))+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        plot.title = element_text(size = 16)
  )

foi_age_step_plot

ggsave(paste0("figures/",modelid,"/foi_age_step_",modelid,".png"),
      foi_age_step_plot,
      height = 6,
      width = 8)

###########################################################
###
### beta0_survival_inf/sus
###
############################################################
# beta0_survival_sus <- out[,grep("beta0_survival_sus",rownames(fit_sum))]
# beta0_survival_inf <- out[,grep("beta0_survival_inf",rownames(fit_sum))]

beta0_survival_sus <- do.call(c,out[,grep("beta0_survival_sus",rownames(fit_sum))])
beta0_survival_inf <- do.call(c,out[,grep("beta0_survival_inf",rownames(fit_sum))])

df_beta0_survival <- data.frame(beta0_survival_sus,beta0_survival_inf) %>% 
                    pivot_longer(cols = everything())
names(df_beta0_survival)[1] <- "cwd_status"
df_beta0_survival$cwd_status <- as.factor(df_beta0_survival$cwd_status)
levels(df_beta0_survival$cwd_status) <- c("Susceptible","Infected")
# df_beta0_survival$cwd_status <- factor(df_beta0_survival$cwd_status, levels = c("Susceptible","Infected"))
# max(df_beta0_survival$value)

beta0_survival_plot <- ggplot(df_beta0_survival)+
    geom_density(aes(x=value,fill = cwd_status,color = cwd_status), alpha = .8)+
    theme_bw() +
    ylab("Density")+xlab("Coefficient")+
    ggtitle("Posterior Intercepts Log Hazard Mortality")+
    labs(fill = "CWD Status",color = "CWD Status")+
    scale_fill_manual(values = troy_pal[c(7, 2)])+
    scale_color_manual(values = troy_pal[c(7, 2)])+
    geom_vline(xintercept = beta0_survival_inf_true,linetype = "dashed") +
    geom_vline(xintercept = beta0_survival_sus_true,linetype = "dotted") +
    theme(axis.text=element_text(size=14),
            axis.title=element_text(size=16),
            title =element_text(size=18),
            strip.text.x = element_text(size = 12),
            legend.title = element_text(size=12))

beta0_survival_plot

ggsave(paste0("figures/",modelid,"/beta0_survival_plot","_",modelid,".png"),
      beta0_survival_plot,
      height = 6,
      width = 10)

fit_sum[grep("beta0_survival_",rownames(fit_sum)),]
write.csv(fit_sum[grep("beta0_survival_",rownames(fit_sum)),],file = paste0("results/",modelid,"/beta0_survival_summary.csv"))


###############################################
###
### Plot of period effects for mortality hazard
###
###############################################

te_indx <- grep("period_effect_surv",rownames(fit_sum))
out_period_effect <- data.frame(weeks = 1:nT_period,
                                period_effect_mean = fit_sum[te_indx,1],
                                period_effect_lower = fit_sum[te_indx,4],
                                period_effect_upper = fit_sum[te_indx,5],
                                truth = period_effect)

period_effect_plot <- ggplot(data =out_period_effect,aes(x = weeks))+
  geom_ribbon(aes(ymin=period_effect_lower,ymax=period_effect_upper),alpha=.2,linetype=0)+
  geom_line(aes(x = weeks,y=period_effect_mean),size=1)+
  geom_line(aes(x = weeks,y=truth),size=1,color=troy_pal[3])+
  ggtitle("Log Mortality Hazard Period Effect ")+xlab("Time")+ylab("Log Mortality Hazard")+
  theme_bw() +
  scale_x_continuous(breaks = seq(0,nT_period,by=52),labels=paste0("Jan ", 2017:2022)) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        axis.text.x = element_text(angle = 45,hjust = 1),
        plot.title = element_text(size = 16)
)
  # scale_x_continuous(breaks = seq(0,inf_nT_period,by=104),labels=seq(0,18,by=2))+
  # theme(axis.text.x = element_text(angle = 45, hjust = 1))

period_effect_plot

# ggsave("figures/period_effect_inf_int.pdf",period_effect_plot_inf_int)
ggsave(paste0("figures/",modelid,"/period_effects_",modelid,".png"),period_effect_plot, height = 6, width = 10)

###############################################
###
### Plot of age effects for mortality hazard
###
###############################################

ae_indx <- grep("age_effect_survival",rownames(fit_sum))

out_age_effect <- data.frame(weeks = 1:nT_age,
                            age_effect_mean = fit_sum[ae_indx,2],
                            age_effect_lower =  fit_sum[ae_indx,4],
                            age_effect_upper = fit_sum[ae_indx,5],
                            truth = age_effect)

age_effect_plot <- ggplot(data =out_age_effect,aes(x = weeks))+
  geom_ribbon(aes(ymin=age_effect_lower,ymax=age_effect_upper),alpha=.2,linetype=0)+
  geom_line(aes(x = weeks,y=age_effect_mean),size=1,alpha = .6)+
  geom_line(aes(x = weeks,y=truth),size=1,linetype = "dotted")+
  ggtitle("Log Mortality Hazard Age Effect ")+xlab("Age") + ylab("Log Mortality Hazard")+
  theme_bw() +
  scale_x_continuous(breaks = seq(104,nT_age,by=104),labels=seq(2,19,by =2)) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        plot.title = element_text(size = 16)
)
age_effect_plot
ggsave(paste0("figures/",modelid,"/age_effect_",modelid,".png"),age_effect_plot, height = 6, width = 8)

##################################################
##################################################
###
### save mcmcout 
###
###################################################
##################################################

save(mcmcout,file = paste0("results/",modelid,"/mcmcout_",modelid,".Rdata"))
save(fit_sum,file = paste0("results/",modelid,"/fit_sum_",modelid,".Rdata"))


#############################
### psi
#############################

# psi_out <- data.frame(psi = fit_sum[grep("psi", rownames(fit_sum)), 1],
#     study_area = rep(c("east", "west"), 20),
#     sex = rep(c("Female", "Female", "Male", "Male"), 10),
#     age = rep(1:10, each = 4)
#     )
# psi_out$age <- as.factor(psi_out$age)
# psi_out <- psi_out[!(psi_out$sex == "Male" & psi_out$age %in% c(8,9,10)),] 


# psi_plot <- ggplot(data = psi_out) + geom_point(aes(x=age,y=psi)) + facet_nested_wrap(~ study_area + sex)+ theme_bw()
# ggsave(paste0("figures/",modelid,"/psi_mn_",modelid,".png"), psi_plot)


# # #############################
# # ### sn_sus
# # #############################
# # nrow(fit_sum[grep("sn_sus", rownames(fit_sum)), ]) == 100

# sn_sus_out <- data.frame(
#     sn_sus = fit_sum[grep("sn_sus", rownames(fit_sum)), 1],
#     sex = rep(c("Female", "Male"), 50),
#     age = rep(rep(1:10, each = 2), n_year),
#     year = rep(2017:2021, each = 20)
# )
# sn_sus_out$age <- as.factor(sn_sus_out$age)
# sn_sus_out <- sn_sus_out[!(sn_sus_out$sex == "Male" & sn_sus_out$age %in% c(8,9,10)),] 

# sn_sus_plot <- ggplot(data = sn_sus_out) +
#                geom_point(aes(x=year,y=sn_sus,color=age)) +
#                facet_wrap(~ sex) +
#                theme_bw()
# ggsave(paste0("figures/",modelid,"/sn_sus_mn_",modelid,".png"), sn_sus_plot,height=4,width=7)


# #############################
# ### sn_inf
# #############################

# sn_inf_out <- data.frame(
#     sn_inf = fit_sum[grep("sn_inf", rownames(fit_sum)), 1],
#     sex = rep(c("Female", "Male"), 50),
#     age = rep(rep(1:10, each = 2), n_year),
#     year = rep(2017:2021, each = 20)
# )
# sn_inf_out$age <- as.factor(sn_inf_out$age)
# sn_inf_out <- sn_inf_out[!(sn_inf_out$sex == "Male" & sn_inf_out$age %in% c(8,9,10)),] 

# sn_inf_plot <- ggplot(data = sn_inf_out) + geom_point(aes(x=year,y=sn_inf,color=age)) + facet_wrap(~ sex)+ theme_bw()
# ggsave(paste0("figures/",modelid,"/sn_inf_mn_",modelid,".png"), sn_inf_plot,height=4,width=7)

# #############################
# ### sh_sus
# #############################

# sh_sus_out <- data.frame(sh_sus = fit_sum[grep("sh_sus", rownames(fit_sum)), 1],
#     sex = rep(c("Female", "Male"), 50),
#     age = rep(rep(1:10, each = 2), n_year),
#     year = rep(2017:2021, each = 20)
# )
# sh_sus_out$age <- as.factor(sh_sus_out$age)
# sh_sus_out <- sh_sus_out[!(sh_sus_out$sex == "Male" & sh_sus_out$age %in% c(8,9,10)),] 


# sh_sus_plot <- ggplot(data = sh_sus_out) + geom_point(aes(x=year,y=sh_sus,color=age)) + facet_wrap(~ sex)+ theme_bw()
# ggsave(paste0("figures/",modelid,"/sh_sus_mn_",modelid,".png"), sh_sus_plot,height=4,width=7)

# #############################
# ### sh_inf
# #############################

# sh_inf_out <- data.frame(sh_inf = fit_sum[grep("sh_inf", rownames(fit_sum)), 1],
# sex = rep(c("Female", "Male"), 50),
#     age = rep(rep(1:10, each = 2), n_year),
#     year = rep(2017:2021, each = 20)
# )
# sh_inf_out$age <- as.factor(sh_inf_out$age)
# sh_inf_out <- sh_inf_out[!(sh_inf_out$sex == "Male" & sh_inf_out$age %in% c(8,9,10)),] 


# sh_inf_plot <- ggplot(data = sh_inf_out) + geom_point(aes(x=year,y=sh_inf,color=age)) + facet_wrap(~ sex)+ theme_bw()
# ggsave(paste0("figures/",modelid,"/sh_inf_mn_",modelid,".png"), sh_inf_plot,height=4,width=7)


# pdf(paste0("figures/",modelid,"/psi_surv_plots_mn_",modelid,".pdf"))
# psi_plot
# sn_sus_plot
# sn_inf_plot
# sh_sus_plot
# sh_inf_plot
# dev.off()

##################################################
##################################################
##################################################
##################################################
##################################################
##################################################
##################################################
##################################################

# modelid <- "N"

# load(paste0("results/",modelid,"/mcmcout_",modelid,".Rdata"))

# # fit_sum <- mcmcout$summary
# fit_sum <- mcmcout$summary$all.chains
# out <- mcmcout$samples

# gelman.diag(out[,grep("beta",rownames(fit_sum))],multivariate = FALSE)
# gelman.diag(out[,grep("tau",rownames(fit_sum))],multivariate = FALSE)
# gelman.diag(out[,grep("sd",rownames(fit_sum))],multivariate = FALSE)
# gelman.diag(out[,grep("f_age",rownames(fit_sum))],multivariate = FALSE)
# gelman.diag(out[,grep("m_age",rownames(fit_sum))],multivariate = FALSE)

# ##########################################################
# #plots of annual survival
# ##########################################################

# age_effect_survival <- do.call(rbind,out[,grep("age_effect",rownames(fit_sum))])
# period_effect_survival <- do.call(rbind,out[,grep("period_effect_surv",rownames(fit_sum))])
# beta0_survival_sus <- do.call(c,out[,grep("beta0_survival_sus",rownames(fit_sum))])
# beta0_survival_inf <- do.call(c,out[,grep("beta0_survival_inf",rownames(fit_sum))])
# tot_iter <- nrow(age_effect_survival)
# sn_sus <-  do.call(rbind,out[,grep("sn_sus",rownames(fit_sum))])
# sn_inf <-  do.call(rbind,out[,grep("sn_inf",rownames(fit_sum))])
# sh_sus <-  do.call(rbind,out[,grep("sh_sus",rownames(fit_sum))])
# sh_inf <-  do.call(rbind,out[,grep("sh_inf",rownames(fit_sum))])
# psi <-  do.call(rbind,out[,grep("psi",rownames(fit_sum))])



# ##################################################################
# ###
# ### Clean plot of annual survival - Susceptibles (sn_sus)
# ###
# ###################################################################

# n_df <- length(fit_sum[grep("sn_sus",rownames(fit_sum)),1])

# df_sn_sus  <- data.frame(survival = fit_sum[grep("sn_sus",rownames(fit_sum)),1],
#                          lower = fit_sum[grep("sn_sus",rownames(fit_sum)),4],                    
#                          upper = fit_sum[grep("sn_sus",rownames(fit_sum)),5],
#                          sex = rep(c("Female","Male"),n_df/2),
#                          age = rep(1:10,each = 2,n_df/20),
#                          year = rep(2017:2021, each = n_df/5)
#                         )

# df_sn_sus <- df_sn_sus %>% filter(!is.na(survival))
# df_sn_sus$age <- as.factor(df_sn_sus$age)
# levels(df_sn_sus$age) <- c("Fawn",as.character(1:8),"9+")

# out_plot_sus <- ggplot(data=df_sn_sus) +
#                 geom_point(aes(x= year,y=survival,color = age),size = 4,alpha = .6) +
#                 facet_wrap(~sex) +
#                 theme_bw()+
#                 ylab("Annual Survival Probability May 15-May14 the following year")+
#                 xlab("Year")+
#                 ylim(0,1)+
#                 scale_color_manual("Age",values=renoir_pal)+
#                 theme(axis.text = element_text(size = 14),
#                       axis.title = element_text(size = 16),
#                       strip.text = element_text(size = 16),
#                       legend.title = element_text(size = 16),
#                       legend.text = element_text(size = 14)
#                 )
# out_plot_sus
# ggsave(paste0("figures/",modelid,"/annual_survival_sn_sus","_",modelid,".png"),out_plot_sus,height = 7, width = 10)

# df_sn_sus_f <- df_sn_sus %>% filter(sex == "Female")
# f_out_plot_sus_ci <- ggplot(data=df_sn_sus_f) +
#                 geom_point(aes(x= year,y=survival,color = age),size = 3,alpha = .6) +
#                 geom_errorbar(aes(x= year,ymin=lower, ymax=upper,color = age), width=.3)+
#                 facet_wrap(~age, nrow  = 2) +
#                 theme_bw()+
#                 ylab("Annual Survival Probability")+
#                 xlab("Year")+
#                 ylim(0,1)+
#                 scale_color_manual("Age",values=renoir_pal)+
#                 theme(axis.text = element_text(size = 12),
#                       axis.title = element_text(size = 16),
#                       strip.text = element_text(size = 16),
#                       legend.title = element_text(size = 16),
#                       legend.text = element_text(size = 14),
#                       axis.text.x = element_text(angle = 45,hjust = 1)
#                 )
# f_out_plot_sus_ci

# ggsave(paste0("figures/",modelid,"/annual_survival_sn_sus_f","_",modelid,".png"),
#               f_out_plot_sus_ci,height = 6, width = 10.5)


# df_sn_sus_m <- df_sn_sus %>% filter(sex == "Male")
# df_sn_sus_m$age <- as.factor(as.character(df_sn_sus_m$age))
# levels(df_sn_sus_m$age) <- c(as.character(1:5),"6+","Fawn")
# df_sn_sus_m$age <- factor(df_sn_sus_m$age,levels = c("Fawn",as.character(1:5),"6+"))

# m_out_plot_sus_ci <- ggplot(data=df_sn_sus_m) +
#                 geom_point(aes(x= year,y=survival,color = age),size = 3,alpha = .6) +
#                 geom_errorbar(aes(x= year,ymin=lower, ymax=upper,color = age), width=.3)+
#                 facet_wrap(~age, nrow  = 2) +
#                 theme_bw()+
#                 ylab("Annual Survival Probability")+
#                 xlab("Year")+
#                 ylim(0,1)+
#                 scale_color_manual("Age",values=renoir_pal)+
#                 theme(axis.text = element_text(size = 12),
#                       axis.title = element_text(size = 16),
#                       strip.text = element_text(size = 16),
#                       legend.title = element_text(size = 16),
#                       legend.text = element_text(size = 14),
#                       axis.text.x = element_text(angle = 45,hjust = 1)
#                 )
# m_out_plot_sus_ci
# ggsave(paste0("figures/",modelid,"/annual_survival_sn_sus_m","_",modelid,".png"),m_out_plot_sus_ci,height = 6, width = 10.5)


# ##################################################################
# ###
# ### Clean plot of annual survival - Infected (sn_inf)
# ###
# ###################################################################

# n_df <- length(fit_sum[grep("sn_inf",rownames(fit_sum)),1])

# df_sn_inf  <- data.frame(survival = fit_sum[grep("sn_inf",rownames(fit_sum)),1],
#                          lower = fit_sum[grep("sn_inf",rownames(fit_sum)),4],                    
#                          upper = fit_sum[grep("sn_inf",rownames(fit_sum)),5],
#                          sex = rep(c("Female","Male"),n_df/2),
#                          age = rep(1:10,each = 2,n_df/20),
#                          year = rep(2017:2021, each = n_df/5)
#                         )

# df_sn_inf <- df_sn_inf %>% filter(!is.na(survival) & age != 1)
# df_sn_inf$age <- as.factor(df_sn_inf$age)

# levels(df_sn_inf$age)  <- c(as.character(1:8),"9+")


# out_plot_inf <- ggplot(data=df_sn_inf) +
#                 geom_point(aes(x= year,y=survival,color = age),size = 4,alpha = .6) +
#                 facet_wrap(~sex) +
#                 theme_bw()+
#                 ylab("Annual Survival Probability May 15-May14 the following year")+
#                 xlab("Year")+
#                 ylim(0,.5)+
#                 scale_color_manual("Age",values=renoir_pal)+
#                 theme(axis.text = element_text(size = 14),
#                       axis.title = element_text(size = 16),
#                       strip.text = element_text(size = 16),
#                       legend.title = element_text(size = 16),
#                       legend.text = element_text(size = 14)
#                 )
# out_plot_inf
# ggsave(paste0("figures/",modelid,"/annual_survival_sn_inf","_",modelid,".png"),out_plot_inf,height = 7, width = 10)

# df_sn_inf_f <- df_sn_inf %>% filter(sex == "Female")
# f_out_plot_inf_ci <- ggplot(data=df_sn_inf_f) +
#                 geom_point(aes(x= year,y=survival,color = age),size = 3,alpha = .6) +
#                 geom_errorbar(aes(x= year,ymin=lower, ymax=upper,color = age), width=.3)+
#                 facet_wrap(~age, nrow  = 2) +
#                 theme_bw()+
#                 ylab("Annual Survival Probability")+
#                 xlab("Year")+
#                 ylim(0,.5)+
#                 scale_color_manual("Age",values=renoir_pal)+
#                 theme(axis.text = element_text(size = 12),
#                       axis.title = element_text(size = 16),
#                       strip.text = element_text(size = 16),
#                       legend.title = element_text(size = 16),
#                       legend.text = element_text(size = 14),
#                       axis.text.x = element_text(angle = 45,hjust = 1)
#                 )
# f_out_plot_inf_ci

# ggsave(paste0("figures/",modelid,"/annual_survival_sn_inf_f","_",modelid,".png"),
#               f_out_plot_inf_ci,height = 6, width = 10.5)


# df_sn_inf_m <- df_sn_inf %>% filter(sex == "Male")
# df_sn_inf_m$age <- as.factor(as.character(df_sn_inf_m$age))
# levels(df_sn_inf_m$age) <- c(as.character(1:5),"6+")
 
# m_out_plot_inf_ci <- ggplot(data=df_sn_inf_m) +
#                 geom_point(aes(x= year,y=survival,color = age),size = 3,alpha = .6) +
#                 geom_errorbar(aes(x= year,ymin=lower, ymax=upper,color = age), width=.3)+
#                 facet_wrap(~age, nrow  = 2) +
#                 theme_bw()+
#                 ylab("Annual Survival Probability")+
#                 xlab("Year")+
#                 ylim(0,.5)+
#                 scale_color_manual("Age",values=renoir_pal)+
#                 theme(axis.text = element_text(size = 12),
#                       axis.title = element_text(size = 16),
#                       strip.text = element_text(size = 16),
#                       legend.title = element_text(size = 16),
#                       legend.text = element_text(size = 14),
#                       axis.text.x = element_text(angle = 45,hjust = 1)
#                 )
# m_out_plot_inf_ci

# ggsave(paste0("figures/",modelid,"/annual_survival_sn_inf_m","_",modelid,".png"),
#               m_out_plot_inf_ci,height = 6, width = 10.5)



# ##############################################
# ###
# ### infection probability sex/age specific
# ###
# ###############################################

# n_df <- length(fit_sum[grep("psi",rownames(fit_sum)),1])
# df_psi  <- data.frame(infection_prob = fit_sum[grep("psi",rownames(fit_sum)),1],
#                          lower = fit_sum[grep("psi",rownames(fit_sum)),4],                    
#                          upper = fit_sum[grep("psi",rownames(fit_sum)),5],
#                          study_area = rep(c("East","West"),n_df/2),
#                          sex = rep(c(rep("Female",2),rep("Male",2)),n_df/4),
#                          age = rep(1:10,each = 4)
#                         )

# df_psi <- df_psi %>% filter(infection_prob != 0)
# df_psi$age <- as.factor(df_psi$age)
# levels(df_psi$age) <- c("Fawn",as.character(1:8),"9+")

# out_plot_psi <- ggplot(data=df_psi) +
#                 geom_point(aes(x= age,y=infection_prob,color = sex),size = 4,alpha = .8) +
#                 geom_errorbar(aes(x= age,ymin=lower, ymax=upper,color = sex), width=.3)+
#                 facet_nested_wrap(~sex + study_area,nrow = 2)+
#                 theme_bw()+
#                 ylab("Incidence")+
#                 xlab("Age Class")+
#                 ylim(0,1)+
#                 scale_color_manual("Age",values=pillement_pal)+
#                 theme(axis.text = element_text(size = 14),
#                       axis.title = element_text(size = 16),
#                       strip.text = element_text(size = 16),
#                       legend.position = "n"
#                 )
# out_plot_psi


# ggsave(paste0("figures/",modelid,"/incidence_psi","_",modelid,".png"),out_plot_psi,height = 7, width = 10)


# df_psi_f  <- df_psi %>% filter(sex == "Female")
# out_plot_psi_f <- ggplot(data=df_psi_f) +
#                 geom_point(aes(x= age,y=infection_prob,color = sex),size = 4,alpha = .8) +
#                 geom_errorbar(aes(x= age,ymin=lower, ymax=upper,color = sex), width=.3)+
#                 facet_wrap(~study_area)+
#                 theme_bw()+
#                 ylab("Incidence")+
#                 xlab("Age Class")+
#                 ylim(0,0.3)+
#                 scale_color_manual("Age",values=pillement_pal)+
#                 theme(axis.text = element_text(size = 14),
#                       axis.title = element_text(size = 16),
#                       strip.text = element_text(size = 16),
#                       legend.position = "n"
#                 )
# out_plot_psi_f
# ggsave(paste0("figures/",modelid,"/incidence_psi_f","_",modelid,".png"),out_plot_psi_f,height = 7, width = 10)


# df_psi_m <- df_psi %>% filter(sex == "Male")
# df_psi_m$age <- as.factor(as.character(df_psi_m$age))
# levels(df_psi_m$age) <- c(as.character(1:5),"6+","Fawn")

# # df_psi_m$age <- c("Fawn",as.character(1:5),"6+")
# df_psi_m$age <- factor(df_psi_m$age,levels = c("Fawn",as.character(1:5),"6+"))
# out_plot_psi_m <- ggplot(data=df_psi_m) +
#                 geom_point(aes(x= age,y=infection_prob,color = sex),size = 4,alpha = .8) +
#                 geom_errorbar(aes(x= age,ymin=lower, ymax=upper,color = sex), width=.3)+
#                 facet_wrap(~study_area)+
#                 theme_bw()+
#                 ylab("Incidence")+
#                 xlab("Age Class")+
#                 ylim(0,1)+
#                 scale_color_manual("Age",values=pillement_pal)+
#                 theme(axis.text = element_text(size = 14),
#                       axis.title = element_text(size = 16),
#                       strip.text = element_text(size = 16),
#                       legend.position = "n"
#                 )
# out_plot_psi_m
# ggsave(paste0("figures/",modelid,"/incidence_psi_m","_",modelid,".png"),out_plot_psi_m,height = 7, width = 10)

# ####################################################################################
# ###
# ### Apparent prevalence over time
# ###
# #######################################################################################

# df_prev <- cwd_df %>%
#               filter(kill_year>2016)%>%
#               group_by(age, sex, ew,kill_year) %>%
#               summarise(n_cases = n(), .groups = 'drop',
#               n_pos = sum(teststatus==1),
#               prevalence = sum(teststatus==1)/n()) %>% 
#               filter(prevalence !=1)

# df_prev$kill_year[df_prev$kill_year==2022] <- 2021

# df_prev$sex <- as.factor(df_prev$sex)
# levels(df_prev$sex) <- c("Female","Male")
# df_prev$sex <- factor(df_prev$sex,levels = c("Female","Male"))

# df_prev$ew <- as.factor(df_prev$ew)
# levels(df_prev$ew) <- c("East","West")

# prev_plot  <- ggplot(df_prev,aes(fill=age,y=prevalence,x=kill_year))+
#   geom_bar(position="dodge",stat="identity")+
#   scale_fill_manual("Age Class",values=met.brewer("Veronese",7),labels=c("Fawn","1.5","2.5","3.5","4.5-5.5","6.5-8.5","9.5+"))+
#   theme_bw()+facet_nested_wrap(~sex+ew)+
#   ylab("Apparent Prevalence")+xlab("Year")+labs("Age Class")+
#   theme(axis.text=element_text(size=12),
#         axis.title=element_text(size=14),
#         title =element_text(size=12),
#         strip.text.x = element_text(size = 12),
#         legend.title = element_text(size=12))#+
#   # ggtitle("Core Area - Apparent Prevalence")
# prev_plot
# ggsave(paste0("figures/",modelid,"/prevalence","_",modelid,".png"),prev_plot,height = 6, width = 10)


# #############################################
# ###
# ### apparent prevalence w/o year 2017-2021
# ###
# #############################################

# df_prev_noyear <- cwd_df %>%
#             #   filter(age == 3) %>%
#               filter(kill_year>2016)%>%
#               group_by(age,sex,ew) %>%
#               summarise(n_cases = n(), .groups = 'drop',
#               n_pos = sum(teststatus==1),
#               prevalence = sum(teststatus==1)/n())

# df_prev_noyear$sex <- as.factor(df_prev_noyear$sex)
# levels(df_prev_noyear$sex) <- c("Female","Male")
# df_prev_noyear$sex <- factor(df_prev_noyear$sex,levels = c("Male","Female"))

# df_prev_noyear$ew <- as.factor(df_prev_noyear$ew)
# levels(df_prev_noyear$ew) <- c("East","West")


# df_prev_noyear$age <- as.factor(df_prev_noyear$age)
# levels(df_prev_noyear$age) <-c("Fawn","1.5","9.5+","2.5","3.5","4.5-5.5","6.5-8.5") 
# df_prev_noyear$age <- factor(df_prev_noyear$age,levels = c("Fawn","1.5","2.5","3.5","4.5-5.5","6.5-8.5","9.5+"))
# levels(df_prev_noyear$age)

# prev_plot_noyear  <- ggplot(df_prev_noyear,aes(y=prevalence,x=age,fill = age))+
#   geom_bar(position="dodge",stat="identity")+
#   scale_fill_manual("Age Class",values=met.brewer("Veronese",7))+
#   theme_bw()+facet_nested_wrap(~sex+ew)+
#   ylim(0,1)+
#   ylab("Apparent Prevalence")+xlab("Age Class")+
#   theme(axis.text=element_text(size=12),
#         axis.title=element_text(size=14),
#         title =element_text(size=12),
#         strip.text.x = element_text(size = 12),
#         legend.title = element_text(size=12))#+
# prev_plot_noyear

# ggsave(paste0("figures/",modelid,"/prevalence_noyear","_",modelid,".png"),prev_plot_noyear
# ,height = 6, width = 10)

#   # ggtitle("Core Area - Apparent Prevalence")

# df_prev_noyear %>% filter(ew =="East" & sex == "Male") %>% 
#     write_csv(path = "results/male_east_prev.csv")

# df_prev_noyear %>% filter(ew =="East" & sex == "Female") %>% 
#     write_csv(path = "results/female_east_prev.csv")

# df_prev_noyear %>% filter(ew =="West" & sex == "Female") %>% 
#     write_csv(path = "results/female_west_prev.csv")

# df_prev_noyear %>% filter(ew =="West" & sex == "Male") %>% 
#     write_csv(path = "results/male_west_prev.csv")


# df_psi_f %>% filter(study_area =="East") %>% 
#     select(age,sex,study_area,infection_prob) %>% 
#     write_csv(path = paste0("results/",modelid,"/","female_east_psi.csv"))

# df_psi_f %>% filter(study_area =="West") %>% 
#     select(age,sex,study_area,infection_prob) %>% 
#     write_csv(path = paste0("results/",modelid,"/","female_west_psi.csv"))

# df_psi_m %>% filter(study_area =="East") %>% 
#     select(age,sex,study_area,infection_prob) %>% 
#     write_csv(path = paste0("results/",modelid,"/","male_east_psi.csv"))

# df_psi_m %>% filter(study_area =="West") %>% 
#     select(age,sex,study_area,infection_prob) %>% 
#     write_csv(path = paste0("results/",modelid,"/","male_west_psi.csv"))



# d_surv_nofawn <- d_surv %>% filter(left_age_e>30)
# names(d_surv_nofawn)

# head(df_cap)

# data.frame(sex="Combined",
#     prev_capture = sum(d_surv_nofawn$cwd_cap)/nrow(d_surv_nofawn),
#            prev_mort = sum(d_surv_nofawn$cwd_mort,na.rm=TRUE)/nrow(d_surv_nofawn))%>% 
#     write_csv(path = paste0("results/prevalence_collars.csv"))


# d_surv_nofawn_f <- d_surv %>% filter(left_age_e>30) %>% filter(sex==0)
# data.frame(sex = "Female",
#     prev_capture = sum(d_surv_nofawn_f$cwd_cap)/nrow(d_surv_nofawn_f),
#            prev_mort = sum(d_surv_nofawn_f$cwd_mort,na.rm=TRUE)/nrow(d_surv_nofawn_f))%>% 
#     write_csv(path = paste0("results/prevalence_collars_f.csv"))

# d_surv_nofawn_m <- d_surv %>% filter(left_age_e>30) %>% filter(sex==1)

# data.frame(sex = "Male",
#     prev_capture = sum(d_surv_nofawn_m$cwd_cap)/nrow(d_surv_nofawn_m),
#            prev_mort = sum(d_surv_nofawn_m$cwd_mort,na.rm=TRUE)/nrow(d_surv_nofawn_m))%>% 
#     write_csv(path = paste0("results/prevalence_collars_m.csv"))    


# rbind(data.frame(sex="Combined",
#     prev_capture = sum(d_surv_nofawn$cwd_cap)/nrow(d_surv_nofawn),
#            prev_mort = sum(d_surv_nofawn$cwd_mort,na.rm=TRUE)/nrow(d_surv_nofawn)),
# data.frame(sex = "Female",
#     prev_capture = sum(d_surv_nofawn_f$cwd_cap)/nrow(d_surv_nofawn_f),
#            prev_mort = sum(d_surv_nofawn_f$cwd_mort,na.rm=TRUE)/nrow(d_surv_nofawn_f)),
# data.frame(sex = "Male",
#     prev_capture = sum(d_surv_nofawn_m$cwd_cap)/nrow(d_surv_nofawn_m),
#            prev_mort = sum(d_surv_nofawn_m$cwd_mort,na.rm=TRUE)/nrow(d_surv_nofawn_m)))%>% 
#     write_csv(path = paste0("results/prevalence_collars_all.csv"))    



# ###prime age males
# table(d_surv_nofawn_m$right_age_r)

# plot(d_surv_nofawn_m$left_age_e)


#  prev_capture = sum(d_surv_nofawn_m$cwd_cap)/nrow(d_surv_nofawn_m),
#            prev_mort = sum(d_surv_nofawn_m$cwd_mort,na.rm=TRUE)/nrow(d_surv_nofawn_m))


# d_surv_hh <- d_surv_nofawn %>% filter(p4==1)

# prev_mort_hh = sum(d_surv_hh$cwd_mort,na.rm=TRUE)/nrow(d_surv_hh)
# prev_mort_hh_ = sum(d_surv_hh$cwd_mort,na.rm=TRUE)/nrow(d_surv_hh)

# d_surv_hh_f <- d_surv_nofawn %>% filter(p4==1) %>% filter(sex==0)
# d_surv_hh_f <- d_surv_nofawn %>% filter(p4==1) %>% filter(sex==0)

