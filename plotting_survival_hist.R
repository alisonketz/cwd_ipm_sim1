head(d_surv)

d_surv_out = data.frame(d_surv)
d_surv_out$censored = as.factor(d_surv_out$censored)
levels(d_surv_out$censored) = c("mortality","right_censored")
right_age_plot = ggplot(data=d_surv_out, aes(x = right_period_r,fill = censored)) +
                    geom_histogram() +
                    theme_bw()
ggsave("figures/right_age_censored.png",right_age_plot)
