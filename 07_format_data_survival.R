#########################################
###
### setting the number of periods
### for collar data,
### overall,
### precollar study periods
###
#########################################

###########################
### Weekly
###########################

#setting the maximum number of periods that are estimated by the collar data
nT_period_overall_ext <- ceiling(interval("1985-05-15","2022-05-14") / weeks(1))
nT_period_precollar_ext <- floor(interval("1985-05-15","2017-01-08") / weeks(1))
nT_period_prestudy_ext <- floor(interval("1985-05-15","1994-05-14") / weeks(1))
nT_period_collar <- nT_period_overall_ext - nT_period_precollar_ext
nT_period_overall <- nT_period_overall_ext - nT_period_prestudy_ext
nT_period_precollar <- nT_period_overall - nT_period_collar

#setting the maximum number of periods that are estimated by the collar data
# nT_period_collar <- interval("2017-01-09", "2022-05-14") %/% weeks(1)
# nT_period_collar_monthly <- interval("2017-01-09", "2022-05-14") %/% months(1)

# #the first birth is in 1985
# # weekly calculation
# #from first birth in 1985 until end of study in 2022
# nT_period_overall_ext <- interval("1985-05-15","2022-05-14") %/% weeks(1)

# #from the start of the population model in May 1994 until the end of the study May 2022
# nT_period_overall <- interval("1994-05-15", "2022-05-14") %/% weeks(1)

# #number of weekly intervals prior to the first collared individual
# #from the start of the pop model in 1994
# nT_period_precollar <- nT_period_overall - nT_period_collar

# #number of weekly intervals prior to the start of the pop model in 1994
# nT_period_prestudy_ext <- nT_period_overall_ext - nT_period_overall

# #number of weekly intervals from the birth of the first deer in 1985
# #until the start of the collar study. 
# nT_period_precollar_ext <- nT_period_overall_ext - nT_period_collar

# nT_period_prestudy_ext + nT_period_overall

###########################
### Monthly
###########################

 

###########################

### Monthly

###########################


nT_period_overall_ext_monthly <- ceiling(interval("1985-05-15","2022-05-14") / months(1))
nT_period_precollar_ext_monthly <- floor(interval("1985-05-15","2017-01-08") / months(1))
nT_period_prestudy_ext_monthly <- floor(interval("1985-05-15","1994-05-14") / months(1))
nT_period_collar_monthly <- nT_period_overall_ext_monthly - nT_period_precollar_ext_monthly
nT_period_overall_monthly <- nT_period_overall_ext_monthly - nT_period_prestudy_ext_monthly
nT_period_precollar_monthly <- nT_period_overall_monthly - nT_period_collar_monthly

# and then every date of interest is defined in terms of weeks as:
# week_period <- <- ceiling(interval("1985-05-15","xxxx-xx-xx") / weeks(1))

# nT_period_overall_ext_monthly <- interval("1985-05-15", "2022-05-14") %/% months(1)
# # nT_period_overall_ext_monthly <- interval("1992-05-15","2022-05-14") %/% months(1)
# nT_period_overall_monthly <- interval("1994-05-15","2022-05-14") %/% months(1)
# nT_period_precollar_ext_monthly <- nT_period_overall_ext_monthly - nT_period_collar

###########################
### Years
###########################

n_year_precollar <- 2017 - 1994
n_year_precollar_ext <- 2017 - 1985
n_year_prestudy_ext <- 1994 - 1985
n_year_ext <- 2022 - 1985

########################################
### spreadsheet to confirm timeline
########################################

#the first birth is in 1985
# weekly calculation
timing_df <- data.frame(var = c("nT_period_collar",
    "nT_period_overall_ext",
    "nT_period_overall ",
    "nT_period_precollar ",
    "nT_period_prestudy_ext ",
    "nT_period_precollar_ext ",
    "nT_period_prestudy_ext + nT_period_overall",
    "n_year",
    "n_year_ext",
    "n_year_precollar",
    "n_year_precollar_ext",
    "n_year_prestudy_ext"),
    value = c(nT_period_collar,
    nT_period_overall_ext,
    nT_period_overall,
    nT_period_precollar,
    nT_period_prestudy_ext,
    nT_period_precollar_ext,
    nT_period_prestudy_ext + nT_period_overall,
    n_year,
    n_year_ext,
    n_year_precollar,
    n_year_precollar_ext,
    n_year_prestudy_ext),
    formula = c("interval(2017-01-09, 2022-05-14) %/% weeks(1)",
    "interval(1985-05-15,2022-05-14) %/% weeks(1)",
    "interval(1994-05-15, 2022-05-14) %/% weeks(1)",
    "nT_period_overall - nT_period_collar",
    "nT_period_overall_ext - nT_period_overall",
    "nT_period_overall_ext - nT_period_collar",
    "nT_period_prestudy_ext + nT_period_overall",
    "2022 - 1994",
    "2022 - 1985",
    "2017 - 1994",
    "2017 - 1985",
    "1994 - 1985")
)
write.csv(timing_df,file="results/timing_params.csv")

########################################################
###
### Formatting/combining of data
###
########################################################
###
### d_surv[,1] = e
### d_surv[,2] = r
### d_surv[,3] = s
### d_surv[,4] = censor
### d_surv[,5] = CWD status at capture
### d_surv[,6] = sex, sex=0=females, sex=1=males
### d_surv[,7] = gun harvest
### d_surv[,8] = bow harvest
### d_surv[,9] = lowtag
### d_surv[,10] = year of study captured
### d_surv[,11] = CWD status at mortality
### d_surv[,12] = e age at capture (weeks)
### d_surv[,13] = r age (weeks)
### d_surv[,14] = s age (weeks)
### d_surv[,15:21] prob of different causes of mort
### d_surv[,22] = e agemonths at capture
### d_surv[,23] = ageright r months
### d_surv[,24] = ageright s months
### d_surv[,25] = cwd status antemortem test at recapture
### d_surv[,26] = age in weeks at recapture
### d_surv[,27] = age in months at recapture
### d_surv[,28] = period week at recapture
### d_surv[,29] = period week at recapture
### d_surv[,30] =  "dsection"
### d_surv[,31] =  "study_area"

n_cap <- nrow(df_cap)
n_mort <- nrow(d_mort)
n_cens <- nrow(d_cens)

#initialize matrix of captures/recaptures/morts
d_surv <- data.frame(matrix(NA, nr = n_cap, nc = 31))
names(d_surv)  <- c("left_period_e","right_period_r","right_period_s",
                "censored",
                "cwd_cap",
                "sex",
                "gun",
                "bow",
                "lowtag",
                "year_cap",
                "cwd_mort",
                "left_age_e",
                "right_age_r",
                "right_age_s",
                "p1","p2","p3","p4","p5","p6","p7",
                "emonth",
                "rmonth",
                "smonth",
                "recap_cwd",
                "ageweek_recap",
                "agemonth_recap",
                "periodweek_recap",
                "periodmonth_recap",
                "dsection",
                "study_area"
                )

#initialize with all animals alive
d_surv$censored <- 1

#initialize with 0's al remaining columns
d_surv[,5:ncol(d_surv)] <- 0

#Filling d_surv with all times
d_surv$left_period_e  <-  df_cap$eweek#e
d_surv$cwd_cap <-  df_cap$cwd_cap
d_surv$sex  <-  df_cap$sex
d_surv$lowtag  <-  df_cap$lowtag
d_surv$year_cap <-  df_cap$year_cap
d_surv$left_age_e <- df_cap$ageweek_cap
d_surv$emonth <- df_cap$emonth
d_surv$recap_cwd <- df_cap$recap_cwd
d_surv$ageweek_recap <- df_cap$recap_disagewk
d_surv$agemonth_recap <- df_cap$recap_disagemth
d_surv$periodweek_recap <- df_cap$recap_disweek
d_surv$periodmonth_recap <- df_cap$recap_dismonth
d_surv$dsection <- df_cap$home_section
d_surv$study_area <- ifelse(df_cap$ew=="east",1,2)

for(i in 1:n_cap) {
  for(j in 1:n_mort) {
    if(df_cap$lowtag[i] == d_mort$lowtag[j]) {
      d_surv$right_period_r[i] <- d_mort$rweek[j]
      d_surv$right_period_s[i] <- d_mort$sweek[j]
      d_surv$censored[i] <- 0
      d_surv$p1[i] <- d_mort$pc1[j]
      d_surv$p2[i] <- d_mort$pc2[j]
      d_surv$p3[i] <- d_mort$pc3[j]
      d_surv$p4[i] <- d_mort$pc4[j]
      d_surv$p5[i] <- d_mort$pc5[j]
      d_surv$p6[i] <- d_mort$pc6[j]
      d_surv$p7[i] <- d_mort$pc7[j]
      d_surv$rmonth[i] <- d_mort$rmonth[j]
      d_surv$smonth[i] <- d_mort$smonth[j]
      if (d_mort$cause1[j] == "hunter harvest") {
        if(d_mort$weapon[j] == "rifle" |
           d_mort$weapon[j] == "shotgun" |
           d_mort$weapon[j] == "muzzleloader") {
          d_surv[i,7] <- 1
        } else {
          d_surv[i,8] <- 1
        }
      }
    }
  }
}
n_surv <- dim(d_surv)[1]
n_surv

###
### censored deer
###

for(i in 1:n_cens){
  d_surv$right_period_r[d_surv$lowtag %in% d_cens$lowtag[i]] <- d_cens$rweek[i]
  d_surv$rmonth[d_surv$lowtag %in% d_cens$lowtag[i]] <- d_cens$rmonth[i]
}

#there's 81 alive after study
low_endlive <- d_surv[which(is.na(d_surv[,2])),9]
nrow(d_surv[which(is.na(d_surv[,2])),])
nrow(df_cap) - (nrow(d_cens) + nrow(d_mort))

#Right censor these
d_surv[which(is.na(d_surv[, 2])), 2] <- nT_period_overall_ext
d_surv$rmonth[which(d_surv$rmonth == 0)] <- nT_period_overall_ext_monthly
d_surv$smonth[which(d_surv$smonth == 0)] <- nT_period_overall_ext_monthly

###
### converting to a data.frame
###

for(i in 1:5) {
  d_surv[, i] <- as.integer(d_surv[, i])
}
for(i in 7:13) {
  d_surv[, i] <- as.integer(d_surv[, i])
}
d_surv$sex <- ifelse(d_surv$sex == "Female", 0, 1)

###
### Setting CWD status at mortality
###

d_surv$cwd_mort <- NA
for (i in 1:n_postcwd) {
    d_surv$cwd_mort[d_surv$lowtag %in% d_post_cwd$lowtag[i]] <- d_post_cwd$cwdresult[i]
}
d_surv$cwd_mort <- as.factor(d_surv$cwd_mort)
levels(d_surv$cwd_mort) <- c("Negative",NA,"Positive")
levels(d_surv$cwd_mort) <- c(0, 1)

d_surv$cwd_mort <- as.numeric(as.character(d_surv$cwd_mort))
d_surv$cwd_mort[which(d_surv$cwd_cap == 1)] <- 1

###
### right age in weeks
###

d_surv$right_age_r <- d_surv$left_age_e + d_surv$right_period_r - d_surv$left_period_e
d_surv$right_age_s <- d_surv$left_age_e + d_surv$right_period_s - d_surv$left_period_e

min_age <- min(d_surv$left_age_e)
max_ager <- max(d_surv$right_age_r, na.rm = TRUE)
max_ages <- max(d_surv$right_age_s, na.rm = TRUE)

###
### left_age/right age in months
###

d_surv$left_age_month <- df_cap$agemonth_cap

#fixing fast mortalities in months
fix_month_indx <- which(d_surv$rmonth < d_surv$emonth)
d_surv$rmonth[fix_month_indx] <- d_surv$smonth[fix_month_indx]
d_surv$smonth[fix_month_indx] <- d_surv$smonth[fix_month_indx] + 1

d_surv$right_age_rmonth <- d_surv$left_age_month + d_surv$rmonth - d_surv$emonth
d_surv$right_age_smonth <- d_surv$left_age_month + d_surv$smonth - d_surv$emonth

#constants
n_surv <- dim(d_surv)[1]
n_mort <- dim(d_mort)[1]
n_cens <- dim(d_cens)[1]

#fixing fast right censored individuals
d_surv$right_age_r[which(d_surv$right_period_r - d_surv$left_period_e < 1 & d_surv$censored == 1)] <- 
  d_surv$right_age_r[which(d_surv$right_period_r - d_surv$left_period_e < 1 & d_surv$censored == 1)] + 1
d_surv$right_period_r[which(d_surv$right_period_r - d_surv$left_period_e < 1 & d_surv$censored == 1)] <- 
  d_surv$right_period_r[which(d_surv$right_period_r - d_surv$left_period_e < 1 & d_surv$censored == 1)] + 1

#fixing fast mortalities individuals
d_surv$right_age_s[which(d_surv$right_period_s - d_surv$left_period_e < 1 & d_surv$censored == 0)] <- 
  d_surv$right_age_s[which(d_surv$right_period_s - d_surv$left_period_e < 1 & d_surv$censored == 0)] + 1
d_surv$right_period_s[which(d_surv$right_period_s - d_surv$left_period_e < 1 & d_surv$censored == 0)] <- 
  d_surv$right_period_s[which(d_surv$right_period_s - d_surv$left_period_e < 1 & d_surv$censored == 0)] + 1

fix_fast_mortalities_indx <- which(d_surv$right_age_r == 0)
d_surv$right_period_r[fix_fast_mortalities_indx] <- d_surv$left_period_e[fix_fast_mortalities_indx]
d_surv$right_age_r[fix_fast_mortalities_indx] <- 1

# there are 2 fast mortalities for non-neonates
# these need to have r fixed too, where r should equal e
# the lowtags for these are 5080, 5219
fix_fast_mortalities_nonneonate_indx <- which(d_surv$lowtag %in% c(5080,5219))
d_surv$right_period_r[fix_fast_mortalities_nonneonate_indx] <- d_surv$left_period_e[fix_fast_mortalities_nonneonate_indx]
d_surv$right_age_r[fix_fast_mortalities_nonneonate_indx] <- d_surv$left_age_e[fix_fast_mortalities_nonneonate_indx] 

##################################################
###
### removing the recapture info for the deer 
### that were recaptured after right censoring 
###
### 6817 6876 5153
###
##################################################

censor_fix_low <- c(6817,6876,5153)
# d_surv[d_surv$lowtag == 6817,]
# d_surv[d_surv$lowtag == 6876,]
# d_surv[d_surv$lowtag == 5153,]

#removing these from recaptured deer
#just using the status at recapture as same as status at right censor
rm_censor_fix <- which(low_recap %in% censor_fix_low)
low_recap <- low_recap[-rm_censor_fix]

rm_censor_fix_neg <- which(low_recap_neg %in% censor_fix_low)
low_recap_neg <- low_recap_neg[-rm_censor_fix_neg]


##################################################
###
### removing the deer that were VERY fast
### right censors. i.e., lost collars within
### first week after capture, these don't contribute
### any information to estimate the effects
###
### 5052 6400 6081 5257 7113 7787
###
##################################################

low_remove_fast_cens <- c(5052, 6400, 6081, 5257, 7113, 7787)

d_surv <- d_surv[!(d_surv$lowtag %in% low_remove_fast_cens), ]
n_surv <- nrow(d_surv)

##########################################################################
###
### calibrating collar study time with start of harvest study time
### rescaling the origin of the  period effects
### from the collar study start to overall
### pop model start
###
##########################################################################

# #recalibrating periods of collar data
# #weekly adjustment

# d_surv$left_period_e <- d_surv$left_period_e 
# d_surv$right_period_r <- d_surv$right_period_r 
# d_surv$right_period_s <- d_surv$right_period_s 
# d_surv$periodweek_recap[d_surv$periodweek_recap != 0] <-
#           d_surv$periodweek_recap[d_surv$periodweek_recap != 0] +
#           nT_period_precollar_ext

# #monthly adjustment
# d_surv$emonth <- d_surv$emonth 
# d_surv$rmonth <- d_surv$rmonth 
# d_surv$smonth <- d_surv$smonth 
# d_surv$periodmonth_recap[d_surv$periodmonth_recap != 0] <-
#           d_surv$periodmonth_recap[d_surv$periodmonth_recap != 0] +
#           nT_period_precollar_ext_monthly
