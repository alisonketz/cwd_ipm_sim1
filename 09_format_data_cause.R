###############################################################
###
### Formatting data to fit the conditional probability that 
### individuals that die during the hunting season die from 
### hunter harvest mortality
###
##############################################################

d_fit_hh <- d_surv[d_surv$censored == 0,]
d_fit_hh[,1:3] <- d_fit_hh[,1:3] - nT_period_precollar_ext

#females are 0
#males are 1
#creating the response for the cause-specific likelihood
#mort_h == 0 means that deer are not hunter harvest
#mort_h == 1 means that deer are hunter harvested 

d_fit_hh$mort_h <- 0
d_fit_hh$mort_h[d_fit_hh$p4 == 1] <- 1
records_cause <- nrow(d_fit_hh)

###########################################
###
### Read data - hunting season dates
###
###########################################

### huntseason is for duration of collar data only
d_huntseason <- read_xlsx(paste0(filepath,"Hunting_SeasonDates.xlsx"),1)
d_huntseason <- d_huntseason %>% filter(Year > 2016)

### season is specified from beginning of 
### the pop model study (May 15, 1994)
### for the full study
d_season <- read_xlsx(paste0(filepath,"Hunting_SeasonDates.xlsx"),1)
d_season <- d_season %>% filter(Year > 1993)


# startdate <- "2017-01-09"
startdate <- "1985-05-15"
n_year <- length(unique(d_season$Year))
n_year_collar <- 5

startng <- c()
endng <- c()
startgun <- c()
endgun <- c()
for (i in 1:n_year_collar) {
    startng[i] <- ceiling(interval(startdate,
        min(d_huntseason$OpenDate[d_huntseason$Year ==
                                  (i + 2016)])) / weeks(1)) -
        nT_period_precollar_ext
    endng[i] <- ceiling(interval(startdate,
        max(d_huntseason$CloseDate[d_huntseason$Year ==
                                  (i + 2016)])) / weeks(1)) -
        nT_period_precollar_ext
    startgun[i] <- ceiling(interval(startdate,
        d_huntseason$OpenDate[d_huntseason$Year ==
                              (i + 2016) &
        d_huntseason$SeasonType == "nineday"]) / weeks(1)) -
        nT_period_precollar_ext
    endgun[i] <- ceiling(interval(startdate,
        d_huntseason$CloseDate[d_huntseason$Year ==
                                (i + 2016) &
        d_huntseason$SeasonType == "nineday"]) / weeks(1)) -
        nT_period_precollar_ext
}

#########################################################################
### these seasons are not quite right and are 1 week too late based on 
### when then hunted deer in the survival collar 
### data were dying from hunter harvest
#########################################################################

# startgun[1] <- startgun[1]-1
# endgun[1] <- endgun[1]-1
# startgun[4] <- startgun[4]-1
# endgun[4] <- endgun[4]-1

#########################################################################
### 
### Season indexing across full study timeline
### 
#########################################################################

# study_start2 <- "1994-05-15"
study_origin <- "1985-05-15"
season_ng_start <- c()
season_ng_end <- c()
for (i in 1:length(unique(d_season$Year))) {
    season_ng_start[i] <- ceiling(interval(study_origin,
        min(d_season$OpenDate[d_season$Year ==
                                  (i + 1993)])) / weeks(1)) - 
        nT_period_prestudy_ext

    season_ng_end[i] <- ceiling(interval(study_origin,
        max(d_season$CloseDate[d_season$Year ==
                                  (i + 1993)])) / weeks(1)) -
        nT_period_prestudy_ext
}

# pulling out approximate 9 day gun season
temp <- d_season %>% filter(SeasonType == "nineday")
nonineday <- c(1994:2021)[which(!(c(1994:2021) %in% temp$Year))]
temp2 <- d_season %>% 
            filter(Year %in% nonineday) %>% 
            filter(SeasonType != "Archery" &
                   SeasonType != "Muzzleloader" &
                   SeasonType != "Youth Gun")
temp2rm  <- c(5,6,8,11,12)
temp2 <- temp2[-temp2rm,]
temp <- rbind(temp,temp2)
temp <- temp[order(temp$Year),]

season_gun_start <- c()
season_gun_end <- c()
for (i in 1:length(unique(temp$Year))) {
    season_gun_start[i] <- ceiling(interval(study_origin,
          temp$OpenDate[temp$Year == (i + 1993)]) / weeks(1)) - 
    nT_period_prestudy_ext

    season_gun_end[i] <-  ceiling(interval(study_origin,
          temp$CloseDate[temp$Year == (i + 1993)]) / weeks(1)) - 
    nT_period_prestudy_ext
}

#########################################################################
### these seasons are not quite right and are 1 week too late based on 
### when then hunted deer in the survival collar 
### data were dying from hunter harvest
#########################################################################

# season_gun_start[n_year_precollar + 1] <-
#         season_gun_start[n_year_precollar + 1] - 1
# season_gun_end[n_year_precollar + 4] <- 
#         season_gun_end[n_year_precollar + 4] - 1

#########################################################################
###
### Setup dataframe for seasons across full study time
###
### need 8 start/end points throughout year
### for the season indexing which corresponds to 
### start of birth pulse period
### then end non-harvest period
### then ng_harvest start
### then gun_harvest start
### gun harvest end
### ng harvest end
### postharvest start
### end of aah annual year, expressed in period effects
###
#########################################################################

end_dates <- paste(1986:2022,"05","14",sep="-")
start_date <- "1985-05-15"
interval_cuts <- floor(interval(start_date,end_dates)/weeks(1))
(interval_step <- c(diff(interval_cuts),52))
length(interval_step)
intvl_step_yearly <- interval_step[(n_year_prestudy_ext + 1):n_year_ext]

d_fit_season <- matrix(NA, nrow = n_year, ncol = 8)
d_fit_season[1, ]  <-  c( 1,
                          season_ng_start[1] - 1,
                          season_ng_start[1],
                          season_gun_start[1],
                          season_gun_end[1],
                          season_ng_end[1],
                          season_ng_end[1] + 1,
                          intvl_step_yearly[1])
for(t in 2:n_year) {
   d_fit_season[t, ] <- c(sum(intvl_step_yearly[1:(t-1)]) + 1,
                          season_ng_start[t] - 1,
                          season_ng_start[t],
                          season_gun_start[t],
                          season_gun_end[t],
                          season_ng_end[t],
                          season_ng_end[t] + 1,
                          sum(intvl_step_yearly[1:(t)]))#this is changed
}

d_fit_season <- data.frame(d_fit_season)
colnames(d_fit_season) <- c("yr_start",
                            "pre_hunt_end",
                            "ng_start",
                            "gun_start",
                            "gun_end",
                            "ng_end",
                            "post_hunt_start",
                            "yr_end")

# d_fit_season$yr_end <- cumsum(intvl_step_yearly)
#saving for aah_disease test
save(d_fit_season, file = paste0(filepath, "d_fit_season.Rdata"))

d_fit_season$year <- 1994:2021

##########################################################################
###
### Preliminaries gun season
###
##########################################################################

Z_cause_ng <- rep(0,nT_period_collar)
Z_cause_gun <- rep(0,nT_period_collar)
for(i in 1:5){
    Z_cause_ng[startng[i]:endng[i]] <- 1
    Z_cause_gun[startgun[i]:endgun[i]] <- 1
}

Z_overall_ng <- rep(0,nT_period_overall_ext)
Z_overall_gun <- rep(0,nT_period_overall_ext)
for(i in 1:n_year){
    Z_overall_ng[(d_fit_season$ng_start[i]:
                  d_fit_season$ng_end[i]) + nT_period_prestudy_ext] <- 1
    Z_overall_gun[(d_fit_season$gun_start[i]:
                   d_fit_season$gun_end[i]) + nT_period_prestudy_ext] <- 1
}


Z_collar_ng <- Z_cause_ng
Z_collar_ng[Z_cause_gun==1] <- 0 
Z_collar_gun <- Z_cause_gun

###
### this should be used in the likelihoods to ensure use
### the time varying hunterharvest only collar data to estimate the
###

Z_overall_gun_collar <- Z_overall_gun
Z_overall_gun_collar[1:nT_period_precollar_ext] <- 0

Z_overall_ng_collar <- Z_overall_ng
Z_overall_ng_collar[1:nT_period_precollar_ext] <- 0

################################################################
###
### preliminaries for hunter harvest Z at a monthly time step
###
################################################################

mort_h <- d_fit_hh$mort_h
interval_cause <- d_fit_hh$right_period_s - 1

### debugging hunting season dates
# table(interval_cause[d_fit_hh$mort_h==1])
# which(Z_cause_gun==1)
# interval_cause[which(Z_cause_gun==1)]
# length(interval_cause[which(interval_cause %in% which(Z_cause_gun==1))])
# sum(mort_h[which(interval_cause %in% which(Z_cause_gun==1))])
# length(interval_cause[which(interval_cause %in% which(Z_cause_ng==1))])
# sum(mort_h[which(interval_cause %in% which(Z_cause_ng==1))])
# length(interval_cause[which(interval_cause %in% which(Z_cause_gun==1))])
# interval_cause[which(interval_cause %in% which(Z_cause_ng==1))]

low_gun <- d_mort$lowtag[d_mort$weapon == "rifle"]
low_gun <- low_gun[!is.na(low_gun)]
low_bow <- d_mort$lowtag[d_mort$weapon == "bow" |
                         d_mort$weapon == "crossbow" |
                         d_mort$weapon == "compound bow" |
                         d_mort$weapon == "long bow" |
                         d_mort$weapon == "shotgun"]
low_bow <- low_bow[!is.na(low_gun)]

d_fit_hh$gun[d_fit_hh$lowtag %in% low_gun] <- 1
d_fit_hh$bow[d_fit_hh$lowtag %in% low_bow] <- 1
d_fit_hh$gun[which(interval_cause %in% which(Z_cause_gun==1))]

# names(d_fit_hh)
# which(apply(d_fit_hh[d_fit_hh$mort_h==1,7:8],1,sum)==0)
# d_fit_hh[d_fit_hh$mort_h==1,][c(5,85),]
# 5736
# 6915
# d_fit_hh[d_fit_hh$mort_h,]
# names(d_fit_hh)
# table(interval_cause[d_fit_hh$gun==1])
# table(interval_cause[d_fit_hh$bow==1])
# which(Z_cause_gun==1)
# startgun
# endgun
 

 ##################################################################
 ###
 ### Figuring out mistakes in the data,
 ### which deer were gun harvest not during the nineday gun season
 ### are these correct? 
 ###
 ##################################################################
# test1 <- which(d_fit_hh$lowtag %in% low_gun)
# test2 <- which(interval_cause %in% which(Z_cause_gun==1))
# test1[which(!(test1 %in% test2))]
# testlow <- d_fit_hh$lowtag[test1[which(!(test1 %in% test2))]]
# testlowdf <- d_mort[which(d_mort$lowtag %in% testlow),c(2,14:16,55)]
# write.csv(testlowdf,file="results/testlowgun.csv")
