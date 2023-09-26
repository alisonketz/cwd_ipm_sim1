####################################################################################
###
### Structuring data to fit into model
###
#####################################################################################


### Number of age classes and sex classes
n_study_area <- 2
n_year_ext <- 2022 - 1885
n_year <- 2022 - 1994
n_ageclass <- 7
n_ageclassm <- 6
n_ageclassf <- 7
n_agem <- 7
n_agef <- 10
n_sex <- 2

#structuring classification data to fit into the model
Cage_sus <- array(NA, c(n_study_area, n_sex, n_ageclassf, n_year))
for(j in 1:(n_year)) {
    Cage_sus[1,1,,j] <- df_age_sus$n[df_age_sus$year == (1993 + j) &
                                    df_age_sus$sex == "Female" &
                                    df_age_sus$study_area == "east"]
    Cage_sus[2,1,,j] <- df_age_sus$n[df_age_sus$year == (1993 + j) &
                                    df_age_sus$sex == "Female" &
                                    df_age_sus$study_area == "west"]
    Cage_sus[1,2,,j] <- df_age_sus$n[df_age_sus$year == (1993 + j) &
                                    df_age_sus$sex == "Male" &
                                    df_age_sus$study_area == "east"]
    Cage_sus[2,2,,j] <- df_age_sus$n[df_age_sus$year == (1993 + j) &
                                    df_age_sus$sex == "Male" &
                                    df_age_sus$study_area == "west"]
}
# Cage_sus[2,2,,]

#structuring classification data to fit into the model
n_yr_surveillance <- length(unique(df_age_inf$year))
Cage_inf <- array(NA,c(n_study_area,n_sex,n_ageclassf,n_yr_surveillance))
for(j in 1:n_yr_surveillance){
    Cage_inf[1,1,,j] <- df_age_inf$n[df_age_inf$year == (2001 + j) &
                                    df_age_inf$sex == "Female" &
                                    df_age_inf$study_area == "east"]
    Cage_inf[2,1,,j] <- df_age_inf$n[df_age_inf$year == (2001 + j) &
                                    df_age_inf$sex == "Female" &
                                    df_age_inf$study_area == "west"]
    Cage_inf[1,2,,j] <- df_age_inf$n[df_age_inf$year == (2001 + j) &
                                    df_age_inf$sex == "Male" &
                                    df_age_inf$study_area == "east"]
    Cage_inf[2,2,,j] <- df_age_inf$n[df_age_inf$year == (2001 + j) &
                                    df_age_inf$sex == "Male" &
                                    df_age_inf$study_area == "west"]
}
Cage <- Cage_sus
Cage[, , , 9:28] <- Cage[, , , 9:28] + Cage_inf

#Aggregating the oldest age class for males into the next oldest age
Cage[, 2, 6, ] <- Cage[, 2, 6, ] + Cage[, 2, 7,]
Cage[, 2, 7, ] <- 0

####################################################################################
###
### Total harvest data
###
#####################################################################################

### if not using data by gun/bow types:
names(df_harv_overall_total)[1] <- "year"

Ototal_east <- df_harv_overall_total[df_harv_overall_total$year > 1993 &
                                     df_harv_overall_total$study_area == "east",]

Ototal_west <-  df_harv_overall_total[df_harv_overall_total$year > 1993 &
                                     df_harv_overall_total$study_area == "west",]

Ototal <- array(NA,c(n_study_area,n_sex,nrow(Ototal_east)))
Ototal[1,,] <- t(Ototal_east[,4:3])
Ototal[2,,] <- t(Ototal_west[,4:3])

####################################################################################
### Loading and cleaning harvest compliance rate data
#####################################################################################

report_df <- suppressWarnings(read_excel("datafiles/ComplianceRate2020.xlsx",sheet=7))
report_hyp_sum <- apply(report_df[,2:3],2,mean)

report_hyp_all <- unlist(beta.moments(report_hyp_sum[1], report_hyp_sum[2]))
report_hyp_y <- matrix(NA, nrow(report_df), 2)

for(i in 1:nrow(report_df)){
    report_hyp_y[i,] <- unlist(beta.moments(report_df$compliance_rate[i],report_df$se[i]))
}
report_hyp_y <- data.frame(report_hyp_y)
names(report_hyp_y) <- c("alpha", "beta")

####################################################################################
###
### Loading and cleaning fawn:doe ratio estimates data
###
#####################################################################################

fawndoe_df <- read.csv(paste0(filepath,"fawndoe_1997_2017.csv"),header = TRUE)

#calculating overall fawn:doe ratios across all three counties
fawndoe_df$overall_doe <- fawndoe_df$dane_num_doe +
                          fawndoe_df$iowa_num_doe +
                          fawndoe_df$grant_num_doe
fawndoe_df$overall_fawn <- fawndoe_df$dane_num_fawn +
                           fawndoe_df$iowa_num_fawn +
                           fawndoe_df$grant_num_fawn
fawndoe_df$overall_fd <- fawndoe_df$overall_fawn / fawndoe_df$overall_doe

#Restricting to the years of the study
fawndoe_df <- fawndoe_df[fawndoe_df$year < 2017, ]

#2017-2021
df_camtrap_fd <- read.csv(paste0(filepath, "Iowa_FDR_2017-2021_with_sd.csv"), header = TRUE)

#reading data from 1992-2015
fd_older_df <- read_excel(paste0(filepath, "SW_FDR_1992-2015.xlsx"), 1)
fd_older_df  <- fd_older_df %>% filter(year>1991 & year < 1997)

names(fd_older_df) <- c("spatial.unit", "year", "overall_fawn", "overall_doe", "overall_fd")
for(j in 1:5){
    fawndoe_df[nrow(fawndoe_df) + 1, ] <- NA
}
indx_add <- which(is.na(fawndoe_df$year))
fawndoe_df$year[indx_add] <- fd_older_df$year
fawndoe_df$overall_doe[indx_add] <- fd_older_df$overall_doe
fawndoe_df$overall_fawn[indx_add] <- fd_older_df$overall_fawn
fawndoe_df$overall_fd[indx_add] <- fd_older_df$overall_fd
fawndoe_df <- fawndoe_df[order(fawndoe_df$year), ]
fawndoe_df <- fawndoe_df[fawndoe_df$year > 1993, ]

####################################################################################
###
### Loading data on occurance of earn-a-buck
###
#####################################################################################

df_eab <- read.csv(paste0(filepath,"eab_present.csv"))

####################################################################################
###
### Loading data SAK estimates of total population size
###
#####################################################################################

df_pop_estimate <-  read.csv(paste0(filepath, "Total_pop_size_UNIT.csv"))
names(df_pop_estimate) <- c("year","east","west","total")
pop_estimate_east <- df_pop_estimate$east[3]
pop_estimate_west <- df_pop_estimate$west[3]

