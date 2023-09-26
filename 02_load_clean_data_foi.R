###################################################################################################################################
# Load Data
###################################################################################################################################

cwd_df_1 <- suppressWarnings(read_excel("datafiles/1999-2021countydeer.xlsx",sheet=1))
cwd_df_2 <- suppressWarnings(read_excel("datafiles/1999-2021countydeer.xlsx",sheet=2))
cwd_df_3 <- suppressWarnings(read_excel("datafiles/1999-2021countydeer.xlsx",sheet=3))

cwd_df <- as.data.frame(rbind(cwd_df_1,cwd_df_2,cwd_df_3))
names(cwd_df) <- gsub(" ","_",names(cwd_df))
names(cwd_df) <- gsub("/","_",names(cwd_df))
names(cwd_df) <- tolower(names(cwd_df))

###################################################################################################################################
# Cleaning Data
###################################################################################################################################

###
### Restricting to counties we care about
###

counties <- c("Iowa","Dane","Grant")
cwd_df <- cwd_df[cwd_df$county %in% counties,]

###
### Cleaning up age classes
###
# "44291" = "4-5"
# "44450" = "9-11"
# "44355" = "6-8"
# "44656" = "4-5"
# "44720" =  "6-8"
# "44815" = "9-11"

cwd_df$age <- as.factor(cwd_df$age)
levels(cwd_df$age) <- c("1","12+","2","3","4-5","6-8","9-11","4-5","6-8","9-11","ADULT","0")
cwd_df <- cwd_df[cwd_df$age != "ADULT", ]
cwd_df$age <- as.factor(as.character(cwd_df$age))
cwd_df$kill_date <- as.Date(cwd_df$kill_date,origin="1899-12-30")

#removing surveillence data with no town, range, range direction, or section
cwd_df <- cwd_df[!is.na(cwd_df$dtrs),]
cwd_df <- cwd_df[!is.na(cwd_df$sect),]
cwd_df <- cwd_df[!is.na(cwd_df$town),]
cwd_df <- cwd_df[!is.na(cwd_df$range),]
cwd_df <- cwd_df[!is.na(cwd_df$range_dir),]

# cwd_df$trs <- paste0(cwd_df$town,"-",cwd_df$range,"-",cwd_df$sect)
# removing deer without a kill date
cwd_df <- cwd_df[order(cwd_df$kill_date),]
cwd_df <- cwd_df[!is.na(cwd_df$kill_date),]

#removing surveillence data with no town, range, range direction, or section
cwd_df$trs <- paste0(cwd_df$town,"-",cwd_df$range,"-",cwd_df$sect)
cwd_df <- cwd_df[order(cwd_df$kill_date),]
cwd_df <- cwd_df[!is.na(cwd_df$kill_date),]


# nrow(cwd_df[cwd_df$negative=="N" & cwd_df$positive=="N" & cwd_df$has_results=="Y",])

source("cleanData_foi.R")

cwd_df <- cleanData(cwd_df)

#######################################
###
### Restricing to the study area
###
#######################################

study_df <- sf::st_read("datafiles/study_df.shp")
study_df <- study_df %>% st_transform(crs = 3071)

#creating sections that account for range_direction first
study_df$dsection <- paste0(study_df$dir,"-",study_df$sectionid)
cwd_df$dsection <- do.call(paste, c(cwd_df[c("range_dir",
                                             "range",
                                             "town",
                                             "sect")], sep = "-"))

cwd_df <- cwd_df[cwd_df$dsection %in% study_df$dsection, ]

cwd_df$year <- lubridate::year(cwd_df$kill_date)

#setting all deer killed in January 2022 to be part of study year 2021
cwd_df$year[cwd_df$year == 2022] <- 2021

#setting up east/west study areas
cwd_df$ew <- rep(NA,nrow(cwd_df))
for(i in 1:nrow(study_df)) {
        j <- which(cwd_df$dsection %in% study_df$dsection[i])
            cwd_df$ew[j] <- study_df$ew[i]
        }
#table(cwd_df$ew)
# head(cwd_df$ew)
cwd_df$ew <- as.numeric(as.factor(cwd_df$ew))#east ==1, west == 2

######################################################
# Aggregate Oldest Age Class of Males
######################################################

male6 <- cwd_df %>% filter(age_num == 6 & sex == 1)
cwd_df$agedays[cwd_df$age_num == 9 & cwd_df$sex == 1] <- max(male6$agedays)

#######################################################
### setting these to the maximum number of 
### possible days, weeks and months
### that a 6.5-8.5 year old could live until
### aging into the next age class
#######################################################

# hist(cwd_df$ageweeks)

# cwd_df$ageweeks[cwd_df$age_num == 9 & cwd_df$sex == 1] <- 467#max(male6$ageweeks)
# cwd_df$agemonths[cwd_df$age_num == 9 & cwd_df$sex == 1] <- 107#max(male6$agemonths)
cwd_df$ageweeks[cwd_df$age_num == 9 & cwd_df$sex == 1] <- 338#max(male6$ageweeks)
cwd_df$agemonths[cwd_df$age_num == 9 & cwd_df$sex == 1] <- 78#max(male6$agemonths)
cwd_df$agemonths[cwd_df$agemonths > 78 & cwd_df$sex == 1] <- 78#max(male6$agemonths)


cwd_df$age_num[cwd_df$age_num == 9 & cwd_df$sex == 1] <- 6
ageclass <- as.numeric(levels(as.factor(cwd_df$age_num)))

# length(c(rep(c(0.5,1.5,2.5,3.5), each = 12),
#     rep(4.5,24),6.5))
#
# 73 weeks is the maximum need to integrate over for individuals that are 6.5+ for males in months
# 72+24+11
# rep(1,2,3,4),5,6
# length(c(rep(c(0.5,1.5,2.5,3.5), each = 52), rep(4.5,104),6.5))
# length(c(rep(c(0.5,1.5,2.5,3.5), each = 52), rep(4.5,104),rep(6.5,103)))
# 6*52+26
# 6*12+6
# table(cwd_df$agemonths)
