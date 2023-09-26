
# rm(list = ls())
# setwd("~/Documents/aah/aah_disease_ew")

# library(ggplot2)
# library(nimble)
# library(Matrix)
# library(coda)
# library(lattice)
# library(splines)
# library(Hmisc)
# library(lubridate)
# library(readxl)
# library(gridExtra)
# library(xtable)
# library(beepr)
# library(RColorBrewer)
# library(reshape2)
# library(viridis)
# library(ggridges)
# library(doParallel)
# library(dplyr)
library(tidyr)
library(data.table)
library(abind)
library(sf)
library(MetBrewer)
library(terra)
library(tidyterra)


######################################################
### Loading county-based polygon shapefile and 
### clip to counties of study area
######################################################

county <- sf::st_read("datafiles/county_bnds.shp")
county <- st_transform(county,3701)
county <- county[county$CTY_NAME %in% c("Dane","Iowa","Grant"),]
grant <- county[county$CTY_NAME %in% c("Grant"),]
iowa <- county[county$CTY_NAME %in% c("Iowa"),]
dane <- county[county$CTY_NAME %in% c("Dane"),]
county_tres <- st_union(county)
county_tres_bound <- sf::st_union(county_tres)

# county <- terra::project(county,"epsg:3071")

######################################################
### Loading study area section-based polygon shapefile 
######################################################

# study_df <- terra::vect("~/Documents/Data/Study_Area/study_df.shp")
# plot(study_df)
# study_df <- terra::project(study_df,"epsg:4326")

study_df <- sf::st_read("datafiles/study_df.shp")
study_df <- st_transform(study_df,3701)

study_df_east_sect <- study_df[study_df$ew == "east",]
study_df_west_sect <- study_df[study_df$ew == "west",]
study_area_bound <- sf::st_union(study_df)
study_east_bound <- sf::st_union(study_df_east_sect)
study_west_bound <- sf::st_union(study_df_west_sect)

# st_write(study_area_bound,"~/Documents/Data/Study_Area/study_area_bound.shp")
# st_write(study_west_bound,"~/Documents/Data/Study_Area/study_west_bound.shp")
# st_write(study_east_bound,"~/Documents/Data/Study_Area/study_east_bound.shp")


grant_w_study <- st_intersection(grant,study_area_bound)
iowa_w_study <- st_intersection(iowa,study_west_bound)
iowa_e_study <- st_intersection(iowa,study_east_bound)
dane_e_study <- st_intersection(dane,study_area_bound)

# st_write(grant_w_study,"~/Documents/Data/Study_Area/grant_w_study.shp")
# st_write(iowa_w_study,"~/Documents/Data/Study_Area/iowa_w_study.shp")
# st_write(iowa_e_study,"~/Documents/Data/Study_Area/iowa_e_study.shp")
# st_write(dane_e_study,"~/Documents/Data/Study_Area/dane_e_study.shp")

#convert from sf object to terra object for extraction
grant_w_study <- vect(grant_w_study)
iowa_w_study <- vect(iowa_w_study)
iowa_e_study <- vect(iowa_e_study)
dane_e_study <- vect(dane_e_study)
county_tres_bound <- vect(county_tres_bound)
grant <- vect(grant)
iowa <- vect(iowa)
dane <- vect(dane)
# plot(study_df)
# plot(grant_w_study)
# plot(iowa_w_study)
# plot(iowa_e_study)
# plot(dane_e_study)
# plot(county_tres_bound)


######################################################
### Loading raster of deer habitat
######################################################

deer_habitat <- terra::rast("/home/aketz/Documents/Data/230523_harvest_data/finalgrid_2017.tif")
# deer_habitat <- terra::rast("/home/aketz/Documents/Data/230523_harvest_data/deer_habitat.tif")
deer_habitat <- project(deer_habitat,"EPSG:3701")

#crop the WI deer habitat to the three relevant counties
hab <- crop(deer_habitat,county_tres_bound)
# plot(hab)
# head(hab)

#setting all na values to 0, which are not deer habitat
hab[is.na(hab)] <- 0
unique(hab$SQ_MI)

######################################################
### Loading county-based polygon shapefile and 
### clip to counties of study area
######################################################

dmu86 <- sf::st_read("datafiles/dmu_1986.shp")
dmu86 <- st_transform(dmu86,3701)
dmu86 <- dmu86[dmu86$UNIT_ID %in% c("70A","70C","70D","73C","73"),]
dmu86_bound <- st_union(dmu86)

dmu99 <- sf::st_read("datafiles/dmu_1999.shp")
dmu99 <- st_transform(dmu99,3701)
dmu99 <- dmu99[dmu99$UNIT_ID %in% c("70A","70C","70D","73E"),]
dmu99_bound <- st_union(dmu99)

dmu13 <- sf::st_read("datafiles/dmu_2013.shp")
dmu13 <- st_transform(dmu13,3701)
dmu13 <- dmu13[dmu13$UNIT_ID %in% c("70A-CWD","70C-CWD","70D-CWD","73E-CWD"),]
dmu13_bound <- st_union(dmu13)


# st_write(dmu86_bound,"~/Documents/Data/Study_Area/dmu86_bound.shp")
# st_write(dmu99_bound,"~/Documents/Data/Study_Area/dmu99_bound.shp")
# st_write(dmu13_bound,"~/Documents/Data/Study_Area/dmu13_bound.shp")

# plot(dmu86)
# table(dmu86$UNIT_ID)

###########################################################################
### creating dmu based intersections with the study area for each year
### 1986
###########################################################################

dmu86_73 <- dmu86[dmu86$UNIT_ID == "73",]
dmu86_73c <- dmu86[dmu86$UNIT_ID == "73C",]
dmu86_70d <- dmu86[dmu86$UNIT_ID == "70D",]
dmu86_70c <- dmu86[dmu86$UNIT_ID == "70C",]
dmu86_70a <- dmu86[dmu86$UNIT_ID == "70A",]

dmu86_73_w_study <- st_intersection(dmu86_73,study_area_bound)
dmu86_73c_w_study <- st_intersection(dmu86_73c,study_area_bound)
dmu86_70d_e_study <- st_intersection(dmu86_70d,study_area_bound)
dmu86_70c_e_study <- st_intersection(dmu86_70c,study_area_bound)
dmu86_70a_e_study <- st_intersection(dmu86_70a,study_area_bound)

### 1999 DMU
dmu99_73e <- dmu99[dmu99$UNIT_ID == "73E",]
dmu99_70d <- dmu99[dmu99$UNIT_ID == "70D",]
dmu99_70c <- dmu99[dmu99$UNIT_ID == "70C",]
dmu99_70a <- dmu99[dmu99$UNIT_ID == "70A",]

dmu99_73e_w_study <- st_intersection(dmu99_73e,study_west_bound)
dmu99_70d_e_study <- st_intersection(dmu99_70d,study_east_bound)
dmu99_70c_e_study <- st_intersection(dmu99_70c,study_area_bound)
dmu99_70a_e_study <- st_intersection(dmu99_70a,study_area_bound)


####2002-2013 DMU boundaries
dmu13_73e_cwd <- dmu13[dmu13$UNIT_ID == "73E-CWD",]
dmu13_70d_cwd <- dmu13[dmu13$UNIT_ID == "70D-CWD",]
dmu13_70c_cwd <- dmu13[dmu13$UNIT_ID == "70C-CWD",]
dmu13_70a_cwd <- dmu13[dmu13$UNIT_ID == "70A-CWD",]

dmu13_73e_w_study <- st_intersection(dmu13_73e_cwd,study_west_bound)
dmu13_70d_e_study <- st_intersection(dmu13_70d_cwd,study_east_bound)
dmu13_70c_e_study <- st_intersection(dmu13_70c_cwd,study_area_bound)
dmu13_70a_e_study <- st_intersection(dmu13_70a_cwd,study_area_bound)

############################
###
### converting dmus to vect 
### objects for extraction
###
############################

dmu86_73 <- vect(dmu86_73)
dmu86_73c <- vect(dmu86_73c)
dmu86_70d <- vect(dmu86_70d)
dmu86_70c <- vect(dmu86_70c)
dmu86_70a <- vect(dmu86_70a)

dmu86_73_w_study <- vect(dmu86_73_w_study)
dmu86_73c_w_study <- vect(dmu86_73c_w_study)
dmu86_70d_e_study <- vect(dmu86_70d_e_study)
dmu86_70c_e_study <- vect(dmu86_70c_e_study)
dmu86_70a_e_study <- vect(dmu86_70a_e_study)

### 1999 DMU
dmu99_73e <- vect(dmu99_73e)
dmu99_70d <- vect(dmu99_70d)
dmu99_70c <- vect(dmu99_70c)
dmu99_70a <- vect(dmu99_70a)

dmu99_73e_w_study <- vect(dmu99_73e_w_study)
dmu99_70d_e_study <- vect(dmu99_70d_e_study)
dmu99_70c_e_study <- vect(dmu99_70c_e_study)
dmu99_70a_e_study <- vect(dmu99_70a_e_study)


####2002-2013 DMU boundaries
dmu13_73e_cwd <- vect(dmu13_73e_cwd)
dmu13_70d_cwd <- vect(dmu13_70d_cwd)
dmu13_70c_cwd <- vect(dmu13_70c_cwd)
dmu13_70a_cwd <- vect(dmu13_70a_cwd)

dmu13_73e_w_study <- vect(dmu13_73e_w_study)
dmu13_70d_e_study <- vect(dmu13_70d_e_study)
dmu13_70c_e_study <- vect(dmu13_70c_e_study)
dmu13_70a_e_study <- vect(dmu13_70a_e_study)



#############################################################################
### extract the sum of the number of raster pixels 
### that are deer habitat for each dmu/count/studyarea intersection part
#############################################################################

grant_all_hab <- extract(hab, grant)
iowa_all_hab <- extract(hab, iowa)
dane_all_hab <- extract(hab, dane)

grant_w_study_hab <- extract(hab, grant_w_study)
iowa_w_study_hab <- extract(hab, iowa_w_study)
iowa_e_study_hab <- extract(hab, iowa_e_study)
dane_e_study_hab <- extract(hab, dane_e_study)

dmu86_73_all_hab <- extract(hab, dmu86_73)
dmu86_73c_all_hab <- extract(hab, dmu86_73c)
dmu86_70d_all_hab <- extract(hab, dmu86_70d)
dmu86_70c_all_hab <- extract(hab, dmu86_70c)
dmu86_70a_all_hab <- extract(hab, dmu86_70a)

dmu86_73_w_study_hab <- extract(hab, dmu86_73_w_study)
dmu86_73c_w_study_hab <- extract(hab, dmu86_73c_w_study)
dmu86_70d_e_study_hab <- extract(hab, dmu86_70d_e_study)
dmu86_70c_e_study_hab <- extract(hab, dmu86_70c_e_study)
dmu86_70a_e_study_hab <- extract(hab, dmu86_70a_e_study)

dmu99_73e_all_hab <- extract(hab, dmu99_73e)
dmu99_70d_all_hab <- extract(hab, dmu99_70d)
dmu99_70c_all_hab <- extract(hab, dmu99_70c)
dmu99_70a_all_hab <- extract(hab, dmu99_70a)

dmu99_73e_w_study_hab <- extract(hab, dmu99_73e_w_study)
dmu99_70d_e_study_hab <- extract(hab, dmu99_70d_e_study)
dmu99_70c_e_study_hab <- extract(hab, dmu99_70c_e_study)
dmu99_70a_e_study_hab <- extract(hab, dmu99_70a_e_study)

dmu13_73e_cwd_all_hab <- extract(hab, dmu13_73e_cwd)
dmu13_70d_cwd_all_hab <- extract(hab, dmu13_70d_cwd)
dmu13_70c_cwd_all_hab <- extract(hab, dmu13_70c_cwd)
dmu13_70a_cwd_all_hab <- extract(hab, dmu13_70a_cwd)

dmu13_73e_w_study_hab <- extract(hab, dmu13_73e_w_study)
dmu13_70d_e_study_hab <- extract(hab, dmu13_70d_e_study)
dmu13_70c_e_study_hab <- extract(hab, dmu13_70c_e_study)
dmu13_70a_e_study_hab <- extract(hab, dmu13_70a_e_study)




######################################################
### Makin' relevant maps
######################################################

# study_area_old_dmus <- ggplot()+geom_sf(data=dmu86)+geom_sf(data=study_area_bound,color="darkred",alpha=.5)+theme_bw()
# ggsave("figures/study_area_old_dmus.png")

######################################################
### Crop raster by county
######################################################

county_deer_habitat <- terra::crop(deer_habitat,county_tres)
study_deer_habitat <- terra::crop(deer_habitat,study_df)

png("figures/county_deer_habitat_rast.png")
plot(county_deer_habitat,main="deer habitat raster counties")
dev.off()
png("figures/study_deer_habitat_rast.png")
plot(study_deer_habitat,main="deer habitat raster study_overall")
dev.off()
# values(county_deer_habitat)[!is.finite(values(county_deer_habitat))] <- 0
# values(study_deer_habitat)[!is.finite(values(study_deer_habitat))] <- 0
study_deer_hab <- mask(study_deer_habitat,vect(study_area_bound))
png("figures/study_deer_habitat_rast.png")
plot(study_deer_hab)
dev.off()



# iowa_all_hab$SQ_MI <- as.factor(iowa_all_hab$SQ_MI)
# levels(iowa_all_hab$SQ_MI) <- c(0,1)
# iowa_all_hab$SQ_MI <- as.numeric(as.character(iowa_all_hab$SQ_MI))
# iowa_all_hab_rel <- table(iowa_all_hab$SQ_MI)
# iowa_all_hab_rel[1]/(iowa_all_hab_rel[1]+iowa_all_hab_rel[2])

# grant_all_hab$SQ_MI <- as.factor(grant_all_hab$SQ_MI)
# levels(grant_all_hab$SQ_MI) <- c(0,1)
# grant_all_hab$SQ_MI <- as.numeric(as.character(grant_all_hab$SQ_MI))
# grant_all_hab_rel <- table(grant_all_hab$SQ_MI)
# grant_all_hab_rel[1]/(grant_all_hab_rel[1]+grant_all_hab_rel[2])

# dane_all_hab$SQ_MI <- as.factor(dane_all_hab$SQ_MI)
# levels(dane_all_hab$SQ_MI) <- c(0,1)
# dane_all_hab$SQ_MI <- as.numeric(as.character(dane_all_hab$SQ_MI))
# dane_all_hab_rel <- table(dane_all_hab$SQ_MI)
# dane_all_hab_rel[1]/(dane_all_hab_rel[1]+dane_all_hab_rel[2])


# iowa_w_study_hab$SQ_MI <- as.factor(iowa_w_study_hab$SQ_MI)
# levels(iowa_w_study_hab$SQ_MI) <- c(0,1)
# iowa_w_study_hab$SQ_MI <- as.numeric(as.character(iowa_w_study_hab$SQ_MI))
# iowa_w_study_hab_rel <- table(iowa_w_study_hab$SQ_MI)
# iowa_w_study_hab_rel[1]/(iowa_w_study_hab_rel[1] + iowa_w_study_hab_rel[2])

# iowa_e_study_hab$SQ_MI <- as.factor(iowa_e_study_hab$SQ_MI)
# levels(iowa_e_study_hab$SQ_MI) <- c(0,1)
# iowa_e_study_hab$SQ_MI <- as.numeric(as.character(iowa_e_study_hab$SQ_MI))
# iowa_e_study_hab_rel <- table(iowa_e_study_hab$SQ_MI)
# iowa_e_study_hab_rel[1]/(iowa_e_study_hab_rel[1] + iowa_e_study_hab_rel[2])


calc_sum_deer_habitat <- function(x){
    x$SQ_MI <- as.factor(x$SQ_MI)
    levels(x$SQ_MI) <- c(0,1)
    x$SQ_MI <- as.numeric(as.character(x$SQ_MI))
    x_hab_rel <- table(x$SQ_MI)
    x_hab_rel[1]
}


grant_all_habsum <- calc_sum_deer_habitat(grant_all_hab)
iowa_all_habsum <- calc_sum_deer_habitat(iowa_all_hab)
dane_all_habsum <- calc_sum_deer_habitat(dane_all_hab)

grant_w_study_habsum <- calc_sum_deer_habitat(grant_w_study_hab)
iowa_w_study_habsum <- calc_sum_deer_habitat(iowa_w_study_hab)
iowa_e_study_habsum <- calc_sum_deer_habitat(iowa_e_study_hab)
dane_e_study_habsum <- calc_sum_deer_habitat(dane_e_study_hab)

dmu86_73_all_habsum <- calc_sum_deer_habitat(dmu86_73_all_hab)
dmu86_73c_all_habsum <- calc_sum_deer_habitat(dmu86_73c_all_hab)
dmu86_70d_all_habsum <- calc_sum_deer_habitat(dmu86_70d_all_hab)
dmu86_70c_all_habsum <- calc_sum_deer_habitat(dmu86_70c_all_hab)
dmu86_70a_all_habsum <- calc_sum_deer_habitat(dmu86_70a_all_hab)

dmu86_73_w_study_habsum <- calc_sum_deer_habitat(dmu86_73_w_study_hab)
dmu86_73c_w_study_habsum <- calc_sum_deer_habitat(dmu86_73c_w_study_hab)
dmu86_70d_e_study_habsum <- calc_sum_deer_habitat(dmu86_70d_e_study_hab)
dmu86_70c_e_study_habsum <- calc_sum_deer_habitat(dmu86_70c_e_study_hab)
dmu86_70a_e_study_habsum <- calc_sum_deer_habitat(dmu86_70a_e_study_hab)

dmu99_73e_all_habsum <- calc_sum_deer_habitat(dmu99_73e_all_hab)
dmu99_70d_all_habsum <- calc_sum_deer_habitat(dmu99_70d_all_hab)
dmu99_70c_all_habsum <- calc_sum_deer_habitat(dmu99_70c_all_hab)
dmu99_70a_all_habsum <- calc_sum_deer_habitat(dmu99_70a_all_hab)

dmu99_73e_w_study_habsum <- calc_sum_deer_habitat(dmu99_73e_w_study_hab)
dmu99_70d_e_study_habsum <- calc_sum_deer_habitat(dmu99_70d_e_study_hab)
dmu99_70c_e_study_habsum <- calc_sum_deer_habitat(dmu99_70c_e_study_hab)
dmu99_70a_e_study_habsum <- calc_sum_deer_habitat(dmu99_70a_e_study_hab)

dmu13_73e_cwd_all_habsum <- calc_sum_deer_habitat(dmu13_73e_cwd_all_hab)
dmu13_70d_cwd_all_habsum <- calc_sum_deer_habitat(dmu13_70d_cwd_all_hab)
dmu13_70c_cwd_all_habsum <- calc_sum_deer_habitat(dmu13_70c_cwd_all_hab)
dmu13_70a_cwd_all_habsum <- calc_sum_deer_habitat(dmu13_70a_cwd_all_hab)

dmu13_73e_w_study_habsum <- calc_sum_deer_habitat(dmu13_73e_w_study_hab)
dmu13_70d_e_study_habsum <- calc_sum_deer_habitat(dmu13_70d_e_study_hab)
dmu13_70c_e_study_habsum <- calc_sum_deer_habitat(dmu13_70c_e_study_hab)
dmu13_70a_e_study_habsum <- calc_sum_deer_habitat(dmu13_70a_e_study_hab)




################################################
### calculating the deer habitat correction
#################################################

grant_correct <- grant_w_study_habsum/grant_all_habsum
iowa_w_correct <- iowa_w_study_habsum/iowa_all_habsum
iowa_e_correct <- iowa_e_study_habsum/iowa_all_habsum
dane_e_correct <- dane_e_study_habsum/dane_all_habsum
dmu86_73_w_correct <- dmu86_73_w_study_habsum/dmu86_73_all_habsum
dmu86_73c_w_correct <- dmu86_73c_w_study_habsum/dmu86_73c_all_habsum
dmu86_70d_e_correct <- dmu86_70d_e_study_habsum/dmu86_70d_all_habsum
dmu86_70c_e_correct <- dmu86_70c_e_study_habsum/dmu86_70c_all_habsum
dmu86_70a_e_correct <- dmu86_70a_e_study_habsum/dmu86_70a_all_habsum
dmu99_73e_w_correct <- dmu99_73e_w_study_habsum/dmu99_73e_all_habsum
dmu99_70d_e_correct <- dmu99_70d_e_study_habsum/dmu99_70d_all_habsum
dmu99_70c_e_correct <- dmu99_70c_e_study_habsum/dmu99_70c_all_habsum
dmu99_70a_e_correct <- dmu99_70a_e_study_habsum/dmu99_70a_all_habsum
dmu13_73e_w_correct <- dmu13_73e_w_study_habsum/dmu13_73e_cwd_all_habsum
dmu13_70d_e_correct <- dmu13_70d_e_study_habsum/dmu13_70d_cwd_all_habsum
dmu13_70c_e_correct <- dmu13_70c_e_study_habsum/dmu13_70c_cwd_all_habsum
dmu13_70a_e_correct <- dmu13_70a_e_study_habsum/dmu13_70a_cwd_all_habsum



################################################
### saving the deer habitat correction
#################################################

# filepath2 <- "~/Documents/Data/Habitat_Calibration/"

# save(grant_correct, file = paste0(filepath2,"grant_correct.Rdata"))
# save(iowa_w_correct, file = paste0(filepath2,"iowa_w_correct.Rdata"))
# save(iowa_e_correct, file = paste0(filepath2,"iowa_e_correct.Rdata"))
# save(dane_e_correct, file = paste0(filepath2,"dane_e_correct.Rdata"))
# save(dmu86_73_w_correct, file = paste0(filepath2,"dmu86_73_w_correct.Rdata"))
# save(dmu86_73c_w_correct, file = paste0(filepath2,"dmu86_73c_w_correct.Rdata"))
# save(dmu86_70d_e_correct, file = paste0(filepath2,"dmu86_70d_e_correct.Rdata"))
# save(dmu86_70c_e_correct, file = paste0(filepath2,"dmu86_70c_e_correct.Rdata"))
# save(dmu86_70a_e_correct, file = paste0(filepath2,"dmu86_70a_e_correct.Rdata"))
# save(dmu99_73e_w_correct, file = paste0(filepath2,"dmu99_73e_w_correct.Rdata"))
# save(dmu99_70d_e_correct, file = paste0(filepath2,"dmu99_70d_e_correct.Rdata"))
# save(dmu99_70c_e_correct, file = paste0(filepath2,"dmu99_70c_e_correct.Rdata"))
# save(dmu99_70a_e_correct, file = paste0(filepath2,"dmu99_70a_e_correct.Rdata"))
# save(dmu13_73e_w_correct, file = paste0(filepath2,"dmu13_73e_w_correct.Rdata"))
# save(dmu13_70d_e_correct, file = paste0(filepath2,"dmu13_70d_e_correct.Rdata"))
# save(dmu13_70c_e_correct, file = paste0(filepath2,"dmu13_70c_e_correct.Rdata"))
# save(dmu13_70a_e_correct, file = paste0(filepath2,"dmu13_70a_e_correct.Rdata"))


