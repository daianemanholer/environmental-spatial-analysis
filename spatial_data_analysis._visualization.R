#Author: Daiane Domingos Manholer

#Repository purpose: Selected R script fragments developed as part of the data analysis for my master's 
#research on the behavioral responses of cetaceans to whale-watching activities in Brazil

#Application context: This repository is shared exclusively in support of my application to the 2027 IAI HARMONIZE Workshop,
#"Harmonizing and Modeling Climate, Environment and Health Data to Support Public Health Decision-Making." 
#NOTE: Only selected script fragments are shared in this repository, as the research is part of a manuscript currently under
#preparation for publication. The complete data, scripts and analyses are therefore not publicly available at this stage

library(dplyr)
library(tidyverse)
library(sf)
library(readxl)
library(lubridate)
library(mapview)
library(REdaS)
library(geosphere)
library(nngeo)
library(GEOmap)
library(openxlsx)
library(hms)
library(ggplot2)
library(dunn.test)
library(mgcv)
library(FSA)
library(MASS)
library(ggpubr)
library(car)
library(betareg)
library(broom)
library(truncnorm)
library(patchwork)

#######################################1########################################
####################################Data input##################################

dbTheod2022 <- readxl::read_xlsx("./data_input/theod_behavior_2022.xlsx") %>%
  dplyr::select(-c(8)) #removing columns that will not be used

dbTheod2023 <- readxl::read_xlsx("./data_input/theod_behavior_2023_corrected.xlsx") %>%
  dplyr::select(-c(8)) #removing columns that will not be used

dbDrone2023_2024 <- readxl::read_xlsx("./data_input/drone_behavior_2023_2024.xlsx")

df2022 <- readxl::read_xlsx("./data_input/location_theod_2022.xlsx") %>%
  dplyr::select(-c(2, 9, 25, 30:36, 38)) #removing columns without GPS location information and other information that will not be used

df2023 <- readxl::read_xlsx("./data_input/location_theod_2023.xlsx") %>%
  dplyr::select(-c(2, 9, 25, 30:36, 38)) #removing columns without GPS location information and other information that will not be used

dlocDrone2023_2024 <- readxl::read_xlsx("./data_input/location_drone_2023_2024.xlsx")


#######################################3########################################
################################Filtrate data###################################

#Exclude focal follows with < 3 locational points (not able to properly evaluate linearity and,
#swim speed, but will keep being analyzed for states and events

#loc theodolite
df2022 <- df2022 %>%
  dplyr::group_by(`ID MPF`, Group) %>%
  dplyr::filter(dplyr::n() >= 3) %>% 
  dplyr::ungroup()

df2023 <- df2023 %>%
  dplyr::group_by(`ID MPF`, Group) %>%
  dplyr::filter(dplyr::n() >= 3) %>% 
  dplyr::ungroup()

#loc drone
dlocDrone2023_2024 <- dlocDrone2023_2024 %>%
  dplyr::group_by(`ID MPF`, Group) %>%
  dplyr::filter(dplyr::n() >= 3) %>% 
  dplyr::ungroup()


#######################################4########################################
################################Obtain metrics##################################

#Swim speed, linearity and reorientation rate for drone data
dlocDroneSpatial <- dlocDrone2023_2024 %>%
  dplyr::mutate(Destlat1 = Lat,
                Destlon1 = Long) %>%
  sf::st_as_sf(coords = c("Destlon1","Destlat1"), crs = 4326) %>%
  sf::st_transform(31982)

mapview::mapview(dlocDroneSpatial)

#Swim speed
DroneSpeed <- 
  dlocDroneSpatial %>%
  dplyr::mutate(Time = lubridate::as_datetime(Time)) %>%
  dplyr::group_by(`ID MPF`, Group) %>%
  dplyr::mutate(
    distFix = as.numeric(st_distance(geometry, lead(geometry),  #calculate distance between consecutive points inside each group
                                     by_element = TRUE))/1000,
    timeDiff = as.numeric(difftime(lead(Time), Time, units = "hours")),
    speed = round(distFix/timeDiff,digits = 1)) %>%
  ungroup()

DroneSpeedMean <- DroneSpeed %>%
  dplyr::filter(!is.na(speed)) %>%
  dplyr::group_by(`ID MPF`, Group) %>%
  dplyr::summarise(meanSpeed = mean(speed), .groups = "drop")


#Linearity 
DroneLinearity <- 
  dlocDroneSpatial %>%
  dplyr::group_by(`ID MPF`, Group) %>%
  dplyr::mutate(
    distFix = as.numeric(st_distance(geometry, lead(geometry), 
                                     by_element = TRUE))/1000 # Consecutive distances
  ) %>%
  dplyr::summarise(
    distSum = sum(distFix, na.rm = TRUE),  # Total path distance
    distInitialFinal = as.numeric(st_distance(first(geometry), 
                                              last(geometry)))/1000, # Start-end distance
    Linearity = distInitialFinal / distSum, # Linearity ratio
    .groups = "drop"
  )


#Reorientation rate

#Filtrating
dlocDroneRemoveRR <- dlocDroneSpatial %>%  #crating data frame to see how many they are
  group_by(`ID MPF`, Group) %>%
  filter(n() <= 6) %>%
  ungroup()

dlocDroneSpatialRR <- dlocDroneSpatial %>%
  group_by(`ID MPF`, Group) %>%
  filter(n() > 6) %>%  #removing focal follow with < 7 location points
  ungroup()

#add total time to behavior data
dlocDroneTotalTime <- dbDrone2023_2024 %>%
  dplyr::distinct(`ID MPF`, Group, Duration_min)

dlocDroneTotalTime <- dlocDroneTotalTime %>%
  dplyr::rename(TotalTime = Duration_min)

dlocDroneSpatialRR <- dlocDroneSpatialRR %>%
  dplyr::left_join(dlocDroneTotalTime, by = c("ID MPF", "Group"))

#Calculating reorientation rate 
dlocDroneSpatialRR <- dlocDroneSpatialRR %>%
  dplyr::mutate(
    Lat = as.numeric(Lat),
    Long = as.numeric(Long)
  ) 

options(scipen = 999)

DroneReorientation <- dlocDroneSpatialRR %>%
  dplyr::mutate(Time = lubridate::as_datetime(Time)) %>%
  dplyr::group_by(`ID MPF`, Group) %>%
  dplyr::mutate(
    bearing = c(NA, sapply(1:(n() - 1), function(i) {     #Bearing calculation
      bear <- geosphere::bearing(c(Long[i], Lat[i]),
                                 c(Long[i + 1], Lat[i + 1]))
      if (bear < 0) bear + 360 else bear     #Transform negative into positive
    })),
    diffAngle = c(NA, diff(bearing)),     #Difference between consecutive bearings 
    diffAngle = abs(diffAngle),
    diffAngleCorr = ifelse(diffAngle > 180, 360 - diffAngle, diffAngle), #Set range from 0 to 180, not considering direction
    reorientationRate = sum(diffAngleCorr, na.rm = TRUE) / unique(TotalTime) #Reorientation rate
  ) %>%
  dplyr::ungroup()

#Remove bearings = 0 that could overestimate reorientation rate
DroneReorientation <- DroneReorientation %>%
  dplyr::filter(bearing != 0 | is.na(bearing))

#Recalculate reorientation rate
DroneReorientation <- DroneReorientation %>%
  dplyr::mutate(Time = lubridate::as_datetime(Time)) %>%
  dplyr::group_by(`ID MPF`, Group) %>%
  dplyr::mutate(
    bearing = c(NA, sapply(1:(n() - 1), function(i) {     #Bearing calculation
      bear <- geosphere::bearing(c(Long[i], Lat[i]),
                                 c(Long[i + 1], Lat[i + 1]))
      if (bear < 0) bear + 360 else bear     #Transform negative into positive
    })),
    diffAngle = c(NA, diff(bearing)),     #Difference between consecutive bearings 
    diffAngle = abs(diffAngle),
    diffAngleCorr = ifelse(diffAngle > 180, 360 - diffAngle, diffAngle), #Set range from 0 to 180, not considering direction
    reorientationRate = sum(diffAngleCorr, na.rm = TRUE) / unique(TotalTime) #Reorientation rate
  ) %>%
  dplyr::ungroup()

#Box plot to see how data looks like

#DroneLinearity
ggplot(
  DroneLinearity %>%
    filter(str_detect(Group, "C|T")) %>%
    mutate(Scenario = case_when(
      str_detect(Group, "C") ~ "Control",
      str_detect(Group, "T") ~ "Treatment"
    )),
  aes(x = Scenario, y = Linearity, fill = Scenario)
) +
  geom_boxplot() +
  labs(
    title = "Drone - Linearity",
    x = "Scenario",
    y = "Linearity"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("Control" = "steelblue", "Treatment" = "salmon"))

#DroneReorientation
ggplot(
  DroneReorientation %>%
    filter(str_detect(Group, "C|T")) %>%
    mutate(Scenario = case_when(
      str_detect(Group, "C") ~ "Control",
      str_detect(Group, "T") ~ "Treatment"
    )) %>%
    distinct(`ID MPF`, Group, Scenario, reorientationRate),  # remove duplicatas
  aes(x = Scenario, y = reorientationRate, fill = Scenario)
) +
  geom_boxplot() +
  labs(
    title = "Drone - Reorientation",
    x = "Scenario",
    y = "Reorientation rate"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("Control" = "steelblue", "Treatment" = "salmon"))


#DroneSpeedMean
ggplot(
  DroneSpeedMean %>%
    filter(str_detect(Group, "C|T")) %>%
    mutate(Scenario = case_when(
      str_detect(Group, "C") ~ "Control",
      str_detect(Group, "T") ~ "Treatment"
    )),
  aes(x = Scenario, y = meanSpeed, fill = Scenario)
) +
  geom_boxplot() +
  labs(
    title = "Drone - Speed Mean",
    x = "Scenario",
    y = "Speed Mean"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("Control" = "steelblue", "Treatment" = "salmon"))
