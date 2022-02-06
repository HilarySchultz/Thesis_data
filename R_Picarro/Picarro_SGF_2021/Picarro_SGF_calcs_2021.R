library(dplyr) 
library(readr) 
library(ggplot2) 
library(tidyverse)
library(lubridate)
options(scipen = 999)

## GET SGF FILES ##
filelist <-list.files(path = getwd(), pattern="*.csv", full.names = TRUE)
SGF_files <- list()
for(f in 1:length(filelist)) {
  SGF_files[[f]] <- read.csv(filelist[f], header = T, sep = ",", as.is = T)
}

# Creating a time vector for LR #

time_vec <- rep(c(5,35,65))

SGF_df <- bind_rows(SGF_files) %>%
  arrange(DATE) %>%
  mutate(Time = rep(time_vec, 287)) %>%
  rename(Reach = Design, Date = DATE, Chamber = chamber) %>%
  select(-peak) 
 
dates <- distinct(SGF_df, Date)
sites <- distinct(SGF_df, Site)
chambers <- distinct(SGF_df, Chamber)

## GET SLOPES FOR GAS CONCENTRATIONS ##
slope_df <- data.frame(Date = character(), Site = character(), Chamber = integer(), CO2_slope = numeric(), CH4_slope = numeric())
for(date in dates$Date){
  for(site in sites$Site) {
    for(chamber in chambers$Chamber){
      data <- SGF_df %>% filter(Site %in% site) %>% filter(Date %in% date) %>% filter(Chamber %in% chamber)
      if (length(data$X > 1)) {
        CO2_slope <- summary(lm(avg_CO2~Time, data))$coefficients[2,1]
        CH4_slope <- summary(lm(avg_CH4~Time, data))$coefficients[2,1]
        R2_CO2 <- summary(lm(avg_CO2~Time, data))$r.squared
        R2_CH4 <- summary(lm(avg_CH4~Time, data))$r.squared
        data_to_add <- data.frame(Date = date, Site = site, Chamber = chamber, 
                                  CO2_slope = CO2_slope, CH4_slope = CH4_slope, 
                                  R2_CO2 = R2_CO2, R2_CH4 = R2_CH4)
        slope_df <- rbind(slope_df, data_to_add)
      }
    }
  }
}
view(slope_df) 

slope_df <- slope_df %>%
  mutate(Reach = if_else(Chamber <= 6, "BDA", "REF"))

hist(slope_df$R2_CO2)
plot(slope_df$R2_CO2, slope_df$R2_CH4)

# Separting methane and carbon dioxide readings that are below a 0.9 R-squared value

CO2_df <- slope_df %>% filter(R2_CO2 >= 0.9)
CH4_df <- slope_df %>% filter(R2_CH4 >= 0.9)


## GET SGF FIELD DATA ##
chamber_data <- read.csv("../Flux_Chamber_Field_Data.csv", header = T, sep = ",")

chamber_data <- chamber_data %>%
  mutate(DATE = format(as.Date(Date, format = "%m/%d/%Y"), "%Y-%m-%d"), Date = NULL) %>%
  rename(Temp_1 = Temp.1...C., Temp_2 = Temp..2...C., Temp_3 = Temp.3...C.,
         Baro_pressure = Barometric.Pressure..hPa., Date = DATE, Reach = Location) %>%
  rowwise() %>%
  mutate(AvgTemp = mean(c(Temp_1, Temp_2, Temp_3)), Temp_1 = NULL, Temp_2 = NULL, Temp_3 = NULL) %>%
  mutate(Baro_pressure_atm = Baro_pressure * 0.0009869233, Baro_pressure = NULL) %>%
  arrange(Date)

## GET HYDRAGO FILES ##
hydrago_data <- read.csv("../HydraGO_2021.csv", header = T, sep = ",")  

names(hydrago_data)

reach_vec <- c(rep("BDA", 6), rep("REF", 6))

hydrago_data <- hydrago_data %>% 
  mutate(Date = format(as.Date(Timestamp, format = "%m/%d/%Y"), "%Y-%m-%d")) %>%
  select(Date, Location, Soil.Moisture..wfv., Soil.Moisture...., Soil.Temperature..C., Bulk.EC.TC...dS.m.) %>%
  rename(Site = Location, Soil_Moisture_wfv = Soil.Moisture..wfv., Soil_Moisture_Percent = Soil.Moisture....,
         Soil_Temp_C = Soil.Temperature..C., EC = Bulk.EC.TC...dS.m.) %>%
  mutate(Site = case_when(Site == "Blackfoot - Lost Prairie" ~ "LP", 
                              Site == "Blackfoot - Fish Creek" ~ "FH", 
                              Site == "Lolo - Teepee Creek" ~ "TP")) %>%
  distinct(Date, Site, Soil_Moisture_wfv, Soil_Moisture_Percent, Soil_Temp_C, EC) %>%
  mutate(Reach = rep(reach_vec, 24)) %>%
  mutate(Chamber = rep(1:12, 24))

chamber_dates <- chamber_data %>% distinct(Date)
hydrago_dates <- hydrago_data %>% distinct(Date)
field_bind <- bind_cols(chamber_dates, hydrago_dates)

# Combining chamber and hydrago data #           
field_data <- inner_join(chamber_data, hydrago_data)

# Checking all dates line up
field_dates <- field_data %>% distinct(Date)
CO2_dates <- CO2_df %>% distinct(Date)
CH4_dates <- CH4_df %>% distinct(Date)
bind <- bind_cols(field_dates, CO2_dates, CH4_dates)

# New dataframes to reference!                       
CO2_data <- inner_join(CO2_df, field_data)                    
map(CO2_data, ~sum(is.na(.)))

CH4_data <- inner_join(CH4_df, field_data)
map(CH4_data, ~sum(is.na(.)))

# FLUX CALCULATIONS #
CO2_fluxes <- CO2_data %>%
  rowwise() %>%
  mutate(CO2am = (CO2_slope*12*Baro_pressure_atm)/(0.0821*(AvgTemp + 273))) %>%
  mutate(CO2_flux_micrograms = (CO2am*12.31*60)/0.062635) %>%
  mutate(CO2_flux_g = CO2_flux_micrograms /1e+6)

CH4_fluxes <- CH4_data %>%
  rowwise() %>%
  mutate(CH4am = (CH4_slope*12*Baro_pressure_atm)/(0.0821*(AvgTemp + 273))) %>%
  mutate(CH4_flux_micrograms = (CH4am*12.31*60)/0.062635) %>%
  mutate(CH4_flux_g = CH4_flux_micrograms*24/1e+6) # g C/m^2/day

write.csv(CO2_fluxes, 
          file = "/Users/hilaryschultz/Desktop/Thesis/R_Picarro/Picarro_Analyses/CO2_fluxes.csv")

write.csv(CH4_fluxes, 
          file = "/Users/hilaryschultz/Desktop/Thesis/R_Picarro/Picarro_Analyses/CH4_fluxes.csv")


