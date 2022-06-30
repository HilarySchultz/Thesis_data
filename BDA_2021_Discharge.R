# library(Rmisc)
# library(conflicted)
library(dplyr)
library(tidyverse)
# library(lme4)
# library(lsmeans)
# library(car)
# library(readtext)
# library(vroom)
library(lubridate)

#### Discharge ####
# Just another way to change jday to date

# time_vec <- pull(Discharge, jday)
# origin <- as.Date("2021-01-01")
# Discharge$Date <- origin + time_vec
Discharge_2019 <- read.csv("BDA_Discharge.csv", header = T, sep = ",") %>%
  filter(year != 2020, year != 2021, creek != "North Fork Howard", 
         creek != "Little Fish", creek != "Lost Horse") %>%
  mutate(creek = case_when(creek == "Teepee" ~"TP", 
                           creek == "Fish" ~ "FH",
                           creek == "Lost Prairie" ~ "LP")) %>%
  add_column(Origin = as.Date("2020-01-01")) %>%
  mutate(Date = Origin + jday, .keep = "unused") %>%
  select(-X) %>%
  rename(Station = us.ds, Site = creek) %>%
  mutate(Date = format(Date, "%m/%d/%Y")) %>% # Makes date format match DOC dates
  na.exclude()

Discharge_2020 <- read.csv("BDA_Discharge.csv", header = T, sep = ",") %>%
  filter(year != 2019, year != 2021, creek != "North Fork Howard", 
         creek != "Little Fish", creek != "Lost Horse") %>%
  mutate(creek = case_when(creek == "Teepee" ~"TP", 
                           creek == "Fish" ~ "FH",
                           creek == "Lost Prairie" ~ "LP")) %>%
  add_column(Origin = as.Date("2020-01-01")) %>%
  mutate(Date = Origin + jday, .keep = "unused") %>%
  select(-X) %>%
  rename(Station = us.ds, Site = creek) %>%
  mutate(Date = format(Date, "%m/%d/%Y")) %>% # Makes date format match DOC dates
  na.exclude()

Discharge_2021 <- read.csv("BDA_Discharge.csv", header = T, sep = ",") %>%
  filter(year != 2019, year != 2020, creek != "North Fork Howard", 
         creek != "Little Fish", creek != "Lost Horse") %>%
  mutate(creek = case_when(creek == "Teepee" ~"TP", 
                           creek == "Fish" ~ "FH",
                           creek == "Lost Prairie" ~ "LP")) %>%
  add_column(Origin = as.Date("2021-01-01")) %>%
  mutate(Date = Origin + jday, .keep = "unused") %>%
  select(-X) %>%
  rename(Station = us.ds, Site = creek) %>%
  mutate(Date = format(Date, "%m/%d/%Y")) %>% # Makes date format match DOC and SPOC dates
  na.exclude()

# csv to transfer into the other script
write.csv(Discharge_2021, "Discharge_2021.csv", row.names = F)

# LP has missing Spring discharge from 2021, so we need to use the data from 2020
# to create a lm to predict Spring discharge for 2021. 

# Below steps are pulling out the discharge for each station at LP and combining into df

# 2019 data
LP_flow_19 <- Discharge_2019 %>%
  filter(Site == "LP") 
  
LP_DS_19 <- Discharge_2019 %>%
  filter(Site == "LP", Station == "DS") %>%
  rename(DS_discharge = mean.discharge.cms) %>%
  select(-Station)

LP_US_19 <- Discharge_2019 %>%
  filter(Site == "LP", Station == "US") %>%
  rename(US_discharge = mean.discharge.cms) %>%
  select(-Station)

LP_Spring_19 <- Discharge_2019 %>%
  filter(Site == "LP", Station == "Spring") %>%
  rename(Spring_discharge = mean.discharge.cms) %>%
  select(-Station)

LP_19_half <- full_join(LP_DS_19, LP_Spring_19) %>%
  mutate(LP_BDA_test = DS_discharge - Spring_discharge) %>%
  na.exclude()

LP_2019 <- full_join(LP_19_half, LP_US_19) %>%
  na.exclude()


## 2020 data ##
LP_flow_20 <- Discharge_2020 %>%
  filter(Site == "LP") 

LP_DS_20 <- Discharge_2020 %>%
  filter(Site == "LP", Station == "DS") %>%
  rename(DS_discharge = mean.discharge.cms) %>%
  select(-Station)

LP_US_20 <- Discharge_2020 %>%
  filter(Site == "LP", Station == "US") %>%
  rename(US_discharge = mean.discharge.cms) %>%
  select(-Station)

LP_Spring_20 <- Discharge_2020 %>%
  filter(Site == "LP", Station == "Spring") %>%
  rename(Spring_discharge = mean.discharge.cms) %>%
  select(-Station)

LP_20_half <- full_join(LP_DS_20, LP_Spring_20) %>%
  mutate(LP_BDA_test = DS_discharge - Spring_discharge) %>%
  na.exclude()

LP_2020 <- full_join(LP_20_half, LP_US_20) %>%
  na.exclude()

LP_spring_test <- LP_2020 %>%
  mutate(spring_test = DS_discharge - US_discharge) %>%
  select(spring_test, Date) %>%
  rename(mean.discharge.cms = spring_test) %>%
  mutate(Site = "LP", Station = "Spring_test")

LP_BDA <- LP_20_half %>%
  select(LP_BDA_test, Date) %>%
  rename(mean.discharge.cms = LP_BDA_test) %>%
  mutate(Site = "LP", Station = "BDA")

LP_flow_20 <- full_join(LP_flow_20, LP_BDA)
LP_flow_20 <- full_join(LP_flow_20, LP_spring_test)

# Model using data from 2020 on estimated BDA discharge
LP_2020 %>%
  ggplot(aes(DS_discharge, LP_BDA_test)) +
  geom_point() +
  geom_smooth(method = 'lm', se = F)

flow_lm_20 <- lm(LP_BDA_test~DS_discharge * US_discharge, data = LP_2020)
flow_summary_20 <- summary(flow_lm_20)
flow_Rsquared_20 <- summary(flow_lm_20)$r.squared # 0.9489

LP_Q_2021_test <- as_tibble(predict(flow_lm_20, newdata = c(LP_DS_2021, LP_US_2021))) %>%
  rename(mean.discharge.cms = value) %>%
  mutate(Station = "BDA", Site = "LP") %>%
  add_column(LP_Date_2021) %>%
  full_join(LP_DSUS) 

LP_Q_2021_test %>% 
  ggplot() +
  geom_line(aes(x = Date, 
                y = mean.discharge.cms, 
                group = Station, 
                color = Station))
# Plotting Spring discharge against US discharge to determine the correlation.
# First plot shows all 3 stations
LP_flow_20 %>% 
  ggplot() +
  geom_line(aes(x = Date, 
                y = mean.discharge.cms, 
                group = Station, 
                color = Station))

Discharge_2019 %>% 
  filter(Site == "LP") %>%
  ggplot() +
  geom_line(aes(x = Date, 
                y = mean.discharge.cms, 
                group = Station, 
                color = Station))
Discharge_2021 %>% 
  filter(Site == "FH") %>%
  ggplot() +
  geom_line(aes(x = Date, 
                y = mean.discharge.cms, 
                group = Station, 
                color = Station))


# Plot shows relationship between US and Spring station
LP_2019 %>%
  ggplot(aes(US_discharge, Spring_discharge)) +
  geom_point() +
  geom_smooth(method = 'lm', se = F)

LP_2020 %>%
  ggplot(aes(US_discharge, Spring_discharge)) +
  geom_point() +
  geom_smooth(method = 'lm', se = F)

# Percent differences between stations to show which ones are least different
flow_test <- LP_2020 %>%
  mutate(spring_test = DS_discharge - US_discharge,
         spring_percent_diff = (abs(spring_test - Spring_discharge))/((spring_test + Spring_discharge)/2) * 100,
         us_percent_diff = (abs(US_discharge - Spring_discharge))/((US_discharge + Spring_discharge)/2) * 100,)

# Creating a linear model to determine the relationship between US and Spring Q
# Model to predict Spring Q from US Q
flow_lm <- lm(Spring_discharge~US_discharge, data = flow_test)
flow_summary <- summary(flow_lm)
flow_Rsquared <- summary(flow_lm)$r.squared # r.squared = 0.92

flow_test$Date <- as.factor(flow_test$Date)
flow_lm1 <- lm(Spring_discharge~US_discharge*DS_discharge, data = flow_test)
flow_summary1 <- summary(flow_lm1)
flow_Rsquared1 <- summary(flow_lm1)$r.squared # 0.9369

flow_lm2 <- lm(Spring_discharge~DS_discharge, data = flow_test)
flow_summary2 <- summary(flow_lm2)
flow_Rsquared2 <- summary(flow_lm2)$r.squared # 0.9368887

cor(flow_test$Spring_discharge, flow_test$US_discharge)
cor(flow_test$US_discharge, flow_test$Spring_discharge)

lm_test <- full_join(LP_2019, LP_2020)


lm_2019 <- lm(Spring_discharge~US_discharge, data = lm_test)
flow_summary <- summary(lm_2019)
flow_Rsquared <- summary(lm_2019)$r.squared #0.77

# Getting Q for just LP
LP_2021 <- Discharge_2021 %>%
  filter(Site == "LP")
# map(LP_2021, ~sum(is.na(.)))

# Plot of 2021 US and DS Q
LP_2021 %>% 
  ggplot() +
  geom_line(aes(x = Date, 
                y = mean.discharge.cms, 
                group = Station, 
                color = Station))

# Q column for lm prediction
LP_US_2021 <- LP_2021 %>%
  filter(Station == "US") %>%
  select(mean.discharge.cms) %>%
  rename(US_discharge = mean.discharge.cms) %>%
  slice(1:67)

LP_DS_2021 <- LP_2021 %>%
  filter(Station == "DS") %>%
  select(mean.discharge.cms) %>%
  rename(DS_discharge = mean.discharge.cms) 

# US and DS discharge df
LP_DSUS <- LP_2021 %>%
  filter(Station != "Spring") %>%
  select(-year)
  
# Date column
LP_Date_2021 <- LP_2021 %>%
  filter(Station == "US") %>%
  select(Date) %>%
  slice(1:67)

# Predicting the spring discharge from the lm
# DS Q only goes to 8/18/21
LP_Q_2021 <- as_tibble(predict(flow_lm1, newdata = c(LP_US_2021, LP_DS_2021))) %>%
  rename(mean.discharge.cms = value) %>%
  mutate(Station = "Spring", Site = "LP") %>%
  add_column(LP_Date_2021) %>%
  full_join(LP_DSUS) 


# Plot for 2021 Q (DS, US, Spring)
LP_Q_2021 %>% 
  ggplot() +
  geom_line(aes(x = Date, 
                y = mean.discharge.cms, 
                group = Station, 
                color = Station))


# Calculating the BDA discharge - need to separate the discharge columns first.
DS <- LP_Q_2021 %>%
  filter(Station == "DS") %>%
  select(mean.discharge.cms, Date) %>%
  rename(DS_Q = mean.discharge.cms)
  
Spring <- LP_Q_2021 %>%
  filter(Station == "Spring") %>%
  select(mean.discharge.cms, Date) %>%
  rename(Spring_Q = mean.discharge.cms)

BDA <- full_join(DS, Spring) %>%
  mutate(BDA_Q = DS_Q - Spring_Q) %>%
  select(Date, BDA_Q) %>%
  rename(mean.discharge.cms = BDA_Q) %>%
  mutate(Station = "BDA", Site = "LP")

LP_Q_2021 <- full_join(LP_Q_2021, BDA) %>%
  na.exclude()

# .csv to transfer into other script
write.csv(LP_Q_2021, "LP_Q_2021.csv", row.names = F)

# Starts 4/25/2021
TP_Q_dates <- Discharge_2021 %>% 
  filter(Site == "TP")

# Starts 6/2/2021
FH_Q_dates <- Discharge_2021 %>% 
  filter(Site == "FH")

# Starts 6/10/2021
LP_Q_dates <- Discharge_2021 %>% 
  filter(Site == "LP")


