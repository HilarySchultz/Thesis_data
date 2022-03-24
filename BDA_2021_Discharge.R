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
# Discharge for only 2022 #
Discharge <- read.csv("BDA_Discharge.csv", header = T, sep = ",") %>%
  filter(year != 2019, year != 2020, creek != "North Fork Howard", 
         creek != "Little Fish", creek != "Lost Horse") %>%
  mutate(creek = case_when(creek == "Teepee" ~"TP", 
                           creek == "Fish" ~ "FH",
                           creek == "Lost Prairie" ~ "LP")) %>%
  add_column(Origin = as.Date("2021-01-01")) %>%
  mutate(Date = Origin + jday, .keep = "unused") %>%
  select(-X, -year) %>%
  rename(Station = us.ds, Site = creek) %>%
  mutate(Date = format(Date, "%m/%d/%Y")) %>% # Makes date format match DOC dates
  na.exclude()

# Getting data for just LP

LP_DS <- Discharge %>%
  filter(Site == "LP", Station == "DS") %>%
  rename(DS_discharge = mean.discharge.cms) %>%
  select(-Station)

LP_US <- Discharge %>%
  filter(Site == "LP", Station == "US") %>%
  rename(US_discharge = mean.discharge.cms) %>%
  select(-Station)

LP_2021 <- full_join(LP_DS, LP_US)

LP_2021 %>% ggplot(aes(Date, US_discharge)) +
  geom_point() + geom_smooth(method = 'lm', se = F)

LP_2021 %>% ggplot(aes(Date, DS_discharge)) +
  geom_point() + geom_smooth(method = 'lm', se = F)

# 2020 discharge for exterpolation
discharge_2020 <- read.csv("allBDAhydro.csv", header = T, sep = ",") %>%
  filter(year == 2020, creek == "Lost Prairie") %>%
  add_column(Origin = as.Date("2020-01-01")) %>%
  mutate(Date = Origin + jday, .keep = "unused") %>%
  select(-X, -year) %>%
  rename(Station = us.ds, Site = creek) %>%
  mutate(Date = format(Date, "%m/%d/%Y"))

DS_LP <- discharge_2020 %>% 
  filter(Station == "DS") %>%
  rename(DS_discharge = mean.discharge.cms) %>%
  select(-Station) %>%
  filter(!is.na(DS_discharge))

US_LP <- discharge_2020 %>% 
  filter(Station == "US") %>%
  rename(US_discharge = mean.discharge.cms) %>%
  select(-Station) %>%
  filter(!is.na(US_discharge))

Spring_LP <- discharge_2020 %>% 
  filter(Station == "Spring") %>%
  rename(Spring_discharge = mean.discharge.cms) %>%
  select(-Station) %>%
  filter(!is.na(Spring_discharge))

# Full dataframe
discharge_LP <- join_all(list(DS_LP, US_LP, Spring_LP)) %>%
  na.exclude()


discharge_LP %>% ggplot(aes(US_discharge, Spring_discharge)) +
  geom_point() + geom_smooth(method = 'lm', se = F)

discharge_LP %>% ggplot(aes(DS_discharge, Spring_discharge)) +
  geom_point() + geom_smooth(method = 'lm', se = F)

# Just another way to change jday to date

# time_vec <- pull(Discharge, jday)
# origin <- as.Date("2021-01-01")
# Discharge$Date <- origin + time_vec

#### DOC LOADS ####

# Site: TP
TP_dates <- cleaned_DOC %>% # pulling a column of distinct dates
  filter(Site == "TP") %>%
  distinct(Date) %>%
  pull(Date)

TP_DOC <- cleaned_DOC %>% 
  filter(Site == "TP") %>%
  mutate(Station = case_when(Reach == "REF" ~ "US", # Matching gauging station to reach location
                             Reach == "BDA" ~ "DS",
                             Reach == "DS" ~ "DS"))

TP_cms <- Discharge %>% # Matching distinct DOC sample dates to discharge dates
  filter(Site == "TP") %>% 
  filter(Date %in% TP_dates) %>%
  arrange(Date)

TP <- full_join(TP_DOC, TP_cms) %>% # Load df
  mutate(Load = Conc_ppm * mean.discharge.cms) 

# Site: FH
FH_dates <- cleaned_DOC %>%
  filter(Site == "FH") %>%
  distinct(Date) %>%
  pull(Date)

FH_DOC <- cleaned_DOC %>%
  filter(Site == "FH") %>%
  mutate(Station = case_when(Reach == "REF" ~ "US", # Matching gauging station to reach location
                             Reach == "BDA" ~ "DS",
                             Reach == "DS" ~ "DS"))

# No discharge data past 9/1/21, so I need to use that date for 9/7/21
# Pulled discharge data from 9/1/2021 for US and DS
FH_new_date <- Discharge %>%
  filter(Site == "FH", Date == "09/01/2021")

FH_cms <- Discharge %>%
  filter(Site == "FH") %>% 
  filter(Date %in% FH_dates) %>%
  bind_rows(FH_new_date) %>%
  arrange(Date) %>%
  mutate(Date = ifelse(Date == "09/01/2021", "09/07/2021", Date)) 
# Changed the discharge date so that it would match with the DOC date.

FH <- full_join(FH_DOC, FH_cms) %>%
  mutate(Load = Conc_ppm * mean.discharge.cms) 



# Site: LP
LP_dates <- cleaned_DOC %>%
  filter(Site == "LP") %>%
  distinct(Date) %>%
  pull(Date)

LP_discharge <- Discharge %>%
  filter(Site == "LP")

LP_DOC <- cleaned_DOC %>%
  filter(Site == "LP")

LP_cms <- Discharge %>%
  filter(Site == "LP") %>% 
  filter(Date %in% LP_dates) %>%
  arrange(Date) 

LP <- full_join(LP_DOC, LP_cms)
get_FH_cms <- Discharge %>%
  filter(Site =="FH")

get_LP_cms <- Discharge %>%
  filter(Site =="LP")

mutate(Load = Conc_ppm * mean.discharge.cms)

join <- full_join(LP, FH)

DOC_loads <- full_join(join, TP) %>%
  arrange(Date) 

# TP_difference <- TP %>%
#   mutate(US_cms = case_when(Station == "US" ~ mean.discharge.cms),
#          DS_cms = case_when(Station == "DS" ~ mean.discharge.cms))
# LP_difference <- LP_cms %>%
#   group_by(Site, Date) %>%
#   mutate(US_cms = case_when(Station == "US" ~ mean.discharge.cms),
#          DS_cms = case_when(Station == "DS" ~ mean.discharge.cms))


mutate(percent_diff = mean.discharge.cms/lag(mean.discharge.cms)*100)

((a-b)/a+b/2)*100
test2 <- cleaned_DOC %>%
  distinct(Date) %>%
  pull(Date)

testdate <- cleaned_DOC %>%
  distinct(Date)

testdate1 <- Discharge %>%
  filter(Date %in% testdate)


test2.1 <- cleaned_DOC %>%
  distinct(Date) %>%
  pull(Site)

test2_cms <- Discharge %>%
  filter(Date %in% test2 & Site %in% test2.1)

test_join <- full_join(cleaned_DOC, test2_cms)