library(dplyr) 
library(readr) 
library(ggplot2) 
library(tidyverse)

## GET FILES ##
filelist <-list.files(path = getwd(), pattern="*.csv", full.names = TRUE)
DIC_files <- list()
for(f in 1:length(filelist)) {
  DIC_files[[f]] <- read.csv(filelist[f], header = T, sep = ",", as.is = T)
}

DIC_df <- bind_rows(DIC_files)

DIC_data <- DIC_df %>%
  rename(Segment = Design, 
         Date = DATE, 
         Reach = Location,
         CO2_ppm = avg_CO2,
         CH4_ppm = avg_CH4) %>%
  mutate(Segment = if_else(Segment == "REF", "Reference", "Treatment")) %>%
  mutate(Reach = case_when(Reach == "LWR" ~ "Lower", 
                           Reach == "MID" ~ "Middle", 
                           Reach == "UPR" ~ "Upper")) %>%
  mutate(Date = case_when(Date == "2021-06-01" ~ "06/01/2021",
                          Date == "2021-06-02" ~ "06/01/2021",
                          Date == "2021-06-03" ~ "06/01/2021",
                          Date == "2021-06-14" ~ "06/14/2021",
                          Date == "2021-06-15" ~ "06/14/2021",
                          Date == "2021-06-17" ~ "06/14/2021",
                          Date == "2021-06-28" ~ "06/28/2021",
                          Date == "2021-06-29" ~ "06/28/2021",
                          Date == "2021-07-01" ~ "06/28/2021",
                          Date == "2021-07-12" ~ "07/12/2021",
                          Date == "2021-07-13" ~ "07/12/2021",
                          Date == "2021-07-15" ~ "07/12/2021",
                          Date == "2021-07-26" ~ "07/26/2021",
                          Date == "2021-07-27" ~ "07/26/2021",
                          Date == "2021-07-29" ~ "07/26/2021", 
                          Date == "2021-08-09" ~ "08/09/2021",
                          Date == "2021-08-10" ~ "08/09/2021",
                          Date == "2021-08-12" ~ "08/09/2021",
                          Date == "2021-08-25" ~ "08/25/2021",
                          Date == "2021-08-26" ~ "08/25/2021",
                          Date == "2021-08-29" ~ "08/25/2021",
                          Date == "2021-09-07" ~ "09/07/2021",
                          Date == "2021-09-08" ~ "09/07/2021",
                          Date == "2021-09-10" ~ "09/07/2021" ))

# Need to add replicates in here but there were some faulty samples so reps aren't 1:3 always
# The month of June has all the replicates
DICreps <- DIC_data %>%
  filter(Date %in% c("06/01/2021","06/14/2021","06/28/2021", 
                     "07/12/2021", "07/26/2021","08/25/2021", "09/07/2021")) %>%
  mutate(Replicate = rep(1:3, 126))

shortrep <- DIC_data %>%
  filter(Date == "08/09/2021",
         Site == "LP",
         Segment == "Treatment",
         Reach == "Upper") %>%
  mutate(Replicate = 1)

DICshortreps <- DIC_data %>%
  filter(Date == "08/09/2021") %>%
  slice(-43) %>%
  mutate(Replicate = rep(1:3, 17))

DICmissingreps <- full_join(shortrep, DICshortreps)

DIC_data <- full_join(DICreps, DICmissingreps)

# Separated into CO2 and CH4 dfs
DIC_CO2 <- DIC_data %>%
  select(-Triplicate_CH4_avg, -avg_deltaCH4, -CH4_ppm, -Triplicate_delta_CH4) 

DIC_CH4 <- DIC_data %>%
  select(-Triplicate_CO2_avg, -avg_deltaCO2, -CO2_ppm, -Triplicate_CO2_avg) 

# Reading alkalinity file
Alkalinity <- read.csv("../Alkalinity_BDA.csv", header = T, sep = ",")

Alkalinity_data <- Alkalinity %>%
  rename(Date = Date.Sample.Taken,
         Segment = Reach,
         alkalinity_µmolperkg. = Alkalinity..µmol.kg.) %>%
  filter(Segment != "DS") %>%
  select(-File.Name, -X, -LA..ml., -pHi, -pHf, -End.T.C, -Analysis.Date) %>%
  mutate(Segment = case_when(Segment == "BDA" ~ "Treatment",
                             Segment == "REF" ~ "Reference")) %>%
  mutate(Date = case_when(Date == "6/1/2021" ~ "06/01/2021",
                          Date == "6/2/2021" ~ "06/01/2021",
                          Date == "6/3/2021" ~ "06/01/2021",
                          Date == "6/14/2021" ~ "06/14/2021",
                          Date == "6/15/2021" ~ "06/14/2021",
                          Date == "6/17/2021" ~ "06/14/2021",
                          Date == "6/28/2021" ~ "06/28/2021",
                          Date == "6/29/2021" ~ "06/28/2021",
                          Date == "7/1/2021" ~ "06/28/2021",
                          Date == "7/12/2021" ~ "07/12/2021",
                          Date == "7/13/2021" ~ "07/12/2021",
                          Date == "7/15/2021" ~ "07/12/2021",
                          Date == "7/26/2021" ~ "07/26/2021",
                          Date == "7/27/2021" ~ "07/26/2021",
                          Date == "7/29/2021" ~ "07/26/2021", 
                          Date == "8/9/2021" ~ "08/09/2021",
                          Date == "8/10/2021" ~ "08/09/2021",
                          Date == "8/12/2021" ~ "08/09/2021",
                          Date == "8/25/2021" ~ "08/25/2021",
                          Date == "8/26/2021" ~ "08/25/2021",
                          Date == "8/29/2021" ~ "08/25/2021",
                          Date == "9/7/2021" ~ "09/07/2021",
                          Date == "9/8/2021" ~ "09/07/2021",
                          Date == "9/10/2021" ~ "09/07/2021" ))