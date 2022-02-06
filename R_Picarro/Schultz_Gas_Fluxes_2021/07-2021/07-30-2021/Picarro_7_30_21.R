library(dplyr) 
library(readr) 
library(ggplot2) 
library(tidyverse)


## GET FILES ##
filelist <-list.files(path = getwd(), pattern="*.dat", full.names = TRUE)
rawdata <- list()
for(f in 1:length(filelist)) {
  rawdata[[f]] <- read.delim(filelist[f], header = T, sep = "", as.is = T)
}

names(rawdata)
rawdata <- bind_rows(rawdata) %>% select(DATE, TIME, EPOCH_TIME, X12CO2_dry, HP_12CH4, 
                                         Delta_Raw_iCO2, HP_Delta_iCH4_Raw)
map(rawdata, ~sum(is.na(.)))
rawdata <- na.omit(rawdata)
map(rawdata, ~sum(is.na(.)))

DIC_data <- rawdata %>% filter(INDEX < 7155)
SGF_data <- rawdata %>% filter(INDEX > 7155)

## PLOTTING RAW DATA##
plot(rawdata$X12CO2_dry, type = "l")
lines(DIC_data$X12CO2_dry, col = "red")
lines(rawdata$HP_12CH4, type = "l", col = 4)

## NEW BASELINE ##
# filtered_data <- dplyr::if_else(rawdata$X12CO2_dry>150, rawdata$X12CO2_dry, 0.0) 
DIC_filtered_data <- if_else(DIC_data$X12CO2_dry>75, DIC_data$X12CO2_dry, 0.0)
SGF_filtered_data <- if_else(SGF_data$X12CO2_dry > 300, SGF_data$X12CO2_dry, 0.0)

DIC_index_stamp <- 1:length(DIC_data$EPOCH_TIME)
SGF_index_stamp <- 1:length(SGF_data$EPOCH_TIME)

DIC_data$INDEX <- DIC_index_stamp
SGF_data$INDEX <- SGF_index_stamp

## GET PEAKS ##
get_peaks <- function(data_vec, index_vec){
  peak_num <- c()
  peak_value <- c()
  current_peak <- 0
  peak_index <- c()
  for(i in 2:length(data_vec)){
    if (data_vec[i-1] == 0 & data_vec[i] > 0){
      current_peak <- current_peak + 1
    }
    if (data_vec[i-1] != 0 & data_vec[i] != 0) {
      peak_value <- append(peak_value, data_vec[i-1])
      peak_num <- append(peak_num, current_peak)
      peak_index <- append(peak_index, index_vec[i-1])
    }
  }
  peak_df <- data.frame(peak_num = peak_num,
                        peak_value = peak_value, peak_index = peak_index
                        
  )
  return(peak_df)
}


## USE GET PEAKS ON FILTERED DATA ##
# This grabs the values > 300 (or whatever is speficied in the get_peaks()) for each peak
DIC_filtered_data_peaks <- get_peaks(DIC_filtered_data, DIC_index_stamp)
SGF_filtered_data_peaks <- get_peaks(SGF_filtered_data, SGF_index_stamp)

## GET INDEXES ##
# get_index takes the range of values obtained for each peak in the previous step and arranges
# them in descending order then takes the highest 30 (or whatever is speficied values)
get_index <- function(peaks_df){
  index_stamp <- c()
  peak_num <- c()
  for (i in 1:max(peaks_df$peak_num)){ # for each peak
    single_peak <- peaks_df %>% filter(peak_num==i)
    single_peak_sort <- single_peak %>% arrange(desc(peak_value)) %>% slice(1:120) %>% select(peak_value, peak_index, peak_num)
    index_stamp <- append(index_stamp, single_peak_sort$peak_index)
    peak_num <- append(peak_num, single_peak_sort$peak_num)
  }
  index_stamp_df <- data.frame(index_stamp = index_stamp, peak_num = peak_num)
  return(index_stamp_df)
}

## USE GET INDEX ON FILTERED DATA PEAK AND INDEX RAWDATA ##
DIC_filtered_index_num <- get_index(DIC_filtered_data_peaks)
SGF_filtered_index_num <- get_index(SGF_filtered_data_peaks)

# Peak tally
DIC_filtered_index_num %>% distinct(peak_num)
SGF_filtered_index_num %>% distinct(peak_num)

###################################################################
#### STOP STOP STOP, DON'T RUN BELOW UNTIL PEAKS ARE EVALUATED ####
# Matches index_stamp and peak_num with rawdata
DIC_filtered_rawdata <- DIC_data %>% 
  filter(INDEX %in% DIC_filtered_index_num$index_stamp) %>%
  add_column(peak = DIC_filtered_index_num$peak_num) %>%
  arrange(INDEX) 

DIC_filtered_rawdata %>% distinct(peak) # 18 peaks

SGF_filtered_rawdata <- SGF_data %>% 
  filter(INDEX %in% SGF_filtered_index_num$index_stamp) %>%
  add_column(peak = SGF_filtered_index_num$peak_num) %>%
  arrange(INDEX) 

SGF_peaks <- SGF_filtered_rawdata %>% distinct(peak) # 36 peaks

# Double check there are the correct number of peaks after slice the erroneous peaks
# filtered_rawdata %>% distinct(peak)

# Making sure plot looks correct without unwanted peaks
plot(DIC_filtered_rawdata$X12CO2_dry, type = "l")
plot(SGF_filtered_rawdata$X12CO2_dry, type = "l")
# plot(filtered_rawdata$HP_12CH4, type = "l")

# SEPARATE CO2 DATA INTO ITS OWN DATAFRAME FOR DIC AND SGF
DIC_carbondioxide_df <- DIC_filtered_rawdata %>%
  select(-HP_12CH4, -HP_Delta_iCH4_Raw) %>%
  arrange(INDEX)

SGF_carbondioxide_df <- SGF_filtered_rawdata %>%
  select(-HP_12CH4, -HP_Delta_iCH4_Raw) %>%
  arrange(INDEX)

# SEPARATE CH4 DATA INTO ITS OWN DATAFRAME FOR DIC AND SGF
DIC_methane_df <- DIC_filtered_rawdata %>% 
  select(-X12CO2_dry, -Delta_Raw_iCO2) %>%
  arrange(INDEX)

SGF_methane_df <- SGF_filtered_rawdata %>% 
  select(-X12CO2_dry, -Delta_Raw_iCO2) %>%
  arrange(INDEX)

## CO2 FUNCTION TO GET THE BEST PEAK WINDOW ##
CO2_clean_peaks <- function(whole_peak){
  # window_size is the adjustable hyper parameter; each reading is a set of 4 (36/4 = 8 data points)
  window_size = 36
  best_df = NA
  min_std = 1e100
  for (v in 1:(length(whole_peak$X12CO2_dry) - window_size)){
    CO2_peak <- whole_peak %>% slice(v: (v + window_size - 1))
    sd_CO2 <- sd(CO2_peak$X12CO2_dry)
    if (sd_CO2 < min_std){
      min_std <- sd_CO2
      best_df <- CO2_peak
    }
  }
  return(best_df)
}


## CH4 FUNCTION TO GET THE BEST PEAK WINDOW ##
CH4_clean_peaks <- function(whole_peak){
  # window_size is the adjustable hyper parameter; each reading is a set of 4 (36/4 = 8 data points)
  window_size = 36
  best_df = NA
  min_std = 1e100
  for (v in 1:(length(whole_peak$HP_12CH4) - window_size)){
    CH4_peak <- whole_peak %>% slice(v: (v + window_size - 1))
    sd_CH4 <- sd(CH4_peak$HP_12CH4)
    if (sd_CH4 < min_std){
      min_std <- sd_CH4
      best_df <- CH4_peak
    }
  }
  return(best_df)
}

## CALLING CLEAN_PEAKS ON THE DAYS DATA FOR CO2 AND CH4 ##
DIC_carbondioxide_cleaned_df <- DIC_carbondioxide_df %>%
  group_by(peak) %>%
  nest() %>%
  mutate(CO2_cleaned_peaks = map(data, CO2_clean_peaks)) %>%
  unnest(CO2_cleaned_peaks) %>%
  select(-data)

SGF_carbondioxide_cleaned_df <- SGF_carbondioxide_df %>%
  group_by(peak) %>%
  nest() %>%
  mutate(CO2_cleaned_peaks = map(data, CO2_clean_peaks)) %>%
  unnest(CO2_cleaned_peaks) %>%
  select(-data)

DIC_methane_cleaned_df <- DIC_methane_df %>%
  group_by(peak) %>%
  nest() %>%
  mutate(CH4_cleaned_peaks = map(data, CH4_clean_peaks)) %>%
  unnest(CH4_cleaned_peaks) %>%
  select(-data)

SGF_methane_cleaned_df <- SGF_methane_df %>%
  group_by(peak) %>%
  nest() %>%
  mutate(CH4_cleaned_peaks = map(data, CH4_clean_peaks)) %>%
  unnest(CH4_cleaned_peaks) %>%
  select(-data)

## SUMMARIZED DATAFRAME FOR CARBON DIOXIDE AND METHANE DATAFRAMES ##
DIC_CO2_final_df <- DIC_carbondioxide_cleaned_df %>%
    group_by(peak, DATE) %>%
    summarise(avg_CO2 = mean(X12CO2_dry), avg_deltaCO2=mean(Delta_Raw_iCO2))

SGF_CO2_final_df <- SGF_carbondioxide_cleaned_df %>%
  group_by(peak, DATE) %>%
  summarise(avg_CO2 = mean(X12CO2_dry), avg_deltaCO2=mean(Delta_Raw_iCO2))

DIC_CH4_final_df <- DIC_methane_cleaned_df %>%
  group_by(peak, DATE) %>%
  summarise(avg_CH4 = mean(HP_12CH4), avg_deltaCH4 = mean(HP_Delta_iCH4_Raw))

SGF_CH4_final_df <- SGF_methane_cleaned_df %>%
  group_by(peak, DATE) %>%
  summarise(avg_CH4 = mean(HP_12CH4), avg_deltaCH4 = mean(HP_Delta_iCH4_Raw))

DIC_Final_gas_df <- inner_join(DIC_CO2_final_df, DIC_CH4_final_df) 

SGF_Final_gas_df <- inner_join(SGF_CO2_final_df, SGF_CH4_final_df)

# location_vec for 18 samples (one site)
location_vec <- rep(c(rep("LWR", 3), rep("MID", 3), rep("UPR", 3)), 2)

# design_vec for 18 samples (one site)
design_vec <- c(rep("BDA", 9), rep("REF", 9))

#######################################################
### MAKE SURE THE SLICES FOR EACH SAMPLE COLLECTION ARE CORRECT ###
## SUBSETTING FINAL_GAS_DF INTO DIC AND SGF ##
July_30_DIC <- as_tibble(DIC_Final_gas_df) %>%
  mutate(DATE = replace(DATE, DATE == "2021-07-30","2021-07-29")) %>%
  mutate(Site = "TP") %>%
  mutate(Location = location_vec) %>%
  mutate(Design = design_vec) %>%
  group_by(Site, Design, Location) %>% 
  mutate(Triplicate_CO2_avg = mean(avg_CO2), 
         Triplicate_CH4_avg = mean(avg_CH4), 
         Triplicate_delta_CO2_avg = mean(avg_deltaCO2),
         Triplicate_delta_CH4 = mean(avg_deltaCH4)) %>%
  ungroup()

# DIC .csv
write.csv(July_30_DIC, 
          file = "/Users/hilaryschultz/Desktop/Thesis/R_Picarro/Picarro_DIC_2021/July_30_DIC.csv")


chamber_vector <- c(rep.int(1,3), rep.int(2,3), rep.int(3,3), rep.int(4,3), rep.int(5,3),
                    rep.int(6,3), rep.int(7,3), rep.int(8,3), rep.int(9,3), rep.int(10,3),
                    rep.int(11,3), rep.int(12,3))

# SGF .csv
July_30_SGF <- as_tibble(SGF_Final_gas_df) %>% 
  mutate(Site = "TP") %>%
  mutate(DATE = replace(DATE, DATE == "2021-07-30","2021-07-29")) %>%
  mutate(chamber = chamber_vector) %>%
  mutate(Design = if_else(chamber <= 6, "BDA", "REF")) 

write.csv(July_30_SGF, 
          file = "/Users/hilaryschultz/Desktop/Thesis/R_Picarro/Picarro_SGF_2021/July_30_SGF.csv")


plot(rawdata$EPOCH_TIME, rawdata$HP_12CH4, type = "l")
points(DIC_methane_cleaned_df$EPOCH_TIME, DIC_methane_cleaned_df$HP_12CH4, col = "green")
points(SGF_methane_cleaned_df$EPOCH_TIME, SGF_methane_cleaned_df$HP_12CH4, col = "green")

# THIS FILE HAS BEEN RUN 11/19/21 