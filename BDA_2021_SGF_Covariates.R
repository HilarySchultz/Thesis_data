library(multcomp) # for doing elegant pairwise comparisons
library(Rmisc)  # for summarySE function to summarize data frame
library(boot) # for diagnostics for GLM
library(lme4)   # for linear, general, and nonlinear mixed models
library(car)    # for doing ANOVA
library(emmeans) # posthoc tests
library(multcompView) # making posthoc tests easier to view and plot
library(AICcmodavg) # for comparing models
library(PerformanceAnalytics) # for making amazing correlation plots
library(tidyverse)  # loads several useful packages, including ggplot2, 
# tidyr, and dplyr
options(scipen = 999)


# Reading in .csv files and removing columns that are not needed to simplify dfs

## SOC data
SOC <- read.csv("SOC_data.csv", header = T, sep = ",") %>%
  select(Site, Reach, Chamber, C_Conc) # C_conc is g C/g soil

## Methane data
Methane <- read.csv("Methane_data.csv", header = T, sep = ",") %>%
  select(Date, Site, Reach, Chamber, Soil_Moisture_wfv, Soil_Moisture_Percent, Soil_Temp_C, 
         EC, CH4_flux_g)


## Carbon dioxide data
CDioxide <- read.csv("CarbonDioxide_data.csv", header = T, sep = ",") %>%
  select(Date, Site, Reach, Chamber, Soil_Moisture_wfv, Soil_Moisture_Percent, Soil_Temp_C, 
         EC, CO2_flux_g)

# Joining SOC into the dataframes 
Methane_data <- full_join(Methane, SOC)
Methane_data$Date <- ordered(Methane_data$Date, 
                        levels = c("6/7/2021", "6/14/2021", "6/28/2021",
                                   "7/12/2021", "7/26/2021", "8/9/2021", 
                                   "8/25/2021", "9/7/2021"))
Methane_data$Site <- as.factor(Methane_data$Site)
Methane_data$Reach <- as.factor(Methane_data$Reach)
Methane_data$Chamber <- ordered(Methane_data$Chamber, 
                           levels = c(1:12))
##
CDioxide_data <- full_join(CDioxide, SOC)
CDioxide_data$Date <- ordered(CDioxide_data$Date, 
                         levels = c("6/7/2021", "6/14/2021", "6/28/2021",
                                    "7/12/2021", "7/26/2021", "8/9/2021", 
                                    "8/25/2021", "9/7/2021"))
CDioxide_data$Site <- as.factor(CDioxide_data$Site)
CDioxide_data$Reach <- as.factor(CDioxide_data$Reach)
CDioxide_data$Chamber <- ordered(CDioxide_data$Chamber, 
                            levels = c(1:12))

# I think I need to change the data types for this model to run. 
CDioxide_data <- CDioxide_data %>%
  mutate()
### Exploratory Models ###
# Full interaction models
cd_model <- glmer(CO2_flux_g ~ Date*Reach*Soil_Moisture_wfv*Soil_Temp_C*EC*C_Conc + (1|Site) + (1|Chamber), 
                       data = CDioxide_data,
                       family = Gamma(link = "log"))
# I may need to make all of these numeric for this to work, or should I do a random forest?
Anova(cd_model)