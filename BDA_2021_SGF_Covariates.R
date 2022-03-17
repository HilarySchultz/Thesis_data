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


# Reading in .csv files and removing columns that are not needed to simplify df
SOC <- read.csv("SOC_data.csv", header = T, sep = ",") %>%
  select(Site, Reach, Chamber, C_Conc) # C_conc is g C/g soil

Methane <- read.csv("Methane_data.csv", header = T, sep = ",") %>%
  select(Date, Site, Reach, Chamber, Soil_Moisture_wfv, Soil_Moisture_Percent, Soil_Temp_C, 
         EC, CH4_flux_g)

CDioxide <- read.csv("CarbonDioxide_data.csv", header = T, sep = ",") %>%
  select(Date, Site, Reach, Chamber, Soil_Moisture_wfv, Soil_Moisture_Percent, Soil_Temp_C, 
         EC, CO2_flux_g)

# Joining SOC into the dataframes 
Methane_data <- full_join(Methane, SOC)
CDioxide_data <- full_join(CDioxide, SOC)

