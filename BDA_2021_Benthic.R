# library(Rmisc)
library(dplyr)
library(tidyverse)
# library(lme4)
# library(lsmeans)
# library(car)
# library(readtext)
# library(vroom)
library(lubridate)

#### BENTHIC ####

# Reading in FBPOC data and doing the calculations to compare to google sheets
# Final calculations: Jaunary 2022
FBPOC_data <- read.csv("BDA_FBPOC_Calc.csv", header = T, sep = ",") %>% 
  select(2:8, 14:15) %>%
  rename(Dry_Weight_Filter = Dry.Weight...Filter..g., 
         Ashed_Weight_Filter = Ash...Filter..g., 
         Average_Depth = Average.Depth..cm., 
         Filtered_Volume = Filtered.Volume..ml.) %>%
  mutate(FBPOC_OM = Dry_Weight_Filter - Ashed_Weight_Filter, 
         Benthinator_Volume = pi*8.75^2*Average_Depth, 
         FBPOC_Total_OM = Benthinator_Volume/Filtered_Volume*FBPOC_OM, 
         FBPOC_Carbon_Mass_per_Area = FBPOC_Total_OM/240.53 * 1000 * 0.52)
# Units are in g C/m^2

# Pulling out the average volume column values
Benthinator_volume <- FBPOC_data %>% 
  select(Site, Reach, Location, Replicate, Date, Benthinator_Volume)

# Reading in BPOC data and doing the calculations to compare to google sheets
# You can pass multiple mutate arguments but I am separating them so I can have check points along the way
# Final calculations: Jaunary 2022
BPOC_data <- read.csv("BDA_CBPOC_Calc.csv", header = T, sep = ",") %>% select(3:10, 12) %>% # selecting certain columns
  rename(BPOC_Dry_Weight = Dry.Weight..g., BPOC_Tin_Weight = Tin.Weight..g., # renaming columns
         BPOC_Subsample_Dry_Tin = Subsample.Dry.Weight...Tin..g., BPOC_Ash_Tin = Tin...Ash..g.) %>%
  # "BPOC_Subsample_Dry_Tin" is the weight of the subsample + weight of the tin
  # "BPOC_Ash_Tin" is the weight of either the subsample or sample ash + weight of the tin
  mutate(BPOC_Subsample_Dry_Weight = ifelse(BPOC_Subsample_Dry_Tin > 0, BPOC_Subsample_Dry_Tin - BPOC_Tin_Weight, 0)) %>% 
  # Not all samples have subsamples so if there is a subsample, the ifelse finds which samples have subsamples
  # then subtracts the weight of the tin from the dry subsample in the tin. If no subsample then entry is 0. 
  mutate(BPOC_Ash = BPOC_Ash_Tin - BPOC_Tin_Weight) %>%
  # Determining the mass of the ash, so what's actually there after all OM burns off, so whats left is inorganic.
  mutate(BPOC_Subsample_OM = ifelse(BPOC_Subsample_Dry_Weight > 0, BPOC_Subsample_Dry_Weight - BPOC_Ash, 0)) %>%
  # Determining the mass of what actually burned off, which is the OM content. 
  # If there is a subsample then the ash is subtracted from the dry weight, 
  # meaning that the inorganic material is subtracted from the dry weight (inorganic + OM) to yield the OM content
  mutate(BPOC_Sample_OM = ifelse(BPOC_Subsample_Dry_Tin != 0, 
                                 BPOC_Dry_Weight/BPOC_Subsample_Dry_Weight*BPOC_Subsample_OM, BPOC_Dry_Weight - BPOC_Ash)) %>%
  # Determining the sample OM mass based on a sample/subsample ratio
  # If there is a subsample, then ((dry_sample/dry_subsample) = (x/subsample_OM)), where x = sample OM
  # If there is no subsample, then the sample ash weight is subtracted from the sample dry weight to yield OM content
  mutate(BPOM_mass_per_area = BPOC_Sample_OM/240.53*1000) %>%
  # Sampler area is 240.53 cm^2 (radius of 8.75 cm), when we change to m^2
  # The units here are (g/m^2)
  mutate(BPOC_Carbon_Mass_per_Area = BPOM_mass_per_area * 0.52)
# 52% carbon 
# Units here are (g C/m^2)

# Need the filtered volume for the suspended CBPOC samples (SBPOC)
Benthic_field_data <- read.csv("Benthic_Field_Data.csv", header = T, sep = ",") %>%
  select(2:6, 11) %>%
  rename(Filtered_sample_volume = Filtered.BCPOC..ml.)

volumes <- full_join(Benthic_field_data, Benthinator_volume)

# Reading in SBPOC data and doing the calculations to compare to google sheets
# Final calculations: Jaunary 2022
SBPOC <- read.csv("BDA_SCBPOC_Calc.csv", header = T, sep = ",") %>%
  select(1, 3:10, 12) # selecting certain columns

SBPOC_data <- full_join(SBPOC, volumes) %>% # Joining the benthinator and sample volumes in with the df
  rename(SBPOC_Dry_Weight = Dry.Weight..g., SBPOC_Tin_Weight = Tin.Weight..g., 
         SBPOC_Subsample_Dry_Tin = Subsample.Dry.Weight...Tin..g., SBPOC_Ash_Tin = Tin...Ash..g.) %>%
  mutate(SBPOC_Subsample_Dry_Weight = ifelse(SBPOC_Subsample_Dry_Tin > 0, SBPOC_Subsample_Dry_Tin - SBPOC_Tin_Weight, 0),
         SBPOC_Ash = SBPOC_Ash_Tin - SBPOC_Tin_Weight, 
         SBPOC_Subsample_OM = ifelse(SBPOC_Subsample_Dry_Weight > 0, SBPOC_Subsample_Dry_Weight - SBPOC_Ash, 0), 
         SBPOC_Sample_OM = ifelse(SBPOC_Subsample_Dry_Tin != 0, 
                                  SBPOC_Dry_Weight/SBPOC_Subsample_Dry_Weight*SBPOC_Subsample_OM, SBPOC_Dry_Weight - SBPOC_Ash)) %>%
  # Above are the same steps as the CBPOC, however, I need to normalize by the volume of water that was taken.
  mutate(SBPOM_Total_Mass = (Benthinator_Volume/Filtered_sample_volume)*SBPOC_Sample_OM) %>%
  # Determining the total amount of OM in the total benthinator volume
  # Units are in g
  mutate(SBPOM_mass_per_area = SBPOM_Total_Mass/240.53*1000) %>%
  mutate(SBPOC_Carbon_Mass_per_Area = SBPOM_mass_per_area * 0.52)
# Units are in g C/m^2  

# Binding benthic CBPOC, suspended CBPOC, and FBPOC into one df then selecting rows of interest
CBPOC_data <- full_join(BPOC_data, SBPOC_data) # Coarse samples df

Benthic_data <- full_join(CBPOC_data, FBPOC_data) %>% # Coarse with fine samples into df
  select("Site", "Reach", "Location", "Replicate", "Date", "BPOC_Carbon_Mass_per_Area", 
         "SBPOC_Carbon_Mass_per_Area", "FBPOC_Carbon_Mass_per_Area") %>%
  mutate(Total_BPOC_Mass_per_Area = BPOC_Carbon_Mass_per_Area + SBPOC_Carbon_Mass_per_Area + FBPOC_Carbon_Mass_per_Area)

# write.csv(Benthic_data, "Benthic_data", row.names = F)

Benthic_data %>% group_by(Site, Date, Reach, Location) %>% summarise(AvgBPOC= mean(Total_BPOC_Mass_per_Area)) %>%
  ggplot(aes(Site, AvgBPOC, fill = Site)) + geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 2, color = "blue", fill = "blue") +
  labs(title = "Benthic Particulate Organic Carbon", x = NULL, y = expression(BPOC~(g~C/m^2))) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_blank(), 
                     axis.ticks = element_blank(), axis.text = element_text(colour = "black")) +
  scale_fill_brewer(palette = "Spectral", name = "Site") + 
  facet_grid(cols = vars((Reach))) +
  scale_y_continuous(limits = c(0,120), breaks = c(0,15,30,45,60,75,90,105,120))

# Looking at the distribution
# Benthic_data %>% ggplot() + geom_histogram(aes(Total_BPOC_Mass_per_Area), binwidth = 2)
# Benthic_data %>% ggplot(aes(Total_BPOC_Mass_per_Area)) + geom_histogram() + facet_grid(Site ~ Reach) + 
#   stat_bin(bins = 50) +
#   theme(axis.text.x = element_text(angle = 45))

## GLM for CBPOC
BPOC_glm1 <- glm(Total_BPOC_Mass_per_Area ~ Site*Reach*Location*Date, data = Benthic_data, family = Gamma(link="log"))
Anova(BPOC_glm1) # No interactions were significant. Only Reach, Location, Replicate, and Date have a p-value < 0.5 
AIC(BPOC_glm1) # 1295.187

BPOC_glm2 <- glm(Total_BPOC_Mass_per_Area ~ Reach + Location + Replicate + Date, data = Benthic_data, family = Gamma(link="log") )
Anova(BPOC_glm2) # Now location is the least significant
AIC(BPOC_glm2) # 1314.862

BPOC_glm3 <- glm(Total_BPOC_Mass_per_Area ~ Reach + Replicate + Date, data = Benthic_data, family = Gamma(link="log") )
Anova(BPOC_glm3) 
AIC(BPOC_glm3) # 1321.015

anova(BPOC_glm1, BPOC_glm2, BPOC_glm3)
