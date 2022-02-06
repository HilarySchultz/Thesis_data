# library(Rmisc)
library(dplyr)
library(tidyverse)
# library(lme4)
# library(lsmeans)
# library(car)
# library(readtext)
# library(vroom)
library(lubridate)

#### SOC ####
# Reading in SOC data, selecting columns of interest, and renaming those columns
SOC_data <- read.csv("BDA_SOC.csv", header = T, sep = ",") %>% select(1:4, 6:9, 11:12, 14:15) %>%
  rename(Plaster_Volume = Plaster.Volume..cm.3., Tin_Weight = Tin.Weight..g.,
         Total_Soil_Wet = Total.Sample.Wet..g., Subsample_Wet_Tin = Subsample.Wet..g., 
         Gravel_Weight = Gravel..g., Subsample_Dry_Tin = Pre.grind.Dry.Subsample...Tin..g.,
         Ground_Subsample_Tin = Re.dried.Subsample..g., 
         Ground_Subsample_Ashed_Tin = Ground.Subsample.Ashed...Tin..g.)

# Mutating data for calculations to check against Excel document
SOC_data <- SOC_data %>% mutate(Gravel_Volume = Gravel_Weight/2.65, # Gravel Volume
                                Subsample_Wet = Subsample_Wet_Tin - Tin_Weight, # Subsample weight without tin
                                Sample_Subsample_Ratio = Total_Soil_Wet/Subsample_Wet, # Sample/Subsample Ratio
                                Subsample_Dry = Subsample_Dry_Tin - Tin_Weight, # Subsample dry without tin
                                Total_Soil_Dry = Sample_Subsample_Ratio*Subsample_Dry, # Sample dry weight
                                Soil_Bulk_Density = Total_Soil_Dry/(Plaster_Volume - Gravel_Volume),
                                Ground_Subsample_Dry = Ground_Subsample_Tin - Tin_Weight, # Ground subsample without tin
                                Ground_Subsample_Ashed = Ground_Subsample_Ashed_Tin - Tin_Weight, # Ground ashed subsample without tin 
                                Subsample_Ashed = Ground_Subsample_Ashed/Ground_Subsample_Dry*Subsample_Dry, # Subsample ashed accounting for the change in mass during grinding
                                Total_Soil_Ashed = Total_Soil_Dry/Subsample_Dry*Subsample_Ashed, # Ashed weight of sample
                                Total_Soil_OM = Total_Soil_Dry - Total_Soil_Ashed, # Sample OM content
                                Total_Soil_C = Total_Soil_OM * 0.52,
                                SOC_Bulk_Density = Total_Soil_C/(Plaster_Volume - Gravel_Volume)) 

write.csv(SOC_data, "SOC_data.csv", row.names = F)
names(SOC_data)
range(SOC_data$SOC_Bulk_Density)

SOC_data %>% 
  ggplot(aes(Site, SOC_Bulk_Density, fill = Site)) + geom_boxplot() +
  labs(title = "Soil Organic Carbon Bulk Density", x = NULL, y = expression(SOC~Bulk~Density~(g~C/cm^3))) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 2, color = "blue", fill = "blue") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_blank(), 
                     axis.ticks = element_blank(), axis.text = element_text(colour = "black")) +
  scale_fill_brewer(palette = "Spectral", name = "Site") + 
  facet_grid(cols = vars(Reach)) +
  scale_y_continuous(limits = c(0,0.08))

# SOC_data <- read.csv("BDA_SOC.csv", header = T, sep = ",") %>% 
#             rename(Plaster_Volume = Plaster.Volume..cm.3.,
#             Gravel = Gravel..g., Sample_OM = Sample.OM..g., 
#             Sample_Carbon = Sample.Carbon.Content..g.,
#             Ratio = Sample.Subsample.Ratio, 
#             Subsample_Dry = Pre.grind.Dry.Subsample..g.) %>%
#             mutate(Gravel_Volume = Gravel/2.65, 
#             Total_Sample_Dry = Ratio*Subsample_Dry,
#             Soil_Volume = Plaster_Volume - Gravel_Volume,
#             SOC_Bulk_Density = Sample_Carbon/Soil_Volume)
View(SOC_data)
names(SOC_data)