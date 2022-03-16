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
SOC_data <- read.csv("BDA_SOC.csv", header = T, sep = ",") %>%
  dplyr::select(1:9, 11:12, 14:15) %>%
  dplyr::rename(Plaster_Volume = Plaster.Volume..cm.3.,
                Plaster_Weight_g = Plaster.Weight..g.,
                Tin_Weight = Tin.Weight..g.,
                Total_Soil_Wet = Total.Sample.Wet..g., 
                Subsample_Wet_Tin = Subsample.Wet..g., 
                Gravel_Weight = Gravel..g., 
                Subsample_Dry_Tin = Pre.grind.Dry.Subsample...Tin..g., 
                Ground_Subsample_Tin = Re.dried.Subsample..g., 
                Ground_Subsample_Ashed_Tin = Ground.Subsample.Ashed...Tin..g.)

# Mutating data for calculations to check against Excel document
# The samples were not actually "wet" they were air-dried enough to sieve but just using
# that term so the different drying levels don't get mixed up. 
# density of plaster of paris (650 to 710 Kg/m^3)
SOC_data <- SOC_data %>% 
  mutate(Subsample_Wet = Subsample_Wet_Tin - Tin_Weight, # Subsample weight
         Sample_Subsample_Ratio = Total_Soil_Wet/Subsample_Wet, # Sample/Subsample Ratio
         # The total soil wet is only the weight of the soil, this exclues the coarse fraction
         Subsample_Dry = Subsample_Dry_Tin - Tin_Weight, # Subsample dry without tin
         Total_Soil_Dry = Sample_Subsample_Ratio*Subsample_Dry, # Sample dry weight
         # The subsamples were ground and weighed again the ashed.
         Ground_Subsample_Dry = Ground_Subsample_Tin - Tin_Weight, # Ground subsample without tin
         Ground_Subsample_Ashed = Ground_Subsample_Ashed_Tin - Tin_Weight, # Ground ashed subsample without tin
         Subsample_Ashed = Ground_Subsample_Ashed/Ground_Subsample_Dry*Subsample_Dry, # Subsample ashed accounting for the change in mass during grinding
         Total_Soil_Ashed = Sample_Subsample_Ratio*Subsample_Ashed, # Ashed weight of sample
         Total_Soil_OM = Total_Soil_Dry - Total_Soil_Ashed, # Sample OM content
         Conc_OM = Total_Soil_OM/Total_Soil_Dry,
         C_Conc = Conc_OM*0.52)  # This is mass C per mass soil
         
  
# We are not         
# Need bulk density of fine fraction g soil per cm^3 (which is the same as g soil per ml)                 
Gravel_Volume = Gravel_Weight/2.65, # Gravel Volume
Soil_Volume = (Plaster_Volume - Gravel_Volume)
Bulk_density = Total_Soil_Dry/Soil_Volume
# # gC/g soil, or g C/area
# SOC_Bulk_Density = Total_Soil_C/(Plaster_Volume - Gravel_Volume) #m^2
                                    
# We want to be able to do this: g C/m^2 and g C/g                                  
                                    
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