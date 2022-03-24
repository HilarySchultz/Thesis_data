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

#### BENTHIC CALCULATIONS ####

# FBPOC calculations
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

# Benthinator volumes df
Benthinator_volume <- FBPOC_data %>% 
  select(Site, Reach, Location, Replicate, Date, Benthinator_Volume)

# BPOC calculations
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

# SBPOC calculations to compare
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

# Need to change dates so that they represent the three sampling trips (June, July, August)
Benthic_data <- Benthic_data %>%
  mutate(Date = case_when(Date == "6/14/2021" ~ "June", 
                          Date == "6/15/2021" ~ "June", 
                          Date == "6/17/2021" ~ "June", 
                          Date == "7/5/2021" ~ "July",
                          Date == "7/6/2021" ~ "July", 
                          Date == "7/8/2021" ~ "July",
                          Date == "8/2/2021" ~ "August",
                          Date == "8/3/2021" ~ "August",
                          Date == "8/5/2021" ~ "August"))

write.csv(Benthic_data, "Final_benthic_df", row.names = F)

benthic_outliers <- Benthic_data %>%
  rstatix::identify_outliers("Total_BPOC_Mass_per_Area") %>%
  filter(is.extreme == TRUE)

Benthic_data <- anti_join(Benthic_data, benthic_outliers)

#### Model for Benthic Data ####

# Histogram to see the distribution of the data
# Right-skewed
Benthic_data %>% ggplot(aes(Total_BPOC_Mass_per_Area)) +
  geom_histogram(binwidth = 10) +
  facet_grid(rows = vars(Reach)) 

# Changing Site, Reach, Location, Date, and Replicate to ordered factors.
Benthic_data[,1:5] <- lapply(Benthic_data[,1:5], as.factor)

Benthic_data$Date <- ordered(Benthic_data$Date, levels = c("June", "July", "August"))
Benthic_data$Location <- ordered(Benthic_data$Location, levels = c("LWR", "MID", "UPR"))
Benthic_data$Replicate <- ordered(Benthic_data$Replicate, levels = c(1:3))

### GLMM ###
# Final model decided 3/17/22
finalbenthicmodel <- glmer(Total_BPOC_Mass_per_Area ~ Date/Site/Reach + (1|Location/Replicate),
              data = Benthic_data,
              family = Gamma(link="log"))
plot(finalbenthicmodel)
Anova(finalbenthicmodel)
AICc(finalbenthicmodel)  # 1301.195
AIC(finalbenthicmodel) # 1294.401

# Checking to see that the model is not overdispersed
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

overdisp_fun(finalbenthicmodel) # Model is not overdispersed

#### Post-hoc tests ####
### Emmeans
benthic_emm <- emmeans(finalbenthicmodel, ~ Reach|Date|Site,
                       type = "response")
benthic_emm_sum <- summary(benthic_emm)

### CLD
benthic_reach_cld <- cld(benthic_emm,
                 by = c("Site", "Date"),
                 alpha = 0.05, 
                 Letters = letters,
                 decreasing = TRUE)
benthic_reach_cld$.group = gsub(" ", "", benthic_reach_cld$.group)
benthic_reach_cld <- arrange(benthic_reach_cld, Reach, Site, Date)

# Asterisks
# reach_cld$.group <- if_else(reach_cld$.group == "b", "*","")

#### PLOTS ####
ggplot() +
  geom_boxplot(data = Benthic_data, aes(x = Reach, y = Total_BPOC_Mass_per_Area, fill = Reach)) +
  geom_point(data = benthic_reach_cld, aes(x = Reach, y = response), size = 1, shape = 19,
             color = "blue") +
  geom_text(data = benthic_reach_cld, aes(x = Reach, y = response, label= .group,
                                       vjust = -2.1, hjust = 0.5),
            size = 6, position = position_dodge(0.5), color = "black") +
  geom_text(aes()) +
  scale_fill_manual(name = "Reach", labels = c("Treatment", "Reference"), values = c("#3399FF", "#CC99FF")) +
  # scale_fill_brewer(palette = "Spectral") +
  labs(title = "Benthic Particulate Organic Carbon Pools", 
       x = NULL, 
       y = expression(BPOC~(g~C~m^-2))) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(colour = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  theme(axis.text = element_text(size = 12)) +
  facet_grid(Date~Site) 

# This code is to change the orientation of the letter display
# rep(-1.8, 5), 2, rep(-1.8, 5)


