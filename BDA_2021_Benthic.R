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
library(ggpubr)
library(MuMIn)
library(gridExtra)
# library(mrfDepth)
# library(robustbase)
library(univOutl)
library(Outliers)
# library(extremevalues)

options(scipen = 999) # Get rid of scientific notation

#### BENTHIC CALCULATIONS ####

# FBPOC calculations
# Final calculations: June 2022
Sampler_radius_cm <- 8.75 # radius in cm

Sampler_area_cm <- pi*(Sampler_radius_cm^2) # 240.5 cm^2

FBPOC_data <- read.csv("BDA_FBPOC_Calc.csv", header = T, sep = ",") %>% 
  select(2:8, 14:15) %>%
  rename(Dry_Weight_Filter_g = Dry.Weight...Filter..g., 
         Ashed_Weight_Filter_g = Ash...Filter..g., 
         Average_Depth_cm = Average.Depth..cm.,
         Filtered_Volume_ml = Filtered.Volume..ml.) %>%
  mutate(FBPOM_Sample = Dry_Weight_Filter_g - Ashed_Weight_Filter_g,
         Benthinator_Volume_ml = pi*(Sampler_radius_cm^2)*Average_Depth_cm,
         Ratio = Benthinator_Volume_ml/Filtered_Volume_ml, 
         # cm^3 = 1 ml, so now the units are in ml 
         FBPOM_Total = Ratio * FBPOM_Sample,
         FBPOC_Carbon_Mass_per_Area = FBPOM_Total/(Sampler_area_cm/100) * 0.52)
         # 52 percent C, need to change cm to m so divide the sampler area by 100
  # Units are in g/m^2

# BPOC calculations
# Final calculations: June 2022
BPOC_data <- read.csv("BDA_CBPOC_Calc.csv", header = T, sep = ",") %>%
  select(3:12) %>%
  select(-Subsample.Dry.Weight..g.) %>%
  rename(BPOC_Dry_Weight = Dry.Weight..g., 
         BPOC_Tin_Weight = Tin.Weight..g., # renaming columns
         BPOC_Subsample_Dry_Tin = Subsample.Dry.Weight...Tin..g., 
         BPOC_Ash_Tin = Tin...Ash..g.) %>%
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
  # meaning that the inorganic material is subtracted from the dry weight to yield the OM content
  mutate(BPOC_Sample= ifelse(BPOC_Subsample_Dry_Tin != 0, 
                             BPOC_Subsample_OM*(BPOC_Dry_Weight/BPOC_Subsample_Dry_Weight),
                             BPOC_Dry_Weight - BPOC_Ash)) %>%
  # Determining the sample OM mass based on a sample/subsample ratio
  # If there is no subsample, then the sample ash weight is subtracted from the sample dry weight to yield OM content
  mutate(BPOC_Carbon_Mass_per_Area = BPOC_Sample/(Sampler_area_cm/100) * 0.52) 
# Units here are (g C/m^2)

# Need the filtered volume for the suspended CBPOC samples (SBPOC)
Benthic_field_data <- read.csv("Benthic_Field_Data.csv", header = T, sep = ",") %>%
  select(2:6, 11) %>%
  rename(Filtered_suspended_volume_ml = Filtered.BCPOC..ml.)

# Benthinator volumes df
Benthinator_volume <- FBPOC_data %>% 
  select(Site, Reach, Location, Replicate, Date, Benthinator_Volume_ml)

volumes <- full_join(Benthic_field_data, Benthinator_volume)

# SBPOC calculations to compare
# Final calculations: June 2022
SBPOC <- read.csv("BDA_SCBPOC_Calc.csv", header = T, sep = ",") %>%
  select(3:12) %>%
  select(-Subsample.Dry.Weight..g.)

SBPOC_data <- full_join(SBPOC, volumes) %>% # Joining the benthinator and sample volumes in with the df
  rename(SBPOC_Dry_Weight = Dry.Weight..g., 
         SBPOC_Tin_Weight = Tin.Weight..g., 
         SBPOC_Subsample_Dry_Tin = Subsample.Dry.Weight...Tin..g., 
         SBPOC_Ash_Tin = Tin...Ash..g.) %>%
  mutate(SBPOC_Subsample_Dry_Weight = ifelse(SBPOC_Subsample_Dry_Tin > 0, SBPOC_Subsample_Dry_Tin - SBPOC_Tin_Weight, 0),
         SBPOC_Ash = SBPOC_Ash_Tin - SBPOC_Tin_Weight, 
         SBPOC_Subsample_OM = ifelse(SBPOC_Subsample_Dry_Weight > 0, SBPOC_Subsample_Dry_Weight - SBPOC_Ash, 0), 
         SBPOC_Sample_OM = ifelse(SBPOC_Subsample_Dry_Tin != 0, 
                                  SBPOC_Subsample_OM * (SBPOC_Dry_Weight/SBPOC_Subsample_Dry_Weight),
                                  SBPOC_Dry_Weight - SBPOC_Ash)) %>%
  mutate(Ratio = Benthinator_Volume_ml/(Filtered_suspended_volume_ml), 
         SBPOC_Carbon_Mass_per_Area = SBPOC_Sample_OM/(Sampler_area_cm/100) * 0.52)
# Units are in g C/m^2  

# Binding benthic CBPOC, suspended CBPOC, and FBPOC into one df then selecting rows of interest
CBPOC_data <- full_join(BPOC_data, SBPOC_data) %>%
  select(Date, Site, Reach, Location, Replicate, SBPOC_Carbon_Mass_per_Area, BPOC_Carbon_Mass_per_Area)

FBPOC_join_data <- FBPOC_data %>%
  select(Date, Site, Reach, Location, Replicate, FBPOC_Carbon_Mass_per_Area)

Benthic_data <- full_join(CBPOC_data, FBPOC_join_data) %>%
  mutate(Total_BPOC_Mass_per_Area = BPOC_Carbon_Mass_per_Area + SBPOC_Carbon_Mass_per_Area + FBPOC_Carbon_Mass_per_Area) %>%
  rename(Segment = Reach, 
         Reach = Location,
         Location = Replicate)

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

# write.csv(Benthic_data, "Final_benthic_df", row.names = F)

##### Identifying outliers ####
qqPlot(Benthic_data$Total_BPOC_Mass_per_Area)

Benthic_data <- Benthic_data %>%
  slice(-127, -92)

#### Model for Benthic Data ####

# Histogram to see the distribution of the data
# Right-skewed
Benthic_data %>% ggplot(aes(Total_BPOC_Mass_per_Area)) +
  geom_histogram(binwidth = 1) +
  facet_grid(rows = vars(Segment)) 

# Changing Site, Reach, Location, Date, and Replicate to ordered factors.
Benthic_data$Date <- ordered(Benthic_data$Date, levels = c("June", "July", "August"))
Benthic_data$Reach <- ordered(Benthic_data$Reach, levels = c("UPR", "MID", "LWR"))
Benthic_data$Segment <- as.factor(Benthic_data$Segment)
Benthic_data$Site <- as.factor(Benthic_data$Site)
Benthic_data$Location <- ordered(Benthic_data$Location, levels = c("1","2","3"))

write.csv(Benthic_data, "Benthic_data.csv", row.names = FALSE)

### GLMM ###
# Final model decided 3/17/22 - Reworked data 6/22
finalbenthicmodel <- glmer(Total_BPOC_Mass_per_Area ~ Date/Site/Segment + (1|Reach) + (1|Location),
              data = Benthic_data,
              control = glmerControl(optimizer = "bobyqa",
                       optCtrl = list(maxfun = 100000)),
              family = Gamma(link = "log"))

qqnorm(resid(finalbenthicmodel))
qqline(resid(finalbenthicmodel))

# We want to look at the trigamma R2c - conditional for the entire model
benthic_fit <- r.squaredGLMM(finalbenthicmodel) 

#### Post-hoc tests ####
### Emmeans
benthic_emm <- emmeans(finalbenthicmodel, pairwise ~ Segment|Date|Site,
                       type = "response",
                       nesting = "Date %in% Site, Segment %in% (Site*Date)")
benthic_emm_sum <- summary(benthic_emm)
benthic_emmeans <- summary(benthic_emm$emmeans)
benthic_emm_contrast <- summary(benthic_emm$contrasts)

### CLD
benthic_cld <- cld(benthic_emm,
                 by = c("Site", "Date"),
                 alpha = 0.05, 
                 Letters = letters,
                 decreasing = TRUE)
benthic_cld$.group = gsub(" ", "", benthic_cld$.group)
benthic_cld <- arrange(benthic_cld, Date, Site, Segment)

# Asterisks
benthic_cld$.group <- if_else(benthic_cld$.group == "b", "*","")

#### Mean-error-plots for Site ####
ggplot() +
  geom_jitter(Benthic_data, 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6),
              mapping = aes(Date, Total_BPOC_Mass_per_Area, color = Segment)) + 
              # colour = "black", 
              # width = 0.1) +
  geom_point(data = benthic_emmeans, 
             mapping = aes(Date, response, fill = Segment), 
             position = position_dodge(0.6),
             # colour = "black",
             shape = 1, 
             size = 3) +
  geom_text(data = benthic_cld, aes(x = Date, y = response, label= .group,
                                    vjust = -1, hjust = 0.5),
            size = 6, position = position_dodge(0.5), color = "black") +
  scale_y_log10(limits = c(0.001,1e2)) +
  scale_color_manual(name = "Segment", 
                     labels = c("Treatment", "Reference"), 
                     values = c("#0072B2", "#009E73")) +
  scale_fill_manual(name = "Segment", 
                    labels = c("Treatment", "Reference"), 
                    values = c("#0072B2", "#009E73")) +
  geom_errorbar(data = benthic_emmeans, width = 0.2, position = position_dodge(0.6),
                mapping = aes(Date, response, ymin = asymp.LCL, ymax = asymp.UCL, fill = Segment)) +
  # scale_color_manual(values = c())
  labs(title = NULL, 
       x = NULL, 
       y = expression(Benthic~Particulate~Organic~Carbon~(g~m^-2))) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(colour = "black", size = 12),
        panel.spacing.x = unit(1, "lines"),
        axis.title = element_text(size = 12),
        axis.title.y = element_text(size = 14)) +
  facet_grid(rows = vars(Site))

#### Stats for benthic model by site ####
benthic_date_sum <- benthic_emmeans %>%
  group_by(Date, Site) %>%
  summarise(Avg = mean(response),
            mean_UCL = mean(asymp.UCL),
            mean_LCL = mean(asymp.LCL))

benthic_reach_sum <- benthic_emmeans %>%
  group_by(Site, Segment) %>%
  summarise(Avg = mean(response),
            mean_UCL = mean(asymp.UCL),
            mean_LCL = mean(asymp.LCL))

benthic_site_sum <- benthic_emmeans %>%
  group_by(Site) %>%
  summarise(Avg = mean(response),
            mean_UCL = mean(asymp.UCL),
            mean_LCL = mean(asymp.LCL))

date_difference <- benthic_emmeans %>%
  group_by(Date, Site, Segment) %>%
  summarise(Avg = mean(response),
            mean_UCL = mean(asymp.UCL),
            mean_LCL = mean(asymp.LCL))

#### SEGMENT BENTHIC CALCULATIONS ####
Wetted_widths <- read.csv("wetted widths_raw.csv", header = T, sep = ",") 

Average_widths <- Wetted_widths %>%
  group_by(Site, Segment, Reach, Date) %>%
  summarise(avg_width = mean(Width)) %>%
  na.exclude()

Average_widths$Date <- ordered(Average_widths$Date, levels = c("June", "July", "August"))

Benthic_data <- full_join(Benthic_data, Average_widths) %>%
  mutate(BPOC_by_100m = Total_BPOC_Mass_per_Area*avg_width*100)
Benthic_data$Site <- as.factor(Benthic_data$Site)
Benthic_data$Segment <- as.factor(Benthic_data$Segment)
Benthic_data$Reach <- as.factor(Benthic_data$Reach)

benthicmodel <- glmer(BPOC_by_100m ~ Site/Segment + (1|Reach) + (1|Location),
                           data = Benthic_data,
                           control = glmerControl(optimizer = "bobyqa",
                           optCtrl = list(maxfun = 100000)),
                           family = Gamma(link = "log"))

plot(benthicmodel)
Anova(benthicmodel)

### Emmeans
BPOC_emm <- emmeans(benthicmodel, pairwise ~ Segment|Site,
                    type = "response")
BPOC_emm_sum <- summary(BPOC_emm)
BPOC_emmeans <- summary(BPOC_emm$emmeans)
BPOC_emm_contrast <- summary(BPOC_emm$contrasts)

### CLD
BPOC_cld <- cld(BPOC_emm,
                   by = "Site",
                   alpha = 0.05, 
                   Letters = letters,
                   decreasing = TRUE)
BPOC_cld$.group = gsub(" ", "", BPOC_cld$.group)
BPOC_cld <- arrange(BPOC_cld, Site, Segment)

# Asterisks
BPOC_cld$.group <- if_else(BPOC_cld$.group == "b", "*","")

benthic_segment_plot <- ggplot() +
  geom_jitter(Benthic_data, 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6),
              mapping = aes(Segment, BPOC_by_100m, color = Segment)) + 
  # colour = "black", 
  # width = 0.1) +
  geom_point(data = BPOC_emmeans, 
             mapping = aes(Segment, response, fill = Segment), 
             position = position_dodge(0.6),
             # colour = "black",
             shape = 1, 
             size = 3) +
  geom_text(data = BPOC_cld, aes(x = Segment, y = response, label= .group,
                                    vjust = -1, hjust = -15),
            size = 6, color = "black") +
  scale_y_log10(limits = c(1,1e4)) +
  scale_color_manual(name = "Segment", 
                     labels = c("Treatment", "Reference"), 
                     values = c("#0072B2", "#009E73")) +
  scale_fill_manual(name = "Segment", 
                    labels = c("Treatment", "Reference"), 
                    values = c("#0072B2", "#009E73")) +
  geom_errorbar(data = BPOC_emmeans, width = 0.2, position = position_dodge(0.6),
                mapping = aes(Segment, response, ymin = asymp.LCL, ymax = asymp.UCL, fill = Segment)) +
  # scale_color_manual(values = c())
  labs(title = NULL, 
       x = NULL, 
       y = expression(Benthic~Particulate~Organic~Carbon~per~100~m~(g))) +
  theme_bw() +
  scale_x_discrete(labels = c("BDA" = "Treatment", "REF" = "Reference")) +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(colour = "black", size = 12),
        panel.spacing.x = unit(1, "lines"),
        axis.title = element_text(size = 12),
        axis.title.y = element_text(size = 14)) +
  facet_grid(rows = vars(Site))

save_plot("benthic_segment_plot.png", benthic_segment_plot, base_height = 6, base_width = 9)

ggarrange(benthic_plot +
            theme(axis.ticks.y = element_blank(),
                  plot.margin = margin(r = 1) ),
          benthic_m_plot +
            theme(
                  plot.margin = margin(r = 1, l = 1) ),
          nrow = 1, ncol = 2,
          common.legend = TRUE,
          legend = "bottom")


### Data for Locations in Treatment #### 
benthic_BDA <- Benthic_data %>%
  filter(Segment == "BDA")

benthic_reach_model <- glmer(Total_BPOC_Mass_per_Area ~ Site/Reach + (1|Location),
                        data = benthic_BDA,
                           # control = glmerControl(optimizer = "bobyqa",
                           # optCtrl = list(maxfun = 100000)),
                           family = Gamma(link = "log"))
plot(benthic_reach_model)
Anova(benthic_reach_model)

benthic_reach_emm <- emmeans(benthic_reach_model, pairwise ~ Reach|Site,
                       type = "response")
                       # nesting = "Site %in% Reach")
benthic_reach_emm_sum <- summary(benthic_reach_emm)
benthic_reach_emmeans <- summary(benthic_reach_emm$emmeans)
benthic_reach_emm_contrast <- summary(benthic_reach_emm$contrasts)

benthic_reach_cld <- cld(benthic_reach_emm,
                         # by = c("Site", "Reach"),
                         alpha = 0.05,
                         Letters = letters,
                         decreasing = F)
benthic_reach_cld$.group = gsub(" ", "", benthic_reach_cld$.group)
benthic_reach_cld <- arrange(benthic_reach_cld, Site, Reach)

#### Mean-error-plots for Reach ####
ggplot() +
  geom_jitter(benthic_BDA, 
              mapping = aes(Reach, Total_BPOC_Mass_per_Area), 
              colour = "black", 
              width = 0.1) +
  geom_point(data = benthic_reach_emmeans, 
             mapping = aes(Reach, response), 
             shape = 1, 
             size = 4) +
  geom_errorbar(data = benthic_reach_emmeans, width = 0.2,
                mapping = aes(Reach, response, ymin = asymp.LCL, ymax = asymp.UCL)) +
  scale_y_log10(limits = c(1,1e3)) +
  labs(title = NULL, 
       x = "Reach", 
       y = expression(Benthic~Particulate~Organic~Carbon~(g~m^-2))) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(colour = "black", size = 12),
        panel.spacing.x = unit(1, "lines"),
        axis.title = element_text(size = 12),
        axis.title.y = element_text(size = 14)) +
  facet_grid(~Site)

#### Reach Stats ####
b_reach <- benthic_reach_emmeans %>%
  group_by(Reach) %>%
  summarise(Avg = mean(response),
            mean_UCL = mean(asymp.UCL),
            mean_LCL = mean(asymp.LCL))

#### Location stats ####

benthic_location_model <- glmer(Total_BPOC_Mass_per_Area ~ Site/Location + (1|Reach),
                       data = benthic_BDA,
                       control = glmerControl(optimizer = "bobyqa",
                                              optCtrl = list(maxfun = 100000)),
                       family = Gamma(link = "log"))
plot(benthic_location_model)
Anova(benthic_location_model)
b_loc_emm <- emmeans(benthic_location_model, pairwise ~ Location|Site,
                          type = "response")

b_loc_emm_sum <- summary(b_loc_emm)
b_loc_emmeans <- summary(b_loc_emm$emmeans)
b_loc_emm_contrast <- summary(b_loc_emm$contrasts)

b_loc_cld <- cld(b_loc_emm,
                         by = c("Site"),
                         alpha = 0.05, 
                         Letters = letters,
                         decreasing = TRUE)
b_loc_cld$.group = gsub(" ", "", b_loc_cld$.group)
b_loc_cld <- arrange(b_loc_cld, Site, Location)

b_reach_cld <- b_reach_cld %>%
  mutate(.group = case_when(.group == "b" ~ "*",
                            .group == "c" ~ "*",
                            .group == "a" ~ ""))

benthic_BDA$Reach <- ordered(benthic_BDA$Reach, levels = c("UPR", "MID", "LWR"))

ggplot() +
  geom_boxplot(data = benthic_BDA, aes(x = , y = Total_BPOC_Mass_per_Area, fill = Reach)) +
  # geom_point(data = b_reach_cld, aes(x = Location, y = response), size = 3, shape = 2,
  #            color = "black") +
  geom_text(data = benthic_reach_cld, aes(x = Reach, y = response, label= .group,
                                          vjust = c(-5,-9,-12, -5, -8, -5, rep(-5,3)), hjust = 0.5),
            size = 6, color = "black", position = position_dodge(0.5)) +
  scale_fill_manual(name = "Reach", labels = c("UPR", "MID", "LWR"), values = c("#edf8b1", "#7fcdbb", "#2c7fb8")) +
  # scale_y_log10() +
  labs(title = NULL, 
       x = "Reach", 
       y = expression(Benthic~Particulate~Carbon~Pools~(g~C~m^-2))) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(colour = "black", size = 12),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 0.5),
        axis.title.y = element_text(size = 14)) +
  facet_grid(~Site) 

bsum <- b_reach_emmeans %>%
  group_by(Location) %>%
  summarise(Avg = mean(response),
            mean_UCL = mean(asymp.UCL),
            mean_LCL = mean(asymp.LCL))


#### SEGMENT BENTHIC CALCULATIONS ####
Wetted_widths <- read.csv("wetted widths_raw.csv", header = T, sep = ",") 

LP_treat <- 60
LP_ref <- 18

FH_treat <- 70
FH_ref <- 20

TP_treat <- 170
TP_ref <- 160

benthic_means <- Benthic_data %>%
  group_by(Site, Segment, Date)


Average_widths <- Wetted_widths %>%
  group_by(Site, Segment, Reach, Date) %>%
  summarise(avg_width = mean(Width)) %>%
  na.exclude()

Area <- Average_widths %>%
  mutate(Segment_length = case_when(Site == "FH" & Segment == "BDA" ~ FH_treat,
                                    Site == "FH" & Segment == "REF" ~ FH_ref,
                                    Site == "LP" & Segment == "BDA" ~ LP_treat,
                                    Site == "LP" & Segment == "REF" ~ LP_ref,
                                    Site == "TP" & Segment == "BDA" ~ TP_treat,
                                    Site == "TP" & Segment == "REF" ~ TP_ref))

Benthic_data <- full_join(Benthic_data, Area) %>%
  mutate(streambed_area = avg_width * Segment_length) %>%
  mutate(benthic_per_m = avg_width * 100 * Total_BPOC_Mass_per_Area)


# Asterisks
benthic_reach_cld$.group <- if_else(benthic_reach_cld$.group == "b", "*","")

# FH_loc <- benthic_means %>%
#   filter(Site == "FH") %>%
#   ggplot(aes(x = Location, y = Total_BPOC_Mass_per_Area)) +
#   geom_boxplot(aes(fill = Replicate), position=position_dodge(.9)) +
#   scale_fill_manual(values = c("#edf8b1", "#7fcdbb", "#2c7fb8")) +
#   stat_summary(fun = mean, aes(group = Replicate), 
#                geom = "point", shape = 20, size = 3, color = "black", fill = "black",
#                position=position_dodge(.9)) +
#   theme_bw() + 
#   theme(plot.title = element_text(hjust = 0.5), 
#         axis.text = element_text(colour = "black", size = 12),
#         axis.title.y = element_text(size = 14)) +
#   labs(title = "Fish Creek", 
#        x = NULL, 
#        y = expression(Benthic~Particulate~Organic~Carbon~(g~m^-2))) +
#   ylim(0,325)
# 
# LP_loc <- benthic_means %>%
#   filter(Site == "LP") %>%
#   ggplot(aes(x = Location, y = Total_BPOC_Mass_per_Area)) +
#   geom_boxplot(aes(fill = Replicate), position=position_dodge(.9)) +
#   scale_fill_manual(values = c("#edf8b1", "#7fcdbb", "#2c7fb8")) +
#   stat_summary(fun = mean, aes(group = Replicate), 
#                geom = "point", shape = 20, size = 3, color = "black", fill = "black",
#                position=position_dodge(.9)) +
#   theme_bw() + 
#   theme(plot.title = element_text(hjust = 0.5), 
#         axis.text = element_text(colour = "black", size = 12),
#         axis.title.y = element_text(size = 14)) +
#   labs(title = "Lost Prairie Creek", 
#        x = NULL, 
#        y = expression(Benthic~Particulate~Organic~Carbon~(g~m^-2))) +
#   ylim(0,325)
# 
# TP_loc <-benthic_means %>%
#   filter(Site == "TP") %>%
#   ggplot(aes(x = Location, y = Total_BPOC_Mass_per_Area)) +
#   geom_boxplot(aes(fill = Replicate), position=position_dodge(.9)) +
#   scale_fill_manual(values = c("#edf8b1", "#7fcdbb", "#2c7fb8")) +
#   stat_summary(fun = mean, aes(group = Replicate), 
#                geom = "point", shape = 20, size = 3, color = "black", fill = "black",
#                position=position_dodge(.9)) +
#   theme_bw() + 
#   theme(plot.title = element_text(hjust = 0.5), 
#         axis.text = element_text(colour = "black", size = 12),
#         axis.title.y = element_text(size = 14)) +
#   labs(title = "Teepee Creek", 
#        x = NULL, 
#        y = expression(Benthic~Particulate~Organic~Carbon~(g~m^-2))) +
#   ylim(0,325) 
#   
# ggarrange(FH_loc +
#             theme(axis.ticks.y = element_blank(),
#                   plot.margin = margin(r = 1) ), 
#           LP_loc + 
#             theme(axis.text.y = element_blank(),
#                   axis.ticks.y = element_blank(),
#                   axis.title.y = element_blank(),
#                   plot.margin = margin(r = 1, l = 1) ), 
#           TP_loc + 
#             theme(axis.text.y = element_blank(),
#                   axis.ticks.y = element_blank(),
#                   axis.title.y = element_blank(),
#                   plot.margin = margin(l = 1)  ),
#           nrow = 1, ncol = 3, 
#           common.legend = TRUE, 
#           legend = "bottom")
# 
# rep_avg <- benthic_means %>%
#   select(-Reach) %>%
#   group_by(Site, Location) %>%
#   summarise(mean_rep = mean(Total_BPOC_Mass_per_Area))
# 
# benthic_means %>%
#   ggplot(aes(x = Location, y = Total_BPOC_Mass_per_Area)) +
#   geom_boxplot(aes(fill = Replicate), position=position_dodge(.9)) +
#   scale_fill_manual(values = c("#edf8b1", "#7fcdbb", "#2c7fb8")) +
#   stat_summary(fun = mean, aes(group = Replicate),
#                geom = "point", shape = 20, size = 3, color = "black", fill = "black",
#                position=position_dodge(.9)) +
#   theme_bw() + 
#   theme(plot.title = element_text(hjust = 0.5), 
#         axis.text = element_text(colour = "black", size = 12),
#         axis.title.y = element_text(size = 14)) +
#   labs(title = "Carbon Pools by Location", 
#        x = NULL, 
#        y = expression(Benthic~Particulate~Organic~Carbon~(g~m^-2))) +
#   ylim(0,325) + 
#   facet_grid(~Site)
# 
# 
# benthic_means %>%
#   ggplot(aes(x = Location, y = Total_BPOC_Mass_per_Area)) +
#   geom_point() +
#   # scale_fill_manual(values = c("#edf8b1", "#7fcdbb", "#2c7fb8")) +
#   # stat_summary(fun = mean, aes(group = Replicate), 
#   #              geom = "point", shape = 20, size = 3, color = "black", fill = "black",
#   #              position=position_dodge(.9)) +
#   theme_bw() + 
#   theme(plot.title = element_text(hjust = 0.5), 
#         axis.text = element_text(colour = "black", size = 12),
#         axis.title.y = element_text(size = 14)) +
#   labs(title = "Carbon Pools by Location", 
#        x = NULL, 
#        y = expression(Benthic~Particulate~Organic~Carbon~(g~m^-2))) +
#   ylim(0,325) + 
#   facet_grid(Site~Replicate)
# 

#### Data for proportions ####
Benthic_sum <- Benthic_data %>%
  mutate(Coarse_POC = SBPOC_Carbon_Mass_per_Area + BPOC_Carbon_Mass_per_Area) %>%
  group_by(Date, Site, Reach) %>%
  summarise(mean_CBPOC = mean(Coarse_POC), mean_FBPOC = mean(FBPOC_Carbon_Mass_per_Area)) %>%
  mutate(Total_mean_BPOC = mean_CBPOC + mean_FBPOC) %>%
  gather(mean_CBPOC, mean_FBPOC, key = "sample", value = "value")

Benthic_sum$Date <- ordered(Benthic_sum$Date, levels = c("June", "July", "August"))
levels(Benthic_sum$Reach) <- c("Treatment", "Reference")

FH_bar <- Benthic_sum %>%
  filter(Site == "FH") %>%
  ggplot(aes(Date, value, fill = sample)) +
  geom_bar(position="dodge", stat="identity") +
  facet_grid(rows = vars(Reach)) +
  geom_point(aes(y = Total_mean_BPOC, color = "Total BPOC")) +
  labs(title = "Fish Creek", 
       x = NULL, 
       y = expression(Benthic~Particulate~Organic~Carbon~(g~C~m^-2))) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(colour = "black")) +
  theme(axis.text = element_text(size = 12)) +
  scale_fill_manual(name = "Particulate \nCarbon", labels = c("Coarse", "Fine"), 
                    values = c("#F8766D", "#00BFC4")) +
  scale_color_manual(name = NULL, values = c("Total BPOC" = "black")) +
  guides(fill = guide_legend(override.aes = list(shape = NA))) +
  ylim(0,80)


LP_bar <- Benthic_sum %>%
  filter(Site == "LP") %>%
  ggplot(aes(Date, value, fill = sample)) +
  geom_bar(position="dodge", stat="identity") +
  facet_grid(rows = vars(Reach)) +
  geom_point(aes(y = Total_mean_BPOC, color = "Total BPOC")) +
  labs(title = "Lost Prairie Creek", 
       x = NULL, 
       y = expression(BPOC~(g~C~m^-2))) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(colour = "black")) +
  theme(axis.text = element_text(size = 12)) +
  scale_fill_manual(name = "Particulate \nCarbon", labels = c("Coarse", "Fine"), 
                    values = c("#F8766D", "#00BFC4")) +
  scale_color_manual(name = NULL, values = c("Total BPOC" = "black")) +
  guides(fill = guide_legend(override.aes = list(shape = NA))) +
  ylim(0,80)


TP_bar <- Benthic_sum %>%
  filter(Site == "TP") %>%
  ggplot(aes(Date, value, fill = sample)) +
  geom_bar(position="dodge", stat="identity") +
  facet_grid(rows = vars(Reach)) +
  geom_point(aes(y = Total_mean_BPOC, color = "Total BPOC")) +
  labs(title = "Teepee Creek", 
       x = NULL, 
       y = expression(BPOC~(g~C~m^-2))) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(colour = "black")) +
  theme(axis.text = element_text(size = 12)) +
  scale_fill_manual(name = "Particulate \nCarbon", labels = c("Coarse", "Fine"), 
                    values = c("#F8766D", "#00BFC4")) +
  scale_color_manual(name = NULL, values = c("Total BPOC" = "black")) +
  guides(fill = guide_legend(override.aes = list(shape = NA))) +
  ylim(0,80)


ggarrange(FH_bar +
            theme(axis.ticks.y = element_blank(),
                  plot.margin = margin(r = 1) ), 
          LP_bar + 
            theme(axis.text.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  axis.title.y = element_blank(),
                  plot.margin = margin(r = 1, l = 1) ), 
          TP_bar + 
            theme(axis.text.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  axis.title.y = element_blank(),
                  plot.margin = margin(l = 1)  ),
          nrow = 1, ncol = 3, 
          common.legend = TRUE, 
          legend = "bottom")

# ggplot(Benthic_data, aes(fill=c(Coarse_POC, FBPOC_Carbon_Mass_per_Area), y=value, x=specie)) + 
#   geom_bar(position="dodge", stat="identity")
