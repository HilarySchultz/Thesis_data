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


#### INITIAL DATA INPUT ####
# Reading in .csv of SGF, chamber field data, and HydraGO data
CO2_fluxes <- read.csv("CO2_fluxes_final.csv", header = T, sep = ",") %>%
  rename(Time = Date)
CH4_fluxes <- read.csv("CH4_fluxes_final.csv", header = T, sep = ",") %>%
  rename(Time = Date)

# Need to change the dates to the week of sampling
CO2_fluxes <- CO2_fluxes %>%
  mutate(Date = case_when(Time == "2021-06-07" ~ "6/7/2021",
                          Time == "2021-06-08" ~ "6/7/2021",
                          Time == "2021-06-10" ~ "6/7/2021",
                          Time == "2021-06-14" ~ "6/14/2021",
                          Time == "2021-06-15" ~ "6/14/2021",
                          Time == "2021-06-17" ~ "6/14/2021",
                          Time == "2021-06-28" ~ "6/28/2021",
                          Time == "2021-06-29" ~ "6/28/2021",
                          Time == "2021-07-01" ~ "6/28/2021",
                          Time == "2021-07-12" ~ "7/12/2021",
                          Time == "2021-07-13" ~ "7/12/2021",
                          Time == "2021-07-15" ~ "7/12/2021",
                          Time == "2021-07-26" ~ "7/26/2021",
                          Time == "2021-07-27" ~ "7/26/2021",
                          Time == "2021-07-29" ~ "7/26/2021", 
                          Time == "2021-08-09" ~ "8/9/2021",
                          Time == "2021-08-10" ~ "8/9/2021",
                          Time == "2021-08-12" ~ "8/9/2021",
                          Time == "2021-08-25" ~ "8/25/2021",
                          Time == "2021-08-26" ~ "8/25/2021",
                          Time == "2021-08-29" ~ "8/25/2021",
                          Time == "2021-09-07" ~ "9/7/2021",
                          Time == "2021-09-08" ~ "9/7/2021",
                          Time == "2021-09-10" ~ "9/7/2021" ))
                          

CH4_fluxes <- CH4_fluxes %>%
  mutate(Date = case_when(Time == "2021-06-07" ~ "6/7/2021",
                          Time == "2021-06-08" ~ "6/7/2021",
                          Time == "2021-06-10" ~ "6/7/2021",
                          Time == "2021-06-14" ~ "6/14/2021",
                          Time == "2021-06-15" ~ "6/14/2021",
                          Time == "2021-06-17" ~ "6/14/2021",
                          Time == "2021-06-28" ~ "6/28/2021",
                          Time == "2021-06-29" ~ "6/28/2021",
                          Time == "2021-07-01" ~ "6/28/2021",
                          Time == "2021-07-12" ~ "7/12/2021",
                          Time == "2021-07-13" ~ "7/12/2021",
                          Time == "2021-07-15" ~ "7/12/2021",
                          Time == "2021-07-26" ~ "7/26/2021",
                          Time == "2021-07-27" ~ "7/26/2021",
                          Time == "2021-07-29" ~ "7/26/2021", 
                          Time == "2021-08-09" ~ "8/9/2021",
                          Time == "2021-08-10" ~ "8/9/2021",
                          Time == "2021-08-12" ~ "8/9/2021",
                          Time == "2021-08-25" ~ "8/25/2021",
                          Time == "2021-08-26" ~ "8/25/2021",
                          Time == "2021-08-29" ~ "8/25/2021",
                          Time == "2021-09-07" ~ "9/7/2021",
                          Time == "2021-09-08" ~ "9/7/2021",
                          Time == "2021-09-10" ~ "9/7/2021"))

# Changing characters to factors
CO2_fluxes$Date <- ordered(CO2_fluxes$Date, 
                              levels = c("6/7/2021", "6/14/2021", "6/28/2021",
                                         "7/12/2021", "7/26/2021", "8/9/2021", 
                                         "8/25/2021", "9/7/2021"))
CO2_fluxes$Site <- as.factor(CO2_fluxes$Site)
CO2_fluxes$Reach <- as.factor(CO2_fluxes$Reach)
CO2_fluxes$Chamber <- ordered(CO2_fluxes$Chamber, 
                              levels = c(1:12))

CH4_fluxes$Date <- ordered(CH4_fluxes$Date, 
                          levels = c("6/7/2021", "6/14/2021", "6/28/2021",
                                     "7/12/2021", "7/26/2021", "8/9/2021", 
                                     "8/25/2021", "9/7/2021"))
CH4_fluxes$Site <- as.factor(CH4_fluxes$Site)
CH4_fluxes$Reach <- as.factor(CH4_fluxes$Reach)
CH4_fluxes$Chamber <- ordered(CH4_fluxes$Chamber, levels=c(1:12))

# None
CO2_outliers <- CO2_fluxes %>%
  rstatix::identify_outliers("CO2_flux_g") %>%
  filter(is.extreme == TRUE)

# Need to see what happens here when I pull all this data out - not done yet
CH4_outliers <- CH4_fluxes %>%
  rstatix::identify_outliers("CH4_flux_g") %>%
  filter(is.extreme == TRUE)

CH4_fluxes <- anti_join(CH4_fluxes, CH4_outliers)

# Writing out .csv for soil gas flux covariates
# These .csv will have outliers removed as well as new ordered and factored variables
# There was just going to be too much going on in this script

write.csv(CH4_fluxes, "Methane_data.csv", row.names = F)
write.csv(CO2_fluxes, "CarbonDioxide_data.csv", row.names = F)

# Covariates: soil moisture, soil temperature, and EC (SOC - in another RScript)
#### Histograms to see the distribution of the data ####
CO2_fluxes %>% 
  ggplot(aes(CO2_flux_g)) +
  geom_histogram(binwidth = 0.1) +
  facet_grid(Reach~Site)
  # The distributions across sites are roughly normal, but the numbers are small
  
CH4_fluxes %>% 
  ggplot(aes(CH4_flux_g)) +
  geom_histogram(binwidth = 0.001) +
  facet_grid(Reach~Site) +
  theme(axis.text.x = element_text(angle = 45))
  # The distributions across sites are left skewed, but the numbers are small and sometimes negative
## The issue is that you cannot log negative numbers...


#### Exploratory CO2 Models ####
cdioxidemodel <- glmer(CO2_flux_g ~ Date + Reach + (1|Site) + (1|Chamber), 
                       data = CO2_fluxes,
                       family = Gamma(link = "log"))
plot(cdioxidemodel)
Anova(cdioxidemodel)
AICc(cdioxidemodel) # 647.1857
AIC(cdioxidemodel) # 645.9574

cdioxidemodel1 <- lmer(CO2_flux_g ~ Date + Reach + (1|Site) + (1|Chamber), 
                       data = CO2_fluxes)
plot(cdioxidemodel1)
Anova(cdioxidemodel1)
AICc(cdioxidemodel1) # 727.0645
AIC(cdioxidemodel1) # 725.8361


#### Post-hoc Test ####
# Reach by Date
# For ribbon plot
cdioxide_reachdate_emm <- emmeans(cdioxidemodel, ~ Reach|Date,
                       type = "response")
cdioxide_reachdate_emm_sum <- summary(cdioxide_reachdate_emm)

##
# Reach by reach
# For boxplot
cdioxide_reach_emm <- emmeans(cdioxidemodel, ~ Reach,
                                  type = "response")

###
cdioxide_reach_cld <- cld(cdioxide_reach_emm,
                         alpha = 0.05, 
                         Letters = letters)
cdioxide_reach_cld$.group <- gsub(" ", "", cdioxide_reach_cld$.group)


#### CO2 PLOTS ####
# Ribbon plot
ggplot(data = cdioxide_reachdate_emm_sum) +
  geom_ribbon(aes(x = Date,
                  ymin = asymp.LCL, 
                  ymax = asymp.UCL,
                  group = Reach,
                  fill = Reach), 
              alpha = 0.40, 
              color = NA) + # opaqueness of the CI
  # fill = "#3984ff") +
  geom_line(aes(x = Date, 
                y = response, 
                group = Reach, 
                color = Reach), 
            lwd = 1) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  scale_color_manual(name = "Reach", labels = c("BDA", "Reference"), values = c("Blue", "Purple")) +
  scale_fill_manual(name = "Reach", labels = c("BDA", "Reference"), values = c("Blue", "Purple")) + 
  labs(title = "Riparian Soil Carbon Dioxide Fluxes", 
       x = "Date", 
       y = expression(CO[2]~Fluxes~(g~C/m^2/d))) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(colour = "black"), 
        axis.text.x = element_text(angle = 45, vjust = 0.5))

# labels = c("6/7/2021", "", "6/28/2021", "", "7/26/2021", "", "8/25/2021", ""))

## Boxplot
ggplot(data = CO2_fluxes, aes(x = Reach, y = CO2_flux_g)) +
  geom_boxplot(aes(fill = Reach)) +
  geom_point(data = cdioxide_reach_cld, aes(x = Reach, y = response), size = 1, shape = 19,
             color = "blue") +
  geom_text(data = cdioxide_reach_cld, aes(x = Reach, y = response, label= .group,
                                           vjust = -2.2, hjust = 0.5),
            size = 5, position = position_dodge(0.5), color = "black") +
  scale_fill_manual(name = "Reach", labels = c("BDA", "Reference"), values = c("#3399FF", "#CC99FF")) +
  # scale_fill_brewer(palette = "Spectral") +
  labs(title = "Riparian Soil Carbon Dioxide Fluxes", 
       x = NULL,
       y = expression(CO[2]~Fluxes~(g~C~m^-2~d^-1))) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(colour = "black", size = 12),
        # axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        legend.position = "none") +
  scale_x_discrete(labels = c("Treatment", "Reference"))

#### CH4 GLMM ####
methanemodel <- lmer(CH4_flux_g ~ Date + Reach + (1|Site) + (1|Chamber),
                           data = CH4_fluxes)
Anova(methanemodel)
AICc(methanemodel)
AIC(methanemodel)

#### Post-hoc ####
### Emmeans
# Reach by Date for ribbon plot
methane_reachdate_emm <- emmeans(methanemodel, ~ Reach|Date,
                       type = "response")
methane_reachdate_emm_sum <- summary(methane_reachdate_emm)

# Reach by reach
# Boxplot 
methane_reach_emm <- emmeans(methanemodel, ~ Reach,
                                 type = "response")

### CLD
methane_reach_cld <- cld(methane_reach_emm,
                          alpha = 0.05, 
                          Letters = letters)
methane_reach_cld$.group <- gsub(" ", "", methane_reach_cld$.group)


#### Methane Plots ####
ggplot(data = methane_reachdate_emm_sum) +
  geom_ribbon(aes(x = Date,
                  ymin = lower.CL, 
                  ymax = upper.CL,
                  group = Reach,
                  fill = Reach), 
              alpha = 0.40, 
              color = NA) + # opaqueness of the CI
  # fill = "#3984ff") +
  geom_line(aes(x = Date, 
                y = emmean, 
                group = Reach, 
                color = Reach), 
            lwd = 1) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  scale_color_manual(name = "Reach", labels = c("BDA", "Reference"), values = c("Blue", "Purple")) +
  scale_fill_manual(name = "Reach", labels = c("BDA", "Reference"), values = c("Blue", "Purple")) + 
  labs(title = "Riparian Soil Methane Fluxes", 
       x = "Date", 
       y = expression(CH[4]~Fluxes~(g~C~m^-2~d^-1))) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(colour = "black"), 
        axis.text.x = element_text(angle = 45, vjust = 0.5)) 

# Boxplot
ggplot(data = CH4_fluxes, aes(x = Reach, y = CH4_flux_g)) +
  geom_boxplot(aes(fill = Reach)) +
  geom_point(data = methane_reach_cld, aes(x = Reach, y = emmean), size = 1, shape = 19,
             color = "blue") +
  geom_text(data = methane_reach_cld, aes(x = Reach, y = emmean, label= .group,
                                           vjust = -2.2, hjust = 0.5),
            size = 5, position = position_dodge(0.5), color = "black") +
  scale_fill_manual(name = "Reach", labels = c("BDA", "Reference"), values = c("#3399FF", "#CC99FF")) +
  # scale_fill_brewer(palette = "Spectral") +
  labs(title = "Riparian Soil Methane Fluxes", 
       x = NULL,
       y = expression(CH[4]~Fluxes~(g~C~m^-2~d^-1))) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(colour = "black", size = 12),
        # axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        legend.position = "none") +
  scale_x_discrete(labels = c("Treatment", "Reference"))

### Figuring out where missing data is 
CH4_fluxes %>%
  filter(Site == "LP", Reach == "BDA")

CH4_fluxes %>%
  filter(Site == "TP", Reach == "REF")

missingdata <- CH4_fluxes %>%
  filter(Site == "FH")

length(CH4_fluxes$X)
length(CO2_fluxes$X)









