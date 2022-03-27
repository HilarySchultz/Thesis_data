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
  mutate(Date = case_when(Time == "2021-06-07" ~ "06/07/2021",
                          Time == "2021-06-08" ~ "06/07/2021",
                          Time == "2021-06-10" ~ "06/07/2021",
                          Time == "2021-06-14" ~ "06/14/2021",
                          Time == "2021-06-15" ~ "06/14/2021",
                          Time == "2021-06-17" ~ "06/14/2021",
                          Time == "2021-06-28" ~ "06/28/2021",
                          Time == "2021-06-29" ~ "06/28/2021",
                          Time == "2021-07-01" ~ "06/28/2021",
                          Time == "2021-07-12" ~ "07/12/2021",
                          Time == "2021-07-13" ~ "07/12/2021",
                          Time == "2021-07-15" ~ "07/12/2021",
                          Time == "2021-07-26" ~ "07/26/2021",
                          Time == "2021-07-27" ~ "07/26/2021",
                          Time == "2021-07-29" ~ "07/26/2021", 
                          Time == "2021-08-09" ~ "08/09/2021",
                          Time == "2021-08-10" ~ "08/09/2021",
                          Time == "2021-08-12" ~ "08/09/2021",
                          Time == "2021-08-25" ~ "08/25/2021",
                          Time == "2021-08-26" ~ "08/25/2021",
                          Time == "2021-08-29" ~ "08/25/2021",
                          Time == "2021-09-07" ~ "09/07/2021",
                          Time == "2021-09-08" ~ "09/07/2021",
                          Time == "2021-09-10" ~ "09/07/2021" ))


CH4_fluxes <- CH4_fluxes %>%
  mutate(Date = case_when(Time == "2021-06-07" ~ "06/07/2021",
                          Time == "2021-06-08" ~ "06/07/2021",
                          Time == "2021-06-10" ~ "06/07/2021",
                          Time == "2021-06-14" ~ "06/14/2021",
                          Time == "2021-06-15" ~ "06/14/2021",
                          Time == "2021-06-17" ~ "06/14/2021",
                          Time == "2021-06-28" ~ "06/28/2021",
                          Time == "2021-06-29" ~ "06/28/2021",
                          Time == "2021-07-01" ~ "06/28/2021",
                          Time == "2021-07-12" ~ "07/12/2021",
                          Time == "2021-07-13" ~ "07/12/2021",
                          Time == "2021-07-15" ~ "07/12/2021",
                          Time == "2021-07-26" ~ "07/26/2021",
                          Time == "2021-07-27" ~ "07/26/2021",
                          Time == "2021-07-29" ~ "07/26/2021", 
                          Time == "2021-08-09" ~ "08/09/2021",
                          Time == "2021-08-10" ~ "08/09/2021",
                          Time == "2021-08-12" ~ "08/09/2021",
                          Time == "2021-08-25" ~ "08/25/2021",
                          Time == "2021-08-26" ~ "08/25/2021",
                          Time == "2021-08-29" ~ "08/25/2021",
                          Time == "2021-09-07" ~ "09/07/2021",
                          Time == "2021-09-08" ~ "09/07/2021",
                          Time == "2021-09-10" ~ "09/07/2021"))

# Changing characters to factors
# CO2_fluxes$Date <- ordered(CO2_fluxes$Date, 
#                               levels = c("6/7/2021", "6/14/2021", "6/28/2021",
#                                          "7/12/2021", "7/26/2021", "8/9/2021", 
#                                          "8/25/2021", "9/7/2021"))
CO2_fluxes$Site <- as.factor(CO2_fluxes$Site)
CO2_fluxes$Reach <- as.factor(CO2_fluxes$Reach)
CO2_fluxes$Chamber <- as.factor(CO2_fluxes$Chamber)

# CH4_fluxes$Date <- ordered(CH4_fluxes$Date, 
#                           levels = c("6/7/2021", "6/14/2021", "6/28/2021",
#                                      "7/12/2021", "7/26/2021", "8/9/2021", 
#                                      "8/25/2021", "9/7/2021"))
CH4_fluxes$Site <- as.factor(CH4_fluxes$Site)
CH4_fluxes$Reach <- as.factor(CH4_fluxes$Reach)
CH4_fluxes$Chamber <- as.factor(CH4_fluxes$Chamber)

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


#### CO2 Model ####
finalcdioxidemodel <- glmer(CO2_flux_g ~ Date/Site/Reach + (1|Chamber), 
                            data = CO2_fluxes,
                            control = glmerControl(optimizer = "bobyqa",
                                                   optCtrl = list(maxfun = 100000)),
                            family = Gamma(link = "log"))
plot(finalcdioxidemodel)
Anova(finalcdioxidemodel)
AICc(finalcdioxidemodel) # 592.5527
AIC(finalcdioxidemodel) # 568.9416
summary(finalcdioxidemodel)



overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}
overdisp_fun(finalcdioxidemodel) # No overdispersion

#### Post-hoc Test ####
### Emmeans
cdioxide_emm <- emmeans(finalcdioxidemodel, ~ Reach|Date|Site,
                        type = "response")
cdioxide_emm_sum <- summary(cdioxide_emm)
co2_pvalues <- pairs(cdioxide_emm) %>%
  summary()

### CLD 
co2_cld <- cld(cdioxide_emm,
               by = c("Site", "Date"),
               alpha = 0.05, 
               Letters = letters,
               decreasing = TRUE)
co2_cld$.group = gsub(" ", "", co2_cld$.group)
co2_cld <- arrange(co2_cld, Date, Site, Reach)



#### CO2 PLOTS ####
# Ribbon plot
ggplot(data = cdioxide_emm_sum) +
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
  scale_x_discrete(expand = c(0,0), guide = guide_axis(angle = 45)) +
  scale_color_manual(name = "Reach", labels = c("Treatment", "Reference"), values = c("Blue", "Purple")) +
  scale_fill_manual(name = "Reach", labels = c("Treatment", "Reference"), values = c("Blue", "Purple")) + 
  labs(title = "Riparian Soil Carbon Dioxide Fluxes", 
       x = "Date", 
       y = expression(CO[2]~Fluxes~(g~C~m^-2~d^-1))) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(colour = "black"),
        panel.spacing.x = unit(1, "lines")) +
  # axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5)) +
  facet_grid(rows = vars(Site))


co2sum <- CO2_fluxes %>%
  filter(Site == "FH") %>%
  group_by(Date, Reach) %>%
  summarise(Avg = mean(CO2_flux_g))
# 7/26/2021:FH:REF 5.564029
# 8/9/2021:FH:REF 3.899408
# 8/25/2021:FH:REF 2.654296
# 9/7/2021:FH:REF 3.227161


co2sum1 <- CO2_fluxes %>%
  group_by(Site, Reach) %>%
  summarise(Avg = mean(CO2_flux_g))
# FH    BDA    1.13
# FH    REF    2.03
## 1.79646
# LP    BDA    1.24
# LP    REF    1.07
## 1.158879
# TP    BDA    2.55
# TP    REF    2.80
## 1.098039

co2sum2 <- CO2_fluxes %>%
  group_by(Reach) %>%
  summarise(Avg = mean(CO2_flux_g))
# BDA    1.68
# REF    1.96
## 1.166667

cdioxiderange <- CO2_fluxes %>%
  group_by(Site) %>%
  summarise(Range = max(CO2_flux_g)- min(CO2_flux_g))
# Range values for each Site across Date and Reach 
# FH     5.58
# LP     2.73
# TP     6.36

cdioxiderange1 <- CO2_fluxes %>%
  group_by(Reach) %>%
  summarise(Range = max(CO2_flux_g)- min(CO2_flux_g))
# Range values across reaches
# BDA    6.30
# REF    5.91

cdioxiderange2 <- CO2_fluxes %>%
  summarise(Range = max(CO2_flux_g)- min(CO2_flux_g))
# Range value across Date, Site, and Reach
# 6.497249

# labels = c("6/7/2021", "", "6/28/2021", "", "7/26/2021", "", "8/25/2021", ""))

#### CH4 GLMM ####
finalmethanemodel <- lmer(CH4_flux_g ~ Date/Site/Reach + (1|Chamber),
                          data = CH4_fluxes)
Anova(finalmethanemodel)
AICc(finalmethanemodel) # -1851.34
AIC(finalmethanemodel) # -1878.613
summary(finalmethanemodel)
#### Post-hoc ####
### Emmeans
methane_emm <- emmeans(finalmethanemodel, ~ Reach|Date,
                       type = "response")
methane_emm_sum <- summary(methane_emm)

### CLD 
ch4_cld <- cld(methane_emm,
               by = c("Site", "Date"),
               alpha = 0.05, 
               Letters = letters,
               decreasing = TRUE)

#### Methane Plots ####
ggplot(data = methane_emm_sum) +
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
        axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  facet_grid(~Site)


ch4sum1 <- CH4_fluxes %>%
  group_by(Site, Reach) %>%
  summarise(Avg = mean(CH4_flux_g))
# FH    BDA   -0.0000478
# FH    REF    0.000550 
## 11.50628 higher in REF
# LP    BDA    0.000798 
# LP    REF    0.000649 
## 1.229584 higher in BDA 
# TP    BDA    0.000546 
# TP    REF    0.000852
## 1.56044 higher in REF

ch4sum1 <- CH4_fluxes %>%
  group_by(Reach) %>%
  summarise(Avg = mean(CH4_flux_g))
# BDA    1.68
# REF    1.96
## 1.166667

methanerange1 <- CH4_fluxes %>%
  group_by(Site) %>%
  summarise(Range = max(CH4_flux_g)- min(CH4_flux_g))
# Range values for each Site across Date and Reach 
# FH    0.00338
# LP    0.00633
# TP    0.00668

methanerange2 <- CH4_fluxes %>%
  group_by(Reach) %>%
  summarise(Range = max(CH4_flux_g)- min(CH4_flux_g))
# Range values across reaches
# BDA    0.00698
# REF    0.00512

methanerange3 <- CH4_fluxes %>%
  summarise(Range = max(CH4_flux_g)- min(CH4_flux_g))
# Range value across Date, Site, and Reach
# 0.006977156

#### Soil Moisture ####

HydraGO_data <- CO2_fluxes %>%
  select(Date, Site, Chamber, Reach, Soil_Moisture_Percent, Soil_Temp_C, Soil_Moisture_wfv, EC)

## Plots
HydraGO_data %>%
  group_by(Date, Site, Reach) %>%
  summarise(Avg = mean(Soil_Moisture_Percent)) %>%
  ggplot()+
  geom_line(aes(Date, Avg, color = Reach, group = Reach), lwd = 1.5) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  scale_color_manual(name = "Reach", labels = c("Treatment", "Reference"), values = c("Blue", "Purple")) +
  scale_fill_manual(name = "Reach", labels = c("Treatment", "Reference"), values = c("Blue", "Purple")) + 
  labs(title = "Riparian Soil Moisture", 
       x = "Date", 
       y = "Soil Moisture (%)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(colour = "black"), 
        axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  facet_grid(~Site)


