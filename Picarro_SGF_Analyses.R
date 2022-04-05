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

# write.csv(CH4_fluxes, "Methane_data.csv", row.names = F)
# write.csv(CO2_fluxes, "CarbonDioxide_data.csv", row.names = F)

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

#### Stats ####


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
  scale_x_discrete(expand = c(0,0), guide = guide_axis(angle = 45)) +
  scale_color_manual(name = "Reach", labels = c("Treatment", "Reference"), values = c("Blue", "Purple")) +
  scale_fill_manual(name = "Reach", labels = c("Treatment", "Reference"), values = c("Blue", "Purple")) + 
  labs(title = "Riparian Soil Methane Fluxes", 
       x = "Date", 
       y = expression(CH[4]~Fluxes~(g~C~m^-2~d^-1))) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(colour = "black")) +
  facet_grid(rows = vars(Site))

#### Stats ####


#### Soil Moisture ####

HydraGO_data <- CO2_fluxes %>%
  select(Date, Site, Chamber, Reach, Soil_Moisture_Percent, Soil_Temp_C, Soil_Moisture_wfv, EC)

hist(HydraGO_data$Soil_Moisture_wfv)

soillm <- lmer(Soil_Moisture_wfv ~ Date/Site/Reach + (1|Chamber),
                                         data = HydraGO_data)
plot(soillm)
Anova(soillm)

soillm_emm <- emmeans(soillm, ~ Reach|Date|Site,
                        type = "response")
soillm_emm_sum <- summary(soillm_emm)

soil_mois_cld <- cld(soillm_emm,
               by = c("Site", "Date"),
               alpha = 0.05, 
               Letters = letters,
               decreasing = TRUE)


## Plots
ggplot(data = soillm_emm_sum) +
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
  scale_x_discrete(expand = c(0,0), guide = guide_axis(angle = 45)) +
  scale_color_manual(name = "Reach", labels = c("Treatment", "Reference"), values = c("Blue", "Purple")) +
  scale_fill_manual(name = "Reach", labels = c("Treatment", "Reference"), values = c("Blue", "Purple")) + 
  labs(title = "Riparian Soil Moisture", 
       x = "Date", 
       y = "Soil Moisture (wfv)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(colour = "black")) +
  facet_grid(rows = vars(Site))

## Only line plot
HydraGO_data %>%
  group_by(Date, Site, Reach) %>%
  summarise(Avg = mean(Soil_Moisture_Percent)) %>%
  ggplot(aes(Date, Avg, color = Reach, group = Reach)) +
  geom_line(lwd = 1) +
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
  facet_grid(rows = vars(Site))

#### Temperature 
hist(HydraGO_data$Soil_Temp_C)

templm <- lmer(Soil_Temp_C ~ Date/Site/Reach + (1|Chamber),
               data = HydraGO_data)
plot(templm)
Anova(templm)

templm_emm <- emmeans(templm, ~ Reach|Date|Site,
                      type = "response")
templm_emm_sum <- summary(templm_emm)

templm_cld <- cld(templm_emm,
                     by = c("Site", "Date"),
                     alpha = 0.05, 
                     Letters = letters,
                     decreasing = TRUE)


## Temperature plot
ggplot(data = templm_emm_sum) +
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
  scale_x_discrete(expand = c(0,0), guide = guide_axis(angle = 45)) +
  scale_color_manual(name = "Reach", labels = c("Treatment", "Reference"), values = c("Blue", "Purple")) +
  scale_fill_manual(name = "Reach", labels = c("Treatment", "Reference"), values = c("Blue", "Purple")) + 
  labs(title = "Riparian Soil Temperature", 
       x = "Date", 
       y = "Soil Temperature (°C)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(colour = "black")) +
  facet_grid(rows = vars(Site)) 
  
### EC 
hist(HydraGO_data$EC)

eclm <- glmer(EC ~ Date/Site/Reach + (1|Chamber),
               data = HydraGO_data,
              control = glmerControl(optimizer = "bobyqa",
                                     optCtrl = list(maxfun = 100000)),
              family = Gamma("log"))
plot(eclm)
Anova(eclm)

eclm_emm <- emmeans(eclm, ~ Reach|Date|Site,
                      type = "response")
eclm_emm_sum <- summary(eclm_emm)

## Ribbon plot
ggplot(data = eclm_emm_sum) +
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
  labs(title = "Riparian Soil EC", 
       x = "Date", 
       y = "Soil EC(S/m)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(colour = "black")) +
  facet_grid(rows = vars(Site)) 

## line plot only
HydraGO_data %>%
  group_by(Date, Site, Reach) %>%
  summarise(Avg = mean(EC)) %>%
  ggplot(aes(Date, Avg, color = Reach, group = Reach)) +
  geom_line(lwd = 1) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  scale_color_manual(name = "Reach", labels = c("Treatment", "Reference"), values = c("Blue", "Purple")) +
  scale_fill_manual(name = "Reach", labels = c("Treatment", "Reference"), values = c("Blue", "Purple")) + 
  labs(title = "Riparian Soil EC", 
       x = "Date", 
       y = "Soil EC (S/m)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(colour = "black"), 
        axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  facet_grid(rows = vars(Site))

