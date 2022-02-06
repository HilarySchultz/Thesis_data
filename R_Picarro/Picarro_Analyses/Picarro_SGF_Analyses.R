# library(multcomp) # for doing elegant pairwise comparisons
library(tidyverse)  # loads several useful packages, including ggplot2, 
# tidyr, and dplyr
# library(Rmisc)  # for summarySE function to summarize data frame
# library(boot) # for diagnostics for GLM
# library(lme4)   # for linear, general, and nonlinear mixed models
# library(car)    # for doing ANOVA
# library(emmeans) # posthoc tests
# library(multcompView) # making posthoc tests easier to view and plot
# library(AICcmodavg) # for comparing models
# library(PerformanceAnalytics) # for making amazing correlation plots
options(scipen = 999)

# Reading in .csv of SGF, chamber field data, and HydraGO data
CO2_fluxes <- read.csv("CO2_fluxes.csv", header = T, sep = ",")
CH4_fluxes <- read.csv("CH4_fluxes.csv", header = T, sep = ",")

# Setting theme


#             theme(axis.line = element_line(linetype = "solid"), 
#                   axis.ticks = element_line(colour = "black", size = 1), 
#                   panel.grid.major = element_line(colour = NA), 
#                   panel.grid.minor = element_line(colour = NA),
#                   axis.title = element_text(size = 22),
#                   axis.text = element_text(size = 20, colour = "black"),
#                   plot.title = element_text(size = 16),
#                   panel.background = element_rect(fill = NA),
#                   plot.background = element_rect(colour = NA),
#                   legend.key = element_blank(),
#                   legend.text = element_text(size = 14),
#                   legend.title = element_text(size = 14)))

# CO2 fluxes by site
CO2_fluxes %>% ggplot(aes(Site, CO2_flux_g, fill = TP))


CO2_fluxes %>% ggplot(aes(Site, CO2_flux_g, fill = Site)) + geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 2, color = "blue", fill = "blue") +
  labs(title = expression(Riparian~Soil~CO[2]~Fluxes), x = NULL, y = expression(g~C /m^2/day)) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_blank(), 
                     axis.ticks = element_blank(), axis.text = element_text(colour = "black")) +
  scale_fill_brewer(palette = "Spectral", name = "Site") + facet_grid(cols = vars(Reach))

# CH4 fluxes by site
CH4_fluxes %>% ggplot(aes(Site, CH4_flux_g, fill = Site)) + geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 2, color = "blue", fill = "blue") +
  labs(title = expression(Riparian~Soil~CH[4]~Fluxes), x = NULL, y = expression(g~C /m^2/day)) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_blank(), 
                     axis.ticks = element_blank(), axis.text = element_text(colour = "black")) +
  scale_fill_brewer(palette = "Spectral", name = "Site") + facet_grid(cols = vars(Reach))


CO2_fluxes %>% ggplot(aes(Site, Soil_Temp_C, fill = Site)) + geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 2, color = "blue", fill = "blue") +
  labs(title = "Riparian Soil Temperature", x = NULL, y = "Soil Temperature (°C)") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_blank(), 
                     axis.ticks = element_blank(), axis.text = element_text(colour = "black")) +
  scale_fill_brewer(palette = "Spectral", name = "Site") + facet_grid(cols = vars(Reach))
