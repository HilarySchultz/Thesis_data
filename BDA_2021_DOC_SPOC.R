# library(Rmisc)
library(dplyr)
library(tidyverse)
# library(lme4)
# library(lsmeans)
# library(car)
# library(readtext)
# library(vroom)
library(lubridate)

#### DOC #### 
DOC <- read.csv("DOC_Data.csv", header = T, sep = ",") %>% 
  filter(Conc_ppm > 0, Reach != "GS") %>% # Filtering out null data and removing GS
  select(-Location)

# Identifying outliers. The function identify_outliers() creates two new columns that ID outliers (1.5)
# and extreme outliers (3), Which are 1.5 and 3 standard deviations away from the mean, respectively. 
# I removed extreme outliers and doubled checked that those matched up with the observational data. 

outliers <- DOC %>%
  group_by(Date, Site) %>%
  identify_outliers("Conc_ppm") %>%
  filter(is.extreme == TRUE)

# Removed outliers from dataframe then calculated the avg concentration, mass, and area between replicates
DOC_data <- anti_join(DOC_data, outliers) %>%
  group_by(Date, Site, Reach, .drop = F) %>%
  summarise_at(c("Avg_area", "Mass_ug", "Conc_ppm"), mean) 

DOC_data %>% group_by(Site, Date, Reach) %>% 
  ggplot(aes(Site, Conc_ppm, fill = Site)) + 
  geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 2, color = "blue", fill = "blue") +
  labs(title = "Dissolved Organic Carbon", x = NULL, y = expression(DOC~(g~C/L))) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_blank(), 
                     axis.ticks = element_blank(), axis.text = element_text(colour = "black")) +
  scale_fill_brewer(palette = "Spectral", name = "Site") + 
  facet_grid(cols = vars(fct_relevel(Reach, 'REF', 'BDA', 'DS'))) +
  scale_y_continuous(limits = c(0, 8), breaks = c(0,2,4,6,8))

#### SPOC #####

# Reading in SPOC calculations and creating a summary table of SE and CI
SPOC_data <- read.csv("BDA_SPOC_Calc.csv", header = T, sep = ",") %>% 
  filter(SPOC > 0)
# SPOC_summary <- SPOC_data %>% summarySE(groupvars = c("Site", "Date", "Reach"), measurevar = "SPOC") %>%
#     rename(stdev = sd, stderror = se)
# names(SPOC_summary)

# Quick plots
SPOC_data %>% group_by(Site, Date, Reach) %>% summarise(AvgSPOC= mean(SPOC)) %>%
  ggplot(aes(Site, AvgSPOC, fill = Site)) + geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 2, color = "blue", fill = "blue") +
  labs(title = "Suspended Particulate Organic Carbon", x = NULL, y = "SPOC (mg C/L)") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_blank(), 
                     axis.ticks = element_blank(), axis.text = element_text(colour = "black")) +
  scale_fill_brewer(palette = "Spectral", name = "Site") + 
  facet_grid(cols = vars(fct_relevel(Reach, 'REF', 'GS', 'BDA', 'DS'))) +
  scale_y_continuous(limits = c(0,3), breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3))



# SPOC_data %>% summarySE(groupvars = c("Site", "Date", "Reach"), measurevar = "SPOC") %>%
#   ggplot(aes(Date, SPOC)) + geom_point() + geom_errorbar(aes(ymin = SPOC - sd, ymax = SPOC + sd)) +
#   facet_grid(Site~fct_relevel(Reach, 'GS', 'REF', 'BDA', 'DS'), space = "free") + 
#   theme(axis.text.x = element_text(angle =45))

# Plot data for those outlier days!!!!!! Not really sure what to do here. 
# SPOC_summary %>% filter(stdev>1) %>% ggplot(aes(Date, SPOC)) + geom_boxplot()

## GLM for SPOC
SPOC_data <- SPOC_data %>% 
  mutate(Date = as.factor(as.Date(mdy(Date))), 
          Site = as.factor(ordered(Site, levels = c("LP", "FH", "TP"))), 
         Reach = as.factor(ordered(Reach, levels = c("DS", "BDA", "REF", "GS"))))
 
SPOC.glm1 <- glm(SPOC ~ Site*Reach*Date, data = SPOC_data, family = Gamma(link = "log"))
Anova(SPOC.glm1) 
AIC(SPOC.glm1) # -97.16237; eliminating Site:Reach, Site:Date, Site:Reach:Date since p-values are high

SPOC.glm2 <- glm(SPOC ~ Site + Reach + Date + Reach*Date, data = SPOC_data, family = Gamma(link = "log") )
Anova(SPOC.glm2)
AIC(SPOC.glm2) # -97.16237; all interactions have an equally low p-value
par(mfrow = c(2,2)) 
plot(SPOC.glm2) # Diagnositc plots that show data points 68, 174, 240 as outliers. Should I slice them out?
summary(SPOC.glm2)

# Stuck here: do I call lsmeans on each parameter and/or on each interaction?
SPOC.glm2.lsmean <- lsmeans(SPOC.glm2, ~ Site, Reach, Date, Reach | Date , type="response")






