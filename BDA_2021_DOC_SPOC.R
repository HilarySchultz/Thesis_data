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

#### DOC #### 
DOC <- read.csv("DOC_Data.csv", header = T, sep = ",") %>% 
  filter(Conc_ppm > 0, Reach != "GS", Reach != "DS") %>% # Filtering out null data and removing GS and DS
  select(-Location)

# Identifying outliers. The function identify_outliers() creates two new columns that ID outliers (1.5)
# and extreme outliers (3), Which are 1.5 and 3 standard deviations away from the mean, respectively. 
# I removed extreme outliers and doubled checked that those matched up with the observational data. 

outliers <- DOC %>%
  group_by(Date, Site) %>%
  rstatix::identify_outliers("Conc_ppm") %>%
  filter(is.extreme == TRUE)

# Removed outliers from dataframe changing dates
DOC_data <- anti_join(DOC, outliers) %>%
  rename(Time = Date) 

# Setting this so the sampling dates are all at the start of the week
DOC_data <- DOC_data %>%
  mutate(Time = if_else(Time == "07/14/2021", "07/15/2021", Time)) %>% # Date was off 
  mutate(Date = case_when(Time == "05/31/2021" ~ "6/1/2021",
                          Time == "06/02/2021" ~ "6/1/2021",
                          Time == "06/03/2021" ~ "6/1/2021",
                          Time == "06/07/2021" ~ "6/7/2021",
                          Time == "06/08/2021" ~ "6/7/2021",
                          Time == "06/10/2021" ~ "6/7/2021",
                          Time == "06/14/2021" ~ "6/14/2021",
                          Time == "06/15/2021" ~ "6/14/2021",
                          Time == "06/17/2021" ~ "6/14/2021",
                          Time == "06/28/2021" ~ "6/28/2021",
                          Time == "06/29/2021" ~ "6/28/2021",
                          Time == "07/01/2021" ~ "6/28/2021",
                          Time == "07/12/2021" ~ "7/12/2021",
                          Time == "07/13/2021" ~ "7/12/2021",
                          Time == "07/15/2021" ~ "7/12/2021",
                          Time == "07/26/2021" ~ "7/26/2021",
                          Time == "07/27/2021" ~ "7/26/2021",
                          Time == "07/29/2021" ~ "7/26/2021", 
                          Time == "08/09/2021" ~ "8/9/2021",
                          Time == "08/10/2021" ~ "8/9/2021",
                          Time == "08/12/2021" ~ "8/9/2021",
                          Time == "08/25/2021" ~ "8/25/2021",
                          Time == "08/26/2021" ~ "8/25/2021",
                          Time == "08/29/2021" ~ "8/25/2021",
                          Time == "09/07/2021" ~ "9/7/2021",
                          Time == "09/08/2021" ~ "9/7/2021",
                          Time == "09/10/2021" ~ "9/7/2021" ))

# Changing Time, Date, Sample, Site, Reach, Replicate to ordered factors
DOC_data[,1:5] <- lapply(DOC_data[,1:5], as.ordered)
DOC_data$Date <- ordered(DOC_data$Date, 
                         levels = c("6/1/2021", "6/7/2021", "6/14/2021", "6/28/2021", 
                                    "7/12/2021", "7/26/2021", "8/9/2021", "8/25/2021", "9/7/2021"))

#### DOC Models ####
# Histogram to see the distribution of the data
# Data is right skewed, and sort of has a bimodal character
DOC_data %>% ggplot(aes(Conc_ppm)) +
  geom_histogram(binwidth = 0.1) 

# Trying to normalize the distribution
DOC_data$transformed <- abs(DOC_data$Conc_ppm - mean(DOC_data$Conc_ppm))
shapiro.test(DOC_data$transformed)
hist(DOC_data$transformed)

#### Exploratory/Bad GLMMs ####

# Models are singular - no good, overfitting.
DOC_fullnested_trans <- glmer(transformed ~ Date/Site/Reach + (1|Replicate),
                        data = DOC_data, 
                        family = Gamma(link = "log"))
isSingular(DOC_fullnested, tol = 0.0001) #TRUE

DOC_fullnested_trans1 <- glmer(transformed ~ Date/Site*Reach + (1|Replicate),
                              data = DOC_data, 
                              family = Gamma(link = "log"))

DOC_fullnested_trans2 <- glmer(transformed ~ Date*Site*Reach + (1|Replicate),
                               data = DOC_data, 
                               family = Gamma(link = "log"))

DOC_test2 <- glmer(Conc_ppm ~ Site/(1|Replicate),
                   data = DOC_data, 
                   family = Gamma(link = "log"))

DOC_test3 <- glmer(Conc_ppm ~ Site + (1|Replicate),
                   data = DOC_data, 
                   family = Gamma(link = "log"))

DOC_fullnested <- glmer(Conc_ppm ~ Date/Site/Reach + (1|Replicate),
                        data = DOC_data, 
                        family = Gamma(link = "log"))

DOC_fullnested1 <- glmer(Conc_ppm ~ Date/Site*Reach + (1|Replicate),
                        data = DOC_data, 
                        family = Gamma(link = "log"))


DOC_fullnested2 <- glmer(Conc_ppm ~ Date*Site*Reach + (1|Replicate),
                         data = DOC_data, 
                         family = Gamma(link = "log"))

DOC_fullnested3 <- glmer(Conc_ppm ~ Date*Site*Reach + (1|Replicate),
                         data = DOC_data, 
                         family = gaussian(link = "log"))

DOC_test <- glmer(Conc_ppm ~ Date/Site/Reach/(1|Replicate),
                  data = DOC_data, 
                  family = Gamma(link = "log"))

# This is without Date. Only site is significant but still singular
DOC_fullnested4 <- glmer(Conc_ppm ~ Site/Reach + (1|Replicate),
                         data = DOC_data, 
                         family = Gamma(link = "log"))
plot(DOC_fullnested4)
Anova(DOC_fullnested4)

DOC_fullnested5 <- glmer(Conc_ppm ~ Site*Reach + (1|Replicate),
                         data = DOC_data, 
                         family = Gamma(link = "log"))
plot(DOC_fullnested5)
Anova(DOC_fullnested5)

DOC_test1 <- glmer(Conc_ppm ~ Site/Reach/(1|Replicate),
                   data = DOC_data, 
                   family = Gamma(link = "log"))
plot(DOC_test1)
Anova(DOC_test1)

anova(DOC_fullnested4, DOC_fullnested5, DOC_test1)

# This is only BDA and REF reach and Dates reordered correctly.
# Trying to get date in the model - model is singular when everything is in date
DOC_datetest1.1 <- glmer(Conc_ppm ~ Site/Reach + Date + (1|Replicate),
                         data = DOC_data, 
                         family = Gamma(link = "log"))
plot(DOC_datetest1.1)
Anova(DOC_datetest1.1)

DOC_datetest2.1 <- glmer(Conc_ppm ~ Site/Reach*Date + (1|Replicate),
                         data = DOC_data, 
                         family = Gamma(link = "log"))
plot(DOC_datetest2.1)
Anova(DOC_datetest2.1)

# Date shouldn't be nested within Site and Reach.
DOC_datetest3 <- glmer(Conc_ppm ~ Site/Reach/Date + (1|Replicate),
                       data = DOC_data, 
                       family = Gamma(link = "log"))
plot(DOC_datetest3)
Anova(DOC_datetest3)

#### DOC GLMs ####
DOC_fullglm <- glm(Conc_ppm ~ Date*Site*Reach*Replicate, 
                   data = DOC_data,
                   family = gaussian(link = "log"))
par(mfrow = c(2,2))
plot(DOC_fullglm) # a lot of points with leverage
Anova(DOC_fullglm)
AICc(DOC_fullglm) # 1150.068
AIC(DOC_fullglm) # 349.8864

# There is a four way interaction but I'm thinking that this is due to the unnested 
# nature of the way replicate is included in the model. 

DOC_fullglm2 <- glm(Conc_ppm ~ Site*Reach*Replicate, 
                    data = DOC_data,
                    family = gaussian(link = "log"))

par(mfrow = c(2,2))
plot(DOC_fullglm2) 
Anova(DOC_fullglm2)
AICc(DOC_fullglm2) # 867.6869
AIC(DOC_fullglm2)  # 864.0505
# Only site is significant

DOC_fullglm3 <- glm(Conc_ppm ~ Site, 
                    data = DOC_data,
                    family = gaussian(link = "log"))
par(mfrow = c(2,2))
plot(DOC_fullglm3) 
Anova(DOC_fullglm3)
AICc(DOC_fullglm3) # 838.7048
AIC(DOC_fullglm3) # 838.5263

DOC_gammma <- glm(Conc_ppm ~ Site*Reach*Replicate, 
                  data = DOC_data,
                  family = Gamma(link = "log"))
par(mfrow = c(2,2))
plot(DOC_gammma) 
Anova(DOC_gammma)
AICc(DOC_gammma) # 787.5597
AIC(DOC_gammma) # 782.9233

DOC_gammma2 <- glm(Conc_ppm ~ Site, 
                   data = DOC_data,
                   family = Gamma(link = "log"))
par(mfrow = c(2,2))
plot(DOC_gammma2) 
Anova(DOC_gammma2)
AICc(DOC_gammma2) # 759.2363
AIC(DOC_gammma2) # 759.0577

# 4-way interaction
DOC_gammma3 <- glm(Conc_ppm ~ Date*Site*Reach*Replicate, 
                   data = DOC_data,
                   family = Gamma(link = "log"))
par(mfrow = c(2,2))
plot(DOC_gammma3) 
Anova(DOC_gammma3)
AICc(DOC_gammma3) # 1184.621
AIC(DOC_gammma3) #384.4397

# Reduce
DOC_gammma4 <- glm(Conc_ppm ~ Date:Site:Reach:Replicate + Date:Site:Replicate +
                     Date:Site:Reach + Site:Reach + Date:Site + Replicate +Site +Date, 
                   data = DOC_data,
                   family = Gamma(link = "log"))
par(mfrow = c(2,2))
plot(DOC_gammma4) 
Anova(DOC_gammma4)
AICc(DOC_gammma4) # 1184.621
AIC(DOC_gammma4) #384.4397

# REduce
DOC_gammma5 <- glm(Conc_ppm ~ Date:Site:Reach:Replicate + 
                     Date:Site:Replicate + Date:Site + Date + Site, 
                   data = DOC_data,
                   family = Gamma(link = "log"))
par(mfrow = c(2,2))
plot(DOC_gammma5) 
Anova(DOC_gammma5)
AICc(DOC_gammma5) # 664.8359
AIC(DOC_gammma5) # 462.3233

DOC_gammma6 <- glm(Conc_ppm ~ Date:Site + Site +Date, 
                   data = DOC_data,
                   family = Gamma(link = "log"))
par(mfrow = c(2,2))
plot(DOC_gammma6) 
Anova(DOC_gammma6)
AICc(DOC_gammma6) # 664.8359
AIC(DOC_gammma6)# 463.3233

#### Final Models ####
doctest <- glm(Conc_ppm ~ Date/Site/Reach,
                              data = DOC_data, 
                              family = Gamma(link = "log"))
par(mfrow = c(2,2))
plot(doctest)
Anova(doctest)
AIC(doctest) # 303.7063
AICc(doctest) # 367.2115

doctest1 <- glm(Conc_ppm ~ Date/Site + Date,
               data = DOC_data, 
               family = Gamma(link = "log"))
par(mfrow = c(2,2))
plot(doctest1)
Anova(doctest1)
AIC(doctest1) # 291.9806
AICc(doctest1) # 305.0774

#### Post-hoc tests ####
doc_emm <- emmeans(doctest, ~ Reach|Site|Date,
                       type = "response")
doc_emm <- summary(doc_emm)

### CLD
doc_reach_cld <- cld(doc_emm,
                 by = c("Site", "Date"),
                 alpha = 0.05, 
                 Letters = letters,
                 decreasing = TRUE)
doc_reach_cld$.group = gsub(" ", "", doc_reach_cld$.group)
doc_reach_cld <- arrange(doc_reach_cld, Reach, Site, Date)


#### DOC Box Plots ####
# Boxplots don't work because you only have 3 points. 
ggplot() +
  geom_boxplot(data = DOC_data, aes(x = Date, y = Conc_ppm, fill = Reach)) +
  geom_point(data = doc_reach_cld, aes(x = Date, y = response), size = 1, shape = 19,
             color = "blue") +
  geom_text(data = doc_reach_cld, aes(x = Reach, y = response, label= .group,
                                  vjust = -2.2, hjust = 0.5),
            size = 5, position = position_dodge(0.5), color = "black") +
  geom_text(aes()) +
  # scale_fill_manual(name = "Reach", labels = c("BDA", "Reference"), values = c("#3399FF", "#CC99FF")) +
  scale_fill_brewer(palette = "Spectral", name = "Reach", labels = c("BDA", "Reference")) +
  labs(title = "Dissolved Organic Carbon Concentrations", 
       x = NULL, 
       y = expression(DOC~(g~C~ml^-1))) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(colour = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  theme(axis.text = element_text(size = 12)) +
  facet_grid(Date~Site)   

#### DOC Ribbon Plot ####
ggplot(data = doc_emm) +
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
  labs(title = "Dissolved Organic Carbon Concentrations", 
       x = "Date", 
       y = "Estimated Marginal Means") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(colour = "black"), 
        axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  facet_grid(vars(rows = Site)) 




#### SPOC #####

# Reading in SPOC calculations and creating a summary table of SE and CI
SPOC_data <- read.csv("BDA_SPOC_Calc.csv", header = T, sep = ",") %>% 
  filter(SPOC > 0, Reach != "GS") %>%
  rename(Time = Date)

# SPOC_summary <- SPOC_data %>% summarySE(groupvars = c("Site", "Date", "Reach"), measurevar = "SPOC") %>%
#     rename(stdev = sd, stderror = se)
# names(SPOC_summary)

SPOC_data <- SPOC_data %>%
  mutate(Time = if_else(Time == "8/9/2021", "08/09/2021", Time),
         Time = if_else(Time == "8/10/2021", "08/10/2021", Time), 
         Time = if_else(Time == "8/10/2021", "08/10/2021", Time),
         Time = if_else(Time == "8/10/2021", "08/10/2021", Time),
         Time = if_else(Time == "8/12/2021", "08/12/2021", Time),
         Time = if_else(Time == "8/25/2021", "08/25/2021", Time),
         Time = if_else(Time == "8/26/2021", "08/26/2021", Time),
         Time = if_else(Time == "8/29/2021", "08/29/2021", Time),
         Time = if_else(Time == "9/7/2021", "09/07/2021", Time),
         Time = if_else(Time == "9/8/2021", "09/08/2021", Time),
         Time = if_else(Time == "9/10/2021", "09/10/2021", Time)) %>%
  mutate(Date = case_when(Time == "06/14/2021" ~ "6/14/2021",
                          Time == "06/15/2021" ~ "6/14/2021",
                          Time == "06/17/2021" ~ "6/14/2021",
                          Time == "06/28/2021" ~ "6/28/2021",
                          Time == "06/29/2021" ~ "6/28/2021",
                          Time == "07/01/2021" ~ "6/28/2021",
                          Time == "07/12/2021" ~ "7/12/2021",
                          Time == "07/13/2021" ~ "7/12/2021",
                          Time == "07/15/2021" ~ "7/12/2021",
                          Time == "07/26/2021" ~ "7/26/2021",
                          Time == "07/27/2021" ~ "7/26/2021",
                          Time == "07/29/2021" ~ "7/26/2021", 
                          Time == "08/09/2021" ~ "8/9/2021",
                          Time == "08/10/2021" ~ "8/9/2021",
                          Time == "08/12/2021" ~ "8/9/2021",
                          Time == "08/25/2021" ~ "8/25/2021",
                          Time == "08/26/2021" ~ "8/25/2021",
                          Time == "08/29/2021" ~ "8/25/2021",
                          Time == "09/07/2021" ~ "9/7/2021",
                          Time == "09/08/2021" ~ "9/7/2021",
                          Time == "09/10/2021" ~ "9/7/2021" ))

SPOC_data[,1:5] <- lapply(SPOC_data[,1:5], as.ordered)
SPOC_data$Date <- as.ordered(SPOC_data$Date)
SPOC_data$Time <- as.ordered(SPOC_data$Time)

hist(SPOC_data$SPOC)

#### GLMMs for SPOC ####
# SPOC_data <- SPOC_data %>% 
#   mutate(Date = as.factor(as.Date(mdy(Date))), 
#           Site = as.factor(ordered(Site, levels = c("LP", "FH", "TP"))), 
#          Reach = as.factor(ordered(Reach, levels = c("DS", "BDA", "REF"))))

# Rank defficient
SPOC_nested <- glmer(SPOC ~ Date/Site/Reach + (1|Replicate),
                              data = SPOC_data, 
                              family = Gamma(link = "log"))

# did not converge
SPOC_nested1 <- glmer(SPOC ~ Site/Reach*(Date) + (1|Replicate),
                         data = SPOC_data, 
                         family = Gamma(link = "log"))
SPOC_full <- glmer(SPOC ~ Site*Reach*Date + (1|Replicate),
                      data = SPOC_data, 
                      family = Gamma(link = "log"))
# Singular
SPOC_GLMM <- glmer(SPOC ~ Site*Reach + (1|Replicate),
                   data = SPOC_data, 
                   family = Gamma(link = "log"))
SPOC_GLMM2 <- glmer(SPOC ~ Reach + (1|Replicate),
                   data = SPOC_data, 
                   family = Gamma(link = "log"))
SPOC_GLMM3 <- glmer(SPOC ~ Site + (1|Replicate),
                    data = SPOC_data, 
                    family = Gamma(link = "log"))
SPOC_GLMM4 <- glmer(SPOC ~ Date/Site/Reach/(1|Replicate),
                   data = SPOC_data, 
                   family = Gamma(link = "log"))

#### GLMs for SPOC ####
SPOC.glm1 <- glm(SPOC ~ Site*Reach*Date, data = SPOC_data, family = Gamma(link = "log"))
par(mfrow = c(2,2))
plot(SPOC.glm1)
Anova(SPOC.glm1) 
AIC(SPOC.glm1) # -99.00951
AICc(SPOC.glm1) # -34.35993
# All interactions are ***

# Model did not converge.
SPOC.glm2 <- glm(SPOC ~ Date*Site*Reach*Replicate, 
                 data = SPOC_data, 
                 family = Gamma(link = "log"))
par(mfrow = c(2,2))
plot(SPOC.glm2)
Anova(SPOC.glm2) 
AIC(SPOC.glm2) # -99.00951
AICc(SPOC.glm2)

SPOC.glm3 <- glm(SPOC ~ Site*Reach*Replicate, 
                 data = SPOC_data, 
                 family = Gamma(link = "log"))
plot(SPOC.glm3)
Anova(SPOC.glm3) 
AIC(SPOC.glm3) # 20.2126
AICc(SPOC.glm3) # 24.99248

SPOC.glm4 <- glm(SPOC ~ Site:Reach + Reach + Site, 
                 data = SPOC_data, 
                 family = Gamma(link = "log"))
plot(SPOC.glm4)
Anova(SPOC.glm4) 
AIC(SPOC.glm4) # 15.34293
AICc(SPOC.glm4) # 16.65246
