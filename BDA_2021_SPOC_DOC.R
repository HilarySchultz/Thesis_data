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
library(lubridate)
library(ggpubr)

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

# # Trying to normalize the distribution
# DOC_data$transformed <- abs(DOC_data$Conc_ppm - mean(DOC_data$Conc_ppm))
# shapiro.test(DOC_data$transformed)
# hist(DOC_data$transformed)

#### DOC GLMMs ####
# nAGQ = 0 is a less exact form of parameter estimation for GLMMs (helps model run)
docmodel <- glmer(Conc_ppm ~ Date + Reach + (1|Site) + (1|Replicate),
                  data = DOC_data, 
                  nAGQ = 0,
                  family = Gamma(link = "log"))
Anova(docmodel)
AIC(docmodel) # 358.2897
AICc(docmodel) # 360.9085


docmodel1 <- glmer(Conc_ppm ~ Date + Reach + (1|Site),
                  data = DOC_data, 
                  family = Gamma(link = "log"))
plot(docmodel1)
Anova(docmodel1)
AIC(docmodel1) # 356.2851
AICc(docmodel1) # 358.5137
summary(docmodel1)

docmodel2 <- glmer(Conc_ppm ~ Date + Reach + Site + (1|Replicate),
                   data = DOC_data,
                   nAGQ = 0,
                   family = Gamma(link = "log"))
plot(docmodel2)
Anova(docmodel2)
AIC(docmodel2) # 357.1709
AICc(docmodel2) # 360.2143
summary(docmodel2)


#### Post-hoc tests ####
doc_emm_reach <- emmeans(docmodel1, ~ Reach|Date,
                       type = "response")
doc_emm_reach_sum <- summary(doc_emm_reach)

##
doc_emm_date <- emmeans(docmodel1, ~ Date|Reach,
                   type = "response")
doc_emm_date_sum <- summary(doc_emm_date)

##
doc_rg <- ref_grid(docmodel1)
doc_emm_reachtest1 <- emmeans(ref_grid(doc_rg), "source")
doc_emm_reachtest <- regrid(emmeans(docmodel1, ~ Reach|Reach,
                         type = "response"))
doc_emm_reachtest_sum <- summary(doc_emm_reachtest)

### CLD
# Reach CLDs
doc_reach_cld <- cld(doc_emm_reach,
                 by = "Date",
                 alpha = 0.05, 
                 Letters = letters,
                 decreasing = TRUE)
doc_reach_cld$.group = gsub(" ", "", doc_reach_cld$.group)
doc_reach_cld <- arrange(doc_reach_cld, Reach, Date)

# Reach by Reach

doc_reach_cldtest <- cld(doc_emm_reachtest,
                     by = "Reach",
                     alpha = 0.05, 
                     Letters = letters,
                     decreasing = TRUE)
doc_reach_cldtest$.group = gsub(" ", "", doc_reach_cldtest$.group)
doc_reach_cldtest <- arrange(doc_reach_cldtest, Reach)

# Date CLD
doc_date_cld <- cld(doc_emm_date,
                     by = "Reach",
                     alpha = 0.05, 
                     Letters = letters,
                     decreasing = TRUE)
doc_date_cld$.group = gsub(" ", "", doc_date_cld$.group)
doc_date_cld <- arrange(doc_date_cld, Reach, Date)

#### DOC Box Plots ####
ggplot(data = DOC_data, aes(x = Reach, y = Conc_ppm)) +
  geom_boxplot(aes(fill = Reach)) +
  geom_point(data = doc_reach_cldtest, aes(x = Reach, y = response), size = 1, shape = 19,
             color = "blue") +
  geom_text(data = doc_reach_cldtest, aes(x = Reach, y = response, label= .group,
                                  vjust = -2.2, hjust = 0.5),
            size = 5, position = position_dodge(0.5), color = "black") +
  scale_fill_manual(name = "Reach", labels = c("BDA", "Reference"), values = c("#3399FF", "#CC99FF")) +
  # scale_fill_brewer(palette = "Spectral") +
  labs(title = "Dissolved Organic Carbon Concentrations", 
       x = NULL,
       y = expression(DOC~(g~C~ml^-1))) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(colour = "black", size = 12),
        # axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        legend.position = "none") +
  scale_x_discrete(labels = c("Treatment", "Reference"))

#### DOC Ribbon Plot ####
ggplot(data = doc_emm_sum) +
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
       y = expression(Concentration~(mg~C~L^-1))) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(colour = "black"), 
        axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  facet_grid(vars(rows = Site)) 

# Concentration by Date - fill = Reach
ggplot(data = doc_emm_reach_sum) +
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
       y = expression(Concentration~(mg~C~L^-1))) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(colour = "black"), 
        axis.text.x = element_text(angle = 45, vjust = 0.5)) 

#### SPOC #####

# Reading in SPOC calculations and creating a summary table of SE and CI
SPOC_data <- read.csv("BDA_SPOC_Calc.csv", header = T, sep = ",") %>% 
  filter(SPOC > 0, Reach != "GS", Reach != "DS") %>%
  select(!Site.ID) %>%
  rename(Time = Date)

outliers <- SPOC_data %>%
  group_by(Time, Site) %>%
  rstatix::identify_outliers("SPOC") %>%
  filter(is.extreme == TRUE)

SPOC_data <- anti_join(SPOC_data, outliers)

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
SPOC_data$Date <- ordered(SPOC_data$Date, 
                         levels = c("6/14/2021", "6/28/2021","7/12/2021", "7/26/2021", 
                                    "8/9/2021", "8/25/2021", "9/7/2021"))
hist(SPOC_data$SPOC)

#### GLMMs for SPOC ####
# SPOC_data <- SPOC_data %>% 
#   mutate(Date = as.factor(as.Date(mdy(Date))), 
#           Site = as.factor(ordered(Site, levels = c("LP", "FH", "TP"))), 
#          Reach = as.factor(ordered(Reach, levels = c("DS", "BDA", "REF"))))

SPOC_nested <- glmer(SPOC ~ Date/Site/Reach + (1|Replicate),
                              data = SPOC_data, 
                              family = Gamma(link = "log"))
Anova(SPOC_nested)
plot(SPOC_nested)
AIC(SPOC_nested) # -100.7296
AICc(SPOC_nested) # -49.301

# Response: SPOC
# Chisq Df Pr(>Chisq)    
# Date            155.04  6  < 2.2e-16 ***
#   Date:Site       472.28 14  < 2.2e-16 ***
#   Date:Site:Reach 241.90 21  < 2.2e-16 ***

SPOC_nested1 <- glmer(SPOC ~ Date/Site + Reach + (1|Replicate),
                     data = SPOC_data, 
                     family = Gamma(link = "log"))
Anova(SPOC_nested1)
plot(SPOC_nested1)
AIC(SPOC_nested1) # -28.56463
AICc(SPOC_nested1) # -16.1935

# Response: SPOC
# Chisq Df Pr(>Chisq)    
# Date       80.035  6  3.515e-15 ***
#   Reach      11.106  1  0.0008606 ***
#   Date:Site 199.946 14  < 2.2e-16 ***


SPOC_nested2 <- glmer(SPOC ~ Date/Site/(1|Replicate) + Reach,
                      data = SPOC_data, 
                      family = Gamma(link = "log"))
Anova(SPOC_nested2)
plot(SPOC_nested2)
AIC(SPOC_nested2) # -28.56463
AICc(SPOC_nested2) # -16.1935

SPOC_GLMM4 <- glmer(SPOC ~ Date/Site/Reach/(1|Replicate),
                    data = SPOC_data, 
                    family = Gamma(link = "log"))
Anova(SPOC_GLMM4)
plot(SPOC_GLMM4)
AIC(SPOC_GLMM4) # -100.7296
AICc(SPOC_GLMM4) # -49.301

# Response: SPOC
# Chisq Df Pr(>Chisq)    
# Date            155.04  6  < 2.2e-16 ***
#   Date:Site       472.28 14  < 2.2e-16 ***
#   Date:Site:Reach 241.90 21  < 2.2e-16 ***

spoctest <- SPOC_nested <- glm(SPOC ~ Date/Site/Reach,
                                 data = SPOC_data, 
                                 family = Gamma(link = "log"))
par(mfrow = c(2,2))
plot(spoctest)
Anova(spoctest)
AIC(spoctest) # -100.7296
AICc(spoctest) # -52.22545

# Response: SPOC
# LR Chisq Df Pr(>Chisq)    
# Date              127.76  6  < 2.2e-16 ***
#   Date:Site         323.29 14  < 2.2e-16 ***
#   Date:Site:Reach   134.76 21  < 2.2e-16 ***

spoctest2 <- SPOC_nested <- glm(SPOC ~ Date/Site + Reach,
                               data = SPOC_data, 
                               family = Gamma(link = "log"))
par(mfrow = c(2,2))
plot(spoctest2)
Anova(spoctest2)
AIC(spoctest2) # -32.33628
AICc(spoctest2) # -21.07098

spoctest3 <- SPOC_nested <- glm(SPOC ~ Date/Site*Reach,
                                data = SPOC_data, 
                                family = Gamma(link = "log"))
par(mfrow = c(2,2))
plot(spoctest3)
Anova(spoctest3)
AIC(spoctest3) # -100.7296
AICc(spoctest3) # -52.22545

# Not helpful
anova(SPOC_nested, SPOC_nested1, SPOC_nested2, SPOC_GLMM4)
anova(spoctest, spoctest2, spoctest3)

#### Post-hoc Tests ####
spoc_emm <- emmeans(SPOC_GLMM4, ~ Reach|Site|Date,
                   type = "response")

spoc_emm_sum <- summary(spoc_emm)

### CLD
spoc_reach_cld <- cld(spoc_emm,
                     by = c("Site", "Date"),
                     alpha = 0.05, 
                     Letters = letters,
                     decreasing = TRUE)
spoc_reach_cld$.group = gsub(" ", "", spoc_reach_cld$.group)
spoc_reach_cld <- arrange(spoc_reach_cld, Reach, Site, Date)

#### SPOC Ribbon Plot ####
ggplot(data = spoc_emm_sum) +
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
  labs(title = "Particulate Organic Carbon Concentrations", 
       x = "Date", 
       y = expression(Concentration~(mg~C~L^-1))) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(colour = "black"), 
        axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  facet_grid(vars(rows = Site))

#### Exploratory SPOC Models ####
SPOC_nested1 <- glmer(SPOC ~ Site/Reach*(Date) + (1|Replicate),
                         data = SPOC_data, 
                         family = Gamma(link = "log"))
SPOC_full <- glmer(SPOC ~ Site*Reach*Date + (1|Replicate),
                      data = SPOC_data, 
                      family = Gamma(link = "log"))
SPOC_GLMM <- glmer(SPOC ~ Site*Reach + (1|Replicate),
                   data = SPOC_data, 
                   family = Gamma(link = "log"))
SPOC_GLMM2 <- glmer(SPOC ~ Reach + (1|Replicate),
                   data = SPOC_data, 
                   family = Gamma(link = "log"))
SPOC_GLMM3 <- glmer(SPOC ~ Site + (1|Replicate),
                    data = SPOC_data, 
                    family = Gamma(link = "log"))
SPOC.glm1 <- glm(SPOC ~ Site*Reach*Date, 
                 data = SPOC_data, 
                 family = Gamma(link = "log"))
SPOC.glm2 <- glm(SPOC ~ Date*Site*Reach*Replicate, 
                 data = SPOC_data, 
                 family = Gamma(link = "log"))
SPOC.glm3 <- glm(SPOC ~ Site*Reach*Replicate, 
                 data = SPOC_data, 
                 family = Gamma(link = "log"))
SPOC.glm4 <- glm(SPOC ~ Site:Reach + Reach + Site, 
                 data = SPOC_data, 
                 family = Gamma(link = "log"))


#### DOC SPOC Ratios ####
DOC <- DOC_data %>%
  select(Site, Reach, Replicate, Date, Conc_ppm) %>%
  rename(DOC_conc = Conc_ppm)

DOC$Date <- as.character(DOC$Date)


SPOC <- SPOC_data %>%
  select(Site, Reach, Replicate, Date, SPOC) %>%
  rename(SPOC_conc = SPOC)
SPOC$Date <- as.character(SPOC$Date)


ratio_data <- inner_join(DOC, SPOC) %>%
  mutate(TOC_conc = DOC_conc + SPOC_conc,
         percent_DOC = DOC_conc/TOC_conc*100, 
         percent_SPOC = SPOC_conc/TOC_conc*100) 
  
ratio_data_sum <- ratio_data %>%
  group_by(Date, Site, Reach) %>%
  summarise(meanSPOCpercent = mean(percent_SPOC), 
            meanDOCpercent = mean(percent_DOC),
            meanSPOCconc = mean(SPOC_conc), 
            meanDOCconc = mean(DOC_conc))

test_DOC <- DOC %>%
  add_column(Sample = "DOC") %>%
  rename(Concentration = DOC_conc)

test_SPOC <- SPOC %>%
  add_column(Sample = "SPOC") %>%
  rename(Concentration = SPOC_conc)

test_OC_df <- full_join(test_DOC, test_SPOC) %>%
  mutate(TOC_conc = DOC_conc + SPOC_conc,
         percent_DOC = DOC_conc/TOC_conc*100, 
         percent_SPOC = SPOC_conc/TOC_conc*100) 

test_OC_df %>%
  filter(Sample == "SPOC") %>%
  ggplot(aes(Concentration)) +
  geom_histogram() + 
  facet_grid(vars(rows = Site), vars(cols = Reach))




