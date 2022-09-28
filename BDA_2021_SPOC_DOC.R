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
# library(ggpubr)
library(rstatix)
library(MuMIn)
library(gridExtra)
library(cowplot)

#### DOC #### 
DOC_data <- read.csv("DOC_Data.csv", header = T, sep = ",") %>%
  filter(Conc_ppm > 0, Reach != "GS", Reach != "DS") %>%
  select(Date, Site, Reach, Replicate, Conc_ppm) %>%
  rename(Segment = Reach, Time = Date) 

#### Outlier removal ####
# Identifying outliers. The function identify_outliers() creates two new columns that ID outliers (1.5)
# and extreme outliers (3), Which are 1.5 and 3 standard deviations away from the mean, respectively. 
# I removed extreme outliers and doubled checked that those matched up with the observational data. 
# Outlier removal due to obvious contaimination
doc_outliers <- DOC_data %>%
  rstatix::identify_outliers("Conc_ppm") %>%
  filter(is.extreme == TRUE)
  
DOC_data <- anti_join(DOC_data, doc_outliers)
  
# Setting this so the sampling dates are all at the start of the week
DOC_data <- DOC_data %>%
  mutate(Segment = if_else(Segment == "REF", "Reference", "Treatment")) %>%
  mutate(Time = if_else(Time == "07/14/2021", "07/15/2021", Time)) %>% # Date was off 
  mutate(Date = case_when(Time == "05/31/2021" ~ "06/01/2021",
                          Time == "06/02/2021" ~ "06/01/2021",
                          Time == "06/03/2021" ~ "06/01/2021",
                          Time == "06/07/2021" ~ "06/07/2021",
                          Time == "06/08/2021" ~ "06/07/2021",
                          Time == "06/10/2021" ~ "06/07/2021",
                          Time == "06/14/2021" ~ "06/14/2021",
                          Time == "06/15/2021" ~ "06/14/2021",
                          Time == "06/17/2021" ~ "06/14/2021",
                          Time == "06/28/2021" ~ "06/28/2021",
                          Time == "06/29/2021" ~ "06/28/2021",
                          Time == "07/01/2021" ~ "06/28/2021",
                          Time == "07/12/2021" ~ "07/12/2021",
                          Time == "07/13/2021" ~ "07/12/2021",
                          Time == "07/15/2021" ~ "07/12/2021",
                          Time == "07/26/2021" ~ "07/26/2021",
                          Time == "07/27/2021" ~ "07/26/2021",
                          Time == "07/29/2021" ~ "07/26/2021", 
                          Time == "08/09/2021" ~ "08/09/2021",
                          Time == "08/10/2021" ~ "08/09/2021",
                          Time == "08/12/2021" ~ "08/09/2021",
                          Time == "08/25/2021" ~ "08/25/2021",
                          Time == "08/26/2021" ~ "08/25/2021",
                          Time == "08/29/2021" ~ "08/25/2021",
                          Time == "09/07/2021" ~ "09/07/2021",
                          Time == "09/08/2021" ~ "09/07/2021",
                          Time == "09/10/2021" ~ "09/07/2021" )) 

# Time and Date as ordered factors; Segment, Site, and Replicate are factors
DOC_data$Time <- as.ordered(DOC_data$Time)
DOC_data$Date <- as.ordered(DOC_data$Date)
DOC_data$Segment <- as.factor(DOC_data$Segment)
DOC_data$Site <- as.factor(DOC_data$Site)
DOC_data$Replicate <- as.factor(DOC_data$Replicate)

# write.csv(DOC_data, "DOC_cleaned_data.csv", row.names = F)

#### DOC Models ####
# Histogram to see the distribution of the data
# Data is right skewed, and sort of has a bimodal character
DOC_data %>% 
ggplot(aes(Conc_ppm)) +
  geom_histogram(binwidth = 0.1) 

#### DOC GLMM ####
finaldocmodel <- glmer(Conc_ppm ~ Date/Site/Segment + (1|Replicate),
                  data = DOC_data, 
                  control = glmerControl(optimizer = "bobyqa",
                                         optCtrl = list(maxfun = 100000)),
                  family = Gamma(link = "log"))

doc_norm <- lmer(Conc_ppm ~ Date/Site/Segment + (1|Replicate),
                 data = DOC_data)

plot(finaldocmodel)
Anova(finaldocmodel)

qqnorm(resid(finaldocmodel))
qqline(resid(finaldocmodel))

# R2m = marginal - the variance explained by the fixed effects
# R2c = conditional - the variance explained by the entire model, including random and fixed
doc_fit <- r.squaredGLMM(finaldocmodel)

#### Post-hoc tests ####
### Emmeans
# Tukey adjusted, this is what I'm using 4/9
doc_emm <- emmeans(finaldocmodel, pairwise ~ Segment|Date|Site, 
                   adjust = "none",
                     type = "response",
                     nesting = "Date %in% Site, Segment %in% (Site*Date)")
doc_emm_sum <- summary(doc_emm)
doc_emmeans <- summary(doc_emm$emmeans)
doc_emm_contrast <- summary(doc_emm$contrasts) # p values

### CLD 
doc_cld <- cld(doc_emm,
                         by = c("Site", "Date"),
                         alpha = 0.05, 
                         Letters = letters,
                         decreasing = TRUE)
doc_cld$.group = gsub(" ", "", doc_cld$.group)
doc_cld <- arrange(doc_cld, Date, Site, Segment)
doc_cld$.group <- if_else(doc_cld$.group == "b", "*","")

#### DOC Ribbon Plot ####
doc_plot <-ggplot(data = doc_emmeans) +
  geom_ribbon(aes(x = Date,
                  ymin = asymp.LCL, 
                  ymax = asymp.UCL,
                  group = Segment,
                  fill = Segment), 
              alpha = 0.40, 
              color = NA) + # opaqueness of the CI
  # fill = "#3984ff") 
  geom_line(aes(x = Date, 
                y = response, 
                group = Segment, 
                color = Segment,
                linetype = Segment), 
            lwd = 1) +
  geom_text(data = doc_cld, aes(x = Date, y = response, label= .group,
                                 vjust = -1.5, hjust = 0.5),
            size = 6, position = position_dodge(0.5), color = "black") +
  ylim(0,10) +
  # scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(guide = guide_axis(angle = 45), 
                   # labels = c("June 1, 2021", "", "June 14,2021", "", "July 12, 2021", "",
                   #            "August 9, 2021", "", "September 7, 2021"),
                   # labels = c("June","June","June", "June", "July","July", "August", "August", "September"),
                   expand = c(0,0)) +
  scale_color_manual(values = c("#0072B2", "#009E73")) +
  scale_fill_manual(values = c("#0072B2", "#009E73")) +
  labs(title = NULL, 
       x = NULL, 
       y = expression(Dissolved~Organic~Carbon~(mg~L^-1))) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(colour = "black", size = 12),
        panel.spacing.x = unit(1, "lines"),
        axis.title = element_text(size = 12),
        axis.title.y = element_text(size = 14)) +
  # axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5)) +
  facet_grid(rows = vars(Site))

save_plot("docfile.png", doc_plot, base_height = 6, base_width = 9)

#### Stats ####
doc_date_sum <- doc_emmeans %>%
  group_by(Date, Site) %>%
  summarise(Avg = mean(response),
            mean_UCL = mean(asymp.UCL),
            mean_LCL = mean(asymp.LCL))

doc_site_sum <- doc_emmeans %>%
  group_by(Site) %>%
  summarise(Avg = mean(response),
            mean_UCL = mean(asymp.UCL),
            mean_LCL = mean(asymp.LCL))

doc_reach_sum <- doc_emmeans %>%
  group_by(Site, Segment) %>%
  summarise(Avg = mean(response),
            mean_UCL = mean(asymp.UCL),
            mean_LCL = mean(asymp.LCL)) 

std_error <- function(x) sd(x)/sqrt(length(x))

#### SPOC #####
SPOC_data <- read.csv("BDA_SPOC_Calc.csv", header = T, sep = ",") %>% 
  filter(SPOC > 0, Reach != "GS", Reach != "DS") %>%
  select(-Site.ID, -Sample, -Location) %>%
  rename(Time = Date, Segment = Reach) 

spoc_outliers <- SPOC_data %>%
  rstatix::identify_outliers("SPOC") %>%
  filter(is.extreme == T)

SPOC_data <- SPOC_data %>%
  anti_join(spoc_outliers)

SPOC_data <- SPOC_data %>%
  mutate(Segment = if_else(Segment == "REF", "Reference", "Treatment")) %>%
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
  mutate(Date = case_when(Time == "06/14/2021" ~ "06/14/2021",
                          Time == "06/15/2021" ~ "06/14/2021",
                          Time == "06/17/2021" ~ "06/14/2021",
                          Time == "06/28/2021" ~ "06/28/2021",
                          Time == "06/29/2021" ~ "06/28/2021",
                          Time == "07/01/2021" ~ "06/28/2021",
                          Time == "07/12/2021" ~ "07/12/2021",
                          Time == "07/13/2021" ~ "07/12/2021",
                          Time == "07/15/2021" ~ "07/12/2021",
                          Time == "07/26/2021" ~ "07/26/2021",
                          Time == "07/27/2021" ~ "07/26/2021",
                          Time == "07/29/2021" ~ "07/26/2021", 
                          Time == "08/09/2021" ~ "08/09/2021",
                          Time == "08/10/2021" ~ "08/09/2021",
                          Time == "08/12/2021" ~ "08/09/2021",
                          Time == "08/25/2021" ~ "08/25/2021",
                          Time == "08/26/2021" ~ "08/25/2021",
                          Time == "08/29/2021" ~ "08/25/2021",
                          Time == "09/07/2021" ~ "09/07/2021",
                          Time == "09/08/2021" ~ "09/07/2021",
                          Time == "09/10/2021" ~ "09/07/2021" ))

# Setting info to factors and numeric values
SPOC_data$Time <- as.ordered(SPOC_data$Time)
SPOC_data$Date <- as.ordered(SPOC_data$Date)
SPOC_data$Segment <- as.factor(SPOC_data$Segment)
SPOC_data$Site <- as.factor(SPOC_data$Site)
SPOC_data$Replicate <- as.factor(SPOC_data$Replicate)

hist(SPOC_data$SPOC)

# write.csv(SPOC_data, "SPOC_data.csv", row.names = F)

#### GLMM for SPOC ####
# SPOC_data <- SPOC_data %>% 
#   mutate(Date = as.factor(as.Date(mdy(Date)))

finalspocmodel <- glmer(SPOC ~ Date/Site/Segment + (1|Replicate),
                   data = SPOC_data, 
                   control = glmerControl(optimizer = "bobyqa",
                                          optCtrl = list(maxfun = 100000)),
                   family = Gamma(link = "log"))

plot(finalspocmodel)
Anova(finalspocmodel)
AIC(finalspocmodel) # -129.5816
AICc(finalspocmodel) # -79.45506
summary(finalspocmodel)

spoc_fit <- r.squaredGLMM(finalspocmodel)

#### Post-hoc Test
spoc_emm <- emmeans(finalspocmodel, pairwise ~ Segment|Date|Site, 
                     type = "response",
                     nesting = "Date %in% Site, Segment %in% (Site*Date)")
spoc_emm_sum <- summary(spoc_emm)
spoc_emmeans <- summary(spoc_emm$emmeans)
spoc_emm_contrast <- summary(spoc_emm$contrasts) # p values

### CLD
spoc_cld <- cld(spoc_emm,
               by = c("Site", "Date"),
               alpha = 0.05, 
               Letters = letters,
               decreasing = TRUE)
spoc_cld$.group = gsub(" ", "", spoc_cld$.group)
spoc_cld <- arrange(spoc_cld, Date, Site, Segment)

spoc_cld$.group <- if_else(spoc_cld$.group == "b", "*","")

#### SPOC Ribbon Plot ####
spoc_plot <- ggplot(data = spoc_emmeans) +
  geom_ribbon(aes(x = Date,
                  ymin = asymp.LCL, 
                  ymax = asymp.UCL,
                  group = Segment,
                  fill = Segment), 
              alpha = 0.40, 
              color = NA) + # opaqueness of the CI
  # fill = "#3984ff") +
  geom_line(aes(x = Date, 
                y = response, 
                group = Segment, 
                color = Segment,
                linetype = Segment),
            lwd = 1) +
  geom_text(data = spoc_cld, aes(x = Date, y = response, label= .group,
                                          vjust = -1, hjust = 0.5),
            size = 6, position = position_dodge(0.5), color = "black") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(guide = guide_axis(angle = 45),
  #                  labels = c("June 14,2021", "", "July 12, 2021", "",
  #                             "August 9, 2021", "", "September 7, 2021"),
                   # labels = c("June","June","June", "July","July", "August", "August", "September"),
                   expand = c(0,0)) +
  scale_color_manual(values = c("#0072B2", "#009E73")) +
  scale_fill_manual(values = c("#0072B2", "#009E73")) + 
  labs(title = NULL, 
       x = NULL, 
       y = expression(Suspended~Particulate~Organic~Carbon~(mg~L^-1))) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(colour = "black", size = 12),
        panel.spacing.x = unit(1, "lines"),
        axis.title = element_text(size = 12),
        axis.title.y = element_text(size = 14)) +
  # axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5)) +
  facet_grid(rows = vars(Site))

save_plot("spoc_plot.png", spoc_plot, base_height =6, base_width = 9)

#### Stats ####
emm_avg <- spoc_emmeans %>%
  group_by(Site) %>%
  summarise(mean_response = mean(response),
            mean_UCL = mean(asymp.UCL),
            mean_LCL = mean(asymp.LCL))

emm_avg1 <- spoc_emmeans %>%
  group_by(Site, Segment) %>%
  summarise(mean_response = mean(response),
            mean_UCL = mean(asymp.UCL),
            mean_LCL = mean(asymp.LCL)) 
pivot <- emm_avg1 %>%
  pivot_wider(names_from = Reach, values_from = c(mean_response, mean_UCL, mean_LCL)) 


FH <- SPOC_data %>%
  filter(Site == "FH") %>%
  pull(SPOC)
std_error(FH) # 0.0372177

LP <- SPOC_data %>%
  filter(Site == "LP") %>%
  pull(SPOC)
std_error(LP) # 0.06365265

TP <- SPOC_data %>%
  filter(Site == "TP") %>%
  pull(SPOC)
std_error(TP) # 0.02063711

sd(LP) # 0.4025747
sd(FH) #  0.2411983
sd(TP) #  0.1337438

LP_bda <- SPOC_data %>%
  filter(Site == "LP", Reach == "BDA") %>%
  pull(SPOC)
sd(LP_bda) # 0.4458614

LP_ref <- SPOC_data %>%
  filter(Site == "LP", Reach == "REF") %>%
  pull(SPOC)
sd(LP_ref) # 0.3567699

FH_bda <- SPOC_data %>%
  filter(Site == "FH", Reach == "BDA") %>%
  pull(SPOC)
sd(FH_bda) # 0.2357238

FH_ref <- SPOC_data %>%
  filter(Site == "FH", Reach == "REF") %>%
  pull(SPOC)
sd(FH_ref) # 0.2068724

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




