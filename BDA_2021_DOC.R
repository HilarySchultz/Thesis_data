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

#### Cleaning DOC data #### 
DOC <- read.csv("DOC_Data.csv", header = T, sep = ",") %>%
  filter(Conc_ppm > 0) %>%
  select(Date, Site, Reach, Replicate, Conc_ppm)

# Identifying outliers. The function identify_outliers() creates two new columns that ID outliers (1.5)
# and extreme outliers (3), Which are 1.5 and 3 standard deviations away from the mean, respectively. 
# I removed extreme outliers and doubled checked that those matched up with the observational data. 

doc_outliers <- DOC %>%
  # group_by(Date, Site) %>%
  rstatix::identify_outliers("Conc_ppm") %>%
  filter(is.extreme == TRUE)

# Removed outliers from dataframe changing dates
# There were two: FH REF 6/15/21, FH DS 6/02/21
DOC_data <- anti_join(DOC, doc_outliers) %>%
  rename(Time = Date) # renaming column so I can mutate the dates into lumped dates

# Setting this so the sampling dates are all at the start of the week
DOC_data <- DOC_data %>%
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

# Changing Time, Date, Sample, Site, Reach, Replicate to ordered factors
DOC_data$Time <- as.ordered(DOC_data$Time)
DOC_data$Date <- as.ordered(DOC_data$Date)

# DOC_data$Date <- ordered(DOC_data$Date, 
#                          levels = c("6/1/2021", "6/7/2021", "6/14/2021", "6/28/2021", 
#                                     "7/12/2021", "7/26/2021", "8/9/2021", "8/25/2021", "9/7/2021"))
# DOC_data$Replicate <- ordered(DOC_data$Replicate, levels = c(1:3))

write.csv(DOC_data, "DOC_cleaned_data.csv", row.names = F)

#### Load Calcs ####
# Reading in discharge data to calculate loads and separating for each site to simplify
Q_2021 <- read.csv("Discharge_2021.csv", header = T, sep = ",")

# Site: TP
TP_dates <- DOC_data %>% # pulling a column of distinct dates
  filter(Site == "TP") %>%
  distinct(Date) %>%
  pull(Date)

TP_DOC <- DOC_data %>% 
  filter(Site == "TP") %>%
  mutate(Station = case_when(Reach == "REF" ~ "US", # Matching gauging station to reach location
                             Reach == "GS" ~ "US",
                             Reach == "BDA" ~ "DS",
                             Reach == "DS" ~ "DS"))

TP_Q <- Q_2021 %>% # Matching distinct DOC sample dates to discharge dates
  filter(Site == "TP") %>% 
  filter(Date %in% TP_dates) %>%
  arrange(Date) %>%
  select(-year)

TP_DOC_loads <- full_join(TP_DOC, TP_Q) %>% # Load df
  mutate(DOC_Load = Conc_ppm * mean.discharge.cms) %>%
  rename(monitoring_location = Reach)

# Site: FH
FH_dates <- DOC_data %>%
  filter(Site == "FH") %>%
  distinct(Date) %>%
  pull(Date)

FH_DOC <- DOC_data %>%
  filter(Site == "FH") %>%
  mutate(Station = case_when(Reach == "REF" ~ "US", # Matching gauging station to reach location
                             Reach == "GS" ~ "US",
                             Reach == "BDA" ~ "DS",
                             Reach == "DS" ~ "DS"))

# No discharge data for 6/1/21 so I will pull from 6/2/2021 - just changed the date
# No discharge data past 9/1/21, so I need to use that date for 9/7/21
# Pulled discharge data from 9/1/2021 for US and DS
FH_new_date <- Q_2021 %>%
  filter(Site == "FH", Date == "09/01/2021")

FH_Q <- Q_2021 %>%
  mutate(Date = ifelse(Date == "06/02/2021", "06/01/2021", Date)) %>%
  filter(Site == "FH", 
         Date %in% FH_dates) %>%
  bind_rows(FH_new_date) %>%
  arrange(Date) %>%
  mutate(Date = ifelse(Date == "09/01/2021", "09/07/2021", Date)) %>%
  select(-year)
# Changed the discharge date so that it would match with the DOC date.

FH_DOC_loads<- full_join(FH_DOC, FH_Q) %>%
  mutate(DOC_Load = Conc_ppm * mean.discharge.cms) %>%
  rename(monitoring_location = Reach)

# Site: LP
LP_dates <- DOC_data %>%
  filter(Site == "LP") %>%
  distinct(Date) %>%
  pull(Date)

LP_DOC <- DOC_data %>%
  filter(Site == "LP") %>%
  mutate(Station = case_when(Reach == "REF" ~ "US", # Matching gauging station to reach location
                             Reach == "GS" ~ "US",
                             Reach == "BDA" ~ "US",
                             Reach == "DS" ~ "DS"))

# Using date 6/10 becasuse info starts there for 6/1 and 6/7 for US and DS
# Filling in the DS info
LP_6.1 <- Q_2021 %>%
  filter(Site == "LP", Date == "06/10/2021") %>%
  mutate(Date = ifelse(Date == "06/10/2021", "06/01/2021", Date))

LP_6.7 <- Q_2021 %>%
  filter(Site == "LP", Date == "06/10/2021") %>%
  mutate(Date = ifelse(Date == "06/10/2021", "06/07/2021", Date))

LP_DS8.9 <- Q_2021 %>%
  filter(Site == "LP", Station == "DS", Date == "08/08/2021") %>%
  mutate(Date = ifelse(Date == "08/08/2021", "08/09/2021", Date))

LP_DS8.25 <- Q_2021 %>%
  filter(Site == "LP", Station == "DS", Date == "08/18/2021") %>%
  mutate(Date = ifelse(Date == "08/18/2021", "08/25/2021", Date))

LP_DS9.7 <- Q_2021 %>%
  filter(Site == "LP", Station == "DS", Date == "08/18/2021") %>%
  mutate(Date = ifelse(Date == "08/18/2021", "09/07/2021", Date))

LP_Q <- Q_2021 %>% 
  filter(Site == "LP", Date %in% LP_dates) %>%
  bind_rows(LP_6.1, LP_6.7, LP_DS8.9, LP_DS8.25, LP_DS9.7) %>%
  arrange(Date) %>%
  select(-year)

LP_DOC_loads<- full_join(LP_DOC, LP_Q) %>%
  mutate(DOC_Load = Conc_ppm * mean.discharge.cms) %>%
  rename(monitoring_location = Reach)

#### LOAD DF ####
half_load <- full_join(LP_DOC_loads, TP_DOC_loads)

# Filtering out GS and DS because we don't need it now. 
filter <- c("BDA", "REF")

Load_data <- full_join(half_load, FH_DOC_loads) %>%
  filter(monitoring_location %in% filter) %>%
  rename(Segment = monitoring_location)

#### MODELS ####
Load_data %>% ggplot(aes(DOC_Load)) +
  geom_histogram(binwidth = 0.1) 

docmodel <- glmer(DOC_Load ~ Date/Site/Segment + (1|Replicate),
                       data = Load_data,
                       family = Gamma(link = "log"))
plot(docmodel)
qqnorm(residuals(docmodel))
Anova(docmodel)
AIC(docmodel) # -1863.199
BIC(docmodel) # -1691.341
AICc(docmodel) # -1800.611 

doc_emm <- emmeans(docmodel, pairwise ~ Segment|Date|Site, 
                   adjust = "none",
                   type = "response",
                   nesting = "Date %in% Site, Segment %in% (Site*Date)")
doc_emm_sum <- summary(doc_emm)
doc_emmeans <- summary(doc_emm$emmeans)
doc_emm_contrast <- summary(doc_emm$contrasts)

doc_cld <- cld(doc_emm,
               by = c("Site", "Date"),
               alpha = 0.05, 
               Letters = letters,
               decreasing = TRUE)
doc_cld$.group = gsub(" ", "", doc_cld$.group)
doc_cld <- arrange(doc_cld, Date, Site, Segment)
doc_cld$.group <- if_else(doc_cld$.group == "b", "*","")

#### DOC Ribbon Plot ####
ggplot(data = doc_emmeans) +
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
                color = Segment), 
            lwd = 1) +
  geom_text(data = doc_cld, aes(x = Date, y = response, label= .group,
                                vjust = -1.5, hjust = 0.5),
            size = 6, position = position_dodge(0.5), color = "black") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(guide = guide_axis(angle = 45), 
                   labels = c("June 1, 2021", "", "June 14,2021", "", "July 12, 2021", "",
                              "August 9, 2021", "", "September 7, 2021"),
                   # labels = c("June","June","June", "June", "July","July", "August", "August", "September"),
                   expand = c(0,0)) +
  scale_color_manual(name = "Segment", labels = c("Treatment", "Reference"), values = c("#0072B2", "#009E73")) +
  scale_fill_manual(name = "Segment", labels = c("Treatment", "Reference"), values = c("#0072B2", "#009E73")) +
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



docmodel1 <- lmer(log(DOC_Load) ~ Date/Site/Segment + (1|Replicate),
                  data = Load_data)
plot(docmodel1)
Anova(docmodel1)
AIC(docmodel1) # 178.4624
BIC(docmodel1) # 350.3211
AICc(docmodel1) # 178.4624

doc_emm1 <- emmeans(docmodel1, pairwise ~ Segment|Date|Site, 
                   adjust = "none",
                   type = "response",
                   nesting = "Date %in% Site, Segment %in% (Site*Date)")
doc_emm_sum1 <- summary(doc_emm1)
doc_emmeans1 <- summary(doc_emm1$emmeans)
doc_emm_contrast1 <- summary(doc_emm1$contrasts)

doc_cld1 <- cld(doc_emm1,
               by = c("Site", "Date"),
               alpha = 0.05, 
               Letters = letters,
               decreasing = TRUE)
doc_cld1$.group = gsub(" ", "", doc_cld1$.group)
doc_cld1 <- arrange(doc_cld1, Date, Site, Segment)
doc_cld1$.group <- if_else(doc_cld1$.group == "b", "*","")

#### DOC Ribbon Plot ####
ggplot(data = doc_emmeans1) +
  geom_ribbon(aes(x = Date,
                  ymin = lower.CL, 
                  ymax = upper.CL,
                  group = Segment,
                  fill = Segment), 
              alpha = 0.40, 
              color = NA) + # opaqueness of the CI
  # fill = "#3984ff") 
  geom_line(aes(x = Date, 
                y = response, 
                group = Segment, 
                color = Segment), 
            lwd = 1) +
  geom_text(data = doc_cld1, aes(x = Date, y = response, label= .group,
                                 vjust = -1.5, hjust = 0.5),
            size = 6, position = position_dodge(0.5), color = "black") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(guide = guide_axis(angle = 45), 
                   labels = c("June 1, 2021", "", "June 14,2021", "", "July 12, 2021", "",
                              "August 9, 2021", "", "September 7, 2021"),
                   # labels = c("June","June","June", "June", "July","July", "August", "August", "September"),
                   expand = c(0,0)) +
  scale_color_manual(name = "Segment", labels = c("Treatment", "Reference"), values = c("#0072B2", "#009E73")) +
  scale_fill_manual(name = "Segment", labels = c("Treatment", "Reference"), values = c("#0072B2", "#009E73")) +
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
ggplot(data = doc_emmeans1) +
  geom_ribbon(aes(x = Date,
                  ymin = lower.CL, 
                  ymax = upper.CL,
                  group = Segment,
                  fill = Segment), 
              alpha = 0.40, 
              color = NA) + # opaqueness of the CI
  # fill = "#3984ff") 
  geom_line(aes(x = Date, 
                y = response, 
                group = Segment, 
                color = Segment), 
            lwd = 1) +
  geom_text(data = doc_cld1, aes(x = Date, y = response, label= .group,
                                vjust = -1.5, hjust = 0.5),
            size = 6, position = position_dodge(0.5), color = "black") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(guide = guide_axis(angle = 45), 
                   labels = c("June 1, 2021", "", "June 14,2021", "", "July 12, 2021", "",
                              "August 9, 2021", "", "September 7, 2021"),
                   # labels = c("June","June","June", "June", "July","July", "August", "August", "September"),
                   expand = c(0,0)) +
  scale_color_manual(name = "Segment", labels = c("Treatment", "Reference"), values = c("#0072B2", "#009E73")) +
  scale_fill_manual(name = "Segment", labels = c("Treatment", "Reference"), values = c("#0072B2", "#009E73")) +
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

control = glmerControl(optimizer = "bobyqa",
                       optCtrl = list(maxfun = 100000))
