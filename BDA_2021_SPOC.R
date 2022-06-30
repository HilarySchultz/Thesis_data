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

# END WEIGHT UNITS ARE IN MG/L
SPOC_data <- read.csv("BDA_SPOC_Calc.csv", header = T, sep = ",") %>% 
  select(Date, Site, Reach, Replicate, Dry_Weight_and_Filter, Ash_and_Filter) %>%
  mutate(AFDW_mg = (Dry_Weight_and_Filter - Ash_and_Filter) * 1000, # change to mg
         SPOC = AFDW_mg * 0.52) %>% # assume 0.52 percent of OM is carbon
  # Sample volume is 1 L, so no unit conversion required since we want L
  # Units are in mg/L
  filter(SPOC > 0, Reach != "GS", Reach != "DS") %>%
  rename(Time = Date, Segment = Reach)

# spoc_outliers <- SPOC_data %>%
#   rstatix::identify_outliers("SPOC") %>%
#   filter(is.extreme == TRUE)

SPOC_data <- anti_join(SPOC_data, spoc_outliers)

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

SPOC_data$Time <- as.factor(SPOC_data$Time)
SPOC_data$Site <- as.factor(SPOC_data$Site)
SPOC_data$Segment <- as.factor(SPOC_data$Segment)
SPOC_data$Replicate <- as.factor(SPOC_data$Replicate)
SPOC_data$Date <- as.factor(SPOC_data$Date)

# SPOC_data$Date <- ordered(SPOC_data$Date,
#                          levels = c("6/14/2021", "6/28/2021","7/12/2021", "7/26/2021",
#                                     "8/9/2021", "8/25/2021", "9/7/2021"))

Q_2021 <- read.csv("Discharge_2021.csv", header = T, sep = ",")

# Site: TP
TP_spoc_dates <- SPOC_data %>% # pulling a column of distinct dates
  filter(Site == "TP") %>%
  distinct(Date) %>%
  pull(Date)

TP_SPOC <- SPOC_data %>% 
  filter(Site == "TP") %>%
  mutate(Station = case_when(Segment == "REF" ~ "US", # Matching gauging station to reach location
                             Segment == "GS" ~ "US",
                             Segment == "BDA" ~ "DS",
                             Segment == "DS" ~ "DS"))

TP_spoc_Q <- Q_2021 %>% # Matching distinct DOC sample dates to discharge dates
  filter(Site == "TP") %>% 
  filter(Date %in% TP_spoc_dates) %>%
  arrange(Date) %>%
  select(-year)

# UNITS ARE IN MG/S - DISCHARGE IS M^3/S
TP_SPOC_loads <- full_join(TP_SPOC, TP_spoc_Q) %>% # Load df
  mutate(SPOC_Load = SPOC * mean.discharge.cms * 1000) 


# Site: FH
FH_spoc_dates <- SPOC_data %>%
  filter(Site == "FH") %>%
  distinct(Date) %>%
  pull(Date)

FH_SPOC <- SPOC_data %>%
  filter(Site == "FH") %>%
  mutate(Station = case_when(Segment == "REF" ~ "US", # Matching gauging station to reach location
                             Segment == "GS" ~ "US",
                             Segment == "BDA" ~ "DS",
                             Segment == "DS" ~ "DS"))

# No discharge data past 9/1/21, so I need to use that date for 9/7/21
# Pulled discharge data from 9/1/2021 for US and DS
FH_date <- Q_2021 %>%
  filter(Site == "FH", Date == "09/01/2021")

FH_spoc_Q <- Q_2021 %>%
  filter(Site == "FH", 
         Date %in% FH_spoc_dates) %>%
  bind_rows(FH_date) %>%
  arrange(Date) %>%
  mutate(Date = ifelse(Date == "09/01/2021", "09/07/2021", Date)) %>%
  select(-year)

# UNITS ARE IN MG/S - DISCHARGE IS M^3/S
FH_SPOC_loads<- full_join(FH_SPOC, FH_spoc_Q) %>%
  mutate(SPOC_Load = SPOC * mean.discharge.cms*1000) 

# Site: LP
LP_spoc_dates <- SPOC_data %>%
  filter(Site == "LP") %>%
  distinct(Date) %>%
  pull(Date)

LP_SPOC <- SPOC_data %>%
  filter(Site == "LP") %>%
  mutate(Station = case_when(Segment == "REF" ~ "US", # Matching gauging station to reach location
                             Segment == "BDA" ~ "US"))

LP_spoc_Q <- Q_2021 %>% 
  filter(Site == "LP", Date %in% LP_spoc_dates, Station == "US") %>%
  arrange(Date) %>%
  select(-year)

# UNITS ARE IN MG/S - DISCHARGE IS M^3/S
LP_SPOC_loads<- full_join(LP_SPOC, LP_spoc_Q) %>%
  mutate(SPOC_Load = SPOC * mean.discharge.cms*1000) 

#### LOAD DF ####
half_spoc_load <- full_join(LP_SPOC_loads, TP_SPOC_loads)

SPOC_load_data <- full_join(half_spoc_load, FH_SPOC_loads) 

#### MODELS ####
spoc_load_outliers <- SPOC_load_data %>%
  group_by(Date, Site) %>%
  rstatix::identify_outliers("SPOC_Load") %>%
  filter(is.extreme == TRUE)

SPOC_data <- anti_join(SPOC_data, spoc_load_outliers)

SPOC_load_data %>% ggplot(aes(SPOC_Load)) +
  geom_histogram(binwidth = 0.5) 

spocmodel <- glmer(SPOC_Load ~ Date/Site/Segment + (1|Replicate),
                        data = SPOC_load_data, 
                        family = Gamma(link = "log"),
                   control = glmerControl(optimizer = "bobyqa",
                                          optCtrl = list(maxfun = 100000)))
plot(spocmodel)
Anova(spocmodel)
AIC(spocmodel) # -287.7
AICc(spocmodel) # -232.5

spoc_emm <- emmeans(spocmodel, pairwise ~ Segment|Date|Site, 
                   adjust = "none",
                   type = "response",
                   nesting = "Date %in% Site, Segment %in% (Site*Date)")
spoc_emm_sum <- summary(spoc_emm)
spoc_emmeans <- summary(spoc_emm$emmeans)
spoc_emm_contrast <- summary(spoc_emm$contrasts)

spoc_cld <- cld(spoc_emm,
               by = c("Site", "Date"),
               alpha = 0.05, 
               Letters = letters,
               decreasing = TRUE)
spoc_cld$.group = gsub(" ", "", spoc_cld$.group)
spoc_cld <- arrange(spoc_cld, Date, Site, Segment)
spoc_cld$.group <- if_else(spoc_cld$.group == "b", "*","")

#### PLOTS ####
ggplot(data = spoc_emmeans) +
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
  geom_text(data = spoc_cld, aes(x = Date, y = response, label= .group,
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
       y = expression(SPOC~Loads~(mg~s^-1))) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(colour = "black", size = 12),
        panel.spacing.x = unit(1, "lines"),
        axis.title = element_text(size = 12),
        axis.title.y = element_text(size = 14)) +
  # axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5)) +
  facet_grid(rows = vars(Site))


spocmodel1 <- lmer(log(SPOC_Load) ~ Date/Site/Segment + (1|Replicate),
                   data = SPOC_load_data)
plot(spocmodel1)
Anova(spocmodel1)
AIC(spocmodel1) # 156.6936
AICc(spocmodel1) # 206.8202
