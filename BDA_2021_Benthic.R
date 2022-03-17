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

#### BENTHIC CALCULATIONS ####

# Reading in FBPOC data and doing the calculations to compare to google sheets
# Final calculations: Jaunary 2022
FBPOC_data <- read.csv("BDA_FBPOC_Calc.csv", header = T, sep = ",") %>% 
  select(2:8, 14:15) %>%
  rename(Dry_Weight_Filter = Dry.Weight...Filter..g., 
         Ashed_Weight_Filter = Ash...Filter..g., 
         Average_Depth = Average.Depth..cm., 
         Filtered_Volume = Filtered.Volume..ml.) %>%
  mutate(FBPOC_OM = Dry_Weight_Filter - Ashed_Weight_Filter, 
         Benthinator_Volume = pi*8.75^2*Average_Depth, 
         FBPOC_Total_OM = Benthinator_Volume/Filtered_Volume*FBPOC_OM, 
         FBPOC_Carbon_Mass_per_Area = FBPOC_Total_OM/240.53 * 1000 * 0.52)
# Units are in g C/m^2

# Pulling out the average volume column values
Benthinator_volume <- FBPOC_data %>% 
  select(Site, Reach, Location, Replicate, Date, Benthinator_Volume)

# Reading in BPOC data and doing the calculations to compare to google sheets
# You can pass multiple mutate arguments but I am separating them so I can have check points along the way
# Final calculations: Jaunary 2022
BPOC_data <- read.csv("BDA_CBPOC_Calc.csv", header = T, sep = ",") %>% select(3:10, 12) %>% # selecting certain columns
  rename(BPOC_Dry_Weight = Dry.Weight..g., BPOC_Tin_Weight = Tin.Weight..g., # renaming columns
         BPOC_Subsample_Dry_Tin = Subsample.Dry.Weight...Tin..g., BPOC_Ash_Tin = Tin...Ash..g.) %>%
  # "BPOC_Subsample_Dry_Tin" is the weight of the subsample + weight of the tin
  # "BPOC_Ash_Tin" is the weight of either the subsample or sample ash + weight of the tin
  mutate(BPOC_Subsample_Dry_Weight = ifelse(BPOC_Subsample_Dry_Tin > 0, BPOC_Subsample_Dry_Tin - BPOC_Tin_Weight, 0)) %>% 
  # Not all samples have subsamples so if there is a subsample, the ifelse finds which samples have subsamples
  # then subtracts the weight of the tin from the dry subsample in the tin. If no subsample then entry is 0. 
  mutate(BPOC_Ash = BPOC_Ash_Tin - BPOC_Tin_Weight) %>%
  # Determining the mass of the ash, so what's actually there after all OM burns off, so whats left is inorganic.
  mutate(BPOC_Subsample_OM = ifelse(BPOC_Subsample_Dry_Weight > 0, BPOC_Subsample_Dry_Weight - BPOC_Ash, 0)) %>%
  # Determining the mass of what actually burned off, which is the OM content. 
  # If there is a subsample then the ash is subtracted from the dry weight, 
  # meaning that the inorganic material is subtracted from the dry weight (inorganic + OM) to yield the OM content
  mutate(BPOC_Sample_OM = ifelse(BPOC_Subsample_Dry_Tin != 0, 
                                 BPOC_Dry_Weight/BPOC_Subsample_Dry_Weight*BPOC_Subsample_OM, BPOC_Dry_Weight - BPOC_Ash)) %>%
  # Determining the sample OM mass based on a sample/subsample ratio
  # If there is a subsample, then ((dry_sample/dry_subsample) = (x/subsample_OM)), where x = sample OM
  # If there is no subsample, then the sample ash weight is subtracted from the sample dry weight to yield OM content
  mutate(BPOM_mass_per_area = BPOC_Sample_OM/240.53*1000) %>%
  # Sampler area is 240.53 cm^2 (radius of 8.75 cm), when we change to m^2
  # The units here are (g/m^2)
  mutate(BPOC_Carbon_Mass_per_Area = BPOM_mass_per_area * 0.52)
# 52% carbon 
# Units here are (g C/m^2)

# Need the filtered volume for the suspended CBPOC samples (SBPOC)
Benthic_field_data <- read.csv("Benthic_Field_Data.csv", header = T, sep = ",") %>%
  select(2:6, 11) %>%
  rename(Filtered_sample_volume = Filtered.BCPOC..ml.)

volumes <- full_join(Benthic_field_data, Benthinator_volume)

# Reading in SBPOC data and doing the calculations to compare to google sheets
# Final calculations: Jaunary 2022
SBPOC <- read.csv("BDA_SCBPOC_Calc.csv", header = T, sep = ",") %>%
  select(1, 3:10, 12) # selecting certain columns

SBPOC_data <- full_join(SBPOC, volumes) %>% # Joining the benthinator and sample volumes in with the df
  rename(SBPOC_Dry_Weight = Dry.Weight..g., SBPOC_Tin_Weight = Tin.Weight..g., 
         SBPOC_Subsample_Dry_Tin = Subsample.Dry.Weight...Tin..g., SBPOC_Ash_Tin = Tin...Ash..g.) %>%
  mutate(SBPOC_Subsample_Dry_Weight = ifelse(SBPOC_Subsample_Dry_Tin > 0, SBPOC_Subsample_Dry_Tin - SBPOC_Tin_Weight, 0),
         SBPOC_Ash = SBPOC_Ash_Tin - SBPOC_Tin_Weight, 
         SBPOC_Subsample_OM = ifelse(SBPOC_Subsample_Dry_Weight > 0, SBPOC_Subsample_Dry_Weight - SBPOC_Ash, 0), 
         SBPOC_Sample_OM = ifelse(SBPOC_Subsample_Dry_Tin != 0, 
                                  SBPOC_Dry_Weight/SBPOC_Subsample_Dry_Weight*SBPOC_Subsample_OM, SBPOC_Dry_Weight - SBPOC_Ash)) %>%
  # Above are the same steps as the CBPOC, however, I need to normalize by the volume of water that was taken.
  mutate(SBPOM_Total_Mass = (Benthinator_Volume/Filtered_sample_volume)*SBPOC_Sample_OM) %>%
  # Determining the total amount of OM in the total benthinator volume
  # Units are in g
  mutate(SBPOM_mass_per_area = SBPOM_Total_Mass/240.53*1000) %>%
  mutate(SBPOC_Carbon_Mass_per_Area = SBPOM_mass_per_area * 0.52)
# Units are in g C/m^2  

# Binding benthic CBPOC, suspended CBPOC, and FBPOC into one df then selecting rows of interest
CBPOC_data <- full_join(BPOC_data, SBPOC_data) # Coarse samples df

Benthic_data <- full_join(CBPOC_data, FBPOC_data) %>% # Coarse with fine samples into df
  select("Site", "Reach", "Location", "Replicate", "Date", "BPOC_Carbon_Mass_per_Area", 
         "SBPOC_Carbon_Mass_per_Area", "FBPOC_Carbon_Mass_per_Area") %>%
  mutate(Total_BPOC_Mass_per_Area = BPOC_Carbon_Mass_per_Area + SBPOC_Carbon_Mass_per_Area + FBPOC_Carbon_Mass_per_Area)

# Need to change dates so that they represent the three sampling trips (June, July, August)
Benthic_data <- Benthic_data %>%
  mutate(Date = case_when(Date == "6/14/2021" ~ "June", 
                          Date == "6/15/2021" ~ "June", 
                          Date == "6/17/2021" ~ "June", 
                          Date == "7/5/2021" ~ "July",
                          Date == "7/6/2021" ~ "July", 
                          Date == "7/8/2021" ~ "July",
                          Date == "8/2/2021" ~ "August",
                          Date == "8/3/2021" ~ "August",
                          Date == "8/5/2021" ~ "August"))

# write.csv(Benthic_data, "Benthic_data", row.names = F)

benthic_outliers <- Benthic_data %>%
  rstatix::identify_outliers("Total_BPOC_Mass_per_Area") %>%
  filter(is.extreme == TRUE)

Benthic_data <- anti_join(Benthic_data, benthic_outliers)

#### Simple Plots ####
# Histogram to see the distribution of the data
Benthic_data %>% ggplot(aes(Total_BPOC_Mass_per_Area)) +
  geom_histogram(binwidth = 10) +
  facet_grid(rows = vars(Reach)) 

# Boxplots for SE presentation
# Benthic_data %>% group_by(Site, Date, Reach, Location) %>% summarise(AvgBPOC= mean(Total_BPOC_Mass_per_Area)) %>%
#   ggplot(aes(Site, AvgBPOC, fill = Site)) + geom_boxplot() +
#   stat_summary(fun = mean, geom = "point", shape = 20, size = 2, color = "blue", fill = "blue") +
#   labs(title = "Benthic Particulate Organic Carbon", x = NULL, y = expression(BPOC~(g~C/m^2))) +
#   theme_bw() + theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_blank(), 
#                      axis.ticks = element_blank(), axis.text = element_text(colour = "black")) +
#   scale_fill_brewer(palette = "Spectral", name = "Site") + 
#   facet_grid(cols = vars((Reach))) +
#   scale_y_continuous(limits = c(0,120), breaks = c(0,15,30,45,60,75,90,105,120))


#### Models for Benthic Data ####
# Changing the class of each variable. Categorical variables enter into R differently than continuous variables
# storing data as factors insures taht the modeling functions will treat such data correctly. 
# Factors in R are stored as a vector of integer values with a corresponding set of character values to use
# when the factor is displayed. Factors levels are always character values. 

# Changing Site, Reach, Location, Date, and Replicate to ordered factors.
Benthic_data[,1:5] <- lapply(Benthic_data[,1:5], as.factor)

Benthic_data$Date <- ordered(Benthic_data$Date, levels = c("June", "July", "August"))
Benthic_data$Location <- ordered(Benthic_data$Location, levels = c("LWR", "MID", "UPR"))
Benthic_data$Replicate <- ordered(Benthic_data$Replicate, levels = c(1:3))
#### Final GLMMs ####

### Date and Reach = Fixed effects, Site = random effect
# Trying to determine the best model by changing up location and replicate. 

# Location = fixed, replicate = random
bmodel <- glmer(Total_BPOC_Mass_per_Area ~ Date + Reach + Location + (1|Site) + (1|Replicate),
                data = Benthic_data,
                family = Gamma(link="log")) 
plot(bmodel) # Bunched in the lower left-hand corner. 
Anova(bmodel) # Date, Reach, Location are ***
AICc(bmodel)  # 1344.603
AIC(bmodel) # 1343.387
summary(bmodel)

# Location = random, replicate = random 
bmodel1 <- glmer(Total_BPOC_Mass_per_Area ~ Date + Reach + (1|Location) + (1|Site) + (1|Replicate),
                data = Benthic_data, 
                family = Gamma(link="log")) 
plot(bmodel1) # Bunched in the lower left-hand corner. 
Anova(bmodel1) # Date and Reach are ***
AICc(bmodel1)  # 1349.593
AIC(bmodel1) # 1348.627
summary(bmodel1)

# Location = random, replicate = NULL
bmodel2 <- glmer(Total_BPOC_Mass_per_Area ~ Date + Reach + (1|Location) + (1|Site),
                 data = Benthic_data, 
                 family = Gamma(link="log")) 
plot(bmodel2) # Bunched in the lower left-hand corner, less so than previous bmodels 
Anova(bmodel2) # Date and Reach are ***
AICc(bmodel2)  # 1355.239
AIC(bmodel2) # 1354.493

# Location = fixed, replicate = NULL
bmodel3 <- glmer(Total_BPOC_Mass_per_Area ~ Date + Reach + Location + (1|Site),
                 data = Benthic_data, 
                 family = Gamma(link="log")) 
plot(bmodel3) # Bunched in the lower left-hand corner. 
Anova(bmodel3) # Date, Reach, Location are ***
AICc(bmodel3)  # 1349.928
AIC(bmodel3) # 1348.961

test <- glmer(Total_BPOC_Mass_per_Area ~ Site/Reach/Location*Date + (1|Replicate),
              data = Benthic_data, 
              nAGQ = 0,
              family = Gamma(link="log"))
plot(test)
Anova(test)
AICc(test)  # 1381.576
AIC(test) # 1318.368

finalbenthicmodel <- glmer(Total_BPOC_Mass_per_Area ~ Date/Site/Reach + (1|Location/Replicate),
              data = Benthic_data,
              family = Gamma(link="log"))
plot(test1)
Anova(test1)
AICc(test1)  # 1301.195
AIC(test1) # 1294.401

test2 <- glmer(Total_BPOC_Mass_per_Area ~ Date*Site*Reach + (1|Location/Replicate),
               data = Benthic_data,
               family = Gamma(link="log"))
plot(test2)
Anova(test2)
AICc(test2)  # 1307.7
AIC(test2) # 1300.906

anova(test1, test2)

bmodelnest <- glmer(Total_BPOC_Mass_per_Area ~ Date + Reach + (1|Location/Replicate) + (1|Site),
                data = Benthic_data,
                family = Gamma(link="log")) 
plot(bmodelnest) # Bunched in the lower left-hand corner. 
Anova(bmodelnest) # Date, Reach, Location are ***
AICc(bmodelnest)  # 1339.944
AIC(bmodelnest) # 1338.978
# Here, with the random effects nesting we are saying:
# The intercept varies among locations and among replicates within locations (nested)
# (1|Location/Replicate) = (1|Location) + (1|Location:Replicate)

overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

overdisp_fun(bmodelnest) # This model has less dispersion
overdisp_fun(bmodel)
overdisp_fun(finalbenthicmodel)


quasi_table <- function(model,ctab=coef(summary(model)),
                        phi=overdisp_fun(model)["ratio"]) {
  qctab <- within(as.data.frame(ctab),
                  {   `Std. Error` <- `Std. Error`*sqrt(phi)
                  `z value` <- Estimate/`Std. Error`
                  `Pr(>|z|)` <- 2*pnorm(abs(`z value`), lower.tail=FALSE)
                  })
  return(qctab)
}
printCoefmat(quasi_table(bmodelnest),digits=3)


#### Post-hoc tests ####

### Emmeans
# Reach by Date and Location
benthic_emm <- emmeans(finalbenthicmodel, ~ Reach|Date|Site,
                       type = "response")
benthic_emm_sum <- summary(benthic_emm)

# Reach by Reach
# benthic_reach_emm <- emmeans(finalbenthicmodel, ~ Reach,
#                              type = "response")
# benthic_reach_emm_sum <- summary(benthic_reachbyreach_emm)

# Reach by Date
# benthic_reachdate_emm <- emmeans(finalbenthicmodel, ~ Reach|Date,
#                              type = "response")
# benthic_reachdate_emm_sum <- summary(benthic_reachdate_emm)

### CLD
# Reach by Date and Location
benthic_reach_cld <- cld(benthic_emm,
                 by = c("Site", "Date"),
                 alpha = 0.05, 
                 Letters = letters,
                 decreasing = TRUE)
benthic_reach_cld$.group = gsub(" ", "", benthic_reach_cld$.group)
benthic_reach_cld <- arrange(benthic_reach_cld, Reach, Site, Date)




# Reach by Reach
benthic_cld <- cld(benthic_reachbyreach_emm,
                         by = c("Site", "Date"),
                         alpha = 0.05, 
                         Letters = letters,
                         decreasing = TRUE)
benthic_cld$.group = gsub(" ", "", benthic_cld$.group)
benthic_cld <- arrange(benthic_cld, Reach)

# Reach by Date
# benthic_reachdate_cld <- cld(benthic_reachdate_emm,
#                                 by = "Date",
#                                 alpha = 0.05, 
#                                 Letters = letters,
#                                 decreasing = TRUE)
# benthic_reachdate_cld$.group = gsub(" ", "", benthic_reachdate_cld$.group)
# benthic_reachdate_cld <- arrange(benthic_reachdate_cld, Reach, Date)


# Asterisks
# reach_cld$.group <- if_else(reach_cld$.group == "b", "*","")


#### PLOTS ####
# Reach by Date 
ggplot() +
  geom_boxplot(data = Benthic_data, aes(x = Reach, y = Total_BPOC_Mass_per_Area, fill = Reach)) +
  geom_point(data = benthic_cld, aes(x = Reach, y = response), size = 1, shape = 19,
             color = "blue") +
  geom_text(data = benthic_cld, aes(x = Reach, y = response, label= .group,
                                       vjust = -2.1, hjust = 0.5),
            size = 6, position = position_dodge(0.5), color = "red") +
  geom_text(aes()) +
  scale_fill_manual(name = "Reach", labels = c("Treatment", "Reference"), values = c("#3399FF", "#CC99FF")) +
  # scale_fill_brewer(palette = "Spectral") +
  labs(title = "Benthic Particulate Organic Carbon Pools", 
       x = NULL, 
       y = expression(BPOC~(g~C~m^-2))) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(colour = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  theme(axis.text = element_text(size = 12)) +
  facet_grid(Date~Site) 

# Reach comparison  
ggplot() +
  geom_boxplot(data = Benthic_data, aes(x = Reach, y = Total_BPOC_Mass_per_Area, fill = Reach)) +
  geom_point(data = benthic_reachbyreach_cld, aes(x = Reach, y = response), size = 1, shape = 19,
             color = "blue") +
  geom_text(data = benthic_reachbyreach_cld, aes(x = Reach, y = response, label= .group,
                                              vjust = -2.1, hjust = 0.5),
            size = 6, position = position_dodge(0.5), color = "red") +
  geom_text(aes()) +
  scale_fill_manual(name = "Reach", labels = c("Treatment", "Reference"), values = c("#3399FF", "#CC99FF")) +
  # scale_fill_brewer(palette = "Spectral") +
  labs(title = "Benthic Particulate Organic Carbon Pools", 
       x = NULL, 
       y = expression(BPOC~(g~C~m^-2))) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(colour = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  theme(axis.text = element_text(size = 12)) 

# Reach by Date and Location
ggplot() +
  geom_boxplot(data = Benthic_data, aes(x = Reach, y = Total_BPOC_Mass_per_Area, fill = Reach)) +
  geom_point(data = benthic_reach_cld, aes(x = Reach, y = response), size = 1, shape = 19,
             color = "blue") +
  geom_text(data = benthic_reach_cld, aes(x = Reach, y = response, label= .group,
                                          vjust = -2.1, hjust = 0.5),
            size = 6, position = position_dodge(0.5), color = "red") +
  geom_text(aes()) +
  scale_fill_manual(name = "Reach", labels = c("Treatment", "Reference"), values = c("#3399FF", "#CC99FF")) +
  # scale_fill_brewer(palette = "Spectral") +
  labs(title = "Benthic Particulate Organic Carbon Pools", 
       x = NULL, 
       y = expression(BPOC~(g~C~m^-2))) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(colour = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  theme(axis.text = element_text(size = 12)) +
  facet_grid(Date~Location)

#rep(-1.8, 5), 2, rep(-1.8, 5)


