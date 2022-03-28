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
library(ggpubr)

BPOC_plot_data <- read.csv("Benthic_data.csv", header = T, sep = ",") %>%
  select(Date, Site, Reach, Location, Replicate,  Total_BPOC_Mass_per_Area) %>%
  rename(BPOC = Total_BPOC_Mass_per_Area)

SPOC_plot_data <- read.csv("SPOC_data.csv", header = T, sep = ",") %>%
  select(Date, Site, Reach, Replicate, SPOC) %>%
  group_by(Date, Site, Reach) %>%
  summarise(avg_SPOC = mean(SPOC))

DOC_plot_data <- read.csv("DOC_data.csv", header = T, sep = ",") %>%
  select(Date, Site, Reach, Conc_ppm) %>%
  rename(DOC = Conc_ppm) %>%
  group_by(Date, Site, Reach) %>%
  summarise(avg_DOC = mean(DOC))

SPOC_DOC_data <- full_join(SPOC_plot_data,DOC_plot_data) %>%
  filter(Date != "06/01/2021", Date != "06/07/2021") %>%
  mutate(Total_OC = avg_DOC + avg_SPOC) %>%
  gather(avg_DOC, avg_SPOC, key = "sample", value = "value") 
  
levels(SPOC_DOC_data$Reach) <- c("Treatment", "Reference")

FH_streamOC <- SPOC_DOC_data %>%
  filter(Site == "FH") %>%
  ggplot(aes(Date, value, fill = sample)) +
  geom_bar(position="dodge", stat="identity") +
  facet_grid(rows = vars(Reach)) +
  geom_point(aes(y = Total_OC, color = "Total Suspended Organic Carbon")) +
  labs(title = "Fish Creek", 
       x = NULL, 
       y = expression(Suspended~Organic~Carbon~(mg~C~L^-1))) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(colour = "black")) +
  theme(axis.text = element_text(size = 12)) +
  scale_fill_manual(name = "Particulate \nCarbon", 
                    values = c("#F8766D", "#00BFC4"),labels = c("DOC", "SPOC")) +
  scale_color_manual(name = NULL, values = c("Total Suspended Organic Carbon" = "black")) +
  scale_x_discrete(expand = c(0,0), guide = guide_axis(angle = 45)) +
  guides(fill = guide_legend(override.aes = list(shape = NA)))
FH_streamOC

# labels = c("Coarse", "Fine")
