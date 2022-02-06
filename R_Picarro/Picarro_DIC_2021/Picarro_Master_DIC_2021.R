library(dplyr) 
library(readr) 
library(ggplot2) 
library(tidyverse)


## GET FILES ##
filelist <-list.files(path = getwd(), pattern="*.csv", full.names = TRUE)
DIC_files <- list()
for(f in 1:length(filelist)) {
  DIC_files[[f]] <- read.csv(filelist[f], header = T, sep = ",", as.is = T)
}


DIC_df <- bind_rows(DIC_files)

DIC_data <- DIC_df %>%
  group_by(DATE, Site, Location, Design) %>%
  distinct(Triplicate_CO2_avg, Triplicate_CH4_avg) %>%
  rename(Reach = Design)

DIC_data %>% 
  ggplot(aes(Site, Triplicate_CO2_avg, fill = Site)) + geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 2, color = "blue", fill = "blue") +
  labs(title = expression(Instream~Dissolved~CO[2]~Concentrations), 
       x = NULL, y = expression(CO[2]~Concentration~(ppm))) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_blank(), 
                     axis.ticks = element_blank(), axis.text = element_text(colour = "black")) +
  scale_fill_brewer(palette = "Spectral", name = "Site") + facet_grid(cols = vars(Location), 
                                                                      rows = vars(Reach))

DIC_data %>% 
  ggplot(aes(Site, Triplicate_CH4_avg, fill = Site)) + geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 2, color = "blue", fill = "blue") +
  labs(title = expression(Instream~Dissolved~CH[4]~Concentrations), 
       x = NULL, y = expression(CH[4]~Concentration~(ppm))) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_blank(), 
                     axis.ticks = element_blank(), axis.text = element_text(colour = "black")) +
  scale_fill_brewer(palette = "Spectral", name = "Site") + facet_grid(cols = vars(Location), 
                                                                      rows = vars(Reach))

DIC_data %>% 
  filter(Site == "FH") %>%
  ggplot(aes(Site, Triplicate_CH4_avg, fill = Site)) + geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 2, color = "blue", fill = "blue") +
  labs(title = expression(Instream~Dissolved~CH[4]~Concentrations), 
       x = NULL, y = expression(CH[4]~Concentration~(ppm))) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_blank(), 
                     axis.ticks = element_blank(), axis.text = element_text(colour = "black")) +
  scale_fill_brewer(palette = "Spectral", name = "Site") + facet_grid(cols = vars(Location), 
                                                                      rows = vars(Reach))

DIC_data %>% 
  ggplot(aes(Reach, Triplicate_CH4_avg, fill = Reach)) + geom_boxplot() +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 2, color = "blue", fill = "blue") +
  labs(title = expression(Instream~Dissolved~CH[4]~Concentrations), 
       x = NULL, y = expression(CH[4]~Concentration~(ppm))) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_blank(), 
                     axis.ticks = element_blank(), axis.text = element_text(colour = "black")) +
  scale_fill_brewer(palette = "Spectral", name = "Reach") + facet_grid(cols = vars(Location))



