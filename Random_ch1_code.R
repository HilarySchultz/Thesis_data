#### Benthic ####

# getOutliers(Benthic_data$Total_BPOC_Mass_per_Area, distribution = "lognormal")
#### Outlier removal ####
# We want to transform the data into a normal distribution 
# Benthic_data_transformations <- Benthic_data %>%
#   mutate(log_trans_BPOC = log(Total_BPOC_Mass_per_Area), 
#          reciprocal_root = -1/(sqrt(Total_BPOC_Mass_per_Area)),
#          reciprocal = -1/Total_BPOC_Mass_per_Area,
#          reciprocal_square = -1/(Total_BPOC_Mass_per_Area^2),
#          square = Total_BPOC_Mass_per_Area^2,
#          square_root = sqrt(Total_BPOC_Mass_per_Area))
# 
# hist(Benthic_data_transformations$log_trans_BPOC)
# hist(Benthic_data_transformations$reciprocal_root)
# hist(Benthic_data_transformations$reciprocal)
# hist(Benthic_data_transformations$reciprocal_square)
# hist(Benthic_data_transformations$square)
# hist(Benthic_data_transformations$square_root)
# 
# 
# skewness(Benthic_data$Total_BPOC_Mass_per_Area)[1] # 2.93
# skewness(Benthic_data_transformations$log_trans_BPOC)[1] # -0.42
# skewness(Benthic_data_transformations$reciprocal_root)[1] # -2.57
# skewness(Benthic_data_transformations$reciprocal)[1] # -10.15
# skewness(Benthic_data_transformations$reciprocal_square)[1] # 6.5 
# skewness(Benthic_data_transformations$square)[1] 
# skewness(Benthic_data_transformations$square_root)[1]
skewness(Benthic_data$Total_BPOC_Mass_per_Area)[1] # measure of skewness 
mc(Benthic_data$Total_BPOC_Mass_per_Area)[1] # medcouple measurement - needs to be less than 0.6 and greater than -0.6
# mc = 0.40
adj_boxplot <- boxB(Benthic_data$Total_BPOC_Mass_per_Area, method = "adjbox")
Benthic_outliers <- adj_boxplot$outliers # row numbers containing outliers
Benthic_data$Total_BPOC_Mass_per_Area[adj_boxplot$outliers] # tells you the values that are excluded

# benthic_outliers <- Benthic_data %>%
#   rstatix::identify_outliers("Total_BPOC_Mass_per_Area")
# 
# 
# mc(Benthic_data_transformations$log_trans_BPOC)[1] # medcouple measurement - needs to be less than 0.6 and greater than -0.6
# # mc = 0.40
# adj_boxplot <- boxB(Benthic_data_transformations$log_trans_BPOC, method = "adjbox")
# Benthic_outliers <- adj_boxplot$outliers # row numbers containing outliers
# Benthic_data_transformations$log_trans_BPOC[adj_boxplot$outliers] # tells you the values that are excluded
# 
# boxplot(log(Benthic_data$Total_BPOC_Mass_per_Area))
# 
# hist(log(Benthic_data$Total_BPOC_Mass_per_Area))

# benthic_outliers <- Benthic_data %>%
#   rstatix::identify_outliers(Total_BPOC_Mass_per_Area) %>%
#   filter(is.extreme == T)
# 
# Benthic_data <- anti_join(Benthic_data, benthic_outliers)


# finalbenthicmodel_log <- lmer(log(Total_BPOC_Mass_per_Area) ~ Date/Site/Segment + (1|Reach) + (1|Location),
#                            data = Benthic_data)
# benthic_fit_log <- r.squaredGLMM(finalbenthicmodel_log) # 0.54
# 
# par(mfrow = c(2,2))
# 
# plot(finalbenthicmodel_log)
# Anova(finalbenthicmodel_log)
# Area <- Average_widths %>%
#   mutate(Segment_length = case_when(Site == "FH" & Segment == "BDA" ~ FH_treat,
#                                     Site == "FH" & Segment == "REF" ~ FH_ref,
#                                     Site == "LP" & Segment == "BDA" ~ LP_treat,
#                                     Site == "LP" & Segment == "REF" ~ LP_ref,
#                                     Site == "TP" & Segment == "BDA" ~ TP_treat,
#                                     Site == "TP" & Segment == "REF" ~ TP_ref))
# 
# Benthic_data <- full_join(Benthic_data, Area) %>%
#   mutate(streambed_area = avg_width * Segment_length) %>%
#   mutate(benthic_per_m = avg_width * 100 * Total_BPOC_Mass_per_Area)
# Checking to see that the model is not overdispersed
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

overdisp_fun(finalbenthicmodel) # Model is not overdispersed
#### BOXPOLT ####
ggplot() +
  geom_boxplot(data = Benthic_data, aes(x = Segment, y = Total_BPOC_Mass_per_Area, fill = Segment)) +
  # geom_point(data = benthic_cld, aes(x = Segment, y = response), size = 1, shape = 8,
  #            color = "darkred") +
  geom_text(data = benthic_cld, aes(x = Segment, y = response, label= .group,
                                    vjust = -3, hjust = -1.5),
            size = 8, position = position_dodge(0.5), color = "black") +
  geom_text(aes()) +
  scale_fill_manual(name = "Segment", labels = c("Treatment", "Reference"), values = c("#0072B2", "#009E73")) +
  # scale_fill_brewer(palette = "Spectral") +
  labs(title = NULL, 
       x = NULL, 
       y = expression(Benthic~Particulate~Organic~Carbon~Pools~(g~C~m^-2))) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(colour = "black", size = 12),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  facet_grid(Site~Date) 

benthic_emmeans$Date <- ordered(benthic_emmeans$Date, levels = c("June", "July", "August"))

#### RIBBON PLOT 
ggplot(data = benthic_emmeans) +
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
  # scale_y_log10(limits = c(1,1e4)) +
  # annotation_logticks(sides = "l") +
  geom_text(data = benthic_cld, aes(x = Date, y = response, label= .group,
                                    vjust = -1.5, hjust = 0.5),
            size = 6, position = position_dodge(0.5), color = "black") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(guide = guide_axis(angle = 45), 
                   # labels = c("June","June","June", "July","July", "August", "August", "September"),
                   expand = c(0,0)) +
  scale_color_manual(name = "Segment", labels = c("Treatment", "Reference"), values = c("#0072B2", "#009E73")) +
  scale_fill_manual(name = "Segment", labels = c("Treatment", "Reference"), values = c("#0072B2", "#009E73")) + 
  labs(title = NULL, 
       x = NULL, 
       y = expression(Benthic~Particulate~Organic~Carbon~(g~m^-2))) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(colour = "black", size = 12),
        panel.spacing.x = unit(1, "lines"),
        axis.title = element_text(size = 12),
        axis.title.y = element_text(size = 14)) +
  # axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5)) +
  facet_grid(rows = vars(Site))
# This code is to change the orientation of the letter display
# rep(-1.8, 5), 2, rep(-1.8, 5)

#### CART 
library(rpart)
library(rpart.plot)
rep_tree <- rpart(Total_BPOC_Mass_per_Area~Site + Reach +Location +Replicate,
                  data = benthic_means, method = "class",
                  control=rpart.control(minsplit=10,cp=0.0001))
plot(rep_tree)
text(rep_tree)
print(rep_tree)
printcp(rep_tree)
plotcp(rep_tree)
rep_tree_prune <- prune(rep_tree, cp = 0.05)
print(rep_tree_prune)
plot(rep_tree_prune, main = "")
text(rep_tree_prune, pretty=0,use.n=T,xpd=T) 
rpart.plot(rep_tree, extra = 9)

#### DOC & SPOC ####
DOC_data_transformations <- DOC_data %>%
  mutate(log_trans = log(Conc_ppm), 
         reciprocal_root = -1/(sqrt(Conc_ppm)),
         reciprocal = -1/Conc_ppm,
         reciprocal_square = -1/(Conc_ppm^2),
         square = Conc_ppm^2,
         square_root = sqrt(Conc_ppm))

hist(DOC_data_transformations$log_trans)
hist(DOC_data_transformations$reciprocal_root)
hist(DOC_data_transformations$reciprocal)

skewness(DOC_data_transformations$log_trans)[1] # 0.54
skewness(DOC_data_transformations$reciprocal_root)[1] # 0.033
skewness(DOC_data_transformations$reciprocal)[1] # -0.50

qqPlot(DOC_data_transformations$log_trans)
qqPlot(DOC_data_transformations$reciprocal_root)
qqPlot(DOC_data_transformations$reciprocal)




log_outliers <- DOC_data_transformations %>%
  rstatix::identify_outliers("log_trans") %>%
  filter(is.extreme == TRUE)

reciprocal_outliers <- DOC_data_transformations %>%
  rstatix::identify_outliers("reciprocal_root") %>%
  filter(is.extreme == TRUE)
DOC_data <- DOC_data %>%
  mutate(log_DOC = log(Conc_ppm))
qqPlot(DOC_data$Conc_ppm, distribution = "gamma")

skewness(DOC_data$log_DOC)[1] # measure of skewness 
mc(DOC_data$log_DOC)[1] # medcouple measurement - needs to be less than 0.6 and greater than -0.6
# mc = 0.40
adj_boxplot <- boxB(DOC_data$log_DOC, method = "adjbox")
doc_outliers <- adj_boxplot$outliers # row numbers containing outliers
DOC_data$log_DOC[adj_boxplot$outliers] # tells you the values that are excluded


mod <- lm(Conc_ppm~., data = DOC_data)
cooksd <- cooks.distance(mod)
plot(cooksd, pch = "*", cex = 2)
abline(h = 4*mean(cooksd, na.rm = T), col = "red")
text(x = 1:length(cooksd) + 1, y = cooksd, labels = ifelse(cooksd > 4 * mean(cooksd, na.rm=T),names(cooksd),""), col="red")

car::outlierTest(mod)
outlierPlot(DOC$Conc_ppm)
skewness(DOC_data$Conc_ppm)[1] # measure of skewness 
mc(DOC_data$Conc_ppm)[1] # medcouple measurement - needs to be less than 0.6 and greater than -0.6
# mc = 0.37
adj_boxplot_doc <- boxB(DOC_data$Conc_ppm, method = "adjbox")
DOC_outliers <- adj_boxplot_doc$outliers # row numbers containing outliers
DOC_data$Conc_ppm[adj_boxplot_doc$outliers] # tells you the values that are excluded
# Seeing if the data meets the normal distribution assumption
# shapiro.test(log(DOC_data$Conc_ppm)) # W = 0.85684, p-value = 0.00000000003756
# qqPlot(log(DOC_data$Conc_ppm)) # This tests/shows if the data is normally distributed
# small p-value implies that the distribution of the data is drastically special
# from the normal distribution so we don't count on normality. 

# It is not so we need to use the Levene test to see if there are differences
# between the tested sample variances. 
# levenetestdoc <- leveneTest(Conc_ppm ~ Date*Site*Segment, 
# data = DOC_data)
# p<0.05 says that the variances are relatively equal so we can use Tukey as post-hoc
# doc_emm_contrast1 <- doc_emm_contrast %>%
#   mutate(p.bon = p.adjust(p.value, method = "bonferroni"),
#          p.holm = p.adjust(p.value, method = "holm"),
#          p.sig = if_else(p.holm < signif(0.05,5), "*", ""))
# doc_emm <- emmeans(finaldocmodel, ~ Reach|Site|Date,
#                        type = "response",
#                    nesting = "Date %in% Site, Reach %in% (Site*Date)")
# doc_emm_sum <- summary(doc_emm)
# doc_pair_sum$p.value <- p.adjust(doc_pair_sum$p.value, 
#                              method = "bonferroni")
# doc_p_adjust <- pairwise.t.test(doc_pair_sum$p.value, doc_pair_sum$contrast, p.adjust.method="bonferroni")

# SPOC_data_transformations <- SPOC_data %>%
#   mutate(log_trans = log(SPOC), 
#          reciprocal_root = -1/(sqrt(SPOC)),
#          reciprocal = -1/SPOC,
#          reciprocal_square = -1/(SPOC^2),
#          square = SPOC^2,
#          square_root = sqrt(SPOC))
# 
# hist(SPOC_data_transformations$log_trans)
# hist(SPOC_data_transformations$reciprocal_root)
# hist(SPOC_data_transformations$reciprocal)
# 
# skewness(SPOC_data_transformations$log_trans)[1] # 0.54
# skewness(SPOC_data_transformations$reciprocal_root)[1] # 0.033
# skewness(SPOC_data_transformations$reciprocal)[1] # -0.50
# 
# qqPlot(SPOC_data_transformations$log_trans)
# qqPlot(SPOC_data_transformations$reciprocal_root)
# qqPlot(SPOC_data_transformations$reciprocal)
# 
# spoc_outliers_trans <- SPOC_data_transformations %>%
#   rstatix::identify_outliers("log_trans") 


# Removing outliers
# skewness(SPOC_data$SPOC) # measure of skewness 
# mc(SPOC_data$SPOC)[1] # medcouple measurement - needs to be less than 0.6 and greater than -0.6
# # mc = 0.37
# adj_boxplot_spoc <- boxB(SPOC_data$SPOC, method = "adjbox", k = 2)
# SPOC_outliers <- adj_boxplot_spoc$outliers # row numbers containing outliers
# SPOC_data$SPOC[adj_boxplot_spoc$outliers] # tells you the values that are excluded
# 
# boxplot(log(SPOC_data$SPOC))

# SPOC_data <- SPOC_data %>%
#   slice(-SPOC_outliers) 
# SPOC_summary <- SPOC_data %>% summarySE(groupvars = c("Site", "Date", "Reach"), measurevar = "SPOC") %>%
#     rename(stdev = sd, stderror = se)
# names(SPOC_summary)

# SPOC_data[,1:5] <- lapply(SPOC_data[,1:5], as.factor)

# SPOC_data$Date <- ordered(SPOC_data$Date,
#                          levels = c("6/14/2021", "6/28/2021","7/12/2021", "7/26/2021",
#                                     "8/9/2021", "8/25/2021", "9/7/2021"))
# shapiro.test(SPOC_data)
# shapiro.test(SPOC_data$SPOC) # W = 0.96985, p-value = 0.007091
# qqPlot(SPOC_data$SPOC)
# Not normally distributed 

# levene_spoc <- SPOC_data %>%
#   group_by(Date, Site) %>%
#   levene_test(SPOC~Reach)
# 
# levenetestspoc <- leveneTest(SPOC ~ Date/Site/Reach, 
#                             data = SPOC_data)
# 
# pwc <- SPOC_data %>%
#   group_by(Date, Site) %>%
#   pairwise_t_test(
#     SPOC ~ Reach, paired = FALSE,
#     p.adjust.method = "bonferroni"
#   ) 
# spoc_emm_contrast <- spoc_emm_contrast %>%
#   mutate(p.bon = p.adjust(p.value, method = "bonferroni"),
#          p.holm = p.adjust(p.value, method = "holm"),
#          p.sig = if_else(p.value < signif(0.05,5), "*", ""))
# 
# spoc_emm1.1 <- spoc_emm$contrasts %>%
#   rbind() %>%
#   summary(infer = TRUE) %>%
#   p.sig = if_else(p.value <= 0.05, "*", "")
# 
# contrast_df <- spoc_emm$contrasts %>%
#   rbind() %>%
#   as.data.frame() %>%
#   mutate(p.sig = if_else(p.value <= 0.05, "*", ""))
