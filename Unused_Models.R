# A respository for all unused models in Thesis

#### DOC ####
docmodeltest1 <- glmer(Conc_ppm ~ Date*Site*Reach + (1|Replicate),
                       data = DOC_data, 
                       control = glmerControl(optimizer = "bobyqa",
                                              optCtrl = list(maxfun = 100000)),
                       family = Gamma(link = "log"))
plot(docmodeltest1)
Anova(docmodeltest1)
AIC(docmodeltest1) # 360.3364
AICc(docmodeltest1) # 422.9247

# Does not converge
docmodeltest2 <- glmer(Conc_ppm ~ Date+Site+Reach + (1|Replicate),
                       data = DOC_data, 
                       family = Gamma(link = "log"))
plot(docmodeltest2)
Anova(docmodeltest2)
AIC(docmodeltest2) # 385.8355
AICc(docmodeltest2) # 388.7521



# Replicate = random
docmodel <- glmer(Conc_ppm ~ Date + Reach + (1|Site) + (1|Replicate),
                  data = DOC_data, 
                  family = Gamma(link = "log"))
plot(docmodel)
Anova(docmodel)
AIC(docmodel) # 386.9168
AICc(docmodel) # 389.4272

# This is the current working model 3/15/22
# Replicate = NULL
docmodel1 <- glmer(Conc_ppm ~ Date + Reach + (1|Site),
                   data = DOC_data, 
                   family = Gamma(link = "log"))
plot(docmodel1)
Anova(docmodel1)
AIC(docmodel1) # 385.1746
AICc(docmodel1) # 387.3116
summary(docmodel1)

DOC_nested <- glmer(Conc_ppm ~ Date/Reach + (1|Site) + (1|Replicate),
                    data = DOC_data, 
                    family = Gamma(link = "log"))
Anova(DOC_nested)
plot(DOC_nested)
AIC(DOC_nested) # 398.9339
AICc(DOC_nested) # 405.6785
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
summary(DOC_fullnested)
Anova(DOC_fullnested)

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

# This is without Date. Only site is significant but all these models are still singular
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

DOC_datetest4 <- glmer(Conc_ppm ~ Site + Reach + Date + (1|Replicate),
                       data = DOC_data, 
                       family = Gamma(link = "log"))
DOC_datetest5 <- glmer(Conc_ppm ~ Site + Reach + (1|Date) + (1|Replicate),
                       data = DOC_data, 
                       family = Gamma(link = "log"))
DOC_datetest6 <- glmer(Conc_ppm ~ Date/Site + Reach + (1|Replicate),
                       data = DOC_data, 
                       family = Gamma(link = "log"))
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

DOC_fullnested <- glmer(Conc_ppm ~ Date/Site/Reach + (1|Replicate),
                        data = DOC_data,
                        nAGQ = 0,
                        family = Gamma(link = "log"))
Anova(DOC_fullnested)
AIC(DOC_fullnested) # 307.7063
AICc(DOC_fullnested) # 374.2063

# These glms work
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

doctest6 <- glm(Conc_ppm ~ Date/Site + Reach,
                data = DOC_data, 
                family = Gamma(link = "log"))
par(mfrow = c(2,2))
plot(doctest6)
Anova(doctest6)
AIC(doctest6) # 291.6785
AICc(doctest6) # 305.8249
# All other variations end up reducing down to Date:Site + Date

#### SPOC ####
### These are GLMs ###
spocmodeltest1 <- glmer(SPOC ~ Date*Site*Reach + (1|Replicate),
                        data = SPOC_data, 
                        control = glmerControl(optimizer = "bobyqa",
                                               optCtrl = list(maxfun = 100000)),
                        family = Gamma(link = "log"))
plot(spocmodeltest1)
Anova(spocmodeltest1)
AIC(spocmodeltest1) # -129.5816
AICc(spocmodeltest1) # -79.45506
summary(spocmodeltest1)



spocmodel <- glmer(SPOC ~ Date + Reach + (1|Site) + (1|Replicate),
                   data = SPOC_data, 
                   family = Gamma(link = "log"))
plot(spocmodel)
Anova(spocmodel)
AIC(spocmodel) # -58.65812
AICc(spocmodel) # -56.30098
summary(spocmodel)

# Replicate = NULL
# This is the same structure as the DOC model 
spocmodel1 <- glmer(SPOC ~ Date + Reach + (1|Site),
                    data = SPOC_data, 
                    family = Gamma(link = "log"))
plot(spocmodel1)
Anova(spocmodel1)
AIC(spocmodel1) # -57.46977
AICc(spocmodel1) # -55.52287
summary(spocmodel1)


SPOC_nested <- glmer(SPOC ~ Date/Reach + (1|Site) + (1|Replicate),
                     data = SPOC_data, 
                     family = Gamma(link = "log"))
Anova(SPOC_nested)
plot(SPOC_nested)
AIC(SPOC_nested) # -75.63007
AICc(SPOC_nested) # -69.85649
summary(SPOC_nested)

spoctest <- SPOC_nested <- glm(SPOC ~ Date/Site/Reach,
                               data = SPOC_data, 
                               family = Gamma(link = "log"))
par(mfrow = c(2,2))
plot(spoctest)
Anova(spoctest)
AIC(spoctest) # -100.7296
AICc(spoctest) # -52.22545
summary(spoctest)
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
summary(SPOC_nested2)

SPOC_GLMM4 <- glmer(SPOC ~ Date/Site/Reach/(1|Replicate),
                    data = SPOC_data, 
                    family = Gamma(link = "log"))
Anova(SPOC_GLMM4)
plot(SPOC_GLMM4)
AIC(SPOC_GLMM4) # -100.7296
AICc(SPOC_GLMM4) # -49.301
summary(SPOC_GLMM4)

# Response: SPOC
# Chisq Df Pr(>Chisq)    
# Date            155.04  6  < 2.2e-16 ***
#   Date:Site       472.28 14  < 2.2e-16 ***
#   Date:Site:Reach 241.90 21  < 2.2e-16 ***
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
#### BENTHIC ####
#### Exploratory GLMMs ####
# Nesting date within the interactions between site, reach, and location
# Replicate is a random factor.
# This plot looks decent. 
Nested_BPOC_GLMM <- glmer(Total_BPOC_Mass_per_Area ~ Date/(Site*Reach*Location) + (1|Replicate),
                          data = Benthic_data,
                          family = Gamma(link="log"))

plot(Nested_BPOC_GLMM)
Anova(Nested_BPOC_GLMM)
AICc(Nested_BPOC_GLMM) # 1381.576
AIC(Nested_BPOC_GLMM) #1318.368

Nested_BPOC_GLMM1 <- glmer(Total_BPOC_Mass_per_Area ~ Date/Site/Reach + 
                             Date/Location + Date/Reach + Date/Site + Date + (1|Replicate),
                           data = Benthic_data,
                           family = Gamma(link="log"))
plot(Nested_BPOC_GLMM1) # bunched
Anova(Nested_BPOC_GLMM1)
AICc(Nested_BPOC_GLMM1) # 1309.32
AIC(Nested_BPOC_GLMM1) #1298.602

# Testing out different ways of nesting 
test <- glmer(Total_BPOC_Mass_per_Area ~ Site/Reach/Location*Date + (1|Replicate),
              data = Benthic_data, 
              family = Gamma(link="log"))
plot(test)
Anova(test)
AICc(test)  # 1381.576
AIC(test) # 1318.368


test2 <- glmer(Total_BPOC_Mass_per_Area ~ Site/Reach*Date + Site/Reach/Location + (1|Replicate),
               data = Benthic_data, 
               family = Gamma(link="log"))
plot(test2) # bunched
Anova(test2)
AICc(test2)  # 1377.408
AIC(test2) # 1361.036

test2.1 <- glmer(Total_BPOC_Mass_per_Area ~ Site:Reach:Location + Site:Date +
                   Site:Reach + Date + Site + (1|Replicate),
                 data = Benthic_data, 
                 family = Gamma(link="log"))
plot(test2.1)
Anova(test2.1)
AICc(test2.1)  # 1369.942
AIC(test2.1) # 1359.542

# No date
# All interactions are significant *** 
# Cannot be reduced further
test1 <- glmer(Total_BPOC_Mass_per_Area ~ Site/Reach/Location + (1|Replicate),
               data = Benthic_data, 
               family = Gamma(link="log"))
plot(test1)
Anova(test1)
AICc(test1)  # 1404.622
AIC(test1) # 1398.664

# This has a 4-way interaction, all interactions are significant but AIC are higher
test4 <- glmer(Total_BPOC_Mass_per_Area ~ Date/Site/Reach/Location + (1|Replicate),
               data = Benthic_data, 
               family = Gamma(link="log"))
plot(test4)
Anova(test4)
AICc(test4)  # 1443.321
AIC(test) # 1382.521

### Full interaction models ###
Full_BPOC_GLMM <- glmer(Total_BPOC_Mass_per_Area ~ Site*Reach*Location*Date + (1|Replicate),
                        data = Benthic_data, 
                        family = Gamma(link="log"))

plot(Full_BPOC_GLMM)
Anova(Full_BPOC_GLMM)
AICc(Full_BPOC_GLMM) # 1443.321
AIC(Full_BPOC_GLMM) #1382.521

# Reducing
Full_BPOC_GLMM1 <- glmer(Total_BPOC_Mass_per_Area ~ Site:Reach:Date + Site:Reach:Location + 
                           Location:Date + Site:Date + Site:Location + Site:Reach + Date + Location + 
                           Site + (1|Replicate),
                         data = Benthic_data, 
                         family = Gamma(link="log"))
plot(Full_BPOC_GLMM1)
Anova(Full_BPOC_GLMM1)
AICc(Full_BPOC_GLMM1) # 1381.962
AIC(Full_BPOC_GLMM1) #1360.65

# Reducing
Full_BPOC_GLMM1.1 <- glmer(Total_BPOC_Mass_per_Area ~ Site:Reach:Location + Location:Date + Site:Date + Site:Location + Site:Reach + Date + Location + Site + (1|Replicate),
                           data = Benthic_data, 
                           family = Gamma(link="log"))
plot(Full_BPOC_GLMM1.1)
Anova(Full_BPOC_GLMM1.1)
AICc(Full_BPOC_GLMM1.1) # 1372.867
AIC(Full_BPOC_GLMM1.1) # 1358.669

# Reducing
Full_BPOC_GLMM1.2 <- glmer(Total_BPOC_Mass_per_Area ~ Site:Reach:Location + Site:Date + Site:Location + Site:Reach + Date + Location + Site + (1|Replicate),
                           data = Benthic_data, 
                           family = Gamma(link="log"))
plot(Full_BPOC_GLMM1.2)
Anova(Full_BPOC_GLMM1.2)
AICc(Full_BPOC_GLMM1.2) # 1369.942
AIC(Full_BPOC_GLMM1.2) # 1359.542

Full_BPOC_GLMM1.3 <- glmer(Total_BPOC_Mass_per_Area ~ Site:Date + Site:Reach + Date + Location + Site + (1|Replicate),
                           data = Benthic_data, 
                           family = Gamma(link="log"))
plot(Full_BPOC_GLMM1.3)
Anova(Full_BPOC_GLMM1.3)
AICc(Full_BPOC_GLMM1.3) # 1366.942
AIC(Full_BPOC_GLMM1.3) # 1362.542


### Main-effects models ###
# Building models up from simple to more complex with no interactions
b_glm <- glmer(Total_BPOC_Mass_per_Area ~ Reach + (1|Replicate),
               data = Benthic_data,
               family = Gamma(link="log"))

plot(b_glm)
Anova(b_glm)
AICc(b_glm) # 1440.214
AIC(b_glm) # 1439.959

b_glm1 <- glmer(Total_BPOC_Mass_per_Area ~ Reach + Site + (1|Replicate),
                data = Benthic_data,
                family = Gamma(link="log"))

plot(b_glm1)
Anova(b_glm1)
AICc(b_glm1) # 1432.186
AIC(b_glm1) # 1431.644

b_glm2 <- glmer(Total_BPOC_Mass_per_Area ~ Reach + Site + Location + Date + (1|Replicate),
                data = Benthic_data,
                family = Gamma(link="log"))

plot(b_glm2)
Anova(b_glm2)
AICc(b_glm2) # 1406.561
AIC(b_glm2) # 1405.104

test4 <- glmer(Total_BPOC_Mass_per_Area ~ Site + Reach + Date + (1|Location) + (1|Replicate),
               data = Benthic_data, 
               family = Gamma(link="log"))
plot(test4)
Anova(test4)
AICc(test4) 
AIC(test4) 

newtest <- glmer(Total_BPOC_Mass_per_Area ~ Site/Reach*Date + (1|Location) + (1|Replicate),
                 data = Benthic_data, 
                 family = Gamma(link="log"))
plot(newtest) # Bunched in the lower left-hand corner. 
Anova(newtest)
AICc(newtest)  # 1307.7
AIC(newtest) # 1300.906

# All of these interactions and main effects are significant.
# Cannot be reduced further
newtest1 <- glmer(Total_BPOC_Mass_per_Area ~ 
                    Site + Date + Site/Reach + Site/Date +
                    (1|Location) + (1|Replicate),
                  data = Benthic_data, 
                  family = Gamma(link="log"))
plot(newtest1) # Bunched in the lower left-hand corner. 
Anova(newtest1) 
AICc(newtest1)  # 1305.625
AIC(newtest1) # 1302.244

btest <- glmer(Total_BPOC_Mass_per_Area ~ Date/Site + Reach + (1|Location) + (1|Replicate),
               data = Benthic_data, 
               family = Gamma(link="log"))
plot(btest) # Bunched in the lower left-hand corner. 
Anova(btest)
AICc(btest) # 1341.965
AIC(btest) # 1339.437

### Full nested models ###
# All interactions are significant. 
# Cannot be reduced further
newtest2 <- glmer(Total_BPOC_Mass_per_Area ~ Date/Site/Reach + (1|Location) + (1|Replicate),
                  data = Benthic_data, 
                  family = Gamma(link="log"))
plot(newtest2) # Bunched in the lower left-hand corner. 
Anova(newtest2)
AICc(newtest2) # 1307.7
AIC(newtest2) # 1300.906

# All interactions are significant ***
# Cannot be reduced further
test3 <- glmer(Total_BPOC_Mass_per_Area ~ Site/Reach/Location + Date + (1|Replicate),
               data = Benthic_data, 
               family = Gamma(link="log"))
plot(test3)
Anova(test3)
AICc(test3)  # 1307.887
AIC(test3) # 1300.39

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
###
# Nesting some models to see if the interaction improves the model

# Location = fixed, replicate = random
bmodelnest <- glmer(Total_BPOC_Mass_per_Area ~ Date/Reach + Location + (1|Site) + (1|Replicate),
                    data = Benthic_data, 
                    family = Gamma(link="log"))
plot(bmodelnest) # Bunched in the lower left-hand corner. 
Anova(bmodelnest) # All ***
AICc(bmodelnest)  # 1348.44
AIC(bmodelnest) # 1346.632

# Location = random, replicate = random
bmodelnest1 <- glmer(Total_BPOC_Mass_per_Area ~ Date/Reach + (1|Location) + (1|Site) + (1|Replicate),
                     data = Benthic_data, 
                     family = Gamma(link="log"))
plot(bmodelnest1) # Bunched in the lower left-hand corner. 
Anova(bmodelnest1) # All ***
AICc(bmodelnest1)  # 1353.229
AIC(bmodelnest1) # 1351.732

# Location = fixed, replicate = NULL
bmodelnest2 <- glmer(Total_BPOC_Mass_per_Area ~ Date/Reach + Location + (1|Site),
                     data = Benthic_data, 
                     family = Gamma(link="log"))
plot(bmodelnest2) # Bunched in the lower left-hand corner. 
Anova(bmodelnest2) # All ***
AICc(bmodelnest2)  # 1353.344
AIC(bmodelnest2) # 1351.847

# Location = random, replicate = NULL
bmodelnest3 <- glmer(Total_BPOC_Mass_per_Area ~ Date/Reach + (1|Location) + (1|Site),
                     data = Benthic_data, 
                     family = Gamma(link="log"))
plot(bmodelnest3) # Bunched in the lower left-hand corner. 
Anova(bmodelnest3) # All ***
AICc(bmodelnest3)  # 1358.365
AIC(bmodelnest3) # 1357.148


#### SGF ####
cdioxidemodel <- glmer(CO2_flux_g ~ Date + Reach + (1|Site) + (1|Chamber), 
                       data = CO2_fluxes,
                       family = Gamma(link = "log"))
plot(cdioxidemodel)
Anova(cdioxidemodel)
AICc(cdioxidemodel) # 647.1857
AIC(cdioxidemodel) # 645.9574

cdioxidemodel1 <- lmer(CO2_flux_g ~ Date + Reach + (1|Site) + (1|Chamber), 
                       data = CO2_fluxes)
plot(cdioxidemodel1)
Anova(cdioxidemodel1)
AICc(cdioxidemodel1) # 727.0645
AIC(cdioxidemodel1) # 725.8361

Nestedco2glm <- glm(CO2_flux_g ~ Date/Site/Reach, 
                    data = CO2_fluxes, 
                    family = gaussian(link = "identity"))
plot(Nestedco2glm)
Anova(Nestedco2glm)
AICc(Nestedco2glm) # 705.5293
AIC(Nestedco2glm) # 682.9487

Nestedco2glm2 <- glm(CO2_flux_g ~ Date/Site/Reach, 
                     data = CO2_fluxes, 
                     family = Gamma(link = "identity"))
plot(Nestedco2glm2)
Anova(Nestedco2glm2)
AICc(Nestedco2glm2) # 652.2665
AIC(Nestedco2glm2) # 629.6859

# I think this is the best model based on p-values and plots
Nestedco2glm3 <- glm(CO2_flux_g ~ Date/Site/Reach, 
                     data = CO2_fluxes, 
                     family = Gamma(link = "log"))
plot(Nestedco2glm3)
Anova(Nestedco2glm3)
AICc(Nestedco2glm3) # 652.2665
AIC(Nestedco2glm3) # 629.6859


Full_cdioxide_glm <- glm(CO2_flux_g ~ Site*Reach*Date, 
                         data = CO2_fluxes,
                         family = gaussian(link = "identity"))
par(mfrow = c(2,2))
plot(Full_cdioxide_glm)
Anova(Full_cdioxide_glm)
AICc(Full_cdioxide_glm) # -935.0952
AIC(Full_cdioxide_glm) # -958.0416

# Reducing
Full_cdioxide_glm1 <- glm(CO2_flux_g ~ Site:Date + Site:Reach + Date + Site + Reach, 
                          data = CO2_fluxes,
                          family = gaussian(link = "identity"))
par(mfrow = c(2,2))
plot(Full_cdioxide_glm1)
Anova(Full_cdioxide_glm1)
AICc(Full_cdioxide_glm1) # -968.4241
AIC(Full_cdioxide_glm1) #-975.6419

##### GLMMs ####
# Need to use lmer on this b/c glmer doesn't like normal?
# Final model
Nestedco2 <- lmer(CO2_flux_g ~ Date/Site/Reach + (1|Chamber), 
                  data = CO2_fluxes)
plot(Nestedco2)
Anova(Nestedco2)
AICc(Nestedco2) # 706.8167
AIC(Nestedco2) # 683.5291

# Doesn't work
Nestedco2 <- glmer(CO2_flux_g ~ Date/Site/Reach + (1|Chamber), 
                   data = CO2_fluxes,
                   family = gaussian("log"))


Nestedco22 <- glmer(CO2_flux_g ~ Date/Site/Reach + (1|Chamber), 
                    data = CO2_fluxes,
                    nAGQ=0,
                    family = Gamma(link = "log"))
plot(Nestedco22)
Anova(Nestedco22)
AICc(Nestedco22) # 592.6799
AIC(Nestedco22) # 569.0688

Nestedco23 <- glmer(CO2_flux_g ~ Date/Site/Reach + (1|Chamber), 
                    data = CO2_fluxes,
                    nAGQ=0,
                    family = Gamma)
plot(Nestedco23)
Anova(Nestedco23)
AICc(Nestedco23) # 623.417
AIC(Nestedco23) # 599.8059

Full_methane_glm1 <- glmer(CH4_flux_g ~ Site*Reach*Date + (1|Chamber),
                           data = CH4_fluxes,
                           family = gaussian(link = "identity"))

Full_methane_GLMM <- glmer(CH4_flux_g ~ Site*Reach*Date + (1|Chamber),
                           data = CH4_fluxes) # Not working

Full_methane_GLMM1 <- glmer(CH4_flux_g ~ Date/(Site*Reach) + (1|Chamber),
                            data = CH4_fluxes) # Not working

#### Methane GLMs ####
Nested_methane_glm <- lmer(CH4_flux_g ~ Date/Site/Reach + (1|Chamber), 
                           data = CH4_fluxes)
plot(Nested_methane_glm)
Anova(Nested_methane_glm)
AIC(Nested_methane_glm) # -1656.464
AICc(Nested_methane_glm)# -1630.836

# This is the best lmer
Nested_methane_glm1.1 <- lmer(CH4_flux_g ~ Date/Site + Date + (1|Chamber), 
                              data = CH4_fluxes)
plot(Nested_methane_glm1.1)
Anova(Nested_methane_glm1.1)
AIC(Nested_methane_glm1.1) # -1975.843
AICc(Nested_methane_glm1.1)# -1969.547

meth_glm <- glm(CH4_flux_g ~ Date/Site/Reach,
                data = CH4_fluxes,
                family = gaussian("identity"))
plot(meth_glm)
Anova(meth_glm)
AIC(meth_glm) # -2219.12
AICc(meth_glm)# -2194.62

# Best glm
meth_glm1.1 <- glm(CH4_flux_g ~ Date/Site + Date,
                   data = CH4_fluxes,
                   family = gaussian("identity"))
plot(meth_glm1.1)
Anova(meth_glm1.1)
AIC(meth_glm1.1) # -2252.736
AICc(meth_glm1.1)# -2246.932


Full_methane <- glm(CH4_flux_g ~ Site*Reach*Date, 
                    data = CH4_fluxes,
                    family = gaussian(link = "identity"))
par(mfrow = c(2,2))
plot(Full_methane_glm)
Anova(Full_methane_glm)
AICc(Full_methane_glm) # -2025.942
AIC(Full_methane_glm) # -2050.332

test_methane_glm <- glm(CH4_flux_g ~ Reach:Date + Site:Date + Site:Reach + Date + Reach + Site,
                        data = CH4_fluxes,
                        family = gaussian(link = "identity"))
par(mfrow = c(2,2))
plot(test_methane_glm)
Anova(test_methane_glm)
AICc(test_methane_glm) # -2025.942
AIC(test_methane_glm)

test_methane_glm2 <- glm(CH4_flux_g ~ Site:Reach + Date + Reach + Site,
                         data = CH4_fluxes,
                         family = gaussian(link = "identity"))
par(mfrow = c(2,2))
plot(test_methane_glm2)
Anova(test_methane_glm2)
AICc(test_methane_glm2) # -2025.942
AIC(test_methane_glm2)

Full_methane_glm1 <- glm(CH4_flux_g ~ Site + Date, 
                         data = CH4_fluxes,
                         family = gaussian(link = "identity"))
par(mfrow = c(2,2))
plot(Full_methane_glm1)
Anova(Full_methane_glm1)
AICc(Full_methane_glm1) # -2092.876
AIC(Full_methane_glm1) # -2094.071
