# Step_5_Modeling.R
# Contains all code used to construct models for the six selected analytes (PFOS, PFNA, PFDA, PFUnA, PFDoA, PFTrDA)

# Written by Peter Martin
# Created December 13, 2024
# Finalized

# Working directory
setwd("~/Desktop/Publications/Leyerle Martin et al., 2025")

# # Packages for linear and mixed model regression code
# library(broom)
# library(fitdistrplus)
# library(emmeans)
# library(DescTools)
# library(ggpp) # recognizes formulations like geom = "table"
# library(ggResidpanel)
# library(lmtest)
# library(devtools)
# 
# library(Matrix)
# library(TMB)
# library(glmmTMB)
# 

## Packages
# Data formatting and combining
# library(stringr)
# library(Hmisc) # capitalize function
# library(plyr)
# library(dplyr)
library(tidyverse)
library(lattice)
library(caret)


# Modeling
library(nlme)
library(mgcv)
library(DHARMa)

# Visualizing model output and calculating results
library(gratia)
library(carData)
library(car)
library(emmeans)


# library(mapview)
# library(leafsync)
# # library(ggplot2)
# # library(rnaturalearth)

# # Imputation
# library(MASS)
# library(survival)
# library(NADA)
# library(truncnorm)
# library(zCompositions)
#
# #
# Assigning watershed designations from the GL watershed shapefile
# library(sf)
# library(s2)
# #


## Functions
source("PFAS_Review_supportingFunctions.R") # Load supporting functions

## Loading in data frame created from Step 1
final_imputed_data <- read.csv("Step_3_full_imputed_df.csv",header = TRUE)


###### Pre-Modeling Formatting ################################################
# Converting columns to proper format (i.e., ordered factors)
final_imputed_data$Waterbody<-factor(final_imputed_data$Waterbody,
                                     levels = c("Lake Superior",
                                                "Lake Michigan",
                                                "Lake Huron",
                                                "Lake Erie",
                                                "Lake Ontario"))

final_imputed_data$Waterbody_Type<-factor(final_imputed_data$Waterbody_Type,
                                          levels = c("Inland waters",
                                                     "Connecting channel",
                                                     "Lake"))

final_imputed_data$Composite<-factor(final_imputed_data$Composite)

final_imputed_data$Trophic_Level<-factor(final_imputed_data$Trophic_Level,
                                         levels = c("Primary Producer",
                                                    "Primary Consumer",
                                                    "Secondary Consumer",
                                                    "Tertiary Consumer",
                                                    "Quaternary Consumer",
                                                    "Piscivorous/Insectivorous Bird",
                                                    "Apex Predator"))

final_imputed_data$Class<-factor(final_imputed_data$Class,
                                 levels = c("Algae","Plantae (Magnoliopsida)",
                                            "Annelida","Bivalvia","Gastropoda",
                                            "Zooplankton",
                                            "Shrimp, water fleas, and allies",
                                            "Insecta","Astacoidea","Amphibia",
                                            "Pisces","Reptilia","Aves","Mammalia"))

final_imputed_data$Revised_Tissue<-factor(final_imputed_data$Revised_Tissue,
                                          levels = c("Muscle",
                                                     "Whole organism homogenate",
                                                     "Misc. Tissue",
                                                     "Liver","Blood","Eggs"))

final_imputed_data <- final_imputed_data %>% group_by(Waterbody) %>% 
  mutate(Water_Level_sc = scale(Water_Level))
# Idea for this line of code from 
# https://stackoverflow.com/questions/41761018/scale-all-values-depending-on-group

final_imputed_data %>% group_by(Waterbody) %>% 
  summarise(mean(Water_Level_sc),sd(Water_Level_sc)) # Check that the scaling was performed correctly (mean should be 0 and standard deviation should be 1)


###### Full Data Frame Modeling ###############################################
# We decided to focus analysis on the following contaminants that all had >80% Detection Frequency:
# PFOS, PFNA, PFDA, PFUnA, PFDoA, PFTrDA

# Iterative modeling led us to choose the following variables for these six GAMs

# Fixed Effects: Waterbody, Waterbody_Type, Trophic_Level, Revised_Tissue, Water_Level and Sampling_Year
# Random Effects: Composite

######################### Full PFOS model ######################################
PFOS<-final_imputed_data[is.na(final_imputed_data$PFOS)==FALSE,]

car::qqPlot(log(PFOS$PFOS,base = 10))

# use multiple threads for fitting
ctrl <- gam.control(nthreads = 4)

# Hold out validation (80:20 method)
# From https://stackoverflow.com/questions/22972854/how-to-implement-a-hold-out-validation-in-r
set.seed(100) 
# 26 or 226 runs in about 51 minutes (significant temporal autocorrelation, marginally significant deviations from 1:1)
# 36 runs in about 46 minutes (significant DHARMA output, potentially significant patterns in residuals)


PFOS_in_train<-createDataPartition(PFOS$PFOS, p = 4/5, list = FALSE)
PFOS_train_set<-PFOS[PFOS_in_train,]
PFOS_test_set<-PFOS[-PFOS_in_train,]

# ########################## Load the model for analysis ########################
 load("full_PFOS_gam.Rdata")
# ###############################################################################

# About 50 minutes
system.time(
  full_PFOS_gam <- gam(log10(PFOS) ~ Waterbody + Waterbody_Type
                       + s(Water_Level_sc,by = Waterbody,k=20)
                       + Trophic_Level
                       + Revised_Tissue
                       + s(Sampling.Year,by = Waterbody,k=20)
                       + s(Longitude,Latitude,k=525)
                       + s(Composite,bs="re"),
                       data = PFOS_train_set,
                       method = "REML",
                       family = scat(link = "identity"),
                       control = ctrl)
)
# When we set k in s(Longitude,Latitude,k=525) higher than 525 (at 642, 
# the maximum), there isn't a noticeable increase in edf, and the model shows
# signs of overfitting (does not accurately predict observed values)
# And computation time is much longer (about 50 minutes)


# See email conversation of July 19th, 2024
# Model diagnostics and output support the inclusion of the 15 data points in 
# our models

# Adding in Class as a random effect (final option for including Class)
# yielded more variable estimates of Trophic_Level levels without changing 
# the magnitude of estimates appreciably (results become less coherent with 
# theory — Apex predator estimates are not significantly different from 
# secondary consumer estimates)
#                         df      AIC
# full_PFOS_gam       107.5607 2050.497
# class_full_PFOS_gam 117.5529 1929.754
# DHARMa diagnostics don't really change (KS test slightly less significant but
# increased quantile deviations), and both QQ plot and Observed vs Fitted plot
# in appraise() get slightly worse

# Thus, I think we can refute the idea that Class should be included in the model

# system.time(
#   class_full_PFOS_gam <- gam(log(PFOS,base = 10) ~ Waterbody + Waterbody_Type
#                        + s(Water_Level_sc,by = Waterbody)
#                        + Trophic_Level
#                        + Revised_Tissue
#                        + s(Sampling.Year,by = Waterbody)
#                        + s(Longitude,Latitude)
#                        + s(Composite,bs="re")
#                        + s(Class,bs="re"),
#                        data = PFOS_train_set,
#                        method = "REML",
#                        family = scat(link = "identity"),
#                        control = ctrl)
# )
# 
# dwplot(list(no_class = full_PFOS_gam, 
#             class = class_full_PFOS_gam), 
#        effects = "fixed")

# Outliers (23 at the two margins (n = 1981), p-value = 0.07543)
# No outliers according to Cook's distance
testOutliers(full_PFOS_gam)
hist(cooks.distance(full_PFOS_gam)) 
# Values greater than 1 indicate a potential problem
# perhaps a better test for "outliers" than testOutliers() functionality
# https://github.com/florianhartig/DHARMa/issues/171
# https://www.rdocumentation.org/packages/DHARMa/versions/0.4.6/topics/testOutliers

residualPlot(full_PFOS_gam)
# residualPlot(na_out_full_PFOS_gam)
# No real heteroscedacity that's observable, pink line doesn't deviate
# substantially from the theoretical horizontal line 

# get data
refit_PFOS <- full_PFOS_gam$model
# make residuals column
refit_PFOS$resid <- residuals(full_PFOS_gam,type="deviance")
# fit a model (same model)
PFOS_resid_fit <- gam(resid ~ Waterbody + Waterbody_Type
                      + s(Water_Level_sc,by = Waterbody,k=20)
                      + Trophic_Level
                      + Revised_Tissue
                      + s(Sampling.Year,by = Waterbody,k=20)
                      + s(Longitude,Latitude,k=525)
                      + s(Composite,bs="re"),
                      family=gaussian(), data=refit_PFOS, method="REML")
summary(PFOS_resid_fit)

ggplot(refit_PFOS, aes(x = Waterbody_Type, y = resid)) + 
  geom_point() + 
  geom_hline(yintercept = 0, color='red', linetype='dashed')


plot(simulateResiduals(full_PFOS_gam),quantreg=T)
# plot(simulateResiduals(na_out_full_PFOS_gam),quantreg=T)
# Some deviation according to KS test, as well as a slight deviation in 
# quantiles, but nothing to really be concerned about

appraise(full_PFOS_gam)
# appraise(na_out_full_PFOS_gam)
# Histogram looks great, QQ plot doesn't reveal any glaring problems (only a bit
# of the left tail falls outside the confidence intervals), and observed vs
# fitted values looks great!

########## Test for Spatial Autocorrelation: N.S. 
# p-value = 0.9534
PFOS_test_set$coords <- paste(PFOS_test_set$Longitude,", ",
                              PFOS_test_set$Latitude)
coords <- c(unique(PFOS_test_set$coords))
x_unique <- c(str_extract(coords, "^.+(?=,)"))
y_unique <- c(str_extract(coords, "(?<=, ).+$"))

res<-recalculateResiduals(simulateResiduals(full_PFOS_gam),
                          group = PFOS_test_set$coords)
testSpatialAutocorrelation(res, x = x_unique, y = y_unique)

# na_out_PFOS_test_set$coords <- paste(na_out_PFOS_test_set$Longitude,", ",
#                                      na_out_PFOS_test_set$Latitude)
# coords <- c(unique(na_out_PFOS_test_set$coords))
# x_unique <- c(str_extract(coords, "^.+(?=,)"))
# y_unique <- c(str_extract(coords, "(?<=, ).+$"))
# 
# res<-recalculateResiduals(simulateResiduals(na_out_full_PFOS_gam),
#                           group = na_out_PFOS_test_set$coords)
# testSpatialAutocorrelation(res, x = x_unique, y = y_unique)

############ Test for Temporal Autocorrelation: N.S.
# p-value = 0.5369
res = recalculateResiduals(simulateResiduals(full_PFOS_gam), 
                           group = PFOS_test_set$Sampling.Year)
testTemporalAutocorrelation(res, time = unique(PFOS_test_set$Sampling.Year))

# res = recalculateResiduals(simulateResiduals(na_out_full_PFOS_gam), 
#                            group = na_out_PFOS$Sampling.Year)
# testTemporalAutocorrelation(res, time = unique(PFOS$Sampling.Year))

############## Testing effectiveness of model in predicting test set
PFOS_results<-cbind(PFOS_test_set,as.data.frame(predict(full_PFOS_gam,
                                                        PFOS_test_set,
                                                        type="response",
                                                        se.fit=TRUE)))

PFOS_line<-lm(log(PFOS, base = 10) ~ (fit),
              data = PFOS_results)
summary(PFOS_line)

# Clever piece of code from 
# https://stackoverflow.com/questions/33060601/test-if-the-slope-in-simple-linear-regression-equals-to-a-given-constant-in-r
# Significance test (using two emmeans functions) that test the 1:1 assumption
# of predicted vs observed
# https://cran.r-project.org/web/packages/emmeans/vignettes/confidence-intervals.html
fit.emt <- emtrends(PFOS_line, ~1, var="fit")
fit.emt
test(fit.emt, null=1)

# OR (per https://stats.stackexchange.com/questions/225956/how-do-i-test-hypothesis-of-slope-1-and-intercept-0-for-observed-vs-predict)
PFOS_results$dif<-(log(PFOS_results$PFOS,base = 10)-PFOS_results$fit)
one_to_one<-lm(dif ~ fit,
               data = PFOS_results)
summary(one_to_one)
# Both validate the assumption of 1:1

# Plot of predicted (X) vs observed (y) on the log scale
ggplot(data=PFOS_results, aes(x=(fit), y=log(PFOS,base = 10),color=Waterbody)) +
  geom_point() +
  geom_pointrange(aes(xmin=fit-se.fit, xmax=fit+se.fit)) +
  geom_abline(slope=1, intercept=0)

##################### Save the model for future work ##########################
save(full_PFOS_gam,file = "c_full_PFOS_gam.Rdata")
###############################################################################


######################### Full PFNA Model #####################################
PFNA<-final_imputed_data[is.na(final_imputed_data$PFNA)==FALSE,]
car::qqPlot(log(PFNA$PFNA,base = 10))

# use multiple threads for fitting
ctrl <- gam.control(nthreads = 4)

# Hold out validation (80:20 method)
# From https://stackoverflow.com/questions/22972854/how-to-implement-a-hold-out-validation-in-r
set.seed(16)
PFNA_in_train<-createDataPartition(PFNA$PFNA, p = 4/5, list = FALSE)
PFNA_train_set<-PFNA[PFNA_in_train,]
PFNA_test_set<-PFNA[-PFNA_in_train,]

########################## Load the model for analysis ########################
load("full_PFNA_gam.Rdata")
###############################################################################

# About 6 minutes and 40 seconds
system.time(
  full_PFNA_gam <- gam(log10(PFNA) ~ Waterbody + Waterbody_Type
                       + s(Water_Level_sc,by = Waterbody)
                       + Trophic_Level
                       + Revised_Tissue
                       + s(Sampling.Year,by = Waterbody,k=20)
                       + s(Longitude,Latitude,k=350)
                       + s(Composite,bs="re"),
                       data = PFNA_train_set,
                       method = "REML",
                       family = scat(link = "identity"),
                       control = ctrl)
)
# s(Longitude,Latitude,k=350) maximizes the goodness of fit (indicated by k-index
# higher p-value, appraise() and DHARMa) while avoiding long computation times
# k > 350 (e.g., 400, 450, 500) did not change appreciably in the edf (increase 
# of 2-3)


# Outliers (11 at the two margins (n = 1684), p-value = 0.6791)
# No outliers according to Cook's distance
testOutliers(full_PFNA_gam)
hist(cooks.distance(full_PFNA_gam)) 

residualPlot(full_PFNA_gam)
# residualPlot(na_out_full_PFOS_gam)
# No real heteroscedacity that's observable, pink line doesn't deviate
# substantially from the theoretical horizontal line 

# get data
refit_PFNA <- full_PFNA_gam$model
# make residuals column
refit_PFNA$resid <- residuals(full_PFNA_gam,type="deviance")
# fit a model (same model)
PFNA_resid_fit <- gam(resid ~ Waterbody + Waterbody_Type
                      + s(Water_Level_sc,by = Waterbody)
                      + Trophic_Level
                      + Revised_Tissue
                      + s(Sampling.Year,by = Waterbody,k=20)
                      + s(Longitude,Latitude,k=350)
                      + s(Composite,bs="re"),
                      family=gaussian(), data=refit_PFNA, method="REML")
summary(PFNA_resid_fit)

ggplot(refit_PFNA, aes(x = Waterbody_Type, y = resid)) + 
  geom_point() + 
  geom_hline(yintercept = 0, color='red', linetype='dashed')


plot(simulateResiduals(full_PFNA_gam),quantreg=T)
# plot(simulateResiduals(na_out_full_PFOS_gam),quantreg=T)
# No outliers and minimal quantile deviations, significance for KS test
# Nothing to be concerned about

appraise(full_PFNA_gam)
# appraise(na_out_full_PFOS_gam)
# Histogram looks great, QQ plot doesn't reveal any glaring problems (only a bit
# of the right tail falls outside the confidence intervals), and observed vs
# fitted values looks great!

########## Test for Spatial Autocorrelation: N.S. 
# p-value = 0.1704
PFNA_test_set$coords <- paste(PFNA_test_set$Longitude,", ",
                              PFNA_test_set$Latitude)
coords <- c(unique(PFNA_test_set$coords))
x_unique <- c(str_extract(coords, "^.+(?=,)"))
y_unique <- c(str_extract(coords, "(?<=, ).+$"))

res<-recalculateResiduals(simulateResiduals(full_PFNA_gam),
                          group = PFNA_test_set$coords)
testSpatialAutocorrelation(res, x = x_unique, y = y_unique)

############ Test for Temporal Autocorrelation: N.S.
# p-value = 0.3425
res = recalculateResiduals(simulateResiduals(full_PFNA_gam), 
                           group = PFNA_test_set$Sampling.Year)
testTemporalAutocorrelation(res, time = unique(PFNA_test_set$Sampling.Year))


############## Testing effectiveness of model in predicting test set
PFNA_results<-cbind(PFNA_test_set,as.data.frame(predict(full_PFNA_gam,
                                                        PFNA_test_set,
                                                        type="response",
                                                        se.fit=TRUE)))

PFNA_line<-lm(log(PFNA, base = 10) ~ (fit),
              data = PFNA_results)
summary(PFNA_line)

# Clever piece of code from 
# https://stackoverflow.com/questions/33060601/test-if-the-slope-in-simple-linear-regression-equals-to-a-given-constant-in-r
# Significance test (using two emmeans functions) that test the 1:1 assumption
# of predicted vs observed
# https://cran.r-project.org/web/packages/emmeans/vignettes/confidence-intervals.html
fit.emt <- emtrends(PFNA_line, ~1, var="fit")
fit.emt
test(fit.emt, null=1)

# OR (per https://stats.stackexchange.com/questions/225956/how-do-i-test-hypothesis-of-slope-1-and-intercept-0-for-observed-vs-predict)
PFNA_results$dif<-(log(PFNA_results$PFNA,base = 10)-PFNA_results$fit)
one_to_one<-lm(dif ~ fit,
               data = PFNA_results)
summary(one_to_one)
# Both validate the assumption of 1:1

# Plot of predicted (X) vs observed (y) on the log scale
ggplot(data=PFNA_results, aes(x=(fit), y=log(PFNA,base = 10),color=Trophic_Level)) +
  geom_point() +
  geom_pointrange(aes(xmin=fit-se.fit, xmax=fit+se.fit)) +
  geom_abline(slope=1, intercept=0)

##################### Save the model for future work ##########################
save(full_PFNA_gam,file = "full_PFNA_gam.Rdata")
###############################################################################


