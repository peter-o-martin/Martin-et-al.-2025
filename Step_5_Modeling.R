# Step_5_Modeling.R
# Contains all code used to construct models for the six selected analytes (PFOS, PFNA, PFDA, PFUnA, PFDoA, PFTrDA)

# Written by Peter Martin
# Created December 13, 2024
# Finalized December 17, 2024

# Working directory
setwd("~/Desktop/Publications/Leyerle Martin et al., 2025")

## Packages
# Data formatting and combining
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
library(ggplot2)
library(dotwhisker)

## Functions
source("PFAS_Review_supportingFunctions.R") # Load supporting functions

## Loading in data frame created from Step 1
final_imputed_data <- read.csv("Step_3_full_imputed_df.csv",header = TRUE)
validation_supp <- read.csv("Step_4_validation_supp.csv",header = TRUE)

# Load the imputation results (07/12/24) whose values were used to construct the models reported in Leyerle Martin et al., 2025
# Done this way since lrDA() will produce slightly different imputed estimates each time the segment of code (i.e., the Step_2 R script) is run
modeling_concentration_values <- read.csv("full_imputed_df.csv")
modeling_concentration_values <- modeling_concentration_values[,32:41]

final_imputed_data[,32:41] <- modeling_concentration_values

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

###### HERE can access the saved models to explore diagnostics ################
load("full_PFOS_gam.Rdata") 
load("full_PFNA_gam.Rdata")
load("full_PFDA_gam.Rdata")
load("full_PFUnA_gam.Rdata")
load("full_PFDoA_gam.Rdata")
load("full_PFTrDA_gam.Rdata")

###### Full PFOS model ########################################################
# Extract the subset of samples that quantified PFOS concentrations (i.e., remove NAs)
PFOS<-subset(final_imputed_data,!is.na(PFOS))
car::qqPlot(log(PFOS$PFOS,base = 10)) # constructs a QQ plot of log10-transformed PFOS data (clear evidence of heavy tails, which supports our choice of the GAM scaled t family)

ctrl <- gam.control(nthreads = 4) # use multiple threads for fitting GAMs

### Hold out validation (80:20 method)
set.seed(166) 
PFOS_in_train<-createDataPartition(PFOS$PFOS, p = 4/5, list = FALSE)
PFOS_train_set<-PFOS[PFOS_in_train,]
PFOS_test_set<-rbind(PFOS[-PFOS_in_train,],validation_supp) # also added the eight samples initially omitted from modeling to reinforce the strenght of the test set


# About 24 minutes
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


appraise(full_PFOS_gam)
# Observed vs fitted values plot does not deviate from an approximate 1:1 trend line 
# QQ plot and histogram of residuals show few departures from normality (a few more points on the left tail than expected from normal data being the major note)
# Deviance residuals vs linear predictor plot shows little observable heteroscedasticity

# The same type of observations can be made about the pearson residuals vs linear predictor plot 
residualPlot(full_PFOS_gam)
# A few "outliers" shown in the plot, but the model (as would be expected of a scaled-t distribution) seems to produce output that is robust to any such points


plot(simulateResiduals(full_PFOS_gam),quantreg=T)
# Some deviation according to KS test, as well as a slight deviation in 
# quantiles, but not much that wouldn't be expected at large sample sizes
# Dispersion and Outlier Tests not significant

# Outliers (23 at the two margins (n = 1981), p-value = 0.07543)
testOutliers(full_PFOS_gam)


### Test for any significant patterns between predictor variable values and residuals
# get data
refit_PFOS <- full_PFOS_gam$model
# make residuals column
refit_PFOS$resid <- residuals(full_PFOS_gam,type="deviance")
# fit the same specified model to the residuals
PFOS_resid_fit <- gam(resid ~ Waterbody + Waterbody_Type
                      + s(Water_Level_sc,by = Waterbody,k=20)
                      + Trophic_Level
                      + Revised_Tissue
                      + s(Sampling.Year,by = Waterbody,k=20)
                      + s(Longitude,Latitude,k=525)
                      + s(Composite,bs="re"),
                      family=gaussian(), data=refit_PFOS, method="REML")
summary(PFOS_resid_fit)
# No real significant patterns detected: 1.98% of the deviance explained with only one significant difference from 0 in one factor level (Liver) of the tissue variable

ggplot(refit_PFOS, aes(x = Waterbody_Type, y = resid)) + 
  geom_point() + 
  geom_hline(yintercept = 0, color='red', linetype='dashed')


### Test for Spatial Autocorrelation: N.S. 
# p-value = 0.6336
PFOS_test_set$coords <- paste(PFOS_test_set$Longitude,", ",
                              PFOS_test_set$Latitude)
coords <- c(unique(PFOS_test_set$coords))
x_unique <- c(str_extract(coords, "^.+(?=,)"))
y_unique <- c(str_extract(coords, "(?<=, ).+$"))

res<-recalculateResiduals(simulateResiduals(full_PFOS_gam),
                          group = PFOS_test_set$coords)
testSpatialAutocorrelation(res, x = x_unique, y = y_unique)


### Test for Temporal Autocorrelation: N.S.
# p-value = 0.4839
res = recalculateResiduals(simulateResiduals(full_PFOS_gam), 
                           group = PFOS_test_set$Sampling.Year)
testTemporalAutocorrelation(res, time = unique(PFOS_test_set$Sampling.Year))


### The effectiveness of the model in predicting the test set
PFOS_results<-cbind(PFOS_test_set,as.data.frame(predict(full_PFOS_gam,
                                                        PFOS_test_set,
                                                        type="response",
                                                        se.fit=TRUE)))

PFOS_line<-lm(log(PFOS, base = 10) ~ (fit),
              data = PFOS_results)
summary(PFOS_line)

# Significance test (using two emmeans functions) that test the 1:1 slope assumption
# of predicted vs observed
# p = 0.4195
fit.emt <- emtrends(PFOS_line, ~1, var="fit")
fit.emt
test(fit.emt, null=1)

# OR 
# (Intercept) p = 0.793
# (fit) p = 0.419
PFOS_results$dif<-(log(PFOS_results$PFOS,base = 10)-PFOS_results$fit)
one_to_one<-lm(dif ~ fit,
               data = PFOS_results)
summary(one_to_one)
# Both validate the assumption of 1:1

# Plot of predicted (X) vs observed (y) on the log scale
ggplot(data=PFOS_results, aes(x=(fit), y=log(PFOS,base = 10),color=Waterbody)) +
  geom_point() +
  geom_pointrange(aes(xmin=fit-se.fit, xmax=fit+se.fit)) +
  geom_abline(slope=1, intercept=0) +
  theme_bw()


# This same process is repeated for the five other contaminants 

###### Full PFNA Model ########################################################
PFNA <- final_imputed_data[is.na(final_imputed_data$PFNA)==FALSE,]
car::qqPlot(log(PFNA$PFNA,base = 10))

ctrl <- gam.control(nthreads = 4) # use multiple threads for fitting

### Hold out validation (80:20 method)
set.seed(166)
PFNA_in_train<-createDataPartition(PFNA$PFNA, p = 4/5, list = FALSE)
PFNA_train_set<-PFNA[PFNA_in_train,]
PFNA_test_set<-rbind(PFNA[-PFNA_in_train,],validation_supp)

PFNA_test_set<-subset(PFNA_test_set,!is.na(PFNA))


# About 7 minutes
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

appraise(full_PFNA_gam)
# Histogram shows approximately normally distributed residuals, as does the QQ plot
# Deviance residuals plotted against a linear predictor show little evidence of heteroscedasticity
# Observed vs fitted values plot demonstrates an approximate 1:1 trend line

residualPlot(full_PFNA_gam)
# Besides a few outliers, the density of Pearson residuals is well distributed about
# the 0 horizontal line with little observable heteroscedasticity


plot(simulateResiduals(full_PFNA_gam),quantreg=T)
# No outliers or significant dispersion, minimal quantile deviations
# Significance for KS test (though again, is to be expected at higher sample sizes)

# Outliers (11 at the two margins (n = 1684), p-value = 0.6791)
testOutliers(full_PFNA_gam)


### Test for any significant patterns between predictor variable values and residuals
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
# No significant patterns detected: explains 1.35% of the deviance with no significant terms

ggplot(refit_PFNA, aes(x = Waterbody_Type, y = resid)) + 
  geom_point() + 
  geom_hline(yintercept = 0, color='red', linetype='dashed')


### Test for Spatial Autocorrelation: N.S. 
# p-value = 0.1704
# p-value = 0.01099 when the eight/seven extra data points are added, so there might
# be some evidence of residual spatial autocorrelation
# But we might expect such findings from data taken in an otherwise completely unsampled
# part of the Great Lakes watersheds, and any autocorrelation in these few point,
# which, as suggested by an observed Moran's I value of 0.06987, is fairly weak,
# doesn't appear to create temporal autocorrelation, nor does it impede the
# predictive power of the model for the test data (see below)
PFNA_test_set$coords <- paste(PFNA_test_set$Longitude,", ",
                              PFNA_test_set$Latitude)
coords <- c(unique(PFNA_test_set$coords))
x_unique <- c(str_extract(coords, "^.+(?=,)"))
y_unique <- c(str_extract(coords, "(?<=, ).+$"))

res<-recalculateResiduals(simulateResiduals(full_PFNA_gam),
                          group = PFNA_test_set$coords)
testSpatialAutocorrelation(res, x = x_unique, y = y_unique)


### Test for Temporal Autocorrelation: N.S.
# p-value = 0.2677
res = recalculateResiduals(simulateResiduals(full_PFNA_gam), 
                           group = PFNA_test_set$Sampling.Year)
testTemporalAutocorrelation(res, time = unique(PFNA_test_set$Sampling.Year))


### The effectiveness of the model in predicting the test set
PFNA_results<-cbind(PFNA_test_set,as.data.frame(predict(full_PFNA_gam,
                                                        PFNA_test_set,
                                                        type="response",
                                                        se.fit=TRUE)))
PFNA_line<-lm(log(PFNA, base = 10) ~ (fit),
              data = PFNA_results)
summary(PFNA_line)

# Significance test (using two emmeans functions) that test the 1:1 slope assumption
# of predicted vs observed
# p = 0.4088
fit.emt <- emtrends(PFNA_line, ~1, var="fit")
fit.emt
test(fit.emt, null=1)

# OR 
# (Intercept) p = 0.380
# (fit) p = 0.409
PFNA_results$dif<-(log(PFNA_results$PFNA,base = 10)-PFNA_results$fit)
one_to_one<-lm(dif ~ fit,
               data = PFNA_results)
summary(one_to_one)
# Both validate the assumption of 1:1

# Plot of predicted (X) vs observed (y) on the log scale
ggplot(data=PFNA_results, aes(x=(fit), y=log(PFNA,base = 10),color=Waterbody)) +
  geom_point() +
  geom_pointrange(aes(xmin=fit-se.fit, xmax=fit+se.fit)) +
  geom_abline(slope=1, intercept=0) +
  theme_bw()


###### Full PFDA Model ######################################
PFDA<-final_imputed_data[is.na(final_imputed_data$PFDA)==FALSE,]
car::qqPlot(log(PFDA$PFDA,base = 10))

ctrl <- gam.control(nthreads = 4) # use multiple threads for fitting

### Hold out validation (80:20 method)
set.seed(166)
PFDA_in_train<-createDataPartition(PFDA$PFDA, p = 4/5, list = FALSE)
PFDA_train_set<-PFDA[PFDA_in_train,]
PFDA_test_set<-rbind(PFDA[-PFDA_in_train,],validation_supp)

PFDA_test_set<-subset(PFDA_test_set,!is.na(PFDA))


# About 6 minutes
system.time(
  full_PFDA_gam <- gam(log10(PFDA) ~ Waterbody + Waterbody_Type
                       + s(Water_Level_sc,by = Waterbody,k=20)
                       + Trophic_Level
                       + Revised_Tissue
                       + s(Sampling.Year,by = Waterbody,k=20)
                       + s(Longitude,Latitude,k=250)
                       + s(Composite,bs="re"),
                       data = PFDA_train_set,
                       method = "REML",
                       family = scat(link = "identity"),
                       control = ctrl)
)

appraise(full_PFDA_gam)
# Histogram indicates approximate normality of residuals, QQ plot doesn't reveal any noticeable problems with the
# distribution of the residuals, and the plot of observed vs fitted values demonstrates a 1:1 fit
# The plot of deviance residuals vs linear predictor shows approximate homoscedasticity

residualPlot(full_PFDA_gam)
# No significant heteroscedacity: pink line doesn't diverge substantially from the
# theoretical 0 horizontal line (homoscedasticity)


plot(simulateResiduals(full_PFDA_gam),quantreg=T)
# No outliers, KS and Dispersion N.S., only slight lower quantile deviations

# Outliers (7 at the two margins (n = 1682), p-value = 0.09703)
testOutliers(full_PFDA_gam)


### Test for any significant patterns between predictor variable values and residuals
# get data
refit_PFDA <- full_PFDA_gam$model
# make residuals column
refit_PFDA$resid <- residuals(full_PFDA_gam,type="deviance")
# fit a model (same model)
PFDA_resid_fit <- gam(resid ~ Waterbody + Waterbody_Type
                      + s(Water_Level_sc,by = Waterbody,k=20)
                      + Trophic_Level
                      + Revised_Tissue
                      + s(Sampling.Year,by = Waterbody,k=20)
                      + s(Longitude,Latitude,k=250)
                      + s(Composite,bs="re"),
                      family=gaussian(), data=refit_PFDA, method="REML")
summary(PFDA_resid_fit)
# No significant patterns detected: explains 0.584% of the deviance with no significant terms

ggplot(refit_PFDA, aes(x = Waterbody_Type, y = resid)) + 
  geom_point() + 
  geom_hline(yintercept = 0, color='red', linetype='dashed')


### Test for Spatial Autocorrelation: N.S. 
# p-value = 0.7708
PFDA_test_set$coords <- paste(PFDA_test_set$Longitude,", ",
                              PFDA_test_set$Latitude)
coords <- c(unique(PFDA_test_set$coords))
x_unique <- c(str_extract(coords, "^.+(?=,)"))
y_unique <- c(str_extract(coords, "(?<=, ).+$"))

res<-recalculateResiduals(simulateResiduals(full_PFDA_gam),
                          group = PFDA_test_set$coords)
testSpatialAutocorrelation(res, x = x_unique, y = y_unique)

### Test for Temporal Autocorrelation: N.S.
# p-value = 0.4974
res = recalculateResiduals(simulateResiduals(full_PFDA_gam), 
                           group = PFDA_test_set$Sampling.Year)
testTemporalAutocorrelation(res, time = unique(PFDA_test_set$Sampling.Year))


### The effectiveness of the model in predicting the test set
PFDA_results<-cbind(PFDA_test_set,as.data.frame(predict(full_PFDA_gam,
                                                        PFDA_test_set,
                                                        type="response",
                                                        se.fit=TRUE)))

PFDA_line<-lm(log(PFDA, base = 10) ~ (fit),
              data = PFDA_results)
summary(PFDA_line)

# Significance test (using two emmeans functions) that test the 1:1 slope assumption
# of predicted vs observed
# p = 0.3266
fit.emt <- emtrends(PFDA_line, ~1, var="fit")
fit.emt
test(fit.emt, null=1)

# OR 
# (Intercept) p = 0.0968
# (fit) p = 0.3266
PFDA_results$dif<-(log(PFDA_results$PFDA,base = 10)-PFDA_results$fit)
one_to_one<-lm(dif ~ fit,
               data = PFDA_results)
summary(one_to_one)
# Both validate the assumption of 1:1

# Plot of predicted (X) vs observed (y) on the log scale
ggplot(data=PFDA_results, aes(x=(fit), y=log(PFDA,base = 10),color=Trophic_Level)) +
  geom_point() +
  geom_pointrange(aes(xmin=fit-se.fit, xmax=fit+se.fit)) +
  geom_abline(slope=1, intercept=0) +
  theme_bw()


###### Full PFUnA Model ######################################
PFUnA<-final_imputed_data[is.na(final_imputed_data$PFUnA)==FALSE,]
car::qqPlot(log(PFUnA$PFUnA,base = 10))

ctrl <- gam.control(nthreads = 4) # use multiple threads for fitting

### Hold out validation (80:20 method)
set.seed(130)
PFUnA_in_train<-createDataPartition(PFUnA$PFUnA, p = 4/5, list = FALSE)
PFUnA_train_set<-PFUnA[PFUnA_in_train,]
PFUnA_test_set<-rbind(PFUnA[-PFUnA_in_train,],validation_supp)

PFUnA_test_set<-subset(PFUnA_test_set,!is.na(PFUnA))


# About 8 minutes
system.time(
  full_PFUnA_gam <- gam(log10(PFUnA) ~ Waterbody + Waterbody_Type
                        + s(Water_Level_sc,by = Waterbody,k=20)
                        + Trophic_Level
                        + Revised_Tissue
                        + s(Sampling.Year,by = Waterbody,k=20)
                        + s(Longitude,Latitude,k=200)
                        + s(Composite,bs="re"),
                        data = PFUnA_train_set,
                        method = "REML",
                        family = scat(link = "identity"),
                        control = ctrl)
)

appraise(full_PFUnA_gam)
# Histogram shows approximate normality of residuals, QQ plot doesn't reveal any noticeable problems with the
# distribution of the residuals, and the plot of observed vs fitted values demonstrates a 1:1 fit
# The plot of deviance residuals vs linear predictor shows approximate homoscedasticity

residualPlot(full_PFUnA_gam)
# Same types of conclusions. No real heteroscedacity, pink line doesn't deviate
# substantially from the theoretical horizontal line


plot(simulateResiduals(full_PFUnA_gam),quantreg=T)
# No outliers, Dispersion test not significant, KS significant, slight quantile deviations

# Outliers (14 at the two margins (n = 1686), p-value = 0.7845)
testOutliers(full_PFUnA_gam)

### Test for any significant patterns between predictor variable values and residuals
# get data
refit_PFUnA <- full_PFUnA_gam$model
# make residuals column
refit_PFUnA$resid <- residuals(full_PFUnA_gam,type="deviance")
# fit a model (same model)
PFUnA_resid_fit <- gam(resid ~ Waterbody + Waterbody_Type
                       + s(Water_Level_sc,by = Waterbody,k=20)
                       + Trophic_Level
                       + Revised_Tissue
                       + s(Sampling.Year,by = Waterbody,k=20)
                       + s(Longitude,Latitude,k=200)
                       + s(Composite,bs="re"),
                       family=gaussian(), data=refit_PFUnA, method="REML")
summary(PFUnA_resid_fit)
# No significant patterns detected: explains 1.04% of the deviance with no significant terms (only one marginally significant (0.0593) factor level of the Revised_Tissue variable, Eggs)

ggplot(refit_PFUnA, aes(x = Waterbody, y = resid)) + 
  geom_point() + 
  geom_hline(yintercept = 0, color='red', linetype='dashed')


### Test for Spatial Autocorrelation: N.S. 
# p-value = 0.8511
PFUnA_test_set$coords <- paste(PFUnA_test_set$Longitude,", ",
                               PFUnA_test_set$Latitude)
coords <- c(unique(PFUnA_test_set$coords))
x_unique <- c(str_extract(coords, "^.+(?=,)"))
y_unique <- c(str_extract(coords, "(?<=, ).+$"))

res<-recalculateResiduals(simulateResiduals(full_PFUnA_gam),
                          group = PFUnA_test_set$coords)
testSpatialAutocorrelation(res, x = x_unique, y = y_unique)

### Test for Temporal Autocorrelation: N.S.
# p-value = 0.4416
res = recalculateResiduals(simulateResiduals(full_PFUnA_gam), 
                           group = PFUnA_test_set$Sampling.Year)
testTemporalAutocorrelation(res, time = unique(PFUnA_test_set$Sampling.Year))


### The effectiveness of the model in predicting the test set
PFUnA_results<-cbind(PFUnA_test_set,as.data.frame(predict(full_PFUnA_gam,
                                                          PFUnA_test_set,
                                                          type="response",
                                                          se.fit=TRUE)))

PFUnA_line<-lm(log(PFUnA, base = 10) ~ (fit),
               data = PFUnA_results)
summary(PFUnA_line)

# Significance test (using two emmeans functions) that test the 1:1 slope assumption
# of predicted vs observed
# p = 0.9431
fit.emt <- emtrends(PFUnA_line, ~1, var="fit")
fit.emt
test(fit.emt, null=1)

# OR 
# (Intercept) p = 0.519
# (fit) p = 0.943
PFUnA_results$dif<-(log(PFUnA_results$PFUnA,base = 10)-PFUnA_results$fit)
one_to_one<-lm(dif ~ fit,
               data = PFUnA_results)
summary(one_to_one)
# Both validate the assumption of 1:1

# Plot of predicted (X) vs observed (y) on the log scale
ggplot(data=PFUnA_results, aes(x=(fit), y=log(PFUnA,base = 10),color=Waterbody)) +
  geom_point() +
  geom_pointrange(aes(xmin=fit-se.fit, xmax=fit+se.fit)) +
  geom_abline(slope=1, intercept=0) +
  theme_bw()


###### Full PFDoA Model #####################################
PFDoA<-final_imputed_data[is.na(final_imputed_data$PFDoA)==FALSE,]
car::qqPlot(log(PFDoA$PFDoA,base = 10))

ctrl <- gam.control(nthreads = 4) # use multiple threads for fitting

### Hold out validation (80:20 method)
set.seed(145)
PFDoA_in_train <- createDataPartition(PFDoA$PFDoA, p = 4/5, list = FALSE)
PFDoA_train_set <- PFDoA[PFDoA_in_train,]
PFDoA_test_set <- rbind(PFDoA[-PFDoA_in_train,],validation_supp)

PFDoA_test_set<-subset(PFDoA_test_set,!is.na(PFDoA))


# About 11 minutes
system.time(
  full_PFDoA_gam <- gam(log10(PFDoA) ~ Waterbody + Waterbody_Type
                        + s(Water_Level_sc,by = Waterbody)
                        + Trophic_Level
                        + Revised_Tissue
                        + s(Sampling.Year,by = Waterbody,k=20)
                        + s(Longitude,Latitude,k=300)
                        + s(Composite,bs="re"),
                        data = PFDoA_train_set,
                        method = "REML",
                        family = scat(link = "identity"),
                        control = ctrl)
)

appraise(full_PFDoA_gam)
# Histogram shows approximately normally distributed residuals, as does the QQ plot
# Deviance residuals plotted against a linear predictor show little evidence of heteroscedasticity
# Observed vs fitted values plot demonstrates an approximate 1:1 trend line

residualPlot(full_PFDoA_gam)
# No heteroscedacity that's readily observable: pink line doesn't deviate
# substantially from the theoretical horizontal line (except marginally at the
# ends)


plot(simulateResiduals(full_PFDoA_gam),quantreg=T)
# No outliers, KS and Dispersion tests not significant, some quantile deviations

# Outliers (7 at the two margins (n = 1668), p-value = 0.09616)
testOutliers(full_PFDoA_gam)


### Test for any significant patterns between predictor variable values and residuals
# get data
refit_PFDoA <- full_PFDoA_gam$model
# make residuals column
refit_PFDoA$resid <- residuals(full_PFDoA_gam,type="deviance")
# fit a model (same model)
PFDoA_resid_fit <- gam(resid ~ Waterbody + Waterbody_Type
                       + s(Water_Level_sc,by = Waterbody)
                       + Trophic_Level
                       + Revised_Tissue
                       + s(Sampling.Year,by = Waterbody,k=20)
                       + s(Longitude,Latitude,k=300)
                       + s(Composite,bs="re"),
                       family=gaussian(), data=refit_PFDoA, method="REML")
summary(PFDoA_resid_fit)
# No significant patterns detected: explains 1.09% of the deviance, and, besides significance in a few of the factor levels of the Revised Tissue variable, no systematically significant relationships between predictor variables and the model residuals

ggplot(refit_PFDoA, aes(x = Revised_Tissue, y = resid)) + 
  geom_point() + 
  geom_hline(yintercept = 0, color='red', linetype='dashed')


### Test for Spatial Autocorrelation: N.S. 
# p-value = 0.8459
PFDoA_test_set$coords <- paste(PFDoA_test_set$Longitude,", ",
                               PFDoA_test_set$Latitude)
coords <- c(unique(PFDoA_test_set$coords))
x_unique <- c(str_extract(coords, "^.+(?=,)"))
y_unique <- c(str_extract(coords, "(?<=, ).+$"))

res<-recalculateResiduals(simulateResiduals(full_PFDoA_gam),
                          group = PFDoA_test_set$coords)
testSpatialAutocorrelation(res, x = x_unique, y = y_unique)

### Test for Temporal Autocorrelation: N.S.
# p-value = 0.8813
res = recalculateResiduals(simulateResiduals(full_PFDoA_gam), 
                           group = PFDoA_test_set$Sampling.Year)
testTemporalAutocorrelation(res, time = unique(PFDoA_test_set$Sampling.Year))


### The effectiveness of the model in predicting the test set
PFDoA_results<-cbind(PFDoA_test_set,as.data.frame(predict(full_PFDoA_gam,
                                                          PFDoA_test_set,
                                                          type="response",
                                                          se.fit=TRUE)))

PFDoA_line<-lm(log(PFDoA, base = 10) ~ (fit),
               data = PFDoA_results)
summary(PFDoA_line)

# Significance test (using two emmeans functions) that test the 1:1 slope assumption
# of predicted vs observed
# p = 0.0639
fit.emt <- emtrends(PFDoA_line, ~1, var="fit")
fit.emt
test(fit.emt, null=1)

# OR 
# (Intercept) p = 0.1463
# (fit) p = 0.0639
PFDoA_results$dif<-(log(PFDoA_results$PFDoA,base = 10)-PFDoA_results$fit)
one_to_one<-lm(dif ~ fit,
               data = PFDoA_results)
summary(one_to_one)
# Both validate the assumption of 1:1

# Plot of predicted (X) vs observed (y) on the log scale
ggplot(data=PFDoA_results, aes(x=(fit), y=log(PFDoA,base = 10),color=Waterbody)) +
  geom_point() +
  geom_pointrange(aes(xmin=fit-se.fit, xmax=fit+se.fit)) +
  geom_abline(slope=1, intercept=0)


###### Full PFTrDA Model ######################################
PFTrDA<-final_imputed_data[is.na(final_imputed_data$PFTrDA)==FALSE,]
car::qqPlot(log(PFTrDA$PFTrDA,base = 10))

ctrl <- gam.control(nthreads = 4) # use multiple threads for fitting

### Hold out validation (80:20 method)
set.seed(145)
PFTrDA_in_train<-createDataPartition(PFTrDA$PFTrDA, p = 4/5, list = FALSE)
PFTrDA_train_set<-PFTrDA[PFTrDA_in_train,]
PFTrDA_test_set<-rbind(PFTrDA[-PFTrDA_in_train,],validation_supp)

PFTrDA_test_set<-subset(PFTrDA_test_set,!is.na(PFTrDA))


# About 11 minutes
system.time(
  full_PFTrDA_gam <- gam(log10(PFTrDA) ~ Waterbody + Waterbody_Type
                         + s(Water_Level_sc,by = Waterbody,k=20)
                         + Trophic_Level
                         + Revised_Tissue
                         + s(Sampling.Year, by = Waterbody)
                         + s(Longitude,Latitude,k=300)
                         + s(Composite,bs="re"),
                         data = PFTrDA_train_set,
                         method = "REML",
                         family = scat(link = "identity"),
                         control = ctrl)
)

appraise(full_PFTrDA_gam)
# Fitted vs observed values plot demonstrates an approximate 1:1 trend line
# Some deviations from normality in the residual distribution along the left tail,
# and some heteroscedasticity in deviance residuals
# But otherwise, residual distributions look fairly normal

residualPlot(full_PFTrDA_gam)
# Pearson residuals support this observation: no significant heteroscedacity that
# would shift the balance of residuals off of the center line (homoscedasticity) 


plot(simulateResiduals(full_PFTrDA_gam),quantreg=T)
# Significant Outlier and KS tests, some quantile deviations
# The other diagnostics (see below) show a good fit to the data that can accurately predict test values

# Outliers (18 at the two margins (n = 1263), p-value = 0.0244)
testOutliers(full_PFTrDA_gam) # Not too many outliers (would expect some in environmental contaminant data, especially since our literature search includes some publications that recorded concentrations after major spill events)
# And we intentionally picked a data distribution family (scaled-t) that is robust to outliers and can still make accurate predictions of test data (see below)


### Test for any significant patterns between predictor variable values and residuals
# get data
refit_PFTrDA <- full_PFTrDA_gam$model
# make residuals column
refit_PFTrDA$resid <- residuals(full_PFTrDA_gam,type="deviance")
# fit a model (same model)
PFTrDA_resid_fit <- gam(resid ~ Waterbody + Waterbody_Type
                        + s(Water_Level_sc,by = Waterbody,k=20)
                        + Trophic_Level
                        + Revised_Tissue
                        + s(Sampling.Year,by = Waterbody)
                        + s(Longitude,Latitude,k=300)
                        + s(Composite,bs="re"),
                        family=gaussian(), data=refit_PFTrDA, method="REML")
summary(PFTrDA_resid_fit)
# No significant patterns detected: explains 2.28% of the deviance with no significant terms (only a marginally significant (0.0517) factor level (Eggs) of the Revised_Tissue variable)

ggplot(refit_PFTrDA, aes(x = Waterbody, y = resid)) + 
  geom_point() + 
  geom_hline(yintercept = 0, color='red', linetype='dashed')


### Test for Spatial Autocorrelation: N.S. 
# p-value = 0.876
PFTrDA_test_set$coords <- paste(PFTrDA_test_set$Longitude,", ",
                                PFTrDA_test_set$Latitude)
coords <- c(unique(PFTrDA_test_set$coords))
x_unique <- c(str_extract(coords, "^.+(?=,)"))
y_unique <- c(str_extract(coords, "(?<=, ).+$"))

res<-recalculateResiduals(simulateResiduals(full_PFTrDA_gam),
                          group = PFTrDA_test_set$coords)
testSpatialAutocorrelation(res, x = x_unique, y = y_unique)

### Test for Temporal Autocorrelation: N.S.
# p-value = 0.3146
res = recalculateResiduals(simulateResiduals(full_PFTrDA_gam), 
                           group = PFTrDA_test_set$Sampling.Year)
testTemporalAutocorrelation(res, time = unique(PFTrDA_test_set$Sampling.Year))


### The effectiveness of the model in predicting the test set
PFTrDA_results<-cbind(PFTrDA_test_set,as.data.frame(predict(full_PFTrDA_gam,
                                                            PFTrDA_test_set,
                                                            type="response",
                                                            se.fit=TRUE)))

PFTrDA_line<-lm(log(PFTrDA, base = 10) ~ (fit),
                data = PFTrDA_results)
summary(PFTrDA_line)

# Significance test (using two emmeans functions) that test the 1:1 slope assumption
# of predicted vs observed
# p = 0.0919
fit.emt <- emtrends(PFTrDA_line, ~1, var="fit")
fit.emt
test(fit.emt, null=1)

# OR 
# (Intercept) p = 0.1885
# (fit) p = 0.0919
PFTrDA_results$dif<-(log(PFTrDA_results$PFTrDA,base = 10)-PFTrDA_results$fit)
one_to_one<-lm(dif ~ fit,
               data = PFTrDA_results)
summary(one_to_one)
# Both validate the assumption of 1:1

# Plot of predicted (X) vs observed (y) on the log scale
ggplot(data=PFTrDA_results, aes(x=(fit), y=log(PFTrDA,base = 10),color=Waterbody)) +
  geom_point() +
  geom_pointrange(aes(xmin=fit-se.fit, xmax=fit+se.fit)) +
  geom_abline(slope=1, intercept=0) +
  theme_bw()

###### Save the model for future work ##########################
save(full_PFOS_gam,file = "full_PFOS_gam.Rdata")
save(full_PFNA_gam,file = "full_PFNA_gam.Rdata")
save(full_PFDA_gam,file = "full_PFDA_gam.Rdata")
save(full_PFUnA_gam,file = "full_PFUnA_gam.Rdata")
save(full_PFDoA_gam,file = "full_PFDoA_gam.Rdata")
save(full_PFTrDA_gam,file = "full_PFTrDA_gam.Rdata")





