# Step_5_Modeling.R
# Contains all code used to construct models for the six selected analytes (PFOS, PFNA,
# PFDA, PFUnA, PFDoA, PFTrDA)

# Written by Peter O. Martin (https://orcid.org/0009-0009-9070-9200)

# Working directory
setwd("~/Desktop/Publications/Martin et al., 2025")

## Packages
# Data formatting and combining
library(tidyverse)
library(lattice)
library(caret)

# Modeling
library(nlme)
library(mgcv)
library(mgcViz)
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
final_imputed_data <- read.csv("Finalized_Imputed_Data_Frame.csv",header = TRUE)
validation_supp <- read.csv("Finalized_Supp_Validation_Data_Frame.csv",header = TRUE)

###### Pre-Modeling Formatting ################################################
# Converting columns in each data frame to proper format (i.e., ordered factors)
# Waterbody variable (5 levels)
final_imputed_data$Waterbody<-factor(final_imputed_data$Waterbody,
                                     levels = c("Lake Superior",
                                                "Lake Michigan",
                                                "Lake Huron",
                                                "Lake Erie",
                                                "Lake Ontario"))
validation_supp$Waterbody<-factor(validation_supp$Waterbody,
                                  levels = c("Lake Superior",
                                             "Lake Michigan",
                                             "Lake Huron",
                                             "Lake Erie",
                                             "Lake Ontario"))

# Waterbody_Type variable (3 levels)
final_imputed_data$Waterbody_Type<-factor(final_imputed_data$Waterbody_Type,
                                          levels = c("Inland waters",
                                                     "Connecting channel",
                                                     "Lake"))
validation_supp$Waterbody_Type<-factor(validation_supp$Waterbody_Type,
                                       levels = c("Inland waters",
                                                  "Connecting channel",
                                                  "Lake"))

# Composite variable (2 levels)
final_imputed_data$Composite<-factor(final_imputed_data$Composite)
validation_supp$Composite<-factor(validation_supp$Composite)

# Trophic_Level variable (7 levels)
final_imputed_data$Trophic_Level<-factor(final_imputed_data$Trophic_Level,
                                         levels = c("Primary Producer",
                                                    "Primary Consumer",
                                                    "Secondary Consumer",
                                                    "Tertiary Consumer",
                                                    "Quaternary Consumer",
                                                    "Piscivorous/Insectivorous Bird",
                                                    "Apex Predator"))
validation_supp$Trophic_Level<-factor(validation_supp$Trophic_Level,
                                      levels = c("Primary Producer",
                                                 "Primary Consumer",
                                                 "Secondary Consumer",
                                                 "Tertiary Consumer",
                                                 "Quaternary Consumer",
                                                 "Piscivorous/Insectivorous Bird",
                                                 "Apex Predator"))

# Class variable (14 levels)
final_imputed_data$Class<-factor(final_imputed_data$Class,
                                 levels = c("Algae","Plantae (Magnoliopsida)",
                                            "Annelida","Bivalvia","Gastropoda",
                                            "Zooplankton",
                                            "Shrimp, water fleas, and allies",
                                            "Insecta","Astacoidea","Amphibia",
                                            "Pisces","Reptilia","Aves","Mammalia"))
validation_supp$Class<-factor(validation_supp$Class,
                              levels = c("Algae","Plantae (Magnoliopsida)",
                                         "Annelida","Bivalvia","Gastropoda",
                                         "Zooplankton",
                                         "Shrimp, water fleas, and allies",
                                         "Insecta","Astacoidea","Amphibia",
                                         "Pisces","Reptilia","Aves","Mammalia"))

# Revised_Tissue variable (6 levels)
final_imputed_data$Revised_Tissue<-factor(final_imputed_data$Revised_Tissue,
                                          levels = c("Muscle",
                                                     "Whole organism homogenate",
                                                     "Misc. Tissue",
                                                     "Liver","Blood","Eggs"))
validation_supp$Revised_Tissue<-factor(validation_supp$Revised_Tissue,
                                       levels = c("Muscle",
                                                  "Whole organism homogenate",
                                                  "Misc. Tissue",
                                                  "Liver","Blood","Eggs"))

str(final_imputed_data)
str(validation_supp)

###### Full Data Frame Modeling ###############################################
# We decided to focus analysis on the following contaminants that all had ≥ 80% 
# detection frequency:
# PFOS, PFNA, PFDA, PFUnA, PFDoA, and PFTrDA

# Iterative modeling led us to choose the following variables for these six GAMs

# Fixed Effects: Waterbody, Waterbody_Type, Trophic_Level, Revised_Tissue, 
# Water_Level_sc, and Sampling_Year
# Random Effects: Composite

###### HERE can access the saved models to explore diagnostics ################
load("Models/full_PFOS_gam.Rdata") 
load("Models/full_PFNA_gam.Rdata")
load("Models/full_PFDA_gam.Rdata")
load("Models/full_PFUnA_gam.Rdata")
load("Models/full_PFDoA_gam.Rdata")
load("Models/full_PFTrDA_gam.Rdata")

###### Full PFOS model ########################################################
# Extract the subset of samples that quantified PFOS concentrations (i.e., remove NAs)
PFOS<-subset(final_imputed_data,!is.na(PFOS))

# constructs a QQ plot of log10-transformed PFOS data (clear evidence of heavy tails,
# which supports our choice of the GAM scaled t family)
car::qqPlot(log(PFOS$PFOS,base = 10))

ctrl <- gam.control(nthreads = 4) # use multiple threads for fitting GAMs

### Hold out validation (80:20 method)
set.seed(166) 
# Create a random 80:20 partition of the PFOS data frame into training and test data 
# sets, using the createDataPartition() function from the caret package
PFOS_in_train<-createDataPartition(PFOS$PFOS, p = 4/5, list = FALSE)
PFOS_train_set<-PFOS[PFOS_in_train,]
PFOS_test_set<-rbind(PFOS[-PFOS_in_train,],validation_supp) # also added the eight
# samples initially omitted from modeling to reinforce the strength of the test set


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
# QQ plot and histogram of residuals show few departures from normality (a few more
# points on the left tail than expected from normal data being the major note)
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
# Extract data from model
refit_PFOS <- full_PFOS_gam$model
# Make a residuals column
refit_PFOS$resid <- residuals(full_PFOS_gam,type="deviance")
# Fit the same specified model to the residuals
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

# We can use this piece of code to examine the distribution of residuals within one 
# predictor variable across its different levels
ggplot(refit_PFOS, aes(x = Revised_Tissue, y = resid)) + 
  geom_point() + 
  geom_hline(yintercept = 0, color='red', linetype='dashed')


### Test for Spatial Autocorrelation: N.S. 
# p-value = 0.6336
# Using the testSpatialAutocorrelation() function from the DHARMa package
# We first create a new column in the test_set with the lat/long coordinates
# bound together
PFOS_test_set$coords <- paste(PFOS_test_set$Longitude,", ",
                              PFOS_test_set$Latitude)
# Then we extract the unique combinations of those coordinates, and separate the
# combinations out into the unique Longitudes (x_unique) and Latitudes (y_unique)
coords <- c(unique(PFOS_test_set$coords))
x_unique <- c(str_extract(coords, "^.+(?=,)"))
y_unique <- c(str_extract(coords, "(?<=, ).+$"))

# We recalculate the residuals (using DHARMa's built-in function created exactly for
# this purpose) over the set of unique coordinates (i.e., with each unique coordinate
# pair treated as a group by the function)
res<-recalculateResiduals(simulateResiduals(full_PFOS_gam),
                          group = PFOS_test_set$coords)
# Finally, we test for spatial autocorrelation, feeding into the function the
# recalculated residuals and the set of unique Longitudes and Latitudes
testSpatialAutocorrelation(res, x = x_unique, y = y_unique)


### Test for Temporal Autocorrelation: N.S.
# p-value = 0.4839
# Using the testTemporalAutocorrelation() function from the DHARMa package
# Again, we recalculate model residuals, this time over the array of sampling years
res <- recalculateResiduals(simulateResiduals(full_PFOS_gam), 
                           group = PFOS_test_set$Sampling.Year)
# We then test for temporal autocorrelation, feeding into the function the recalculated
# residuals as well as the unique set of sampling years included in the test set
testTemporalAutocorrelation(res, time = unique(PFOS_test_set$Sampling.Year))


### The effectiveness of the model in predicting the test set
# We first create two new columns in the test_set data frame: one column, fit, contains
# the model's predicted concentration values from the test set, and the other column,
# se.fit, contains the standard error value associated with each model prediction
PFOS_results<-cbind(PFOS_test_set,as.data.frame(predict(full_PFOS_gam,
                                                        PFOS_test_set,
                                                        type="response",
                                                        se.fit=TRUE)))

# We then run a linear model to test whether we see a 1:1 observed:predicted value
# relationship. We would expect that, if the model has good predictive accuracy,
# that a linear regression of observed vs predicted concentration values would yield
# a slope estimate of 1 and an intercept estimate of 0. In other words, there would be
# no significant differences between observed PFOS concentrations and the values that  
# the model predicts from the set of predictor variables and predictor-response 
# relationships that we have specified.
PFOS_line<-lm(log(PFOS, base = 10) ~ (fit),
              data = PFOS_results)
summary(PFOS_line)
# We find, from this regression, an intercept value (-0.01553) that is not significantly
# different from 0 (p = 0.793) and a slope value that is highly significant (p < 2e-16)
# and approximately 1 (0.97108). These findings provide evidence for the predictive
# accuracy of the model.


# We can also calculate residuals (dif), and perform a linear regression of residuals 
# against predicted values. For this linear regression, we would expect there to be no
# significant estimates, either for the intercept or for the slope. If the line of
# observed vs predicted values is truly 1:1, then we would expect non-significant
# deviations from a horizontal line through 0.
PFOS_results$dif<-(log(PFOS_results$PFOS,base = 10)-PFOS_results$fit)
one_to_one<-lm(dif ~ fit,
               data = PFOS_results)
summary(one_to_one)
# Intercept = -0.01553; p = 0.793
# fit = -0.02892; p = 0.419
# This test also validates the assumption of 1:1

# We can perform the exact same residuals test using two emmeans functions:
# p = 0.4195 (notice that the p value here is identical to that calculated for the slope
# estimate of the other residual test above)
fit.emt <- emtrends(PFOS_line, ~1, var="fit")
fit.emt
test(fit.emt, null=1)

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
refit_PFNA <- full_PFNA_gam$model
refit_PFNA$resid <- residuals(full_PFNA_gam,type="deviance")

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

ggplot(refit_PFNA, aes(x = Trophic_Level, y = resid)) + 
  geom_point() + 
  geom_hline(yintercept = 0, color='red', linetype='dashed')


### Test for Spatial Autocorrelation: N.S. 
# p-value = 0.1704
# p-value = 0.01099 when the eight/seven extra data points are added, so there might
# be some evidence of residual spatial autocorrelation
# But we might expect such findings from data that was sampled from an otherwise completely unsampled
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
# Intercept = -0.02075; p = 0.38
# fit = 0.96793; p < 2e-16

# OR 
PFNA_results$dif<-(log(PFNA_results$PFNA,base = 10)-PFNA_results$fit)
one_to_one<-lm(dif ~ fit,
               data = PFNA_results)
summary(one_to_one)
# Intercept = -0.02075; p = 0.380
# fit = -0.03207; p = 0.409
# Both validate the assumption of 1:1

# Significant residuals test (using two emmeans functions)
# p = 0.4088
fit.emt <- emtrends(PFNA_line, ~1, var="fit")
fit.emt
test(fit.emt, null=1)

# Plot of predicted (X) vs observed (y) on the log scale
ggplot(data=PFNA_results, aes(x=(fit), y=log(PFNA,base = 10),color=Waterbody)) +
  geom_point() +
  geom_pointrange(aes(xmin=fit-se.fit, xmax=fit+se.fit)) +
  geom_abline(slope=1, intercept=0) +
  theme_bw()


###### Full PFDA Model ########################################################
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
# Histogram indicates approximate normality of residuals, QQ plot doesn't reveal any 
# noticeable problems with the distribution of the residuals, and the plot of observed 
# vs fitted values demonstrates a 1:1 fit. The plot of deviance residuals vs linear
# predictor shows approximate homoscedasticity

residualPlot(full_PFDA_gam)
# No significant heteroscedacity: pink line doesn't diverge substantially from the
# theoretical 0 horizontal line (homoscedasticity)


plot(simulateResiduals(full_PFDA_gam),quantreg=T)
# No outliers, KS and Dispersion N.S., only slight lower quantile deviations

# Outliers (7 at the two margins (n = 1682), p-value = 0.09703)
testOutliers(full_PFDA_gam)


### Test for any significant patterns between predictor variable values and residuals
refit_PFDA <- full_PFDA_gam$model
refit_PFDA$resid <- residuals(full_PFDA_gam,type="deviance")

PFDA_resid_fit <- gam(resid ~ Waterbody + Waterbody_Type
                      + s(Water_Level_sc,by = Waterbody,k=20)
                      + Trophic_Level
                      + Revised_Tissue
                      + s(Sampling.Year,by = Waterbody,k=20)
                      + s(Longitude,Latitude,k=250)
                      + s(Composite,bs="re"),
                      family=gaussian(), data=refit_PFDA, method="REML")
summary(PFDA_resid_fit)
# No significant patterns detected: explains 0.584% of the deviance with no 
# significant terms

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
# Intercept = -0.03329; p = 0.0968
# fit = 1.03789; p < 2e-16

# OR 
PFDA_results$dif<-(log(PFDA_results$PFDA,base = 10)-PFDA_results$fit)
one_to_one<-lm(dif ~ fit,
               data = PFDA_results)
summary(one_to_one)
# Intercept = -0.03329; p = 0.0968
# fit = 0.03789; p = 0.3266
# Both validate the assumption of 1:1

# Significant residuals test (using two emmeans functions)
# p = 0.3266
fit.emt <- emtrends(PFDA_line, ~1, var="fit")
fit.emt
test(fit.emt, null=1)

# Plot of predicted (X) vs observed (y) on the log scale
ggplot(data=PFDA_results, aes(x=(fit), y=log(PFDA,base = 10),color=Trophic_Level)) +
  geom_point() +
  geom_pointrange(aes(xmin=fit-se.fit, xmax=fit+se.fit)) +
  geom_abline(slope=1, intercept=0) +
  theme_bw()


###### Full PFUnA Model #######################################################
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
# Histogram shows approximate normality of residuals, QQ plot doesn't reveal any 
# noticeable problems with the distribution of the residuals, and the plot of observed 
# vs fitted values demonstrates a 1:1 fit
# The plot of deviance residuals vs linear predictor shows approximate homoscedasticity

residualPlot(full_PFUnA_gam)
# Same types of conclusions. No real heteroscedacity, pink line doesn't deviate
# substantially from the theoretical horizontal line


plot(simulateResiduals(full_PFUnA_gam),quantreg=T)
# No outliers, Dispersion test not significant, KS significant, slight quantile deviations

# Outliers (14 at the two margins (n = 1686), p-value = 0.7845)
testOutliers(full_PFUnA_gam)

### Test for any significant patterns between predictor variable values and residuals
refit_PFUnA <- full_PFUnA_gam$model
refit_PFUnA$resid <- residuals(full_PFUnA_gam,type="deviance")

PFUnA_resid_fit <- gam(resid ~ Waterbody + Waterbody_Type
                       + s(Water_Level_sc,by = Waterbody,k=20)
                       + Trophic_Level
                       + Revised_Tissue
                       + s(Sampling.Year,by = Waterbody,k=20)
                       + s(Longitude,Latitude,k=200)
                       + s(Composite,bs="re"),
                       family=gaussian(), data=refit_PFUnA, method="REML")
summary(PFUnA_resid_fit)
# No significant patterns detected: explains 1.04% of the deviance with no significant
# terms (only one marginally significant (0.0593) factor level of the Revised_Tissue 
# variable, Eggs)

ggplot(refit_PFUnA, aes(x = Revised_Tissue, y = resid)) + 
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
# Intercept = 0.01290; p = 0.519
# fit = 1.00212; p < 2e-16

# OR 
PFUnA_results$dif<-(log(PFUnA_results$PFUnA,base = 10)-PFUnA_results$fit)
one_to_one<-lm(dif ~ fit,
               data = PFUnA_results)
summary(one_to_one)
# Intercept = 0.012900; p = 0.519
# fit = 0.002123; p = 0.943
# Both validate the assumption of 1:1

# Significant residuals test (using two emmeans functions)
# p = 0.9431
fit.emt <- emtrends(PFUnA_line, ~1, var="fit")
fit.emt
test(fit.emt, null=1)

# Plot of predicted (X) vs observed (y) on the log scale
ggplot(data=PFUnA_results, aes(x=(fit), y=log(PFUnA,base = 10),color=Waterbody)) +
  geom_point() +
  geom_pointrange(aes(xmin=fit-se.fit, xmax=fit+se.fit)) +
  geom_abline(slope=1, intercept=0) +
  theme_bw()


###### Full PFDoA Model #######################################################
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
refit_PFDoA <- full_PFDoA_gam$model
refit_PFDoA$resid <- residuals(full_PFDoA_gam,type="deviance")

PFDoA_resid_fit <- gam(resid ~ Waterbody + Waterbody_Type
                       + s(Water_Level_sc,by = Waterbody)
                       + Trophic_Level
                       + Revised_Tissue
                       + s(Sampling.Year,by = Waterbody,k=20)
                       + s(Longitude,Latitude,k=300)
                       + s(Composite,bs="re"),
                       family=gaussian(), data=refit_PFDoA, method="REML")
summary(PFDoA_resid_fit)
# No significant patterns detected: explains 1.09% of the deviance, and, besides
# significance in a few of the factor levels of the Revised Tissue variable, no
# systematically significant relationships between predictor variables and the model 
# residuals

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
# Intercept = -0.03050; p = 0.146
# fit = 0.92938; p < 2e-16

# OR 
PFDoA_results$dif<-(log(PFDoA_results$PFDoA,base = 10)-PFDoA_results$fit)
one_to_one<-lm(dif ~ fit,
               data = PFDoA_results)
summary(one_to_one)
# Intercept = -0.03050; p = 0.1463
# fit = -0.07062; p = 0.0639
# Both validate the assumption of 1:1

# Significant residuals test (using two emmeans functions)
# p = 0.0639
fit.emt <- emtrends(PFDoA_line, ~1, var="fit")
fit.emt
test(fit.emt, null=1)

# Plot of predicted (X) vs observed (y) on the log scale
ggplot(data=PFDoA_results, aes(x=(fit), y=log(PFDoA,base = 10),color=Waterbody)) +
  geom_point() +
  geom_pointrange(aes(xmin=fit-se.fit, xmax=fit+se.fit)) +
  geom_abline(slope=1, intercept=0) +
  theme_bw()


###### Full PFTrDA Model ######################################################
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
testOutliers(full_PFTrDA_gam) 
# Not too many outliers (would expect some in environmental contaminant data, especially
# since our literature search includes some publications that recorded concentrations
# after major spill events)
# And we intentionally picked a data distribution family (scaled-t) that is robust to
# outliers and can still make accurate predictions of test data (see below)


### Test for any significant patterns between predictor variable values and residuals
refit_PFTrDA <- full_PFTrDA_gam$model
refit_PFTrDA$resid <- residuals(full_PFTrDA_gam,type="deviance")

PFTrDA_resid_fit <- gam(resid ~ Waterbody + Waterbody_Type
                        + s(Water_Level_sc,by = Waterbody,k=20)
                        + Trophic_Level
                        + Revised_Tissue
                        + s(Sampling.Year,by = Waterbody)
                        + s(Longitude,Latitude,k=300)
                        + s(Composite,bs="re"),
                        family=gaussian(), data=refit_PFTrDA, method="REML")
summary(PFTrDA_resid_fit)
# No significant patterns detected: explains 2.28% of the deviance with no significant
# terms (only a marginally significant (0.0517) factor level (Eggs) of the Revised_Tissue
# variable)

ggplot(refit_PFTrDA, aes(x = Revised_Tissue, y = resid)) + 
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
# Intercept = -0.03506; p = 0.188
# fit = 0.93669; p < 2e-16

# OR 
PFTrDA_results$dif<-(log(PFTrDA_results$PFTrDA,base = 10)-PFTrDA_results$fit)
one_to_one<-lm(dif ~ fit,
               data = PFTrDA_results)
summary(one_to_one)
# Intercept = -0.03506; p = 0.1885
# fit = -0.06331; p = 0.0919
# Both validate the assumption of 1:1

# Significant residuals test (using two emmeans functions)
# p = 0.0919
fit.emt <- emtrends(PFTrDA_line, ~1, var="fit")
fit.emt
test(fit.emt, null=1)

# Plot of predicted (X) vs observed (y) on the log scale
ggplot(data=PFTrDA_results, aes(x=(fit), y=log(PFTrDA,base = 10),color=Waterbody)) +
  geom_point() +
  geom_pointrange(aes(xmin=fit-se.fit, xmax=fit+se.fit)) +
  geom_abline(slope=1, intercept=0) +
  theme_bw()

###### Save the model for table and figure generation #########################
save(full_PFOS_gam,file = "Models/full_PFOS_gam.Rdata")
save(full_PFNA_gam,file = "Models/full_PFNA_gam.Rdata")
save(full_PFDA_gam,file = "Models/full_PFDA_gam.Rdata")
save(full_PFUnA_gam,file = "Models/full_PFUnA_gam.Rdata")
save(full_PFDoA_gam,file = "Models/full_PFDoA_gam.Rdata")
save(full_PFTrDA_gam,file = "Models/full_PFTrDA_gam.Rdata")

# Delete excess variables
rm(coords,x_unique,y_unique,res,
   refit_PFOS,refit_PFNA,refit_PFDA,refit_PFUnA,refit_PFDoA,refit_PFTrDA,
   PFOS_resid_fit,PFNA_resid_fit,PFDA_resid_fit,PFUnA_resid_fit,
   PFDoA_resid_fit,PFTrDA_resid_fit,
   PFOS_test_set,PFNA_test_set,PFDA_test_set,PFUnA_test_set,PFDoA_test_set,PFTrDA_test_set,
   PFOS_train_set,PFNA_train_set,PFDA_train_set,PFUnA_train_set,
   PFDoA_train_set,PFTrDA_train_set,
   PFOS_results,PFNA_results,PFDA_results,PFUnA_results,PFDoA_results,PFTrDA_results,
   PFOS_line,PFNA_line,PFDA_line,PFUnA_line,PFDoA_line,PFTrDA_line,
   PFOS_in_train,PFNA_in_train,PFDA_in_train,PFUnA_in_train,PFDoA_in_train,PFTrDA_in_train,
   fit.emt,one_to_one,ctrl)


