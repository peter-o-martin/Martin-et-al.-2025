# Step_4_Supplementary_Validation_Dataset.R
# Contains all code used to format the supplementary validation data frame so that
# it can be used in the modeling script (Step 5) to help test the predictive accuracy
# of our models

# Working directory
setwd("~/Desktop/Publications/Leyerle Martin et al., 2025")

## Packages
# Data formatting and combining
library(plyr)
library(tidyverse)

## Functions
source("PFAS_Review_supportingFunctions.R") # Load supporting functions

## Loading in the finalized data frame created from Step 3, as well as the subset 
# of data points to be used as a supplementary test/validation set for our models
final_imputed_data <- read.csv("Finalized_Imputed_Data_Frame.csv",header = TRUE)
validation_supp<-read.csv("Supplementary_Validation_Data_Frame.csv")

# Select out only the six contaminants in the validation_supp data frame that will
# be modeled in Step 5
validation_supp<-cbind(validation_supp[,1:21],
                       validation_supp[,c("PFOS","PFNA","PFDA",
                                          "PFUnA","PFDoA","PFTrDA")])

# In the following sections, we execute the exact same procedure found in Step 3 to 
# add in several variables (i.e., Water_Level and Water_Level_sc) to the validation_supp 
# data frame. We also create some formatting code for ordering several variables' levels.
# For extended explanations of the code, we ask that the reader refer back to the
# relevant sections in Step 3
#---------------- Additional Variables ----------------------------------------
################# Water Level #################################################
water_level<-read.csv("~/Desktop/Publications/Leyerle Martin et al., 2025/Great Lakes Water Level Data/GL Water Levels (NOAA).csv",
                      header = TRUE)

validation_supp$Water_Level<-NA

for (i in 1:nrow(validation_supp)){
  if(validation_supp$Waterbody[i]=="Lake Superior"){
    for (j in 1:nrow(water_level)){
      if(validation_supp$Sampling.Year[i]==water_level$Year[j]){
        validation_supp$Water_Level[i]<-water_level$Lake.Superior[j]
      }
    }
  }
  else if(validation_supp$Waterbody[i]=="Lake Michigan" 
          || validation_supp$Waterbody[i]=="Lake Huron"){
    for (j in 1:nrow(water_level)){
      if(validation_supp$Sampling.Year[i]==water_level$Year[j]){
        validation_supp$Water_Level[i]<-water_level$Lake.Michigan.Huron[j]
      }
    }
  }
  else if(validation_supp$Waterbody[i]=="Lake Erie"){
    for (j in 1:nrow(water_level)){
      if(validation_supp$Sampling.Year[i]==water_level$Year[j]){
        validation_supp$Water_Level[i]<-water_level$Lake.Erie[j]
      }
    }
  }
  else if(validation_supp$Waterbody[i]=="Lake Ontario"){
    for (j in 1:nrow(water_level)){
      if(validation_supp$Sampling.Year[i]==water_level$Year[j]){
        validation_supp$Water_Level[i]<-water_level$Lake.Ontario[j]
      }
    }
  }
}

unique(validation_supp$Water_Level)

################# Water_Level_sc ##############################################
# To ensure that the scaling procedure for the Water_Level variable in validation_supp
# produces values as close to those of the finalized imputed data frame 
# (final_imputed_data) as possible, we bind final_imputed_data to validation_supp. Thus,
# the group-level calculations of means and standard deviations for water levels will be
# preformed over 2497 samples instead of over 8 samples, and the scaling procedure will
# thereby produce values that are more reflective of the distribution of water levels
# across the larger dataset used in modeling
validation_supp<-rbind.fill(validation_supp,final_imputed_data)

Water_Level_sc <- validation_supp %>% group_by(Waterbody) %>% 
  mutate(Water_Level_sc = scale(Water_Level)) %>% subset(select=Water_Level_sc) %>%
  as.vector() %>% unlist() %>% unname()

validation_supp$Water_Level_sc <- Water_Level_sc

validation_supp %>% group_by(Waterbody) %>% 
  summarise(mean(Water_Level_sc),sd(Water_Level_sc)) 

names(validation_supp)[36]

#---------------- Final Data Frame Editing ------------------------------------
################# Variable Formatting #########################################
# Converting columns to proper format (i.e., ordered factors)
validation_supp$Waterbody<-factor(validation_supp$Waterbody,
                                     levels = c("Lake Superior",
                                                "Lake Michigan",
                                                "Lake Huron",
                                                "Lake Erie",
                                                "Lake Ontario"))

validation_supp$Waterbody_Type<-factor(validation_supp$Waterbody_Type,
                                          levels = c("Inland waters",
                                                     "Connecting channel",
                                                     "Lake"))

validation_supp$Composite<-factor(validation_supp$Composite)

validation_supp$Trophic_Level<-factor(validation_supp$Trophic_Level,
                                         levels = c("Primary Producer",
                                                    "Primary Consumer",
                                                    "Secondary Consumer",
                                                    "Tertiary Consumer",
                                                    "Quaternary Consumer",
                                                    "Piscivorous/Insectivorous Bird",
                                                    "Apex Predator"))

validation_supp$Class<-factor(validation_supp$Class,
                                 levels = c("Algae","Plantae (Magnoliopsida)",
                                            "Annelida","Bivalvia","Gastropoda",
                                            "Zooplankton",
                                            "Shrimp, water fleas, and allies",
                                            "Insecta","Astacoidea","Amphibia",
                                            "Pisces","Reptilia","Aves","Mammalia"))

validation_supp$Revised_Tissue<-factor(validation_supp$Revised_Tissue,
                                          levels = c("Muscle",
                                                     "Whole organism homogenate",
                                                     "Misc. Tissue",
                                                     "Liver","Blood","Eggs"))

################# Data Frame Editing and Reordering ###########################
str(validation_supp)
# Remove extraneous Region variable
validation_supp <- subset(validation_supp,select = -Region)

# Reorder columns to mirror final_imputed_data
validation_supp <- validation_supp[,names(final_imputed_data)]

# Remove final_imputed_data from validation_supp (so that only the 8 supplementary data
# points are left)
validation_supp <- validation_supp[1:8,]

################# Save data frame and delete excess variables #################
write.table(validation_supp,file = "Finalized_Supp_Validation_Data_Frame.csv",
            sep = ",",
            row.names = FALSE)

rm(i,j,water_level,Water_Level_sc)




