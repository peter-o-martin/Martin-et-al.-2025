# Step_3_Additional_Variables.R
# Contains all code used to add the additional variables

# Written by Peter Martin
# Created December 13, 2024
# Finalized December 13, 2024

# Working directory
setwd("~/Desktop/Publications/Leyerle Martin et al., 2025")

## Packages
# Data formatting and combining
library(stringr)
library(Hmisc) # capitalize function
library(tidyverse)

# Assigning watershed designations from the GL watershed shapefile
library(sf)
library(s2)

# Visualizing data using maps
library(mapview)
library(leafsync)

## Functions
source("PFAS_Review_supportingFunctions.R") # Load supporting functions

final_imputed_data <- read.csv("Step_3_final_imputed_data.csv",header = TRUE)

validation_supp<-read.csv("Supplementary Validation Data.csv")
supp_samples<-unique(validation_supp$Sample.ID)

validation_supp<-cbind(validation_supp[,1:21],
                       validation_supp[,c("PFOS","PFNA","PFDA",
                                          "PFUnA","PFDoA","PFTrDA")])

water_level<-read.csv("~/Desktop/Publications/Leyerle Martin et al., 2025/Great Lakes Water Level Data/GL Water Levels (NOAA).csv",header = TRUE)

validation_supp<-rbind.fill(validation_supp,
                            subset(final_imputed_data,Waterbody=="Lake Ontario"))

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

validation_supp <- validation_supp %>% group_by(Waterbody) %>% 
  mutate(Water_Level_sc = scale(Water_Level))

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

validation_supp<-validation_supp[1:8,]

write.table(validation_supp,file = "Step_4_validation_supp.csv",sep = ",",
            row.names = FALSE)

rm(i,j,supp_samples,water_level)




