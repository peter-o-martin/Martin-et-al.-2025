# Step_3_Additional_Variables.R
# Contains all code used to add the additional variables

# Written by Peter Martin
# Created December 13, 2024
# Finalized December 13, 2024

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
# library(DHARMa)
# library(Matrix)
# library(TMB)
# library(glmmTMB)
# 

## Packages
# Data formatting and combining
library(stringr)
library(Hmisc) # capitalize function
# # library(plyr)
# # library(dplyr)
# # library(tidyverse)
# 
# # Imputation
# library(MASS)
# library(survival)
# library(NADA)
# library(truncnorm)
# library(zCompositions)
# 
# # 
# Assigning watershed designations from the GL watershed shapefile
library(sf)
library(s2)
# # 
# Visualizing data using maps
library(mapview)
library(leafsync)
# # library(ggplot2)
# # library(rnaturalearth)

## Functions
source("PFAS_Review_supportingFunctions.R") # Load supporting functions

## Loading in data frame created from Step 1
final_imputed_data <- read.csv("Step_2_final_imputed_data.csv",header = TRUE)

################# Additional Variables ########################################
# Add in a seasonality column: Summer (May-Oct) and Winter (Nov-Apr)
# as defined by Beletsky, Saylor and Schwab, 1999 (the two "major dynamical
# regimes: stratified (May through October), and isothermal (November through 
# April)")
final_imputed_data$Season<-format(as.Date(final_imputed_data$Collection.Date,
                                          "%m/%d/%Y"), "%m")
final_imputed_data$Season<-sapply(final_imputed_data$Season,as.numeric)
for (i in 1:nrow(final_imputed_data)){
  if(is.na(final_imputed_data$Collection.Date[i])==FALSE){
    if(any(final_imputed_data$Season[i]==5:10)){
      final_imputed_data$Season[i]<-"Summer"
    }
    else{
      final_imputed_data$Season[i]<-"Winter"
    }
  }
}

unique(final_imputed_data$Season)

# Add in a revised tissue column: the composite classification is removed, and
# overly detailed tissue descriptions are generalized
final_imputed_data$Revised_Tissue<-final_imputed_data$Tissue
for (i in 1:nrow(final_imputed_data)){
  if(grepl("Composite",final_imputed_data$Tissue[i])==TRUE){
    final_imputed_data$Revised_Tissue[i]<-substr(final_imputed_data$Tissue[i],
                                                 11,str_length(final_imputed_data$Tissue[i]))
    final_imputed_data$Revised_Tissue[i]<-capitalize(final_imputed_data$Revised_Tissue[i])
  }
}
for (i in 1:nrow(final_imputed_data)){
  if(grepl("fillet",final_imputed_data$Tissue[i])==TRUE || 
     grepl("Fillet",final_imputed_data$Tissue[i])==TRUE){
    final_imputed_data$Revised_Tissue[i]<-"Muscle"
  }
}
for (i in 1:nrow(final_imputed_data)){
  if(any(final_imputed_data$Tissue[i]==c("Whole fish homogenate","Biofilm",
                                         "Bulk composite",
                                         "Plant sample homogenate",
                                         "Composite whole fish homogenate"))){
    final_imputed_data$Revised_Tissue[i]<-"Whole organism homogenate"
  }
}
for (i in 1:nrow(final_imputed_data)){
  if(grepl("uscle",final_imputed_data$Tissue[i])==TRUE){
    final_imputed_data$Revised_Tissue[i]<-"Muscle"
  }
}
for (i in 1:nrow(final_imputed_data)){
  if(grepl("Egg",final_imputed_data$Revised_Tissue[i])==TRUE){
    final_imputed_data$Revised_Tissue[i]<-"Eggs"
  }
}
for (i in 1:nrow(final_imputed_data)){
  if(any(final_imputed_data$Revised_Tissue[i]==c("Plasma","Red blood cells",
                                                 "Whole blood","Serum"))){
    final_imputed_data$Revised_Tissue[i]<-"Blood"
  }
}
for (i in 1:nrow(final_imputed_data)){
  if(any(final_imputed_data$Revised_Tissue[i]==c("Diet pool","Ovary","Adipose tissue",
                                                 "Gall bladder","Kidney","Brain",
                                                 "Carcass","Testes","Hepatopancreas",
                                                 "GI-tract"))){
    final_imputed_data$Revised_Tissue[i]<-"Misc. Tissue"
  }
}

unique(final_imputed_data$Revised_Tissue)

# Add in a composite column: Yes or No depending on whether the Tissue 
# variable entry indicates a composite/biofilm/bulk composite
final_imputed_data$Composite<-NA
for (i in 1:nrow(final_imputed_data)){
  if(grepl("Composite",final_imputed_data$Tissue[i])==TRUE || 
     grepl("composite",final_imputed_data$Tissue[i])==TRUE || 
     grepl("Biofilm",final_imputed_data$Tissue[i])==TRUE){
    final_imputed_data$Composite[i]<-"Yes"
  }
  else{
    final_imputed_data$Composite[i]<-"No"
  }
}

unique(final_imputed_data$Composite)

# Add in a trophic level column: 
# Primary Producer (e.g., aquatic macrophytes and algae)

# Primary Consumer (e.g., insects, shrimp, amphipods, filter feeders and other
# species that rely on primary production)

# Secondary Consumer (e.g., minnows, goby, dace, redhorse, green frog, smaller-
# bodied bottom feeders, fish that are primarily insectivorous/planktivores/
# detritivores)

# Tertiary Consumer (e.g., most smaller sunfish and perch, as well as smelt 
# and perch, that are smaller bodied piscivorous fish that prey upon a mix of 
# detritivorous/herbivorous primary and secondary consumers. Also included are 
# larger bodied, more omnivorous bottom-feeders)

# Quaternary Consumer (e.g., the piscivorous fish, like salmon, trout, bass,
# pike, pickerel, burbot and walleye)

# Piscivorous/Insectivorous Bird (e.g., gulls, cormorants, herons and swallows, birds that
# are either piscivorous or that represent a synthesis of terrestrial and 
# aquatic sources of PFAS in the watershed while still being beneath the apex
# predators in the foodweb)

# Apex Predator (e.g., bald eagles and mink, the top of the food chain and
# high rate of omnivory. These animals are a synthesis of terrestrial and 
# aquatic sources of PFAS in the watershed)

final_imputed_data$Trophic_Level<-NA

for (i in 1:nrow(final_imputed_data)){
  if(any(final_imputed_data$Class[i]==c("Plantae (Magnoliopsida)","Algae"))){
    final_imputed_data$Trophic_Level[i]<-"Primary Producer"
  }
  else if(any(final_imputed_data$Class[i]==c("Insecta","Zooplankton",
                                             "Shrimp, water fleas, and allies",
                                             "Bivalvia","Annelida","Amphipoda",
                                             "Gastropoda"))){
    final_imputed_data$Trophic_Level[i]<-"Primary Consumer"
  }
  else if(final_imputed_data$Class[i]=="Astacoidea" || 
          grepl("shiner",final_imputed_data$Species[i]) ||
          any(final_imputed_data$Species[i]==c("General forage fish larvae",
                                               "Round goby","Green frog",
                                               "Gizzard shad","Blacknose dace",
                                               "Trout-perch","Freshwater drum",
                                               "Sicklefin redhorse",
                                               "Silver redhorse",
                                               "Bluntnose minnow","Bloater",
                                               "Tench","Cisco")) || 
          grepl("sucker",final_imputed_data$Species[i]) ||
          grepl("sculpin",final_imputed_data$Species[i])){
    final_imputed_data$Trophic_Level[i]<-"Secondary Consumer"
  }
  else if(any(final_imputed_data$Species[i]==c("Common carp",
                                               "Channel catfish",
                                               "Rainbow smelt","Alewife",
                                               "Sunfish","Pumpkinseed",
                                               "Bluegill",
                                               "Yellow perch","White perch")) ||
          grepl("whitefish",final_imputed_data$Species[i]) ||
          grepl("ullhead",final_imputed_data$Species[i]) ||
          grepl("crappie",final_imputed_data$Species[i])){
    final_imputed_data$Trophic_Level[i]<-"Tertiary Consumer"
  }
  else if(any(final_imputed_data$Species[i]==c("Walleye","Splake",
                                               "Chain pickerel","Burbot",
                                               "Northern pike")) ||
          grepl("bass",final_imputed_data$Species[i]) ||
          grepl(" salmon",final_imputed_data$Species[i]) || 
          grepl("trout",final_imputed_data$Species[i])){
    final_imputed_data$Trophic_Level[i]<-"Quaternary Consumer"
  }
  else if(any(final_imputed_data$Species[i]==c("Great blue heron",
                                               "Tree swallow",
                                               "Double-crested cormorant",
                                               "Caspian tern")) ||
          grepl("gull",final_imputed_data$Species[i])){
    final_imputed_data$Trophic_Level[i]<-"Piscivorous/Insectivorous Bird"
  }
  else if(any(final_imputed_data$Species[i]==c("Mink","Snapping turtle",
                                               "Bald eagle","Peregrine falcon"))){
    final_imputed_data$Trophic_Level[i]<-"Apex Predator"
  }
}

unique(final_imputed_data$Trophic_Level)

# Add a Waterbody_Type column
# Lake: Sample that is physically within one of the lake basins

# Inland waters: Sample that is located anywhere within the bounds of a 
# tributary, inland lake, stream or creek that is not one of the Laurentian 
# Great Lakes

# Connecting channel: for example, Niagara River, St. Lawrence River, Lake St.
# Clair, Detroit River, or Saint Marys River, bodies of water/rivers that 
# connect one of the lake basins to another or to the ocean (i.e., St. Lawrence
# River)

final_imputed_data$Waterbody_Type<-NA

Basins_Channels <- st_read(
  "~/Desktop/Publications/Leyerle Martin et al., 2025/Great Lakes Shapefiles/Custom Shapefiles/Individual Waterbodies (Articulated)/Finalized/(Edited) Finalized Lake Basins and Channels.shp")
st_is_valid(Basins_Channels)

CRS<-st_crs(Basins_Channels) # Save the coordinate reference system (CRS) from the shapefile so that the data frame can inherit the same CRS

# Convert data frame to a sf object (using CRS inherited from the Great_Lakes_watershed shapefile)
geo_data <- st_as_sf(final_imputed_data, coords=c("Longitude","Latitude"), crs=CRS, remove=FALSE)
setdiff(geo_data$Longitude,final_imputed_data$Longitude) # Check that exact coordinate values were not lost in converting the original data frame to a sf object

sf_use_s2(TRUE) # Turns on the s2 package to construct spherical geometries from the polygons

# Perform the spatial join
geo_data<-st_join(geo_data, left = T, Basins_Channels["NAME"])

# Remove duplicates generated by the spatial join procedure
geo_data<-subset(geo_data,duplicated(Sample.ID)==FALSE)
length(unique(geo_data$Sample.ID))
unique(geo_data$NAME)

final_imputed_data<-st_drop_geometry(geo_data) # Removes the sf designation

# The for loop used for Waterbody_Type classification of data points
for(i in 1:nrow(final_imputed_data)){
  if(grepl("Lake",final_imputed_data$NAME[i]) && 
     final_imputed_data$NAME[i]!="Lake Saint Clair"){
    final_imputed_data$Waterbody_Type[i]<-"Lake"
  }
  else if(is.na(final_imputed_data$NAME[i])){
    final_imputed_data$Waterbody_Type[i]<-"Inland waters"
  }
  else {
    final_imputed_data$Waterbody_Type[i]<-"Connecting channel"
  }
}

final_imputed_data<-subset(final_imputed_data, select = -NAME)

unique(final_imputed_data$Waterbody_Type)

# Tool for examining how samples in a given publication were classified
mapview(final_imputed_data,xcol = "Longitude", ycol = "Latitude", 
        zcol="Waterbody_Type",
        crs = 4326, grid = FALSE)


# Add in a revised species column: overly detailed tissue descriptions are 
# generalized (e.g., "Faxonius rusticus" and "Astacoidea crayfish" are 
# simplified to one "Crayfish" label)

final_imputed_data$Revised_Species<-final_imputed_data$Species

for (i in 1:nrow(final_imputed_data)){
  if(grepl("mussel",final_imputed_data$Species[i])==TRUE){
    final_imputed_data$Revised_Species[i]<-"Dreissenid mussel"
  }
  else if(grepl("Fax",final_imputed_data$Species[i])==TRUE ||
          grepl("rayfish",final_imputed_data$Species[i])==TRUE){
    final_imputed_data$Revised_Species[i]<-"Crayfish"
  }
  else if(grepl("Diporeia",final_imputed_data$Species[i])==TRUE){
    final_imputed_data$Revised_Species[i]<-"Diporeia spp."
  }
  else if(grepl("ullhead",final_imputed_data$Species[i])==TRUE){
    final_imputed_data$Revised_Species[i]<-"Bullhead (Ameiurus sp.)"
  }
  # Gammaridae and Hyalellidae species compressed into "Amphipods" designation
  else if(grepl("Amphipods",final_imputed_data$Species[i])==TRUE){
    final_imputed_data$Revised_Species[i]<-"Amphipods (principally Gammarus and Hyalella)"
  }
  else if(any(final_imputed_data$Species[i]==c("Mayflies","Damselfly"))){
    final_imputed_data$Revised_Species[i]<-"Aquatic insects"
  }
  else if(any(final_imputed_data$Species[i]==c("Caridea shrimp","Mysis relicta"))){
    final_imputed_data$Revised_Species[i]<-"Freshwater shrimp (Caridea and Mysida)"
  }
}

unique(final_imputed_data$Revised_Species)


# Add in a revised maturity column: overly detailed maturity descriptions are 
# generalized, and three categories are formed: Adult (Reproductive), 
# Adult (Non-reproductive) and Juvenile/Immature

final_imputed_data$Revised_Maturity<-final_imputed_data$Maturity

# Gives Maturity category percentages and counts
final_imputed_data %>% group_by(Maturity) %>% 
  summarise(Percentage=(n()/nrow(.))*100,Count=n())

for (i in 1:nrow(final_imputed_data)){
  if(!is.na(final_imputed_data$Maturity[i])){
    if(any(final_imputed_data$Maturity[i]==c("Juvenile (YOY)","Larvae",
                                             "Nestling","Juvenile"))){
      final_imputed_data$Revised_Maturity[i]<-"Juvenile/Immature"
    }
    else if(final_imputed_data$Maturity[i]=="Adult (Sexually Mature)"){
      final_imputed_data$Revised_Maturity[i]<-"Adult (Reproductive)"
    }
    else if(final_imputed_data$Maturity[i]=="YAO"){
      final_imputed_data$Revised_Maturity[i]<-NA
    }
  }
}

unique(final_imputed_data$Revised_Maturity)


# Add in a Water_Level variable, with values that are specific to each lake in
# each year
# Source: https://www.lre.usace.army.mil/Missions/Great-Lakes-Information/Great-Lakes-Information-2/Water-Level-Data/
# and https://www.glerl.noaa.gov/data/wlevels/dashboard/#longterm

water_level<-read.csv("~/Desktop/Publications/Leyerle Martin et al., 2025/Great Lakes Water Level Data/GL Water Levels (NOAA).csv",header = TRUE)
final_imputed_data$Water_Level<-NA

for (i in 1:nrow(final_imputed_data)){
  if(final_imputed_data$Waterbody[i]=="Lake Superior"){
    for (j in 1:nrow(water_level)){
      if(final_imputed_data$Sampling.Year[i]==water_level$Year[j]){
        final_imputed_data$Water_Level[i]<-water_level$Lake.Superior[j]
      }
    }
  }
  else if(final_imputed_data$Waterbody[i]=="Lake Michigan" 
          || final_imputed_data$Waterbody[i]=="Lake Huron"){
    for (j in 1:nrow(water_level)){
      if(final_imputed_data$Sampling.Year[i]==water_level$Year[j]){
        final_imputed_data$Water_Level[i]<-water_level$Lake.Michigan.Huron[j]
      }
    }
  }
  else if(final_imputed_data$Waterbody[i]=="Lake Erie"){
    for (j in 1:nrow(water_level)){
      if(final_imputed_data$Sampling.Year[i]==water_level$Year[j]){
        final_imputed_data$Water_Level[i]<-water_level$Lake.Erie[j]
      }
    }
  }
  else if(final_imputed_data$Waterbody[i]=="Lake Ontario"){
    for (j in 1:nrow(water_level)){
      if(final_imputed_data$Sampling.Year[i]==water_level$Year[j]){
        final_imputed_data$Water_Level[i]<-water_level$Lake.Ontario[j]
      }
    }
  }
}

unique(final_imputed_data$Water_Level)


# Reorder columns so that contaminants come after all other variables
final_imputed_data<-
  final_imputed_data[,c(colnames(final_imputed_data[,1:6]),"Revised_Maturity",
                        colnames(final_imputed_data[,7:10]),"Water_Level",
                        "Waterbody_Type",
                        colnames(final_imputed_data[,11:18]),"Revised_Species",
                        "Trophic_Level","Tissue","Revised_Tissue","Composite",
                        colnames(final_imputed_data[,20:22]),"Season",
                        colnames(final_imputed_data[,23:33]))]

# Finally, make sure that all columns register properly as characters or
# as numeric
final_imputed_data[,9:10]<-lapply(final_imputed_data[,9:10],as.numeric)


write.table(final_imputed_data,file = "Step_3_full_imputed_df.csv",sep = ",",
            row.names = FALSE)



