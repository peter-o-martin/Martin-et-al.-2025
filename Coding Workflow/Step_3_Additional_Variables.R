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

# The for loop used for classification (some segments had to get very explicitly
# coded since Location values varied between publications and the Location
# variable was the only effective way, that I found, to make the necessary
# distinctions for Waterbody_Type)
for (i in 1:nrow(final_imputed_data)){
  #1 Su et al (2017), #12 Letcher et al (2015), #18 Gebbink et al (2009), 
  #27 Gebbink and Letcher (2010), and #33 Gewurtz et al (2016) are a mix of 
  # Connecting channel and Lake sites
  if(grepl(paste0("#",c(1,12,18,27,33)," ",collapse="|"),
           final_imputed_data$Author..Citation.[i])){
    if(any(final_imputed_data$Location[i]==c("Two Tree Island",
                                             "Pipe Island Twin",
                                             "Five Mile Island",
                                             "Niagara River",
                                             "Fighting Island",
                                             "Turkey Island",
                                             "Weseloh Rocks",
                                             "Strachan Island",
                                             "Swinburn Island",
                                             "Ile Deslauriers",
                                             "Ile Bellechasse")) ||
       final_imputed_data$Sample.ID[i]=="33GWHG 31"){
      final_imputed_data$Waterbody_Type[i]<-"Connecting channel"
    }
    else{
      final_imputed_data$Waterbody_Type[i]<-"Lake"
    }
  }
  #15 Kannan et al (2005) and #37 De Silva et al (2016) are a mix of Connecting 
  # channel, Lake, and Inland water sites
  else if(grepl(paste0("#",c(15,37)," ",collapse="|"),
                final_imputed_data$Author..Citation.[i])){
    if(any(final_imputed_data$Location[i]==c("Saginaw Bay","Thunder Bay",
                                             "Mackinac","Cecil Bay",
                                             "Calumet River, Calumet Park",
                                             "Calumet River, Calumet Harbor",
                                             "Scotch Bonnet Island",
                                             "Hamilton Harbor","Mohawk Island"))){
      final_imputed_data$Waterbody_Type[i]<-"Lake"
    }
    else if(any(final_imputed_data$Location[i]==c("Macomb County, along Lake St. Clair",
                                                  "St. Clair River, Marine City",
                                                  "St. Clair River",
                                                  "Bergin Island",
                                                  "Lac des Deux Montagnes",
                                                  "Îlet Vert"))){
      final_imputed_data$Waterbody_Type[i]<-"Connecting channel"
    }
    else{
      final_imputed_data$Waterbody_Type[i]<-"Inland waters"
    }
  }
  #22 Guo et al (2012), #29 Wu et al (2019), and #56 Hopkins et al (2023) are 
  # a mix of Lake and Inland water sites
  else if(grepl(paste0("#",c(22,29,56)," ",collapse="|"),
                final_imputed_data$Author..Citation.[i])){
    if(any(final_imputed_data$Location[i]==c("Lake Nipigon",
                                             "Inland, Upper Peninsula",
                                             "Mountsberg Conservation Area (Reference)"))){
      final_imputed_data$Waterbody_Type[i]<-"Inland waters"
    }
    else{
      final_imputed_data$Waterbody_Type[i]<-"Lake"
    }
  }
  #25 Dykstra et al (2021) and #64 Route et al (Unpublished) are a mix of Lake 
  # and Inland water sites 
  else if(grepl(paste0("#",c(25,64)," ",collapse="|"),
                final_imputed_data$Author..Citation.[i])){
    if(any(final_imputed_data$Location[i]==c("Little Pokegama Bay",
                                             "RED RIVER, DOUGLAS COUNTY")) ||
       any(final_imputed_data$Region[i]==c("L-SACN","U-SACN")) ||
       final_imputed_data$Sample.ID[i]=="62936500"){
      final_imputed_data$Waterbody_Type[i]<-"Inland waters"
    }
    else{
      final_imputed_data$Waterbody_Type[i]<-"Lake"
    }
  }
  #36 Custer et al (2016) is a mix of Connecting channel, Lake, and Inland
  # water sites
  else if(grepl(paste0("#",36," "),final_imputed_data$Author..Citation.[i])){
    if(any(final_imputed_data$Location[i]==c("MillerCreek","HogIsland",
                                             "TorchLake","MenomineeRiver",
                                             "LittleTailPoint","BayBeach",
                                             "WorkersPark","LakeshorePark",
                                             "Waukegan","OttawaNWRnorth",
                                             "OttawaNWRsouth","BackBay",
                                             "PresqueIsleStPark","PresqueIsleWaterworks",
                                             "KatherineStreet","ChaumontRiver",
                                             "CapeVincent")) ||
       final_imputed_data$Region[i]=="ManistiqueRiver, MI"){
      final_imputed_data$Waterbody_Type[i]<-"Lake"
    }
    else if(any(final_imputed_data$Location[i]==c("AshmunBay")) ||
            any(final_imputed_data$Region[i]==c("StClairRiver, MI",
                                                "DetroitRiver, MI",
                                                "ClintonRiver, MI",
                                                "RougeRiver, MI",
                                                "NiagaraRiver, NY"))){
      final_imputed_data$Waterbody_Type[i]<-"Connecting channel"
    }
    else{
      final_imputed_data$Waterbody_Type[i]<-"Inland waters"
    }
  }
  #65 Kannan et al (2001) is a mix of Connecting channel, Lake, and Inland
  # water sites
  else if(grepl(paste0("#",65," "),final_imputed_data$Author..Citation.[i])){
    if(any(final_imputed_data$Location[i]==c("Hyrn Island, Lake Superior",
                                             "Devil's Is., Lake Superior, WI",
                                             "Otter Island, Lake Superior",
                                             "Rabbit Bay, Houghton County, MI",
                                             "Huron Is., Lake Superior",
                                             "Marquette, MI, Lake Superior",
                                             "Marquette, MI, Lake Superior",
                                             "St. Martin Island",
                                             "Gull Is., Geo Bay",
                                             "Swan Lake, Presque Isle County, Great Lakes, MI",
                                             "Sulphur Is., Thunder Bay, Lake Huron",
                                             "Scarecrow Island, Thunder Bay",
                                             "Little Charity Is., Lake Huron",
                                             "Pt. Movillee, Monroe County, MI",
                                             "White's Landing, OH",
                                             "Sucker Creek, Mackinac, MI"))){
      final_imputed_data$Waterbody_Type[i]<-"Lake"
    }
    else if(any(final_imputed_data$Location[i]==c("Roach Point, Chippewa, MI"))){
      final_imputed_data$Waterbody_Type[i]<-"Connecting channel"
    }
    else{
      final_imputed_data$Waterbody_Type[i]<-"Inland waters"
    }
  }
  #2,5,6,7,9,10,11,17,19,20,26,28,30,34 GLHHFTS,38,47,50,53 are all Lake sites, 
  # Nesting sites adjacent to Lakes, or Nesting sites on islands within 
  # the Lakes (as well as TEST GLHHFTS (2020) and TEST Hopkins et al (2023))
  else if(any(grepl(paste0("#",c(2,5,6,7,9,10,11,17,19,20,26,28,30,38,47,50,53)
                           ," ",collapse="|"),
                    final_imputed_data$Author..Citation.[i])) || 
          final_imputed_data$Author..Citation.[i]=="TEST Hopkins et al (2023)" ||
          final_imputed_data$Author..Citation.[i]=="TEST GLHHFTS (2020)" ||
          any(grepl("#34 GLHHFTS",final_imputed_data$Author..Citation.[i]))){
    final_imputed_data$Waterbody_Type[i]<-"Lake"
  }
  #8,48 are all Connecting channel sites (e.g., St. Lawrence River, 
  # Niagara River)
  else if(any(grepl(paste0("#",c(8,48)," ",collapse="|"),
                    final_imputed_data$Author..Citation.[i]))){
    final_imputed_data$Waterbody_Type[i]<-"Connecting channel"
  }
  #16,21,42,44,46,49,57,58,63 are all Inland waters sites (e.g., tributaries, 
  # inland lakes, inland nesting areas) (as well as TEST Custer et al (2024))
  else if(any(grepl(paste0("#",c(16,21,42,44,46,49,57,58,63)," ",collapse="|"),
                    final_imputed_data$Author..Citation.[i]))){
    final_imputed_data$Waterbody_Type[i]<-"Inland waters"
  }
  # TEST Custer et al (2024) is a mix of Lake and Inland waters sites
  else if(final_imputed_data$Author..Citation.[i]=="TEST Custer et al (2024)"){
    if(final_imputed_data$Location[i]=="Lakeshore Park"){
      final_imputed_data$Waterbody_Type[i]<-"Lake"
    }
    else{
      final_imputed_data$Waterbody_Type[i]<-"Inland waters"
    }
  }
  #34 NRSA FTS dataframes are a mix of Inland waters and Connecting 
  # channel sites
  else if(any(grepl("#34 NRSA",final_imputed_data$Author..Citation.[i]))){
    if(any(final_imputed_data$Location[i]==c("Detroit River",
                                             "Saint Clair River"))){
      final_imputed_data$Waterbody_Type[i]<-"Connecting channel"
    }
    else{
      final_imputed_data$Waterbody_Type[i]<-"Inland waters"
    }
  }
}



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


write.table(final_imputed_data,file = "Step_3_final_imputed_data.csv",sep = ",",
            row.names = FALSE)

rm(water_level,i,j)

