# Step_1_Data_Import_and_Cleaning.R
# Contains all data importing and formatting code for the full data frame

# Written by Peter Martin


# Working directory
setwd("~/Desktop/Rfiles/Leyerle Martin et al., 2025")

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
# library(Hmisc) # capitalize function

## Packages
# Data formatting and combining
library(plyr)
library(dplyr)
library(tidyverse)

# Assigning watershed designations from the GL watershed shapefile
library(sf)
library(s2)

# Visualizing data using maps
library(mapview)
library(leafsync)
library(ggplot2)

############### Formatting and checking for errors ############################
# Compile all .csv files into one data frame
# List of file names used in concatenation
filenames <- 
  list.files(path="~/Desktop/Rfiles/Leyerle Martin et al., 2025/Data Sets", 
                        pattern="#+.*csv") 

source("PFAS_Review_supportingFunctions.R") # Load supporting functions

original_data<-collect.frames(file_names = filenames) # Compiles all data
original_data<-subset(original_data,!is.na(Latitude)) # Removes any NA rows


# write.table(original_data,file = "original_data.csv",sep = ",",row.names = FALSE)



Great_Lakes_watershed <- st_read(
  "~/Desktop/Rfiles/Leyerle Martin et al., 2025/Watershed Shapefile/GL_Watershed_shapefile.shp")

Great_Lakes_watershed$merge[Great_Lakes_watershed$layer==
                              "Gulf of Saint Lawrence"]<-"lk_ont"
Great_Lakes_watershed$merge[Great_Lakes_watershed$layer==
                              "Saint Lawrence Watershed Full"]<-"lk_ont"
Great_Lakes_watershed$merge[Great_Lakes_watershed$merge==
                              "lk_ont"]<-"Lake Ontario"
Great_Lakes_watershed$merge[Great_Lakes_watershed$merge==
                              "lk_erie"]<-"Lake Erie"
Great_Lakes_watershed$merge[Great_Lakes_watershed$merge==
                              "lk_mich"]<-"Lake Michigan"
Great_Lakes_watershed$merge[Great_Lakes_watershed$merge==
                              "lk_huron"]<-"Lake Huron"
Great_Lakes_watershed$merge[Great_Lakes_watershed$merge==
                              "lk_sup"]<-"Lake Superior"

st_geometry_type(Great_Lakes_watershed)
st_is_valid(Great_Lakes_watershed)

Basins_Channels <- st_read(
  "~/Desktop/Publications/Leyerle Martin et al., 2025/Great Lakes Shapefiles/Custom Shapefiles/Individual Waterbodies (Articulated)/Finalized/(Edited) Finalized Lake Basins and Channels.shp")

st_is_valid(Basins_Channels)
LakeB_C<-st_make_valid(Basins_Channels)
st_is_valid(Basins_Channels)

st_crs(Great_Lakes_watershed)

StMarys <- st_read(
  "~/Desktop/Publications/Leyerle Martin et al., 2025/Great Lakes Shapefiles/Custom Shapefiles/Individual Waterbodies (Articulated)/Finalized/St Marys River/Finalized St Marys River (Vector Map, AOC and Natural Earth).shp")
st_is_valid(StMarys)
StMarys[3,]

CRS<-st_crs(Great_Lakes_watershed)
watershed_data <- as.data.frame(original_data) %>% 
  st_as_sf(coords=c("Longitude","Latitude"), crs=4326, remove=FALSE)

ggplot() +
  geom_sf(data = Basins_Channels) +
  geom_sf(data = watershed_data) +
  coord_sf(xlim = c(-75.0,-74.5),ylim = c(44.75,45.25))


# st_crs(Great_Lakes_watershed) <- "EPSG:4326"
sf_use_s2(TRUE)


# Great_Lakes_watershed$geometry <- Great_Lakes_watershed$geometry %>%
#   sf::st_transform(crs = 3857) %>%
#   sf::st_transform(crs = 4326)

a_joined<-st_join(watershed_data, left = T, Basins_Channels["NAME"])
tail(a_joined)


final_a<-subset(a_joined,duplicated(Sample.ID)==FALSE)
setdiff(final_a$Longitude,original_data$Longitude)

final_a<-subset(final_a,is.na(merge)==FALSE)


tail(final_a)
tail(labeled_data)

final_a$Sample.ID==labeled_data$Sample.ID

length(unique(final_a$Sample.ID))
length(unique(labeled_data$Sample.ID))
setdiff(final_a$Sample.ID,labeled_data$Sample.ID)

ID <- final_a$Sample.ID %in% unique(labeled_data$Sample.ID)
excluded <- final_a[!ID,]


excluded<-subset(final_a,Sample.ID==any(samples))
mapview(excluded[,1:20], xcol = "Longitude", ycol = "Latitude",
        crs = 4269, grid = FALSE)


ggplot() +
  geom_sf(data=Great_Lakes_watershed) +
  geom_jitter(data = excluded,aes(x=Longitude,y=Latitude))  +
  coord_sf(xlim = c(-76.0,-74),ylim = c(44,46))
  

sum(is.na(a_joined$merge))



labeled_data<-subset(labeled_data,is.na(Waterbody)==TRUE)
labeled_data<-subset(labeled_data,duplicated(Sample.ID)==FALSE)

unique(final_a$Sample.ID)==unique(labeled_data$Sample.ID)

final_a<-as.data.frame(final_a)

# 
final_a<-subset(final_a,duplicated(Sample.ID)==FALSE)

str(final_a)
colnames(final_a)[20:21]<-c("Sample Size","n")

mapview(final_a, xcol = "Longitude", ycol = "Latitude",
        crs = 4269, grid = FALSE)

nrow(subset(final_a,is.na(merge)))
nrow(subset(labeled_data,Final.Waterbody=="  "))

b<-final_a[final_a$merge!=labeled_data$Waterbody,]





final_a<-subset(final_a,is.na(merge)==FALSE)
final_a<-subset(final_a,duplicated(Sample.ID)==FALSE)
unique(labeled_data$Waterbody)
labeled_data<-labeled_data[,c(1:102)]

# Checking that the functions performed correctly and that all labels 
# were properly transferred over
sum(original_data$Waterbody=="Lake Huron")
sum(labeled_data$Final.Waterbody=="  lk_huron")
sum(labeled_data$Waterbody=="Lake Huron")

# Found an annoying "hyphen instead of a dash" negative sign problem in the
# code, so this for loop just corrects that mistake for the four entries
# Rows 1134, 1135, 1136, and 1139
for (j in 8:9){
  final_a[j]<-lapply(final_a[j], gsub, pattern = "–", 
                          replacement = "-")
}

# Collect test data frames together and bind them to the labeled_data variable
test_filenames <- list.files(path="~/Desktop/Rfiles/PFAS Review Paper", 
                             pattern="TEST+.*csv") 
test_data<-collect.frames(file_names = test_filenames) # Compiles all data
test_data<-subset(test_data,is.na(Latitude)==FALSE) # Removes any NA rows

for (j in 7:8){
  test_data[j]<-lapply(test_data[j], gsub, pattern = "–", 
                       replacement = "-")
}
final_a[,7:8]<-lapply(final_a[,7:8],as.numeric)
colnames(test_data)[12:13]<-c("N","n_1")

labeled_data<-rbind.fill(labeled_data,test_data)
rownames(labeled_data)<-NULL

b<-cbind(final_a$Waterbody,final_a$merge,labeled_data$Waterbody)


# Mapping of watershed divisions in R (mapview function)




################# Post-QGIS watershed assignments, load resulting .csv file ####
# Load in data that has had watershed splits from QGIS applied
QGIS_data<-read.csv("QGIS_data.csv",header = TRUE)

# Reorder data file so it matches the original (original_data), and remove any 
# duplicates produced by QGIS
QGIS_data<-QGIS_data[order(match(QGIS_data[,2],original_data[,2])),]
QGIS_data<-subset(QGIS_data,duplicated(Sample.ID)==FALSE)

# Keep variables generated by QGIS that describe the watershed designations, and
# then merge the contents of these variables into one column: Final.Waterbody
labeled_data<-cbind(QGIS_data[,1:102],QGIS_data$name,QGIS_data$CAN_MDA_EN,
                    QGIS_data$merge,deparse.level = 0)
labeled_data$Final.Waterbody<-paste(labeled_data$`QGIS_data$name`,
                                    labeled_data$`QGIS_data$CAN_MDA_EN`,
                                    labeled_data$`QGIS_data$merge`)
unique(labeled_data$Final.Waterbody)

# Remove QGIS variables and keep only Final.Waterbody column
labeled_data<-labeled_data[,c(1:102,106)]
labeled_data$Final.Waterbody

# Change QGIS-derived classifications to the watershed/waterbody designations
# that we have decided on
for(i in 1:nrow(labeled_data)){
  if(labeled_data$Final.Waterbody[i]=="  lk_huron"){
    labeled_data$Waterbody[i]<-"Lake Huron"
  }
  else if(labeled_data$Final.Waterbody[i]=="  lk_erie"){
    labeled_data$Waterbody[i]<-"Lake Erie"
  }
  else if(labeled_data$Final.Waterbody[i]=="  lk_sup"){
    labeled_data$Waterbody[i]<-"Lake Superior"
  }
  else if(labeled_data$Final.Waterbody[i]=="  lk_ont"){
    labeled_data$Waterbody[i]<-"Lake Ontario"
  }
  else if(labeled_data$Final.Waterbody[i]=="  lk_mich"){
    labeled_data$Waterbody[i]<-"Lake Michigan"
  }
  else if(labeled_data$Final.Waterbody[i]==" St. Lawrence Drainage Area " || 
          labeled_data$Final.Waterbody[i]=="Gulf of St. Lawrence  "){
    labeled_data$Waterbody[i]<-"Lake Ontario"
  }
  else if(labeled_data$Final.Waterbody[i]=="  "){
    labeled_data$Waterbody[i]<-NA
  }
}

# Remove all points that fall outside the shapefile in QGIS (i.e., those points
# for which Final.Waterbody==" " and for which is.na(Waterbody)==TRUE)
# Also remove duplicates, and delete the Final.Waterbody column (superfluous now)
labeled_data<-subset(labeled_data,is.na(Waterbody)==FALSE)
labeled_data<-subset(labeled_data,duplicated(Sample.ID)==FALSE)
unique(labeled_data$Waterbody)
labeled_data<-labeled_data[,c(1:102)]

# Checking that the functions performed correctly and that all labels 
# were properly transferred over
sum(original_data$Waterbody=="Lake Huron")
sum(labeled_data$Final.Waterbody=="  lk_huron")
sum(labeled_data$Waterbody=="Lake Huron")

# Mapping of watershed divisions in R (mapview function)
mapview(labeled_data, xcol = "Longitude", ycol = "Latitude", zcol="Waterbody",
        crs = 4269, grid = FALSE)

# Found an annoying "hyphen instead of a dash" negative sign problem in the
# code, so this for loop just corrects that mistake for the four entries
# Rows 1134, 1135, 1136, and 1139
for (j in 8:9){
  labeled_data[j]<-lapply(labeled_data[j], gsub, pattern = "–", 
                          replacement = "-")
}

# Collect test data frames together and bind them to the labeled_data variable
test_filenames <- list.files(path="~/Desktop/Rfiles/PFAS Review Paper", 
                             pattern="TEST+.*csv") 
test_data<-collect.frames(file_names = test_filenames) # Compiles all data
test_data<-subset(test_data,is.na(Latitude)==FALSE) # Removes any NA rows

for (j in 7:8){
  test_data[j]<-lapply(test_data[j], gsub, pattern = "–", 
                       replacement = "-")
}
test_data[,7:8]<-lapply(test_data[,7:8],as.numeric)
colnames(test_data)[12:13]<-c("N","n_1")

labeled_data<-rbind.fill(labeled_data,test_data)
rownames(labeled_data)<-NULL

mapview(labeled_data, xcol = "Longitude", ycol = "Latitude", zcol="Waterbody",
        crs = 4269, grid = FALSE)






