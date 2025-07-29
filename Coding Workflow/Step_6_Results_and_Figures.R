# Step_6_Results_and_Figures.R
# Contains all code used to extract results from models and construct tables and figures

# Written by Peter O. Martin (https://orcid.org/0009-0009-9070-9200)

# Working directory
setwd("~/Desktop/Publications/Leyerle Martin et al., 2025")

## Packages
# Data formatting and combining
library(tidyverse)
library(gt)
library(webshot2)
library(caret)
 
# Modeling
library(mgcv)

# Visualizing model output and calculating results
library(emmeans)
library(ggrepel)
library(ggpubr)
library(rnaturalearth)
library(sf)
library(terra)
library(ggspatial)
library(RColorBrewer)
library(ggpattern)
library(patchwork)
library(scales)
library(ggpmisc)

## Functions
source("PFAS_Review_supportingFunctions.R") # Load supporting functions

## Loading in finalized data frame from Step 5
final_imputed_data <- read.csv("Finalized_Imputed_Data_Frame.csv",header = TRUE)

PFOS<-final_imputed_data[is.na(final_imputed_data$PFOS)==FALSE,]

## Loading in the models from Step 5
load("Models/full_PFOS_gam.Rdata") 
load("Models/full_PFNA_gam.Rdata")
load("Models/full_PFDA_gam.Rdata")
load("Models/full_PFUnA_gam.Rdata")
load("Models/full_PFDoA_gam.Rdata")
load("Models/full_PFTrDA_gam.Rdata")


model_list <- list(full_PFOS_gam,full_PFNA_gam,full_PFDA_gam,full_PFUnA_gam,
                   full_PFDoA_gam,full_PFTrDA_gam)

WB_level_order<-c("Lake Superior","Lake Michigan","Lake Huron",
                  "Lake Erie","Lake Ontario")

# color palette used
display.brewer.pal(n=8,"RdYlBu")

Great_Lakes_region <- ne_states(country=c("canada","united states of america"),
                                returnclass = "sf")
Great_Lakes_watershed <-
  read_sf("~/Desktop/Publications/Leyerle Martin et al., 2025/Great Lakes Shapefiles/Custom Shapefiles/Full Watershed Great Lakes",
                                 "GL_Watershed_shapefile")

######## Tables ###############################################################
######## Table 1 ##############################################################
# Code to generate the estimated marginal means and pairwise contrasts reported
# in Table 1, where full_PFOS_gam, full_PFNA_gam, full_PFDA_gam, full_PFUnA_gam,
# full_PFDoA_gam, and full_PFTrDA_gam are plugged into the emmeans() function
emmeans(full_PFOS_gam, 
        specs = pairwise ~ Waterbody,
        type = "response",tran = "log10",adjust="tukey")

# -----------------------------------------------------------------------------
######## Table 2 ##############################################################
# Code for generating model-estimated mean concentrations and pairwise 
# contrasts reported in Table 2, where full_PFOS_gam, full_PFNA_gam, 
# full_PFDA_gam, full_PFUnA_gam, full_PFDoA_gam, and full_PFTrDA_gam are plugged
# into the emmeans() function
emmeans(full_PFOS_gam, 
        specs = pairwise ~ Revised_Tissue,
        type = "response",tran = "log10",adjust="tukey")

# Code for calculating the relative distributions (%) of each PFAS across
# four different tissue groups (Eggs, Blood, Liver, and Combined Tissue and
# Blood [i.e., all levels of the variable Revised_Tissue except Eggs])

# Create the matrix that will be filled with the calculated percentage values
Table_2_percent <- matrix(nrow = 6, ncol = 6,
                  dimnames = list(
                    c("Eggs","Blood","Liver","Combined Tissue and Blood",
                      "Muscle","Whole Organism Homogenate"),
                    c("PFOS (C8)","PFNA (C9)","PFDA (C10)","PFUnA (C11)",
                     "PFDoA (C12)","PFTrDA (C13)")
                    ))

# for loop to make calculations for each of the six models in the model_list object
for(i in 1:length(model_list)){
  # Generate output from emmeans
  tissue_estimates <- emmeans(model_list[[i]],
                              specs = ~ Revised_Tissue,
                              type='response',tran = "log10") |>
    as_tibble()
  
  # Calculate percentages and store them in a new column labeled as t_abundance
  tissue_estimates <- tissue_estimates %>% select(Revised_Tissue,response) %>% 
    mutate(t_abundance = 100*round(response/sum(response),digits=5))
  
  # Fill in the correct column of the matrix with the calculated percentages
  Table_2_percent[1:3,i] <- tissue_estimates$t_abundance[6:4]
  Table_2_percent[4,i] <- sum(tissue_estimates$t_abundance[1:5])
  Table_2_percent[5,i] <- tissue_estimates$t_abundance[1]
  Table_2_percent[6,i] <- tissue_estimates$t_abundance[2]
}

# Convert to tibble format
Table_2_percent <- as_tibble(Table_2_percent)

# Label rows with the correct tissue type and create a new column called group
# that can be used by the function gt() to label our final table
Table_2_percent$rowname <- c("Eggs","Blood","Liver","Combined Tissue and Blood",
                             "Muscle","Whole Organism Homogenate")
Table_2_percent$group <- "Percent Composition"

# Produce the finalized section of Table 2 using the function gt() and save
# this output as a .rtf file with the function gtsave()
Table_2_final<-gt(Table_2_percent,
   rowname_col = "rowname",
   groupname_col = "group",
   row_group_as_column = F) |> 
  tab_stubhead(label = "Tissue Type") |>
  cols_align(
    align = "center",
    columns = everything()
  ) |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) |>
  opt_table_font(
    font = "Arial",
    weight = 350,
    size = 13) |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_row_groups()
  ) %>%
  fmt_markdown(columns = everything())

Table_2_final
gtsave(Table_2_final,"Tables and Figures/Table_2.rtf")
# -----------------------------------------------------------------------------
######## Table S4 #############################################################
# Code to generate the model summary statistics and Wald test of significance
# reported in Table S4, where full_PFOS_gam, full_PFNA_gam, full_PFDA_gam,
# full_PFUnA_gam, full_PFDoA_gam, and full_PFTrDA_gam are plugged into the 
# summary() and anova() functions
summary(full_PFOS_gam)
anova(full_PFOS_gam)

# -----------------------------------------------------------------------------
######## Table S5 #############################################################
# Code for generating model-estimated mean concentrations and pairwise 
# contrasts reported in Table S5, where full_PFOS_gam, full_PFNA_gam, 
# full_PFDA_gam, full_PFUnA_gam, full_PFDoA_gam, and full_PFTrDA_gam are plugged
# into the emmeans() function
emmeans(full_PFOS_gam, 
        specs = pairwise ~ Waterbody_Type,
        type = "response",tran = "log10",adjust="tukey")

# Code for calculating the relative distributions (%) of each PFAS across
# three different waterbody type groups (Inland waters, Lake, and Connecting 
# channel)

# Create the matrix that will be filled with the calculated percentage values
Table_S5_percent <- matrix(nrow = 3, ncol = 6,
                          dimnames = list(
                            c("Inland waters","Lake","Connecting channel"),
                            c("PFOS (C8)","PFNA (C9)","PFDA (C10)","PFUnA (C11)",
                              "PFDoA (C12)","PFTrDA (C13)")
                          ))

# for loop to make calculations for each of the six models in the model_list object
for(i in 1:length(model_list)){
  # Generate output from emmeans
  wt_estimates <- emmeans(model_list[[i]],
                          specs = ~ Waterbody_Type,
                          type='response',tran = "log10") |>
    as_tibble()
  
  # Calculate percentages and store them in a new column labeled as wt_abundance
  wt_estimates <- wt_estimates %>% select(Waterbody_Type,response) %>% 
    mutate(wt_abundance = 100*round(response/sum(response),digits=5))
  
  # Fill in the correct column of the matrix with the calculated percentages
  Table_S5_percent[1:3,i] <- wt_estimates$wt_abundance[c(1,3,2)]
}

# Convert to tibble format
Table_S5_percent <- as_tibble(Table_S5_percent)

# Label rows with the correct tissue type and create a new column called group
# that can be used by the function gt() to label our final table
Table_S5_percent$rowname <- c("Inland waters","Lake","Connecting channel")
Table_S5_percent$group <- "Percent Composition"

# Produce the finalized section of Table S5 using the function gt() and save
# this output as a .rtf file with the function gtsave()
Table_S5_final<-gt(Table_S5_percent,
                  rowname_col = "rowname",
                  groupname_col = "group",
                  row_group_as_column = F) |> 
  tab_stubhead(label = "Waterbody Type") |>
  cols_align(
    align = "center",
    columns = everything()
  ) |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) |>
  opt_table_font(
    font = "Arial",
    weight = 350,
    size = 13) |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_row_groups()
  ) %>%
  fmt_markdown(columns = everything())

Table_S5_final
gtsave(Table_S5_final,"Tables and Figures/Table_S5.rtf")
# -----------------------------------------------------------------------------
######## Table S6 #############################################################
# Code to generate the estimated marginal means and pairwise contrasts reported
# in Table S6, where full_PFOS_gam, full_PFNA_gam, full_PFDA_gam, full_PFUnA_gam,
# full_PFDoA_gam, and full_PFTrDA_gam are plugged into the emmeans() function
emmeans(full_PFOS_gam, 
        specs = pairwise ~ Trophic_Level,
        type = "response",tran = "log10",adjust="tukey")

# Code for calculating the relative distributions (%) of each PFAS across
# three different waterbody type groups (Inland waters, Lake, and Connecting 
# channel)

# Create the matrix that will be filled with the calculated percentage values
Table_S6_percent <- matrix(nrow = 7, ncol = 6,
                           dimnames = list(
                             c("Primary Producer","Primary Consumer",
                               "Secondary Consumer",
                               "Tertiary Consumer","Quaternary Consumer",
                               "Piscivorous/Insectivorous Bird",
                               "Apex Predator"),
                             c("PFOS (C8)","PFNA (C9)","PFDA (C10)","PFUnA (C11)",
                               "PFDoA (C12)","PFTrDA (C13)")
                           ))

# for loop to make calculations for each of the six models in the model_list object
for(i in 1:length(model_list)){
  # Generate output from emmeans
  tl_estimates <- emmeans(model_list[[i]],
                          specs = ~ Trophic_Level,
                          type='response',tran = "log10") |>
    as_tibble()
  
  # Calculate percentages and store them in a new column labeled as wt_abundance
  tl_estimates <- tl_estimates %>% select(Trophic_Level,response) %>% 
    mutate(tl_abundance = 100*round(response/sum(response),digits=5))
  
  # Fill in the correct column of the matrix with the calculated percentages
  Table_S6_percent[1:7,i] <- tl_estimates$tl_abundance[1:7]
}

# Convert to tibble format
Table_S6_percent <- as_tibble(Table_S6_percent)

# Label rows with the correct tissue type and create a new column called group
# that can be used by the function gt() to label our final table
Table_S6_percent$rowname <- c("Primary Producer","Primary Consumer",
                              "Secondary Consumer",
                              "Tertiary Consumer","Quaternary Consumer",
                              "Piscivorous/Insectivorous Bird",
                              "Apex Predator")
Table_S6_percent$group <- "Percent Composition"

# Produce the finalized section of Table S6 using the function gt() and save
# this output as a .rtf file with the function gtsave()
Table_S6_final<-gt(Table_S6_percent,
                   rowname_col = "rowname",
                   groupname_col = "group",
                   row_group_as_column = F) |> 
  tab_stubhead(label = "Trophic Level") |>
  cols_align(
    align = "center",
    columns = everything()
  ) |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) |>
  opt_table_font(
    font = "Arial",
    weight = 350,
    size = 13) |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_row_groups()
  ) %>%
  fmt_markdown(columns = everything())

Table_S6_final
gtsave(Table_S6_final,"Tables and Figures/Table_S6.rtf")
# -----------------------------------------------------------------------------

######## Figures ##############################################################
######## Conceptual Figure ####################################################
# Generate spatial boundaries of the figure
max_lat <- max(final_imputed_data$Latitude)
min_lat <- min(final_imputed_data$Latitude)
max_lon <- max(final_imputed_data$Longitude)-5
min_lon <- min(final_imputed_data$Longitude)
# Store boundaries in a single extent object
geographic_extent <- ext(x = c(min_lon, max_lon, min_lat, max_lat))

# Crop the Great_Lakes_watershed object to this spatial extent
cropped_Great_Lakes_watershed <- st_crop(
  st_union(Great_Lakes_watershed,by_feature=F),
  geographic_extent)

# Remove internal borders from the sf object to create one seamless graph
intersection_Great_Lakes_watershed <- st_union(st_crop(
  st_intersection(Great_Lakes_region, 
                  st_union(Great_Lakes_watershed,by_feature=F)),
  geographic_extent*2),by_feature=F)

# Plot results
Conceptual_Figure <- ggplot() +
  geom_sf_pattern(data=cropped_Great_Lakes_watershed,
                  pattern='gradient',
                  pattern_fill="blue",pattern_fill2="red",
                  pattern_orientation = 'horizontal') +
  geom_sf(data=Great_Lakes_region,
          fill="white",color="white") +
  geom_sf_pattern(data=intersection_Great_Lakes_watershed,
                  pattern='gradient',
                  pattern_fill="lightblue",pattern_fill2="salmon",
                  pattern_orientation = 'horizontal') +
  geom_point(data=PFOS, 
             aes(x=Longitude,y=Latitude),
             shape=21,size=1.5,fill = "black")+
  coord_sf(xlim = c(min_lon,max_lon),ylim = c(min_lat-0.5,max_lat+0.5)) +
  annotation_north_arrow(location = "tl",which_north = "true", 
                         pad_x = unit(0.2, "cm"), pad_y = unit(0.2, "cm"),
                         style = north_arrow_nautical,width = unit(1.5, "cm"), 
                         height = unit(1.5, "cm")) +
  xlab("") +
  ylab("") +
  guides(fill= "none") +
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())
  
Conceptual_Figure
ggsave("Tables and Figures/Conceptual_Figure.png", plot = Conceptual_Figure,
       width = 12, height = 7, 
       units = "in", dpi = 300)

# -----------------------------------------------------------------------------
######## Figure 1 #############################################################
# Code to generate a sample distribution plot across the entire extent of the Great
# Lakes watersheds for the samples included in this meta-analysis (n = 2,489)
Figure_1 <- 
  ggplot() +
  geom_sf(data=Great_Lakes_region) +
  geom_sf(data=Great_Lakes_watershed,fill="cadetblue1",alpha=0.4) +
  coord_sf(xlim = c(-93.0,-70.0),ylim = c(40.5,51)) +
  annotation_north_arrow(location = "tl",which_north = "true", 
                         pad_x = unit(0.2, "cm"), pad_y = unit(0.2, "cm"),
                         style = north_arrow_nautical,width = unit(1.5, "cm"), 
                         height = unit(1.5, "cm")) +
  annotation_scale() +
  geom_point(data=PFOS, aes(x=Longitude,y=Latitude,
                            fill=factor(Waterbody,levels = WB_level_order)),
             shape = 21,size=2)+
  scale_fill_manual(values = c("#4575B4","#91BFDB","#FEE090","#FC8D59",
                               "#D73027")) +
  guides(fill= guide_legend(title = "Waterbody")) +
  theme_classic(base_size = 14) +
  theme(
    axis.title.x = element_text(size=14, face="bold", colour = "black"),    
    axis.title.y = element_text(size=14, face="bold", colour = "black"),
    legend.title = element_text(size=14, face="bold", colour = "black"),
    legend.text = element_text(size=12, colour = "black"),
    legend.position = c(.93,.888)
  )

Figure_1
ggsave("Tables and Figures/Figure_1S.png", plot = Figure_1, 
       width = 13, height = 9, units = "in",
      dpi = 300)
# -----------------------------------------------------------------------------
######## Figure 2 #############################################################
# Code to generate a unified temporal trends plot
# Create a vector of all years with data in the PFOS data frame over the entire
# 42 year sampling period
Sampling_Years<-sort(unique(PFOS$Sampling.Year))

# Now, for each lake, we 1) subset the PFOS data frame to isolate the data points
# for that lake, 2) isolate the vector of sampling years (starting point set as 
# the first year of sampling in that lake for which there is sufficient data) that
# we'll 3) feed in as an input to emmeans() to calculate model estimated
# concentrations for each of those years, 4) create a new column in the data frame of 
# emmeans() concentration estimates that contains labels for those estimates of the
# Lake watershed to which they belong, and 5) fit a spline to those concentration 
# estimates to help visualize the trend through time (this spline is saved as a 
# data frame object in R)
## Lake Erie
LE_PFOS<-subset(PFOS,Waterbody=="Lake Erie")
LE_Sampling_Years<-subset(Sampling_Years,
                          Sampling_Years>=min(LE_PFOS$Sampling.Year))
PFOS_LE_SY_means <- emmeans(full_PFOS_gam,
                            specs = ~ Sampling.Year,
                            type='response',
                            at = list(Sampling.Year = LE_Sampling_Years,
                                      Waterbody = "Lake Erie"),
                            tran="log10")  |>
  as_tibble()
PFOS_LE_SY_means$Waterbody<-"Lake Erie"

# Convert to a data frame for easy plotting in ggplot()
LE_spline <- as.data.frame(spline(PFOS_LE_SY_means$Sampling.Year, 
                                  PFOS_LE_SY_means$response))

## Lake Ontario
LO_PFOS<-subset(PFOS,Waterbody=="Lake Ontario")
LO_Sampling_Years<-subset(Sampling_Years,Sampling_Years>=1990)
PFOS_LO_SY_means <- emmeans(full_PFOS_gam,
                            specs = ~ Sampling.Year,
                            type='response',
                            at = list(Sampling.Year = LO_Sampling_Years,
                                      Waterbody = "Lake Ontario"),
                            tran="log10")  |>
  as_tibble()
PFOS_LO_SY_means$Waterbody<-"Lake Ontario"

LO_spline <- as.data.frame(spline(PFOS_LO_SY_means$Sampling.Year, 
                                  PFOS_LO_SY_means$response))

## Lake Superior
LS_PFOS<-subset(PFOS,Waterbody=="Lake Superior")
LS_Sampling_Years<-subset(Sampling_Years,
                          Sampling_Years>=min(LS_PFOS$Sampling.Year))
PFOS_LS_SY_means <- emmeans(full_PFOS_gam,
                            specs = ~ Sampling.Year,
                            type='response',
                            at = list(Sampling.Year = LS_Sampling_Years,
                                      Waterbody = "Lake Superior"),
                            tran="log10")  |>
  as_tibble()
PFOS_LS_SY_means$Waterbody<-"Lake Superior"

LS_spline <- as.data.frame(spline(PFOS_LS_SY_means$Sampling.Year, 
                                  PFOS_LS_SY_means$response))

## Lake Michigan
LM_PFOS<-subset(PFOS,Waterbody=="Lake Michigan")
LM_Sampling_Years<-subset(Sampling_Years,
                          Sampling_Years>=min(LM_PFOS$Sampling.Year))
PFOS_LM_SY_means <- emmeans(full_PFOS_gam,
                            specs = ~ Sampling.Year,
                            type='response',
                            at = list(Sampling.Year = LM_Sampling_Years,
                                      Waterbody = "Lake Michigan"),
                            tran="log10")  |>
  as_tibble()
PFOS_LM_SY_means$Waterbody<-"Lake Michigan"

LM_spline <- as.data.frame(spline(PFOS_LM_SY_means$Sampling.Year, 
                                  PFOS_LM_SY_means$response))

## Lake Huron
LH_PFOS<-subset(PFOS,Waterbody=="Lake Huron")
LH_Sampling_Years<-subset(Sampling_Years,
                          Sampling_Years>=min(LH_PFOS$Sampling.Year))
PFOS_LH_SY_means <- emmeans(full_PFOS_gam,
                            specs = ~ Sampling.Year,
                            type='response',
                            at = list(Sampling.Year = LH_Sampling_Years,
                                      Waterbody = "Lake Huron"),
                            tran="log10")  |>
  as_tibble()
PFOS_LH_SY_means$Waterbody<-"Lake Huron"

LH_spline <- as.data.frame(spline(PFOS_LH_SY_means$Sampling.Year, 
                                  PFOS_LH_SY_means$response))

# Bind all five concentration estimate data frames together into a single object
# that can be plotted by ggplot with the group argument set to the Waterbody
# variable
PFOS_SY_means <- rbind(PFOS_LS_SY_means,PFOS_LM_SY_means,PFOS_LH_SY_means,
                       PFOS_LE_SY_means,PFOS_LO_SY_means)

# Plot results
Figure_2 <- ggplot(PFOS_SY_means, aes(x = Sampling.Year, y=(response))) +
  annotate("rect", xmin=2000, xmax=2002, ymin=0,
           ymax=Inf,
           alpha=0.5,fill="orange") +
  geom_jitter(data = PFOS,
              aes(x = Sampling.Year,
                  y = (PFOS)),
              width = .2,
              alpha = .2,color="black") +
  theme_classic(base_size = 14) +
  ylab("PFOS concentration (ng/g w.w.)") +
  xlab("Sampling Year (1979-2021)") +
  coord_cartesian(ylim = c(0,310)) +
  geom_line(data = LM_spline, aes(x = x, y = y),color="#91BFDB",linewidth=1) +
  geom_line(data = LH_spline, aes(x = x, y = y),color="#FEE090",linewidth=1) +
  geom_line(data = LS_spline, aes(x = x, y = y),color="#4575B4",linewidth=1) +
  geom_line(data = LE_spline, aes(x = x, y = y),color="#FC8D59",linewidth=1) +
  geom_line(data = LO_spline, aes(x = x, y = y),color="#D73027",linewidth=1) +
  geom_pointrange(aes(ymin = 0, ymax = 0,
                      group = factor(Waterbody,levels = WB_level_order),
                      fill=factor(Waterbody,levels = WB_level_order)),
                  shape=21,color="black",size=0.5) +
  scale_fill_manual(values = c("#4575B4","#91BFDB","#FEE090","#FC8D59",
                               "#D73027")) +
  guides(fill = guide_legend(title = "Waterbody")) +
  theme(
    axis.title.x = element_text(size=14, face="bold", colour = "black"),    
    axis.title.y = element_text(size=14, face="bold", colour = "black"),
    legend.title = element_text(size=14, face="bold", colour = "black"),
    legend.text = element_text(size=12, colour = "black"),
    legend.position = c(.13,.8)
    )

Figure_2
ggsave("Tables and Figures/Figure_2.png", plot = Figure_2, 
       width = 13, height = 8, 
       units = "in", dpi = 300)
# -----------------------------------------------------------------------------
######## Figure 3 #############################################################
# Code calculating estimated marginal means in each Great Lakes watershed for 
# each of the six PFAS included in this meta-analysis. We 1) use the emmeans()
# function to obtain model-estimated concentrations, converting the output to a
# tibble object, and 2) create a new column in the tibble that labels these
# concentration estimates with the given contaminant
## PFOS
PFOS_WB_means <- emmeans(full_PFOS_gam,
                         specs = ~ Waterbody,
                         type='response',tran = "log10")  |>
  as_tibble()
PFOS_WB_means$contaminant<-"PFOS"

## PFNA
PFNA_WB_means <- emmeans(full_PFNA_gam,
                         specs = ~ Waterbody,
                         type='response',tran = "log10")  |>
  as_tibble()
PFNA_WB_means$contaminant<-"PFNA"

## PFDA
PFDA_WB_means <- emmeans(full_PFDA_gam,
                         specs = ~ Waterbody,
                         type='response',tran = "log10")  |>
  as_tibble()
PFDA_WB_means$contaminant<-"PFDA"

## PFUnA
PFUnA_WB_means <- emmeans(full_PFUnA_gam,
                          specs = ~ Waterbody,
                          type='response',tran = "log10")  |>
  as_tibble()
PFUnA_WB_means$contaminant<-"PFUnA"

## PFDoA
PFDoA_WB_means <- emmeans(full_PFDoA_gam,
                          specs = ~ Waterbody,
                          type='response',tran = "log10")  |>
  as_tibble()
PFDoA_WB_means$contaminant<-"PFDoA"

## PFTrDA
PFTrDA_WB_means <- emmeans(full_PFTrDA_gam,
                           specs = ~ Waterbody,
                           type='response',tran = "log10")  |>
  as_tibble()
PFTrDA_WB_means$contaminant<-"PFTrDA"

# Bind concentration estimates from each model into one object, and calculate
# the concentration sum of all six PFAS for each watershed
PFAS_WB_means<-rbind(PFOS_WB_means,PFNA_WB_means,PFDA_WB_means,PFUnA_WB_means,
                     PFDoA_WB_means,PFTrDA_WB_means)
PFAS_WB_means$sum<-ave(PFAS_WB_means$response, PFAS_WB_means$Waterbody, FUN=sum)

# Order the levels of the Waterbody and contaminant variables in this unified tibble
PFAS_WB_means <- PFAS_WB_means %>%
  arrange(factor(Waterbody, levels = rev(WB_level_order)))
PFAS_WB_means$contaminant <- factor(PFAS_WB_means$contaminant, 
                                    levels=c('PFOS', 'PFNA', 'PFDA', 
                                             'PFUnA','PFDoA','PFTrDA'))

# Plot results as a bar chart
Figure_3 <- 
  ggplot(PFAS_WB_means, aes(y = Waterbody, 
                            x = response)) +
    geom_col(aes(fill = contaminant)) +
    scale_fill_manual(values = c("#4575B4","#91BFDB","#E0F3F8",
                                 "#FEE090","#FC8D59",
                                 "#D73027")) +
    theme_classic(base_size = 14) +
    ylab(NULL) +
    xlab("Concentration (ng/g w.w.)") +
    geom_label(aes(x = sum,y = Waterbody,
                   label = paste0(round(sum,2)," ng/g w.w."),
                   group = factor(contaminant)),
               hjust = -0.25,fontface='bold',size=4) +
    coord_cartesian(xlim = c(0,225)) +
    guides(fill = guide_legend(title = "PFAS Compound")) +
    theme(legend.title=element_text(size=14,face = "bold"),
          axis.title.x = element_text(size=14, face="bold", colour = "black"),
          axis.text.y = element_text(size=12, colour = "black"),
          legend.text = element_text(size=12, colour = "black"),
          legend.position = c(.87,.3)
          )

Figure_3
ggsave("Tables and Figures/Figure_3.png", plot = Figure_3, width = 13, 
       height = 7,  units = "in", dpi = 300)
# -----------------------------------------------------------------------------
######## Figure 4 #############################################################
# Code calculating model-estimated concentrations of PFOS across the seven trophic
# levels of the Great Lakes food web
PFOS_TL_means <- 
  emmeans(full_PFOS_gam,
        specs = ~ Trophic_Level,
        type='response',tran = "log10") |>
  as_tibble()

# Reorder the factor levels of the Trophic_Level variable from the output of
# the emmeans() function
PFOS_TL_means$Trophic_Level<-factor(PFOS_TL_means$Trophic_Level,
                                    levels = c("Primary Producer",
                                               "Primary Consumer",
                                               "Secondary Consumer",
                                               "Piscivorous/Insectivorous Bird",
                                               "Tertiary Consumer",
                                               "Quaternary Consumer",
                                               "Apex Predator"))

# Plot results as a scatter plot with model estimates and 95% CI
Figure_4 <- 
  ggplot(PFOS_TL_means, aes(x = Trophic_Level,y=response))+
  geom_jitter(data = PFOS,
              aes(x = Trophic_Level,
                  y = PFOS),
              width = .2,
              alpha = .2,
              size = 2)+
  geom_pointrange(aes(ymin = lower.CL, ymax = upper.CL),color="red",
                  alpha=1,show.legend = FALSE,
                  size = 0.8)+
  theme_classic(base_size = 14)+
  scale_y_continuous(transform = "log10",
                     labels = format_format(scientific=FALSE),
                     limits = c(0.0005,5000)) +
  scale_x_discrete(labels= c("PP","PC","SC","PIB","TC","QC","AP")) +
  ylab("PFOS concentration (ng/g w.w.)") +
  xlab("") +
  theme(axis.title.y = element_text(size=14, face="bold", colour = "black"))

Figure_4
ggsave("Tables and Figures/Figure_4.png", plot = Figure_4, 
       width = 10, height = 7, 
       units = "in", dpi = 300)
# -----------------------------------------------------------------------------
######## Figure S3 #############################################################
# Code to show sampling efforts for the various taxonomic "classes" through time
# within the five watersheds of the Great Lakes
# For simplicity, a number of taxonomic classes are merged into one group called
# the "Invertebrates"
for (i in 1:nrow(final_imputed_data)){
  if(final_imputed_data$Class[i]=="Bivalvia"
     || final_imputed_data$Class[i]=="Annelida"
     || final_imputed_data$Class[i]=="Amphipoda"
     || final_imputed_data$Class[i]=="Gastropoda"
     || final_imputed_data$Class[i]=="Shrimp, water fleas, and allies"
     || final_imputed_data$Class[i]=="Astacoidea"
     || final_imputed_data$Class[i]=="Zooplankton"
     || final_imputed_data$Class[i]=="Insecta"){
    final_imputed_data$Class[i]<-"Invertebrates"
  }
}

# Create a new summary data frame with a variable called sample_count, that counts
# the number of samples taken for each taxonomic class in each watershed each year
waterbody_class_sampling_df <- final_imputed_data %>% 
  group_by(Sampling.Year,Class,Waterbody) %>% 
  summarise(sample_count = sum(n_samples))

# Define order of taxonomic classes
Class_order<-c("Algae","Plantae (Magnoliopsida)","Invertebrates","Amphibia",
               "Pisces","Reptilia","Mammalia","Aves")

# Plot results in column chart format
Figure_S3<-
  ggplot(waterbody_class_sampling_df)+ 
  geom_bar(aes(x = Sampling.Year, y = sample_count, 
               fill = factor(Class,
                             levels = Class_order)), 
           stat = "identity") +
  scale_fill_brewer(palette = "RdBu") +
  geom_bar(aes(x = Sampling.Year, y = sample_count),
           colour = "black", linewidth = 0.3,
           stat = "summary", fun = sum,fill="transparent") +
  ylab("Sample Count") +
  xlab("Sampling Year (1979-2021)") +
  guides(fill= guide_legend(title = "Taxonomic Class")) +
  theme_bw(base_size = 14) +
  theme(legend.position = c(0.85, 0.25),
        panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank(),   # Remove minor gridlines
        legend.title=element_text(size=14,face = "bold"),
        axis.title.x = element_text(size=14, face="bold", colour = "black"),
        axis.title.y = element_text(size=14, face="bold", colour = "black"),
        legend.text = element_text(size=12, colour = "black")) +
  facet_wrap(~factor(Waterbody,levels = WB_level_order))

Figure_S3
ggsave("Tables and Figures/Figure_S3.png", plot = Figure_S3, 
       width = 10, height = 7, 
       units = "in", dpi = 300)
# -----------------------------------------------------------------------------
######## Figure S4 ############################################################
# Code to show sampling efforts for the various trophic levels through time
# within the five watersheds of the Great Lakes
# Create a new summary data frame with a variable called sample_count, that counts
# the number of samples taken for each trophic level in each watershed each year
waterbody_tl_sampling_df <- final_imputed_data %>% 
  group_by(Sampling.Year,Trophic_Level,Waterbody) %>% 
  summarise(sample_count = sum(n_samples))

# Define order of trophic levels
Trophic_Level_order<-c("Primary Producer","Primary Consumer",
                       "Secondary Consumer","Tertiary Consumer",
                       "Quaternary Consumer","Piscivorous/Insectivorous Bird",
                       "Apex Predator")

# Plot results in column chart format
Figure_S4<-
  ggplot(waterbody_tl_sampling_df)+ 
  geom_bar(aes(x = Sampling.Year, y = sample_count, 
               fill = factor(Trophic_Level,
                             levels = Trophic_Level_order)), 
           stat = "identity") +
  scale_fill_brewer(palette = "RdBu") +
  geom_bar(aes(x = Sampling.Year, y = sample_count),
           colour = "black", linewidth = 0.3,
           stat = "summary", fun = sum,fill="transparent") +
  ylab("Sample Count") +
  xlab("Sampling Year (1979-2021)") +
  guides(fill= guide_legend(title = "Trophic Level")) +
  theme_bw(base_size = 14) +
  theme(legend.position = c(0.85, 0.25),
        panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank(),   # Remove minor gridlines
        legend.title=element_text(size=14,face = "bold"),
        axis.title.x = element_text(size=14, face="bold", colour = "black"),
        axis.title.y = element_text(size=14, face="bold", colour = "black"),
        legend.text = element_text(size=12, colour = "black"),
        ) +
  facet_wrap(~factor(Waterbody,levels = WB_level_order))

Figure_S4
ggsave("Tables and Figures/Figure_S4.png", plot = Figure_S4, 
       width = 10, height = 7, 
       units = "in", dpi = 300)
# -----------------------------------------------------------------------------
######## Figure S5 ############################################################
# Code to generate a multi-panel plot showing temporal trends and distributions of
# PFOS concentration values for the five watersheds of the Great Lakes. Instead of
# generating separate spline and concentration estimates for each watershed and
# then merging them all together to be plotted with one ggplot() command, we 
# instead create five separate plots (one for each watershed of the sampling 
# distribution, spline and model-estimated yearly concentrations) and link those
# five plots together using the ggarrange() function
## Lake Superior Plot
LS_SY_plot<-ggplot(PFOS_LS_SY_means,aes(x = Sampling.Year, y=response))+
  annotate("rect", xmin=2000, xmax=2002, ymin=0,
           ymax=Inf,
           alpha=0.5,fill="orange") +
  geom_jitter(data = LS_PFOS,
              aes(x = Sampling.Year,y = PFOS),
              width = .2,
              alpha = .2,color="black")+
  geom_line(data = LS_spline, aes(x = x, y = y),color="#4575B4",linewidth=1)+
  geom_pointrange(aes(ymin = lower.CL, ymax = upper.CL,
                      group = factor(Waterbody,levels = WB_level_order),
                      fill=factor(Waterbody,levels = WB_level_order)),
                  shape = 21,color="black",size=0.5,
                  show.legend = T)+
  theme_classic(base_size = 14)+
  ylab("PFOS concentration (ng/g w.w.)") +
  xlab("Sampling Year (1991-2021)") +
  ggtitle("Lake Superior") +
  theme(
    axis.title.x = element_text(size=14, face="bold", colour = "black"),    
    axis.title.y = element_text(size=14, face="bold", colour = "black"),
    plot.title = element_text(size=14, face="bold", colour = "black",
                              hjust = 0.5)) +
  scale_fill_manual(values = c("#4575B4","#91BFDB","#FEE090","#FC8D59",
                               "#D73027")) +
  scale_y_continuous(transform = "log10",
                     labels = format_format(scientific=FALSE),
                     limits = c(0.01,73000)) +
  guides(fill = guide_legend(title = "Waterbody"))

## Lake Michigan Plot
LM_SY_plot<-ggplot(PFOS_LM_SY_means, aes(x = Sampling.Year,y=response))+
  annotate("rect", xmin=2000, xmax=2002, ymin=0,
           ymax=Inf,
           alpha=0.5,fill="orange") +
  geom_jitter(data = LM_PFOS,
              aes(x = Sampling.Year,
                  y = PFOS),
              width = .2,
              alpha = .2,color="black")+
  geom_line(data = LM_spline, aes(x = x, y = y),color="#91BFDB",linewidth=1)+
  geom_pointrange(aes(ymin = lower.CL, ymax = upper.CL,
                      group = factor(Waterbody,levels = WB_level_order),
                      fill=factor(Waterbody,levels = WB_level_order)),
                  shape = 21,color="black",size=0.5,
                  show.legend = T)+
  theme_classic(base_size = 14)+
  ylab("") +
  xlab("Sampling Year (1990-2021)")+
  ggtitle("Lake Michigan") +
  theme(
    axis.title.x = element_text(size=14, face="bold", colour = "black"),    
    axis.title.y = element_text(size=14, face="bold", colour = "black"),
    plot.title = element_text(size=14, face="bold", colour = "black",
                              hjust = 0.5)) +
  scale_fill_manual(values = c("#91BFDB","#4575B4","#FEE090","#FC8D59",
                               "#D73027")) +
  scale_y_continuous(transform = "log10",
                     labels = format_format(scientific=FALSE),
                     limits = c(0.01,73000)) +
  guides(fill = guide_legend(title = "Waterbody"))

# Lake Huron Plot
LH_SY_plot<-ggplot(PFOS_LH_SY_means,aes(x = Sampling.Year, y=response))+
  annotate("rect", xmin=2000, xmax=2002, ymin=0,
           ymax=Inf,
           alpha=0.5,fill="orange") +
  geom_jitter(data = LH_PFOS,
              aes(x = Sampling.Year,
                  y = PFOS),
              width = .2,
              alpha = .2,color="black")+
  geom_line(data = LH_spline, aes(x = x, y = y),color="#FEE090",linewidth=1)+
  geom_pointrange(aes(ymin = lower.CL, ymax = upper.CL,
                      group = factor(Waterbody,levels = WB_level_order),
                      fill=factor(Waterbody,levels = WB_level_order)),
                  shape = 21,color="black",size=0.5,
                  show.legend = T)+
  theme_classic(base_size = 14)+
  ylab("") +
  xlab("Sampling Year (1991-2021)")+
  ggtitle("Lake Huron") +
  theme(
    axis.title.x = element_text(size=14, face="bold", colour = "black"),    
    axis.title.y = element_text(size=14, face="bold", colour = "black"),
    plot.title = element_text(size=14, face="bold", colour = "black",
                              hjust = 0.5)) +
  scale_fill_manual(values = c("#FEE090","#91BFDB","#4575B4","#FC8D59",
                               "#D73027")) +
  scale_y_continuous(transform = "log10",
                     labels = format_format(scientific=FALSE),
                     limits = c(0.01,73000)) +
  guides(fill = guide_legend(title = "Waterbody"))

# Lake Erie Plot
LE_SY_plot<-ggplot(PFOS_LE_SY_means, aes(x = Sampling.Year, y=response)) +
  annotate("rect", xmin=2000, xmax=2002, ymin=0,
           ymax=Inf,
           alpha=0.5,fill="orange") +
  geom_jitter(data = LE_PFOS,
              aes(x = Sampling.Year,
                  y = PFOS),
              width = .2,
              alpha = .2,color="black") +
  geom_line(data = LE_spline, aes(x = x, y = y),color="#FC8D59",linewidth=1)+
  geom_pointrange(aes(ymin = lower.CL, ymax = upper.CL,
                      group = factor(Waterbody,levels = WB_level_order),
                      fill=factor(Waterbody,levels = WB_level_order)),
                  shape = 21,color="black",size=0.5,
                  show.legend = T)+
  theme_classic(base_size = 14)+
  ylab("PFOS concentration (ng/g w.w.)") +
  xlab("Sampling Year (1992-2021)")+
  ggtitle("Lake Erie") +
  theme(
    axis.title.x = element_text(size=14, face="bold", colour = "black"),    
    axis.title.y = element_text(size=14, face="bold", colour = "black"),
    plot.title = element_text(size=14, face="bold", colour = "black",
                              hjust = 0.5)) +
  scale_fill_manual(values = c("#FC8D59","#FEE090","#91BFDB","#4575B4",
                               "#D73027")) +
  scale_y_continuous(transform = "log10",
                     labels = format_format(scientific=FALSE),
                     limits = c(0.01,73000)) +
  guides(fill = guide_legend(title = "Waterbody"))

# Lake Ontario Plot
LO_SY_plot<-ggplot(PFOS_LO_SY_means, aes(x = Sampling.Year,y=response))+
  annotate("rect", xmin=2000, xmax=2002, ymin=0,
           ymax=Inf,
           alpha=0.5,fill="orange") +
  geom_jitter(data = LO_PFOS,
              aes(x = Sampling.Year,
                  y = PFOS),
              width = .2,
              alpha = .2,color="black")+
  geom_line(data = LO_spline, aes(x = x, y = y),color="#D73027",linewidth=1)+
  geom_pointrange(aes(ymin = lower.CL, ymax = upper.CL,
                      group = factor(Waterbody,levels = WB_level_order),
                      fill=factor(Waterbody,levels = WB_level_order)),
                  shape = 21,color="black",size=0.5,
                  show.legend = TRUE)+
  theme_classic(base_size = 14)+
  ggtitle("Lake Ontario") +
  theme(
    axis.title.x = element_text(size=14, face="bold", colour = "black"),    
    axis.title.y = element_text(size=14, face="bold", colour = "black"),
    plot.title = element_text(size=14, face="bold", colour = "black",
                              hjust = 0.5)) +
  ylab("") +
  xlab("Sampling Year (1979-2021)")+
  scale_fill_manual(values = c("#D73027","#FC8D59","#FEE090","#91BFDB",
                               "#4575B4")) +
  scale_y_continuous(transform = "log10",
                     labels = format_format(scientific=FALSE),
                     limits = c(0.01,73000)) +
  guides(fill = guide_legend(title = "Waterbody"))

Figure_S5 <- ggarrange(LS_SY_plot, LM_SY_plot,LH_SY_plot,LE_SY_plot,LO_SY_plot,
                       ncol = 3,nrow = 2,legend = "none")

Figure_S5
ggsave("Tables and Figures/Figure_S5.png", plot = Figure_S5, width = 14, 
       height = 10,  units = "in", dpi = 300)
# -----------------------------------------------------------------------------
######## Figure S6 ############################################################
# Code to plot model-estimated concentrations (95% CI) of PFOS across the five 
# watersheds of the Great Lakes, along with a scatter plot of samples in each watershed
Figure_S6 <-
  ggplot(PFOS_WB_means, aes(x = factor(Waterbody, level = WB_level_order),
                          y=response))+
  geom_jitter(data = PFOS,
              aes(x = factor(Waterbody, level = WB_level_order),
                  y = PFOS),
              width = .2,
              alpha = .2,
              size = 2)+
  geom_pointrange(aes(ymin = lower.CL, ymax = upper.CL),color="red",
                  alpha=1,show.legend = FALSE,
                  size = 0.8)+
  theme_classic(base_size = 14)+
  scale_y_continuous(transform = "log10",
                     labels = format_format(scientific=FALSE),
                     limits = c(0.0005,5000)) +
  ylab("PFOS concentration (ng/g w.w.)") +
  xlab("") +
  theme(axis.title.y = element_text(size=14, face="bold", colour = "black"))

Figure_S6
ggsave("Tables and Figures/Figure_S6.png", plot = Figure_S6, 
       width = 10, height = 7, 
       units = "in", dpi = 300)
# -----------------------------------------------------------------------------
######## Figure S7 ############################################################
# Code to plot model-estimated concentrations (95% CI) of the six PFAS across the five 
# watersheds of the Great Lakes, as a way to visualize the contrasting spatial patterns
# of these contaminants
Figure_S7 <- 
  ggplot(PFAS_WB_means, aes(x = factor(Waterbody,level = WB_level_order), 
                                       y = (response), 
                                       fill = contaminant)) +
  geom_pointrange(aes(ymin = lower.CL, ymax = upper.CL),
                  pch=21,alpha=1,position = position_dodge(width = 0.5),
                  color="black") +
  scale_fill_manual(values = c("#4575B4","#91BFDB","#E0F3F8","#FEE090",
                                "#FC8D59","#D73027")) +
  scale_y_continuous(transform = "log10",
                     labels = format_format(scientific=FALSE)) +
  xlab(NULL) +
  ylab("Concentration (ng/g w.w.)") +
  guides(fill = guide_legend(title = "PFAS Compound")) +
  theme_classic(base_size = 14) +
  theme(
    axis.title.x = element_text(size=14, face="bold", colour = "black"),    
    axis.title.y = element_text(size=14, face="bold", colour = "black"),
    legend.title = element_text(size=14, face="bold", colour = "black")
  )

Figure_S7
ggsave("Tables and Figures/Figure_S7.png", plot = Figure_S7, width = 10, 
       height = 7,  units = "in", dpi = 300)
# -----------------------------------------------------------------------------
# Delete excess variables
rm(WB_level_order,Trophic_Level_order,Sampling_Years,Class_order,
   min_lon,min_lat,max_lat,max_lon,geographic_extent,
   LS_Sampling_Years,LO_Sampling_Years,LM_Sampling_Years,LH_Sampling_Years,
   LE_Sampling_Years,
   i,wt_estimates,tl_estimates,tissue_estimates,Table_2_percent,Table_S5_percent,
   Table_S6_percent,waterbody_tl_sampling_df,waterbody_class_sampling_df,
   PFOS_WB_means,PFNA_WB_means,PFDA_WB_means,PFUnA_WB_means,
   PFDoA_WB_means,PFTrDA_WB_means,PFAS_WB_means,
   PFOS_TL_means,PFOS_LE_SY_means,PFOS_LH_SY_means,PFOS_LM_SY_means,PFOS_LO_SY_means,
   PFOS_LS_SY_means,PFOS_SY_means,model_list,
   LS_SY_plot,LE_SY_plot,LO_SY_plot,LM_SY_plot,LH_SY_plot,
   LS_spline,LE_spline,LO_spline,LM_spline,LH_spline,
   LS_PFOS,LE_PFOS,LO_PFOS,LM_PFOS,LH_PFOS,
   intersection_Great_Lakes_watershed,cropped_Great_Lakes_watershed,
   Great_Lakes_watershed, Great_Lakes_region)



