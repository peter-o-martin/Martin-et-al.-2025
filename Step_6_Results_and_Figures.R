# Step_6_Results_and_Figures.R
# Contains all code used to extract results from models and construct figures

# Written by Peter Martin
# Created December 13, 2024
# Finalized

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
final_imputed_data <- read.csv("Step_5_final_imputed_data.csv",header = TRUE)

PFOS<-final_imputed_data[is.na(final_imputed_data$PFOS)==FALSE,]
set.seed(166)
PFOS_in_train<-createDataPartition(PFOS$PFOS, p = 4/5, list = FALSE)
PFOS_train_set<-PFOS[PFOS_in_train,]
PFOS_test_set<-PFOS[-PFOS_in_train,]

## Loading in the models from Step 5
load("full_PFOS_gam.Rdata") 
load("full_PFNA_gam.Rdata")
load("full_PFDA_gam.Rdata")
load("full_PFUnA_gam.Rdata")
load("full_PFDoA_gam.Rdata")
load("full_PFTrDA_gam.Rdata")


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


######## Table 2 ##############################################################
# Code for generating model-estimated mean concentrations and pairwise 
# contrasts
Table_2_me <- emmeans(full_PFOS_gam, 
           specs = pairwise ~ Revised_Tissue,
           type = "response",tran = "log10",adjust="tukey")
Table_2_me

# Code for calculating the relative distributions (%) of each PFAS across
# four different tissue groups (Eggs, Blood, Liver, and Combined Tissue and
# Blood [i.e., all levels of the variable Revised_Tissue except Eggs])

# Create the matrix that will be filled with the calculated percentage values
Table_2_percent <- matrix(nrow = 4, ncol = 6,
                  dimnames = list(
                    c("Eggs","Blood","Liver","Combined Tissue and Blood"),
                    c("PFOS (C8)","PFNA (C9)","PFDA (C10)","PFUnA (C11)",
                     "PFDoA (C12)","PFTrDA (C13)")
                    ))

# for loop to make calculations for each of the six models in the model_list 
# variable
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
}

# Convert to tibble format
Table_2_percent <- as_tibble(Table_2_percent)

# Label rows with the correct tissue type and create a new column called group
# that can be used by the function gt() to label our final table
Table_2_percent$rowname <- c("Eggs","Blood","Liver","Combined Tissue and Blood")
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

######## Figures ##############################################################
######## Conceptual Figure ####################################################
max_lat <- max(final_imputed_data$Latitude)
min_lat <- min(final_imputed_data$Latitude)
max_lon <- max(final_imputed_data$Longitude)-5
min_lon <- min(final_imputed_data$Longitude)
# Store boundaries in a single extent object
geographic_extent <- ext(x = c(min_lon, max_lon, min_lat, max_lat))

cropped_Great_Lakes_watershed <- st_crop(
  st_union(Great_Lakes_watershed,by_feature=F),
  geographic_extent)

intersection_Great_Lakes_watershed <- st_union(st_crop(
  st_intersection(Great_Lakes_region, 
                  st_union(Great_Lakes_watershed,by_feature=F)),
  geographic_extent*2),by_feature=F)

PFOS$Trophic_Level <- factor(PFOS$Trophic_Level,
                             levels = c("Primary Producer","Primary Consumer",
                                        "Secondary Consumer","Tertiary Consumer",
                                        "Quaternary Consumer",
                                        "Piscivorous/Insectivorous Bird",
                                        "Apex Predator")) 
levels(PFOS$Trophic_Level)
levels(PFOS$Trophic_Level) <- 0:6

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
# https://bookdown.org/nicohahn/making_maps_with_r5/docs/ggplot2.html
Figure_1 <- ggplot() +
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
ggsave("Tables and Figures/Figure_1.png", plot = Figure_1, 
       width = 13, height = 9, units = "in",
      dpi = 300)
# -----------------------------------------------------------------------------
######## Figure 2 #############################################################
Sampling_Years<-sort(unique(PFOS$Sampling.Year))

LE_PFOS<-subset(PFOS,Waterbody=="Lake Erie")
LE_Sampling_Years<-subset(Sampling_Years,Sampling_Years>=min(LE_PFOS$Sampling.Year))
PFOS_LE_SY_means <- emmeans(full_PFOS_gam,
                            specs = ~ Sampling.Year,
                            type='response',
                            at = list(Sampling.Year = LE_Sampling_Years,
                                      Waterbody = "Lake Erie"),
                            tran="log10")  |>
  as_tibble()
PFOS_LE_SY_means$Waterbody<-"Lake Erie"

LE_spline <- as.data.frame(spline(PFOS_LE_SY_means$Sampling.Year, 
                                  PFOS_LE_SY_means$response))


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

LS_PFOS<-subset(PFOS,Waterbody=="Lake Superior")
LS_Sampling_Years<-subset(Sampling_Years,Sampling_Years>=min(LS_PFOS$Sampling.Year))
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

LM_PFOS<-subset(PFOS,Waterbody=="Lake Michigan")
LM_Sampling_Years<-subset(Sampling_Years,Sampling_Years>=min(LM_PFOS$Sampling.Year))
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

LH_PFOS<-subset(PFOS,Waterbody=="Lake Huron")
LH_Sampling_Years<-subset(Sampling_Years,Sampling_Years>=min(LH_PFOS$Sampling.Year))
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

PFOS_SY_means <- rbind(PFOS_LS_SY_means,PFOS_LM_SY_means,PFOS_LH_SY_means,
                       PFOS_LE_SY_means,PFOS_LO_SY_means)

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

# sms<-smooths(full_PFOS_gam)
# ds <- full_PFOS_gam |> 
#   data_slice(Sampling.Year=evenly(Sampling.Year),Waterbody = evenly(Waterbody))
# fv <- full_PFOS_gam |>
#   fitted_values(data = ds, exclude = sms[c(1:5,11,12)])
# 
# fv |>
#   ggplot(aes(x=Sampling.Year, y=10^.fitted, group = Waterbody)) +
#   geom_ribbon(aes(ymax=10^.upper_ci,ymin=10^.lower_ci,fill=Waterbody),
#               alpha = 0.2) +
#   geom_line(aes(colour = Waterbody),linewidth=1) +
#   coord_cartesian(ylim = c(0,300)) +
#   theme_bw(base_size = 14) +
#   theme(
#     axis.title.x = element_text(size=14, face="bold", colour = "black"),    
#     axis.title.y = element_text(size=14, face="bold", colour = "black"),
#     legend.title = element_text(size=14, face="bold", colour = "black")
#   )
# 
# full_PFOS_gam |>
#   conditional_values(
#     condition = c("Sampling.Year","Waterbody"),
#     exclude = sms[c(1:5,11,12)]
#   ) |>
#   draw()

# -----------------------------------------------------------------------------
######## Figure 3 #############################################################
# Loading emmeans for Watershed differences #
# PFOS concentrations in different Waterbodies
PFOS_WB_means <- emmeans(full_PFOS_gam,
                         specs = ~ Waterbody,
                         type='response',tran = "log10")  |>
  as_tibble()
PFOS_WB_means$contaminant<-"PFOS"
# PFNA concentrations in different Waterbodies
PFNA_WB_means <- emmeans(full_PFNA_gam,
                         specs = ~ Waterbody,
                         type='response',tran = "log10")  |>
  as_tibble()
PFNA_WB_means$contaminant<-"PFNA"
# PFDA concentrations in different Waterbodies
PFDA_WB_means <- emmeans(full_PFDA_gam,
                         specs = ~ Waterbody,
                         type='response',tran = "log10")  |>
  as_tibble()
PFDA_WB_means$contaminant<-"PFDA"
# PFUnA concentrations in different Waterbodies
PFUnA_WB_means <- emmeans(full_PFUnA_gam,
                          specs = ~ Waterbody,
                          type='response',tran = "log10")  |>
  as_tibble()
PFUnA_WB_means$contaminant<-"PFUnA"
# PFDoA concentrations in different Waterbodies
PFDoA_WB_means <- emmeans(full_PFDoA_gam,
                          specs = ~ Waterbody,
                          type='response',tran = "log10")  |>
  as_tibble()
PFDoA_WB_means$contaminant<-"PFDoA"
# PFTrDA concentrations in different Waterbodies
PFTrDA_WB_means <- emmeans(full_PFTrDA_gam,
                           specs = ~ Waterbody,
                           type='response',tran = "log10")  |>
  as_tibble()
PFTrDA_WB_means$contaminant<-"PFTrDA"


PFAS_WB_means<-rbind(PFOS_WB_means,PFNA_WB_means,PFDA_WB_means,PFUnA_WB_means,
                     PFDoA_WB_means,PFTrDA_WB_means)
PFAS_WB_means$sum<-ave(PFAS_WB_means$response, PFAS_WB_means$Waterbody, FUN=sum)

PFAS_WB_means$percent<-(PFAS_WB_means$response/PFAS_WB_means$sum)

PFOS_WB_means <- PFOS_WB_means %>%
  arrange(factor(Waterbody, levels = rev(WB_level_order)))
PFAS_WB_means <- PFAS_WB_means %>%
  arrange(factor(Waterbody, levels = rev(WB_level_order)))
PFAS_WB_means$contaminant <- factor(PFAS_WB_means$contaminant, 
                                    levels=c('PFOS', 'PFNA', 'PFDA', 
                                             'PFUnA','PFDoA','PFTrDA'))

# Figure_3 <- ggarrange(
#   ggplot(PFAS_WB_means, aes(y = factor(Waterbody, level = WB_level_order), 
#                                     x = response)) +
#     geom_col(aes(fill = contaminant)) +
#     scale_fill_manual(values = c("#4575B4","#91BFDB","#E0F3F8",
#                                  "#FEE090","#FC8D59",
#                                  "#D73027")) +
#     theme_bw(base_size = 14) +
#     ylab(NULL) +
#     xlab("Concentration (ng/g w.w.)") +
#     geom_label(aes(x = sum,y = factor(Waterbody, level = WB_level_order),
#                    label = paste0(round(sum,2)," ng/g w.w."),
#                    group = factor(contaminant)),
#                hjust = -0.25,fontface='bold',size=4) +
#     coord_cartesian(xlim = c(0,225)) +
#     guides(fill = guide_legend(title = "PFAS Compound")) +
#     theme(legend.title=element_text(size=14,face = "bold"),
#           axis.title.x = element_text(size=14, face="bold", colour = "black")),
#   
#   ggplot(PFAS_WB_means, aes(y = factor(Waterbody, level = WB_level_order), 
#                             x = percent*100, fill = contaminant)) +
#     geom_col() +
#     scale_fill_manual(values = c("#4575B4","#91BFDB","#E0F3F8",
#                                  "#FEE090","#FC8D59",
#                                  "#D73027")) +
#     geom_label_repel(aes(label=paste0(sprintf("%1.1f", percent*100),"%")),
#                      position=position_stack(vjust=0.5),
#                      force=0.01,force_pull = 190,
#                      max.overlaps = 8,direction = "both",show.legend = FALSE,
#                      fontface='bold',size=4) +
#     ylab(NULL) +
#     xlab("Contaminant Composition (%)") +
#     theme_bw(base_size = 14) +
#     guides(fill = guide_legend(title = "PFAS Compound")) +
#     theme(legend.title=element_text(size=14,face = "bold"),
#           axis.title.x = element_text(size=14, face="bold", colour = "black")),
#   
#   ncol = 1,
#   nrow = 2,
#   legend="right",
#   common.legend = T
#   )

levels(PFAS_WB_means$Waterbody)

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
          legend.position = c(.87,.3),
    )

# Figure_3 <-
#   ggplot(PFAS_WB_means, aes(y = Waterbody, 
#                             x = response)) +
#     geom_col(aes(fill = contaminant)) +
#     scale_fill_manual(values = c("#4575B4","#91BFDB","#E0F3F8",
#                                  "#FEE090","#FC8D59",
#                                  "#D73027")) +
#     theme_bw(base_size = 14) +
#     ylab(NULL) +
#     xlab("Concentration (ng/g w.w.)") +
#     coord_cartesian(xlim = c(0,210)) +
#     guides(fill = guide_legend(title = "PFAS Compound")) +
#     theme(legend.title=element_text(size=14,face = "bold"),
#           axis.title.x = element_text(size=14, face="bold", colour = "black")) +
#   theme(legend.position = c(.95, .1),
#         legend.justification = c("right", "bottom"),
#         legend.box.just = "right",
#         legend.margin = margin(6, 6, 6, 6)
#           ) +
#   theme(axis.text.x= element_text(size = 14),
#         axis.text.y= element_text(size = 14),
#         legend.text = element_text(size = 14))


Figure_3
ggsave("Tables and Figures/Figure_3.png", plot = Figure_3, width = 13, 
       height = 7,  units = "in", dpi = 300)

# -----------------------------------------------------------------------------
######## Figure 4 #############################################################

PFOS_TL_means <- 
  emmeans(full_PFOS_gam,
        specs = ~ Trophic_Level,
        type='response',tran = "log10") |>
  as_tibble()

PFOS_TL_means$Trophic_Level<-factor(PFOS_TL_means$Trophic_Level,levels = c("Primary Producer","Primary Consumer","Secondary Consumer","Piscivorous/Insectivorous Bird","Tertiary Consumer","Quaternary Consumer","Apex Predator"))

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
  theme(axis.title.y = element_text(size=14, face="bold", colour = "black")
  )

Figure_4
ggsave("Tables and Figures/Figure_4.png", plot = Figure_4, 
       width = 10, height = 7, 
       units = "in", dpi = 300)
# -----------------------------------------------------------------------------
######## Figure S3 #############################################################
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

waterbody_class_sampling_df <- final_imputed_data %>% 
  group_by(Sampling.Year,Class,Waterbody) %>% 
  summarise(sample_count = sum(n_samples))

Class_order<-c("Algae","Plantae (Magnoliopsida)","Invertebrates","Amphibia",
               "Pisces","Reptilia","Mammalia","Aves")

## Sample Count
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
        legend.text = element_text(size=12, colour = "black"),
  ) +
  facet_wrap(~factor(Waterbody,levels = WB_level_order))

Figure_S3
ggsave("Tables and Figures/Figure_S3.png", plot = Figure_S3, 
       width = 10, height = 7, 
       units = "in", dpi = 300)
# -----------------------------------------------------------------------------
######## Figure S4 ############################################################
waterbody_sampling_df <- final_imputed_data %>% group_by(Sampling.Year,
                                                         Trophic_Level,
                                                         Waterbody) %>% 
  summarise(sample_count = sum(n_samples))

Trophic_Level_order<-c("Primary Producer","Primary Consumer",
                       "Secondary Consumer","Tertiary Consumer",
                       "Quaternary Consumer","Piscivorous/Insectivorous Bird",
                       "Apex Predator")
WB_level_order<-c("Lake Superior","Lake Michigan","Lake Huron",
                  "Lake Erie","Lake Ontario")
## Sample Count
Figure_S4<-
  ggplot(waterbody_sampling_df)+ 
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
Figure_S5 <-
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
  theme(axis.title.y = element_text(size=14, face="bold", colour = "black")
  )

Figure_S5
ggsave("Tables and Figures/Figure_S5.png", plot = Figure_S5, 
       width = 10, height = 7, 
       units = "in", dpi = 300)
# -----------------------------------------------------------------------------
######## Figure S9 ############################################################
Figure_S9 <- 
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

Figure_S9
ggsave("Tables and Figures/Figure_S9.png", plot = Figure_S9, width = 10, 
       height = 7,  units = "in", dpi = 300)

# -----------------------------------------------------------------------------
######## Figure S7 ############################################################
# Lake Superior Plot
Superior <- subset(PFOS_SY_means,Waterbody=="Lake Superior")
Superior_data <- subset(PFOS,Waterbody=="Lake Superior")

LS_SY_plot<-ggplot(Superior,aes(x = Sampling.Year, y=response))+
  annotate("rect", xmin=2000, xmax=2002, ymin=0,
           ymax=Inf,
           alpha=0.5,fill="orange") +
  geom_jitter(data = Superior_data,
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

# Lake Michigan Plot
Michigan <- subset(PFOS_SY_means,Waterbody=="Lake Michigan")
Michigan_data <- subset(PFOS,Waterbody=="Lake Michigan") 

LM_SY_plot<-ggplot(Michigan, aes(x = Sampling.Year,y=response))+
  annotate("rect", xmin=2000, xmax=2002, ymin=0,
           ymax=Inf,
           alpha=0.5,fill="orange") +
  geom_jitter(data = subset(PFOS,Waterbody=="Lake Michigan"),
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
Huron <- subset(PFOS_SY_means,Waterbody=="Lake Huron")
Huron_data <- subset(PFOS,Waterbody=="Lake Huron")

LH_SY_plot<-ggplot(Huron,aes(x = Sampling.Year, y=response))+
  annotate("rect", xmin=2000, xmax=2002, ymin=0,
           ymax=Inf,
           alpha=0.5,fill="orange") +
  geom_jitter(data = Huron_data,
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
Erie <- subset(PFOS_SY_means,Waterbody=="Lake Erie")
Erie_data <- subset(PFOS,Waterbody=="Lake Erie")

LE_SY_plot<-ggplot(Erie, aes(x = Sampling.Year, y=response)) +
  annotate("rect", xmin=2000, xmax=2002, ymin=0,
           ymax=Inf,
           alpha=0.5,fill="orange") +
  geom_jitter(data = Erie_data,
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
Ontario <- subset(PFOS_SY_means,Waterbody=="Lake Ontario")
Ontario_data <- subset(PFOS,Waterbody=="Lake Ontario")

LO_SY_plot<-ggplot(Ontario, aes(x = Sampling.Year,y=response))+
  annotate("rect", xmin=2000, xmax=2002, ymin=0,
           ymax=Inf,
           alpha=0.5,fill="orange") +
  geom_jitter(data = Ontario_data,
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

Figure_S7 <- ggarrange(LS_SY_plot, LM_SY_plot,LH_SY_plot,LE_SY_plot,LO_SY_plot,
          ncol = 3,nrow = 2,legend = "none")

Figure_S7
ggsave("Tables and Figures/Figure_S7.png", plot = Figure_S7, width = 14, 
       height = 10,  units = "in", dpi = 300)

# -----------------------------------------------------------------------------




# Summary statistics for the different variables
final_imputed_data %>% group_by(Waterbody) %>% summarise(sample_count=length(Sample.ID))

final_imputed_data %>% group_by(Trophic_Level) %>% summarise(sample_count=length(Sample.ID))

final_imputed_data %>% group_by(Class) %>% summarise(samp=length(Sample.ID),species = length(unique(Revised_Species)))

final_imputed_data %>% group_by(Waterbody_Type) %>% summarise(sample_count=length(Sample.ID))

emmeans(full_PFOS_gam, 
        specs = ~ Waterbody,
        type = "response",tran = "log10") |>
  as_tibble()

a<-emmeans(full_PFOS_gam, 
        specs = pairwise ~ Revised_Tissue,
        type = "response",tran = "log10",adjust="tukey") 
as_tibble(a$contrasts)
