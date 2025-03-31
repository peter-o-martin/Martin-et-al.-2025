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
# library(lattice)
library(caret)
# 
# # Modeling
# library(nlme)
library(mgcv)
# library(DHARMa)
# 
# # Visualizing model output and calculating results
# library(gratia)
# library(carData)
# library(car)
library(emmeans)
library(ggrepel)
library(ggpubr)
library(rnaturalearth)
library(sf)
library(ggspatial)
library(RColorBrewer)
# library(ggplot2)
# library(dotwhisker)

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

######## Table 2 ##############################################################
Table_2a <- matrix(nrow = 4, ncol = 6,
                  dimnames = list(
                    c("Eggs","Blood","Liver","Combined Tissue and Blood"),
                    c("PFOS (C8)","PFNA (C9)","PFDA (C10)","PFUnA (C11)",
                     "PFDoA (C12)","PFTrDA (C13)")
                    ))
Table_2b <- Table_2a
Table_2b_PFOS <- Table_2a

##### Table 2a
for(i in 1:length(model_list)){
  tissue_estimates <- emmeans(model_list[[i]],
                              specs = ~ Revised_Tissue,
                              type='response',tran = "log10") |>
    as_tibble()
  tissue_estimates <- tissue_estimates %>% select(Revised_Tissue,response) %>% 
    mutate(t_abundance = 100*round(response/sum(response),digits=5))
  
  Table_2a[1:3,i] <- tissue_estimates$t_abundance[6:4]
  Table_2a[4,i] <- sum(tissue_estimates$t_abundance[1:5])
}

Table_2a_final <- Table_2a |> as_tibble()
Table_2a_final$rowname <- c("Eggs","Blood","Liver","Combined Tissue and Blood")
Table_2a_final$group <- "(A) Inter-Tissue PFAS Abundance (%)"

##### Table 2b
composition_estimates <- vector("list",5)
for(i in 2:length(model_list)){
  composition_estimates[[i-1]] <- emmeans(model_list[[i]],
                              specs = ~ Revised_Tissue,
                              type='response',tran = "log10") |>
    as_tibble()
}
composition_sums <- rowSums(data.frame(sapply(composition_estimates, `[[`, 2)))
composition_sums[length(composition_sums) + 1] <- 
  sum(composition_sums[1:5])

composition_estimates <- data.frame(
  sapply(composition_estimates, `[[`, 2)) 

composition_estimates[nrow(composition_estimates) + 1,] <- 
  colSums(composition_estimates[1:5,])
composition_estimates <- 
  round(composition_estimates / composition_sums, digits = 5) * 100

for(i in 2:length(model_list)){
  Table_2b[1:3,i] <- composition_estimates[6:4,i-1]
  Table_2b[4,i] <- sum(composition_estimates[7,i-1])
}

# Table 2b (PFOS included)
PFOS_composition_estimates <- vector("list",6)
for(i in 1:length(model_list)){
  PFOS_composition_estimates[[i]] <- emmeans(model_list[[i]],
                                        specs = ~ Revised_Tissue,
                                        type='response',tran = "log10") |>
    as_tibble()
}
PFOS_composition_sums <- rowSums(
  data.frame(sapply(PFOS_composition_estimates, `[[`, 2)))
PFOS_composition_sums[length(PFOS_composition_sums) + 1] <- 
  sum(PFOS_composition_sums[1:6])

PFOS_composition_estimates <- data.frame(
  sapply(PFOS_composition_estimates, `[[`, 2)) 

PFOS_composition_estimates[nrow(PFOS_composition_estimates) + 1,] <- 
  colSums(PFOS_composition_estimates[1:6,])
PFOS_composition_estimates <- 
  round(PFOS_composition_estimates / PFOS_composition_sums, digits = 5) * 100

for(i in 1:length(model_list)){
  Table_2b_PFOS[1:3,i] <- PFOS_composition_estimates[6:4,i]
  Table_2b_PFOS[4,i] <- sum(PFOS_composition_estimates[7,i])
}

Table_2b_final <- Table_2b
Table_2b_final[]<-paste0(Table_2b_final," ", "<br>","(",Table_2b_PFOS,")")
Table_2b_final <- Table_2b_final |> as_tibble()
Table_2b_final$rowname <- c("Eggs","Blood","Liver","Combined Tissue and Blood")
Table_2b_final$group <- "(B) Intra-Tissue PFAS Composition (%)"

##### Table 2c
Table_2c <- round(rbind(Table_2b[2,] / Table_2b[1,],
                  Table_2b[3,] / Table_2b[1,],
                  Table_2b[4,] / Table_2b[1,]),digits = 3)

Table_2c <- matrix(Table_2c,nrow = 3, ncol = 6,
                   dimnames = list(
                     c("Blood:Eggs","Liver:Eggs","Combined Tissue:Eggs"),
                     c("PFOS (C8)","PFNA (C9)","PFDA (C10)","PFUnA (C11)",
                       "PFDoA (C12)","PFTrDA (C13)")
                   ))

Table_2c_final <- Table_2c |> as_tibble()
Table_2c_final$rowname <- c("Blood:Eggs","Liver:Eggs","Combined Tissue:Eggs")
Table_2c_final$group <- "(C) Egg-Standardized Intra-Tissue PFAS Composition"

Table_2_final<-gt(rbind(Table_2a_final,Table_2b_final,Table_2c_final),
   rowname_col = "rowname",
   groupname_col = "group",
   row_group_as_column = F) |> 
  tab_source_note(source_note = md("Values are derived from tissue-specific estimated marginal means, which were calculated over the reference grid of all possible factor combinations using the **emmeans::emmeans()** function (Lenth, 2024)")) |>
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
  fmt_markdown(columns = TRUE)

Table_2_final
gtsave(Table_2_final,"Tables and Figures/Table_2.rtf")

# -----------------------------------------------------------------------------

######## Figures ##############################################################
######## Conceptual Figure ####################################################
Great_Lakes_shp <- ne_states(country=c("canada","united states of america"),
                             returnclass = "sf")

Great_Lakes_watershed <- read_sf("~/Desktop/Rfiles/PFAS Review Paper/Full Watershed Shape",
                                 "Full_Watershed_PFAS_Paper")

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

PFOS$predict<-predict.gam(full_PFOS_gam,newdata=PFOS)
WB_level_order<-c("Lake Superior","Lake Michigan","Lake Huron",
                  "Lake Erie","Lake Ontario")
# color palette used
display.brewer.pal(n=8,"RdYlBu")


ggplot() +
  geom_sf(data=Great_Lakes_shp,color="grey") +
  geom_point(data=PFOS, aes(x=Longitude,y=Latitude,
                            fill = factor(Waterbody,levels = WB_level_order)),
             shape=21,size=0.75)+
  geom_sf(data=Great_Lakes_watershed,
          aes(fill=factor(merge,levels = WB_level_order)),color=NA,alpha=0.70) +
  scale_fill_manual(values = c("#4575B4","#91BFDB","#FEE090","#FC8D59",
                               "#D73027")) +
  coord_sf(xlim = c(-93.0,-74.5),ylim = c(40.5,51)) +
  annotation_north_arrow(location = "tl",which_north = "true", 
                         pad_x = unit(0.2, "cm"), pad_y = unit(0.2, "cm"),
                         style = north_arrow_nautical,width = unit(1.5, "cm"), 
                         height = unit(1.5, "cm")) +
  xlab("") +
  ylab("") +
  guides(fill= "none") +
  theme_classic(base_size = 14) +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())
# -----------------------------------------------------------------------------
######## Figure 1 #############################################################
Great_Lakes_shp <- ne_states(country=c("canada","united states of america"),
                             returnclass = "sf")
Great_Lakes_watershed <- read_sf("~/Desktop/Publications/Leyerle Martin et al., 2025/Great Lakes Shapefiles/Custom Shapefiles/Full Watershed Great Lakes",
                                 "GL_Watershed_shapefile")

PFOS$predict<-predict.gam(full_PFOS_gam,newdata=PFOS)
WB_level_order<-c("Lake Superior","Lake Michigan","Lake Huron",
                  "Lake Erie","Lake Ontario")

# https://bookdown.org/nicohahn/making_maps_with_r5/docs/ggplot2.html
Figure_1 <- ggplot() +
  geom_sf(data=Great_Lakes_shp) +
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
    legend.title = element_text(size=14, face="bold", colour = "black")
  )

Figure_1
ggsave("Figures/Figure_1.png", plot = Figure_1, width = 15, height = 10, 
       units = "in", dpi = 300)
# -----------------------------------------------------------------------------
######## Figure 2 #############################################################
WB_level_order<-c("Lake Superior","Lake Michigan","Lake Huron",
                  "Lake Erie","Lake Ontario")
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
LO_Sampling_Years<-subset(Sampling_Years,Sampling_Years>=min(LO_PFOS$Sampling.Year))
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
  theme_bw(base_size = 14) +
  ylab("PFOS concentration (ng/g w.w.)") +
  xlab("Sampling Year (1979-2021)") +
  coord_cartesian(ylim = c(0,300)) +
  geom_line(data = LM_spline, aes(x = x, y = y),color="#91BFDB",linewidth=1) +
  geom_line(data = LH_spline, aes(x = x, y = y),color="#FEE090",linewidth=1) +
  geom_line(data = LS_spline, aes(x = x, y = y),color="#4575B4",linewidth=1) +
  geom_line(data = LE_spline, aes(x = x, y = y),color="#FC8D59",linewidth=1) +
  geom_line(data = LO_spline, aes(x = x, y = y),color="#D73027",linewidth=1) +
  geom_line(aes(colour = Waterbody)) +
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
    legend.title = element_text(size=14, face="bold", colour = "black")
    )

Figure_2
ggsave("Figures/Figure_2.png", plot = Figure_2, width = 13, height = 8, 
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

Figure_3 <- ggarrange(
  ggplot(PFAS_WB_means, aes(y = factor(Waterbody, level = WB_level_order), 
                                    x = response)) +
    geom_col(aes(fill = contaminant)) +
    scale_fill_manual(values = c("#4575B4","#91BFDB","#E0F3F8",
                                 "#FEE090","#FC8D59",
                                 "#D73027")) +
    theme_bw(base_size = 14) +
    ylab(NULL) +
    xlab("Concentration (ng/g w.w.)") +
    geom_label(aes(x = sum,y = factor(Waterbody, level = WB_level_order),
                   label = paste0(round(sum,2)," ng/g w.w."),
                   group = factor(contaminant)),
               hjust = -0.25,fontface='bold',size=4) +
    coord_cartesian(xlim = c(0,225)) +
    guides(fill = guide_legend(title = "PFAS Compound")) +
    theme(legend.title=element_text(size=14,face = "bold"),
          axis.title.x = element_text(size=14, face="bold", colour = "black")),
  
  ggplot(PFAS_WB_means, aes(y = factor(Waterbody, level = WB_level_order), 
                            x = percent*100, fill = contaminant)) +
    geom_col() +
    scale_fill_manual(values = c("#4575B4","#91BFDB","#E0F3F8",
                                 "#FEE090","#FC8D59",
                                 "#D73027")) +
    geom_label_repel(aes(label=paste0(sprintf("%1.1f", percent*100),"%")),
                     position=position_stack(vjust=0.5),
                     force=0.01,force_pull = 190,
                     max.overlaps = 8,direction = "both",show.legend = FALSE,
                     fontface='bold',size=4) +
    ylab(NULL) +
    xlab("Contaminant Composition (%)") +
    theme_bw(base_size = 14) +
    guides(fill = guide_legend(title = "PFAS Compound")) +
    theme(legend.title=element_text(size=14,face = "bold"),
          axis.title.x = element_text(size=14, face="bold", colour = "black")),
  
  ncol = 1,
  nrow = 2,
  legend="right",
  common.legend = T
  )

Figure_3
ggsave("Figures/Figure_3.png", plot = Figure_3, width = 18, height = 12, 
       units = "in", dpi = 300)

# -----------------------------------------------------------------------------



emmeans(full_PFOS_gam,
        specs = ~ Waterbody,
        type='response',tran = "log10") |>
  as_tibble()

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
