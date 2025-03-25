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
# library(lattice)
# library(caret)
# 
# # Modeling
# library(nlme)
# library(mgcv)
# library(DHARMa)
# 
# # Visualizing model output and calculating results
# library(gratia)
# library(carData)
# library(car)
library(emmeans)
# library(ggplot2)
# library(dotwhisker)

## Functions
source("PFAS_Review_supportingFunctions.R") # Load supporting functions

## Loading in finalized data frame from Step 5
final_imputed_data <- read.csv("Step_5_final_imputed_data.csv",header = TRUE)

## Loading in the models from Step 5
load("full_PFOS_gam.Rdata") 
load("full_PFNA_gam.Rdata")
load("full_PFDA_gam.Rdata")
load("full_PFUnA_gam.Rdata")
load("full_PFDoA_gam.Rdata")
load("full_PFTrDA_gam.Rdata")


model_list <- list(full_PFOS_gam,full_PFNA_gam,full_PFDA_gam,full_PFUnA_gam,
                   full_PFDoA_gam,full_PFTrDA_gam)

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
Table_2a_final$group <- "Inter-Tissue PFAS Abundance (%)"

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
Table_2b_final[]<-paste0(Table_2b_final,"
                         (",Table_2b_PFOS,")")
Table_2b_final <- Table_2b_final |> as_tibble()
Table_2b_final$rowname <- c("Eggs","Blood","Liver","Combined Tissue and Blood")
Table_2b_final$group <- "Intra-Tissue PFAS Composition (%)"

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
Table_2c_final$group <- "Egg-Standardized Intra-Tissue PFAS Composition"

gt(rbind(Table_2a_final,Table_2b_final,Table_2c_final),
   rowname_col = "rowname",
   groupname_col = "group",
   row_group_as_column = TRUE) |> 
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
  )

gtsave(Table_2_final,"tab_1.rtf")

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
