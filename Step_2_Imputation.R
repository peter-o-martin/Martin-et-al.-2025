# Step_2_Imputation.R
# Contains all imputation code for the PFAS analytes chosen at the beginning of this script

# Written by Peter Martin
# Created December 13, 2024
# Finalized December 13, 2024

# Working directory
setwd("~/Desktop/Publications/Leyerle Martin et al., 2025")

## Packages
# Data formatting and combining
library(stringr)

# Imputation
library(MASS)
library(survival)
library(NADA)
library(truncnorm)
library(zCompositions)

# Visualizing data using maps
library(mapview)
library(leafsync)

## Functions
source("PFAS_Review_supportingFunctions.R") # Load supporting functions

## Loading in data frame created from Step 1
labeled_data<-read.csv("Step_1_labeled_data.csv",header = TRUE)

###############################################################################
#################### All biota (no separation of Classes) #####################
###############################################################################
###### Selection of PFAS Compounds that will be imputed and modeled ###########
# Create a data frame (X) that will hold values for each contaminant. Namely, (a) the number of <LOD, ND, NQ samples, (b) the total number of samples, (c) the detection frequency, the ratio of (a):(b), expressed as a percentage, and (d) the number of waterbodies represented in these samples
X<-matrix(data=0,nrow = 4,ncol = ncol(labeled_data)-23)
rownames(X)<-c("# of Unquantified Samples","# of Samples","% Quantified","# of Waterbodies")

# for loops which will calculate these values
for (i in 1:nrow(labeled_data)){
  for (j in 24:ncol(labeled_data)){
    if(grepl("<",labeled_data[i,j])==TRUE || grepl("N",labeled_data[i,j])==TRUE){
      X[1,j-23]<-X[1,j-23]+1
    }
    X[2,j-23]<-sum(!is.na(labeled_data[,j]))
  }
}
for (j in 24:ncol(labeled_data)){
  wtrbdy_count<-labeled_data[is.na(labeled_data[,j])==FALSE,]
  X[4,j-23]<-length(unique(wtrbdy_count$Waterbody))
}
rm(wtrbdy_count)
X[3,]<-((X[2,]-X[1,])/X[2,])*100

colnames(X)<-colnames(labeled_data[,24:ncol(labeled_data)])
View(X)

# Paste the contaminant with the % Quantified value and the total number of 
# samples
paste(colnames(X),"——",X[3,],"% ","——",X[2,])

# Use our criterion (% Quantified >= 80%) to select out contaminants 
# for further analysis
choice_X<-X[3,]
choice_X<-choice_X[choice_X>=60]
choice_X

# Our second criterion (minimum of 3 quantified samples for each waterbody, and at least two waterbodies represented in the samples), is not satisfied for the last two contaminants (6:6 PFPIA, 6:8 PFPIA). We therefore remove these two from the subset that we have pulled out for further analysis

choice_X<-choice_X[1:10]
choice_X

X<-X[,names(choice_X)]
paste(names(choice_X),"——",choice_X,"% ","——",X[2,],"——",X[4,])
X

# Subset of labeled_data with only analytes that meet the threshold
final_data<-labeled_data[,c(colnames(labeled_data[,1:23]),names(choice_X))]

###### Pre-Imputation Formatting ##############################################
# Prep data frame (dl_analyte) that will be used to hold 
# detection/quantification limits during imputation
dl_analyte<-final_data
final_data<-subset(final_data,is.na(Waterbody)==FALSE)

# For loop for dl_analyte that, when the cell entry in final_data is a 
# quantified value or is entered as NQ/NA/ND (and therefore cannot be imputed),
# enters a corresponding value of 0 in dl_analyte. If the value in final_data
# is a value in the form of "<[number]" (i.e., a value that can and needs to
# be imputed), the number (i.e., the limit) is entered into the corresponding 
# cell in dl_analyte
for (i in 1:nrow(dl_analyte)){
  for (j in 24:ncol(dl_analyte)){
    if(grepl("<",dl_analyte[i,j])==FALSE && is.na(dl_analyte[i,j])==FALSE &&
       grepl("N",dl_analyte[i,j])==FALSE){
      dl_analyte[i,j]<-0
    }
    else if(grepl("<",dl_analyte[i,j])==TRUE){
      dl_analyte[i,j]<-substr(final_data[i,j],2,str_length(final_data[i,j]))
    }
    else if (is.na(dl_analyte[i,j])==TRUE || grepl("N",dl_analyte[i,j])==TRUE){
      dl_analyte[i,j]<-0
    }
  }
}

# Calculate geometric means for each column so that the NA values in each
# column of final_data can be replaced with the corresponding geometric mean
# (since the functions in the imputation package zCompositions cannot handle
# NA values)
final_data[,24:ncol(final_data)]<-lapply(final_data[,24:ncol(final_data)],as.numeric)
geo_mean<-matrix(data=0,nrow = 1,ncol = ncol(final_data)-23)

# Geometric means are calculated as e^mean(log scale concentration values),
# where all unquantified values have been removed (na.omit) prior to calculation
for (j in 24:ncol(final_data)){
  value_set<-na.omit(final_data[,j])
  value_set<-as.numeric(value_set)
  geo_mean[j-23]<-exp(mean(log(value_set)))
}
rm(value_set)

final_data<-labeled_data[,c(colnames(labeled_data[,1:23]),names(choice_X))]

# Now we prep the data frame final_data for imputation: for any values that need
# to be imputed (i.e., "<[number]"), we assign the label 0, and for NA/ND/NQ
# values, we standardize the entry in final_data as "NA"
for (i in 1:nrow(final_data)){
  for (j in 24:ncol(final_data)){
    if(grepl("<",final_data[i,j])==TRUE){
      final_data[i,j]<-0
    }
    else if (is.na(final_data[i,j])==TRUE || grepl("N",final_data[i,j])==TRUE){
      final_data[i,j]<-NA
    }
  }
}

# Make sure that R recognizes all concentration and limits data as numeric
# values (i.e., apply the as.numeric function)
final_data[,24:ncol(final_data)]<-lapply(final_data[,24:ncol(final_data)],
                                         as.numeric)
dl_analyte[,24:ncol(dl_analyte)]<-lapply(dl_analyte[,24:ncol(dl_analyte)],
                                         as.numeric)

# Remove one row that has all but one value below the limit of 
# detection/quantification
# final_data<-final_data[c(1:413,415:nrow(final_data)),]
# dl_analyte<-dl_analyte[c(1:413,415:nrow(dl_analyte)),] 

# Look at the patterns of censored data in the data frame using one of 
# zCompositions in-built functions, zPatterns
# This function gives the percentage of left-censored values at the top,
# and, in the grid, all the different patterns of left-censoring in the data,
# as well the frequency of each pattern (e.g., for 66.14% of the data frame, 
# all contaminants except PFPeDA were quantified)
Data.pattern.ID<-zPatterns(final_data[,c(24:ncol(final_data))],label=0,
                           bar.colors=c("#CAFF70","#A2CD5A"),
                           bar.labels=TRUE,cell.colors=c("green4","white"),
                           cell.labels=c("Nondetected","Observed"),
                           cex.axis=0.8)

###### Imputation Code ########################################################
# Set all values labeled as "NA" in the data frame final_data to the
# corresponding geometric mean of that column
for (i in 1:nrow(final_data)){
  for (j in 24:ncol(final_data)){
    if (is.na(final_data[i,j])==TRUE){
      final_data[i,j]<-geo_mean[j-23]
    }
  }
}

# Multiplicative simple replacement
data_multRepl<-multRepl(final_data[,24:ncol(final_data)],
                        label=0,
                        frac=0.5,
                        dl=dl_analyte[,24:ncol(dl_analyte)],
                        z.warning=1,
                        z.delete = FALSE,
                        closure = 10^9)

# Multiplicative lognormal replacement with maximum-likelihood estimates
data_multLN<-multLN(final_data[,24:ncol(final_data)],
                    label=0,
                    dl=dl_analyte[,24:ncol(dl_analyte)],
                    z.warning=1,
                    z.delete = FALSE)

# Multiplicative lognormal replacement with robust estimates
data_multLNrob<-multLN(final_data[,24:ncol(final_data)],
                       label=0,
                       dl=dl_analyte[,24:ncol(dl_analyte)],
                       rob = TRUE,
                       z.warning=1,
                       z.delete = FALSE)

# Multiplicative lognormal replacement with random values below the limit of
# detection
data_multLNrand<-multLN(final_data[,24:ncol(final_data)],
                        label=0,
                        dl=dl_analyte[,24:ncol(dl_analyte)],
                        random = TRUE,
                        z.warning=1,
                        z.delete = FALSE)

# Because the two advanced imputation methods need one column with no values
# that need to be imputed, we used the contaminant with the lowest
# censoring frequency: PFOS (see zPatterns function output). First,
# the complete PFOS values generated from multiplicative simple replacement 
# are bound to the contaminant columns where imputation is needed: the data
# frame adv_impute_PFOS
non_PFOS_analytes<-colnames(subset(final_data[,c(24:ncol(final_data))],select = -PFOS))

adv_impute_PFOS<-cbind(final_data[,1:23],
                       data_multRepl$PFOS,
                       subset(final_data[,-c(1:23)],select = -PFOS))

colnames(adv_impute_PFOS)<-colnames(final_data)
PFOS_dl_analyte<-dl_analyte
PFOS_dl_analyte$PFOS<-0

# Log-ratio Expectation-Maximization algorithm
# adv_impute_PFOS is run through lrEM to generate imputed values for the rest
# of the contaminants, and then these values are bound to the original PFOS
# column so that lrEM imputation can be performed for the PFOS values as well
data_lrEM<-lrEM(adv_impute_PFOS[,24:ncol(adv_impute_PFOS)],
                label=0,
                dl=PFOS_dl_analyte[,24:ncol(PFOS_dl_analyte)],
                ini.cov = "complete.obs",
                max.iter = 100,
                z.warning = 0.9,
                z.delete=FALSE)

final_data_lrEM<-cbind(final_data$PFOS,subset(data_lrEM,select = -PFOS))
PFOS_dl_analyte<-dl_analyte
PFOS_dl_analyte[,non_PFOS_analytes]<-0

final_data_lrEM<-lrEM(final_data_lrEM[,1:ncol(final_data_lrEM)],
                      label=0,
                      dl=PFOS_dl_analyte[,24:ncol(PFOS_dl_analyte)],
                      ini.cov = "complete.obs",
                      max.iter = 100,
                      z.warning = 0.9,
                      z.delete=FALSE)

colnames(final_data_lrEM)<-names(data_multRepl)

# For the final imputation method (our imputation method of choice), we will 
# now use the PFOS values, completed with the EM-algorithm, as our complete 
# column, and we will impute the other contaminants with the lrDA function
adv_impute_PFOS<-cbind(final_data[,1:23],
                       final_data_lrEM$PFOS,
                       subset(final_data[,-c(1:23)],select = -PFOS))
colnames(adv_impute_PFOS)<-colnames(final_data)

PFOS_dl_analyte<-dl_analyte
PFOS_dl_analyte$PFOS<-0

# Log-ratio Data Augmentation algorithm
data_lrDA<-lrDA(adv_impute_PFOS[,24:ncol(adv_impute_PFOS)],
                label=0,
                dl=PFOS_dl_analyte[,24:ncol(PFOS_dl_analyte)],
                ini.cov="lrEM",
                m=5,
                closure = 10^9,
                z.warning = 0.9,
                z.delete=FALSE)

final_data_lrDA<-cbind(final_data$PFOS,subset(data_lrDA,select = -PFOS))
PFOS_dl_analyte<-dl_analyte
PFOS_dl_analyte[,non_PFOS_analytes]<-0

final_data_lrDA<-lrDA(final_data_lrDA[,1:ncol(final_data_lrDA)],
                      label=0,
                      dl=PFOS_dl_analyte[,24:ncol(PFOS_dl_analyte)],
                      ini.cov="lrEM",
                      m=5,
                      closure = 10^9,
                      z.warning = 0.9,
                      z.delete=FALSE)

colnames(final_data_lrDA)<-names(data_multRepl)

###### Post-Imputation Formatting #############################################
# Replace geometric mean entries in each data frame with NA (restoring to
# pre-imputation form)
for (i in 1:nrow(final_data)){
  for (j in 24:ncol(final_data)){
    if (identical(final_data[i,j],geo_mean[j-23])){
      final_data[i,j]<-NA
    }
  }
}
for (i in 1:nrow(data_lrDA)){
  for (j in 1:ncol(data_lrDA)){
    if (identical(data_lrDA[i,j],geo_mean[j])){
      data_lrDA[i,j]<-NA
    }
  }
}
for (i in 1:nrow(data_multLNrand)){
  for (j in 1:ncol(data_multLNrand)){
    if (identical(data_multLNrand[i,j],geo_mean[j])){
      data_multLNrand[i,j]<-NA
    }
  }
}
for (i in 1:nrow(data_multRepl)){
  for (j in 1:ncol(data_multRepl)){
    if (identical(data_multRepl[i,j],geo_mean[j])){
      data_multRepl[i,j]<-NA
    }
  }
}
for (i in 1:nrow(data_multLN)){
  for (j in 1:ncol(data_multLN)){
    if (identical(data_multLN[i,j],geo_mean[j])){
      data_multLN[i,j]<-NA
    }
  }
}
for (i in 1:nrow(data_multLNrob)){
  for (j in 1:ncol(data_multLNrob)){
    if (identical(data_multLNrob[i,j],geo_mean[j])){
      data_multLNrob[i,j]<-NA
    }
  }
}
for (i in 1:nrow(final_data_lrEM)){
  for (j in 1:ncol(final_data_lrEM)){
    if (identical(final_data_lrEM[i,j],geo_mean[j])){
      final_data_lrEM[i,j]<-NA
    }
  }
}
for (i in 1:nrow(final_data_lrDA)){
  for (j in 1:ncol(final_data_lrDA)){
    if (identical(final_data_lrDA[i,j],geo_mean[j])){
      final_data_lrDA[i,j]<-NA
    }
  }
}

# Paste results together from each imputation method for one of the patterns of censored values (zPatterns output contained in the Data.pattern.ID variable)
# Allows for comparisons of different imputation methods
rbind(final_data[Data.pattern.ID==94,c(24:33)],
      dl_analyte[Data.pattern.ID==94,c(24:33)],
      data_multRepl[Data.pattern.ID==94,],
      data_multLN[Data.pattern.ID==94,],
      data_multLNrob[Data.pattern.ID==94,],
      data_multLNrand[Data.pattern.ID==94,],
      final_data_lrEM[Data.pattern.ID==94,],
      data_lrDA[Data.pattern.ID==94,],
      final_data_lrDA[Data.pattern.ID==94,])


###### Save data frame and delete excess variables ############################
# Remove data frames produced by imputation methods that we are not going to use
rm(data_multRepl,
   data_multLN,
   data_multLNrob,
   data_multLNrand,
   adv_impute_PFOS,
   data_lrEM,
   final_data_lrEM,
   data_lrDA,
   PFOS_dl_analyte,
   geo_mean,
   dl_analyte,
   final_data)

# Delete other now extraneous variables
rm(i,j,
   choice_X,
   Data.pattern.ID,
   labeled_data,
   non_PFOS_analytes)

# Create final data frame for analysis, using lrDA() derived values and save
# this data frame for all downstream analysis
final_imputed_data<-cbind(final_data[,1:23],final_data_lrDA)

setdiff(final_imputed_data$PFOS,final_data$PFOS) # Check that the imputed data was properly assigned to the final_imputed_data variable
rownames(final_imputed_data)<-NULL

write.table(final_imputed_data,file = "Step_2_final_imputed_data.csv",sep = ",",row.names = FALSE)
write.table(final_data_lrDA,file = "Step_2_lrDA_imputed_values.csv",sep = ",",row.names = FALSE)



