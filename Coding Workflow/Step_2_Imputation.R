# Step_2_Imputation.R
# Contains all imputation code for the PFAS analytes chosen at the beginning of 
# this script

# Working directory
setwd("~/Desktop/Publications/Leyerle Martin et al., 2025")

## Packages
# Loading the correct version
library(remotes)

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

###### Selection of PFAS Compounds that will be imputed and modeled ###########
# Create a data frame (X) that will track key values for each contaminant. Namely, 
# for each contaminant, this data frame will record 
# (i) the number of <LOD (less than the limit of detection), <LOQ (less than the 
# limit of quantification), ND (not detected), and NQ (not quantified) samples,
# (ii) the total number of samples, 
# (iii) the detection frequency, the ratio of (i):(ii), expressed as a percentage,
# and (iv) the number of Great Lakes watersheds (maximum of five) represented in
# these samples
# Setting up the matrix
X<-matrix(data=0,nrow = 4,ncol = ncol(labeled_data)-23)
rownames(X)<-c("# of Unquantified Samples","# of Samples","% Quantified","# of Waterbodies")

# for loops which will calculate the values described above
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

# Label the columns of this matrix with the correct contaminant name
colnames(X)<-colnames(labeled_data[,24:ncol(labeled_data)])
View(X)

# Display the matrix with each contaminant's detection frequency and the total
# number of samples
paste(colnames(X),"——",X[3,],"% ","——",X[2,])

# Use our criterion (detection frequency ≥ 60%) to select out contaminants 
# for further analysis
# N.B. While our final criterion, a detection frequency of 80% or greater, selected
# out 6 contaminants for modeling (PFOS, PFNA, PFDA, PFUnA, PFDoA, and PFTrDA),
# we used the more lenient criterion of ≥ 60% for imputation (which led to the
# inclusion of PFDS, PFEtCHxS, PFTeDA, and PFPeDA in this part of the analysis).
# This decision allowed us to provide a larger dataset to the model, with more
# information available on the covariance structure between different PFAS
# contaminants that the algorithms can use in estimating <LOD/<LOQ concentrations.
choice_X<-X[3,]
choice_X<-choice_X[choice_X>=60]
choice_X

# Our second criterion (minimum of 3 quantified samples for each watershed, and at
# least two watersheds represented in the samples), is not satisfied for the last
# two contaminants (6:6 PFPIA, 6:8 PFPIA). We therefore remove these two PFAS from
# the subset that we have pulled out for imputation.
choice_X<-choice_X[1:10]
choice_X

X<-X[,names(choice_X)]
paste(names(choice_X),"——",choice_X,"% ","——",X[2,],"——",X[4,])
X

# We now subset the labeled_data data frame from Step 1 to select only the 
# analytes that met our criteria
final_data<-labeled_data[,c(colnames(labeled_data[,1:23]),names(choice_X))]

###### Pre-Imputation Formatting ##############################################
# Prep the data frame (dl_analyte) that will be used to store 
# detection/quantification limits during imputation (the limit values will be used
# by the algorithms)
dl_analyte<-final_data

# A for loop for dl_analyte that, when the cell entry in final_data is a 
# quantified value or is entered as NQ/NA/ND (and therefore cannot be imputed),
# enters a value of 0 in the corresponding cell of dl_analyte. If the value in
# final_data is a value in the form of "<[number]" (i.e., a value that can be
# imputed), that numerical value is entered into the corresponding cell in
# dl_analyte
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

# Because the imputation algorithms we are using (the lrDA and lrEM functions from
# the zCompositions package) cannot handle NA values in the data frame, we calculate
# the geometric mean concentration for each column (i.e., for each of the 10 PFAS)
# and then replace all NA values in each column with the corresponding mean
# We first remove all non-numeric entries in the PFAS columns by applying the 
# function as.numeric() with lapply() to those columns. We then create the vector
# that will store the geometric mean values for each PFAS column
final_data[,24:ncol(final_data)]<-lapply(final_data[,24:ncol(final_data)],
                                         as.numeric)
geo_mean<-matrix(data=0,nrow = 1,ncol = ncol(final_data)-23)

# Geometric means are calculated as e^mean(ln[concentration values]), where all NA
# values have been removed (using the na.omit() function) prior to calculation
for (j in 24:ncol(final_data)){
  value_set<-na.omit(final_data[,j])
  value_set<-as.numeric(value_set)
  geo_mean[j-23]<-exp(mean(log(value_set)))
}
rm(value_set)

# Undo the actions of as.numeric() on our data frame and restore the entries like
# "<[concentration]" that are interpreted by R as character strings
final_data<-labeled_data[,c(colnames(labeled_data[,1:23]),names(choice_X))]

# Now we prep the data frame final_data for imputation: for any values that need
# to be imputed (i.e., "<[concentration]"), we assign the label 0, and for NA/ND/NQ
# values that cannot be imputed, we standardize the entry in final_data as "NA"
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
# values (i.e., apply the as.numeric function again)
final_data[,24:ncol(final_data)]<-lapply(final_data[,24:ncol(final_data)],
                                         as.numeric)
dl_analyte[,24:ncol(dl_analyte)]<-lapply(dl_analyte[,24:ncol(dl_analyte)],
                                         as.numeric)

# Look at the patterns of censored values in our data frame, using one of 
# zCompositions built-in functions, zPatterns()
# This function gives the percentage of left-censored values (labeled as 0) at the 
# top, and, in the grid, all the different patterns of left-censoring in the data,
# as well the frequency of each pattern
Data.pattern.ID<-zPatterns(final_data[,c(24:ncol(final_data))],label=0,
                           bar.colors=c("#CAFF70","#A2CD5A"),
                           bar.labels=TRUE,cell.colors=c("green4","white"),
                           cell.labels=c("Nondetected","Observed"),
                           cex.axis=0.8)

###### Imputation Code ########################################################
# Set all values labeled as "NA" in the columns of the data frame final_data to the
# corresponding geometric mean of that column
for (i in 1:nrow(final_data)){
  for (j in 24:ncol(final_data)){
    if (is.na(final_data[i,j])==TRUE){
      final_data[i,j]<-geo_mean[j-23]
    }
  }
}

# First we conduct, multiplicative simple replacement of missing values in all six 
# columns, using the function multRepl(), where all left-censored values are
# replaced with half of their analytical limit (i.e., 0.5*LOD or 0.5*LOQ)
data_multRepl<-multRepl(final_data[,24:ncol(final_data)],
                        label=0,
                        frac=0.5,
                        dl=dl_analyte[,24:ncol(dl_analyte)],
                        z.warning=1,
                        z.delete = FALSE,
                        closure = 10^9)

# Because the two advanced imputation methods need one "reference" column with no 
# values that need to be imputed, we used the contaminant with the lowest
# censoring frequency, PFOS, as this reference. To set up PFOS as the 
# reference contaminant, we bind the results of multiplicative simple replacement 
# (multRepl() function) for the PFOS column to the remaining contaminant columns
# that need to be imputated: we call this new data frame adv_impute_PFOS
# First we create a vector with the names of all non-PFOS contaminants in our data
# frame that are going to be imputed
non_PFOS_analytes<-colnames(subset(final_data[,c(24:ncol(final_data))],
                                   select = -PFOS))

# Then we create the new data frame that has PFOS set up as a reference column
adv_impute_PFOS<-cbind(final_data[,1:23],
                       data_multRepl$PFOS,
                       final_data[,non_PFOS_analytes])
colnames(adv_impute_PFOS)<-colnames(final_data)

# We then create a new version of dl_analyte for adv_impute_PFOS, which has all
# non-zero entries removed from the PFOS column (i.e., set to 0) so that the 
# algorithms don't try to impute any values from this reference column
PFOS_dl_analyte<-dl_analyte
PFOS_dl_analyte$PFOS<-0

# Log-ratio Expectation-Maximization algorithm
# adv_impute_PFOS is run through lrEM() to generate imputed values for the rest
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

# Now reset PFOS_dl_analyte and then remove all non-zero entries from the non-PFOS 
# columns so that values in those columns won't be accidentally imputed by the 
# algorithms. Only the PFOS column will be imputed, and the other columns will 
# serve as references for the algorithms.
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

# For the final imputation method, we will now use the lrEM()-imputed PFOS column 
# as our reference column, and we will impute the other nine contaminants using a 
# Log-Ratio Data Augmentation algorithm from the lrDA() function
adv_impute_PFOS<-cbind(final_data[,1:23],
                       final_data_lrEM$PFOS,
                       final_data[,non_PFOS_analytes])
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

colnames(data_lrDA)<-names(data_multRepl)

###### Post-Imputation Formatting #############################################
# Replace geometric mean entries in each data frame with NA (thus restoring the 
# data frames to pre-imputation form)
## Data frame that was given as input to the multRepl() function
for (i in 1:nrow(final_data)){
  for (j in 24:ncol(final_data)){
    if (identical(final_data[i,j],geo_mean[j-23])){
      final_data[i,j]<-NA
    }
  }
}
## multRepl() output
for (i in 1:nrow(data_multRepl)){
  for (j in 1:ncol(data_multRepl)){
    if (identical(data_multRepl[i,j],geo_mean[j])){
      data_multRepl[i,j]<-NA
    }
  }
}
## lrEM() output
for (i in 1:nrow(final_data_lrEM)){
  for (j in 1:ncol(final_data_lrEM)){
    if (identical(final_data_lrEM[i,j],geo_mean[j])){
      final_data_lrEM[i,j]<-NA
    }
  }
}
## lrDA() output
for (i in 1:nrow(data_lrDA)){
  for (j in 1:ncol(data_lrDA)){
    if (identical(data_lrDA[i,j],geo_mean[j])){
      data_lrDA[i,j]<-NA
    }
  }
}


# Paste results together from each imputation method for one of the patterns of
# censored values (patterns are given in the zPatterns() output, contained in the 
# Data.pattern.ID object). This table allows us to compare the estimates derived 
# from different imputation methods
estimates_table <- 
  rbind(final_data[Data.pattern.ID==94,c(24:33)],
      dl_analyte[Data.pattern.ID==94,c(24:33)],
      data_multRepl[Data.pattern.ID==94,],
      final_data_lrEM[Data.pattern.ID==94,],
      data_lrDA[Data.pattern.ID==94,])

rownames(estimates_table)<-c("Entry in the Data Frame (0s to be imputed)",
                             "Analytical Limit (when applicable)",
                             "multRepl() results",
                             "lrEM() results",
                             "lrDA() results")
estimates_table

###### Save data frame and delete excess variables ############################
# Create final data frame for analysis, using the output of the lrDA() function,
# and save this data frame for all downstream analysis
final_imputed_data<-cbind(final_data[,1:23],data_lrDA)

setdiff(final_imputed_data$PFOS,final_data$PFOS) # Check that the imputed data 
# was properly assigned to the final_imputed_data variable
rownames(final_imputed_data)<-NULL

write.table(final_imputed_data,file = "Step_2_final_imputed_data.csv",
            sep = ",",row.names = FALSE)

# Also save the output of lrDA() directly, in case this output needs to be
# referenced later on
write.table(data_lrDA,file = "Step_2_lrDA_imputed_values.csv",
            sep = ",",row.names = FALSE)


# Remove data frames produced by imputation methods that we are not going to use
rm(data_multRepl,
   data_lrEM,
   final_data_lrEM,
   adv_impute_PFOS,
   PFOS_dl_analyte,
   geo_mean,
   dl_analyte,
   final_data)

# Delete other now extraneous variables
rm(i,j,
   choice_X,
   Data.pattern.ID,
   labeled_data,
   non_PFOS_analytes,
   estimates_table)

