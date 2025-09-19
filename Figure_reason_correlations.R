salloc -A def-grouleau --time=0-4 -c 1 --mem=200g
salloc -A def-sfarhan --time=0-8 -c 1 --mem=100g
salloc -A def-tdurcan --time=0-8 -c 1 --mem=40g

def-grouleau

module load StdEnv/2023
module load r/4.4.0
R

## Load libraries
library(MAST)
library(Seurat, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") ##good
library(withr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") ##good
library(ggplot2, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
library(backports, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
library(ggpubr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
library(tidyr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
library(stringr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
library(dplyr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
library(ggrepel, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
library(rstatix, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
library(remotes, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
library(labeling, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
library(tidyverse, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
library(tidyverse, lib="/home/fiorini9/scracth/R")
#install.packages("RColorBrewer", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
library(RColorBrewer, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
library(readr)
library(reshape2)

library(tidyverse, lib="/home/fiorini9/scracth/R")
#install.packages("RColorBrewer", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
library(RColorBrewer, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
library(readr)
library(reshape2)
library(ggh4x)
library(ggsankey)
library(ggalluvial)
library(purrr)

The final model performance dataframes are produced in Figure_model_accuracy4_generalizable.R

######################
## Read HVG in final model performance -- All cells. 
######################
## code
    # Read in file
    performance <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_final_performance_cell_level_top.csv')
    
    # Only retain overall performance
    performance <- subset(performance, group == "Overall")

    # Retain columns of interest
    performance <- performance %>% dplyr::select(input, celltype, mean_accuracy)
    nrow(performance)

    performance_HVG <- performance
##

######################
## Read NMF in final model performance -- All cells. 
######################
## code
    # Read in file
    performance_NMF <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/NMF_final_performance_cell_level_top.csv')
    
    # Only retain overall performance
    performance_NMF <- subset(performance_NMF, group == "Overall")

    # Retain columns of interest
    performance_NMF <- performance_NMF %>% dplyr::select(input, celltype, mean_accuracy)
    nrow(performance_NMF)
##

######################
## Read HVG BA4 in final model performance -- All cells. 
######################
## code
    # Read in file
    performance <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_BA4_final_performance_cell_level_top.csv')
    
    # Only retain overall performance
    performance <- subset(performance, group == "Overall")

    # Retain columns of interest
    performance <- performance %>% dplyr::select(input, celltype, mean_accuracy)
    nrow(performance)

    performance_HVG_BA4 <- performance
##

######################
## Read HVG BA9 in final model performance -- All cells. 
######################
## code
    # Read in file
    performance <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_BA9_final_performance_cell_level_top.csv')
    
    # Only retain overall performance
    performance <- subset(performance, group == "Overall")

    # Retain columns of interest
    performance <- performance %>% dplyr::select(input, celltype, mean_accuracy)
    nrow(performance)

    performance_HVG_BA9 <- performance
##

######################
## Read NMF BA4 in final model performance -- All cells. 
######################
## code
    # Read in file
    performance <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/NMF_BA4_final_performance_cell_level_top.csv')
    
    # Only retain overall performance
    performance <- subset(performance, group == "Overall")

    # Retain columns of interest
    performance <- performance %>% dplyr::select(input, celltype, mean_accuracy)
    nrow(performance)

    performance_NMF_BA4 <- performance
##

######################
## Read NMF BA9 in final model performance -- All cells. 
######################
## code
    # Read in file
    performance <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/NMF_BA9_final_performance_cell_level_top.csv')
    
    # Only retain overall performance
    performance <- subset(performance, group == "Overall")

    # Retain columns of interest
    performance <- performance %>% dplyr::select(input, celltype, mean_accuracy)
    nrow(performance)

    performance_NMF_BA9 <- performance
##

######################
## Merge performance
######################
## code 
    performance_HVG$type <- "All cells"
    performance_HVG$input_type <- "HVGs"

    performance_NMF$type <- "All cells"
    performance_NMF$input_type <- "NMF"

    performance_HVG_BA4$type <- "BA4"
    performance_HVG_BA4$input_type <- "HVGs"

    performance_HVG_BA9$type <- "BA9"
    performance_HVG_BA9$input_type <- "HVGs"

    performance_NMF_BA4$type <- "BA4"
    performance_NMF_BA4$input_type <- "NMF"

    performance_NMF_BA9$type <- "BA9"
    performance_NMF_BA9$input_type <- "NMF"

    performance_total <- rbind(performance_HVG, performance_NMF, performance_HVG_BA4, performance_HVG_BA9, performance_NMF_BA4, performance_NMF_BA9)
    nrow(performance_total)
    performance_total$celltype[performance_total$celltype == "Rose"] <- "Rosehip"
##

######################
## Read in meta data info and QC
######################
## code both regions
    # Read in file
    meta_QC <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/metadata_QC_python.csv')
    meta_QC$type <- "All cells"  
##

## code BA4
    # Read in file
    meta_QC_BA4 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/metadata_QC_python_BA4.csv')
    meta_QC_BA4$type <- "BA4"
##

## code BA9
    # Read in file
    meta_QC_BA9 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/metadata_QC_python_BA9.csv')
    meta_QC_BA9$type <- "BA9"
##

## bind and merge
    meta_QC_total <- rbind(meta_QC, meta_QC_BA4, meta_QC_BA9)

    # merge with model performance dataframe
    performance_total <- merge(performance_total, meta_QC_total, by = c("celltype", "type"), all.x = T) 
    nrow(performance_total)
    head(performance_total)
##

######################
## Read TxD -- All cells
######################
## code All cells diagnosis
    # Read in file
    TxD_1 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/TxD_ALS_FTLD_control_python.csv')
    
    # Retain columns of interest
    TxD_1 <- TxD_1 %>% dplyr::select(CellType, Total_divergence)

    # Remove duplicated rows
    TxD_1 <- unique(TxD_1)

    # Add input info
    TxD_1$input <- "Diagnosis"

    # Add type info
    TxD_1$type <- "All cells"


    # change column names
    colnames(TxD_1) <- c("celltype", "Total_divergence", "input", "type" )
##

## code All cells genetic
    # Read in file
    TxD_2 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/TxD_C9_sporadic_control_python.csv')
    
    # Retain columns of interest
    TxD_2 <- TxD_2 %>% dplyr::select(CellType, Total_divergence)

    # Remove duplicated rows
    TxD_2 <- unique(TxD_2)

    # Add input info
    TxD_2$input <- "Genetic group"

    # Add type info
    TxD_2$type <- "All cells"

    # change column names
    colnames(TxD_2) <- c("celltype", "Total_divergence", "input", "type" )
##

## code All cells disease
    # Read in file
    TxD_3 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/TxD_subtypes_control_python.csv')
    
    # Retain columns of interest
    TxD_3 <- TxD_3 %>% dplyr::select(CellType, Total_divergence)

    # Remove duplicated rows
    TxD_3 <- unique(TxD_3)

    # Add input info
    TxD_3$input <- "Disease group"

    # Add type info
    TxD_3$type <- "All cells"

    # change column names
    colnames(TxD_3) <- c("celltype", "Total_divergence", "input", "type" )
##

## code BA4 diagnosis
    # Read in file
    TxD_1_BA4 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/TxD_ALS_FTLD_control_python_BA4.csv')
    
    # Retain columns of interest
    TxD_1_BA4 <- TxD_1_BA4 %>% dplyr::select(CellType, Total_divergence)

    # Remove duplicated rows
    TxD_1_BA4 <- unique(TxD_1_BA4)

    # Add input info
    TxD_1_BA4$input <- "Diagnosis"

    # Add type info
    TxD_1_BA4$type <- "BA4"


    # change column names
    colnames(TxD_1_BA4) <- c("celltype", "Total_divergence", "input", "type" )
##

## code BA4 genetic
    # Read in file
    TxD_2_BA4 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/TxD_sporadic_C9_control_python_BA4.csv')
    
    # Retain columns of interest
    TxD_2_BA4 <- TxD_2_BA4 %>% dplyr::select(CellType, Total_divergence)

    # Remove duplicated rows
    TxD_2_BA4 <- unique(TxD_2_BA4)

    # Add input info
    TxD_2_BA4$input <- "Genetic group"

    # Add type info
    TxD_2_BA4$type <- "BA4"

    # change column names
    colnames(TxD_2_BA4) <- c("celltype", "Total_divergence", "input", "type" )
##

## code BA4 disease
    # Read in file
    TxD_3_BA4 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/TxD_subtypes_control_python_BA4.csv')
    
    # Retain columns of interest
    TxD_3_BA4 <- TxD_3_BA4 %>% dplyr::select(CellType, Total_divergence)

    # Remove duplicated rows
    TxD_3_BA4 <- unique(TxD_3_BA4)

    # Add input info
    TxD_3_BA4$input <- "Disease group"

    # Add type info
    TxD_3_BA4$type <- "BA4"

    # change column names
    colnames(TxD_3_BA4) <- c("celltype", "Total_divergence", "input", "type" )
##

## code BA9 diagnosis
    # Read in file
    TxD_1_BA9 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/TxD_ALS_FTLD_control_python_BA9.csv')
    
    # Retain columns of interest
    TxD_1_BA9 <- TxD_1_BA9 %>% dplyr::select(CellType, Total_divergence)

    # Remove duplicated rows
    TxD_1_BA9 <- unique(TxD_1_BA9)

    # Add input info
    TxD_1_BA9$input <- "Diagnosis"

    # Add type info
    TxD_1_BA9$type <- "BA9"


    # change column names
    colnames(TxD_1_BA9) <- c("celltype", "Total_divergence", "input", "type" )
##

## code BA9 genetic
    # Read in file
    TxD_2_BA9 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/TxD_sporadic_C9_control_python_BA9.csv')
    
    # Retain columns of interest
    TxD_2_BA9 <- TxD_2_BA9 %>% dplyr::select(CellType, Total_divergence)

    # Remove duplicated rows
    TxD_2_BA9 <- unique(TxD_2_BA9)

    # Add input info
    TxD_2_BA9$input <- "Genetic group"

    # Add type info
    TxD_2_BA9$type <- "BA9"

    # change column names
    colnames(TxD_2_BA9) <- c("celltype", "Total_divergence", "input", "type" )
##

## code BA9 disease
    # Read in file
    TxD_3_BA9 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/TxD_subtypes_control_python_BA9.csv')
    
    # Retain columns of interest
    TxD_3_BA9 <- TxD_3_BA9 %>% dplyr::select(CellType, Total_divergence)

    # Remove duplicated rows
    TxD_3_BA9 <- unique(TxD_3_BA9)

    # Add input info
    TxD_3_BA9$input <- "Disease group"

    # Add type info
    TxD_3_BA9$type <- "BA9"

    # change column names
    colnames(TxD_3_BA9) <- c("celltype", "Total_divergence", "input", "type" )
##

## Code: bind and merge
    # bind
    TxD_all <- rbind(TxD_1, TxD_2, TxD_3,
                    TxD_1_BA4, TxD_2_BA4, TxD_3_BA4,
                    TxD_1_BA9, TxD_2_BA9, TxD_3_BA9)
    
    # merge with model performance dataframe
    performance_total <- merge(performance_total, TxD_all, by = c("celltype", "input", "type"), all.x = T) 
    nrow(performance_total)
    head(performance_total)
    


##

######################
## Read class imbalance
######################
## code all cells diagnosis
    # Read in file
    imbalance_1 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/class_imbalance_ALS_FTLD_control_python.csv')
    
    # Retain columns of interest
    imbalance_1 <- imbalance_1 %>% dplyr::select(celltype, norm_shannon_entropy)

    # Add input info
    imbalance_1$input <- "Diagnosis"

    # Add type info
    imbalance_1$type <- "All cells"

    # change column names
    colnames(imbalance_1) <- c("celltype", "norm_shannon_entropy", "input", "type" )
##

## code all cells genetic
    # Read in file
    imbalance_2 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/class_imbalance_sporadic_C9_control_python.csv')
    
    # Retain columns of interest
    imbalance_2 <- imbalance_2 %>% dplyr::select(celltype, norm_shannon_entropy)

    # Add input info
    imbalance_2$input <- "Genetic group"

    # Add type info
    imbalance_2$type <- "All cells"

    # change column names
    colnames(imbalance_2) <- c("celltype", "norm_shannon_entropy", "input", "type" )
##

## code all cells disease
    # Read in file
    imbalance_3 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/class_imbalance_subtypes_control_python.csv')
    
    # Retain columns of interest
    imbalance_3 <- imbalance_3 %>% dplyr::select(celltype, norm_shannon_entropy)

    # Add input info
    imbalance_3$input <- "Disease group"

    # Add type info
    imbalance_3$type <- "All cells"

    # change column names
    colnames(imbalance_3) <- c("celltype", "norm_shannon_entropy", "input", "type" )
##

## code BA4 diagnosis
    # Read in file
    imbalance_1_BA4 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/class_imbalance_ALS_FTLD_control_python_BA4.csv')
    
    # Retain columns of interest
    imbalance_1_BA4 <- imbalance_1_BA4 %>% dplyr::select(celltype, norm_shannon_entropy)

    # Add input info
    imbalance_1_BA4$input <- "Diagnosis"

    # Add type info
    imbalance_1_BA4$type <- "BA4"

    # change column names
    colnames(imbalance_1_BA4) <- c("celltype", "norm_shannon_entropy", "input", "type" )
##

## code BA4 genetic
    # Read in file
    imbalance_2_BA4 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/class_imbalance_sporadic_C9_control_python_BA4.csv')
    
    # Retain columns of interest
    imbalance_2_BA4 <- imbalance_2_BA4 %>% dplyr::select(celltype, norm_shannon_entropy)

    # Add input info
    imbalance_2_BA4$input <- "Genetic group"

    # Add type info
    imbalance_2_BA4$type <- "BA4"

    # change column names
    colnames(imbalance_2_BA4) <- c("celltype", "norm_shannon_entropy", "input", "type" )
##

## code BA4 disease
    # Read in file
    imbalance_3_BA4 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/class_imbalance_subtypes_control_python_BA4.csv')
    
    # Retain columns of interest
    imbalance_3_BA4 <- imbalance_3_BA4 %>% dplyr::select(celltype, norm_shannon_entropy)

    # Add input info
    imbalance_3_BA4$input <- "Disease group"

    # Add type info
    imbalance_3_BA4$type <- "BA4"

    # change column names
    colnames(imbalance_3_BA4) <- c("celltype", "norm_shannon_entropy", "input", "type" )
##

## code BA9 diagnosis
    # Read in file
    imbalance_1_BA9 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/class_imbalance_ALS_FTLD_control_python_BA9.csv')
    
    # Retain columns of interest
    imbalance_1_BA9 <- imbalance_1_BA9 %>% dplyr::select(celltype, norm_shannon_entropy)

    # Add input info
    imbalance_1_BA9$input <- "Diagnosis"

    # Add type info
    imbalance_1_BA9$type <- "BA9"

    # change column names
    colnames(imbalance_1_BA9) <- c("celltype", "norm_shannon_entropy", "input", "type" )
##

## code BA9 genetic
    # Read in file
    imbalance_2_BA9 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/class_imbalance_sporadic_C9_control_python_BA9.csv')
    
    # Retain columns of interest
    imbalance_2_BA9 <- imbalance_2_BA9 %>% dplyr::select(celltype, norm_shannon_entropy)

    # Add input info
    imbalance_2_BA9$input <- "Genetic group"

    # Add type info
    imbalance_2_BA9$type <- "BA9"

    # change column names
    colnames(imbalance_2_BA9) <- c("celltype", "norm_shannon_entropy", "input", "type" )
##

## code BA9 disease
    # Read in file
    imbalance_3_BA9 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/class_imbalance_subtypes_control_python_BA9.csv')
    
    # Retain columns of interest
    imbalance_3_BA9 <- imbalance_3_BA9 %>% dplyr::select(celltype, norm_shannon_entropy)

    # Add input info
    imbalance_3_BA9$input <- "Disease group"

    # Add type info
    imbalance_3_BA9$type <- "BA9"

    # change column names
    colnames(imbalance_3_BA9) <- c("celltype", "norm_shannon_entropy", "input", "type" )
##

## Code: bind and merge
    # bind
    imbalance_all <- rbind(imbalance_1, imbalance_2, imbalance_3,
                    imbalance_1_BA4, imbalance_2_BA4, imbalance_3_BA4,
                    imbalance_1_BA9, imbalance_2_BA9, imbalance_3_BA9)
    
    # merge with model performance dataframe
    performance_total <- merge(performance_total, imbalance_all, by = c("celltype", "input", "type"), all.x = T) 
    nrow(performance_total)
    head(performance_total)
##


######################
## Drop-out-sparsity
######################
## code all cells 
    # Read in file
    drop_out <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/Drop_out_sparsity_python.csv')

    # Add type info
    drop_out$type <- "All cells"
##

## code BA4
    # Read in file
    drop_out_BA4 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/Drop_out_sparsity_python_BA4.csv')

    # Add type info
    drop_out_BA4$type <- "BA4"
##

## code BA9
    # Read in file
    drop_out_BA9 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/Drop_out_sparsity_python_BA9.csv')

    # Add type info
    drop_out_BA9$type <- "BA9"
##

## Code: bind and merge
    # bind
    drop_out_all <- rbind(drop_out, drop_out_BA4, drop_out_BA9)
    
    # merge with model performance dataframe
    performance_total <- merge(performance_total, drop_out_all, by = c("celltype", "type"), all.x = T) 
    nrow(performance_total)
    head(performance_total)
##


######################
## Signal-to_noise ratio
######################
## code All cells diagnosis
    # Read in file
    SNR_1 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/SNR_diagnosis_python.csv')
    
    # Retain columns of interest
    SNR_1 <- SNR_1 %>% dplyr::select(celltype, global_snr)

    # Add input info
    SNR_1$input <- "Diagnosis"

    # Add type info
    SNR_1$type <- "All cells"

    # change column names
    colnames(SNR_1) <- c("celltype", "global_snr", "input", "type" )

##

## code All cells genetic
    # Read in file
    SNR_2 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/SNR_genetic_python.csv')
    
    # Retain columns of interest
    SNR_2 <- SNR_2 %>% dplyr::select(celltype, global_snr)

    # Add input info
    SNR_2$input <- "Genetic group"

    # Add type info
    SNR_2$type <- "All cells"

    # change column names
    colnames(SNR_2) <- c("celltype", "global_snr", "input", "type" )

##

## code All cells disease
    # Read in file
    SNR_3 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/SNR_disease_python.csv')
    
    # Retain columns of interest
    SNR_3 <- SNR_3 %>% dplyr::select(celltype, global_snr)

    # Add input info
    SNR_3$input <- "Disease group"

    # change column names
    colnames(SNR_3) <- c("celltype", "global_snr", "input" )

    # Add type info
    SNR_3$type <- "All cells"

    # change column names
    colnames(SNR_3) <- c("celltype", "global_snr", "input", "type" )
##

## code BA4 diagnosis
    # Read in file
    SNR_1_BA4 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/SNR_diagnosis_python_BA4.csv')
    
    # Retain columns of interest
    SNR_1_BA4 <- SNR_1_BA4 %>% dplyr::select(celltype, global_snr)

    # Add input info
    SNR_1_BA4$input <- "Diagnosis"

    # Add type info
    SNR_1_BA4$type <- "BA4"

    # change column names
    colnames(SNR_1_BA4) <- c("celltype", "global_snr", "input", "type" )

##

## code BA4 genetic
    # Read in file
    SNR_2_BA4 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/SNR_genetic_python_BA4.csv')
    
    # Retain columns of interest
    SNR_2_BA4 <- SNR_2_BA4 %>% dplyr::select(celltype, global_snr)

    # Add input info
    SNR_2_BA4$input <- "Genetic group"

    # Add type info
    SNR_2_BA4$type <- "BA4"

    # change column names
    colnames(SNR_2_BA4) <- c("celltype", "global_snr", "input", "type" )

##

## code BA4 disease
    # Read in file
    SNR_3_BA4 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/SNR_disease_python_BA4.csv')
    
    # Retain columns of interest
    SNR_3_BA4 <- SNR_3_BA4 %>% dplyr::select(celltype, global_snr)

    # Add input info
    SNR_3_BA4$input <- "Disease group"

    # change column names
    colnames(SNR_3_BA4) <- c("celltype", "global_snr", "input" )

    # Add type info
    SNR_3_BA4$type <- "BA4"

    # change column names
    colnames(SNR_3_BA4) <- c("celltype", "global_snr", "input", "type" )
##

## code BA9 diagnosis
    # Read in file
    SNR_1_BA9 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/SNR_diagnosis_python_BA9.csv')
    
    # Retain columns of interest
    SNR_1_BA9 <- SNR_1_BA9 %>% dplyr::select(celltype, global_snr)

    # Add input info
    SNR_1_BA9$input <- "Diagnosis"

    # Add type info
    SNR_1_BA9$type <- "BA9"

    # change column names
    colnames(SNR_1_BA9) <- c("celltype", "global_snr", "input", "type" )

##

## code BA9 genetic
    # Read in file
    SNR_2_BA9 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/SNR_genetic_python_BA9.csv')
    
    # Retain columns of interest
    SNR_2_BA9 <- SNR_2_BA9 %>% dplyr::select(celltype, global_snr)

    # Add input info
    SNR_2_BA9$input <- "Genetic group"

    # Add type info
    SNR_2_BA9$type <- "BA9"

    # change column names
    colnames(SNR_2_BA9) <- c("celltype", "global_snr", "input", "type" )

##

## code BA9 disease
    # Read in file
    SNR_3_BA9 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/SNR_disease_python_BA9.csv')
    
    # Retain columns of interest
    SNR_3_BA9 <- SNR_3_BA9 %>% dplyr::select(celltype, global_snr)

    # Add input info
    SNR_3_BA9$input <- "Disease group"

    # change column names
    colnames(SNR_3_BA9) <- c("celltype", "global_snr", "input" )

    # Add type info
    SNR_3_BA9$type <- "BA9"

    # change column names
    colnames(SNR_3_BA9) <- c("celltype", "global_snr", "input", "type" )
##

## Code: bind and merge
    # bind
    SNR_all <- rbind(SNR_1, SNR_2, SNR_3,
                    SNR_1_BA4, SNR_2_BA4, SNR_3_BA4,
                    SNR_1_BA9, SNR_2_BA9, SNR_3_BA9)
    
    # merge with model performance dataframe
    performance_total <- merge(performance_total, SNR_all, by = c("celltype", "input", "type"), all.x = T) 
    nrow(performance_total)
    head(performance_total)
##

######################
## Donor-level variability v. within-donor variance: ICC
######################
## code all cells
    # Read in file
    ICC <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/ICC_python.csv')

    # Add type info
    ICC$type <- "All cells"
    
##

## code BA4
    # Read in file
    ICC_BA4 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/ICC_python_BA4.csv')

    # Add type info
    ICC_BA4$type <- "BA4"
    
##

## code BA9
    # Read in file
    ICC_BA9 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/ICC_python_BA9.csv')

    # Add type info
    ICC_BA9$type <- "BA9"
    
##

## Code: bind and merge
    # bind
    ICC_total <- rbind(ICC, ICC_BA4, ICC_BA9)
    
    # merge with model performance dataframe
    performance_total <- merge(performance_total, ICC_total, by = c("celltype", "type"), all.x = T) 
    nrow(performance_total)
    head(performance_total)
##

######################
## Cell type subtype
######################
## code all cells
    # Read in file
    cell_subtype <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/Cell_type_subtype.csv')

    # Add type info
    cell_subtype$type <- "All cells"
    
##

## code BA4
    # Read in file
    cell_subtype_BA4 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/Cell_type_subtype_BA4.csv')

    # Add type info
    cell_subtype_BA4$type <- "BA4"
    
##

## code BA9
    # Read in file
    cell_subtype_BA9 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/Cell_type_subtype_BA9.csv')

    # Add type info
    cell_subtype_BA9$type <- "BA9"
    
##


## Code: bind and merge
    # bind
    cell_subtype_all <- rbind(cell_subtype, cell_subtype_BA4, cell_subtype_BA9)
    
    # merge with model performance dataframe
    performance_total <- merge(performance_total, cell_subtype_all, by = c("celltype", "type"), all.x = T) 
    nrow(performance_total)
    head(performance_total)
##


head(performance_total)


######################
## Compute correlation and plot
######################
## code

    ## Grouped by model type
    cor_results_long_grouped <- performance_total %>%
    group_by(input, type, input_type) %>%
    summarise(
    across(c(n_HVGs:AdjustedDivergence),
            ~ {
                test <- cor.test(.x, mean_accuracy, method = "pearson")
                tibble(cor = unname(test$estimate), pval = test$p.value)
            },
            .names = "{.col}")
    ) 

    cor_results_long_grouped_flat <- cor_results_long_grouped %>%
    unnest_wider(n_HVGs, names_sep = "_") %>%
    unnest_wider(n_cells, names_sep = "_") %>%
    unnest_wider(mean_nCount, names_sep = "_") %>%
    unnest_wider(mean_nFeature, names_sep = "_") %>%
    unnest_wider(mean_percent_mt, names_sep = "_") %>%
    unnest_wider(Total_divergence, names_sep = "_") %>%
    unnest_wider(norm_shannon_entropy, names_sep = "_") %>%
    unnest_wider(sparsity, names_sep = "_") %>%
    unnest_wider(global_snr, names_sep = "_") %>%
    unnest_wider(global_icc, names_sep = "_") %>%
    unnest_wider(NumSubtypes, names_sep = "_") %>%
    unnest_wider(ShannonEntropy, names_sep = "_") %>%
    unnest_wider(AvgSubtypeDivergence, names_sep = "_") %>%
    unnest_wider(EffectiveSubtypes, names_sep = "_") %>%
    unnest_wider(AdjustedDivergence, names_sep = "_")

    cor_results_long_grouped_flat <- data.frame(cor_results_long_grouped_flat)


    df_long <- cor_results_long_grouped_flat %>%
    pivot_longer(
        cols = matches("(_cor|_pval)$"),   # all correlation + pval columns
        names_to = c("metric", ".value"),  # split names into metric + value
        names_pattern = "(.*)_(cor|pval)"  # regex captures: metric + suffix
    )

    df_long$sig <- NA
    df_long$sig[df_long$pval <= 0.05] <- "*"
    df_long$sig[df_long$pval <= 0.01] <- "**"
    df_long$sig[df_long$pval <= 0.001] <- "***"


    df_long <- data.frame(df_long)

    df_long_specific <- df_long


    ##############################################################################################################################

    performance_total$n_HVGs[performance_total$input_type == "NMF" ] <- 100
    performance_total$n_HVGs <- as.numeric(performance_total$n_HVGs)

    ## Not Grouped by model type
    cor_results_long_grouped <- performance_total %>%
    #group_by(type) %>%
    summarise(
    across(c(n_HVGs:AdjustedDivergence),
            ~ {
                test <- cor.test(.x, mean_accuracy)
                tibble(cor = unname(test$estimate), pval = test$p.value)
            },
            .names = "{.col}")
    ) 

    cor_results_long_grouped_flat <- cor_results_long_grouped %>%
    unnest_wider(n_HVGs, names_sep = "_") %>%
    unnest_wider(n_cells, names_sep = "_") %>%
    unnest_wider(mean_nCount, names_sep = "_") %>%
    unnest_wider(mean_nFeature, names_sep = "_") %>%
    unnest_wider(mean_percent_mt, names_sep = "_") %>%
    unnest_wider(Total_divergence, names_sep = "_") %>%
    unnest_wider(norm_shannon_entropy, names_sep = "_") %>%
    unnest_wider(sparsity, names_sep = "_") %>%
    unnest_wider(global_snr, names_sep = "_") %>%
    unnest_wider(global_icc, names_sep = "_") %>%
    unnest_wider(NumSubtypes, names_sep = "_") %>%
    unnest_wider(ShannonEntropy, names_sep = "_") %>%
    unnest_wider(AvgSubtypeDivergence, names_sep = "_") %>%
    unnest_wider(EffectiveSubtypes, names_sep = "_") %>%
    unnest_wider(AdjustedDivergence, names_sep = "_")

    cor_results_long_grouped_flat <- data.frame(cor_results_long_grouped_flat)


    df_long <- cor_results_long_grouped_flat %>%
    pivot_longer(
        cols = matches("(_cor|_pval)$"),   # all correlation + pval columns
        names_to = c("metric", ".value"),  # split names into metric + value
        names_pattern = "(.*)_(cor|pval)"  # regex captures: metric + suffix
    )

    df_long$sig <- NA
    df_long$sig[df_long$pval <= 0.05] <- "*"
    df_long$sig[df_long$pval <= 0.01] <- "**"
    df_long$sig[df_long$pval <= 0.001] <- "***"

    df_long <- data.frame(df_long)

    df_long_all <- df_long
    df_long_all$input <- "total"
    df_long_all$type <- "total"
    df_long_all$input_type <- "total"

    df_long_all <- df_long_all %>% dplyr::select(colnames(df_long_specific))

    #####################
    ## Bind and plot
    #####################
    df_long_specific$cor[df_long_specific$input_type == "NMF" & df_long_specific$metric == "n_HVGs" ] <- NA

    df_long_specific$type[df_long_specific$type == "BA4"] <- "MCx"
    df_long_specific$type[df_long_specific$type == "BA9"] <- "FCx"

    df_long_total <- rbind(df_long_specific, df_long_all)

    ## Fix factor levels
    df_long_total$input[df_long_total$input == "total"] <- "All models"
    df_long_total$input <- factor(df_long_total$input, levels = c("Diagnosis", "Genetic group", "Disease group", "All models" ) )

    df_long_total$type <- factor(df_long_total$type, levels = c("All cells", "MCx", "FCx") )


    ## remove intermediate metrics that we dont want
    df_long_total <- subset(df_long_total, metric != "AdjustedDivergence")
    df_long_total <- subset(df_long_total, metric != "EffectiveSubtypes")
    df_long_total <- subset(df_long_total, metric != "NumSubtypes")
    df_long_total <- subset(df_long_total, metric != "ShannonEntropy")
    df_long_total <- subset(df_long_total, metric != "global_snr")

    #df_long_total$metric_clarified[df_long_total$metric == "n_HVGs"] <- "n HVGs used as input"
    #df_long_total$metric_clarified[df_long_total$metric == "n_cells"] <- "n cells used used for LOSO test-train split"
    #df_long_total$metric_clarified[df_long_total$metric == "mean_nCount"] <- "Total transcript counts per cell"
    #df_long_total$metric_clarified[df_long_total$metric == "mean_nFeature"] <- "Number of unique genes expressed per cell"
    #df_long_total$metric_clarified[df_long_total$metric == "mean_percent_mt"] <- "Percentage of mitochondrial-encoded genes"
    #df_long_total$metric_clarified[df_long_total$metric == "sparsity"] <- "Percentage of genes with zero counts (dropout rate)"
    #df_long_total$metric_clarified[df_long_total$metric == "norm_shannon_entropy"] <- "Class balance (normalized Shannon entropy)"
    #df_long_total$metric_clarified[df_long_total$metric == "AvgSubtypeDivergence"] <- "Transcriptional divergence between cell type-specific subtypes"
    #df_long_total$metric_clarified[df_long_total$metric == "global_snr"] <- "Signal-to-noise ratio (ratio of between-class to within-class variance)"
    #df_long_total$metric_clarified[df_long_total$metric == "global_icc"] <- "Transcriptional consistency within subjects (Intraclass Correlation Coefficient)"
    #df_long_total$metric_clarified[df_long_total$metric == "Total_divergence"] <- "Transcriptional divergence between subject classes"

    df_long_total$metric_clarified[df_long_total$metric == "n_HVGs"] <- "n HVGs used as input to DNN"
    df_long_total$metric_clarified[df_long_total$metric == "n_cells"] <- "n cells used used for LOSO test-train split"
    df_long_total$metric_clarified[df_long_total$metric == "mean_nCount"] <- "Total transcript counts per cell"
    df_long_total$metric_clarified[df_long_total$metric == "mean_nFeature"] <- "Number of unique genes expressed per cell"
    df_long_total$metric_clarified[df_long_total$metric == "mean_percent_mt"] <- "Percentage of mitochondrial-encoded genes"
    df_long_total$metric_clarified[df_long_total$metric == "sparsity"] <- "Percentage of genes with zero counts (dropout rate)"
    df_long_total$metric_clarified[df_long_total$metric == "norm_shannon_entropy"] <- "Class balance (normalized Shannon entropy)"
    df_long_total$metric_clarified[df_long_total$metric == "AvgSubtypeDivergence"] <- "Transcriptional divergence between cellular subtypes"
    df_long_total$metric_clarified[df_long_total$metric == "global_icc"] <- "Transcriptional consistency within subjects"
    df_long_total$metric_clarified[df_long_total$metric == "Total_divergence"] <- "Transcriptional divergence between subject classes"


    df_long_total$metric_clarified <- factor(df_long_total$metric_clarified, 
    levels = rev(c("n HVGs used as input to DNN",
                "n cells used used for LOSO test-train split",
                "Total transcript counts per cell",
                "Number of unique genes expressed per cell",
                "Percentage of mitochondrial-encoded genes",
                "Percentage of genes with zero counts (dropout rate)",
                "Class balance (normalized Shannon entropy)",
                "Transcriptional divergence between cellular subtypes",
                "Transcriptional consistency within subjects",
                "Transcriptional divergence between subject classes")))

    write.csv(df_long_total, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/reasons_correlation_total_dataframe.csv' )


    ## Plot
    ggplot(df_long_total, aes(x = input, y = metric_clarified, fill = cor)) +
            theme_bw() + 
            geom_tile() +
            geom_text(aes(label = sig), size = 3, colour = "black") +
            theme(
                legend.position = "none",
                panel.grid = element_blank(),
                axis.text.x = element_text(
                colour = "black",
                size = 8,
                angle = 90,
                hjust = 1,   # aligns the end of text with tick
                vjust = 0.5    # pushes text closer to axis
                ),
                #axis.ticks.x = element_blank(),
                axis.text.y = element_text(colour = "black", size = 8),
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                strip.background = element_rect(fill="lightgrey", colour = "white"),
                strip.text = element_text(size= 6.5, face="bold", colour = "black"),
                strip.placement = "outside"
            ) +
            #facet_nested(. ~ type, scales = "free", space = "free") +
            facet_nested(.~ input_type + type, scales = "free_x", space = "free_x") +
            ylab("Accuracy") +
            #scale_size_continuous(range = c(1, 5))+
            scale_x_discrete(expand = c(0,0)) +
            scale_y_discrete(expand = c(0,0)) +
            scale_fill_gradientn(colors = rev(brewer.pal(11, "RdBu")), limits = c(-1, 1), na.value = "darkgrey" )

            ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 3, width = 6.5)
##



