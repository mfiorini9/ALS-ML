#run in Narval
salloc -A def-sfarhan --time=0-6 -c 1 --mem=40g

module load StdEnv/2023
module load r/4.4.0
R


## load libraries
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
library(labeling, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")



############################################################################ Print files containing genes with Z-score > 1 --> this will go into the permutation scripts. --> WE ARE NOT USING THESE BECAUSE WE RUN THIS SAME CODE IN THE PERMUTATION PYTHON SCRIPT NOW
############################################################################
############################################################################
############################################################################ 

## code
    ## code SALS BA4
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA4"
        par_status = "SALS"
        par_prep = "CombatSeq"
        
        ################################################
        ## LIME outputs not weighted by percent expression
        ## This function subsets the unweighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_across_donors = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_across_donors, celltype)

        for(celltype2 in unique(celltype_list)){
            

            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_unweighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval_2.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_across_donors <- as.numeric(df$Mean_z_score_across_donors)
            df <- df[order( -(df$Mean_z_score_across_donors)),]

            df <- subset(df, Mean_z_score_across_donors >= 1)
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_across_donors)
            df$celltype <- celltype2

            write.csv(df, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_unweighted_Z1_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval_2.csv'), sep = ",")

            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")

        BA4_SALS_unweighted_df <- fill

        ################################################
        ## LIME outputs weighted by percent expression
        ## This function subsets the weighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_weighted_across_donor = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_weighted_across_donor, celltype)

        for(celltype2 in unique(celltype_list)){
            
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval_2.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_weighted_across_donor <- as.numeric(df$Mean_z_score_weighted_across_donor)
            df <- df[order( -(df$Mean_z_score_weighted_across_donor)),]

            df <- subset(df, Mean_z_score_weighted_across_donor >= 1)
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_weighted_across_donor)
            df$celltype <- celltype2

            write.csv(df, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_Z1_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval_2.csv'), sep = ",")

            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")
        BA4_SALS_weighted_df <- fill

    ##

    ## code C9ALS BA4
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA4"
        par_status = "C9ALS"
        par_prep = "CombatSeq"
        
        ################################################
        ## LIME outputs not weighted by percent expression
        ## This function subsets the unweighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_across_donors = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_across_donors, celltype)

        for(celltype2 in unique(celltype_list)){
            

            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_unweighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval_2.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_across_donors <- as.numeric(df$Mean_z_score_across_donors)
            df <- df[order( -(df$Mean_z_score_across_donors)),]

            df <- subset(df, Mean_z_score_across_donors >= 1)
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_across_donors)
            df$celltype <- celltype2

            write.csv(df, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_unweighted_Z1_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval_2.csv'), sep = ",")

            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")

        BA4_C9ALS_unweighted_df <- fill

        ################################################
        ## LIME outputs weighted by percent expression
        ## This function subsets the weighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_weighted_across_donor = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_weighted_across_donor, celltype)

        for(celltype2 in unique(celltype_list)){
            
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval_2.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_weighted_across_donor <- as.numeric(df$Mean_z_score_weighted_across_donor)
            df <- df[order( -(df$Mean_z_score_weighted_across_donor)),]

            df <- subset(df, Mean_z_score_weighted_across_donor >= 1)
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_weighted_across_donor)
            df$celltype <- celltype2

            write.csv(df, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_Z1_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval_2.csv'), sep = ",")

            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")
        BA4_C9ALS_weighted_df <- fill

    ##

    ## code SFTLD BA4
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA4"
        par_status = "SFTLD"
        par_prep = "CombatSeq"
        
        ################################################
        ## LIME outputs not weighted by percent expression
        ## This function subsets the unweighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_across_donors = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_across_donors, celltype)

        for(celltype2 in unique(celltype_list)){
            

            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_unweighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval_2.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_across_donors <- as.numeric(df$Mean_z_score_across_donors)
            df <- df[order( -(df$Mean_z_score_across_donors)),]

            df <- subset(df, Mean_z_score_across_donors >= 1)
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_across_donors)
            df$celltype <- celltype2

            write.csv(df, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_unweighted_Z1_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval_2.csv'), sep = ",")

            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")

        BA4_SFTLD_unweighted_df <- fill

        ################################################
        ## LIME outputs weighted by percent expression
        ## This function subsets the weighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_weighted_across_donor = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_weighted_across_donor, celltype)

        for(celltype2 in unique(celltype_list)){
            
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval_2.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_weighted_across_donor <- as.numeric(df$Mean_z_score_weighted_across_donor)
            df <- df[order( -(df$Mean_z_score_weighted_across_donor)),]

            df <- subset(df, Mean_z_score_weighted_across_donor >= 1)
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_weighted_across_donor)
            df$celltype <- celltype2

            write.csv(df, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_Z1_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval_2.csv'), sep = ",")

            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")
        BA4_SFTLD_weighted_df <- fill

    ##

    ## code C9FTLD BA4
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA4"
        par_status = "C9FTLD"
        par_prep = "CombatSeq"
        
        ################################################
        ## LIME outputs not weighted by percent expression
        ## This function subsets the unweighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_across_donors = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_across_donors, celltype)

        for(celltype2 in unique(celltype_list)){
            

            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_unweighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval_2.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_across_donors <- as.numeric(df$Mean_z_score_across_donors)
            df <- df[order( -(df$Mean_z_score_across_donors)),]

            df <- subset(df, Mean_z_score_across_donors >= 1)
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_across_donors)
            df$celltype <- celltype2

            write.csv(df, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_unweighted_Z1_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval_2.csv'), sep = ",")

            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")

        BA4_C9FTLD_unweighted_df <- fill

        ################################################
        ## LIME outputs weighted by percent expression
        ## This function subsets the weighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_weighted_across_donor = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_weighted_across_donor, celltype)

        for(celltype2 in unique(celltype_list)){
            
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval_2.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_weighted_across_donor <- as.numeric(df$Mean_z_score_weighted_across_donor)
            df <- df[order( -(df$Mean_z_score_weighted_across_donor)),]

            df <- subset(df, Mean_z_score_weighted_across_donor >= 1)
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_weighted_across_donor)
            df$celltype <- celltype2

            write.csv(df, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_Z1_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval_2.csv'), sep = ",")

            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")
        BA4_C9FTLD_weighted_df <- fill

    ##

    ## code SALS BA9
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA9"
        par_status = "SALS"
        par_prep = "CombatSeq"
        
        ################################################
        ## LIME outputs not weighted by percent expression
        ## This function subsets the unweighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_across_donors = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_across_donors, celltype)

        for(celltype2 in unique(celltype_list)){
            

            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_unweighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval_2.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_across_donors <- as.numeric(df$Mean_z_score_across_donors)
            df <- df[order( -(df$Mean_z_score_across_donors)),]

            df <- subset(df, Mean_z_score_across_donors >= 1)
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_across_donors)
            df$celltype <- celltype2

            write.csv(df, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_unweighted_Z1_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval_2.csv'), sep = ",")

            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")

        BA9_SALS_unweighted_df <- fill

        ################################################
        ## LIME outputs weighted by percent expression
        ## This function subsets the weighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_weighted_across_donor = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_weighted_across_donor, celltype)

        for(celltype2 in unique(celltype_list)){
            
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval_2.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_weighted_across_donor <- as.numeric(df$Mean_z_score_weighted_across_donor)
            df <- df[order( -(df$Mean_z_score_weighted_across_donor)),]

            df <- subset(df, Mean_z_score_weighted_across_donor >= 1)
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_weighted_across_donor)
            df$celltype <- celltype2

            write.csv(df, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_Z1_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval_2.csv'), sep = ",")

            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")
        BA9_SALS_weighted_df <- fill

    ##

    ## code C9ALS BA9
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA9"
        par_status = "C9ALS"
        par_prep = "CombatSeq"
        
        ################################################
        ## LIME outputs not weighted by percent expression
        ## This function subsets the unweighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_across_donors = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_across_donors, celltype)

        for(celltype2 in unique(celltype_list)){
            

            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_unweighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval_2.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_across_donors <- as.numeric(df$Mean_z_score_across_donors)
            df <- df[order( -(df$Mean_z_score_across_donors)),]

            df <- subset(df, Mean_z_score_across_donors >= 1)
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_across_donors)
            df$celltype <- celltype2

            write.csv(df, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_unweighted_Z1_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval_2.csv'), sep = ",")

            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")

        BA9_C9ALS_unweighted_df <- fill

        ################################################
        ## LIME outputs weighted by percent expression
        ## This function subsets the weighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_weighted_across_donor = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_weighted_across_donor, celltype)

        for(celltype2 in unique(celltype_list)){
            
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval_2.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_weighted_across_donor <- as.numeric(df$Mean_z_score_weighted_across_donor)
            df <- df[order( -(df$Mean_z_score_weighted_across_donor)),]

            df <- subset(df, Mean_z_score_weighted_across_donor >= 1)
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_weighted_across_donor)
            df$celltype <- celltype2

            write.csv(df, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_Z1_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval_2.csv'), sep = ",")

            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")
        BA9_C9ALS_weighted_df <- fill

    ##

    ## code SFTLD BA9
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA9"
        par_status = "SFTLD"
        par_prep = "CombatSeq"
        
        ################################################
        ## LIME outputs not weighted by percent expression
        ## This function subsets the unweighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_across_donors = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_across_donors, celltype)

        for(celltype2 in unique(celltype_list)){
            

            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_unweighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval_2.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_across_donors <- as.numeric(df$Mean_z_score_across_donors)
            df <- df[order( -(df$Mean_z_score_across_donors)),]

            df <- subset(df, Mean_z_score_across_donors >= 1)
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_across_donors)
            df$celltype <- celltype2

            write.csv(df, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_unweighted_Z1_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval_2.csv'), sep = ",")

            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")

        BA9_SFTLD_unweighted_df <- fill

        ################################################
        ## LIME outputs weighted by percent expression
        ## This function subsets the weighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_weighted_across_donor = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_weighted_across_donor, celltype)

        for(celltype2 in unique(celltype_list)){
            
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval_2.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_weighted_across_donor <- as.numeric(df$Mean_z_score_weighted_across_donor)
            df <- df[order( -(df$Mean_z_score_weighted_across_donor)),]

            df <- subset(df, Mean_z_score_weighted_across_donor >= 1)
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_weighted_across_donor)
            df$celltype <- celltype2

            write.csv(df, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_Z1_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval_2.csv'), sep = ",")

            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")
        BA9_SFTLD_weighted_df <- fill

    ##

    ## code C9FTLD BA9 -- missing TCELL for some reason -- look into it
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA9"
        par_status = "C9FTLD"
        par_prep = "CombatSeq"
        
        ################################################
        ## LIME outputs not weighted by percent expression
        ## This function subsets the unweighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_across_donors = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_across_donors, celltype)

        for(celltype2 in unique(celltype_list)){
            

            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_unweighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval_2.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_across_donors <- as.numeric(df$Mean_z_score_across_donors)
            df <- df[order( -(df$Mean_z_score_across_donors)),]

            df <- subset(df, Mean_z_score_across_donors >= 1)
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_across_donors)
            df$celltype <- celltype2

            write.csv(df, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_unweighted_Z1_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval_2.csv'), sep = ",")

            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")

        BA9_C9FTLD_unweighted_df <- fill

        ################################################
        ## LIME outputs weighted by percent expression
        ## This function subsets the weighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_weighted_across_donor = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_weighted_across_donor, celltype)

        for(celltype2 in unique(celltype_list)){
            
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval_2.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_weighted_across_donor <- as.numeric(df$Mean_z_score_weighted_across_donor)
            df <- df[order( -(df$Mean_z_score_weighted_across_donor)),]

            df <- subset(df, Mean_z_score_weighted_across_donor >= 1)
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_weighted_across_donor)
            df$celltype <- celltype2

            write.csv(df, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_Z1_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval_2.csv'), sep = ",")

            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")
        BA9_C9FTLD_weighted_df <- fill

    ##
## 


## RUN PERMUTATION SCRIPTS.


############################################################################ Z-score distributions -- weighted by percent expression (potential Supplemental Figure)
############################################################################
############################################################################
############################################################################ 

## code
    ###########################
    ## SALS BA4
    ###########################
    ## code
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA4"
        par_status = "SALS"
        par_prep = "CombatSeq"
        
        ################################################
        ## LIME outputs weighted by percent expression
        ## This function subsets the weighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_weighted_across_donor = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_weighted_across_donor, celltype)

        
        for(celltype2 in unique(celltype_list)){
            
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",") ## CHECK IF THEY ARE THE SAME
       
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_weighted_across_donor <- as.numeric(df$Mean_z_score_weighted_across_donor)
            df <- df[order( -(df$Mean_z_score_weighted_across_donor)),]
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_weighted_across_donor)
            df$celltype <- celltype2
            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")
        BA4_SALS_weighted_df <- fill
    ##

    ###########################
    ## C9ALS BA4
    ###########################
    ## code
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA4"
        par_status = "C9ALS"
        par_prep = "CombatSeq"
        
        ################################################
        ## LIME outputs weighted by percent expression
        ## This function subsets the weighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_weighted_across_donor = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_weighted_across_donor, celltype)

        for(celltype2 in unique(celltype_list)){
            
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_weighted_across_donor <- as.numeric(df$Mean_z_score_weighted_across_donor)
            df <- df[order( -(df$Mean_z_score_weighted_across_donor)),]
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_weighted_across_donor)
            df$celltype <- celltype2
            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")
        BA4_C9ALS_weighted_df <- fill

    ##

    ###########################
    ## FTLD BA4
    ###########################
    ## code
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA4"
        par_status = "SFTLD"
        par_prep = "CombatSeq"
        
        ################################################
        ## LIME outputs weighted by percent expression
        ## This function subsets the weighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_weighted_across_donor = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_weighted_across_donor, celltype)

        for(celltype2 in unique(celltype_list)){
            
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_weighted_across_donor <- as.numeric(df$Mean_z_score_weighted_across_donor)
            df <- df[order( -(df$Mean_z_score_weighted_across_donor)),]
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_weighted_across_donor)
            df$celltype <- celltype2
            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")
        BA4_SFTLD_weighted_df <- fill

    ##
    
    ###########################
    ## C9FTLD BA4
    ###########################
    ## code
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA4"
        par_status = "C9FTLD"
        par_prep = "CombatSeq"

        ################################################
        ## LIME outputs weighted by percent expression
        ## This function subsets the weighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_weighted_across_donor = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_weighted_across_donor, celltype)

        for(celltype2 in unique(celltype_list)){
            
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_weighted_across_donor <- as.numeric(df$Mean_z_score_weighted_across_donor)
            df <- df[order( -(df$Mean_z_score_weighted_across_donor)),]
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_weighted_across_donor)
            df$celltype <- celltype2
            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")
        BA4_C9FTLD_weighted_df <- fill

    ##
    
    ###########################
    ## SALS BA9
    ###########################
    ## code 
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA9"
        par_status = "SALS"
        par_prep = "CombatSeq"
        
        ################################################
        ## LIME outputs weighted by percent expression
        ## This function subsets the weighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_weighted_across_donor = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_weighted_across_donor, celltype)

        for(celltype2 in unique(celltype_list)){
            
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_weighted_across_donor <- as.numeric(df$Mean_z_score_weighted_across_donor)
            df <- df[order( -(df$Mean_z_score_weighted_across_donor)),]
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_weighted_across_donor)
            df$celltype <- celltype2
            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")
        BA9_SALS_weighted_df <- fill

    ##
    
    ###########################
    ## C9ALS BA9
    ###########################
    ## code 
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA9"
        par_status = "C9ALS"
        par_prep = "CombatSeq"

        ################################################
        ## LIME outputs weighted by percent expression
        ## This function subsets the weighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_weighted_across_donor = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_weighted_across_donor, celltype)

        for(celltype2 in unique(celltype_list)){
            
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_weighted_across_donor <- as.numeric(df$Mean_z_score_weighted_across_donor)
            df <- df[order( -(df$Mean_z_score_weighted_across_donor)),]
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_weighted_across_donor)
            df$celltype <- celltype2
            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")
        BA9_C9ALS_weighted_df <- fill

    ##

    ###########################
    ## SFTLD BA9
    ###########################
    ## code 
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA9"
        par_status = "SFTLD"
        par_prep = "CombatSeq"
        
        ################################################
        ## LIME outputs weighted by percent expression
        ## This function subsets the weighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_weighted_across_donor = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_weighted_across_donor, celltype)

        for(celltype2 in unique(celltype_list)){
            
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_weighted_across_donor <- as.numeric(df$Mean_z_score_weighted_across_donor)
            df <- df[order( -(df$Mean_z_score_weighted_across_donor)),]
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_weighted_across_donor)
            df$celltype <- celltype2
            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")
        BA9_SFTLD_weighted_df <- fill

    ##

    ###########################
    ## C9FTLD BA9 -- missing TCELL for some reason -- look into it
    ###########################
    ## code 

        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA9"
        par_status = "C9FTLD"
        par_prep = "CombatSeq"
        
        ################################################
        ## LIME outputs weighted by percent expression
        ## This function subsets the weighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_weighted_across_donor = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_weighted_across_donor, celltype)

        for(celltype2 in unique(celltype_list)){
            
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_weighted_across_donor <- as.numeric(df$Mean_z_score_weighted_across_donor)
            df <- df[order( -(df$Mean_z_score_weighted_across_donor)),]
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_weighted_across_donor)
            df$celltype <- celltype2
            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")
        BA9_C9FTLD_weighted_df <- fill

    ##

    ###########################
    ## All conditions BA4
    ###########################
    ## code
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA4"
        par_status = "all_conditions"
        par_prep = "CombatSeq"

        ################################################
        ## LIME outputs weighted by percent expression
        ## This function subsets the weighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_weighted_across_donor = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_weighted_across_donor, celltype)

        for(celltype2 in unique(celltype_list)){
            
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_weighted_across_donor <- as.numeric(df$Mean_z_score_weighted_across_donor)
            df <- df[order( -(df$Mean_z_score_weighted_across_donor)),]
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_weighted_across_donor)
            df$celltype <- celltype2
            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")
        BA4_all_conditions_weighted_df <- fill

    ##

    ###########################
    ## All conditions BA9
    ###########################
    ## code
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA9"
        par_status = "all_conditions"
        par_prep = "CombatSeq"

        ################################################
        ## LIME outputs weighted by percent expression
        ## This function subsets the weighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_weighted_across_donor = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_weighted_across_donor, celltype)

        for(celltype2 in unique(celltype_list)){
            
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_weighted_across_donor <- as.numeric(df$Mean_z_score_weighted_across_donor)
            df <- df[order( -(df$Mean_z_score_weighted_across_donor)),]
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_weighted_across_donor)
            df$celltype <- celltype2
            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")
        BA9_all_conditions_weighted_df <- fill

    ##

    ###########################
    ## Plot BA4
    ###########################
    ## code
        BA4_SALS_weighted_df$disease_cat <- "SALS"
        BA4_C9ALS_weighted_df$disease_cat <- "C9ALS"
        BA4_SFTLD_weighted_df$disease_cat <- "SFTLD"
        BA4_C9FTLD_weighted_df$disease_cat <- "C9FTLD"
        BA4_all_conditions_weighted_df$disease_cat <- "All"

        BA4_weighted_distribution <- rbind(BA4_SALS_weighted_df, BA4_C9ALS_weighted_df, BA4_SFTLD_weighted_df, BA4_C9FTLD_weighted_df, BA4_all_conditions_weighted_df)

        BA4_weighted_distribution$disease_cat <- factor(BA4_weighted_distribution$disease_cat, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All"))
        BA4_weighted_distribution$celltype <- factor(BA4_weighted_distribution$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))


        ggplot(BA4_weighted_distribution, aes(x = Mean_z_score_weighted_across_donor, colour = disease_cat)) + 
        theme_bw() + 
        geom_vline(xintercept = 1.0, linetype = "dashed", colour = "black") +
        geom_density() +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(colour = "black", angle = 90, hjust = 1, vjust = 0.5),
            axis.text.y = element_text(colour = "black"),
            axis.title.y = element_text(face = "bold"),
            axis.title.x = element_text(face = "bold"),
        ) +
        facet_grid(disease_cat~celltype, scales = "free_x") +
        ylab("Density") + xlab("LIME feature importance Z-score") +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "black")) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 4, width = 15)
    ##

    ###########################
    ## Plot BA9
    ###########################
    ## code
        BA9_SALS_weighted_df$disease_cat <- "SALS"
        BA9_C9ALS_weighted_df$disease_cat <- "C9ALS"
        BA9_SFTLD_weighted_df$disease_cat <- "SFTLD"
        BA9_C9FTLD_weighted_df$disease_cat <- "C9FTLD"
        BA9_all_conditions_weighted_df$disease_cat <- "All"

        BA9_weighted_distribution <- rbind(BA9_SALS_weighted_df, BA9_C9ALS_weighted_df, BA9_SFTLD_weighted_df, BA9_C9FTLD_weighted_df, BA9_all_conditions_weighted_df)

        BA9_weighted_distribution$disease_cat <- factor(BA9_weighted_distribution$disease_cat, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All"))
        BA9_weighted_distribution$celltype <- factor(BA9_weighted_distribution$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))


        ggplot(BA9_weighted_distribution, aes(x = Mean_z_score_weighted_across_donor, colour = disease_cat)) + 
        theme_bw() + 
        geom_vline(xintercept = 1.0, linetype = "dashed", colour = "black") +
        geom_density() +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(colour = "black", angle = 90, hjust = 1, vjust = 0.5),
            axis.text.y = element_text(colour = "black"),
            axis.title.y = element_text(face = "bold"),
            axis.title.x = element_text(face = "bold"),
        ) +
        facet_grid(disease_cat~celltype, scales = "free_x") +
        ylab("Density") + xlab("LIME feature importance Z-score") +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "black")) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 4, width = 15)
    ##

##

############################################################################ N features with Z-score > 1 -- weighted by percent expression (potential Supplemental Figure)
############################################################################
############################################################################
############################################################################ 

## code
    ###########################
    ## SALS BA4
    ###########################
    ## code
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA4"
        par_status = "SALS"
        par_prep = "CombatSeq"
        
        ################################################
        ## LIME outputs weighted by percent expression
        ## This function subsets the weighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_weighted_across_donor = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_weighted_across_donor, celltype)

        
        for(celltype2 in unique(celltype_list)){
            
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_weighted_across_donor <- as.numeric(df$Mean_z_score_weighted_across_donor)
            df <- df[order( -(df$Mean_z_score_weighted_across_donor)),]
            df <- subset(df, Mean_z_score_weighted_across_donor >= 1)
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_weighted_across_donor)
            df$celltype <- celltype2
            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")
        BA4_SALS_weighted_df <- fill
    ##

    ###########################
    ## C9ALS BA4
    ###########################
    ## code
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA4"
        par_status = "C9ALS"
        par_prep = "CombatSeq"
        
        ################################################
        ## LIME outputs weighted by percent expression
        ## This function subsets the weighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_weighted_across_donor = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_weighted_across_donor, celltype)

        for(celltype2 in unique(celltype_list)){
            
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_weighted_across_donor <- as.numeric(df$Mean_z_score_weighted_across_donor)
            df <- df[order( -(df$Mean_z_score_weighted_across_donor)),]
            df <- subset(df, Mean_z_score_weighted_across_donor >= 1)
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_weighted_across_donor)
            df$celltype <- celltype2
            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")
        BA4_C9ALS_weighted_df <- fill

    ##

    ###########################
    ## SFTLD BA4
    ###########################
    ## code
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA4"
        par_status = "SFTLD"
        par_prep = "CombatSeq"
        
        ################################################
        ## LIME outputs weighted by percent expression
        ## This function subsets the weighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_weighted_across_donor = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_weighted_across_donor, celltype)

        for(celltype2 in unique(celltype_list)){
            
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_weighted_across_donor <- as.numeric(df$Mean_z_score_weighted_across_donor)
            df <- df[order( -(df$Mean_z_score_weighted_across_donor)),]
            df <- subset(df, Mean_z_score_weighted_across_donor >= 1)
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_weighted_across_donor)
            df$celltype <- celltype2
            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")
        BA4_SFTLD_weighted_df <- fill

    ##
    
    ###########################
    ## C9FTLD BA4
    ###########################
    ## code
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA4"
        par_status = "C9FTLD"
        par_prep = "CombatSeq"

        ################################################
        ## LIME outputs weighted by percent expression
        ## This function subsets the weighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_weighted_across_donor = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_weighted_across_donor, celltype)

        for(celltype2 in unique(celltype_list)){
            
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_weighted_across_donor <- as.numeric(df$Mean_z_score_weighted_across_donor)
            df <- df[order( -(df$Mean_z_score_weighted_across_donor)),]
            df <- subset(df, Mean_z_score_weighted_across_donor >= 1)
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_weighted_across_donor)
            df$celltype <- celltype2
            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")
        BA4_C9FTLD_weighted_df <- fill

    ##
    
    ###########################
    ## SALS BA9
    ###########################
    ## code 
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA9"
        par_status = "SALS"
        par_prep = "CombatSeq"
        
        ################################################
        ## LIME outputs weighted by percent expression
        ## This function subsets the weighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_weighted_across_donor = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_weighted_across_donor, celltype)

        for(celltype2 in unique(celltype_list)){
            
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_weighted_across_donor <- as.numeric(df$Mean_z_score_weighted_across_donor)
            df <- df[order( -(df$Mean_z_score_weighted_across_donor)),]
            df <- subset(df, Mean_z_score_weighted_across_donor >= 1)
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_weighted_across_donor)
            df$celltype <- celltype2
            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")
        BA9_SALS_weighted_df <- fill

    ##
    
    ###########################
    ## C9ALS BA9
    ###########################
    ## code 
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA9"
        par_status = "C9ALS"
        par_prep = "CombatSeq"

        ################################################
        ## LIME outputs weighted by percent expression
        ## This function subsets the weighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_weighted_across_donor = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_weighted_across_donor, celltype)

        for(celltype2 in unique(celltype_list)){
            
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_weighted_across_donor <- as.numeric(df$Mean_z_score_weighted_across_donor)
            df <- df[order( -(df$Mean_z_score_weighted_across_donor)),]
            df <- subset(df, Mean_z_score_weighted_across_donor >= 1)
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_weighted_across_donor)
            df$celltype <- celltype2
            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")
        BA9_C9ALS_weighted_df <- fill

    ##

    ###########################
    ## SFTLD BA9
    ###########################
    ## code 
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA9"
        par_status = "SFTLD"
        par_prep = "CombatSeq"
        
        ################################################
        ## LIME outputs weighted by percent expression
        ## This function subsets the weighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_weighted_across_donor = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_weighted_across_donor, celltype)

        for(celltype2 in unique(celltype_list)){
            
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_weighted_across_donor <- as.numeric(df$Mean_z_score_weighted_across_donor)
            df <- df[order( -(df$Mean_z_score_weighted_across_donor)),]
            df <- subset(df, Mean_z_score_weighted_across_donor >= 1)
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_weighted_across_donor)
            df$celltype <- celltype2
            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")
        BA9_SFTLD_weighted_df <- fill

    ##

    ###########################
    ## C9FTLD BA9 -- missing TCELL for some reason -- look into it
    ###########################
    ## code 

        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA9"
        par_status = "C9FTLD"
        par_prep = "CombatSeq"
        
        ################################################
        ## LIME outputs weighted by percent expression
        ## This function subsets the weighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_weighted_across_donor = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_weighted_across_donor, celltype)

        for(celltype2 in unique(celltype_list)){
            
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_weighted_across_donor <- as.numeric(df$Mean_z_score_weighted_across_donor)
            df <- df[order( -(df$Mean_z_score_weighted_across_donor)),]
            df <- subset(df, Mean_z_score_weighted_across_donor >= 1)
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_weighted_across_donor)
            df$celltype <- celltype2
            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")
        BA9_C9FTLD_weighted_df <- fill

    ##

    ###########################
    ## All conditions BA4
    ###########################
    ## code
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA4"
        par_status = "all_conditions"
        par_prep = "CombatSeq"

        ################################################
        ## LIME outputs weighted by percent expression
        ## This function subsets the weighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_weighted_across_donor = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_weighted_across_donor, celltype)

        for(celltype2 in unique(celltype_list)){
            
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_weighted_across_donor <- as.numeric(df$Mean_z_score_weighted_across_donor)
            df <- df[order( -(df$Mean_z_score_weighted_across_donor)),]
            df <- subset(df, Mean_z_score_weighted_across_donor >= 1)
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_weighted_across_donor)
            df$celltype <- celltype2
            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")
        BA4_all_conditions_weighted_df <- fill

    ##

    ###########################
    ## All conditions BA9
    ###########################
    ## code
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA9"
        par_status = "all_conditions"
        par_prep = "CombatSeq"

        ################################################
        ## LIME outputs weighted by percent expression
        ## This function subsets the weighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_weighted_across_donor = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_weighted_across_donor, celltype)

        for(celltype2 in unique(celltype_list)){
            
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_weighted_across_donor <- as.numeric(df$Mean_z_score_weighted_across_donor)
            df <- df[order( -(df$Mean_z_score_weighted_across_donor)),]
            df <- subset(df, Mean_z_score_weighted_across_donor >= 1)
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_weighted_across_donor)
            df$celltype <- celltype2
            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")
        BA9_all_conditions_weighted_df <- fill

    ##

    ###########################
    ## Plot BA4
    ###########################
    ## code
        BA4_SALS_weighted_df$disease_cat <- "SALS"
        BA4_C9ALS_weighted_df$disease_cat <- "C9ALS"
        BA4_SFTLD_weighted_df$disease_cat <- "SFTLD"
        BA4_C9FTLD_weighted_df$disease_cat <- "C9FTLD"
        BA4_all_conditions_weighted_df$disease_cat <- "All"

        BA4_weighted_counts <- rbind(BA4_SALS_weighted_df, BA4_C9ALS_weighted_df, BA4_SFTLD_weighted_df, BA4_C9FTLD_weighted_df, BA4_all_conditions_weighted_df)

        counts_df <- data.frame(table(BA4_weighted_counts$celltype, BA4_weighted_counts$disease_cat ))
        
        counts_df$Var2 <- factor(counts_df$Var2, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All"))
        counts_df$Var1 <- factor(counts_df$Var1, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))

        ggplot(counts_df, aes(x = Var2, y = Freq, fill = Var2)) + 
        theme_bw() + 
        geom_bar(stat = 'identity') +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(colour = "black", angle = 90, hjust = 1, vjust = 0.5),
            axis.text.y = element_text(colour = "black"),
            axis.title.y = element_text(face = "bold"),
            axis.title.x = element_blank(),
        ) +
        facet_grid(~Var1, scales = "free_x") +
        ylab("N genes\n(Z-score > 1)") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 2, width = 15)
    ##

    ###########################
    ## Plot BA9
    ###########################
    ## code
        BA9_SALS_weighted_df$disease_cat <- "SALS"
        BA9_C9ALS_weighted_df$disease_cat <- "C9ALS"
        BA9_SFTLD_weighted_df$disease_cat <- "SFTLD"
        BA9_C9FTLD_weighted_df$disease_cat <- "C9FTLD"
        BA9_all_conditions_weighted_df$disease_cat <- "All"

        BA9_weighted_counts <- rbind(BA9_SALS_weighted_df, BA9_C9ALS_weighted_df, BA9_SFTLD_weighted_df, BA9_C9FTLD_weighted_df, BA9_all_conditions_weighted_df)

        counts_df <- data.frame(table(BA9_weighted_counts$celltype, BA9_weighted_counts$disease_cat ))
        
        counts_df$Var2 <- factor(counts_df$Var2, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All"))
        counts_df$Var1 <- factor(counts_df$Var1, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))

        ggplot(counts_df, aes(x = Var2, y = Freq, fill = Var2)) + 
        theme_bw() + 
        geom_bar(stat = 'identity') +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(colour = "black", angle = 90, hjust = 1, vjust = 0.5),
            axis.text.y = element_text(colour = "black"),
            axis.title.y = element_text(face = "bold"),
            axis.title.x = element_blank(),
        ) +
        facet_grid(~Var1, scales = "free_x") +
        ylab("N genes\n(Z-score > 1)") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 2, width = 15)
    ##

##

############################################################################ Weighted permutation optimal threshold test -- lambda
############################################################################
############################################################################
############################################################################ 

## code
    ###########################
    ## SALS BA4
    ###########################
    ## code
        #celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",    "OPC",  "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA4"
        par_status = "SALS"
        par_prep = "CombatSeq"
        
        ################################################
        ## LIME outputs weighted by percent expression
        ## This function subsets the weighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        region = "fill"
        group = "fill"
        celltype = "fill"
        n_genes = 0
        test_accuracy = 0
        mean_accuracy = 0
        method = "fill"
        fill <- data.frame(region,group, celltype, n_genes, test_accuracy, mean_accuracy, method)
        

        for(celltype2 in unique(celltype_list)){
        
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/permutation_weighted_report_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df$five_fold_accuracies <- strsplit(as.character(df$five_fold_accuracies), "; ")
            df$five_fold_accuracies <- lapply(df$five_fold_accuracies, as.numeric)
            df <- cbind(df, do.call(rbind, df$five_fold_accuracies))
            colnames(df)[13:17] <- paste0("fold_", 1:5)
            df$mean_accuracy <- rowMeans(df[, c("test_accuracy", "fold_1", "fold_2", "fold_3", "fold_4", "fold_5")], na.rm = TRUE)
            df <- df %>% dplyr::select(region, group, celltype, n_genes, test_accuracy, mean_accuracy, method)
            fill <- rbind(fill, df)

        }

        fill <- subset(fill, region != "fill")
        BA4_SALS_weighted_df <- fill
    ##

    ###########################
    ## C9ALS BA4
    ###########################
    ## code
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA4"
        par_status = "C9ALS"
        par_prep = "CombatSeq"
        
        ################################################
        ## LIME outputs weighted by percent expression
        ## This function subsets the weighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        region = "fill"
        group = "fill"
        celltype = "fill"
        n_genes = 0
        test_accuracy = 0
        mean_accuracy = 0
        method = "fill"
        fill <- data.frame(region,group, celltype, n_genes, test_accuracy, mean_accuracy, method)
        
        for(celltype2 in unique(celltype_list)){
        
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/permutation_weighted_report_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df$five_fold_accuracies <- strsplit(as.character(df$five_fold_accuracies), "; ")
            df$five_fold_accuracies <- lapply(df$five_fold_accuracies, as.numeric)
            df <- cbind(df, do.call(rbind, df$five_fold_accuracies))
            colnames(df)[13:17] <- paste0("fold_", 1:5)
            df$mean_accuracy <- rowMeans(df[, c("test_accuracy", "fold_1", "fold_2", "fold_3", "fold_4", "fold_5")], na.rm = TRUE)
            df <- df %>% dplyr::select(region, group, celltype, n_genes, test_accuracy, mean_accuracy, method)
            fill <- rbind(fill, df)

        }

        fill <- subset(fill, region != "fill")
        BA4_C9ALS_weighted_df <- fill
    ##

    ###########################
    ## SFTLD BA4
    ###########################
    ## code
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA4"
        par_status = "SFTLD"
        par_prep = "CombatSeq"
        
        ################################################
        ## LIME outputs weighted by percent expression
        ## This function subsets the weighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        region = "fill"
        group = "fill"
        celltype = "fill"
        n_genes = 0
        test_accuracy = 0
        mean_accuracy = 0
        method = "fill"
        fill <- data.frame(region,group, celltype, n_genes, test_accuracy, mean_accuracy, method)
        
        for(celltype2 in unique(celltype_list)){
        
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/permutation_weighted_report_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df$five_fold_accuracies <- strsplit(as.character(df$five_fold_accuracies), "; ")
            df$five_fold_accuracies <- lapply(df$five_fold_accuracies, as.numeric)
            df <- cbind(df, do.call(rbind, df$five_fold_accuracies))
            colnames(df)[13:17] <- paste0("fold_", 1:5)
            df$mean_accuracy <- rowMeans(df[, c("test_accuracy", "fold_1", "fold_2", "fold_3", "fold_4", "fold_5")], na.rm = TRUE)
            df <- df %>% dplyr::select(region, group, celltype, n_genes, test_accuracy, mean_accuracy, method)
            fill <- rbind(fill, df)

        }

        fill <- subset(fill, region != "fill")
        BA4_SFTLD_weighted_df <- fill
    ##

    ###########################
    ## C9FTLD BA4 
    ###########################
    ## code
        
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA4"
        par_status = "C9FTLD"
        par_prep = "CombatSeq"
        
        ################################################
        ## LIME outputs weighted by percent expression
        ## This function subsets the weighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        region = "fill"
        group = "fill"
        celltype = "fill"
        n_genes = 0
        test_accuracy = 0
        mean_accuracy = 0
        method = "fill"
        fill <- data.frame(region,group, celltype, n_genes, test_accuracy, mean_accuracy, method)
        
        for(celltype2 in unique(celltype_list)){
        
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/permutation_weighted_report_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df$five_fold_accuracies <- strsplit(as.character(df$five_fold_accuracies), "; ")
            df$five_fold_accuracies <- lapply(df$five_fold_accuracies, as.numeric)
            df <- cbind(df, do.call(rbind, df$five_fold_accuracies))
            colnames(df)[13:17] <- paste0("fold_", 1:5)
            df$mean_accuracy <- rowMeans(df[, c("test_accuracy", "fold_1", "fold_2", "fold_3", "fold_4", "fold_5")], na.rm = TRUE)
            df <- df %>% dplyr::select(region, group, celltype, n_genes, test_accuracy, mean_accuracy, method)
            fill <- rbind(fill, df)

        }

        fill <- subset(fill, region != "fill")
        BA4_C9FTLD_weighted_df <- fill
    ##

    ###########################
    ## SALS BA9
    ###########################
    ## code
        #celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "OPC",  "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA9"
        par_status = "SALS"
        par_prep = "CombatSeq"
        
        ################################################
        ## LIME outputs weighted by percent expression
        ## This function subsets the weighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        region = "fill"
        group = "fill"
        celltype = "fill"
        n_genes = 0
        test_accuracy = 0
        mean_accuracy = 0
        method = "fill"
        fill <- data.frame(region,group, celltype, n_genes, test_accuracy, mean_accuracy, method)
        

        for(celltype2 in unique(celltype_list)){
        
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/permutation_weighted_report_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df$five_fold_accuracies <- strsplit(as.character(df$five_fold_accuracies), "; ")
            df$five_fold_accuracies <- lapply(df$five_fold_accuracies, as.numeric)
            df <- cbind(df, do.call(rbind, df$five_fold_accuracies))
            colnames(df)[13:17] <- paste0("fold_", 1:5)
            df$mean_accuracy <- rowMeans(df[, c("test_accuracy", "fold_1", "fold_2", "fold_3", "fold_4", "fold_5")], na.rm = TRUE)
            df <- df %>% dplyr::select(region, group, celltype, n_genes, test_accuracy, mean_accuracy, method)
            fill <- rbind(fill, df)

        }

        fill <- subset(fill, region != "fill")
        BA9_SALS_weighted_df <- fill
    ##

    ###########################
    ## C9ALS BA9
    ###########################
    ## code
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA9"
        par_status = "C9ALS"
        par_prep = "CombatSeq"
        
        ################################################
        ## LIME outputs weighted by percent expression
        ## This function subsets the weighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        region = "fill"
        group = "fill"
        celltype = "fill"
        n_genes = 0
        test_accuracy = 0
        mean_accuracy = 0
        method = "fill"
        fill <- data.frame(region,group, celltype, n_genes, test_accuracy, mean_accuracy, method)
        
        for(celltype2 in unique(celltype_list)){
        
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/permutation_weighted_report_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df$five_fold_accuracies <- strsplit(as.character(df$five_fold_accuracies), "; ")
            df$five_fold_accuracies <- lapply(df$five_fold_accuracies, as.numeric)
            df <- cbind(df, do.call(rbind, df$five_fold_accuracies))
            colnames(df)[13:17] <- paste0("fold_", 1:5)
            df$mean_accuracy <- rowMeans(df[, c("test_accuracy", "fold_1", "fold_2", "fold_3", "fold_4", "fold_5")], na.rm = TRUE)
            df <- df %>% dplyr::select(region, group, celltype, n_genes, test_accuracy, mean_accuracy, method)
            fill <- rbind(fill, df)

        }

        fill <- subset(fill, region != "fill")
        BA9_C9ALS_weighted_df <- fill
    ##

    ###########################
    ## SFTLD BA9
    ###########################
    ## code
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA9"
        par_status = "SFTLD"
        par_prep = "CombatSeq"
        
        ################################################
        ## LIME outputs weighted by percent expression
        ## This function subsets the weighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        region = "fill"
        group = "fill"
        celltype = "fill"
        n_genes = 0
        test_accuracy = 0
        mean_accuracy = 0
        method = "fill"
        fill <- data.frame(region,group, celltype, n_genes, test_accuracy, mean_accuracy, method)
        
        for(celltype2 in unique(celltype_list)){
        
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/permutation_weighted_report_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df$five_fold_accuracies <- strsplit(as.character(df$five_fold_accuracies), "; ")
            df$five_fold_accuracies <- lapply(df$five_fold_accuracies, as.numeric)
            df <- cbind(df, do.call(rbind, df$five_fold_accuracies))
            colnames(df)[13:17] <- paste0("fold_", 1:5)
            df$mean_accuracy <- rowMeans(df[, c("test_accuracy", "fold_1", "fold_2", "fold_3", "fold_4", "fold_5")], na.rm = TRUE)
            df <- df %>% dplyr::select(region, group, celltype, n_genes, test_accuracy, mean_accuracy, method)
            fill <- rbind(fill, df)

        }

        fill <- subset(fill, region != "fill")
        BA9_SFTLD_weighted_df <- fill
    ##

    ###########################
    ## C9FTLD BA9 
    ###########################
    ## code
        
        #celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",    "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",    "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA9"
        par_status = "C9FTLD"
        par_prep = "CombatSeq"
        
        ################################################
        ## LIME outputs weighted by percent expression
        ## This function subsets the weighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        region = "fill"
        group = "fill"
        celltype = "fill"
        n_genes = 0
        test_accuracy = 0
        mean_accuracy = 0
        method = "fill"
        fill <- data.frame(region,group, celltype, n_genes, test_accuracy, mean_accuracy, method)
        
        for(celltype2 in unique(celltype_list)){
        
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/permutation_weighted_report_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df$five_fold_accuracies <- strsplit(as.character(df$five_fold_accuracies), "; ")
            df$five_fold_accuracies <- lapply(df$five_fold_accuracies, as.numeric)
            df <- cbind(df, do.call(rbind, df$five_fold_accuracies))
            colnames(df)[13:17] <- paste0("fold_", 1:5)
            df$mean_accuracy <- rowMeans(df[, c("test_accuracy", "fold_1", "fold_2", "fold_3", "fold_4", "fold_5")], na.rm = TRUE)
            df <- df %>% dplyr::select(region, group, celltype, n_genes, test_accuracy, mean_accuracy, method)
            fill <- rbind(fill, df)

        }

        fill <- subset(fill, region != "fill")
        BA9_C9FTLD_weighted_df <- fill
    ##

    ###########################
    ## SALS BA4 -- Lambda approach and plot (for now we will include the randoms for Sali but I think we should remove this eventually)
    ###########################
    ## code
        ## sumarize SALS
        BA4_SALS_weighted_df$disease_cat <- "SALS"
        
        BA4_SALS_weighted_df_summary <- BA4_SALS_weighted_df %>%
            dplyr::group_by(group, celltype, n_genes, method) %>%
            dplyr::summarise(mean_accuracy_by_group = mean(mean_accuracy, na.rm = TRUE))
        
        BA4_SALS_weighted_df_summary <- data.frame(BA4_SALS_weighted_df_summary)

        BA4_SALS_weighted_df_summary <- subset(BA4_SALS_weighted_df_summary, method == "LIME" )
        
        BA4_SALS_weighted_df_summary <- subset(BA4_SALS_weighted_df_summary, n_genes >= 25 )


        # Function to calculate the score with regularization
        calculate_score <- function(df, lambda) {
        df %>%
            mutate(score = mean_accuracy_by_group - lambda * n_genes)
        }

        lambda_values <- seq(0, 0.1, by = 0.00005)  # Try different values of lambda
        result <- data.frame()

        for (lambda in lambda_values) {
        df_with_score <- calculate_score(BA4_SALS_weighted_df_summary, lambda)
        
        # Find the optimal row for each celltype
        optimal_for_lambda <- df_with_score %>%
            group_by(celltype) %>%
            slice_max(score, n = 1)
        
        # Store the results for each lambda value
        optimal_for_lambda$lambda <- lambda
        result <- bind_rows(result, optimal_for_lambda)
        }

        head(result, n = 500)

        # View the result
        print(result)

        ## Plot
        BA4_SALS_weighted_df_summary <- BA4_SALS_weighted_df %>%
            dplyr::group_by(group, celltype, n_genes, method) %>%
            dplyr::summarise(mean_accuracy_by_group = mean(mean_accuracy, na.rm = TRUE))

        BA4_SALS_weighted_df_summary <- subset(BA4_SALS_weighted_df_summary, n_genes >= 25)


        result <- subset(result, lambda == 0.00025)
        result <- result %>% dplyr::select(celltype, n_genes)
        
        colnames(result) <- c("celltype", "optimal_genes")
        
        write.csv(result, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_lambda_0.00025_gene_set_BA4_SALS_narval.csv'))


        BA4_SALS_weighted_df_summary <- merge(BA4_SALS_weighted_df_summary, result, by = "celltype")
        BA4_SALS_weighted_df_summary$method <- factor(BA4_SALS_weighted_df_summary$method, levels = c("LIME", "random_HVG", "random_gene"))
        BA4_SALS_weighted_df_summary$celltype <- factor(BA4_SALS_weighted_df_summary$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))
        
        
        ## custom label
        BA4_SALS_weighted_df_summary <- BA4_SALS_weighted_df_summary %>%
            mutate(label = paste(celltype, "\n(", optimal_genes, " genes)", sep = ""))

        BA4_SALS_weighted_df_summary <- BA4_SALS_weighted_df_summary %>%
            arrange(celltype)  # Arrange by the 'celltype' factor

        BA4_SALS_weighted_df_summary$label <- factor(BA4_SALS_weighted_df_summary$label, levels = unique(BA4_SALS_weighted_df_summary$label))
        BA4_SALS_weighted_df_summary <- subset(BA4_SALS_weighted_df_summary, method == "LIME")
        
        ggplot(BA4_SALS_weighted_df_summary, aes(x = n_genes, y = mean_accuracy_by_group, colour = method, group = method)) + 
        theme_bw() + 
        geom_vline(aes(xintercept = optimal_genes), linetype = "dashed", color = "black") +
        geom_point(size = 1) +
        geom_line() +
        scale_x_reverse() +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(colour = "black", angle = 90, hjust = 1, vjust = 0.5),
            axis.text.y = element_text(colour = "black"),
            axis.title.y = element_text(face = "bold"),
            axis.title.x = element_text(face = "bold"),
        ) +
        facet_grid(group ~ label, scales = "free_x") + 
        ylab("Accuracy") +
        xlab("N genes") +
        scale_colour_manual(values = c("#F18F01","#048BA8" , "#2E4057"))
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 2, width = 15)
    ##

    ###########################
    ## C9ALS BA4 -- Lambda approach and plot
    ###########################
    ## code
        ## sumarize C9ALS
        BA4_C9ALS_weighted_df$disease_cat <- "C9ALS"
        BA4_C9ALS_weighted_df_summary <- BA4_C9ALS_weighted_df %>%
            dplyr::group_by(group, celltype, n_genes, method) %>%
            dplyr::summarise(mean_accuracy_by_group = mean(mean_accuracy, na.rm = TRUE))
        BA4_C9ALS_weighted_df_summary <- data.frame(BA4_C9ALS_weighted_df_summary)

        BA4_C9ALS_weighted_df_summary <- subset(BA4_C9ALS_weighted_df_summary, method == "LIME" )

        ## get rid of the fibro 12 because it is an outlier
        #BA4_C9ALS_weighted_df_summary <- BA4_C9ALS_weighted_df_summary %>%
        #    filter(!(celltype == "Fibro" & n_genes == 12))
        
        BA4_C9ALS_weighted_df_summary <- subset(BA4_C9ALS_weighted_df_summary, n_genes >= 25 )


        # Function to calculate the score with regularization
        calculate_score <- function(df, lambda) {
        df %>%
            mutate(score = mean_accuracy_by_group - lambda * n_genes)
        }

        lambda_values <- seq(0, 0.1, by = 0.00005)  # Try different values of lambda
        result <- data.frame()

        for (lambda in lambda_values) {
        df_with_score <- calculate_score(BA4_C9ALS_weighted_df_summary, lambda)
        
        # Find the optimal row for each celltype
        optimal_for_lambda <- df_with_score %>%
            group_by(celltype) %>%
            slice_max(score, n = 1)
        
        # Store the results for each lambda value
        optimal_for_lambda$lambda <- lambda
        result <- bind_rows(result, optimal_for_lambda)
        }

        head(result, n = 500)

        # View the result
        print(result)

        ## Plot
        BA4_C9ALS_weighted_df_summary <- BA4_C9ALS_weighted_df %>%
            dplyr::group_by(group, celltype, n_genes, method) %>%
            dplyr::summarise(mean_accuracy_by_group = mean(mean_accuracy, na.rm = TRUE))

        BA4_C9ALS_weighted_df_summary <- subset(BA4_C9ALS_weighted_df_summary, n_genes >= 25)


        result <- subset(result, lambda == 0.00025)
        result <- result %>% dplyr::select(celltype, n_genes)
        
        colnames(result) <- c("celltype", "optimal_genes")
        write.csv(result, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_lambda_0.00025_gene_set_BA4_C9ALS_narval.csv'))


        BA4_C9ALS_weighted_df_summary <- merge(BA4_C9ALS_weighted_df_summary, result, by = "celltype")


        BA4_C9ALS_weighted_df_summary$method <- factor(BA4_C9ALS_weighted_df_summary$method, levels = c("LIME", "random_HVG", "random_gene"))
        BA4_C9ALS_weighted_df_summary$celltype <- factor(BA4_C9ALS_weighted_df_summary$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))
        
        
        ## custom label
        BA4_C9ALS_weighted_df_summary <- BA4_C9ALS_weighted_df_summary %>%
            mutate(label = paste(celltype, "\n(", optimal_genes, " genes)", sep = ""))

        BA4_C9ALS_weighted_df_summary <- BA4_C9ALS_weighted_df_summary %>%
            arrange(celltype)  # Arrange by the 'celltype' factor

        BA4_C9ALS_weighted_df_summary$label <- factor(BA4_C9ALS_weighted_df_summary$label, levels = unique(BA4_C9ALS_weighted_df_summary$label))
        BA4_C9ALS_weighted_df_summary <- subset(BA4_C9ALS_weighted_df_summary, method == "LIME")


        ggplot(BA4_C9ALS_weighted_df_summary, aes(x = n_genes, y = mean_accuracy_by_group, colour = method, group = method)) + 
        theme_bw() + 
        geom_vline(aes(xintercept = optimal_genes), linetype = "dashed", color = "black") +
        geom_point(size = 1) +
        geom_line() +
        scale_x_reverse() +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(colour = "black", angle = 90, hjust = 1, vjust = 0.5),
            axis.text.y = element_text(colour = "black"),
            axis.title.y = element_text(face = "bold"),
            axis.title.x = element_text(face = "bold"),
        ) +
        facet_grid(group ~ label, scales = "free_x") + 
        ylab("Accuracy") +
        xlab("N genes") +
        scale_colour_manual(values = c("#F18F01","#048BA8" , "#2E4057"))
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp2.pdf', height = 2, width = 15)
    ##
    
    ###########################
    ## SFTLD BA4 -- Lambda approach and plot
    ###########################
    ## code
        ## sumarize SALS
        BA4_SFTLD_weighted_df$disease_cat <- "SFTLD"
        BA4_SFTLD_weighted_df_summary <- BA4_SFTLD_weighted_df %>%
            dplyr::group_by(group, celltype, n_genes, method) %>%
            dplyr::summarise(mean_accuracy_by_group = mean(mean_accuracy, na.rm = TRUE))
        BA4_SFTLD_weighted_df_summary <- data.frame(BA4_SFTLD_weighted_df_summary)

        BA4_SFTLD_weighted_df_summary <- subset(BA4_SFTLD_weighted_df_summary, method == "LIME" )
        
        BA4_SFTLD_weighted_df_summary <- subset(BA4_SFTLD_weighted_df_summary, n_genes >= 25 )


        # Function to calculate the score with regularization
        calculate_score <- function(df, lambda) {
        df %>%
            mutate(score = mean_accuracy_by_group - lambda * n_genes)
        }

        lambda_values <- seq(0, 0.1, by = 0.00005)  # Try different values of lambda
        result <- data.frame()

        for (lambda in lambda_values) {
        df_with_score <- calculate_score(BA4_SFTLD_weighted_df_summary, lambda)
        
        # Find the optimal row for each celltype
        optimal_for_lambda <- df_with_score %>%
            group_by(celltype) %>%
            slice_max(score, n = 1)
        
        # Store the results for each lambda value
        optimal_for_lambda$lambda <- lambda
        result <- bind_rows(result, optimal_for_lambda)
        }

        head(result, n = 500)

        # View the result
        print(result)

        ## Plot
        BA4_SFTLD_weighted_df_summary <- BA4_SFTLD_weighted_df %>%
            dplyr::group_by(group, celltype, n_genes, method) %>%
            dplyr::summarise(mean_accuracy_by_group = mean(mean_accuracy, na.rm = TRUE))

        BA4_SFTLD_weighted_df_summary <- subset(BA4_SFTLD_weighted_df_summary, n_genes >= 25)


        result <- subset(result, lambda == 0.00025)
        result <- result %>% dplyr::select(celltype, n_genes)
        
        colnames(result) <- c("celltype", "optimal_genes")
        write.csv(result, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_lambda_0.00025_gene_set_BA4_SFTLD_narval.csv'))


        BA4_SFTLD_weighted_df_summary <- merge(BA4_SFTLD_weighted_df_summary, result, by = "celltype")


        BA4_SFTLD_weighted_df_summary$method <- factor(BA4_SFTLD_weighted_df_summary$method, levels = c("LIME", "random_HVG", "random_gene"))
        BA4_SFTLD_weighted_df_summary$celltype <- factor(BA4_SFTLD_weighted_df_summary$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))
        
        
        ## custom label
        BA4_SFTLD_weighted_df_summary <- BA4_SFTLD_weighted_df_summary %>%
            mutate(label = paste(celltype, "\n(", optimal_genes, " genes)", sep = ""))

        BA4_SFTLD_weighted_df_summary <- BA4_SFTLD_weighted_df_summary %>%
            arrange(celltype)  # Arrange by the 'celltype' factor

        BA4_SFTLD_weighted_df_summary$label <- factor(BA4_SFTLD_weighted_df_summary$label, levels = unique(BA4_SFTLD_weighted_df_summary$label))
        BA4_SFTLD_weighted_df_summary <- subset(BA4_SFTLD_weighted_df_summary, method == "LIME")

        ggplot(BA4_SFTLD_weighted_df_summary, aes(x = n_genes, y = mean_accuracy_by_group, colour = method, group = method)) + 
        theme_bw() + 
        geom_vline(aes(xintercept = optimal_genes), linetype = "dashed", color = "black") +
        geom_point(size = 1) +
        geom_line() +
        scale_x_reverse() +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(colour = "black", angle = 90, hjust = 1, vjust = 0.5),
            axis.text.y = element_text(colour = "black"),
            axis.title.y = element_text(face = "bold"),
            axis.title.x = element_text(face = "bold"),
        ) +
        facet_grid(group ~ label, scales = "free_x") + 
        ylab("Accuracy") +
        xlab("N genes") +
        scale_colour_manual(values = c("#F18F01","#048BA8" , "#2E4057"))
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp3.pdf', height = 2, width = 15)
    ##

    ###########################
    ## C9FTLD BA4 -- Lambda approach and plot
    ###########################
    ## code
        ## sumarize C9FTLD
        BA4_C9FTLD_weighted_df$disease_cat <- "C9FTLD"
        BA4_C9FTLD_weighted_df_summary <- BA4_C9FTLD_weighted_df %>%
            dplyr::group_by(group, celltype, n_genes, method) %>%
            dplyr::summarise(mean_accuracy_by_group = mean(mean_accuracy, na.rm = TRUE))
        BA4_C9FTLD_weighted_df_summary <- data.frame(BA4_C9FTLD_weighted_df_summary)

        BA4_C9FTLD_weighted_df_summary <- subset(BA4_C9FTLD_weighted_df_summary, method == "LIME" )
        
        BA4_C9FTLD_weighted_df_summary <- subset(BA4_C9FTLD_weighted_df_summary, n_genes >= 25 )


        # Function to calculate the score with regularization
        calculate_score <- function(df, lambda) {
        df %>%
            mutate(score = mean_accuracy_by_group - lambda * n_genes)
        }

        lambda_values <- seq(0, 0.1, by = 0.00005)  # Try different values of lambda
        result <- data.frame()

        for (lambda in lambda_values) {
        df_with_score <- calculate_score(BA4_C9FTLD_weighted_df_summary, lambda)
        
        # Find the optimal row for each celltype
        optimal_for_lambda <- df_with_score %>%
            group_by(celltype) %>%
            slice_max(score, n = 1)
        
        # Store the results for each lambda value
        optimal_for_lambda$lambda <- lambda
        result <- bind_rows(result, optimal_for_lambda)
        }

        head(result, n = 500)

        # View the result
        print(result)

        ## Plot
        BA4_C9FTLD_weighted_df_summary <- BA4_C9FTLD_weighted_df %>%
            dplyr::group_by(group, celltype, n_genes, method) %>%
            dplyr::summarise(mean_accuracy_by_group = mean(mean_accuracy, na.rm = TRUE))

        BA4_C9FTLD_weighted_df_summary <- subset(BA4_C9FTLD_weighted_df_summary, n_genes >= 25)


        result <- subset(result, lambda == 0.00025)
        result <- result %>% dplyr::select(celltype, n_genes)
        
        colnames(result) <- c("celltype", "optimal_genes")
        write.csv(result, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_lambda_0.00025_gene_set_BA4_C9FTLD_narval.csv'))

        BA4_C9FTLD_weighted_df_summary <- merge(BA4_C9FTLD_weighted_df_summary, result, by = "celltype")


        BA4_C9FTLD_weighted_df_summary$method <- factor(BA4_C9FTLD_weighted_df_summary$method, levels = c("LIME", "random_HVG", "random_gene"))
        BA4_C9FTLD_weighted_df_summary$celltype <- factor(BA4_C9FTLD_weighted_df_summary$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))
        
        
        ## custom label
        BA4_C9FTLD_weighted_df_summary <- BA4_C9FTLD_weighted_df_summary %>%
            mutate(label = paste(celltype, "\n(", optimal_genes, " genes)", sep = ""))

        BA4_C9FTLD_weighted_df_summary <- BA4_C9FTLD_weighted_df_summary %>%
            arrange(celltype)  # Arrange by the 'celltype' factor

        BA4_C9FTLD_weighted_df_summary$label <- factor(BA4_C9FTLD_weighted_df_summary$label, levels = unique(BA4_C9FTLD_weighted_df_summary$label))
        BA4_C9FTLD_weighted_df_summary <- subset(BA4_C9FTLD_weighted_df_summary, method == "LIME")

        ggplot(BA4_C9FTLD_weighted_df_summary, aes(x = n_genes, y = mean_accuracy_by_group, colour = method, group = method)) + 
        theme_bw() + 
        geom_vline(aes(xintercept = optimal_genes), linetype = "dashed", color = "black") +
        geom_point(size = 1) +
        geom_line() +
        scale_x_reverse() +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(colour = "black", angle = 90, hjust = 1, vjust = 0.5),
            axis.text.y = element_text(colour = "black"),
            axis.title.y = element_text(face = "bold"),
            axis.title.x = element_text(face = "bold"),
        ) +
        facet_grid(group ~ label, scales = "free_x") + 
        ylab("Accuracy") +
        xlab("N genes") +
        scale_colour_manual(values = c("#F18F01","#048BA8" , "#2E4057"))
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp4.pdf', height = 2, width = 15)
    ##

    ###########################
    ## SALS BA9 -- Lambda approach and plot (for now we will include the randoms for Sali but I think we should remove this eventually)
    ###########################
    ## code
        ## sumarize SALS
        BA9_SALS_weighted_df$disease_cat <- "SALS"
        
        BA9_SALS_weighted_df_summary <- BA9_SALS_weighted_df %>%
            dplyr::group_by(group, celltype, n_genes, method) %>%
            dplyr::summarise(mean_accuracy_by_group = mean(mean_accuracy, na.rm = TRUE))
        
        BA9_SALS_weighted_df_summary <- data.frame(BA9_SALS_weighted_df_summary)

        BA9_SALS_weighted_df_summary <- subset(BA9_SALS_weighted_df_summary, method == "LIME" )
        
        BA9_SALS_weighted_df_summary <- subset(BA9_SALS_weighted_df_summary, n_genes >= 25 )


        # Function to calculate the score with regularization
        calculate_score <- function(df, lambda) {
        df %>%
            mutate(score = mean_accuracy_by_group - lambda * n_genes)
        }

        lambda_values <- seq(0, 0.1, by = 0.00005)  # Try different values of lambda
        result <- data.frame()

        for (lambda in lambda_values) {
        df_with_score <- calculate_score(BA9_SALS_weighted_df_summary, lambda)
        
        # Find the optimal row for each celltype
        optimal_for_lambda <- df_with_score %>%
            group_by(celltype) %>%
            slice_max(score, n = 1)
        
        # Store the results for each lambda value
        optimal_for_lambda$lambda <- lambda
        result <- bind_rows(result, optimal_for_lambda)
        }

        head(result, n = 500)

        # View the result
        print(result)

        ## Plot
        BA9_SALS_weighted_df_summary <- BA9_SALS_weighted_df %>%
            dplyr::group_by(group, celltype, n_genes, method) %>%
            dplyr::summarise(mean_accuracy_by_group = mean(mean_accuracy, na.rm = TRUE))

        BA9_SALS_weighted_df_summary <- subset(BA9_SALS_weighted_df_summary, n_genes >= 25)


        result <- subset(result, lambda == 0.00025)
        result <- result %>% dplyr::select(celltype, n_genes)
        
        colnames(result) <- c("celltype", "optimal_genes")
        
        write.csv(result, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_lambda_0.00025_gene_set_BA9_SALS_narval.csv'))


        BA9_SALS_weighted_df_summary <- merge(BA9_SALS_weighted_df_summary, result, by = "celltype")
        BA9_SALS_weighted_df_summary$method <- factor(BA9_SALS_weighted_df_summary$method, levels = c("LIME", "random_HVG", "random_gene"))
        BA9_SALS_weighted_df_summary$celltype <- factor(BA9_SALS_weighted_df_summary$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))
        
        
        ## custom label
        BA9_SALS_weighted_df_summary <- BA9_SALS_weighted_df_summary %>%
            mutate(label = paste(celltype, "\n(", optimal_genes, " genes)", sep = ""))

        BA9_SALS_weighted_df_summary <- BA9_SALS_weighted_df_summary %>%
            arrange(celltype)  # Arrange by the 'celltype' factor

        BA9_SALS_weighted_df_summary$label <- factor(BA9_SALS_weighted_df_summary$label, levels = unique(BA9_SALS_weighted_df_summary$label))
        BA9_SALS_weighted_df_summary <- subset(BA9_SALS_weighted_df_summary, method == "LIME")

        ggplot(BA9_SALS_weighted_df_summary, aes(x = n_genes, y = mean_accuracy_by_group, colour = method, group = method)) + 
        theme_bw() + 
        geom_vline(aes(xintercept = optimal_genes), linetype = "dashed", color = "black") +
        geom_point(size = 1) +
        geom_line() +
        scale_x_reverse() +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(colour = "black", angle = 90, hjust = 1, vjust = 0.5),
            axis.text.y = element_text(colour = "black"),
            axis.title.y = element_text(face = "bold"),
            axis.title.x = element_text(face = "bold"),
        ) +
        facet_grid(group ~ label, scales = "free_x") + 
        ylab("Accuracy") +
        xlab("N genes") +
        scale_colour_manual(values = c("#F18F01","#048BA8" , "#2E4057"))
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 2, width = 15)
    ##

    ###########################
    ## C9ALS BA9 -- Lambda approach and plot
    ###########################
    ## code
        ## sumarize C9ALS
        BA9_C9ALS_weighted_df$disease_cat <- "C9ALS"
        BA9_C9ALS_weighted_df_summary <- BA9_C9ALS_weighted_df %>%
            dplyr::group_by(group, celltype, n_genes, method) %>%
            dplyr::summarise(mean_accuracy_by_group = mean(mean_accuracy, na.rm = TRUE))
        BA9_C9ALS_weighted_df_summary <- data.frame(BA9_C9ALS_weighted_df_summary)

        BA9_C9ALS_weighted_df_summary <- subset(BA9_C9ALS_weighted_df_summary, method == "LIME" )

        ## get rid of the fibro 12 because it is an outlier
        #BA9_C9ALS_weighted_df_summary <- BA9_C9ALS_weighted_df_summary %>%
        #    filter(!(celltype == "Fibro" & n_genes == 12))
        
        BA9_C9ALS_weighted_df_summary <- subset(BA9_C9ALS_weighted_df_summary, n_genes >= 25 )


        # Function to calculate the score with regularization
        calculate_score <- function(df, lambda) {
        df %>%
            mutate(score = mean_accuracy_by_group - lambda * n_genes)
        }

        lambda_values <- seq(0, 0.1, by = 0.00005)  # Try different values of lambda
        result <- data.frame()

        for (lambda in lambda_values) {
        df_with_score <- calculate_score(BA9_C9ALS_weighted_df_summary, lambda)
        
        # Find the optimal row for each celltype
        optimal_for_lambda <- df_with_score %>%
            group_by(celltype) %>%
            slice_max(score, n = 1)
        
        # Store the results for each lambda value
        optimal_for_lambda$lambda <- lambda
        result <- bind_rows(result, optimal_for_lambda)
        }

        head(result, n = 500)

        # View the result
        print(result)

        ## Plot
        BA9_C9ALS_weighted_df_summary <- BA9_C9ALS_weighted_df %>%
            dplyr::group_by(group, celltype, n_genes, method) %>%
            dplyr::summarise(mean_accuracy_by_group = mean(mean_accuracy, na.rm = TRUE))

        BA9_C9ALS_weighted_df_summary <- subset(BA9_C9ALS_weighted_df_summary, n_genes >= 25)


        result <- subset(result, lambda == 0.00025)
        result <- result %>% dplyr::select(celltype, n_genes)
        
        colnames(result) <- c("celltype", "optimal_genes")
        write.csv(result, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_lambda_0.00025_gene_set_BA9_C9ALS_narval.csv'))


        BA9_C9ALS_weighted_df_summary <- merge(BA9_C9ALS_weighted_df_summary, result, by = "celltype")


        BA9_C9ALS_weighted_df_summary$method <- factor(BA9_C9ALS_weighted_df_summary$method, levels = c("LIME", "random_HVG", "random_gene"))
        BA9_C9ALS_weighted_df_summary$celltype <- factor(BA9_C9ALS_weighted_df_summary$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))
        
        
        ## custom label
        BA9_C9ALS_weighted_df_summary <- BA9_C9ALS_weighted_df_summary %>%
            mutate(label = paste(celltype, "\n(", optimal_genes, " genes)", sep = ""))

        BA9_C9ALS_weighted_df_summary <- BA9_C9ALS_weighted_df_summary %>%
            arrange(celltype)  # Arrange by the 'celltype' factor

        BA9_C9ALS_weighted_df_summary$label <- factor(BA9_C9ALS_weighted_df_summary$label, levels = unique(BA9_C9ALS_weighted_df_summary$label))
        BA9_C9ALS_weighted_df_summary <- subset(BA9_C9ALS_weighted_df_summary, method == "LIME")

        ggplot(BA9_C9ALS_weighted_df_summary, aes(x = n_genes, y = mean_accuracy_by_group, colour = method, group = method)) + 
        theme_bw() + 
        geom_vline(aes(xintercept = optimal_genes), linetype = "dashed", color = "black") +
        geom_point(size = 1) +
        geom_line() +
        scale_x_reverse() +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(colour = "black", angle = 90, hjust = 1, vjust = 0.5),
            axis.text.y = element_text(colour = "black"),
            axis.title.y = element_text(face = "bold"),
            axis.title.x = element_text(face = "bold"),
        ) +
        facet_grid(group ~ label, scales = "free_x") + 
        ylab("Accuracy") +
        xlab("N genes") +
        scale_colour_manual(values = c("#F18F01","#048BA8" , "#2E4057"))
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp2.pdf', height = 2, width = 15)
    ##
    
    ###########################
    ## SFTLD BA9 -- Lambda approach and plot
    ###########################
    ## code
        ## sumarize SALS
        BA9_SFTLD_weighted_df$disease_cat <- "SFTLD"
        BA9_SFTLD_weighted_df_summary <- BA9_SFTLD_weighted_df %>%
            dplyr::group_by(group, celltype, n_genes, method) %>%
            dplyr::summarise(mean_accuracy_by_group = mean(mean_accuracy, na.rm = TRUE))
        BA9_SFTLD_weighted_df_summary <- data.frame(BA9_SFTLD_weighted_df_summary)

        BA9_SFTLD_weighted_df_summary <- subset(BA9_SFTLD_weighted_df_summary, method == "LIME" )
        
        BA9_SFTLD_weighted_df_summary <- subset(BA9_SFTLD_weighted_df_summary, n_genes >= 25 )


        # Function to calculate the score with regularization
        calculate_score <- function(df, lambda) {
        df %>%
            mutate(score = mean_accuracy_by_group - lambda * n_genes)
        }

        lambda_values <- seq(0, 0.1, by = 0.00005)  # Try different values of lambda
        result <- data.frame()

        for (lambda in lambda_values) {
        df_with_score <- calculate_score(BA9_SFTLD_weighted_df_summary, lambda)
        
        # Find the optimal row for each celltype
        optimal_for_lambda <- df_with_score %>%
            group_by(celltype) %>%
            slice_max(score, n = 1)
        
        # Store the results for each lambda value
        optimal_for_lambda$lambda <- lambda
        result <- bind_rows(result, optimal_for_lambda)
        }

        head(result, n = 500)

        # View the result
        print(result)

        ## Plot
        BA9_SFTLD_weighted_df_summary <- BA9_SFTLD_weighted_df %>%
            dplyr::group_by(group, celltype, n_genes, method) %>%
            dplyr::summarise(mean_accuracy_by_group = mean(mean_accuracy, na.rm = TRUE))

        BA9_SFTLD_weighted_df_summary <- subset(BA9_SFTLD_weighted_df_summary, n_genes >= 25)


        result <- subset(result, lambda == 0.00025)
        result <- result %>% dplyr::select(celltype, n_genes)
        
        colnames(result) <- c("celltype", "optimal_genes")
        write.csv(result, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_lambda_0.00025_gene_set_BA9_SFTLD_narval.csv'))


        BA9_SFTLD_weighted_df_summary <- merge(BA9_SFTLD_weighted_df_summary, result, by = "celltype")


        BA9_SFTLD_weighted_df_summary$method <- factor(BA9_SFTLD_weighted_df_summary$method, levels = c("LIME", "random_HVG", "random_gene"))
        BA9_SFTLD_weighted_df_summary$celltype <- factor(BA9_SFTLD_weighted_df_summary$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))
        
        
        ## custom label
        BA9_SFTLD_weighted_df_summary <- BA9_SFTLD_weighted_df_summary %>%
            mutate(label = paste(celltype, "\n(", optimal_genes, " genes)", sep = ""))

        BA9_SFTLD_weighted_df_summary <- BA9_SFTLD_weighted_df_summary %>%
            arrange(celltype)  # Arrange by the 'celltype' factor

        BA9_SFTLD_weighted_df_summary$label <- factor(BA9_SFTLD_weighted_df_summary$label, levels = unique(BA9_SFTLD_weighted_df_summary$label))
        BA9_SFTLD_weighted_df_summary <- subset(BA9_SFTLD_weighted_df_summary, method == "LIME")

        ggplot(BA9_SFTLD_weighted_df_summary, aes(x = n_genes, y = mean_accuracy_by_group, colour = method, group = method)) + 
        theme_bw() + 
        geom_vline(aes(xintercept = optimal_genes), linetype = "dashed", color = "black") +
        geom_point(size = 1) +
        geom_line() +
        scale_x_reverse() +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(colour = "black", angle = 90, hjust = 1, vjust = 0.5),
            axis.text.y = element_text(colour = "black"),
            axis.title.y = element_text(face = "bold"),
            axis.title.x = element_text(face = "bold"),
        ) +
        facet_grid(group ~ label, scales = "free_x") + 
        ylab("Accuracy") +
        xlab("N genes") +
        scale_colour_manual(values = c("#F18F01","#048BA8" , "#2E4057"))
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp3.pdf', height = 2, width = 15)
    ##

    ###########################
    ## C9FTLD BA9 -- Lambda approach and plot
    ###########################
    ## code
        ## sumarize C9FTLD
        BA9_C9FTLD_weighted_df$disease_cat <- "C9FTLD"
        BA9_C9FTLD_weighted_df_summary <- BA9_C9FTLD_weighted_df %>%
            dplyr::group_by(group, celltype, n_genes, method) %>%
            dplyr::summarise(mean_accuracy_by_group = mean(mean_accuracy, na.rm = TRUE))
        BA9_C9FTLD_weighted_df_summary <- data.frame(BA9_C9FTLD_weighted_df_summary)

        BA9_C9FTLD_weighted_df_summary <- subset(BA9_C9FTLD_weighted_df_summary, method == "LIME" )
        
        BA9_C9FTLD_weighted_df_summary <- subset(BA9_C9FTLD_weighted_df_summary, n_genes >= 25 )


        # Function to calculate the score with regularization
        calculate_score <- function(df, lambda) {
        df %>%
            mutate(score = mean_accuracy_by_group - lambda * n_genes)
        }

        lambda_values <- seq(0, 0.1, by = 0.00005)  # Try different values of lambda
        result <- data.frame()

        for (lambda in lambda_values) {
        df_with_score <- calculate_score(BA9_C9FTLD_weighted_df_summary, lambda)
        
        # Find the optimal row for each celltype
        optimal_for_lambda <- df_with_score %>%
            group_by(celltype) %>%
            slice_max(score, n = 1)
        
        # Store the results for each lambda value
        optimal_for_lambda$lambda <- lambda
        result <- bind_rows(result, optimal_for_lambda)
        }

        head(result, n = 500)

        # View the result
        print(result)

        ## Plot
        BA9_C9FTLD_weighted_df_summary <- BA9_C9FTLD_weighted_df %>%
            dplyr::group_by(group, celltype, n_genes, method) %>%
            dplyr::summarise(mean_accuracy_by_group = mean(mean_accuracy, na.rm = TRUE))

        BA9_C9FTLD_weighted_df_summary <- subset(BA9_C9FTLD_weighted_df_summary, n_genes >= 25)


        result <- subset(result, lambda == 0.00025)
        result <- result %>% dplyr::select(celltype, n_genes)
        
        colnames(result) <- c("celltype", "optimal_genes")
        write.csv(result, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_lambda_0.00025_gene_set_BA9_C9FTLD_narval.csv'))

        BA9_C9FTLD_weighted_df_summary <- merge(BA9_C9FTLD_weighted_df_summary, result, by = "celltype")


        BA9_C9FTLD_weighted_df_summary$method <- factor(BA9_C9FTLD_weighted_df_summary$method, levels = c("LIME", "random_HVG", "random_gene"))
        BA9_C9FTLD_weighted_df_summary$celltype <- factor(BA9_C9FTLD_weighted_df_summary$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))
        
        
        ## custom label
        BA9_C9FTLD_weighted_df_summary <- BA9_C9FTLD_weighted_df_summary %>%
            mutate(label = paste(celltype, "\n(", optimal_genes, " genes)", sep = ""))

        BA9_C9FTLD_weighted_df_summary <- BA9_C9FTLD_weighted_df_summary %>%
            arrange(celltype)  # Arrange by the 'celltype' factor

        BA9_C9FTLD_weighted_df_summary$label <- factor(BA9_C9FTLD_weighted_df_summary$label, levels = unique(BA9_C9FTLD_weighted_df_summary$label))
        BA9_C9FTLD_weighted_df_summary <- subset(BA9_C9FTLD_weighted_df_summary, method == "LIME")

        ggplot(BA9_C9FTLD_weighted_df_summary, aes(x = n_genes, y = mean_accuracy_by_group, colour = method, group = method)) + 
        theme_bw() + 
        geom_vline(aes(xintercept = optimal_genes), linetype = "dashed", color = "black") +
        geom_point(size = 1) +
        geom_line() +
        scale_x_reverse() +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(colour = "black", angle = 90, hjust = 1, vjust = 0.5),
            axis.text.y = element_text(colour = "black"),
            axis.title.y = element_text(face = "bold"),
            axis.title.x = element_text(face = "bold"),
        ) +
        facet_grid(group ~ label, scales = "free_x") + 
        ylab("Accuracy") +
        xlab("N genes") +
        scale_colour_manual(values = c("#F18F01","#048BA8" , "#2E4057"))
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp4.pdf', height = 2, width = 15)
    ##



##

############################################################################ Print file containing optimal gene set genes for each DNN group --> we will use this in downstream analyses. 
############################################################################
############################################################################
############################################################################ 

## code
    ###########################
    ## SALS BA4
    ###########################
    ## code
        #celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",  "OPC",  "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA4"
        par_status = "SALS"
        par_prep = "CombatSeq"

        for(celltype2 in unique(celltype_list)){
            
            ## read in optimal number of genes, obtained from the above code
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_lambda_0.00025_gene_set_',par_brain_region,'_',par_status,'_narval.csv'), sep = ",")
            df_lim <- subset(df, celltype == celltype2 )

            ## read in the LIME weighted percentile rank file which was created by and used in the permutation scripts
            df_rank = read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_percentile_rank_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")

            # Extract the optimal_genes value from df_lim
            n <- df_lim$optimal_genes[1]

            # Subset the first n columns of df_rank
            df_rank_subset <- df_rank[ 1:n,]

            # View the result
            print(nrow(df_rank_subset) == n)

            write.csv(df_rank_subset, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_lambda_0.00025_optimal_genes_',par_brain_region,'_',par_status,'_',celltype2,'_narval.csv'))
            
        }
    ##

    ###########################
    ## C9ALS BA4
    ###########################
    ## code
        #celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",   "Astro",   "OPC",  "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA4"
        par_status = "C9ALS"
        par_prep = "CombatSeq"

        for(celltype2 in unique(celltype_list)){
            
            ## read in optimal number of genes, obtained from the above code
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_lambda_0.00025_gene_set_',par_brain_region,'_',par_status,'_narval.csv'), sep = ",")
            df_lim <- subset(df, celltype == celltype2 )

            ## read in the LIME weighted percentile rank file which was created by and used in the permutation scripts
            df_rank = read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_percentile_rank_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")

            # Extract the optimal_genes value from df_lim
            n <- df_lim$optimal_genes[1]

            # Subset the first n columns of df_rank
            df_rank_subset <- df_rank[ 1:n,]

            # View the result
            print(nrow(df_rank_subset) == n)

            write.csv(df_rank_subset, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_lambda_0.00025_optimal_genes_',par_brain_region,'_',par_status,'_',celltype2,'_narval.csv'))
            
        }
    ##

    ###########################
    ## SFTLD BA4
    ###########################
    ## code
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA4"
        par_status = "SFTLD"
        par_prep = "CombatSeq"

        for(celltype2 in unique(celltype_list)){
            
            ## read in optimal number of genes, obtained from the above code
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_lambda_0.00025_gene_set_',par_brain_region,'_',par_status,'_narval.csv'), sep = ",")
            df_lim <- subset(df, celltype == celltype2 )

            ## read in the LIME weighted percentile rank file which was created by and used in the permutation scripts
            df_rank = read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_percentile_rank_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")

            # Extract the optimal_genes value from df_lim
            n <- df_lim$optimal_genes[1]

            # Subset the first n columns of df_rank
            df_rank_subset <- df_rank[ 1:n,]

            # View the result
            print(nrow(df_rank_subset) == n)

            write.csv(df_rank_subset, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_lambda_0.00025_optimal_genes_',par_brain_region,'_',par_status,'_',celltype2,'_narval.csv'))
            
        }
    ##

    ###########################
    ## C9FTLD BA4
    ###########################
    ## code
        #celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",  "Astro",   "OPC",  "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA4"
        par_status = "C9FTLD"
        par_prep = "CombatSeq"

        for(celltype2 in unique(celltype_list)){
            
            ## read in optimal number of genes, obtained from the above code
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_lambda_0.00025_gene_set_',par_brain_region,'_',par_status,'_narval.csv'), sep = ",")
            df_lim <- subset(df, celltype == celltype2 )

            ## read in the LIME weighted percentile rank file which was created by and used in the permutation scripts
            df_rank = read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_percentile_rank_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")

            # Extract the optimal_genes value from df_lim
            n <- df_lim$optimal_genes[1]

            # Subset the first n columns of df_rank
            df_rank_subset <- df_rank[ 1:n,]

            # View the result
            print(nrow(df_rank_subset) == n)

            write.csv(df_rank_subset, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_lambda_0.00025_optimal_genes_',par_brain_region,'_',par_status,'_',celltype2,'_narval.csv'))
            
        }
    ##

    ###########################
    ## SALS BA9
    ###########################
    ## code
        #celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")
        celltype_list = c("L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo", "OPC",  "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA9"
        par_status = "SALS"
        par_prep = "CombatSeq"

        for(celltype2 in unique(celltype_list)){
            
            ## read in optimal number of genes, obtained from the above code
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_lambda_0.00025_gene_set_',par_brain_region,'_',par_status,'_narval.csv'), sep = ",")
            df_lim <- subset(df, celltype == celltype2 )

            ## read in the LIME weighted percentile rank file which was created by and used in the permutation scripts
            df_rank = read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_percentile_rank_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")

            # Extract the optimal_genes value from df_lim
            n <- df_lim$optimal_genes[1]

            # Subset the first n columns of df_rank
            df_rank_subset <- df_rank[ 1:n,]

            # View the result
            print(nrow(df_rank_subset) == n)

            write.csv(df_rank_subset, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_lambda_0.00025_optimal_genes_',par_brain_region,'_',par_status,'_',celltype2,'_narval.csv'))
            
        }
    ##

    ###########################
    ## C9ALS BA9
    ###########################
    ## code
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA9"
        par_status = "C9ALS"
        par_prep = "CombatSeq"

        for(celltype2 in unique(celltype_list)){
            
            ## read in optimal number of genes, obtained from the above code
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_lambda_0.00025_gene_set_',par_brain_region,'_',par_status,'_narval.csv'), sep = ",")
            df_lim <- subset(df, celltype == celltype2 )

            ## read in the LIME weighted percentile rank file which was created by and used in the permutation scripts
            df_rank = read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_percentile_rank_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")

            # Extract the optimal_genes value from df_lim
            n <- df_lim$optimal_genes[1]

            # Subset the first n columns of df_rank
            df_rank_subset <- df_rank[ 1:n,]

            # View the result
            print(nrow(df_rank_subset) == n)

            write.csv(df_rank_subset, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_lambda_0.00025_optimal_genes_',par_brain_region,'_',par_status,'_',celltype2,'_narval.csv'))
            
        }
    ##

    ###########################
    ## SFTLD BA9
    ###########################
    ## code
        #celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",  "Astro",   "OPC",  "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA9"
        par_status = "SFTLD"
        par_prep = "CombatSeq"

        for(celltype2 in unique(celltype_list)){
            
            ## read in optimal number of genes, obtained from the above code
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_lambda_0.00025_gene_set_',par_brain_region,'_',par_status,'_narval.csv'), sep = ",")
            df_lim <- subset(df, celltype == celltype2 )

            ## read in the LIME weighted percentile rank file which was created by and used in the permutation scripts
            df_rank = read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_percentile_rank_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")

            # Extract the optimal_genes value from df_lim
            n <- df_lim$optimal_genes[1]

            # Subset the first n columns of df_rank
            df_rank_subset <- df_rank[ 1:n,]

            # View the result
            print(nrow(df_rank_subset) == n)

            write.csv(df_rank_subset, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_lambda_0.00025_optimal_genes_',par_brain_region,'_',par_status,'_',celltype2,'_narval.csv'))
            
        }
    ##

    ###########################
    ## C9FTLD BA9
    ###########################
    ## code
        #celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",   "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA9"
        par_status = "C9FTLD"
        par_prep = "CombatSeq"

        for(celltype2 in unique(celltype_list)){
            
            ## read in optimal number of genes, obtained from the above code
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_lambda_0.00025_gene_set_',par_brain_region,'_',par_status,'_narval.csv'), sep = ",")
            df_lim <- subset(df, celltype == celltype2 )

            ## read in the LIME weighted percentile rank file which was created by and used in the permutation scripts
            df_rank = read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_percentile_rank_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")

            # Extract the optimal_genes value from df_lim
            n <- df_lim$optimal_genes[1]

            # Subset the first n columns of df_rank
            df_rank_subset <- df_rank[ 1:n,]

            # View the result
            print(nrow(df_rank_subset) == n)

            write.csv(df_rank_subset, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_lambda_0.00025_optimal_genes_',par_brain_region,'_',par_status,'_',celltype2,'_narval.csv'))
            
        }
    ##
##
    
############################################################################ Figure 2 -- example workflow
############################################################################
############################################################################
############################################################################ 

## code Z-score distribution: BA4 SALS L3_5

    ###########################
    ## SALS BA4
    ###########################
    ## code
        celltype_list = c("L3_L5")

        par_brain_region = "BA4"
        par_status = "SALS"
        par_prep = "CombatSeq"
        
        ################################################
        ## LIME outputs weighted by percent expression
        ## This function subsets the weighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_weighted_across_donor = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_weighted_across_donor, celltype)

        
        for(celltype2 in unique(celltype_list)){
            
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_weighted_across_donor <- as.numeric(df$Mean_z_score_weighted_across_donor)
            df <- df[order( -(df$Mean_z_score_weighted_across_donor)),]
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_weighted_across_donor)
            df$celltype <- celltype2
            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")
        BA4_SALS_weighted_df <- fill

        BA4_SALS_weighted_df$disease_cat <- "SALS"

        BA4_weighted_distribution <- rbind(BA4_SALS_weighted_df)

        ggplot(BA4_weighted_distribution, aes(x = Mean_z_score_weighted_across_donor, colour = disease_cat)) + 
        geom_vline(xintercept = 1.0, linetype = "dashed", colour = "darkgrey") +
        theme_bw() + 
        geom_density() +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(colour = "black"),
            axis.text.y = element_text(colour = "black"),
            axis.title.y = element_text(face = "bold"),
            axis.title.x = element_text(face = "bold"),
        ) +
        ylab("Density") + xlab("LIME feature importance\nZ-score") +
        scale_colour_manual(values = c("black")) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 2, width = 2.5)
    ##

##

## code permutation: BA4 SALS L3_5
    ###########################
    ## SALS BA4
    ###########################
    ## code
        celltype_list = c("L3_L5")

        par_brain_region = "BA4"
        par_status = "SALS"
        par_prep = "CombatSeq"
        
        ################################################
        ## LIME outputs weighted by percent expression
        ## This function subsets the weighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        region = "fill"
        group = "fill"
        celltype = "fill"
        n_genes = 0
        test_accuracy = 0
        mean_accuracy = 0
        method = "fill"
        fill <- data.frame(region,group, celltype, n_genes, test_accuracy, mean_accuracy, method)
        

        for(celltype2 in unique(celltype_list)){
        
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/permutation_weighted_report_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df$five_fold_accuracies <- strsplit(as.character(df$five_fold_accuracies), "; ")
            df$five_fold_accuracies <- lapply(df$five_fold_accuracies, as.numeric)
            df <- cbind(df, do.call(rbind, df$five_fold_accuracies))
            colnames(df)[13:17] <- paste0("fold_", 1:5)
            df$mean_accuracy <- rowMeans(df[, c("test_accuracy", "fold_1", "fold_2", "fold_3", "fold_4", "fold_5")], na.rm = TRUE)
            df <- df %>% dplyr::select(region, group, celltype, n_genes, test_accuracy, mean_accuracy, method)
            fill <- rbind(fill, df)

        }

        fill <- subset(fill, region != "fill")
        BA4_SALS_weighted_df <- fill


        ## sumarize SALS
        BA4_SALS_weighted_df$disease_cat <- "SALS"
        
        BA4_SALS_weighted_df_summary <- BA4_SALS_weighted_df %>%
            dplyr::group_by(group, celltype, n_genes, method) %>%
            dplyr::summarise(mean_accuracy_by_group = mean(mean_accuracy, na.rm = TRUE))
        
        BA4_SALS_weighted_df_summary <- data.frame(BA4_SALS_weighted_df_summary)

        BA4_SALS_weighted_df_summary <- subset(BA4_SALS_weighted_df_summary, method == "LIME" )
        
        BA4_SALS_weighted_df_summary <- subset(BA4_SALS_weighted_df_summary, n_genes >= 25 )


        # Function to calculate the score with regularization
        calculate_score <- function(df, lambda) {
        df %>%
            mutate(score = mean_accuracy_by_group - lambda * n_genes)
        }

        lambda_values <- seq(0, 0.1, by = 0.00005)  # Try different values of lambda
        result <- data.frame()

        for (lambda in lambda_values) {
        df_with_score <- calculate_score(BA4_SALS_weighted_df_summary, lambda)
        
        # Find the optimal row for each celltype
        optimal_for_lambda <- df_with_score %>%
            group_by(celltype) %>%
            slice_max(score, n = 1)
        
        # Store the results for each lambda value
        optimal_for_lambda$lambda <- lambda
        result <- bind_rows(result, optimal_for_lambda)
        }

        head(result, n = 500)

        # View the result
        print(result)

        ## Plot
        BA4_SALS_weighted_df_summary <- BA4_SALS_weighted_df %>%
            dplyr::group_by(group, celltype, n_genes, method) %>%
            dplyr::summarise(mean_accuracy_by_group = mean(mean_accuracy, na.rm = TRUE))

        BA4_SALS_weighted_df_summary <- subset(BA4_SALS_weighted_df_summary, n_genes >= 25)


        result <- subset(result, lambda == 0.00025)
        result <- result %>% dplyr::select(celltype, n_genes)
        
        colnames(result) <- c("celltype", "optimal_genes")
        
        BA4_SALS_weighted_df_summary <- merge(BA4_SALS_weighted_df_summary, result, by = "celltype")
        BA4_SALS_weighted_df_summary$method <- factor(BA4_SALS_weighted_df_summary$method, levels = c("LIME", "random_HVG", "random_gene"))
        BA4_SALS_weighted_df_summary$celltype <- factor(BA4_SALS_weighted_df_summary$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))
        
        
        ## custom label
        BA4_SALS_weighted_df_summary <- BA4_SALS_weighted_df_summary %>%
            mutate(label = paste(celltype, "\n(", optimal_genes, " genes)", sep = ""))

        BA4_SALS_weighted_df_summary <- BA4_SALS_weighted_df_summary %>%
            arrange(celltype)  # Arrange by the 'celltype' factor

        BA4_SALS_weighted_df_summary$label <- factor(BA4_SALS_weighted_df_summary$label, levels = unique(BA4_SALS_weighted_df_summary$label))

        BA4_SALS_weighted_df_summary <- subset(BA4_SALS_weighted_df_summary, method != "random_HVG")
        
        ggplot(BA4_SALS_weighted_df_summary, aes(x = n_genes, y = mean_accuracy_by_group, colour = method, group = method)) + 
        theme_bw() + 
        geom_vline(aes(xintercept = optimal_genes), linetype = "dashed", color = "black") +
        geom_point() +
        geom_line() +
        scale_x_reverse() +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(colour = "black"),
            axis.text.y = element_text(colour = "black"),
            axis.title.y = element_text(face = "bold"),
            axis.title.x = element_text(face = "bold"),
        ) +
        ylab("Accuracy") +
        xlab("N genes") +
        scale_colour_manual(values = c("#F18F01", "#2E4057"))
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 2, width = 2.5)


    ##
##

## I think we will use L4_L5 for now.

## code Z-score distribution: BA4 SALS L4_L5

    ###########################
    ## SALS BA4
    ###########################
    ## code
        celltype_list = c("L4_L5")

        par_brain_region = "BA4"
        par_status = "SALS"
        par_prep = "CombatSeq"
        
        ################################################
        ## LIME outputs weighted by percent expression
        ## This function subsets the weighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_weighted_across_donor = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_weighted_across_donor, celltype)

        
        for(celltype2 in unique(celltype_list)){
            
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_weighted_across_donor <- as.numeric(df$Mean_z_score_weighted_across_donor)
            df <- df[order( -(df$Mean_z_score_weighted_across_donor)),]
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_weighted_across_donor)
            df$celltype <- celltype2
            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")
        BA4_SALS_weighted_df <- fill

        BA4_SALS_weighted_df$disease_cat <- "SALS"

        BA4_weighted_distribution <- rbind(BA4_SALS_weighted_df)

        ggplot(BA4_weighted_distribution, aes(x = Mean_z_score_weighted_across_donor, colour = disease_cat)) + 
        geom_vline(xintercept = 1.0, linetype = "dashed", colour = "darkgrey") +
        theme_bw() + 
        geom_density() +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(colour = "black"),
            axis.text.y = element_text(colour = "black"),
            axis.title.y = element_text(face = "bold"),
            axis.title.x = element_text(face = "bold"),
        ) +
        ylab("Density") + xlab("LIME feature importance\nZ-score") +
        scale_colour_manual(values = c("black")) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 2, width = 2.5)
    ##

##

## code permutation: BA4 SALS L4_L5
    ###########################
    ## SALS BA4
    ###########################
    ## code
        celltype_list = c("L4_L5")

        par_brain_region = "BA4"
        par_status = "SALS"
        par_prep = "CombatSeq"
        
        ################################################
        ## LIME outputs weighted by percent expression
        ## This function subsets the weighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        region = "fill"
        group = "fill"
        celltype = "fill"
        n_genes = 0
        test_accuracy = 0
        mean_accuracy = 0
        method = "fill"
        fill <- data.frame(region,group, celltype, n_genes, test_accuracy, mean_accuracy, method)
        

        for(celltype2 in unique(celltype_list)){
        
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/permutation_weighted_report_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df$five_fold_accuracies <- strsplit(as.character(df$five_fold_accuracies), "; ")
            df$five_fold_accuracies <- lapply(df$five_fold_accuracies, as.numeric)
            df <- cbind(df, do.call(rbind, df$five_fold_accuracies))
            colnames(df)[13:17] <- paste0("fold_", 1:5)
            df$mean_accuracy <- rowMeans(df[, c("test_accuracy", "fold_1", "fold_2", "fold_3", "fold_4", "fold_5")], na.rm = TRUE)
            df <- df %>% dplyr::select(region, group, celltype, n_genes, test_accuracy, mean_accuracy, method)
            fill <- rbind(fill, df)

        }

        fill <- subset(fill, region != "fill")
        BA4_SALS_weighted_df <- fill


        ## sumarize SALS
        BA4_SALS_weighted_df$disease_cat <- "SALS"
        
        BA4_SALS_weighted_df_summary <- BA4_SALS_weighted_df %>%
            dplyr::group_by(group, celltype, n_genes, method) %>%
            dplyr::summarise(mean_accuracy_by_group = mean(mean_accuracy, na.rm = TRUE))
        
        BA4_SALS_weighted_df_summary <- data.frame(BA4_SALS_weighted_df_summary)

        BA4_SALS_weighted_df_summary <- subset(BA4_SALS_weighted_df_summary, method == "LIME" )
        
        BA4_SALS_weighted_df_summary <- subset(BA4_SALS_weighted_df_summary, n_genes >= 25 )


        # Function to calculate the score with regularization
        calculate_score <- function(df, lambda) {
        df %>%
            mutate(score = mean_accuracy_by_group - lambda * n_genes)
        }

        lambda_values <- seq(0, 0.1, by = 0.00005)  # Try different values of lambda
        result <- data.frame()

        for (lambda in lambda_values) {
        df_with_score <- calculate_score(BA4_SALS_weighted_df_summary, lambda)
        
        # Find the optimal row for each celltype
        optimal_for_lambda <- df_with_score %>%
            group_by(celltype) %>%
            slice_max(score, n = 1)
        
        # Store the results for each lambda value
        optimal_for_lambda$lambda <- lambda
        result <- bind_rows(result, optimal_for_lambda)
        }

        head(result, n = 500)

        # View the result
        print(result)

        ## Plot
        BA4_SALS_weighted_df_summary <- BA4_SALS_weighted_df %>%
            dplyr::group_by(group, celltype, n_genes, method) %>%
            dplyr::summarise(mean_accuracy_by_group = mean(mean_accuracy, na.rm = TRUE))

        BA4_SALS_weighted_df_summary <- subset(BA4_SALS_weighted_df_summary, n_genes >= 25)


        result <- subset(result, lambda == 0.00025)
        result <- result %>% dplyr::select(celltype, n_genes)
        
        colnames(result) <- c("celltype", "optimal_genes")
        
        BA4_SALS_weighted_df_summary <- merge(BA4_SALS_weighted_df_summary, result, by = "celltype")
        BA4_SALS_weighted_df_summary$method <- factor(BA4_SALS_weighted_df_summary$method, levels = c("LIME", "random_HVG", "random_gene"))
        BA4_SALS_weighted_df_summary$celltype <- factor(BA4_SALS_weighted_df_summary$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))
        
        
        ## custom label
        BA4_SALS_weighted_df_summary <- BA4_SALS_weighted_df_summary %>%
            mutate(label = paste(celltype, "\n(", optimal_genes, " genes)", sep = ""))

        BA4_SALS_weighted_df_summary <- BA4_SALS_weighted_df_summary %>%
            arrange(celltype)  # Arrange by the 'celltype' factor

        BA4_SALS_weighted_df_summary$label <- factor(BA4_SALS_weighted_df_summary$label, levels = unique(BA4_SALS_weighted_df_summary$label))

        BA4_SALS_weighted_df_summary <- subset(BA4_SALS_weighted_df_summary, method == "LIME")
        
        ggplot(BA4_SALS_weighted_df_summary, aes(x = n_genes, y = mean_accuracy_by_group, colour = method, group = method)) + 
        theme_bw() + 
        geom_vline(aes(xintercept = optimal_genes), linetype = "dashed", color = "black") +
        geom_point() +
        geom_line() +
        scale_x_reverse() +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(colour = "black"),
            axis.text.y = element_text(colour = "black"),
            axis.title.y = element_text(face = "bold"),
            axis.title.x = element_text(face = "bold"),
        ) +
        #facet_grid(~label)+
        ylab("Accuracy") +
        xlab("N genes") +
        scale_y_continuous(limits = c(0.7,0.95)) +
        scale_colour_manual(values = c("#F18F01", "#2E4057"))
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 2, width = 2.5)


    ##
##


############################################################################ Figure 2 -- Optimal LIME vs Random gene 5CV (need to modify the code)
############################################################################
############################################################################
############################################################################ 


## code optimal LIME vs random gene accuracy
    ###########################
    ## SALS BA4
    ###########################
    ## code
        #celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM", "OPC",  "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA4"
        par_status = "SALS"
        par_prep = "CombatSeq"
        
        ################################################
        ## LIME outputs weighted by percent expression
        ## This function subsets the weighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        region = "fill"
        group = "fill"
        celltype = "fill"
        n_genes = 0
        #mean_accuracy = 0
        fold_1 = 0
        fold_2 = 0
        fold_3 = 0
        fold_4 = 0
        fold_5 = 0
        method = "fill"
        fill <- data.frame(region,group, celltype, n_genes, fold_1, fold_2, fold_3, fold_4, fold_5, method)
        
        for(celltype2 in unique(celltype_list)){
        
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_report_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df$method[(nrow(df)-4):nrow(df)] <- "random_gene"
            df$five_fold_accuracies <- strsplit(as.character(df$five_fold_accuracies), "; ")
            df$five_fold_accuracies <- lapply(df$five_fold_accuracies, as.numeric)
            df <- cbind(df, do.call(rbind, df$five_fold_accuracies))
            colnames(df)[12:16] <- paste0("fold_", 1:5)
            #df$mean_accuracy <- rowMeans(df[, c("test_accuracy", "fold_1", "fold_2", "fold_3", "fold_4", "fold_5")], na.rm = TRUE)
            df <- df %>% dplyr::select(region, group, celltype, n_genes, fold_1, fold_2, fold_3, fold_4, fold_5, method)
            fill <- rbind(fill, df)

        }

        fill <- subset(fill, region != "fill")
        BA4_SALS_weighted_df <- fill

        ## only keep LIME and random gene
        BA4_SALS_weighted_df_lim <- subset(BA4_SALS_weighted_df, method %in% c("LIME", "random_gene") )

        ## read in the optimal file
        df <- read.delim('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_lambda_0.00025_gene_set_BA4_SALS_narval.csv', sep = ",")

        ## Only keep optimal
        BA4_SALS_weighted_df_lim <- BA4_SALS_weighted_df_lim %>%
            left_join(df, by = "celltype")

        ## Filter the rows where n_genes equals optimal_genes
        BA4_SALS_weighted_df_lim <- BA4_SALS_weighted_df_lim %>%
            filter(n_genes == optimal_genes)
        
        ## place in long format
        df_long <- BA4_SALS_weighted_df_lim %>%
            pivot_longer(cols = c(fold_1, fold_2, fold_3, fold_4, fold_5),
                        names_to = "fold",  # New column to hold the fold names
                        values_to = "accuracy") %>% as.data.frame()  # New column to hold the corresponding values
        
        ## Calculate mean and SD
        df_stats <- df_long %>%
            group_by(celltype, method) %>%  # Group by celltype and method
            summarize(
                mean_accuracy = mean(accuracy, na.rm = TRUE),  # Calculate mean accuracy
                sd_accuracy = sd(accuracy, na.rm = TRUE)      # Calculate standard deviation of accuracy
            ) %>% as.data.frame()
        
        df_stats <- df_stats %>%
            left_join(df, by = "celltype")

        ## set factor level
        df_stats$celltype <- factor(df_stats$celltype, levels = rev(c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  )))
        df_stats$optimal_genes[df_stats$method == "random_gene"] <- NA


        ## plot
        ggplot(df_stats, aes(x = mean_accuracy, y = celltype, fill = method, label = optimal_genes)) + 
        theme_bw() +
        geom_text(aes(x = 1.15), size = 3)+
        geom_vline(xintercept = 1.05) +
        geom_vline(xintercept = 0.90, colour = "lightgrey", linetype = 'dotted') +
        geom_bar(stat = 'identity', position = "dodge") +
        geom_errorbar(aes(xmin=mean_accuracy-sd_accuracy, xmax=mean_accuracy+sd_accuracy), width=.2,
                position=position_dodge(.9)) +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(colour = "black"),
            axis.text.y = element_text(colour = "black"),
            axis.title.y = element_blank(),
            axis.title.x = element_text(face = "bold")
        ) +
        xlab("Accuracy") +
        scale_fill_manual(values = c("#F18F01", "#2E4057")) +
        coord_cartesian(xlim=c(0.5,1.2)) +
        scale_x_continuous(
            breaks = c(0.5, 0.75, 1.00),          # Set the tick positions
            labels = c("0.5", "0.75", "1.0")      # Optional: Customize the tick labels
        )
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 3, width = 2.25)
    ##

    ###########################
    ## C9ALS BA4
    ###########################
    ## code
        #celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM", "Astro",   "OPC",  "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA4"
        par_status = "C9ALS"
        par_prep = "CombatSeq"
        
        ################################################
        ## LIME outputs weighted by percent expression
        ## This function subsets the weighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        region = "fill"
        group = "fill"
        celltype = "fill"
        n_genes = 0
        #mean_accuracy = 0
        fold_1 = 0
        fold_2 = 0
        fold_3 = 0
        fold_4 = 0
        fold_5 = 0
        method = "fill"
        fill <- data.frame(region,group, celltype, n_genes, fold_1, fold_2, fold_3, fold_4, fold_5, method)
        
        for(celltype2 in unique(celltype_list)){
        
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_report_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df$method[(nrow(df)-4):nrow(df)] <- "random_gene"
            df$five_fold_accuracies <- strsplit(as.character(df$five_fold_accuracies), "; ")
            df$five_fold_accuracies <- lapply(df$five_fold_accuracies, as.numeric)
            df <- cbind(df, do.call(rbind, df$five_fold_accuracies))
            colnames(df)[12:16] <- paste0("fold_", 1:5)
            #df$mean_accuracy <- rowMeans(df[, c("test_accuracy", "fold_1", "fold_2", "fold_3", "fold_4", "fold_5")], na.rm = TRUE)
            df <- df %>% dplyr::select(region, group, celltype, n_genes, fold_1, fold_2, fold_3, fold_4, fold_5, method)
            fill <- rbind(fill, df)

        }

        fill <- subset(fill, region != "fill")
        BA4_C9ALS_weighted_df <- fill

        ## only keep LIME and random gene
        BA4_C9ALS_weighted_df_lim <- subset(BA4_C9ALS_weighted_df, method %in% c("LIME", "random_gene") )

        ## read in the optimal file
        df <- read.delim('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_lambda_0.00025_gene_set_BA4_C9ALS_narval.csv', sep = ",")

        ## Only keep optimal
        BA4_C9ALS_weighted_df_lim <- BA4_C9ALS_weighted_df_lim %>%
            left_join(df, by = "celltype")

        ## Filter the rows where n_genes equals optimal_genes
        BA4_C9ALS_weighted_df_lim <- BA4_C9ALS_weighted_df_lim %>%
            filter(n_genes == optimal_genes)
        
        ## place in long format
        df_long <- BA4_C9ALS_weighted_df_lim %>%
            pivot_longer(cols = c(fold_1, fold_2, fold_3, fold_4, fold_5),
                        names_to = "fold",  # New column to hold the fold names
                        values_to = "accuracy") %>% as.data.frame()  # New column to hold the corresponding values
        
        ## Calculate mean and SD
        df_stats <- df_long %>%
            group_by(celltype, method) %>%  # Group by celltype and method
            summarize(
                mean_accuracy = mean(accuracy, na.rm = TRUE),  # Calculate mean accuracy
                sd_accuracy = sd(accuracy, na.rm = TRUE)      # Calculate standard deviation of accuracy
            ) %>% as.data.frame()
        
        df_stats <- df_stats %>%
            left_join(df, by = "celltype")

        ## set factor level
        df_stats$celltype <- factor(df_stats$celltype, levels = rev(c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  )))
        df_stats$optimal_genes[df_stats$method == "random_gene"] <- NA


        ## plot
        ggplot(df_stats, aes(x = mean_accuracy, y = celltype, fill = method, label = optimal_genes)) + 
        theme_bw() +
        geom_text(aes(x = 1.15), size = 3)+
        geom_vline(xintercept = 1.05) +
        geom_vline(xintercept = 0.90, colour = "lightgrey", linetype = 'dotted') +
        geom_bar(stat = 'identity', position = "dodge") +
        geom_errorbar(aes(xmin=mean_accuracy-sd_accuracy, xmax=mean_accuracy+sd_accuracy), width=.2,
                position=position_dodge(.9)) +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(colour = "black"),
            axis.text.y = element_text(colour = "black"),
            axis.title.y = element_blank(),
            axis.title.x = element_text(face = "bold")
        ) +
        xlab("Accuracy") +
        scale_fill_manual(values = c("#F18F01", "#2E4057")) +
        coord_cartesian(xlim=c(0.5,1.2)) +
        scale_x_continuous(
            breaks = c(0.5, 0.75, 1.00),          # Set the tick positions
            labels = c("0.5", "0.75", "1.0")      # Optional: Customize the tick labels
        )
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 3, width = 2.25)
    ##

    ###########################
    ## SFTLD BA4
    ###########################
    ## code
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA4"
        par_status = "SFTLD"
        par_prep = "CombatSeq"
        
        ################################################
        ## LIME outputs weighted by percent expression
        ## This function subsets the weighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        region = "fill"
        group = "fill"
        celltype = "fill"
        n_genes = 0
        #mean_accuracy = 0
        fold_1 = 0
        fold_2 = 0
        fold_3 = 0
        fold_4 = 0
        fold_5 = 0
        method = "fill"
        fill <- data.frame(region,group, celltype, n_genes, fold_1, fold_2, fold_3, fold_4, fold_5, method)
        
        for(celltype2 in unique(celltype_list)){
        
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_report_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df$method[(nrow(df)-4):nrow(df)] <- "random_gene"
            df$five_fold_accuracies <- strsplit(as.character(df$five_fold_accuracies), "; ")
            df$five_fold_accuracies <- lapply(df$five_fold_accuracies, as.numeric)
            df <- cbind(df, do.call(rbind, df$five_fold_accuracies))
            colnames(df)[12:16] <- paste0("fold_", 1:5)
            #df$mean_accuracy <- rowMeans(df[, c("test_accuracy", "fold_1", "fold_2", "fold_3", "fold_4", "fold_5")], na.rm = TRUE)
            df <- df %>% dplyr::select(region, group, celltype, n_genes, fold_1, fold_2, fold_3, fold_4, fold_5, method)
            fill <- rbind(fill, df)

        }

        fill <- subset(fill, region != "fill")
        BA4_SFTLD_weighted_df <- fill

        ## only keep LIME and random gene
        BA4_SFTLD_weighted_df_lim <- subset(BA4_SFTLD_weighted_df, method %in% c("LIME", "random_gene") )

        ## read in the optimal file
        df <- read.delim('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_lambda_0.00025_gene_set_BA4_SFTLD_narval.csv', sep = ",")

        ## Only keep optimal
        BA4_SFTLD_weighted_df_lim <- BA4_SFTLD_weighted_df_lim %>%
            left_join(df, by = "celltype")

        ## Filter the rows where n_genes equals optimal_genes
        BA4_SFTLD_weighted_df_lim <- BA4_SFTLD_weighted_df_lim %>%
            filter(n_genes == optimal_genes)
        
        ## place in long format
        df_long <- BA4_SFTLD_weighted_df_lim %>%
            pivot_longer(cols = c(fold_1, fold_2, fold_3, fold_4, fold_5),
                        names_to = "fold",  # New column to hold the fold names
                        values_to = "accuracy") %>% as.data.frame()  # New column to hold the corresponding values
        
        ## Calculate mean and SD
        df_stats <- df_long %>%
            group_by(celltype, method) %>%  # Group by celltype and method
            summarize(
                mean_accuracy = mean(accuracy, na.rm = TRUE),  # Calculate mean accuracy
                sd_accuracy = sd(accuracy, na.rm = TRUE)      # Calculate standard deviation of accuracy
            ) %>% as.data.frame()
        
        df_stats <- df_stats %>%
            left_join(df, by = "celltype")

        ## set factor level
        df_stats$celltype <- factor(df_stats$celltype, levels = rev(c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  )))
        df_stats$optimal_genes[df_stats$method == "random_gene"] <- NA


        ## plot
        ggplot(df_stats, aes(x = mean_accuracy, y = celltype, fill = method, label = optimal_genes)) + 
        theme_bw() +
        geom_text(aes(x = 1.15), size = 3)+
        geom_vline(xintercept = 1.05) +
        geom_vline(xintercept = 0.90, colour = "lightgrey", linetype = 'dotted') +
        geom_bar(stat = 'identity', position = "dodge") +
        geom_errorbar(aes(xmin=mean_accuracy-sd_accuracy, xmax=mean_accuracy+sd_accuracy), width=.2,
                position=position_dodge(.9)) +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(colour = "black"),
            axis.text.y = element_text(colour = "black"),
            axis.title.y = element_blank(),
            axis.title.x = element_text(face = "bold")
        ) +
        xlab("Accuracy") +
        scale_fill_manual(values = c("#F18F01", "#2E4057")) +
        coord_cartesian(xlim=c(0.5,1.2)) +
        scale_x_continuous(
            breaks = c(0.5, 0.75, 1.00),          # Set the tick positions
            labels = c("0.5", "0.75", "1.0")      # Optional: Customize the tick labels
        )
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 3, width = 2.25)
    ##

    ###########################
    ## C9FTLD BA4
    ###########################
    ## code
        #celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",   "Astro",   "OPC",  "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA4"
        par_status = "C9FTLD"
        par_prep = "CombatSeq"
        
        ################################################
        ## LIME outputs weighted by percent expression
        ## This function subsets the weighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        region = "fill"
        group = "fill"
        celltype = "fill"
        n_genes = 0
        #mean_accuracy = 0
        fold_1 = 0
        fold_2 = 0
        fold_3 = 0
        fold_4 = 0
        fold_5 = 0
        method = "fill"
        fill <- data.frame(region,group, celltype, n_genes, fold_1, fold_2, fold_3, fold_4, fold_5, method)
        
        for(celltype2 in unique(celltype_list)){
        
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_report_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df$method[(nrow(df)-4):nrow(df)] <- "random_gene"
            df$five_fold_accuracies <- strsplit(as.character(df$five_fold_accuracies), "; ")
            df$five_fold_accuracies <- lapply(df$five_fold_accuracies, as.numeric)
            df <- cbind(df, do.call(rbind, df$five_fold_accuracies))
            colnames(df)[12:16] <- paste0("fold_", 1:5)
            #df$mean_accuracy <- rowMeans(df[, c("test_accuracy", "fold_1", "fold_2", "fold_3", "fold_4", "fold_5")], na.rm = TRUE)
            df <- df %>% dplyr::select(region, group, celltype, n_genes, fold_1, fold_2, fold_3, fold_4, fold_5, method)
            fill <- rbind(fill, df)

        }

        fill <- subset(fill, region != "fill")
        BA4_C9FTLD_weighted_df <- fill

        ## only keep LIME and random gene
        BA4_C9FTLD_weighted_df_lim <- subset(BA4_C9FTLD_weighted_df, method %in% c("LIME", "random_gene") )

        ## read in the optimal file
        df <- read.delim('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_lambda_0.00025_gene_set_BA4_C9FTLD_narval.csv', sep = ",")

        ## Only keep optimal
        BA4_C9FTLD_weighted_df_lim <- BA4_C9FTLD_weighted_df_lim %>%
            left_join(df, by = "celltype")

        ## Filter the rows where n_genes equals optimal_genes
        BA4_C9FTLD_weighted_df_lim <- BA4_C9FTLD_weighted_df_lim %>%
            filter(n_genes == optimal_genes)
        
        ## place in long format
        df_long <- BA4_C9FTLD_weighted_df_lim %>%
            pivot_longer(cols = c(fold_1, fold_2, fold_3, fold_4, fold_5),
                        names_to = "fold",  # New column to hold the fold names
                        values_to = "accuracy") %>% as.data.frame()  # New column to hold the corresponding values
        
        ## Calculate mean and SD
        df_stats <- df_long %>%
            group_by(celltype, method) %>%  # Group by celltype and method
            summarize(
                mean_accuracy = mean(accuracy, na.rm = TRUE),  # Calculate mean accuracy
                sd_accuracy = sd(accuracy, na.rm = TRUE)      # Calculate standard deviation of accuracy
            ) %>% as.data.frame()
        
        df_stats <- df_stats %>%
            left_join(df, by = "celltype")

        ## set factor level
        df_stats$celltype <- factor(df_stats$celltype, levels = rev(c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  )))
        df_stats$optimal_genes[df_stats$method == "random_gene"] <- NA


        ## plot
        ggplot(df_stats, aes(x = mean_accuracy, y = celltype, fill = method, label = optimal_genes)) + 
        theme_bw() +
        geom_text(aes(x = 1.15), size = 3)+
        geom_vline(xintercept = 1.05) +
        geom_vline(xintercept = 0.90, colour = "lightgrey", linetype = 'dotted') +
        geom_bar(stat = 'identity', position = "dodge") +
        geom_errorbar(aes(xmin=mean_accuracy-sd_accuracy, xmax=mean_accuracy+sd_accuracy), width=.2,
                position=position_dodge(.9)) +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(colour = "black"),
            axis.text.y = element_text(colour = "black"),
            axis.title.y = element_blank(),
            axis.title.x = element_text(face = "bold")
        ) +
        xlab("Accuracy") +
        scale_fill_manual(values = c("#F18F01", "#2E4057")) +
        coord_cartesian(xlim=c(0.5,1.2)) +
        scale_x_continuous(
            breaks = c(0.5, 0.75, 1.00),          # Set the tick positions
            labels = c("0.5", "0.75", "1.0")      # Optional: Customize the tick labels
        )
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 3, width = 2.25)
    ##

    ###########################
    ## SALS BA9
    ###########################
    ## code
        #celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")
        celltype_list = c( "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM", "OPC",  "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA9"
        par_status = "SALS"
        par_prep = "CombatSeq"
        
        ################################################
        ## LIME outputs weighted by percent expression
        ## This function subsets the weighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        region = "fill"
        group = "fill"
        celltype = "fill"
        n_genes = 0
        #mean_accuracy = 0
        fold_1 = 0
        fold_2 = 0
        fold_3 = 0
        fold_4 = 0
        fold_5 = 0
        method = "fill"
        fill <- data.frame(region,group, celltype, n_genes, fold_1, fold_2, fold_3, fold_4, fold_5, method)
        
        for(celltype2 in unique(celltype_list)){
        
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_report_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df$method[(nrow(df)-4):nrow(df)] <- "random_gene"
            df$five_fold_accuracies <- strsplit(as.character(df$five_fold_accuracies), "; ")
            df$five_fold_accuracies <- lapply(df$five_fold_accuracies, as.numeric)
            df <- cbind(df, do.call(rbind, df$five_fold_accuracies))
            colnames(df)[12:16] <- paste0("fold_", 1:5)
            #df$mean_accuracy <- rowMeans(df[, c("test_accuracy", "fold_1", "fold_2", "fold_3", "fold_4", "fold_5")], na.rm = TRUE)
            df <- df %>% dplyr::select(region, group, celltype, n_genes, fold_1, fold_2, fold_3, fold_4, fold_5, method)
            fill <- rbind(fill, df)

        }

        fill <- subset(fill, region != "fill")
        BA9_SALS_weighted_df <- fill

        ## only keep LIME and random gene
        BA9_SALS_weighted_df_lim <- subset(BA9_SALS_weighted_df, method %in% c("LIME", "random_gene") )

        ## read in the optimal file
        df <- read.delim('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_lambda_0.00025_gene_set_BA9_SALS_narval.csv', sep = ",")

        ## Only keep optimal
        BA9_SALS_weighted_df_lim <- BA9_SALS_weighted_df_lim %>%
            left_join(df, by = "celltype")

        ## Filter the rows where n_genes equals optimal_genes
        BA9_SALS_weighted_df_lim <- BA9_SALS_weighted_df_lim %>%
            filter(n_genes == optimal_genes)
        
        ## place in long format
        df_long <- BA9_SALS_weighted_df_lim %>%
            pivot_longer(cols = c(fold_1, fold_2, fold_3, fold_4, fold_5),
                        names_to = "fold",  # New column to hold the fold names
                        values_to = "accuracy") %>% as.data.frame()  # New column to hold the corresponding values
        
        ## Calculate mean and SD
        df_stats <- df_long %>%
            group_by(celltype, method) %>%  # Group by celltype and method
            summarize(
                mean_accuracy = mean(accuracy, na.rm = TRUE),  # Calculate mean accuracy
                sd_accuracy = sd(accuracy, na.rm = TRUE)      # Calculate standard deviation of accuracy
            ) %>% as.data.frame()
        
        df_stats <- df_stats %>%
            left_join(df, by = "celltype")

        ## set factor level
        df_stats$celltype <- factor(df_stats$celltype, levels = rev(c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  )))
        df_stats$optimal_genes[df_stats$method == "random_gene"] <- NA


        ## plot
        ggplot(df_stats, aes(x = mean_accuracy, y = celltype, fill = method, label = optimal_genes)) + 
        theme_bw() +
        geom_text(aes(x = 1.15), size = 3)+
        geom_vline(xintercept = 1.05) +
        geom_vline(xintercept = 0.90, colour = "lightgrey", linetype = 'dotted') +
        geom_bar(stat = 'identity', position = "dodge") +
        geom_errorbar(aes(xmin=mean_accuracy-sd_accuracy, xmax=mean_accuracy+sd_accuracy), width=.2,
                position=position_dodge(.9)) +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(colour = "black"),
            axis.text.y = element_text(colour = "black"),
            axis.title.y = element_blank(),
            axis.title.x = element_text(face = "bold")
        ) +
        xlab("Accuracy") +
        scale_fill_manual(values = c("#F18F01", "#2E4057")) +
        coord_cartesian(xlim=c(0.5,1.2)) +
        scale_x_continuous(
            breaks = c(0.5, 0.75, 1.00),          # Set the tick positions
            labels = c("0.5", "0.75", "1.0")      # Optional: Customize the tick labels
        )
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 3, width = 2.25)
    ##

##




############################################################################ Figure 2 -- Optimal LIME vs Random gene (Dont think we are going to use this. )
############################################################################
############################################################################
############################################################################ 


## code optimal LIME vs random gene accuracy
    ###########################
    ## SALS BA4
    ###########################
    ## code
        #celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")
        celltype_list = c("L3_L5",  "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV", "5HT3aR",  "Rosehip", "SOM",    "OPC",  "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA4"
        par_status = "SALS"
        par_prep = "CombatSeq"
        
        ################################################
        ## LIME outputs weighted by percent expression
        ## This function subsets the weighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        region = "fill"
        group = "fill"
        celltype = "fill"
        n_genes = 0
        test_accuracy = 0
        #mean_accuracy = 0
        fold_1 = 0
        fold_2 = 0
        fold_3 = 0
        fold_4 = 0
        fold_5 = 0
        method = "fill"
        fill <- data.frame(region,group, celltype, n_genes, fold_1, fold_2, fold_3, fold_4, fold_5, test_accuracy, method)
        

        for(celltype2 in unique(celltype_list)){
        
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/permutation_weighted_report_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df$five_fold_accuracies <- strsplit(as.character(df$five_fold_accuracies), "; ")
            df$five_fold_accuracies <- lapply(df$five_fold_accuracies, as.numeric)
            df <- cbind(df, do.call(rbind, df$five_fold_accuracies))
            colnames(df)[13:17] <- paste0("fold_", 1:5)
            #df$mean_accuracy <- rowMeans(df[, c("test_accuracy", "fold_1", "fold_2", "fold_3", "fold_4", "fold_5")], na.rm = TRUE)
            df <- df %>% dplyr::select(region, group, celltype, n_genes, fold_1, fold_2, fold_3, fold_4, fold_5, test_accuracy, method)
            fill <- rbind(fill, df)

        }

        fill <- subset(fill, region != "fill")
        BA4_SALS_weighted_df <- fill

        ## only keep LIME and random gene
        BA4_SALS_weighted_df_lim <- subset(BA4_SALS_weighted_df, method %in% c("LIME", "random_gene") )

        ## read in the optimal file
        df <- read.delim('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_lambda_0.00025_gene_set_BA4_SALS_narval.csv', sep = ",")

        ## Only keep optimal
        BA4_SALS_weighted_df_lim <- BA4_SALS_weighted_df_lim %>%
            left_join(df, by = "celltype")

        ## Filter the rows where n_genes equals optimal_genes
        BA4_SALS_weighted_df_lim <- BA4_SALS_weighted_df_lim %>%
            filter(n_genes == optimal_genes)

        ## place in long format
        df_long <- BA4_SALS_weighted_df_lim %>%
            pivot_longer(cols = c(fold_1, fold_2, fold_3, fold_4, fold_5, test_accuracy),
                        names_to = "fold",  # New column to hold the fold names
                        values_to = "accuracy") %>% as.data.frame()  # New column to hold the corresponding values
        
        ## Calculate mean and SD
        df_stats <- df_long %>%
            group_by(celltype, method) %>%  # Group by celltype and method
            summarize(
                mean_accuracy = mean(accuracy, na.rm = TRUE),  # Calculate mean accuracy
                sd_accuracy = sd(accuracy, na.rm = TRUE)      # Calculate standard deviation of accuracy
            ) %>% as.data.frame()
        
        df_stats <- df_stats %>%
            left_join(df, by = "celltype")

        ## set factor level
        df_stats$celltype <- factor(df_stats$celltype, levels = rev(c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  )))
        df_stats$optimal_genes[df_stats$method == "random_gene"] <- NA


        ## plot
        ggplot(df_stats, aes(x = mean_accuracy, y = celltype, fill = method, label = optimal_genes)) + 
        theme_bw() +
        geom_text(aes(x = 1.15), size = 3)+
        geom_vline(xintercept = 1.05) +
        geom_bar(stat = 'identity', position = "dodge") +
        geom_errorbar(aes(xmin=mean_accuracy-sd_accuracy, xmax=mean_accuracy+sd_accuracy), width=.2,
                position=position_dodge(.9)) +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(colour = "black"),
            axis.text.y = element_text(colour = "black"),
            axis.title.y = element_blank(),
            axis.title.x = element_text(face = "bold")
        ) +
        xlab("Accuracy") +
        scale_fill_manual(values = c("#F18F01", "#2E4057")) +
        coord_cartesian(xlim=c(0.5,1.2)) +
        scale_x_continuous(
            breaks = c(0.5, 0.75, 1.00),          # Set the tick positions
            labels = c("0.5", "0.75", "1.0")      # Optional: Customize the tick labels
        )
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 3, width = 2.25)
    ##

    ###########################
    ## C9ALS BA4
    ###########################
    ## code
        #celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",  "OPC",  "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA4"
        par_status = "C9ALS"
        par_prep = "CombatSeq"
        
        ################################################
        ## LIME outputs weighted by percent expression
        ## This function subsets the weighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        region = "fill"
        group = "fill"
        celltype = "fill"
        n_genes = 0
        test_accuracy = 0
        #mean_accuracy = 0
        fold_1 = 0
        fold_2 = 0
        fold_3 = 0
        fold_4 = 0
        fold_5 = 0
        method = "fill"
        fill <- data.frame(region,group, celltype, n_genes, fold_1, fold_2, fold_3, fold_4, fold_5, test_accuracy, method)
        

        for(celltype2 in unique(celltype_list)){
        
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/permutation_weighted_report_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df$five_fold_accuracies <- strsplit(as.character(df$five_fold_accuracies), "; ")
            df$five_fold_accuracies <- lapply(df$five_fold_accuracies, as.numeric)
            df <- cbind(df, do.call(rbind, df$five_fold_accuracies))
            colnames(df)[13:17] <- paste0("fold_", 1:5)
            #df$mean_accuracy <- rowMeans(df[, c("test_accuracy", "fold_1", "fold_2", "fold_3", "fold_4", "fold_5")], na.rm = TRUE)
            df <- df %>% dplyr::select(region, group, celltype, n_genes, fold_1, fold_2, fold_3, fold_4, fold_5, test_accuracy, method)
            fill <- rbind(fill, df)

        }

        fill <- subset(fill, region != "fill")
        BA4_C9ALS_weighted_df <- fill

        ## only keep LIME and random gene
        BA4_C9ALS_weighted_df_lim <- subset(BA4_C9ALS_weighted_df, method %in% c("LIME", "random_gene") )

        ## read in the optimal file
        df <- read.delim('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_lambda_0.00025_gene_set_BA4_C9ALS_narval.csv', sep = ",")

        ## Only keep optimal
        BA4_C9ALS_weighted_df_lim <- BA4_C9ALS_weighted_df_lim %>%
            left_join(df, by = "celltype")

        ## Filter the rows where n_genes equals optimal_genes
        BA4_C9ALS_weighted_df_lim <- BA4_C9ALS_weighted_df_lim %>%
            filter(n_genes == optimal_genes)

        ## place in long format
        df_long <- BA4_C9ALS_weighted_df_lim %>%
            pivot_longer(cols = c(fold_1, fold_2, fold_3, fold_4, fold_5, test_accuracy),
                        names_to = "fold",  # New column to hold the fold names
                        values_to = "accuracy") %>% as.data.frame()  # New column to hold the corresponding values
        
        ## Calculate mean and SD
        df_stats <- df_long %>%
            group_by(celltype, method) %>%  # Group by celltype and method
            summarize(
                mean_accuracy = mean(accuracy, na.rm = TRUE),  # Calculate mean accuracy
                sd_accuracy = sd(accuracy, na.rm = TRUE)      # Calculate standard deviation of accuracy
            ) %>% as.data.frame()
        
        df_stats <- df_stats %>%
            left_join(df, by = "celltype")

        ## set factor level
        df_stats$celltype <- factor(df_stats$celltype, levels = rev(c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  )))
        df_stats$optimal_genes[df_stats$method == "random_gene"] <- NA

        ## plot
        ggplot(df_stats, aes(x = mean_accuracy, y = celltype, fill = method, label = optimal_genes)) + 
        theme_bw() +
        geom_text(aes(x = 1.15), size = 3)+
        geom_vline(xintercept = 1.05) +
        geom_bar(stat = 'identity', position = "dodge") +
        geom_errorbar(aes(xmin=mean_accuracy-sd_accuracy, xmax=mean_accuracy+sd_accuracy), width=.2,
                position=position_dodge(.9)) +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(colour = "black"),
            axis.text.y = element_text(colour = "black"),
            axis.title.y = element_blank(),
            axis.title.x = element_text(face = "bold")
        ) +
        xlab("Accuracy") +
        scale_fill_manual(values = c("#F18F01", "#2E4057")) +
        coord_cartesian(xlim=c(0.5,1.2)) +
        scale_x_continuous(
            breaks = c(0.5, 0.75, 1.00),          # Set the tick positions
            labels = c("0.5", "0.75", "1.0")      # Optional: Customize the tick labels
        )
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 3, width = 2.25)
    ##

    ###########################
    ## C9FTLD BA4
    ###########################
    ## code
        #celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM", "Astro",   "OPC",  "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA4"
        par_status = "C9FTLD"
        par_prep = "CombatSeq"
        
        ################################################
        ## LIME outputs weighted by percent expression
        ## This function subsets the weighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        region = "fill"
        group = "fill"
        celltype = "fill"
        n_genes = 0
        test_accuracy = 0
        #mean_accuracy = 0
        fold_1 = 0
        fold_2 = 0
        fold_3 = 0
        fold_4 = 0
        fold_5 = 0
        method = "fill"
        fill <- data.frame(region,group, celltype, n_genes, fold_1, fold_2, fold_3, fold_4, fold_5, test_accuracy, method)
        

        for(celltype2 in unique(celltype_list)){
        
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/permutation_weighted_report_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df$five_fold_accuracies <- strsplit(as.character(df$five_fold_accuracies), "; ")
            df$five_fold_accuracies <- lapply(df$five_fold_accuracies, as.numeric)
            df <- cbind(df, do.call(rbind, df$five_fold_accuracies))
            colnames(df)[13:17] <- paste0("fold_", 1:5)
            #df$mean_accuracy <- rowMeans(df[, c("test_accuracy", "fold_1", "fold_2", "fold_3", "fold_4", "fold_5")], na.rm = TRUE)
            df <- df %>% dplyr::select(region, group, celltype, n_genes, fold_1, fold_2, fold_3, fold_4, fold_5, test_accuracy, method)
            fill <- rbind(fill, df)

        }

        fill <- subset(fill, region != "fill")
        BA4_C9FTLD_weighted_df <- fill

        ## only keep LIME and random gene
        BA4_C9FTLD_weighted_df_lim <- subset(BA4_C9FTLD_weighted_df, method %in% c("LIME", "random_gene") )

        ## read in the optimal file
        df <- read.delim('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_lambda_0.00025_gene_set_BA4_C9FTLD_narval.csv', sep = ",")

        ## Only keep optimal
        BA4_C9FTLD_weighted_df_lim <- BA4_C9FTLD_weighted_df_lim %>%
            left_join(df, by = "celltype")

        ## Filter the rows where n_genes equals optimal_genes
        BA4_C9FTLD_weighted_df_lim <- BA4_C9FTLD_weighted_df_lim %>%
            filter(n_genes == optimal_genes)

        ## place in long format
        df_long <- BA4_C9FTLD_weighted_df_lim %>%
            pivot_longer(cols = c(fold_1, fold_2, fold_3, fold_4, fold_5, test_accuracy),
                        names_to = "fold",  # New column to hold the fold names
                        values_to = "accuracy") %>% as.data.frame()  # New column to hold the corresponding values
        
        ## Calculate mean and SD
        df_stats <- df_long %>%
            group_by(celltype, method) %>%  # Group by celltype and method
            summarize(
                mean_accuracy = mean(accuracy, na.rm = TRUE),  # Calculate mean accuracy
                sd_accuracy = sd(accuracy, na.rm = TRUE)      # Calculate standard deviation of accuracy
            ) %>% as.data.frame()
        
        df_stats <- df_stats %>%
            left_join(df, by = "celltype")

        ## set factor level
        df_stats$celltype <- factor(df_stats$celltype, levels = rev(c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  )))
        df_stats$optimal_genes[df_stats$method == "random_gene"] <- NA

        ## plot
        ggplot(df_stats, aes(x = mean_accuracy, y = celltype, fill = method, label = optimal_genes)) + 
        theme_bw() +
        geom_text(aes(x = 1.15), size = 3)+
        geom_vline(xintercept = 1.05) +
        geom_bar(stat = 'identity', position = "dodge") +
        geom_errorbar(aes(xmin=mean_accuracy-sd_accuracy, xmax=mean_accuracy+sd_accuracy), width=.2,
                position=position_dodge(.9)) +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(colour = "black"),
            axis.text.y = element_text(colour = "black"),
            axis.title.y = element_blank(),
            axis.title.x = element_text(face = "bold")
        ) +
        xlab("Accuracy") +
        scale_fill_manual(values = c("#F18F01", "#2E4057")) +
        coord_cartesian(xlim=c(0.5,1.2)) +
        scale_x_continuous(
            breaks = c(0.5, 0.75, 1.00),          # Set the tick positions
            labels = c("0.5", "0.75", "1.0")      # Optional: Customize the tick labels
        )
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 3, width = 2.25)
    ##

##


## code optimal LIME vs random gene accuracy
    ###########################
    ## SALS BA4
    ###########################
    ## code
        #celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")
        celltype_list = c("L3_L5",  "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV", "5HT3aR",  "Rosehip", "SOM",    "OPC",  "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA4"
        par_status = "SALS"
        par_prep = "CombatSeq"
        
        ################################################
        ## LIME outputs weighted by percent expression
        ## This function subsets the weighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        region = "fill"
        group = "fill"
        celltype = "fill"
        n_genes = 0
        test_accuracy = 0
        #mean_accuracy = 0
        fold_1 = 0
        fold_2 = 0
        fold_3 = 0
        fold_4 = 0
        fold_5 = 0
        method = "fill"
        fill <- data.frame(region,group, celltype, n_genes, fold_1, fold_2, fold_3, fold_4, fold_5, test_accuracy, method)
        

        for(celltype2 in unique(celltype_list)){
        
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/permutation_weighted_report_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df$five_fold_accuracies <- strsplit(as.character(df$five_fold_accuracies), "; ")
            df$five_fold_accuracies <- lapply(df$five_fold_accuracies, as.numeric)
            df <- cbind(df, do.call(rbind, df$five_fold_accuracies))
            colnames(df)[13:17] <- paste0("fold_", 1:5)
            #df$mean_accuracy <- rowMeans(df[, c("test_accuracy", "fold_1", "fold_2", "fold_3", "fold_4", "fold_5")], na.rm = TRUE)
            df <- df %>% dplyr::select(region, group, celltype, n_genes, fold_1, fold_2, fold_3, fold_4, fold_5, test_accuracy, method)
            fill <- rbind(fill, df)

        }

        fill <- subset(fill, region != "fill")
        BA4_SALS_weighted_df <- fill

        ## only keep LIME and random gene
        BA4_SALS_weighted_df_lim <- subset(BA4_SALS_weighted_df, method %in% c("LIME", "random_gene") )

        ## read in the optimal file
        df <- read.delim('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_lambda_0.00025_gene_set_BA4_SALS_narval.csv', sep = ",")

        ## Only keep optimal
        BA4_SALS_weighted_df_lim <- BA4_SALS_weighted_df_lim %>%
            left_join(df, by = "celltype")

        ## Filter the rows where n_genes equals optimal_genes
        BA4_SALS_weighted_df_lim <- BA4_SALS_weighted_df_lim %>%
            filter(n_genes == optimal_genes)

        ## place in long format
        df_long <- BA4_SALS_weighted_df_lim %>%
            pivot_longer(cols = c(fold_1, fold_2, fold_3, fold_4, fold_5, test_accuracy),
                        names_to = "fold",  # New column to hold the fold names
                        values_to = "accuracy") %>% as.data.frame()  # New column to hold the corresponding values
        
        ## Calculate mean and SD
        df_stats <- df_long %>%
            group_by(celltype, method) %>%  # Group by celltype and method
            summarize(
                mean_accuracy = mean(accuracy, na.rm = TRUE),  # Calculate mean accuracy
                sd_accuracy = sd(accuracy, na.rm = TRUE)      # Calculate standard deviation of accuracy
            ) %>% as.data.frame()
        
        df_stats <- df_stats %>%
            left_join(df, by = "celltype")

        ## set factor level
        df_stats$celltype <- factor(df_stats$celltype, levels = rev(c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  )))
        df_stats$optimal_genes[df_stats$method == "random_gene"] <- NA


        ## plot
        ggplot(df_stats, aes(x = mean_accuracy, y = celltype, fill = method, label = optimal_genes)) + 
        theme_bw() +
        geom_text(aes(x = 1.15), size = 3)+
        geom_vline(xintercept = 1.05) +
        geom_bar(stat = 'identity', position = "dodge") +
        geom_errorbar(aes(xmin=mean_accuracy-sd_accuracy, xmax=mean_accuracy+sd_accuracy), width=.2,
                position=position_dodge(.9)) +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(colour = "black"),
            axis.text.y = element_text(colour = "black"),
            axis.title.y = element_blank(),
            axis.title.x = element_text(face = "bold")
        ) +
        xlab("Accuracy") +
        scale_fill_manual(values = c("#F18F01", "#2E4057")) +
        coord_cartesian(xlim=c(0.5,1.2)) +
        scale_x_continuous(
            breaks = c(0.5, 0.75, 1.00),          # Set the tick positions
            labels = c("0.5", "0.75", "1.0")      # Optional: Customize the tick labels
        )
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 3, width = 2.25)
    ##

    ###########################
    ## C9ALS BA4
    ###########################
    ## code
        #celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",  "OPC",  "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA4"
        par_status = "C9ALS"
        par_prep = "CombatSeq"
        
        ################################################
        ## LIME outputs weighted by percent expression
        ## This function subsets the weighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        region = "fill"
        group = "fill"
        celltype = "fill"
        n_genes = 0
        test_accuracy = 0
        #mean_accuracy = 0
        fold_1 = 0
        fold_2 = 0
        fold_3 = 0
        fold_4 = 0
        fold_5 = 0
        method = "fill"
        fill <- data.frame(region,group, celltype, n_genes, fold_1, fold_2, fold_3, fold_4, fold_5, test_accuracy, method)
        

        for(celltype2 in unique(celltype_list)){
        
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/permutation_weighted_report_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df$five_fold_accuracies <- strsplit(as.character(df$five_fold_accuracies), "; ")
            df$five_fold_accuracies <- lapply(df$five_fold_accuracies, as.numeric)
            df <- cbind(df, do.call(rbind, df$five_fold_accuracies))
            colnames(df)[13:17] <- paste0("fold_", 1:5)
            #df$mean_accuracy <- rowMeans(df[, c("test_accuracy", "fold_1", "fold_2", "fold_3", "fold_4", "fold_5")], na.rm = TRUE)
            df <- df %>% dplyr::select(region, group, celltype, n_genes, fold_1, fold_2, fold_3, fold_4, fold_5, test_accuracy, method)
            fill <- rbind(fill, df)

        }

        fill <- subset(fill, region != "fill")
        BA4_C9ALS_weighted_df <- fill

        ## only keep LIME and random gene
        BA4_C9ALS_weighted_df_lim <- subset(BA4_C9ALS_weighted_df, method %in% c("LIME", "random_gene") )

        ## read in the optimal file
        df <- read.delim('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_lambda_0.00025_gene_set_BA4_C9ALS_narval.csv', sep = ",")

        ## Only keep optimal
        BA4_C9ALS_weighted_df_lim <- BA4_C9ALS_weighted_df_lim %>%
            left_join(df, by = "celltype")

        ## Filter the rows where n_genes equals optimal_genes
        BA4_C9ALS_weighted_df_lim <- BA4_C9ALS_weighted_df_lim %>%
            filter(n_genes == optimal_genes)

        ## place in long format
        df_long <- BA4_C9ALS_weighted_df_lim %>%
            pivot_longer(cols = c(fold_1, fold_2, fold_3, fold_4, fold_5, test_accuracy),
                        names_to = "fold",  # New column to hold the fold names
                        values_to = "accuracy") %>% as.data.frame()  # New column to hold the corresponding values
        
        ## Calculate mean and SD
        df_stats <- df_long %>%
            group_by(celltype, method) %>%  # Group by celltype and method
            summarize(
                mean_accuracy = mean(accuracy, na.rm = TRUE),  # Calculate mean accuracy
                sd_accuracy = sd(accuracy, na.rm = TRUE)      # Calculate standard deviation of accuracy
            ) %>% as.data.frame()
        
        df_stats <- df_stats %>%
            left_join(df, by = "celltype")

        ## set factor level
        df_stats$celltype <- factor(df_stats$celltype, levels = rev(c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  )))
        df_stats$optimal_genes[df_stats$method == "random_gene"] <- NA

        ## plot
        ggplot(df_stats, aes(x = mean_accuracy, y = celltype, fill = method, label = optimal_genes)) + 
        theme_bw() +
        geom_text(aes(x = 1.15), size = 3)+
        geom_vline(xintercept = 1.05) +
        geom_bar(stat = 'identity', position = "dodge") +
        geom_errorbar(aes(xmin=mean_accuracy-sd_accuracy, xmax=mean_accuracy+sd_accuracy), width=.2,
                position=position_dodge(.9)) +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(colour = "black"),
            axis.text.y = element_text(colour = "black"),
            axis.title.y = element_blank(),
            axis.title.x = element_text(face = "bold")
        ) +
        xlab("Accuracy") +
        scale_fill_manual(values = c("#F18F01", "#2E4057")) +
        coord_cartesian(xlim=c(0.5,1.2)) +
        scale_x_continuous(
            breaks = c(0.5, 0.75, 1.00),          # Set the tick positions
            labels = c("0.5", "0.75", "1.0")      # Optional: Customize the tick labels
        )
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 3, width = 2.25)
    ##

    ###########################
    ## C9FTLD BA4
    ###########################
    ## code
        #celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM", "Astro",   "OPC",  "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA4"
        par_status = "C9FTLD"
        par_prep = "CombatSeq"
        
        ################################################
        ## LIME outputs weighted by percent expression
        ## This function subsets the weighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        region = "fill"
        group = "fill"
        celltype = "fill"
        n_genes = 0
        test_accuracy = 0
        #mean_accuracy = 0
        fold_1 = 0
        fold_2 = 0
        fold_3 = 0
        fold_4 = 0
        fold_5 = 0
        method = "fill"
        fill <- data.frame(region,group, celltype, n_genes, fold_1, fold_2, fold_3, fold_4, fold_5, test_accuracy, method)
        

        for(celltype2 in unique(celltype_list)){
        
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/permutation_weighted_report_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df$five_fold_accuracies <- strsplit(as.character(df$five_fold_accuracies), "; ")
            df$five_fold_accuracies <- lapply(df$five_fold_accuracies, as.numeric)
            df <- cbind(df, do.call(rbind, df$five_fold_accuracies))
            colnames(df)[13:17] <- paste0("fold_", 1:5)
            #df$mean_accuracy <- rowMeans(df[, c("test_accuracy", "fold_1", "fold_2", "fold_3", "fold_4", "fold_5")], na.rm = TRUE)
            df <- df %>% dplyr::select(region, group, celltype, n_genes, fold_1, fold_2, fold_3, fold_4, fold_5, test_accuracy, method)
            fill <- rbind(fill, df)

        }

        fill <- subset(fill, region != "fill")
        BA4_C9FTLD_weighted_df <- fill

        ## only keep LIME and random gene
        BA4_C9FTLD_weighted_df_lim <- subset(BA4_C9FTLD_weighted_df, method %in% c("LIME", "random_gene") )

        ## read in the optimal file
        df <- read.delim('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_lambda_0.00025_gene_set_BA4_C9FTLD_narval.csv', sep = ",")

        ## Only keep optimal
        BA4_C9FTLD_weighted_df_lim <- BA4_C9FTLD_weighted_df_lim %>%
            left_join(df, by = "celltype")

        ## Filter the rows where n_genes equals optimal_genes
        BA4_C9FTLD_weighted_df_lim <- BA4_C9FTLD_weighted_df_lim %>%
            filter(n_genes == optimal_genes)

        ## place in long format
        df_long <- BA4_C9FTLD_weighted_df_lim %>%
            pivot_longer(cols = c(fold_1, fold_2, fold_3, fold_4, fold_5, test_accuracy),
                        names_to = "fold",  # New column to hold the fold names
                        values_to = "accuracy") %>% as.data.frame()  # New column to hold the corresponding values
        
        ## Calculate mean and SD
        df_stats <- df_long %>%
            group_by(celltype, method) %>%  # Group by celltype and method
            summarize(
                mean_accuracy = mean(accuracy, na.rm = TRUE),  # Calculate mean accuracy
                sd_accuracy = sd(accuracy, na.rm = TRUE)      # Calculate standard deviation of accuracy
            ) %>% as.data.frame()
        
        df_stats <- df_stats %>%
            left_join(df, by = "celltype")

        ## set factor level
        df_stats$celltype <- factor(df_stats$celltype, levels = rev(c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  )))
        df_stats$optimal_genes[df_stats$method == "random_gene"] <- NA

        ## plot
        ggplot(df_stats, aes(x = mean_accuracy, y = celltype, fill = method, label = optimal_genes)) + 
        theme_bw() +
        geom_text(aes(x = 1.15), size = 3)+
        geom_vline(xintercept = 1.05) +
        geom_bar(stat = 'identity', position = "dodge") +
        geom_errorbar(aes(xmin=mean_accuracy-sd_accuracy, xmax=mean_accuracy+sd_accuracy), width=.2,
                position=position_dodge(.9)) +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(colour = "black"),
            axis.text.y = element_text(colour = "black"),
            axis.title.y = element_blank(),
            axis.title.x = element_text(face = "bold")
        ) +
        xlab("Accuracy") +
        scale_fill_manual(values = c("#F18F01", "#2E4057")) +
        coord_cartesian(xlim=c(0.5,1.2)) +
        scale_x_continuous(
            breaks = c(0.5, 0.75, 1.00),          # Set the tick positions
            labels = c("0.5", "0.75", "1.0")      # Optional: Customize the tick labels
        )
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 3, width = 2.25)
    ##

##




















## old code

############################################################################ Z-score distributions -- unweighted by percent expression (potential Supplemental Figure)
############################################################################
############################################################################
############################################################################ 

## code
    ###########################
    ## SALS BA4
    ###########################
    ## code
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA4"
        par_status = "SALS"
        par_prep = "CombatSeq"
        
        ################################################
        ## LIME outputs unweighted by percent expression
        ## This function subsets the unweighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_across_donors = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_across_donors, celltype)

        
        for(celltype2 in unique(celltype_list)){
            
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_unweighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_across_donors <- as.numeric(df$Mean_z_score_across_donors)
            df <- df[order( -(df$Mean_z_score_across_donors)),]
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_across_donors)
            df$celltype <- celltype2
            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")
        BA4_SALS_unweighted_df <- fill
    ##

    ###########################
    ## C9ALS BA4
    ###########################
    ## code
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA4"
        par_status = "C9ALS"
        par_prep = "CombatSeq"
        
        ################################################
        ## LIME outputs unweighted by percent expression
        ## This function subsets the unweighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_across_donors = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_across_donors, celltype)

        for(celltype2 in unique(celltype_list)){
            
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_unweighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_across_donors <- as.numeric(df$Mean_z_score_across_donors)
            df <- df[order( -(df$Mean_z_score_across_donors)),]
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_across_donors)
            df$celltype <- celltype2
            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")
        BA4_C9ALS_unweighted_df <- fill

    ##

    ###########################
    ## FTLD BA4
    ###########################
    ## code
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA4"
        par_status = "SFTLD"
        par_prep = "CombatSeq"
        
        ################################################
        ## LIME outputs unweighted by percent expression
        ## This function subsets the unweighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_across_donors = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_across_donors, celltype)

        for(celltype2 in unique(celltype_list)){
            
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_unweighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_across_donors <- as.numeric(df$Mean_z_score_across_donors)
            df <- df[order( -(df$Mean_z_score_across_donors)),]
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_across_donors)
            df$celltype <- celltype2
            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")
        BA4_SFTLD_unweighted_df <- fill

    ##
    
    ###########################
    ## C9FTLD BA4
    ###########################
    ## code
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA4"
        par_status = "C9FTLD"
        par_prep = "CombatSeq"

        ################################################
        ## LIME outputs unweighted by percent expression
        ## This function subsets the unweighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_across_donors = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_across_donors, celltype)

        for(celltype2 in unique(celltype_list)){
            
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_unweighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_across_donors <- as.numeric(df$Mean_z_score_across_donors)
            df <- df[order( -(df$Mean_z_score_across_donors)),]
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_across_donors)
            df$celltype <- celltype2
            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")
        BA4_C9FTLD_unweighted_df <- fill

    ##
    
    ###########################
    ## SALS BA9
    ###########################
    ## code 
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA9"
        par_status = "SALS"
        par_prep = "CombatSeq"
        
        ################################################
        ## LIME outputs unweighted by percent expression
        ## This function subsets the unweighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_across_donors = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_across_donors, celltype)

        for(celltype2 in unique(celltype_list)){
            
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_unweighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_across_donors <- as.numeric(df$Mean_z_score_across_donors)
            df <- df[order( -(df$Mean_z_score_across_donors)),]
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_across_donors)
            df$celltype <- celltype2
            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")
        BA9_SALS_unweighted_df <- fill

    ##
    
    ###########################
    ## C9ALS BA9
    ###########################
    ## code 
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA9"
        par_status = "C9ALS"
        par_prep = "CombatSeq"

        ################################################
        ## LIME outputs unweighted by percent expression
        ## This function subsets the unweighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_across_donors = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_across_donors, celltype)

        for(celltype2 in unique(celltype_list)){
            
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_unweighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_across_donors <- as.numeric(df$Mean_z_score_across_donors)
            df <- df[order( -(df$Mean_z_score_across_donors)),]
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_across_donors)
            df$celltype <- celltype2
            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")
        BA9_C9ALS_unweighted_df <- fill

    ##

    ###########################
    ## SFTLD BA9
    ###########################
    ## code 
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA9"
        par_status = "SFTLD"
        par_prep = "CombatSeq"
        
        ################################################
        ## LIME outputs unweighted by percent expression
        ## This function subsets the unweighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_across_donors = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_across_donors, celltype)

        for(celltype2 in unique(celltype_list)){
            
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_unweighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_across_donors <- as.numeric(df$Mean_z_score_across_donors)
            df <- df[order( -(df$Mean_z_score_across_donors)),]
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_across_donors)
            df$celltype <- celltype2
            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")
        BA9_SFTLD_unweighted_df <- fill

    ##

    ###########################
    ## C9FTLD BA9 -- missing TCELL for some reason -- look into it
    ###########################
    ## code 

        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA9"
        par_status = "C9FTLD"
        par_prep = "CombatSeq"
        
        ################################################
        ## LIME outputs unweighted by percent expression
        ## This function subsets the unweighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_across_donors = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_across_donors, celltype)

        for(celltype2 in unique(celltype_list)){
            
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_unweighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_across_donors <- as.numeric(df$Mean_z_score_across_donors)
            df <- df[order( -(df$Mean_z_score_across_donors)),]
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_across_donors)
            df$celltype <- celltype2
            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")
        BA9_C9FTLD_unweighted_df <- fill

    ##

    ###########################
    ## All conditions BA4
    ###########################
    ## code
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA4"
        par_status = "all_conditions"
        par_prep = "CombatSeq"

        ################################################
        ## LIME outputs unweighted by percent expression
        ## This function subsets the unweighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_across_donors = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_across_donors, celltype)

        for(celltype2 in unique(celltype_list)){
            
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_unweighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_across_donors <- as.numeric(df$Mean_z_score_across_donors)
            df <- df[order( -(df$Mean_z_score_across_donors)),]
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_across_donors)
            df$celltype <- celltype2
            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")
        BA4_all_conditions_unweighted_df <- fill

    ##

    ###########################
    ## All conditions BA9
    ###########################
    ## code
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA9"
        par_status = "all_conditions"
        par_prep = "CombatSeq"

        ################################################
        ## LIME outputs unweighted by percent expression
        ## This function subsets the unweighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_across_donors = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_across_donors, celltype)

        for(celltype2 in unique(celltype_list)){
            
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_unweighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_across_donors <- as.numeric(df$Mean_z_score_across_donors)
            df <- df[order( -(df$Mean_z_score_across_donors)),]
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_across_donors)
            df$celltype <- celltype2
            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")
        BA9_all_conditions_unweighted_df <- fill

    ##

    ###########################
    ## Plot BA4
    ###########################
    ## code
        BA4_SALS_unweighted_df$disease_cat <- "SALS"
        BA4_C9ALS_unweighted_df$disease_cat <- "C9ALS"
        BA4_SFTLD_unweighted_df$disease_cat <- "SFTLD"
        BA4_C9FTLD_unweighted_df$disease_cat <- "C9FTLD"
        BA4_all_conditions_unweighted_df$disease_cat <- "All"

        BA4_unweighted_distribution <- rbind(BA4_SALS_unweighted_df, BA4_C9ALS_unweighted_df, BA4_SFTLD_unweighted_df, BA4_C9FTLD_unweighted_df, BA4_all_conditions_unweighted_df)

        BA4_unweighted_distribution$disease_cat <- factor(BA4_unweighted_distribution$disease_cat, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All"))
        BA4_unweighted_distribution$celltype <- factor(BA4_unweighted_distribution$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))


        ggplot(BA4_unweighted_distribution, aes(x = Mean_z_score_across_donors, colour = disease_cat)) + 
        theme_bw() + 
        geom_vline(xintercept = 1.0, linetype = "dashed", colour = "black") +
        geom_density() +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(colour = "black", angle = 90, hjust = 1, vjust = 0.5),
            axis.text.y = element_text(colour = "black"),
            axis.title.y = element_text(face = "bold"),
            axis.title.x = element_text(face = "bold"),
        ) +
        facet_grid(disease_cat~celltype, scales = "free_x") +
        ylab("Density") + xlab("LIME feature importance Z-score") +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "black")) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 4, width = 15)
    ##

    ###########################
    ## Plot BA9
    ###########################
    ## code
        BA9_SALS_unweighted_df$disease_cat <- "SALS"
        BA9_C9ALS_unweighted_df$disease_cat <- "C9ALS"
        BA9_SFTLD_unweighted_df$disease_cat <- "SFTLD"
        BA9_C9FTLD_unweighted_df$disease_cat <- "C9FTLD"
        BA9_all_conditions_unweighted_df$disease_cat <- "All"

        BA9_unweighted_distribution <- rbind(BA9_SALS_unweighted_df, BA9_C9ALS_unweighted_df, BA9_SFTLD_unweighted_df, BA9_C9FTLD_unweighted_df, BA9_all_conditions_unweighted_df)

        BA9_unweighted_distribution$disease_cat <- factor(BA9_unweighted_distribution$disease_cat, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All"))
        BA9_unweighted_distribution$celltype <- factor(BA9_unweighted_distribution$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))


        ggplot(BA9_unweighted_distribution, aes(x = Mean_z_score_across_donors, colour = disease_cat)) + 
        theme_bw() + 
        geom_vline(xintercept = 1.0, linetype = "dashed", colour = "black") +
        geom_density() +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(colour = "black", angle = 90, hjust = 1, vjust = 0.5),
            axis.text.y = element_text(colour = "black"),
            axis.title.y = element_text(face = "bold"),
            axis.title.x = element_text(face = "bold"),
        ) +
        facet_grid(disease_cat~celltype, scales = "free_x") +
        ylab("Density") + xlab("LIME feature importance Z-score") +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "black")) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 4, width = 15)
    ##

##

############################################################################ N features with Z-score > 1 -- unweighted by percent expression (potential Supplemental Figure)
############################################################################
############################################################################
############################################################################ 

## code
    ###########################
    ## SALS BA4
    ###########################
    ## code
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA4"
        par_status = "SALS"
        par_prep = "CombatSeq"
        
        ################################################
        ## LIME outputs unweighted by percent expression
        ## This function subsets the unweighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_across_donor = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_across_donor, celltype)

        
        for(celltype2 in unique(celltype_list)){
            
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_unweighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_across_donor <- as.numeric(df$Mean_z_score_across_donor)
            df <- df[order( -(df$Mean_z_score_across_donor)),]
            df <- subset(df, Mean_z_score_across_donor >= 1)
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_across_donor)
            df$celltype <- celltype2
            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")
        BA4_SALS_unweighted_df <- fill
    ##

    ###########################
    ## C9ALS BA4
    ###########################
    ## code
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA4"
        par_status = "C9ALS"
        par_prep = "CombatSeq"
        
        ################################################
        ## LIME outputs unweighted by percent expression
        ## This function subsets the unweighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_across_donor = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_across_donor, celltype)

        for(celltype2 in unique(celltype_list)){
            
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_unweighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_across_donor <- as.numeric(df$Mean_z_score_across_donor)
            df <- df[order( -(df$Mean_z_score_across_donor)),]
            df <- subset(df, Mean_z_score_across_donor >= 1)
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_across_donor)
            df$celltype <- celltype2
            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")
        BA4_C9ALS_unweighted_df <- fill

    ##

    ###########################
    ## SFTLD BA4
    ###########################
    ## code
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA4"
        par_status = "SFTLD"
        par_prep = "CombatSeq"
        
        ################################################
        ## LIME outputs unweighted by percent expression
        ## This function subsets the unweighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_across_donor = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_across_donor, celltype)

        for(celltype2 in unique(celltype_list)){
            
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_unweighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_across_donor <- as.numeric(df$Mean_z_score_across_donor)
            df <- df[order( -(df$Mean_z_score_across_donor)),]
            df <- subset(df, Mean_z_score_across_donor >= 1)
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_across_donor)
            df$celltype <- celltype2
            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")
        BA4_SFTLD_unweighted_df <- fill

    ##
    
    ###########################
    ## C9FTLD BA4
    ###########################
    ## code
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA4"
        par_status = "C9FTLD"
        par_prep = "CombatSeq"

        ################################################
        ## LIME outputs unweighted by percent expression
        ## This function subsets the unweighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_across_donor = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_across_donor, celltype)

        for(celltype2 in unique(celltype_list)){
            
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_unweighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_across_donor <- as.numeric(df$Mean_z_score_across_donor)
            df <- df[order( -(df$Mean_z_score_across_donor)),]
            df <- subset(df, Mean_z_score_across_donor >= 1)
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_across_donor)
            df$celltype <- celltype2
            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")
        BA4_C9FTLD_unweighted_df <- fill

    ##
    
    ###########################
    ## SALS BA9
    ###########################
    ## code 
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA9"
        par_status = "SALS"
        par_prep = "CombatSeq"
        
        ################################################
        ## LIME outputs unweighted by percent expression
        ## This function subsets the unweighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_across_donor = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_across_donor, celltype)

        for(celltype2 in unique(celltype_list)){
            
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_unweighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_across_donor <- as.numeric(df$Mean_z_score_across_donor)
            df <- df[order( -(df$Mean_z_score_across_donor)),]
            df <- subset(df, Mean_z_score_across_donor >= 1)
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_across_donor)
            df$celltype <- celltype2
            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")
        BA9_SALS_unweighted_df <- fill

    ##
    
    ###########################
    ## C9ALS BA9
    ###########################
    ## code 
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA9"
        par_status = "C9ALS"
        par_prep = "CombatSeq"

        ################################################
        ## LIME outputs unweighted by percent expression
        ## This function subsets the unweighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_across_donor = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_across_donor, celltype)

        for(celltype2 in unique(celltype_list)){
            
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_unweighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_across_donor <- as.numeric(df$Mean_z_score_across_donor)
            df <- df[order( -(df$Mean_z_score_across_donor)),]
            df <- subset(df, Mean_z_score_across_donor >= 1)
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_across_donor)
            df$celltype <- celltype2
            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")
        BA9_C9ALS_unweighted_df <- fill

    ##

    ###########################
    ## SFTLD BA9
    ###########################
    ## code 
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA9"
        par_status = "SFTLD"
        par_prep = "CombatSeq"
        
        ################################################
        ## LIME outputs unweighted by percent expression
        ## This function subsets the unweighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_across_donor = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_across_donor, celltype)

        for(celltype2 in unique(celltype_list)){
            
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_unweighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_across_donor <- as.numeric(df$Mean_z_score_across_donor)
            df <- df[order( -(df$Mean_z_score_across_donor)),]
            df <- subset(df, Mean_z_score_across_donor >= 1)
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_across_donor)
            df$celltype <- celltype2
            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")
        BA9_SFTLD_unweighted_df <- fill

    ##

    ###########################
    ## C9FTLD BA9 -- missing TCELL for some reason -- look into it
    ###########################
    ## code 

        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA9"
        par_status = "C9FTLD"
        par_prep = "CombatSeq"
        
        ################################################
        ## LIME outputs unweighted by percent expression
        ## This function subsets the unweighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_across_donor = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_across_donor, celltype)

        for(celltype2 in unique(celltype_list)){
            
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_unweighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_across_donor <- as.numeric(df$Mean_z_score_across_donor)
            df <- df[order( -(df$Mean_z_score_across_donor)),]
            df <- subset(df, Mean_z_score_across_donor >= 1)
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_across_donor)
            df$celltype <- celltype2
            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")
        BA9_C9FTLD_unweighted_df <- fill

    ##

    ###########################
    ## All conditions BA4
    ###########################
    ## code
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA4"
        par_status = "all_conditions"
        par_prep = "CombatSeq"

        ################################################
        ## LIME outputs unweighted by percent expression
        ## This function subsets the unweighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_across_donor = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_across_donor, celltype)

        for(celltype2 in unique(celltype_list)){
            
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_unweighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_across_donor <- as.numeric(df$Mean_z_score_across_donor)
            df <- df[order( -(df$Mean_z_score_across_donor)),]
            df <- subset(df, Mean_z_score_across_donor >= 1)
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_across_donor)
            df$celltype <- celltype2
            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")
        BA4_all_conditions_unweighted_df <- fill

    ##

    ###########################
    ## All conditions BA9
    ###########################
    ## code
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA9"
        par_status = "all_conditions"
        par_prep = "CombatSeq"

        ################################################
        ## LIME outputs unweighted by percent expression
        ## This function subsets the unweighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_across_donor = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_across_donor, celltype)

        for(celltype2 in unique(celltype_list)){
            
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_unweighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_across_donor <- as.numeric(df$Mean_z_score_across_donor)
            df <- df[order( -(df$Mean_z_score_across_donor)),]
            df <- subset(df, Mean_z_score_across_donor >= 1)
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_across_donor)
            df$celltype <- celltype2
            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")
        BA9_all_conditions_unweighted_df <- fill

    ##

    ###########################
    ## Plot BA4
    ###########################
    ## code
        BA4_SALS_unweighted_df$disease_cat <- "SALS"
        BA4_C9ALS_unweighted_df$disease_cat <- "C9ALS"
        BA4_SFTLD_unweighted_df$disease_cat <- "SFTLD"
        BA4_C9FTLD_unweighted_df$disease_cat <- "C9FTLD"
        BA4_all_conditions_unweighted_df$disease_cat <- "All"

        BA4_unweighted_counts <- rbind(BA4_SALS_unweighted_df, BA4_C9ALS_unweighted_df, BA4_SFTLD_unweighted_df, BA4_C9FTLD_unweighted_df, BA4_all_conditions_unweighted_df)

        counts_df <- data.frame(table(BA4_unweighted_counts$celltype, BA4_unweighted_counts$disease_cat ))
        
        counts_df$Var2 <- factor(counts_df$Var2, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All"))
        counts_df$Var1 <- factor(counts_df$Var1, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))

        ggplot(counts_df, aes(x = Var2, y = Freq, fill = Var2)) + 
        theme_bw() + 
        geom_bar(stat = 'identity') +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(colour = "black", angle = 90, hjust = 1, vjust = 0.5),
            axis.text.y = element_text(colour = "black"),
            axis.title.y = element_text(face = "bold"),
            axis.title.x = element_blank(),
        ) +
        facet_grid(~Var1, scales = "free_x") +
        ylab("N genes\n(Z-score > 1)") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 2, width = 15)
    ##

    ###########################
    ## Plot BA9
    ###########################
    ## code
        BA9_SALS_unweighted_df$disease_cat <- "SALS"
        BA9_C9ALS_unweighted_df$disease_cat <- "C9ALS"
        BA9_SFTLD_unweighted_df$disease_cat <- "SFTLD"
        BA9_C9FTLD_unweighted_df$disease_cat <- "C9FTLD"
        BA9_all_conditions_unweighted_df$disease_cat <- "All"

        BA9_unweighted_counts <- rbind(BA9_SALS_unweighted_df, BA9_C9ALS_unweighted_df, BA9_SFTLD_unweighted_df, BA9_C9FTLD_unweighted_df, BA9_all_conditions_unweighted_df)

        counts_df <- data.frame(table(BA9_unweighted_counts$celltype, BA9_unweighted_counts$disease_cat ))
        
        counts_df$Var2 <- factor(counts_df$Var2, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All"))
        counts_df$Var1 <- factor(counts_df$Var1, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))

        ggplot(counts_df, aes(x = Var2, y = Freq, fill = Var2)) + 
        theme_bw() + 
        geom_bar(stat = 'identity') +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(colour = "black", angle = 90, hjust = 1, vjust = 0.5),
            axis.text.y = element_text(colour = "black"),
            axis.title.y = element_text(face = "bold"),
            axis.title.x = element_blank(),
        ) +
        facet_grid(~Var1, scales = "free_x") +
        ylab("N genes\n(Z-score > 1)") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 2, width = 15)
    ##

##

############################################################################ Unweighted permutation optimal threshold test
############################################################################
############################################################################
############################################################################ 

## code
    ###########################
    ## SALS BA4
    ###########################
    ## code
        #celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")
        celltype_list = c("L3_L5", "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",  "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA4"
        par_status = "SALS"
        par_prep = "CombatSeq"
        
        ################################################
        ## LIME outputs unweighted by percent expression
        ## This function subsets the unweighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        region = "fill"
        group = "fill"
        celltype = "fill"
        n_genes = 0
        test_accuracy = 0
        mean_accuracy = 0
        method = "fill"
        fill <- data.frame(region,group, celltype, n_genes, test_accuracy, mean_accuracy, method)
        
        for(celltype2 in unique(celltype_list)){
        
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/permutation_unweighted_report_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df$five_fold_accuracies <- strsplit(as.character(df$five_fold_accuracies), "; ")
            df$five_fold_accuracies <- lapply(df$five_fold_accuracies, as.numeric)
            df <- cbind(df, do.call(rbind, df$five_fold_accuracies))
            colnames(df)[13:17] <- paste0("fold_", 1:5)
            df$mean_accuracy <- rowMeans(df[, c("test_accuracy", "fold_1", "fold_2", "fold_3", "fold_4", "fold_5")], na.rm = TRUE)
            df <- df %>% dplyr::select(region, group, celltype, n_genes, test_accuracy, mean_accuracy, method)
            fill <- rbind(fill, df)

        }

        fill <- subset(fill, region != "fill")
        BA4_SALS_unweighted_df <- fill
    ##

    ###########################
    ## C9ALS BA4
    ###########################
    ## code
        #celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")
        celltype_list = c("L3_L5",    "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",   "OPC",  "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA4"
        par_status = "C9ALS"
        par_prep = "CombatSeq"
        
        ################################################
        ## LIME outputs unweighted by percent expression
        ## This function subsets the unweighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        region = "fill"
        group = "fill"
        celltype = "fill"
        n_genes = 0
        test_accuracy = 0
        mean_accuracy = 0
        method = "fill"
        fill <- data.frame(region,group, celltype, n_genes, test_accuracy, mean_accuracy, method)
        
        for(celltype2 in unique(celltype_list)){
        
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/permutation_unweighted_report_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df$five_fold_accuracies <- strsplit(as.character(df$five_fold_accuracies), "; ")
            df$five_fold_accuracies <- lapply(df$five_fold_accuracies, as.numeric)
            df <- cbind(df, do.call(rbind, df$five_fold_accuracies))
            colnames(df)[13:17] <- paste0("fold_", 1:5)
            df$mean_accuracy <- rowMeans(df[, c("test_accuracy", "fold_1", "fold_2", "fold_3", "fold_4", "fold_5")], na.rm = TRUE)
            df <- df %>% dplyr::select(region, group, celltype, n_genes, test_accuracy, mean_accuracy, method)
            fill <- rbind(fill, df)

        }

        fill <- subset(fill, region != "fill")
        BA4_C9ALS_unweighted_df <- fill
    ##

    ###########################
    ## SFTLD BA4
    ###########################
    ## code
        #celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")
        celltype_list = c("L3_L5", "L2_L3",  "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",  "Astro",   "OPC",  "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA4"
        par_status = "SFTLD"
        par_prep = "CombatSeq"
        
        ################################################
        ## LIME outputs unweighted by percent expression
        ## This function subsets the unweighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        region = "fill"
        group = "fill"
        celltype = "fill"
        n_genes = 0
        test_accuracy = 0
        mean_accuracy = 0
        method = "fill"
        fill <- data.frame(region,group, celltype, n_genes, test_accuracy, mean_accuracy, method)
        
        for(celltype2 in unique(celltype_list)){
        
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/permutation_unweighted_report_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df$five_fold_accuracies <- strsplit(as.character(df$five_fold_accuracies), "; ")
            df$five_fold_accuracies <- lapply(df$five_fold_accuracies, as.numeric)
            df <- cbind(df, do.call(rbind, df$five_fold_accuracies))
            colnames(df)[13:17] <- paste0("fold_", 1:5)
            df$mean_accuracy <- rowMeans(df[, c("test_accuracy", "fold_1", "fold_2", "fold_3", "fold_4", "fold_5")], na.rm = TRUE)
            df <- df %>% dplyr::select(region, group, celltype, n_genes, test_accuracy, mean_accuracy, method)
            fill <- rbind(fill, df)

        }

        fill <- subset(fill, region != "fill")
        BA4_SFTLD_unweighted_df <- fill
    ##

    ###########################
    ## C9FTLD BA4 -- did not run yet
    ###########################
    ## code
        #celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")
        #celltype_list = c("L3_L5", "L2_L3",  "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",  "Astro",   "OPC",  "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA4"
        par_status = "C9FTLD"
        par_prep = "CombatSeq"
        
        ################################################
        ## LIME outputs unweighted by percent expression
        ## This function subsets the unweighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        region = "fill"
        group = "fill"
        celltype = "fill"
        n_genes = 0
        test_accuracy = 0
        mean_accuracy = 0
        method = "fill"
        fill <- data.frame(region,group, celltype, n_genes, test_accuracy, mean_accuracy, method)
        
        for(celltype2 in unique(celltype_list)){
        
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/permutation_unweighted_report_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df$five_fold_accuracies <- strsplit(as.character(df$five_fold_accuracies), "; ")
            df$five_fold_accuracies <- lapply(df$five_fold_accuracies, as.numeric)
            df <- cbind(df, do.call(rbind, df$five_fold_accuracies))
            colnames(df)[13:17] <- paste0("fold_", 1:5)
            df$mean_accuracy <- rowMeans(df[, c("test_accuracy", "fold_1", "fold_2", "fold_3", "fold_4", "fold_5")], na.rm = TRUE)
            df <- df %>% dplyr::select(region, group, celltype, n_genes, test_accuracy, mean_accuracy, method)
            fill <- rbind(fill, df)

        }

        fill <- subset(fill, region != "fill")
        BA4_C9FTLD_unweighted_df <- fill
    ##

    ###########################
    ## Plot BA4
    ###########################
    ## code
        ## sumarize SALS
        BA4_SALS_unweighted_df$disease_cat <- "SALS"
        BA4_SALS_unweighted_df_summary <- BA4_SALS_unweighted_df %>%
            dplyr::group_by(group, celltype, n_genes, method) %>%
            dplyr::summarise(mean_accuracy_by_group = mean(mean_accuracy, na.rm = TRUE))
        BA4_SALS_unweighted_df_summary <- data.frame(BA4_SALS_unweighted_df_summary)

        ## sumarize C9ALS
        BA4_C9ALS_unweighted_df$disease_cat <- "C9ALS"
        BA4_C9ALS_unweighted_df_summary <- BA4_C9ALS_unweighted_df %>%
            dplyr::group_by(group, celltype, n_genes, method) %>%
            dplyr::summarise(mean_accuracy_by_group = mean(mean_accuracy, na.rm = TRUE))
        BA4_C9ALS_unweighted_df_summary <- data.frame(BA4_C9ALS_unweighted_df_summary)

        ## sumarize SFTLD
        BA4_SFTLD_unweighted_df$disease_cat <- "SFTLD"
        BA4_SFTLD_unweighted_df_summary <- BA4_SFTLD_unweighted_df %>%
            dplyr::group_by(group, celltype, n_genes, method) %>%
            dplyr::summarise(mean_accuracy_by_group = mean(mean_accuracy, na.rm = TRUE))
        BA4_SFTLD_unweighted_df_summary <- data.frame(BA4_SFTLD_unweighted_df_summary)

        ## summarize C9FTLD
        #BA4_C9FTLD_unweighted_df$disease_cat <- "C9FTLD"
        #BA4_C9FTLD_unweighted_df_summary <- BA4_C9FTLD_unweighted_df %>%
        #    dplyr::group_by(group, celltype, n_genes, method) %>%
        #    dplyr::summarise(mean_accuracy_by_group = mean(mean_accuracy, na.rm = TRUE))
        #BA4_C9FTLD_unweighted_df_summary <- data.frame(BA4_C9FTLD_unweighted_df_summary)

        ## summarize all conditions
        #BA4_all_conditions_unweighted_df$disease_cat <- "All"
        #BA4_all_conditions_unweighted_df_summary <- BA4_all_conditions_unweighted_df %>%
        #    dplyr::group_by(group, celltype, n_genes, method) %>%
        #    dplyr::summarise(mean_accuracy_by_group = mean(mean_accuracy, na.rm = TRUE))
        #BA4_all_conditions_unweighted_df_summary <- data.frame(BA4_all_conditions_unweighted_df_summary)

        
        #BA4_unweighted_counts <- rbind(BA4_SALS_unweighted_df_summary, BA4_C9ALS_unweighted_df_summary, BA4_SFTLD_unweighted_df_summary, BA4_C9FTLD_unweighted_df_summary, BA4_all_conditions_unweighted_df_summary)
        BA4_unweighted_permutation <- rbind(BA4_SALS_unweighted_df_summary, BA4_C9ALS_unweighted_df_summary, BA4_SFTLD_unweighted_df_summary)
        
        BA4_unweighted_permutation$group <- factor(BA4_unweighted_permutation$group, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All"))
        BA4_unweighted_permutation$method <- factor(BA4_unweighted_permutation$method, levels = c("LIME", "random_HVG", "random_gene"))
        BA4_unweighted_permutation$celltype <- factor(BA4_unweighted_permutation$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))
        
        ggplot(BA4_unweighted_permutation, aes(x = n_genes, y = mean_accuracy_by_group, colour = method, group = method)) + 
        theme_bw() + 
        geom_point() +
        geom_line() +
        scale_x_reverse() +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(colour = "black", angle = 90, hjust = 1, vjust = 0.5),
            axis.text.y = element_text(colour = "black"),
            axis.title.y = element_text(face = "bold"),
            axis.title.x = element_blank(),
        ) +
        facet_grid(group~celltype, scales = "free_x") +
        ylab("N genes\n(Z-score > 1)") +
        scale_colour_manual(values = c("orange", "seagreen", "darkblue")) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 5, width = 17)
    ##
##


    ###########################
    ## SALS BA4 -- Compute optimal and Plot 
    ###########################
    ## code
        ## sumarize SALS
        BA4_SALS_weighted_df$disease_cat <- "SALS"
        BA4_SALS_weighted_df_summary <- BA4_SALS_weighted_df %>%
            dplyr::group_by(group, celltype, n_genes, method) %>%
            dplyr::summarise(mean_accuracy_by_group = mean(mean_accuracy, na.rm = TRUE))
        BA4_SALS_weighted_df_summary <- data.frame(BA4_SALS_weighted_df_summary)

        spread_data <- BA4_SALS_weighted_df_summary %>%
            spread(key = method, value = mean_accuracy_by_group)

        spread_data <- spread_data %>%
            mutate(difference = LIME - random_HVG)

        max_diff_data <- spread_data %>%
            group_by(celltype) %>%
            filter(abs(difference) == max(abs(difference))) %>%
            select(celltype, n_genes, difference) %>% as.data.frame()

        colnames(max_diff_data) <- c("celltype", "optimal_genes", "difference") 
        nrow(BA4_SALS_weighted_df_summary)
        BA4_SALS_weighted_df_summary <- merge(BA4_SALS_weighted_df_summary, max_diff_data, by = "celltype")

        ## plot
        BA4_SALS_weighted_df_summary$method <- factor(BA4_SALS_weighted_df_summary$method, levels = c("LIME", "random_HVG", "random_gene"))
        BA4_SALS_weighted_df_summary$celltype <- factor(BA4_SALS_weighted_df_summary$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))
        
        
        ## custom label
        BA4_SALS_weighted_df_summary <- BA4_SALS_weighted_df_summary %>%
            mutate(label = paste(celltype, "\n(", optimal_genes, " genes)", sep = ""))

        BA4_SALS_weighted_df_summary <- BA4_SALS_weighted_df_summary %>%
            arrange(celltype)  # Arrange by the 'celltype' factor

        BA4_SALS_weighted_df_summary$label <- factor(BA4_SALS_weighted_df_summary$label, levels = unique(BA4_SALS_weighted_df_summary$label))

        ggplot(BA4_SALS_weighted_df_summary, aes(x = n_genes, y = mean_accuracy_by_group, colour = method, group = method)) + 
        theme_bw() + 
        geom_vline(aes(xintercept = optimal_genes), linetype = "dashed", color = "black") +
        geom_point() +
        geom_line() +
        scale_x_reverse() +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(colour = "black", angle = 90, hjust = 1, vjust = 0.5),
            axis.text.y = element_text(colour = "black"),
            axis.title.y = element_text(face = "bold"),
            axis.title.x = element_text(face = "bold"),
        ) +
        facet_grid(group ~ label, scales = "free_x", labeller = labeller(celltype = custom_labeller)) + 
        ylab("Accuracy") +
        xlab("N genes") +
        scale_colour_manual(values = c("#F18F01","#048BA8" , "#2E4057"))
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 2.5, width = 15)
    ##

    ###########################
    ## C9ALS BA4 -- Compute optimal and Plot 
    ###########################
    ## code
        ## sumarize C9ALS
        BA4_C9ALS_weighted_df$disease_cat <- "C9ALS"
        BA4_C9ALS_weighted_df_summary <- BA4_C9ALS_weighted_df %>%
            dplyr::group_by(group, celltype, n_genes, method) %>%
            dplyr::summarise(mean_accuracy_by_group = mean(mean_accuracy, na.rm = TRUE))
        BA4_C9ALS_weighted_df_summary <- data.frame(BA4_C9ALS_weighted_df_summary)

        spread_data <- BA4_C9ALS_weighted_df_summary %>%
            spread(key = method, value = mean_accuracy_by_group)

        spread_data <- spread_data %>%
            mutate(difference = LIME - random_HVG)

        max_diff_data <- spread_data %>%
            group_by(celltype) %>%
            filter(abs(difference) == max(abs(difference))) %>%
            select(celltype, n_genes, difference) %>% as.data.frame()

        colnames(max_diff_data) <- c("celltype", "optimal_genes", "difference") 
        nrow(BA4_C9ALS_weighted_df_summary)
        BA4_C9ALS_weighted_df_summary <- merge(BA4_C9ALS_weighted_df_summary, max_diff_data, by = "celltype")

        ## plot
        BA4_C9ALS_weighted_df_summary$method <- factor(BA4_C9ALS_weighted_df_summary$method, levels = c("LIME", "random_HVG", "random_gene"))
        BA4_C9ALS_weighted_df_summary$celltype <- factor(BA4_C9ALS_weighted_df_summary$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))
        
        
        ## custom label
        BA4_C9ALS_weighted_df_summary <- BA4_C9ALS_weighted_df_summary %>%
            mutate(label = paste(celltype, "\n(", optimal_genes, " genes)", sep = ""))

        BA4_C9ALS_weighted_df_summary <- BA4_C9ALS_weighted_df_summary %>%
            arrange(celltype)  # Arrange by the 'celltype' factor

        BA4_C9ALS_weighted_df_summary$label <- factor(BA4_C9ALS_weighted_df_summary$label, levels = unique(BA4_C9ALS_weighted_df_summary$label))

        ggplot(BA4_C9ALS_weighted_df_summary, aes(x = n_genes, y = mean_accuracy_by_group, colour = method, group = method)) + 
        theme_bw() + 
        geom_vline(aes(xintercept = optimal_genes), linetype = "dashed", color = "black") +
        geom_point() +
        geom_line() +
        scale_x_reverse() +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(colour = "black", angle = 90, hjust = 1, vjust = 0.5),
            axis.text.y = element_text(colour = "black"),
            axis.title.y = element_text(face = "bold"),
            axis.title.x = element_text(face = "bold"),
        ) +
        facet_grid(group ~ label, scales = "free_x", labeller = labeller(celltype = custom_labeller)) + 
        ylab("Accuracy") +
        xlab("N genes") +
        scale_colour_manual(values = c("#F18F01","#048BA8" , "#2E4057"))
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 2.5, width = 15)
    ##

    ###########################
    ## SFTLD BA4 -- Compute optimal and Plot 
    ###########################
    ## code
        ## sumarize SFTLD
        BA4_SFTLD_weighted_df$disease_cat <- "SFTLD"
        BA4_SFTLD_weighted_df_summary <- BA4_SFTLD_weighted_df %>%
            dplyr::group_by(group, celltype, n_genes, method) %>%
            dplyr::summarise(mean_accuracy_by_group = mean(mean_accuracy, na.rm = TRUE))
        BA4_SFTLD_weighted_df_summary <- data.frame(BA4_SFTLD_weighted_df_summary)

        spread_data <- BA4_SFTLD_weighted_df_summary %>%
            spread(key = method, value = mean_accuracy_by_group)

        spread_data <- spread_data %>%
            mutate(difference = LIME - random_HVG)

        max_diff_data <- spread_data %>%
            group_by(celltype) %>%
            filter(abs(difference) == max(abs(difference))) %>%
            select(celltype, n_genes, difference) %>% as.data.frame()

        colnames(max_diff_data) <- c("celltype", "optimal_genes", "difference") 
        nrow(BA4_SFTLD_weighted_df_summary)
        BA4_SFTLD_weighted_df_summary <- merge(BA4_SFTLD_weighted_df_summary, max_diff_data, by = "celltype")

        ## plot
        BA4_SFTLD_weighted_df_summary$method <- factor(BA4_SFTLD_weighted_df_summary$method, levels = c("LIME", "random_HVG", "random_gene"))
        BA4_SFTLD_weighted_df_summary$celltype <- factor(BA4_SFTLD_weighted_df_summary$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))
        
        
        ## custom label
        BA4_SFTLD_weighted_df_summary <- BA4_SFTLD_weighted_df_summary %>%
            mutate(label = paste(celltype, "\n(", optimal_genes, " genes)", sep = ""))

        BA4_SFTLD_weighted_df_summary <- BA4_SFTLD_weighted_df_summary %>%
            arrange(celltype)  # Arrange by the 'celltype' factor

        BA4_SFTLD_weighted_df_summary$label <- factor(BA4_SFTLD_weighted_df_summary$label, levels = unique(BA4_SFTLD_weighted_df_summary$label))

        ggplot(BA4_SFTLD_weighted_df_summary, aes(x = n_genes, y = mean_accuracy_by_group, colour = method, group = method)) + 
        theme_bw() + 
        geom_vline(aes(xintercept = optimal_genes), linetype = "dashed", color = "black") +
        geom_point() +
        geom_line() +
        scale_x_reverse() +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(colour = "black", angle = 90, hjust = 1, vjust = 0.5),
            axis.text.y = element_text(colour = "black"),
            axis.title.y = element_text(face = "bold"),
            axis.title.x = element_text(face = "bold"),
        ) +
        facet_grid(group ~ label, scales = "free_x", labeller = labeller(celltype = custom_labeller)) + 
        ylab("Accuracy") +
        xlab("N genes") +
        scale_colour_manual(values = c("#F18F01","#048BA8" , "#2E4057"))
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 2.5, width = 15)
    ##

    ###########################
    ## C9FTLD BA4 -- Compute optimal and Plot 
    ###########################
    ## code
        ## sumarize C9FTLD
        BA4_C9FTLD_weighted_df$disease_cat <- "C9FTLD"
        BA4_C9FTLD_weighted_df_summary <- BA4_C9FTLD_weighted_df %>%
            dplyr::group_by(group, celltype, n_genes, method) %>%
            dplyr::summarise(mean_accuracy_by_group = mean(mean_accuracy, na.rm = TRUE))
        BA4_C9FTLD_weighted_df_summary <- data.frame(BA4_C9FTLD_weighted_df_summary)

        spread_data <- BA4_C9FTLD_weighted_df_summary %>%
            spread(key = method, value = mean_accuracy_by_group)

        spread_data <- spread_data %>%
            mutate(difference = LIME - random_HVG)

        max_diff_data <- spread_data %>%
            group_by(celltype) %>%
            filter(abs(difference) == max(abs(difference))) %>%
            select(celltype, n_genes, difference) %>% as.data.frame()

        colnames(max_diff_data) <- c("celltype", "optimal_genes", "difference") 
        nrow(BA4_C9FTLD_weighted_df_summary)
        BA4_C9FTLD_weighted_df_summary <- merge(BA4_C9FTLD_weighted_df_summary, max_diff_data, by = "celltype")

        ## plot
        BA4_C9FTLD_weighted_df_summary$method <- factor(BA4_C9FTLD_weighted_df_summary$method, levels = c("LIME", "random_HVG", "random_gene"))
        BA4_C9FTLD_weighted_df_summary$celltype <- factor(BA4_C9FTLD_weighted_df_summary$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))
        
        
        ## custom label
        BA4_C9FTLD_weighted_df_summary <- BA4_C9FTLD_weighted_df_summary %>%
            mutate(label = paste(celltype, "\n(", optimal_genes, " genes)", sep = ""))

        BA4_C9FTLD_weighted_df_summary <- BA4_C9FTLD_weighted_df_summary %>%
            arrange(celltype)  # Arrange by the 'celltype' factor

        BA4_C9FTLD_weighted_df_summary$label <- factor(BA4_C9FTLD_weighted_df_summary$label, levels = unique(BA4_C9FTLD_weighted_df_summary$label))

        ggplot(BA4_C9FTLD_weighted_df_summary, aes(x = n_genes, y = mean_accuracy_by_group, colour = method, group = method)) + 
        theme_bw() + 
        geom_vline(aes(xintercept = optimal_genes), linetype = "dashed", color = "black") +
        geom_point() +
        geom_line() +
        scale_x_reverse() +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(colour = "black", angle = 90, hjust = 1, vjust = 0.5),
            axis.text.y = element_text(colour = "black"),
            axis.title.y = element_text(face = "bold"),
            axis.title.x = element_text(face = "bold"),
        ) +
        facet_grid(group ~ label, scales = "free_x", labeller = labeller(celltype = custom_labeller)) + 
        ylab("Accuracy") +
        xlab("N genes") +
        scale_colour_manual(values = c("#F18F01","#048BA8" , "#2E4057"))
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 2.5, width = 15)
    ##


    ###########################
    ## SALS BA4 -- Compute optimal and Plot -- test to maximize the difference and the LIME accuracy. 

    #1. Normalize the Columns: Standardize the LIME and difference columns so they are on the same scale (e.g., range from 0 to 1). You can use z-scores or min-max normalization.
    #2. Combine the Columns: Create a new score that combines LIME and difference. For example, you could sum the normalized LIME and difference values or compute a weighted average.
    #3. Find the Optimal n_genes: For each celltype, identify the n_genes value that maximizes the combined score.
    ###########################
    ## code
        ## sumarize SALS
        BA4_SALS_weighted_df$disease_cat <- "SALS"
        BA4_SALS_weighted_df_summary <- BA4_SALS_weighted_df %>%
            dplyr::group_by(group, celltype, n_genes, method) %>%
            dplyr::summarise(mean_accuracy_by_group = mean(mean_accuracy, na.rm = TRUE))
        BA4_SALS_weighted_df_summary <- data.frame(BA4_SALS_weighted_df_summary)

        spread_data <- BA4_SALS_weighted_df_summary %>%
            spread(key = method, value = mean_accuracy_by_group)

        spread_data <- spread_data %>%
            mutate(difference = LIME - random_HVG)

        
        # Normalize the LIME and difference columns
        spread_data <- spread_data %>%
        mutate(
            LIME_normalized = (LIME - min(LIME)) / (max(LIME) - min(LIME)),  # Min-max normalization for LIME
            difference_normalized = (difference - min(difference)) / (max(difference) - min(difference))  # Min-max normalization for difference
        )

        # Create a new column that combines LIME and difference
        # You can adjust the weight for each column (e.g., 0.5 for equal weighting)
        spread_data <- spread_data %>%
        mutate(
            combined_score = 1.25*LIME_normalized + difference_normalized  # Equal weight, you can adjust this
        )

        max_diff_data <- spread_data %>%
            group_by(celltype) %>%
            filter(combined_score == max(combined_score)) %>%
            select(celltype, n_genes, combined_score) %>% as.data.frame()

        colnames(max_diff_data) <- c("celltype", "optimal_genes", "combined_score") 
        
        nrow(BA4_SALS_weighted_df_summary)
        BA4_SALS_weighted_df_summary <- merge(BA4_SALS_weighted_df_summary, max_diff_data, by = "celltype")

        ## plot
        BA4_SALS_weighted_df_summary$method <- factor(BA4_SALS_weighted_df_summary$method, levels = c("LIME", "random_HVG", "random_gene"))
        BA4_SALS_weighted_df_summary$celltype <- factor(BA4_SALS_weighted_df_summary$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))
        
        
        ## custom label
        BA4_SALS_weighted_df_summary <- BA4_SALS_weighted_df_summary %>%
            mutate(label = paste(celltype, "\n(", optimal_genes, " genes)", sep = ""))

        BA4_SALS_weighted_df_summary <- BA4_SALS_weighted_df_summary %>%
            arrange(celltype)  # Arrange by the 'celltype' factor

        BA4_SALS_weighted_df_summary$label <- factor(BA4_SALS_weighted_df_summary$label, levels = unique(BA4_SALS_weighted_df_summary$label))

        ## plot
        ggplot(BA4_SALS_weighted_df_summary, aes(x = n_genes, y = mean_accuracy_by_group, colour = method, group = method)) + 
        theme_bw() + 
        geom_vline(aes(xintercept = optimal_genes), linetype = "dashed", color = "black") +
        geom_point() +
        geom_line() +
        scale_x_reverse() +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(colour = "black", angle = 90, hjust = 1, vjust = 0.5),
            axis.text.y = element_text(colour = "black"),
            axis.title.y = element_text(face = "bold"),
            axis.title.x = element_text(face = "bold"),
        ) +
        facet_grid(group ~ label, scales = "free_x", labeller = labeller(celltype = custom_labeller)) + 
        ylab("Accuracy") +
        xlab("N genes") +
        scale_colour_manual(values = c("#F18F01","#048BA8" , "#2E4057"))
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 2.5, width = 15)
    ##

    ###########################
    ## C9ALS BA4 -- Compute optimal and Plot -- test to maximize the difference and the LIME accuracy. 

    #1. Normalize the Columns: Standardize the LIME and difference columns so they are on the same scale (e.g., range from 0 to 1). You can use z-scores or min-max normalization.
    #2. Combine the Columns: Create a new score that combines LIME and difference. For example, you could sum the normalized LIME and difference values or compute a weighted average.
    #3. Find the Optimal n_genes: For each celltype, identify the n_genes value that maximizes the combined score.
    ###########################
    ## code
        ## sumarize C9ALS
        BA4_C9ALS_weighted_df$disease_cat <- "C9ALS"
        BA4_C9ALS_weighted_df_summary <- BA4_C9ALS_weighted_df %>%
            dplyr::group_by(group, celltype, n_genes, method) %>%
            dplyr::summarise(mean_accuracy_by_group = mean(mean_accuracy, na.rm = TRUE))
        BA4_C9ALS_weighted_df_summary <- data.frame(BA4_C9ALS_weighted_df_summary)

        spread_data <- BA4_C9ALS_weighted_df_summary %>%
            spread(key = method, value = mean_accuracy_by_group)

        spread_data <- spread_data %>%
            mutate(difference = LIME - random_HVG)

        
        # Normalize the LIME and difference columns
        spread_data <- spread_data %>%
        mutate(
            LIME_normalized = (LIME - min(LIME)) / (max(LIME) - min(LIME)),  # Min-max normalization for LIME
            difference_normalized = (difference - min(difference)) / (max(difference) - min(difference))  # Min-max normalization for difference
        )

        # Create a new column that combines LIME and difference
        # You can adjust the weight for each column (e.g., 0.5 for equal weighting)
        spread_data <- spread_data %>%
        mutate(
            combined_score = 2*LIME_normalized + difference_normalized  # Equal weight, you can adjust this
        )

        max_diff_data <- spread_data %>%
            group_by(celltype) %>%
            filter(combined_score == max(combined_score)) %>%
            select(celltype, n_genes, combined_score) %>% as.data.frame()

        colnames(max_diff_data) <- c("celltype", "optimal_genes", "combined_score") 
        
        nrow(BA4_C9ALS_weighted_df_summary)
        BA4_C9ALS_weighted_df_summary <- merge(BA4_C9ALS_weighted_df_summary, max_diff_data, by = "celltype")

        ## plot
        BA4_C9ALS_weighted_df_summary$method <- factor(BA4_C9ALS_weighted_df_summary$method, levels = c("LIME", "random_HVG", "random_gene"))
        BA4_C9ALS_weighted_df_summary$celltype <- factor(BA4_C9ALS_weighted_df_summary$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))
        
        
        ## custom label
        BA4_C9ALS_weighted_df_summary <- BA4_C9ALS_weighted_df_summary %>%
            mutate(label = paste(celltype, "\n(", optimal_genes, " genes)", sep = ""))

        BA4_C9ALS_weighted_df_summary <- BA4_C9ALS_weighted_df_summary %>%
            arrange(celltype)  # Arrange by the 'celltype' factor

        BA4_C9ALS_weighted_df_summary$label <- factor(BA4_C9ALS_weighted_df_summary$label, levels = unique(BA4_C9ALS_weighted_df_summary$label))

        ## plot
        ggplot(BA4_C9ALS_weighted_df_summary, aes(x = n_genes, y = mean_accuracy_by_group, colour = method, group = method)) + 
        theme_bw() + 
        geom_vline(aes(xintercept = optimal_genes), linetype = "dashed", color = "black") +
        geom_point() +
        geom_line() +
        scale_x_reverse() +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(colour = "black", angle = 90, hjust = 1, vjust = 0.5),
            axis.text.y = element_text(colour = "black"),
            axis.title.y = element_text(face = "bold"),
            axis.title.x = element_text(face = "bold"),
        ) +
        facet_grid(group ~ label, scales = "free_x", labeller = labeller(celltype = custom_labeller)) + 
        ylab("Accuracy") +
        xlab("N genes") +
        scale_colour_manual(values = c("#F18F01","#048BA8" , "#2E4057"))
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 2.5, width = 15)
    ##


    ###########################
    ## SFTLD BA4 -- Compute optimal and Plot -- test to maximize the difference and the LIME accuracy. 

    #1. Normalize the Columns: Standardize the LIME and difference columns so they are on the same scale (e.g., range from 0 to 1). You can use z-scores or min-max normalization.
    #2. Combine the Columns: Create a new score that combines LIME and difference. For example, you could sum the normalized LIME and difference values or compute a weighted average.
    #3. Find the Optimal n_genes: For each celltype, identify the n_genes value that maximizes the combined score.
    ###########################
    ## code
        ## sumarize SFTLD
        BA4_SFTLD_weighted_df$disease_cat <- "SFTLD"
        BA4_SFTLD_weighted_df_summary <- BA4_SFTLD_weighted_df %>%
            dplyr::group_by(group, celltype, n_genes, method) %>%
            dplyr::summarise(mean_accuracy_by_group = mean(mean_accuracy, na.rm = TRUE))
        BA4_SFTLD_weighted_df_summary <- data.frame(BA4_SFTLD_weighted_df_summary)

        spread_data <- BA4_SFTLD_weighted_df_summary %>%
            spread(key = method, value = mean_accuracy_by_group)

        spread_data <- spread_data %>%
            mutate(difference = LIME - random_HVG)

        
        # Normalize the LIME and difference columns
        spread_data <- spread_data %>%
        mutate(
            LIME_normalized = (LIME - min(LIME)) / (max(LIME) - min(LIME)),  # Min-max normalization for LIME
            difference_normalized = (difference - min(difference)) / (max(difference) - min(difference))  # Min-max normalization for difference
        )

        # Create a new column that combines LIME and difference
        # You can adjust the weight for each column (e.g., 0.5 for equal weighting)
        spread_data <- spread_data %>%
        mutate(
            combined_score = LIME_normalized + 2*difference_normalized  # Equal weight, you can adjust this
        )

        max_diff_data <- spread_data %>%
            group_by(celltype) %>%
            filter(combined_score == max(combined_score)) %>%
            select(celltype, n_genes, combined_score) %>% as.data.frame()

        colnames(max_diff_data) <- c("celltype", "optimal_genes", "combined_score") 
        
        nrow(BA4_SFTLD_weighted_df_summary)
        BA4_SFTLD_weighted_df_summary <- merge(BA4_SFTLD_weighted_df_summary, max_diff_data, by = "celltype")

        ## plot
        BA4_SFTLD_weighted_df_summary$method <- factor(BA4_SFTLD_weighted_df_summary$method, levels = c("LIME", "random_HVG", "random_gene"))
        BA4_SFTLD_weighted_df_summary$celltype <- factor(BA4_SFTLD_weighted_df_summary$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))
        
        
        ## custom label
        BA4_SFTLD_weighted_df_summary <- BA4_SFTLD_weighted_df_summary %>%
            mutate(label = paste(celltype, "\n(", optimal_genes, " genes)", sep = ""))

        BA4_SFTLD_weighted_df_summary <- BA4_SFTLD_weighted_df_summary %>%
            arrange(celltype)  # Arrange by the 'celltype' factor

        BA4_SFTLD_weighted_df_summary$label <- factor(BA4_SFTLD_weighted_df_summary$label, levels = unique(BA4_SFTLD_weighted_df_summary$label))

        ## plot
        ggplot(BA4_SFTLD_weighted_df_summary, aes(x = n_genes, y = mean_accuracy_by_group, colour = method, group = method)) + 
        theme_bw() + 
        geom_vline(aes(xintercept = optimal_genes), linetype = "dashed", color = "black") +
        geom_point() +
        geom_line() +
        scale_x_reverse() +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(colour = "black", angle = 90, hjust = 1, vjust = 0.5),
            axis.text.y = element_text(colour = "black"),
            axis.title.y = element_text(face = "bold"),
            axis.title.x = element_text(face = "bold"),
        ) +
        facet_grid(group ~ label, scales = "free_x", labeller = labeller(celltype = custom_labeller)) + 
        ylab("Accuracy") +
        xlab("N genes") +
        scale_colour_manual(values = c("#F18F01","#048BA8" , "#2E4057"))
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 2.5, width = 15)
    ##
        ## sumarize C9ALS
        BA4_C9ALS_weighted_df$disease_cat <- "C9ALS"
        BA4_C9ALS_weighted_df_summary <- BA4_C9ALS_weighted_df %>%
            dplyr::group_by(group, celltype, n_genes, method) %>%
            dplyr::summarise(mean_accuracy_by_group = mean(mean_accuracy, na.rm = TRUE))
        BA4_C9ALS_weighted_df_summary <- data.frame(BA4_C9ALS_weighted_df_summary)

        ## sumarize SFTLD
        BA4_SFTLD_weighted_df$disease_cat <- "SFTLD"
        BA4_SFTLD_weighted_df_summary <- BA4_SFTLD_weighted_df %>%
            dplyr::group_by(group, celltype, n_genes, method) %>%
            dplyr::summarise(mean_accuracy_by_group = mean(mean_accuracy, na.rm = TRUE))
        BA4_SFTLD_weighted_df_summary <- data.frame(BA4_SFTLD_weighted_df_summary)

        ## summarize C9FTLD
        BA4_C9FTLD_weighted_df$disease_cat <- "C9FTLD"
        BA4_C9FTLD_weighted_df_summary <- BA4_C9FTLD_weighted_df %>%
            dplyr::group_by(group, celltype, n_genes, method) %>%
            dplyr::summarise(mean_accuracy_by_group = mean(mean_accuracy, na.rm = TRUE))
        BA4_C9FTLD_weighted_df_summary <- data.frame(BA4_C9FTLD_weighted_df_summary)

        ## summarize all conditions
        #BA4_all_conditions_weighted_df$disease_cat <- "All"
        #BA4_all_conditions_weighted_df_summary <- BA4_all_conditions_weighted_df %>%
        #    dplyr::group_by(group, celltype, n_genes, method) %>%
        #    dplyr::summarise(mean_accuracy_by_group = mean(mean_accuracy, na.rm = TRUE))
        #BA4_all_conditions_weighted_df_summary <- data.frame(BA4_all_conditions_weighted_df_summary)

        
        #BA4_weighted_counts <- rbind(BA4_SALS_weighted_df_summary, BA4_C9ALS_weighted_df_summary, BA4_SFTLD_weighted_df_summary, BA4_C9FTLD_weighted_df_summary, BA4_all_conditions_weighted_df_summary)
        BA4_weighted_permutation <- rbind(BA4_SALS_weighted_df_summary, BA4_C9ALS_weighted_df_summary, BA4_SFTLD_weighted_df_summary, BA4_C9FTLD_weighted_df_summary)
        
        BA4_weighted_permutation$group <- factor(BA4_weighted_permutation$group, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All"))
        BA4_weighted_permutation$method <- factor(BA4_weighted_permutation$method, levels = c("LIME", "random_HVG", "random_gene"))
        BA4_weighted_permutation$celltype <- factor(BA4_weighted_permutation$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))
        
        ggplot(BA4_weighted_permutation, aes(x = n_genes, y = mean_accuracy_by_group, colour = method, group = method)) + 
        theme_bw() + 
        geom_point() +
        geom_line() +
        scale_x_reverse() +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(colour = "black", angle = 90, hjust = 1, vjust = 0.5),
            axis.text.y = element_text(colour = "black"),
            axis.title.y = element_text(face = "bold"),
            axis.title.x = element_blank(),
        ) +
        facet_grid(group~celltype, scales = "free_x") +
        ylab("N genes\n(Z-score > 1)") +
        scale_colour_manual(values = c("orange", "seagreen", "darkblue")) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 5, width = 17)
    ##


##































######### OLD CODE BELOW HERE


############################################################################ Explore genes with LIME Z score > 1
############################################################################
############################################################################
############################################################################ Condition specific models

## code
    ########################
    ## Identify genes with LIME Z.score > 1
    ########################

    ## code SALS BA4
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA4"
        par_status = "SALS"
        par_prep = "CombatSeq"
        
        ################################################
        ## LIME outputs not weighted by percent expression
        ## This function subsets the unweighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_across_donors = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_across_donors, celltype)

        for(celltype2 in unique(celltype_list)){
            

            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_unweighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_across_donors <- as.numeric(df$Mean_z_score_across_donors)
            df <- df[order( -(df$Mean_z_score_across_donors)),]

            df <- subset(df, Mean_z_score_across_donors >= 1)
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_across_donors)
            df$celltype <- celltype2

            write.csv(df, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_unweighted_Z1_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")

            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")

        BA4_SALS_unweighted_df <- fill

        ################################################
        ## LIME outputs weighted by percent expression
        ## This function subsets the weighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_weighted_across_donor = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_weighted_across_donor, celltype)

        for(celltype2 in unique(celltype_list)){
            
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_weighted_across_donor <- as.numeric(df$Mean_z_score_weighted_across_donor)
            df <- df[order( -(df$Mean_z_score_weighted_across_donor)),]

            df <- subset(df, Mean_z_score_weighted_across_donor >= 1)
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_weighted_across_donor)
            df$celltype <- celltype2

            write.csv(df, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_Z1_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")

            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")
        BA4_SALS_weighted_df <- fill

    ##

    ## code C9ALS BA4
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA4"
        par_status = "C9ALS"
        par_prep = "CombatSeq"
        
        ################################################
        ## LIME outputs not weighted by percent expression
        ## This function subsets the unweighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_across_donors = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_across_donors, celltype)

        for(celltype2 in unique(celltype_list)){
            

            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_unweighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_across_donors <- as.numeric(df$Mean_z_score_across_donors)
            df <- df[order( -(df$Mean_z_score_across_donors)),]

            df <- subset(df, Mean_z_score_across_donors >= 1)
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_across_donors)
            df$celltype <- celltype2

            write.csv(df, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_unweighted_Z1_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")

            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")

        BA4_C9ALS_unweighted_df <- fill

        ################################################
        ## LIME outputs weighted by percent expression
        ## This function subsets the weighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_weighted_across_donor = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_weighted_across_donor, celltype)

        for(celltype2 in unique(celltype_list)){
            
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_weighted_across_donor <- as.numeric(df$Mean_z_score_weighted_across_donor)
            df <- df[order( -(df$Mean_z_score_weighted_across_donor)),]

            df <- subset(df, Mean_z_score_weighted_across_donor >= 1)
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_weighted_across_donor)
            df$celltype <- celltype2

            write.csv(df, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_Z1_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")

            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")
        BA4_C9ALS_weighted_df <- fill

    ##

    ## code SFTLD BA4
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA4"
        par_status = "SFTLD"
        par_prep = "CombatSeq"
        
        ################################################
        ## LIME outputs not weighted by percent expression
        ## This function subsets the unweighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_across_donors = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_across_donors, celltype)

        for(celltype2 in unique(celltype_list)){
            

            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_unweighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_across_donors <- as.numeric(df$Mean_z_score_across_donors)
            df <- df[order( -(df$Mean_z_score_across_donors)),]

            df <- subset(df, Mean_z_score_across_donors >= 1)
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_across_donors)
            df$celltype <- celltype2

            write.csv(df, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_unweighted_Z1_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")

            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")

        BA4_SFTLD_unweighted_df <- fill

        ################################################
        ## LIME outputs weighted by percent expression
        ## This function subsets the weighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_weighted_across_donor = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_weighted_across_donor, celltype)

        for(celltype2 in unique(celltype_list)){
            
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_weighted_across_donor <- as.numeric(df$Mean_z_score_weighted_across_donor)
            df <- df[order( -(df$Mean_z_score_weighted_across_donor)),]

            df <- subset(df, Mean_z_score_weighted_across_donor >= 1)
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_weighted_across_donor)
            df$celltype <- celltype2

            write.csv(df, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_Z1_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")

            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")
        BA4_SFTLD_weighted_df <- fill

    ##

    ## code C9FTLD BA4
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA4"
        par_status = "C9FTLD"
        par_prep = "CombatSeq"
        
        ################################################
        ## LIME outputs not weighted by percent expression
        ## This function subsets the unweighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_across_donors = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_across_donors, celltype)

        for(celltype2 in unique(celltype_list)){
            

            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_unweighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_across_donors <- as.numeric(df$Mean_z_score_across_donors)
            df <- df[order( -(df$Mean_z_score_across_donors)),]

            df <- subset(df, Mean_z_score_across_donors >= 1)
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_across_donors)
            df$celltype <- celltype2

            write.csv(df, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_unweighted_Z1_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")

            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")

        BA4_C9FTLD_unweighted_df <- fill

        ################################################
        ## LIME outputs weighted by percent expression
        ## This function subsets the weighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_weighted_across_donor = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_weighted_across_donor, celltype)

        for(celltype2 in unique(celltype_list)){
            
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_weighted_across_donor <- as.numeric(df$Mean_z_score_weighted_across_donor)
            df <- df[order( -(df$Mean_z_score_weighted_across_donor)),]

            df <- subset(df, Mean_z_score_weighted_across_donor >= 1)
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_weighted_across_donor)
            df$celltype <- celltype2

            write.csv(df, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_Z1_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")

            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")
        BA4_C9FTLD_weighted_df <- fill

    ##

    ## code SALS BA9
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA9"
        par_status = "SALS"
        par_prep = "CombatSeq"
        
        ################################################
        ## LIME outputs not weighted by percent expression
        ## This function subsets the unweighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_across_donors = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_across_donors, celltype)

        for(celltype2 in unique(celltype_list)){
            

            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_unweighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_across_donors <- as.numeric(df$Mean_z_score_across_donors)
            df <- df[order( -(df$Mean_z_score_across_donors)),]

            df <- subset(df, Mean_z_score_across_donors >= 1)
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_across_donors)
            df$celltype <- celltype2

            write.csv(df, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_unweighted_Z1_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")

            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")

        BA9_SALS_unweighted_df <- fill

        ################################################
        ## LIME outputs weighted by percent expression
        ## This function subsets the weighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_weighted_across_donor = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_weighted_across_donor, celltype)

        for(celltype2 in unique(celltype_list)){
            
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_weighted_across_donor <- as.numeric(df$Mean_z_score_weighted_across_donor)
            df <- df[order( -(df$Mean_z_score_weighted_across_donor)),]

            df <- subset(df, Mean_z_score_weighted_across_donor >= 1)
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_weighted_across_donor)
            df$celltype <- celltype2

            write.csv(df, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_Z1_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")

            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")
        BA9_SALS_weighted_df <- fill

    ##

    ## code C9ALS BA9
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA9"
        par_status = "C9ALS"
        par_prep = "CombatSeq"
        
        ################################################
        ## LIME outputs not weighted by percent expression
        ## This function subsets the unweighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_across_donors = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_across_donors, celltype)

        for(celltype2 in unique(celltype_list)){
            

            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_unweighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_across_donors <- as.numeric(df$Mean_z_score_across_donors)
            df <- df[order( -(df$Mean_z_score_across_donors)),]

            df <- subset(df, Mean_z_score_across_donors >= 1)
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_across_donors)
            df$celltype <- celltype2

            write.csv(df, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_unweighted_Z1_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")

            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")

        BA9_C9ALS_unweighted_df <- fill

        ################################################
        ## LIME outputs weighted by percent expression
        ## This function subsets the weighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_weighted_across_donor = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_weighted_across_donor, celltype)

        for(celltype2 in unique(celltype_list)){
            
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_weighted_across_donor <- as.numeric(df$Mean_z_score_weighted_across_donor)
            df <- df[order( -(df$Mean_z_score_weighted_across_donor)),]

            df <- subset(df, Mean_z_score_weighted_across_donor >= 1)
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_weighted_across_donor)
            df$celltype <- celltype2

            write.csv(df, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_Z1_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")

            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")
        BA9_C9ALS_weighted_df <- fill

    ##

    ## code SFTLD BA9
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA9"
        par_status = "SFTLD"
        par_prep = "CombatSeq"
        
        ################################################
        ## LIME outputs not weighted by percent expression
        ## This function subsets the unweighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_across_donors = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_across_donors, celltype)

        for(celltype2 in unique(celltype_list)){
            

            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_unweighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_across_donors <- as.numeric(df$Mean_z_score_across_donors)
            df <- df[order( -(df$Mean_z_score_across_donors)),]

            df <- subset(df, Mean_z_score_across_donors >= 1)
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_across_donors)
            df$celltype <- celltype2

            write.csv(df, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_unweighted_Z1_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")

            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")

        BA9_SFTLD_unweighted_df <- fill

        ################################################
        ## LIME outputs weighted by percent expression
        ## This function subsets the weighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_weighted_across_donor = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_weighted_across_donor, celltype)

        for(celltype2 in unique(celltype_list)){
            
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_weighted_across_donor <- as.numeric(df$Mean_z_score_weighted_across_donor)
            df <- df[order( -(df$Mean_z_score_weighted_across_donor)),]

            df <- subset(df, Mean_z_score_weighted_across_donor >= 1)
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_weighted_across_donor)
            df$celltype <- celltype2

            write.csv(df, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_Z1_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")

            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")
        BA9_SFTLD_weighted_df <- fill

    ##

    ## code C9FTLD BA9 -- missing TCELL for some reason -- look into it
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA9"
        par_status = "C9FTLD"
        par_prep = "CombatSeq"
        
        ################################################
        ## LIME outputs not weighted by percent expression
        ## This function subsets the unweighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_across_donors = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_across_donors, celltype)

        for(celltype2 in unique(celltype_list)){
            

            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_unweighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_across_donors <- as.numeric(df$Mean_z_score_across_donors)
            df <- df[order( -(df$Mean_z_score_across_donors)),]

            df <- subset(df, Mean_z_score_across_donors >= 1)
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_across_donors)
            df$celltype <- celltype2

            write.csv(df, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_unweighted_Z1_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")

            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")

        BA9_C9FTLD_unweighted_df <- fill

        ################################################
        ## LIME outputs weighted by percent expression
        ## This function subsets the weighted LIME outputs to only include genes with Z score > 1. 
        ################################################
        ## Create a filler function
        feature = "fill"
        Mean_z_score_weighted_across_donor = 0
        celltype = "fill"
        fill <- data.frame(feature, Mean_z_score_weighted_across_donor, celltype)

        for(celltype2 in unique(celltype_list)){
            
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_weighted_across_donor <- as.numeric(df$Mean_z_score_weighted_across_donor)
            df <- df[order( -(df$Mean_z_score_weighted_across_donor)),]

            df <- subset(df, Mean_z_score_weighted_across_donor >= 1)
            nrow(df)
            df <- df %>% dplyr::select(feature, Mean_z_score_weighted_across_donor)
            df$celltype <- celltype2

            write.csv(df, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_Z1_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")

            fill <- rbind(fill, df)

        }

        fill <- subset(fill, feature != "fill")
        BA9_C9FTLD_weighted_df <- fill

    ##

    ########################
    ## Plot number of genes above Z > 1 for each 
    ########################

    ## BA4 unweighted
        BA4_SALS_unweighted_df$disease_cat <- "SALS"
        BA4_C9ALS_unweighted_df$disease_cat <- "C9ALS"
        BA4_SFTLD_unweighted_df$disease_cat <- "SFTLD"
        BA4_C9FTLD_unweighted_df$disease_cat <- "C9FTLD"

        BA4_unweighted <- rbind(BA4_SALS_unweighted_df, BA4_C9ALS_unweighted_df, BA4_SFTLD_unweighted_df, BA4_C9FTLD_unweighted_df)

        counts_df <- data.frame(table(BA4_unweighted$celltype, BA4_unweighted$disease_cat ))

        counts_df$Var2 <- factor(counts_df$Var2, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD"))
        counts_df$Var1 <- factor(counts_df$Var1, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))


        ggplot(counts_df, aes(x = Var2, y = Freq, fill = Var2)) + 
        theme_bw() + 
        geom_bar(stat = 'identity') +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~Var1, scales = "free_x") +
        ylab("N genes\n(Z-score > 1)") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) +
        ggtitle('BA4: n genes with LIME Z-score > 1 (Disease-specific models; unweighted LIME outputs)')
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 2, width = 17)
    ##

    ## BA9 unweighted
        BA9_SALS_unweighted_df$disease_cat <- "SALS"
        BA9_C9ALS_unweighted_df$disease_cat <- "C9ALS"
        BA9_SFTLD_unweighted_df$disease_cat <- "SFTLD"
        BA9_C9FTLD_unweighted_df$disease_cat <- "C9FTLD"

        BA9_unweighted <- rbind(BA9_SALS_unweighted_df, BA9_C9ALS_unweighted_df, BA9_SFTLD_unweighted_df, BA9_C9FTLD_unweighted_df)

        counts_df <- data.frame(table(BA9_unweighted$celltype, BA9_unweighted$disease_cat ))

        counts_df$Var2 <- factor(counts_df$Var2, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD"))
        counts_df$Var1 <- factor(counts_df$Var1, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))


        ggplot(counts_df, aes(x = Var2, y = Freq, fill = Var2)) + 
        theme_bw() + 
        geom_bar(stat = 'identity') +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~Var1, scales = "free_x") +
        ylab("N genes\n(Z-score > 1)") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) +
        ggtitle('BA9: n genes with LIME Z-score > 1 (Disease-specific models; unweighted LIME outputs)')
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 2, width = 17)
    ##

    ## BA4 weighted
        BA4_SALS_weighted_df$disease_cat <- "SALS"
        BA4_C9ALS_weighted_df$disease_cat <- "C9ALS"
        BA4_SFTLD_weighted_df$disease_cat <- "SFTLD"
        BA4_C9FTLD_weighted_df$disease_cat <- "C9FTLD"

        BA4_weighted <- rbind(BA4_SALS_weighted_df, BA4_C9ALS_weighted_df, BA4_SFTLD_weighted_df, BA4_C9FTLD_weighted_df)

        counts_df <- data.frame(table(BA4_weighted$celltype, BA4_weighted$disease_cat ))

        counts_df$Var2 <- factor(counts_df$Var2, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD"))
        counts_df$Var1 <- factor(counts_df$Var1, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))


        ggplot(counts_df, aes(x = Var2, y = Freq, fill = Var2)) + 
        theme_bw() + 
        geom_bar(stat = 'identity') +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~Var1, scales = "free_x") +
        ylab("N genes\n(Z-score > 1)") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) +
        ggtitle('BA4: n genes with LIME Z-score > 1 (Disease-specific models; weighted LIME outputs)')
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 2, width = 17)
    ##

    ## BA9 weighted
        BA9_SALS_weighted_df$disease_cat <- "SALS"
        BA9_C9ALS_weighted_df$disease_cat <- "C9ALS"
        BA9_SFTLD_weighted_df$disease_cat <- "SFTLD"
        BA9_C9FTLD_weighted_df$disease_cat <- "C9FTLD"

        BA9_weighted <- rbind(BA9_SALS_weighted_df, BA9_C9ALS_weighted_df, BA9_SFTLD_weighted_df, BA9_C9FTLD_weighted_df)

        counts_df <- data.frame(table(BA9_weighted$celltype, BA9_weighted$disease_cat ))

        counts_df$Var2 <- factor(counts_df$Var2, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD"))
        counts_df$Var1 <- factor(counts_df$Var1, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))


        ggplot(counts_df, aes(x = Var2, y = Freq, fill = Var2)) + 
        theme_bw() + 
        geom_bar(stat = 'identity') +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~Var1, scales = "free_x") +
        ylab("N genes\n(Z-score > 1)") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) +
        ggtitle('BA9: n genes with LIME Z-score > 1 (Disease-specific models; weighted LIME outputs)')
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 2, width = 17)
    ##

    ########################
    ## Plot number of genes above Z > 1 overlaping with NDKP HuGE genes
    ########################

    ## BA4 unweighted
        ## NDKP ALS HuGE -- gene
        NDKP_ALS_HuGE <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/gene_lists/NDKP_HuGE.csv'), sep = ",")
        nrow(NDKP_ALS_HuGE)

        BA4_SALS_unweighted_df$disease_cat <- "SALS"
        BA4_C9ALS_unweighted_df$disease_cat <- "C9ALS"
        BA4_SFTLD_unweighted_df$disease_cat <- "SFTLD"
        BA4_C9FTLD_unweighted_df$disease_cat <- "C9FTLD"

        BA4_unweighted <- rbind(BA4_SALS_unweighted_df, BA4_C9ALS_unweighted_df, BA4_SFTLD_unweighted_df, BA4_C9FTLD_unweighted_df)
        BA4_unweighted <- subset(BA4_unweighted, feature %in% NDKP_ALS_HuGE$gene)

        counts_df <- data.frame(table(BA4_unweighted$celltype, BA4_unweighted$disease_cat ))

        counts_df$Var2 <- factor(counts_df$Var2, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD"))
        counts_df$Var1 <- factor(counts_df$Var1, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))


        ggplot(counts_df, aes(x = Var2, y = Freq, fill = Var2)) + 
        theme_bw() + 
        geom_bar(stat = 'identity') +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~Var1, scales = "free_x") +
        ylab("N genes\n(Z-score > 1)") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) +
        ggtitle('BA4: n genes with LIME Z-score > 1 intersecting with NDKP ALS HuGE (Disease-specific models; unweighted LIME outputs)')
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 2, width = 17)
    ##

    ## BA9 unweighted
        ## NDKP ALS HuGE -- gene
        NDKP_ALS_HuGE <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/gene_lists/NDKP_HuGE.csv'), sep = ",")
        nrow(NDKP_ALS_HuGE)

        BA9_SALS_unweighted_df$disease_cat <- "SALS"
        BA9_C9ALS_unweighted_df$disease_cat <- "C9ALS"
        BA9_SFTLD_unweighted_df$disease_cat <- "SFTLD"
        BA9_C9FTLD_unweighted_df$disease_cat <- "C9FTLD"

        BA9_unweighted <- rbind(BA9_SALS_unweighted_df, BA9_C9ALS_unweighted_df, BA9_SFTLD_unweighted_df, BA9_C9FTLD_unweighted_df)
        BA9_unweighted <- subset(BA9_unweighted, feature %in% NDKP_ALS_HuGE$gene)

        counts_df <- data.frame(table(BA9_unweighted$celltype, BA9_unweighted$disease_cat ))

        counts_df$Var2 <- factor(counts_df$Var2, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD"))
        counts_df$Var1 <- factor(counts_df$Var1, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))


        ggplot(counts_df, aes(x = Var2, y = Freq, fill = Var2)) + 
        theme_bw() + 
        geom_bar(stat = 'identity') +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~Var1, scales = "free_x") +
        ylab("N genes\n(Z-score > 1)") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) +
        ggtitle('BA9: n genes with LIME Z-score > 1 intersecting with NDKP ALS HuGE (Disease-specific models; unweighted LIME outputs)')
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 2, width = 17)
    ##

    ## BA4 weighted
        ## NDKP ALS HuGE -- gene
        NDKP_ALS_HuGE <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/gene_lists/NDKP_HuGE.csv'), sep = ",")
        nrow(NDKP_ALS_HuGE)

        BA4_SALS_weighted_df$disease_cat <- "SALS"
        BA4_C9ALS_weighted_df$disease_cat <- "C9ALS"
        BA4_SFTLD_weighted_df$disease_cat <- "SFTLD"
        BA4_C9FTLD_weighted_df$disease_cat <- "C9FTLD"

        BA4_weighted <- rbind(BA4_SALS_weighted_df, BA4_C9ALS_weighted_df, BA4_SFTLD_weighted_df, BA4_C9FTLD_weighted_df)
        BA4_weighted <- subset(BA4_weighted, feature %in% NDKP_ALS_HuGE$gene)

        counts_df <- data.frame(table(BA4_weighted$celltype, BA4_weighted$disease_cat ))

        counts_df$Var2 <- factor(counts_df$Var2, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD"))
        counts_df$Var1 <- factor(counts_df$Var1, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))


        ggplot(counts_df, aes(x = Var2, y = Freq, fill = Var2)) + 
        theme_bw() + 
        geom_bar(stat = 'identity') +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~Var1, scales = "free_x") +
        ylab("N genes\n(Z-score > 1)") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) +
        ggtitle('BA4: n genes with LIME Z-score > 1 intersecting with NDKP ALS HuGE (Disease-specific models; weighted LIME outputs)')
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 2, width = 17)
    ##

    ## BA9 weighted
        ## NDKP ALS HuGE -- gene
        NDKP_ALS_HuGE <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/gene_lists/NDKP_HuGE.csv'), sep = ",")
        nrow(NDKP_ALS_HuGE)

        BA9_SALS_weighted_df$disease_cat <- "SALS"
        BA9_C9ALS_weighted_df$disease_cat <- "C9ALS"
        BA9_SFTLD_weighted_df$disease_cat <- "SFTLD"
        BA9_C9FTLD_weighted_df$disease_cat <- "C9FTLD"

        BA9_weighted <- rbind(BA9_SALS_weighted_df, BA9_C9ALS_weighted_df, BA9_SFTLD_weighted_df, BA9_C9FTLD_weighted_df)
        BA9_weighted <- subset(BA9_weighted, feature %in% NDKP_ALS_HuGE$gene)

        counts_df <- data.frame(table(BA9_weighted$celltype, BA9_weighted$disease_cat ))

        counts_df$Var2 <- factor(counts_df$Var2, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD"))
        counts_df$Var1 <- factor(counts_df$Var1, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))


        ggplot(counts_df, aes(x = Var2, y = Freq, fill = Var2)) + 
        theme_bw() + 
        geom_bar(stat = 'identity') +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~Var1, scales = "free_x") +
        ylab("N genes\n(Z-score > 1)") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) +
        ggtitle('BA9: n genes with LIME Z-score > 1 intersecting with NDKP ALS HuGE (Disease-specific models; weighted LIME outputs)')
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 2, width = 17)
    ##

    ########################
    ## Plot number of genes above Z > 1 overlaping with NDKP common 
    ########################

    ## BA4 unweighted
        ## NDKP ALS Common -- gene
        gene_list <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/gene_lists/NDKP_common.csv'), sep = ",")
        nrow(gene_list)

        BA4_SALS_unweighted_df$disease_cat <- "SALS"
        BA4_C9ALS_unweighted_df$disease_cat <- "C9ALS"
        BA4_SFTLD_unweighted_df$disease_cat <- "SFTLD"
        BA4_C9FTLD_unweighted_df$disease_cat <- "C9FTLD"

        BA4_unweighted <- rbind(BA4_SALS_unweighted_df, BA4_C9ALS_unweighted_df, BA4_SFTLD_unweighted_df, BA4_C9FTLD_unweighted_df)
        BA4_unweighted <- subset(BA4_unweighted, feature %in% gene_list$gene)

        counts_df <- data.frame(table(BA4_unweighted$celltype, BA4_unweighted$disease_cat ))

        counts_df$Var2 <- factor(counts_df$Var2, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD"))
        counts_df$Var1 <- factor(counts_df$Var1, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))


        ggplot(counts_df, aes(x = Var2, y = Freq, fill = Var2)) + 
        theme_bw() + 
        geom_bar(stat = 'identity') +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~Var1, scales = "free_x") +
        ylab("N genes\n(Z-score > 1)") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) +
        ggtitle('BA4: n genes with LIME Z-score > 1 intersecting with NDKP ALS common (Disease-specific models; unweighted LIME outputs)')
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 2, width = 17)
    ##

    ## BA9 unweighted
        ## NDKP ALS Common -- gene
        gene_list <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/gene_lists/NDKP_common.csv'), sep = ",")
        nrow(gene_list)

        BA9_SALS_unweighted_df$disease_cat <- "SALS"
        BA9_C9ALS_unweighted_df$disease_cat <- "C9ALS"
        BA9_SFTLD_unweighted_df$disease_cat <- "SFTLD"
        BA9_C9FTLD_unweighted_df$disease_cat <- "C9FTLD"

        BA9_unweighted <- rbind(BA9_SALS_unweighted_df, BA9_C9ALS_unweighted_df, BA9_SFTLD_unweighted_df, BA9_C9FTLD_unweighted_df)
        BA9_unweighted <- subset(BA9_unweighted, feature %in% gene_list$gene)

        counts_df <- data.frame(table(BA9_unweighted$celltype, BA9_unweighted$disease_cat ))

        counts_df$Var2 <- factor(counts_df$Var2, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD"))
        counts_df$Var1 <- factor(counts_df$Var1, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))


        ggplot(counts_df, aes(x = Var2, y = Freq, fill = Var2)) + 
        theme_bw() + 
        geom_bar(stat = 'identity') +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~Var1, scales = "free_x") +
        ylab("N genes\n(Z-score > 1)") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) +
        ggtitle('BA9: n genes with LIME Z-score > 1 intersecting with NDKP ALS common (Disease-specific models; unweighted LIME outputs)')
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 2, width = 17)
    ##

    ## BA4 weighted
        ## NDKP ALS Common -- gene
        gene_list <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/gene_lists/NDKP_common.csv'), sep = ",")
        nrow(gene_list)

        BA4_SALS_weighted_df$disease_cat <- "SALS"
        BA4_C9ALS_weighted_df$disease_cat <- "C9ALS"
        BA4_SFTLD_weighted_df$disease_cat <- "SFTLD"
        BA4_C9FTLD_weighted_df$disease_cat <- "C9FTLD"

        BA4_weighted <- rbind(BA4_SALS_weighted_df, BA4_C9ALS_weighted_df, BA4_SFTLD_weighted_df, BA4_C9FTLD_weighted_df)
        BA4_weighted <- subset(BA4_weighted, feature %in% gene_list$gene)

        counts_df <- data.frame(table(BA4_weighted$celltype, BA4_weighted$disease_cat ))

        counts_df$Var2 <- factor(counts_df$Var2, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD"))
        counts_df$Var1 <- factor(counts_df$Var1, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))


        ggplot(counts_df, aes(x = Var2, y = Freq, fill = Var2)) + 
        theme_bw() + 
        geom_bar(stat = 'identity') +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~Var1, scales = "free_x") +
        ylab("N genes\n(Z-score > 1)") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) +
        ggtitle('BA4: n genes with LIME Z-score > 1 intersecting with NDKP ALS common (Disease-specific models; weighted LIME outputs)')
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 2, width = 17)
    ##

    ## BA9 weighted
        ## NDKP ALS Common -- gene
        gene_list <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/gene_lists/NDKP_common.csv'), sep = ",")
        nrow(gene_list)

        BA9_SALS_weighted_df$disease_cat <- "SALS"
        BA9_C9ALS_weighted_df$disease_cat <- "C9ALS"
        BA9_SFTLD_weighted_df$disease_cat <- "SFTLD"
        BA9_C9FTLD_weighted_df$disease_cat <- "C9FTLD"

        BA9_weighted <- rbind(BA9_SALS_weighted_df, BA9_C9ALS_weighted_df, BA9_SFTLD_weighted_df, BA9_C9FTLD_weighted_df)
        BA9_weighted <- subset(BA9_weighted, feature %in% gene_list$gene)

        counts_df <- data.frame(table(BA9_weighted$celltype, BA9_weighted$disease_cat ))

        counts_df$Var2 <- factor(counts_df$Var2, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD"))
        counts_df$Var1 <- factor(counts_df$Var1, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))


        ggplot(counts_df, aes(x = Var2, y = Freq, fill = Var2)) + 
        theme_bw() + 
        geom_bar(stat = 'identity') +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~Var1, scales = "free_x") +
        ylab("N genes\n(Z-score > 1)") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) +
        ggtitle('BA9: n genes with LIME Z-score > 1 intersecting with NDKP ALS common (Disease-specific models; weighted LIME outputs)')
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 2, width = 17)
    ##

    ########################
    ## Plot number of genes above Z > 1 overlaping with NDKP rare 
    ########################

    ## BA4 unweighted
        ## NDKP ALS rare -- gene
        gene_list <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/gene_lists/NDKP_rare.csv'), sep = ",")
        nrow(gene_list)

        BA4_SALS_unweighted_df$disease_cat <- "SALS"
        BA4_C9ALS_unweighted_df$disease_cat <- "C9ALS"
        BA4_SFTLD_unweighted_df$disease_cat <- "SFTLD"
        BA4_C9FTLD_unweighted_df$disease_cat <- "C9FTLD"

        BA4_unweighted <- rbind(BA4_SALS_unweighted_df, BA4_C9ALS_unweighted_df, BA4_SFTLD_unweighted_df, BA4_C9FTLD_unweighted_df)
        BA4_unweighted <- subset(BA4_unweighted, feature %in% gene_list$gene)

        counts_df <- data.frame(table(BA4_unweighted$celltype, BA4_unweighted$disease_cat ))

        counts_df$Var2 <- factor(counts_df$Var2, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD"))
        counts_df$Var1 <- factor(counts_df$Var1, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))


        ggplot(counts_df, aes(x = Var2, y = Freq, fill = Var2)) + 
        theme_bw() + 
        geom_bar(stat = 'identity') +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~Var1, scales = "free_x") +
        ylab("N genes\n(Z-score > 1)") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) +
        ggtitle('BA4: n genes with LIME Z-score > 1 intersecting with NDKP ALS rare (Disease-specific models; unweighted LIME outputs)')
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 2, width = 17)
    ##

    ## BA9 unweighted
        ## NDKP ALS rare -- gene
        gene_list <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/gene_lists/NDKP_rare.csv'), sep = ",")
        nrow(gene_list)

        BA9_SALS_unweighted_df$disease_cat <- "SALS"
        BA9_C9ALS_unweighted_df$disease_cat <- "C9ALS"
        BA9_SFTLD_unweighted_df$disease_cat <- "SFTLD"
        BA9_C9FTLD_unweighted_df$disease_cat <- "C9FTLD"

        BA9_unweighted <- rbind(BA9_SALS_unweighted_df, BA9_C9ALS_unweighted_df, BA9_SFTLD_unweighted_df, BA9_C9FTLD_unweighted_df)
        BA9_unweighted <- subset(BA9_unweighted, feature %in% gene_list$gene)

        counts_df <- data.frame(table(BA9_unweighted$celltype, BA9_unweighted$disease_cat ))

        counts_df$Var2 <- factor(counts_df$Var2, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD"))
        counts_df$Var1 <- factor(counts_df$Var1, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))


        ggplot(counts_df, aes(x = Var2, y = Freq, fill = Var2)) + 
        theme_bw() + 
        geom_bar(stat = 'identity') +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~Var1, scales = "free_x") +
        ylab("N genes\n(Z-score > 1)") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) +
        ggtitle('BA9: n genes with LIME Z-score > 1 intersecting with NDKP ALS rare (Disease-specific models; unweighted LIME outputs)')
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 2, width = 17)
    ##

    ## BA4 weighted
        ## NDKP ALS rare -- gene
        gene_list <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/gene_lists/NDKP_rare.csv'), sep = ",")
        nrow(gene_list)

        BA4_SALS_weighted_df$disease_cat <- "SALS"
        BA4_C9ALS_weighted_df$disease_cat <- "C9ALS"
        BA4_SFTLD_weighted_df$disease_cat <- "SFTLD"
        BA4_C9FTLD_weighted_df$disease_cat <- "C9FTLD"

        BA4_weighted <- rbind(BA4_SALS_weighted_df, BA4_C9ALS_weighted_df, BA4_SFTLD_weighted_df, BA4_C9FTLD_weighted_df)
        BA4_weighted <- subset(BA4_weighted, feature %in% gene_list$gene)

        counts_df <- data.frame(table(BA4_weighted$celltype, BA4_weighted$disease_cat ))

        counts_df$Var2 <- factor(counts_df$Var2, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD"))
        counts_df$Var1 <- factor(counts_df$Var1, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))


        ggplot(counts_df, aes(x = Var2, y = Freq, fill = Var2)) + 
        theme_bw() + 
        geom_bar(stat = 'identity') +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~Var1, scales = "free_x") +
        ylab("N genes\n(Z-score > 1)") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) +
        ggtitle('BA4: n genes with LIME Z-score > 1 intersecting with NDKP ALS rare (Disease-specific models; weighted LIME outputs)')
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 2, width = 17)
    ##

    ## BA9 weighted
        ## NDKP ALS rare -- gene
        gene_list <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/gene_lists/NDKP_rare.csv'), sep = ",")
        nrow(gene_list)

        BA9_SALS_weighted_df$disease_cat <- "SALS"
        BA9_C9ALS_weighted_df$disease_cat <- "C9ALS"
        BA9_SFTLD_weighted_df$disease_cat <- "SFTLD"
        BA9_C9FTLD_weighted_df$disease_cat <- "C9FTLD"

        BA9_weighted <- rbind(BA9_SALS_weighted_df, BA9_C9ALS_weighted_df, BA9_SFTLD_weighted_df, BA9_C9FTLD_weighted_df)
        BA9_weighted <- subset(BA9_weighted, feature %in% gene_list$gene)

        counts_df <- data.frame(table(BA9_weighted$celltype, BA9_weighted$disease_cat ))

        counts_df$Var2 <- factor(counts_df$Var2, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD"))
        counts_df$Var1 <- factor(counts_df$Var1, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))


        ggplot(counts_df, aes(x = Var2, y = Freq, fill = Var2)) + 
        theme_bw() + 
        geom_bar(stat = 'identity') +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~Var1, scales = "free_x") +
        ylab("N genes\n(Z-score > 1)") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) +
        ggtitle('BA9: n genes with LIME Z-score > 1 intersecting with NDKP ALS rare (Disease-specific models; weighted LIME outputs)')
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 2, width = 17)
    ##

    ########################
    ## Plot number of genes above Z > 1 overlaping with Zhang genes
    ########################

    ## BA4 unweighted
        ## Zhang genes -- gene
        gene_list <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/gene_lists/Zhang_neuron.csv'), sep = ",")
        nrow(gene_list)

        BA4_SALS_unweighted_df$disease_cat <- "SALS"
        BA4_C9ALS_unweighted_df$disease_cat <- "C9ALS"
        BA4_SFTLD_unweighted_df$disease_cat <- "SFTLD"
        BA4_C9FTLD_unweighted_df$disease_cat <- "C9FTLD"

        BA4_unweighted <- rbind(BA4_SALS_unweighted_df, BA4_C9ALS_unweighted_df, BA4_SFTLD_unweighted_df, BA4_C9FTLD_unweighted_df)
        BA4_unweighted <- subset(BA4_unweighted, feature %in% gene_list$gene)

        counts_df <- data.frame(table(BA4_unweighted$celltype, BA4_unweighted$disease_cat ))

        counts_df$Var2 <- factor(counts_df$Var2, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD"))
        counts_df$Var1 <- factor(counts_df$Var1, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))


        ggplot(counts_df, aes(x = Var2, y = Freq, fill = Var2)) + 
        theme_bw() + 
        geom_bar(stat = 'identity') +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~Var1, scales = "free_x") +
        ylab("N genes\n(Z-score > 1)") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) +
        ggtitle('BA4: n genes with LIME Z-score > 1 intersecting with Zhang et al. (Disease-specific models; unweighted LIME outputs)')
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 2, width = 17)
    ##

    ## BA9 unweighted
        ## Zhang genes -- gene
        gene_list <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/gene_lists/Zhang_neuron.csv'), sep = ",")
        nrow(gene_list)

        BA9_SALS_unweighted_df$disease_cat <- "SALS"
        BA9_C9ALS_unweighted_df$disease_cat <- "C9ALS"
        BA9_SFTLD_unweighted_df$disease_cat <- "SFTLD"
        BA9_C9FTLD_unweighted_df$disease_cat <- "C9FTLD"

        BA9_unweighted <- rbind(BA9_SALS_unweighted_df, BA9_C9ALS_unweighted_df, BA9_SFTLD_unweighted_df, BA9_C9FTLD_unweighted_df)
        BA9_unweighted <- subset(BA9_unweighted, feature %in% gene_list$gene)

        counts_df <- data.frame(table(BA9_unweighted$celltype, BA9_unweighted$disease_cat ))

        counts_df$Var2 <- factor(counts_df$Var2, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD"))
        counts_df$Var1 <- factor(counts_df$Var1, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))


        ggplot(counts_df, aes(x = Var2, y = Freq, fill = Var2)) + 
        theme_bw() + 
        geom_bar(stat = 'identity') +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~Var1, scales = "free_x") +
        ylab("N genes\n(Z-score > 1)") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) +
        ggtitle('BA9: n genes with LIME Z-score > 1 intersecting with Zhang et al. (Disease-specific models; unweighted LIME outputs)')
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 2, width = 17)
    ##

    ## BA4 weighted
        ## Zhang genes -- gene
        gene_list <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/gene_lists/Zhang_neuron.csv'), sep = ",")
        nrow(gene_list)

        BA4_SALS_weighted_df$disease_cat <- "SALS"
        BA4_C9ALS_weighted_df$disease_cat <- "C9ALS"
        BA4_SFTLD_weighted_df$disease_cat <- "SFTLD"
        BA4_C9FTLD_weighted_df$disease_cat <- "C9FTLD"

        BA4_weighted <- rbind(BA4_SALS_weighted_df, BA4_C9ALS_weighted_df, BA4_SFTLD_weighted_df, BA4_C9FTLD_weighted_df)
        BA4_weighted <- subset(BA4_weighted, feature %in% gene_list$gene)

        counts_df <- data.frame(table(BA4_weighted$celltype, BA4_weighted$disease_cat ))

        counts_df$Var2 <- factor(counts_df$Var2, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD"))
        counts_df$Var1 <- factor(counts_df$Var1, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))


        ggplot(counts_df, aes(x = Var2, y = Freq, fill = Var2)) + 
        theme_bw() + 
        geom_bar(stat = 'identity') +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~Var1, scales = "free_x") +
        ylab("N genes\n(Z-score > 1)") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) +
        ggtitle('BA4: n genes with LIME Z-score > 1 intersecting with Zhang et al. (Disease-specific models; weighted LIME outputs)')
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 2, width = 17)
    ##

    ## BA9 weighted
        ## Zhang genes -- gene
        gene_list <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/gene_lists/Zhang_neuron.csv'), sep = ",")
        nrow(gene_list)

        BA9_SALS_weighted_df$disease_cat <- "SALS"
        BA9_C9ALS_weighted_df$disease_cat <- "C9ALS"
        BA9_SFTLD_weighted_df$disease_cat <- "SFTLD"
        BA9_C9FTLD_weighted_df$disease_cat <- "C9FTLD"

        BA9_weighted <- rbind(BA9_SALS_weighted_df, BA9_C9ALS_weighted_df, BA9_SFTLD_weighted_df, BA9_C9FTLD_weighted_df)
        BA9_weighted <- subset(BA9_weighted, feature %in% gene_list$gene)

        counts_df <- data.frame(table(BA9_weighted$celltype, BA9_weighted$disease_cat ))

        counts_df$Var2 <- factor(counts_df$Var2, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD"))
        counts_df$Var1 <- factor(counts_df$Var1, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))


        ggplot(counts_df, aes(x = Var2, y = Freq, fill = Var2)) + 
        theme_bw() + 
        geom_bar(stat = 'identity') +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~Var1, scales = "free_x") +
        ylab("N genes\n(Z-score > 1)") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) +
        ggtitle('BA9: n genes with LIME Z-score > 1 intersecting with Zhang et al. (Disease-specific models; weighted LIME outputs)')
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 2, width = 17)
    ##


    ########################
    ## Plot number of genes above Z > 1 overlaping with ND genes from UK program
    ########################

    ## BA4 unweighted
        ## ND genes -- gene
        gene_list <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/gene_lists/ND_genes_July_29.tsv'), sep = "\t")
        nrow(gene_list)

        BA4_SALS_unweighted_df$disease_cat <- "SALS"
        BA4_C9ALS_unweighted_df$disease_cat <- "C9ALS"
        BA4_SFTLD_unweighted_df$disease_cat <- "SFTLD"
        BA4_C9FTLD_unweighted_df$disease_cat <- "C9FTLD"

        BA4_unweighted <- rbind(BA4_SALS_unweighted_df, BA4_C9ALS_unweighted_df, BA4_SFTLD_unweighted_df, BA4_C9FTLD_unweighted_df)
        BA4_unweighted <- subset(BA4_unweighted, feature %in% gene_list$Gene.Symbol)

        counts_df <- data.frame(table(BA4_unweighted$celltype, BA4_unweighted$disease_cat ))

        counts_df$Var2 <- factor(counts_df$Var2, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD"))
        counts_df$Var1 <- factor(counts_df$Var1, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))


        ggplot(counts_df, aes(x = Var2, y = Freq, fill = Var2)) + 
        theme_bw() + 
        geom_bar(stat = 'identity') +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~Var1, scales = "free_x") +
        ylab("N genes\n(Z-score > 1)") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) +
        ggtitle('BA4: n genes with LIME Z-score > 1 intersecting with Neurodegenerative Disease genes (Disease-specific models; unweighted LIME outputs)')
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 2, width = 17)
    ##

    ## BA9 unweighted
        ## ND genes -- gene
        gene_list <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/gene_lists/ND_genes_July_29.tsv'), sep = "\t")
        nrow(gene_list)

        BA9_SALS_unweighted_df$disease_cat <- "SALS"
        BA9_C9ALS_unweighted_df$disease_cat <- "C9ALS"
        BA9_SFTLD_unweighted_df$disease_cat <- "SFTLD"
        BA9_C9FTLD_unweighted_df$disease_cat <- "C9FTLD"

        BA9_unweighted <- rbind(BA9_SALS_unweighted_df, BA9_C9ALS_unweighted_df, BA9_SFTLD_unweighted_df, BA9_C9FTLD_unweighted_df)
        BA9_unweighted <- subset(BA9_unweighted, feature %in% gene_list$Gene.Symbol)

        counts_df <- data.frame(table(BA9_unweighted$celltype, BA9_unweighted$disease_cat ))

        counts_df$Var2 <- factor(counts_df$Var2, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD"))
        counts_df$Var1 <- factor(counts_df$Var1, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))


        ggplot(counts_df, aes(x = Var2, y = Freq, fill = Var2)) + 
        theme_bw() + 
        geom_bar(stat = 'identity') +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~Var1, scales = "free_x") +
        ylab("N genes\n(Z-score > 1)") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) +
        ggtitle('BA9: n genes with LIME Z-score > 1 intersecting with Neurodegenerative Disease genes (Disease-specific models; unweighted LIME outputs)')
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 2, width = 17)
    ##

    ## BA4 weighted
        ## ND genes -- gene
        gene_list <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/gene_lists/ND_genes_July_29.tsv'), sep = "\t")
        nrow(gene_list)

        BA4_SALS_weighted_df$disease_cat <- "SALS"
        BA4_C9ALS_weighted_df$disease_cat <- "C9ALS"
        BA4_SFTLD_weighted_df$disease_cat <- "SFTLD"
        BA4_C9FTLD_weighted_df$disease_cat <- "C9FTLD"

        BA4_weighted <- rbind(BA4_SALS_weighted_df, BA4_C9ALS_weighted_df, BA4_SFTLD_weighted_df, BA4_C9FTLD_weighted_df)
        BA4_weighted <- subset(BA4_weighted, feature %in% gene_list$Gene.Symbol)

        counts_df <- data.frame(table(BA4_weighted$celltype, BA4_weighted$disease_cat ))

        counts_df$Var2 <- factor(counts_df$Var2, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD"))
        counts_df$Var1 <- factor(counts_df$Var1, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))


        ggplot(counts_df, aes(x = Var2, y = Freq, fill = Var2)) + 
        theme_bw() + 
        geom_bar(stat = 'identity') +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~Var1, scales = "free_x") +
        ylab("N genes\n(Z-score > 1)") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) +
        ggtitle('BA4: n genes with LIME Z-score > 1 intersecting with Neurodegenerative Disease genes (Disease-specific models; weighted LIME outputs)')
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 2, width = 17)
    ##

    ## BA9 weighted
        ## ND genes -- gene
        gene_list <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/gene_lists/ND_genes_July_29.tsv'), sep = "\t")
        nrow(gene_list)

        BA9_SALS_weighted_df$disease_cat <- "SALS"
        BA9_C9ALS_weighted_df$disease_cat <- "C9ALS"
        BA9_SFTLD_weighted_df$disease_cat <- "SFTLD"
        BA9_C9FTLD_weighted_df$disease_cat <- "C9FTLD"

        BA9_weighted <- rbind(BA9_SALS_weighted_df, BA9_C9ALS_weighted_df, BA9_SFTLD_weighted_df, BA9_C9FTLD_weighted_df)
        BA9_weighted <- subset(BA9_weighted, feature %in% gene_list$Gene.Symbol)

        counts_df <- data.frame(table(BA9_weighted$celltype, BA9_weighted$disease_cat ))

        counts_df$Var2 <- factor(counts_df$Var2, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD"))
        counts_df$Var1 <- factor(counts_df$Var1, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))


        ggplot(counts_df, aes(x = Var2, y = Freq, fill = Var2)) + 
        theme_bw() + 
        geom_bar(stat = 'identity') +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~Var1, scales = "free_x") +
        ylab("N genes\n(Z-score > 1)") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) +
        ggtitle('BA9: n genes with LIME Z-score > 1 intersecting with Neurodegenerative Disease genes (Disease-specific models; weighted LIME outputs)')
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 2, width = 17)
    ##

    ########################
    ## Plot number of genes above Z > 1 overlaping with ALS ClinGen genes
    ########################

    ## BA4 unweighted
        ## ClinGen genes
        gene_list <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/gene_lists/ClinGen_ALS_Dec_30_2024.csv'), sep = ",")
        gene_list$Gene <- gsub("HGNC.*", "", gene_list$Gene)
        nrow(gene_list)

        BA4_SALS_unweighted_df$disease_cat <- "SALS"
        BA4_C9ALS_unweighted_df$disease_cat <- "C9ALS"
        BA4_SFTLD_unweighted_df$disease_cat <- "SFTLD"
        BA4_C9FTLD_unweighted_df$disease_cat <- "C9FTLD"

        BA4_unweighted <- rbind(BA4_SALS_unweighted_df, BA4_C9ALS_unweighted_df, BA4_SFTLD_unweighted_df, BA4_C9FTLD_unweighted_df)
        BA4_unweighted <- subset(BA4_unweighted, feature %in% gene_list$Gene)

        counts_df <- data.frame(table(BA4_unweighted$celltype, BA4_unweighted$disease_cat ))

        counts_df$Var2 <- factor(counts_df$Var2, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD"))
        counts_df$Var1 <- factor(counts_df$Var1, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))


        ggplot(counts_df, aes(x = Var2, y = Freq, fill = Var2)) + 
        theme_bw() + 
        geom_bar(stat = 'identity') +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~Var1, scales = "free_x") +
        ylab("N genes\n(Z-score > 1)") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) +
        ggtitle('BA4: n genes with LIME Z-score > 1 intersecting with ALS ClinGen genes (Disease-specific models; unweighted LIME outputs)')
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 2, width = 17)
    ##

    ## BA9 unweighted
        ## ClinGen genes
        gene_list <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/gene_lists/ClinGen_ALS_Dec_30_2024.csv'), sep = ",")
        gene_list$Gene <- gsub("HGNC.*", "", gene_list$Gene)
        nrow(gene_list)

        BA9_SALS_unweighted_df$disease_cat <- "SALS"
        BA9_C9ALS_unweighted_df$disease_cat <- "C9ALS"
        BA9_SFTLD_unweighted_df$disease_cat <- "SFTLD"
        BA9_C9FTLD_unweighted_df$disease_cat <- "C9FTLD"

        BA9_unweighted <- rbind(BA9_SALS_unweighted_df, BA9_C9ALS_unweighted_df, BA9_SFTLD_unweighted_df, BA9_C9FTLD_unweighted_df)
        BA9_unweighted <- subset(BA9_unweighted, feature %in% gene_list$Gene)

        counts_df <- data.frame(table(BA9_unweighted$celltype, BA9_unweighted$disease_cat ))

        counts_df$Var2 <- factor(counts_df$Var2, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD"))
        counts_df$Var1 <- factor(counts_df$Var1, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))


        ggplot(counts_df, aes(x = Var2, y = Freq, fill = Var2)) + 
        theme_bw() + 
        geom_bar(stat = 'identity') +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~Var1, scales = "free_x") +
        ylab("N genes\n(Z-score > 1)") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) +
        ggtitle('BA9: n genes with LIME Z-score > 1 intersecting with ALS ClinGen genes (Disease-specific models; unweighted LIME outputs)')
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 2, width = 17)
    ##

    ## BA4 weighted
        ## ClinGen genes
        gene_list <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/gene_lists/ClinGen_ALS_Dec_30_2024.csv'), sep = ",")
        gene_list$Gene <- gsub("HGNC.*", "", gene_list$Gene)
        nrow(gene_list)

        BA4_SALS_weighted_df$disease_cat <- "SALS"
        BA4_C9ALS_weighted_df$disease_cat <- "C9ALS"
        BA4_SFTLD_weighted_df$disease_cat <- "SFTLD"
        BA4_C9FTLD_weighted_df$disease_cat <- "C9FTLD"

        BA4_weighted <- rbind(BA4_SALS_weighted_df, BA4_C9ALS_weighted_df, BA4_SFTLD_weighted_df, BA4_C9FTLD_weighted_df)
        BA4_weighted <- subset(BA4_weighted, feature %in% gene_list$Gene)

        counts_df <- data.frame(table(BA4_weighted$celltype, BA4_weighted$disease_cat ))

        counts_df$Var2 <- factor(counts_df$Var2, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD"))
        counts_df$Var1 <- factor(counts_df$Var1, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))


        ggplot(counts_df, aes(x = Var2, y = Freq, fill = Var2)) + 
        theme_bw() + 
        geom_bar(stat = 'identity') +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~Var1, scales = "free_x") +
        ylab("N genes\n(Z-score > 1)") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) +
        ggtitle('BA4: n genes with LIME Z-score > 1 intersecting with ALS ClinGen genes (Disease-specific models; weighted LIME outputs)')
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 2, width = 17)
    ##

    ## BA9 weighted
        ## ClinGen genes
        gene_list <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/gene_lists/ClinGen_ALS_Dec_30_2024.csv'), sep = ",")
        gene_list$Gene <- gsub("HGNC.*", "", gene_list$Gene)
        nrow(gene_list)

        BA9_SALS_weighted_df$disease_cat <- "SALS"
        BA9_C9ALS_weighted_df$disease_cat <- "C9ALS"
        BA9_SFTLD_weighted_df$disease_cat <- "SFTLD"
        BA9_C9FTLD_weighted_df$disease_cat <- "C9FTLD"

        BA9_weighted <- rbind(BA9_SALS_weighted_df, BA9_C9ALS_weighted_df, BA9_SFTLD_weighted_df, BA9_C9FTLD_weighted_df)
        BA9_weighted <- subset(BA9_weighted, feature %in% gene_list$Gene)

        counts_df <- data.frame(table(BA9_weighted$celltype, BA9_weighted$disease_cat ))

        counts_df$Var2 <- factor(counts_df$Var2, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD"))
        counts_df$Var1 <- factor(counts_df$Var1, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))


        ggplot(counts_df, aes(x = Var2, y = Freq, fill = Var2)) + 
        theme_bw() + 
        geom_bar(stat = 'identity') +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~Var1, scales = "free_x") +
        ylab("N genes\n(Z-score > 1)") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) +
        ggtitle('BA9: n genes with LIME Z-score > 1 intersecting with ALS ClinGen genes (Disease-specific models; weighted LIME outputs)')
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 2, width = 17)
    ##
##

## MAKE A Z-score > 1 overlap heatmap







############################################################################ Explore genes with LIME Z score > 1
############################################################################
############################################################################
############################################################################ Condition specific models














########################
## Check intersects with gene lists
########################

## code 
    ##############################
    ## read in gene lists
    ##############################


    ## NDKP ALS Common -- gene
    NDKP_ALS_common <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/gene_lists/NDKP_common.csv'), sep = ",")

    ## NDKP ALS Rare -- gene
    NDKP_ALS_rare <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/gene_lists/NDKP_rare.csv'), sep = ",")

    ## Zhang Neuron paper genes -- gene
    zhang_neuron <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/gene_lists/Zhang_neuron.csv'), sep = ",")

    ## ND Genes from UK panel -- Gene.Symbol
    ND_genes <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/gene_lists/ND_genes_July_29.tsv'), sep = "\t")

    ## ClinGen
    ClinGen_genes <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/gene_lists/ClinGen_ALS_Dec_30_2024.csv'), sep = ",")
    ClinGen_genes$Gene <- gsub("HGNC.*", "", ClinGen_genes$Gene)

    ##############################
    ## read in gene lists
    ##############################



########################
## All conditions together DNN model LIME outputs
########################







