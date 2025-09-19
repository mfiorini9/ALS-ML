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


## RUN THIS CODE AFTER THE PERMUTATIONS.

############################################################################ N features with Z-score > 1 -- weighted by percent expression (Supplemental Figure)
############################################################################
############################################################################
############################################################################ 
## code
    ###########################
    ## SALS BA4
    ###########################
    ## code
        #celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
        #        "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
        #        "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")
        
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
        "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
        "Micro", "Mural",   "Endo",    "Fibro", "L5")

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
        #celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
        #        "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
        #        "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")
        
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
        "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
        "Micro", "Mural",   "Endo",    "Fibro", "L5")

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
        #celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
        #        "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
        #        "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")
        
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
        "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
        "Micro", "Mural",   "Endo",    "Fibro", "L5")

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
        #celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
        #        "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
        #        "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")
        
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
        "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
        "Micro", "Mural",   "Endo",    "Fibro", "L5")

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
        #celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
        #        "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
        #        "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")
        
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
        "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
        "Micro", "Mural",   "Endo",    "Fibro", "L5")

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
        #celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
        #        "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
        #        "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")
        
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
        "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
        "Micro", "Mural",   "Endo",    "Fibro", "L5")

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
        #celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
        #        "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
        #        "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")
        
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
        "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
        "Micro", "Mural",   "Endo",    "Fibro", "L5")

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
    ## C9FTLD BA9 
    ###########################
    ## code 

        #celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
        #        "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
        #        "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")
        
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
        "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
        "Micro", "Mural",   "Endo",    "Fibro", "L5")

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
        #celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
        #        "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
        #        "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")
        
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
        "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
        "Micro", "Mural",   "Endo",    "Fibro", "L5")

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
        #celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
        #        "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
        #        "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")
        
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
        "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
        "Micro", "Mural",   "Endo",    "Fibro", "L5")

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
        BA4_all_conditions_weighted_df$disease_cat <- "All conditions"

        BA4_weighted_counts <- rbind(BA4_SALS_weighted_df, BA4_C9ALS_weighted_df, BA4_SFTLD_weighted_df, BA4_C9FTLD_weighted_df, BA4_all_conditions_weighted_df)

        counts_df <- data.frame(table(BA4_weighted_counts$celltype, BA4_weighted_counts$disease_cat ))
        
        counts_df$Var2 <- factor(counts_df$Var2, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All conditions"))
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
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA4_nGenes_Zscore1.pdf', height = 2, width = 13)
    ##

    ###########################
    ## Plot BA9
    ###########################
    ## code
        BA9_SALS_weighted_df$disease_cat <- "SALS"
        BA9_C9ALS_weighted_df$disease_cat <- "C9ALS"
        BA9_SFTLD_weighted_df$disease_cat <- "SFTLD"
        BA9_C9FTLD_weighted_df$disease_cat <- "C9FTLD"
        BA9_all_conditions_weighted_df$disease_cat <- "All conditions"

        BA9_weighted_counts <- rbind(BA9_SALS_weighted_df, BA9_C9ALS_weighted_df, BA9_SFTLD_weighted_df, BA9_C9FTLD_weighted_df, BA9_all_conditions_weighted_df)

        counts_df <- data.frame(table(BA9_weighted_counts$celltype, BA9_weighted_counts$disease_cat ))
        
        counts_df$Var2 <- factor(counts_df$Var2, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All conditions"))
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
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA9_nGenes_Zscore1.pdf', height = 2, width = 13)
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
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        #celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",  "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")

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
        
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/permutation_weighted_report_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval_2.csv'), sep = ",")
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
        #celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        #celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",  "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")


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
        
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/permutation_weighted_report_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval_2.csv'), sep = ",")
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
        #celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro", "Mural",   "Endo",    "Fibro", "L5")

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
        
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/permutation_weighted_report_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval_2.csv'), sep = ",")
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
        
        #celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro", "Mural",   "Endo",    "Fibro", "L5")

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
        
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/permutation_weighted_report_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval_2.csv'), sep = ",")
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
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro", "Mural",   "Endo",    "Fibro", "L5")

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
        
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/permutation_weighted_report_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval_2.csv'), sep = ",")
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
        #celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro", "Mural",   "Endo",    "Fibro", "L5")

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
        
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/permutation_weighted_report_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval_2.csv'), sep = ",")
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
        #celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro", "Mural",   "Endo",    "Fibro", "L5")

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
        
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/permutation_weighted_report_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval_2.csv'), sep = ",")
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
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")

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
        
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/permutation_weighted_report_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval_2.csv'), sep = ",")
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
    ## SALS BA4 -- Lambda approach and plot -- DONE
    ###########################
    ## code
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        #celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",  "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        
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

        
        BA4_SALS_weighted_df_summary <- BA4_SALS_weighted_df %>%
            dplyr::group_by(group, celltype, n_genes, method) %>%
            dplyr::summarise(mean_accuracy_by_group = mean(mean_accuracy, na.rm = TRUE))

        BA4_SALS_weighted_df_summary <- subset(BA4_SALS_weighted_df_summary, n_genes >= 25)

        
        result <- subset(result, lambda == 0.00025)
        result <- result %>% dplyr::select(celltype, n_genes)
        
        colnames(result) <- c("celltype", "optimal_genes")
        
        BA4_SALS_weighted_df_summary <- merge(BA4_SALS_weighted_df_summary, result, by = "celltype")

        ## modify so that if the threshold value is below we take the number of genes corresponding to 80. 
        for(celltype2 in celltype_list){
        if (BA4_SALS_weighted_df_summary$mean_accuracy_by_group[BA4_SALS_weighted_df_summary$celltype == celltype2 &  BA4_SALS_weighted_df_summary$method == "LIME" & BA4_SALS_weighted_df_summary$n_genes == BA4_SALS_weighted_df_summary$optimal_genes] < 0.80 & max(BA4_SALS_weighted_df_summary$mean_accuracy_by_group[BA4_SALS_weighted_df_summary$celltype == celltype2 &  BA4_SALS_weighted_df_summary$method == "LIME"]) >= 0.8 ) {
    
            # Find the smallest n_genes where mean_accuracy_by_group >= 0.80
            new_optimal_genes <- min(BA4_SALS_weighted_df_summary$n_genes[BA4_SALS_weighted_df_summary$mean_accuracy_by_group >= 0.80 & BA4_SALS_weighted_df_summary$celltype == celltype2 & BA4_SALS_weighted_df_summary$method == "LIME"])
            
            # Update the optimal_genes column with the new value
            BA4_SALS_weighted_df_summary$optimal_genes[BA4_SALS_weighted_df_summary$celltype == celltype2] <- new_optimal_genes
            }
        }

        subset(BA4_SALS_weighted_df_summary, celltype == "L5")
        
        ## custom label
        BA4_SALS_weighted_df_summary <- BA4_SALS_weighted_df_summary %>%
            mutate(label = paste(celltype, "\n(", optimal_genes, " genes)", sep = ""))

        BA4_SALS_weighted_df_summary$method <- factor(BA4_SALS_weighted_df_summary$method, levels = c("LIME", "random_HVG", "random_gene"))
        BA4_SALS_weighted_df_summary$celltype <- factor(BA4_SALS_weighted_df_summary$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))
        
        BA4_SALS_weighted_df_summary <- BA4_SALS_weighted_df_summary %>%
            arrange(celltype)  # Arrange by the 'celltype' factor

        BA4_SALS_weighted_df_summary$label <- factor(BA4_SALS_weighted_df_summary$label, levels = unique(BA4_SALS_weighted_df_summary$label))
        
    
        ## prin the number of genes corresponding to the optimal for each cell type
        result <- BA4_SALS_weighted_df_summary %>% dplyr::select(celltype, optimal_genes)
        result <- result %>% distinct()
        
        write.csv(result, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_lambda_0.00025_gene_set_BA4_SALS_narval_2.csv'))

        ## Plot
        ggplot(BA4_SALS_weighted_df_summary, aes(x = n_genes, y = mean_accuracy_by_group, colour = method, group = method)) + 
        theme_bw() + 
        geom_vline(aes(xintercept = optimal_genes), linetype = "dashed", color = "black") +
        geom_point(size = 1) +
        geom_line() +
        geom_hline(yintercept = 0.8, colour = "red") +
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
        ylab("Balanced accuracy") +
        xlab("N genes") +
        scale_y_continuous(limits = c(0.5, 1)) +
        scale_colour_manual(values = c("#F18F01", "#2E4057"))
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA4_SALS_permutation_dotplot.pdf', height = 2, width = 15)
    ##



    ###########################
    ## C9ALS BA4 -- Lambda approach and plot -- DONE
    ###########################
    ## code
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        #celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",  "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
    
        ## sumarize C9ALS
        BA4_C9ALS_weighted_df$disease_cat <- "C9ALS"
        
        BA4_C9ALS_weighted_df_summary <- BA4_C9ALS_weighted_df %>%
            dplyr::group_by(group, celltype, n_genes, method) %>%
            dplyr::summarise(mean_accuracy_by_group = mean(mean_accuracy, na.rm = TRUE))
        
        BA4_C9ALS_weighted_df_summary <- data.frame(BA4_C9ALS_weighted_df_summary)

        BA4_C9ALS_weighted_df_summary <- subset(BA4_C9ALS_weighted_df_summary, method == "LIME" )
        
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
        
        BA4_C9ALS_weighted_df_summary <- BA4_C9ALS_weighted_df %>%
            dplyr::group_by(group, celltype, n_genes, method) %>%
            dplyr::summarise(mean_accuracy_by_group = mean(mean_accuracy, na.rm = TRUE))

        BA4_C9ALS_weighted_df_summary <- subset(BA4_C9ALS_weighted_df_summary, n_genes >= 25)

        
        result <- subset(result, lambda == 0.00025)
        result <- result %>% dplyr::select(celltype, n_genes)
        
        colnames(result) <- c("celltype", "optimal_genes")
        
        BA4_C9ALS_weighted_df_summary <- merge(BA4_C9ALS_weighted_df_summary, result, by = "celltype")

        ## modify so that if the threshold value is below we take the number of genes corresponding to 80. 
        for(celltype2 in celltype_list){
        if (BA4_C9ALS_weighted_df_summary$mean_accuracy_by_group[BA4_C9ALS_weighted_df_summary$celltype == celltype2 &  BA4_C9ALS_weighted_df_summary$method == "LIME" & BA4_C9ALS_weighted_df_summary$n_genes == BA4_C9ALS_weighted_df_summary$optimal_genes] < 0.80 & max(BA4_C9ALS_weighted_df_summary$mean_accuracy_by_group[BA4_C9ALS_weighted_df_summary$celltype == celltype2 &  BA4_C9ALS_weighted_df_summary$method == "LIME"]) >= 0.8 ) {
    
            # Find the smallest n_genes where mean_accuracy_by_group >= 0.80
            new_optimal_genes <- min(BA4_C9ALS_weighted_df_summary$n_genes[BA4_C9ALS_weighted_df_summary$mean_accuracy_by_group >= 0.80 & BA4_C9ALS_weighted_df_summary$celltype == celltype2 & BA4_C9ALS_weighted_df_summary$method == "LIME"])
            
            # Update the optimal_genes column with the new value
            BA4_C9ALS_weighted_df_summary$optimal_genes[BA4_C9ALS_weighted_df_summary$celltype == celltype2] <- new_optimal_genes
            }
        }

        ## custom label
        BA4_C9ALS_weighted_df_summary <- BA4_C9ALS_weighted_df_summary %>%
            mutate(label = paste(celltype, "\n(", optimal_genes, " genes)", sep = ""))

        BA4_C9ALS_weighted_df_summary$method <- factor(BA4_C9ALS_weighted_df_summary$method, levels = c("LIME", "random_HVG", "random_gene"))
        BA4_C9ALS_weighted_df_summary$celltype <- factor(BA4_C9ALS_weighted_df_summary$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))
        
        BA4_C9ALS_weighted_df_summary <- BA4_C9ALS_weighted_df_summary %>%
            arrange(celltype)  # Arrange by the 'celltype' factor

        BA4_C9ALS_weighted_df_summary$label <- factor(BA4_C9ALS_weighted_df_summary$label, levels = unique(BA4_C9ALS_weighted_df_summary$label))
        
    
        ## prin the number of genes corresponding to the optimal for each cell type
        result <- BA4_C9ALS_weighted_df_summary %>% dplyr::select(celltype, optimal_genes)
        result <- result %>% distinct()
        
        write.csv(result, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_lambda_0.00025_gene_set_BA4_C9ALS_narval_2.csv'))

        ## Plot
        ggplot(BA4_C9ALS_weighted_df_summary, aes(x = n_genes, y = mean_accuracy_by_group, colour = method, group = method)) + 
        theme_bw() + 
        geom_vline(aes(xintercept = optimal_genes), linetype = "dashed", color = "black") +
        geom_point(size = 1) +
        geom_line() +
        geom_hline(yintercept = 0.8, colour = "red") +
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
        ylab("Balanced accuracy") +
        xlab("N genes") +
        scale_y_continuous(limits = c(0.5, 1)) +
        scale_colour_manual(values = c("#F18F01", "#2E4057"))
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA4_C9ALS_permutation_dotplot.pdf', height = 2, width = 15)
    ##
    
    ###########################
    ## SFTLD BA4 -- Lambda approach and plot -- DONE
    ###########################
    ## code
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
    
        ## sumarize SFTLD
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
        
        BA4_SFTLD_weighted_df_summary <- BA4_SFTLD_weighted_df %>%
            dplyr::group_by(group, celltype, n_genes, method) %>%
            dplyr::summarise(mean_accuracy_by_group = mean(mean_accuracy, na.rm = TRUE))

        BA4_SFTLD_weighted_df_summary <- subset(BA4_SFTLD_weighted_df_summary, n_genes >= 25)

        
        result <- subset(result, lambda == 0.00025)
        result <- result %>% dplyr::select(celltype, n_genes)
        
        colnames(result) <- c("celltype", "optimal_genes")
        
        BA4_SFTLD_weighted_df_summary <- merge(BA4_SFTLD_weighted_df_summary, result, by = "celltype")

        ## modify so that if the threshold value is below we take the number of genes corresponding to 80. 
        for(celltype2 in celltype_list){
        if (BA4_SFTLD_weighted_df_summary$mean_accuracy_by_group[BA4_SFTLD_weighted_df_summary$celltype == celltype2 &  BA4_SFTLD_weighted_df_summary$method == "LIME" & BA4_SFTLD_weighted_df_summary$n_genes == BA4_SFTLD_weighted_df_summary$optimal_genes] < 0.80 & max(BA4_SFTLD_weighted_df_summary$mean_accuracy_by_group[BA4_SFTLD_weighted_df_summary$celltype == celltype2 &  BA4_SFTLD_weighted_df_summary$method == "LIME"]) >= 0.8 ) {
    
            # Find the smallest n_genes where mean_accuracy_by_group >= 0.80
            new_optimal_genes <- min(BA4_SFTLD_weighted_df_summary$n_genes[BA4_SFTLD_weighted_df_summary$mean_accuracy_by_group >= 0.80 & BA4_SFTLD_weighted_df_summary$celltype == celltype2 & BA4_SFTLD_weighted_df_summary$method == "LIME"])
            
            # Update the optimal_genes column with the new value
            BA4_SFTLD_weighted_df_summary$optimal_genes[BA4_SFTLD_weighted_df_summary$celltype == celltype2] <- new_optimal_genes
            }
        }

        ## custom label
        BA4_SFTLD_weighted_df_summary <- BA4_SFTLD_weighted_df_summary %>%
            mutate(label = paste(celltype, "\n(", optimal_genes, " genes)", sep = ""))

        BA4_SFTLD_weighted_df_summary$method <- factor(BA4_SFTLD_weighted_df_summary$method, levels = c("LIME", "random_HVG", "random_gene"))
        BA4_SFTLD_weighted_df_summary$celltype <- factor(BA4_SFTLD_weighted_df_summary$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))
        
        BA4_SFTLD_weighted_df_summary <- BA4_SFTLD_weighted_df_summary %>%
            arrange(celltype)  # Arrange by the 'celltype' factor

        BA4_SFTLD_weighted_df_summary$label <- factor(BA4_SFTLD_weighted_df_summary$label, levels = unique(BA4_SFTLD_weighted_df_summary$label))
        
    
        ## prin the number of genes corresponding to the optimal for each cell type
        result <- BA4_SFTLD_weighted_df_summary %>% dplyr::select(celltype, optimal_genes)
        result <- result %>% distinct()
        
        write.csv(result, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_lambda_0.00025_gene_set_BA4_SFTLD_narval_2.csv'))

        ## Plot
        ggplot(BA4_SFTLD_weighted_df_summary, aes(x = n_genes, y = mean_accuracy_by_group, colour = method, group = method)) + 
        theme_bw() + 
        geom_vline(aes(xintercept = optimal_genes), linetype = "dashed", color = "black") +
        geom_point(size = 1) +
        geom_line() +
        geom_hline(yintercept = 0.8, colour = "red") +
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
        ylab("Balanced accuracy") +
        xlab("N genes") +
        scale_y_continuous(limits = c(0.5, 1)) +
        scale_colour_manual(values = c("#F18F01", "#2E4057"))
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA4_SFTLD_permutation_dotplot.pdf', height = 2, width = 15)
    ##

    ###########################
    ## C9FTLD BA4 -- Lambda approach and plot -- DONE
    ###########################
    ## code
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
    
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
        
        BA4_C9FTLD_weighted_df_summary <- BA4_C9FTLD_weighted_df %>%
            dplyr::group_by(group, celltype, n_genes, method) %>%
            dplyr::summarise(mean_accuracy_by_group = mean(mean_accuracy, na.rm = TRUE))

        BA4_C9FTLD_weighted_df_summary <- subset(BA4_C9FTLD_weighted_df_summary, n_genes >= 25)

        
        result <- subset(result, lambda == 0.00025)
        result <- result %>% dplyr::select(celltype, n_genes)
        
        colnames(result) <- c("celltype", "optimal_genes")
        
        BA4_C9FTLD_weighted_df_summary <- merge(BA4_C9FTLD_weighted_df_summary, result, by = "celltype")

        ## modify so that if the threshold value is below we take the number of genes corresponding to 80. 
        for(celltype2 in celltype_list){
        if (BA4_C9FTLD_weighted_df_summary$mean_accuracy_by_group[BA4_C9FTLD_weighted_df_summary$celltype == celltype2 &  BA4_C9FTLD_weighted_df_summary$method == "LIME" & BA4_C9FTLD_weighted_df_summary$n_genes == BA4_C9FTLD_weighted_df_summary$optimal_genes] < 0.80 & max(BA4_C9FTLD_weighted_df_summary$mean_accuracy_by_group[BA4_C9FTLD_weighted_df_summary$celltype == celltype2 &  BA4_C9FTLD_weighted_df_summary$method == "LIME"]) >= 0.8 ) {
    
            # Find the smallest n_genes where mean_accuracy_by_group >= 0.80
            new_optimal_genes <- min(BA4_C9FTLD_weighted_df_summary$n_genes[BA4_C9FTLD_weighted_df_summary$mean_accuracy_by_group >= 0.80 & BA4_C9FTLD_weighted_df_summary$celltype == celltype2 & BA4_C9FTLD_weighted_df_summary$method == "LIME"])
            
            # Update the optimal_genes column with the new value
            BA4_C9FTLD_weighted_df_summary$optimal_genes[BA4_C9FTLD_weighted_df_summary$celltype == celltype2] <- new_optimal_genes
            }
        }

        ## custom label
        BA4_C9FTLD_weighted_df_summary <- BA4_C9FTLD_weighted_df_summary %>%
            mutate(label = paste(celltype, "\n(", optimal_genes, " genes)", sep = ""))

        BA4_C9FTLD_weighted_df_summary$method <- factor(BA4_C9FTLD_weighted_df_summary$method, levels = c("LIME", "random_HVG", "random_gene"))
        BA4_C9FTLD_weighted_df_summary$celltype <- factor(BA4_C9FTLD_weighted_df_summary$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))
        
        BA4_C9FTLD_weighted_df_summary <- BA4_C9FTLD_weighted_df_summary %>%
            arrange(celltype)  # Arrange by the 'celltype' factor

        BA4_C9FTLD_weighted_df_summary$label <- factor(BA4_C9FTLD_weighted_df_summary$label, levels = unique(BA4_C9FTLD_weighted_df_summary$label))
        
    
        ## prin the number of genes corresponding to the optimal for each cell type
        result <- BA4_C9FTLD_weighted_df_summary %>% dplyr::select(celltype, optimal_genes)
        result <- result %>% distinct()
        
        write.csv(result, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_lambda_0.00025_gene_set_BA4_C9FTLD_narval_2.csv'))

        ## Plot
        ggplot(BA4_C9FTLD_weighted_df_summary, aes(x = n_genes, y = mean_accuracy_by_group, colour = method, group = method)) + 
        theme_bw() + 
        geom_vline(aes(xintercept = optimal_genes), linetype = "dashed", color = "black") +
        geom_point(size = 1) +
        geom_line() +
        geom_hline(yintercept = 0.8, colour = "red") +
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
        ylab("Balanced accuracy") +
        xlab("N genes") +
        scale_y_continuous(limits = c(0.5, 1)) +
        scale_colour_manual(values = c("#F18F01", "#2E4057"))
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA4_C9FTLD_permutation_dotplot.pdf', height = 2, width = 15)
    ##

    ###########################
    ## SALS BA9 -- Lambda approach and plot -- DONE
    ###########################
    ## code
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        
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

        
        BA9_SALS_weighted_df_summary <- BA9_SALS_weighted_df %>%
            dplyr::group_by(group, celltype, n_genes, method) %>%
            dplyr::summarise(mean_accuracy_by_group = mean(mean_accuracy, na.rm = TRUE))

        BA9_SALS_weighted_df_summary <- subset(BA9_SALS_weighted_df_summary, n_genes >= 25)

        
        result <- subset(result, lambda == 0.00025)
        result <- result %>% dplyr::select(celltype, n_genes)
        
        colnames(result) <- c("celltype", "optimal_genes")
        
        BA9_SALS_weighted_df_summary <- merge(BA9_SALS_weighted_df_summary, result, by = "celltype")

        ## modify so that if the threshold value is below we take the number of genes corresponding to 80. 
        for(celltype2 in celltype_list){
        if (BA9_SALS_weighted_df_summary$mean_accuracy_by_group[BA9_SALS_weighted_df_summary$celltype == celltype2 &  BA9_SALS_weighted_df_summary$method == "LIME" & BA9_SALS_weighted_df_summary$n_genes == BA9_SALS_weighted_df_summary$optimal_genes] < 0.80 & max(BA9_SALS_weighted_df_summary$mean_accuracy_by_group[BA9_SALS_weighted_df_summary$celltype == celltype2 &  BA9_SALS_weighted_df_summary$method == "LIME"]) >= 0.8 ) {
    
            # Find the smallest n_genes where mean_accuracy_by_group >= 0.80
            new_optimal_genes <- min(BA9_SALS_weighted_df_summary$n_genes[BA9_SALS_weighted_df_summary$mean_accuracy_by_group >= 0.80 & BA9_SALS_weighted_df_summary$celltype == celltype2 & BA9_SALS_weighted_df_summary$method == "LIME"])
            
            # Update the optimal_genes column with the new value
            BA9_SALS_weighted_df_summary$optimal_genes[BA9_SALS_weighted_df_summary$celltype == celltype2] <- new_optimal_genes
            }
        }

        ## custom label
        BA9_SALS_weighted_df_summary <- BA9_SALS_weighted_df_summary %>%
            mutate(label = paste(celltype, "\n(", optimal_genes, " genes)", sep = ""))

        BA9_SALS_weighted_df_summary$method <- factor(BA9_SALS_weighted_df_summary$method, levels = c("LIME", "random_HVG", "random_gene"))
        BA9_SALS_weighted_df_summary$celltype <- factor(BA9_SALS_weighted_df_summary$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))
        
        BA9_SALS_weighted_df_summary <- BA9_SALS_weighted_df_summary %>%
            arrange(celltype)  # Arrange by the 'celltype' factor

        BA9_SALS_weighted_df_summary$label <- factor(BA9_SALS_weighted_df_summary$label, levels = unique(BA9_SALS_weighted_df_summary$label))
        
    
        ## prin the number of genes corresponding to the optimal for each cell type
        result <- BA9_SALS_weighted_df_summary %>% dplyr::select(celltype, optimal_genes)
        result <- result %>% distinct()
        
        write.csv(result, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_lambda_0.00025_gene_set_BA9_SALS_narval_2.csv'))

        ## Plot
        ggplot(BA9_SALS_weighted_df_summary, aes(x = n_genes, y = mean_accuracy_by_group, colour = method, group = method)) + 
        theme_bw() + 
        geom_vline(aes(xintercept = optimal_genes), linetype = "dashed", color = "black") +
        geom_point(size = 1) +
        geom_line() +
        geom_hline(yintercept = 0.8, colour = "red") +
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
        ylab("Balanced accuracy") +
        xlab("N genes") +
        scale_y_continuous(limits = c(0.5, 1)) +
        scale_colour_manual(values = c("#F18F01", "#2E4057"))
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA9_SALS_permutation_dotplot.pdf', height = 2, width = 15)
    ##

    ###########################
    ## C9ALS BA9 -- Lambda approach and plot -- DONE
    ###########################
    ## code
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
    
        ## sumarize C9ALS
        BA9_C9ALS_weighted_df$disease_cat <- "C9ALS"
        
        BA9_C9ALS_weighted_df_summary <- BA9_C9ALS_weighted_df %>%
            dplyr::group_by(group, celltype, n_genes, method) %>%
            dplyr::summarise(mean_accuracy_by_group = mean(mean_accuracy, na.rm = TRUE))
        
        BA9_C9ALS_weighted_df_summary <- data.frame(BA9_C9ALS_weighted_df_summary)

        BA9_C9ALS_weighted_df_summary <- subset(BA9_C9ALS_weighted_df_summary, method == "LIME" )
        
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
        
        BA9_C9ALS_weighted_df_summary <- BA9_C9ALS_weighted_df %>%
            dplyr::group_by(group, celltype, n_genes, method) %>%
            dplyr::summarise(mean_accuracy_by_group = mean(mean_accuracy, na.rm = TRUE))

        BA9_C9ALS_weighted_df_summary <- subset(BA9_C9ALS_weighted_df_summary, n_genes >= 25)

        
        result <- subset(result, lambda == 0.00025)
        result <- result %>% dplyr::select(celltype, n_genes)
        
        colnames(result) <- c("celltype", "optimal_genes")
        
        BA9_C9ALS_weighted_df_summary <- merge(BA9_C9ALS_weighted_df_summary, result, by = "celltype")

        ## modify so that if the threshold value is below we take the number of genes corresponding to 80. 
        for(celltype2 in celltype_list){
        if (BA9_C9ALS_weighted_df_summary$mean_accuracy_by_group[BA9_C9ALS_weighted_df_summary$celltype == celltype2 &  BA9_C9ALS_weighted_df_summary$method == "LIME" & BA9_C9ALS_weighted_df_summary$n_genes == BA9_C9ALS_weighted_df_summary$optimal_genes] < 0.80 & max(BA9_C9ALS_weighted_df_summary$mean_accuracy_by_group[BA9_C9ALS_weighted_df_summary$celltype == celltype2 &  BA9_C9ALS_weighted_df_summary$method == "LIME"]) >= 0.8 ) {
    
            # Find the smallest n_genes where mean_accuracy_by_group >= 0.80
            new_optimal_genes <- min(BA9_C9ALS_weighted_df_summary$n_genes[BA9_C9ALS_weighted_df_summary$mean_accuracy_by_group >= 0.80 & BA9_C9ALS_weighted_df_summary$celltype == celltype2 & BA9_C9ALS_weighted_df_summary$method == "LIME"])
            
            # Update the optimal_genes column with the new value
            BA9_C9ALS_weighted_df_summary$optimal_genes[BA9_C9ALS_weighted_df_summary$celltype == celltype2] <- new_optimal_genes
            }
        }

        ## custom label
        BA9_C9ALS_weighted_df_summary <- BA9_C9ALS_weighted_df_summary %>%
            mutate(label = paste(celltype, "\n(", optimal_genes, " genes)", sep = ""))

        BA9_C9ALS_weighted_df_summary$method <- factor(BA9_C9ALS_weighted_df_summary$method, levels = c("LIME", "random_HVG", "random_gene"))
        BA9_C9ALS_weighted_df_summary$celltype <- factor(BA9_C9ALS_weighted_df_summary$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))
        
        BA9_C9ALS_weighted_df_summary <- BA9_C9ALS_weighted_df_summary %>%
            arrange(celltype)  # Arrange by the 'celltype' factor

        BA9_C9ALS_weighted_df_summary$label <- factor(BA9_C9ALS_weighted_df_summary$label, levels = unique(BA9_C9ALS_weighted_df_summary$label))
        
    
        ## prin the number of genes corresponding to the optimal for each cell type
        result <- BA9_C9ALS_weighted_df_summary %>% dplyr::select(celltype, optimal_genes)
        result <- result %>% distinct()
        
        write.csv(result, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_lambda_0.00025_gene_set_BA9_C9ALS_narval_2.csv'))

        ## Plot
        ggplot(BA9_C9ALS_weighted_df_summary, aes(x = n_genes, y = mean_accuracy_by_group, colour = method, group = method)) + 
        theme_bw() + 
        geom_vline(aes(xintercept = optimal_genes), linetype = "dashed", color = "black") +
        geom_point(size = 1) +
        geom_line() +
        geom_hline(yintercept = 0.8, colour = "red") +
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
        ylab("Balanced accuracy") +
        xlab("N genes") +
        scale_y_continuous(limits = c(0.5, 1)) +
        scale_colour_manual(values = c("#F18F01", "#2E4057"))
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA9_C9ALS_permutation_dotplot.pdf', height = 2, width = 15)
    ##
    
    ###########################
    ## SFTLD BA9 -- Lambda approach and plot -- DONE
    ###########################
    ## code
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
    
        ## sumarize SFTLD
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
        
        BA9_SFTLD_weighted_df_summary <- BA9_SFTLD_weighted_df %>%
            dplyr::group_by(group, celltype, n_genes, method) %>%
            dplyr::summarise(mean_accuracy_by_group = mean(mean_accuracy, na.rm = TRUE))

        BA9_SFTLD_weighted_df_summary <- subset(BA9_SFTLD_weighted_df_summary, n_genes >= 25)

        
        result <- subset(result, lambda == 0.00025)
        result <- result %>% dplyr::select(celltype, n_genes)
        
        colnames(result) <- c("celltype", "optimal_genes")
        
        BA9_SFTLD_weighted_df_summary <- merge(BA9_SFTLD_weighted_df_summary, result, by = "celltype")

        ## modify so that if the threshold value is below we take the number of genes corresponding to 80. 
        for(celltype2 in celltype_list){
        if (BA9_SFTLD_weighted_df_summary$mean_accuracy_by_group[BA9_SFTLD_weighted_df_summary$celltype == celltype2 &  BA9_SFTLD_weighted_df_summary$method == "LIME" & BA9_SFTLD_weighted_df_summary$n_genes == BA9_SFTLD_weighted_df_summary$optimal_genes] < 0.80 & max(BA9_SFTLD_weighted_df_summary$mean_accuracy_by_group[BA9_SFTLD_weighted_df_summary$celltype == celltype2 &  BA9_SFTLD_weighted_df_summary$method == "LIME"]) >= 0.8 ) {
    
            # Find the smallest n_genes where mean_accuracy_by_group >= 0.80
            new_optimal_genes <- min(BA9_SFTLD_weighted_df_summary$n_genes[BA9_SFTLD_weighted_df_summary$mean_accuracy_by_group >= 0.80 & BA9_SFTLD_weighted_df_summary$celltype == celltype2 & BA9_SFTLD_weighted_df_summary$method == "LIME"])
            
            # Update the optimal_genes column with the new value
            BA9_SFTLD_weighted_df_summary$optimal_genes[BA9_SFTLD_weighted_df_summary$celltype == celltype2] <- new_optimal_genes
            }
        }

        ## custom label
        BA9_SFTLD_weighted_df_summary <- BA9_SFTLD_weighted_df_summary %>%
            mutate(label = paste(celltype, "\n(", optimal_genes, " genes)", sep = ""))

        BA9_SFTLD_weighted_df_summary$method <- factor(BA9_SFTLD_weighted_df_summary$method, levels = c("LIME", "random_HVG", "random_gene"))
        BA9_SFTLD_weighted_df_summary$celltype <- factor(BA9_SFTLD_weighted_df_summary$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))
        
        BA9_SFTLD_weighted_df_summary <- BA9_SFTLD_weighted_df_summary %>%
            arrange(celltype)  # Arrange by the 'celltype' factor

        BA9_SFTLD_weighted_df_summary$label <- factor(BA9_SFTLD_weighted_df_summary$label, levels = unique(BA9_SFTLD_weighted_df_summary$label))
        
    
        ## prin the number of genes corresponding to the optimal for each cell type
        result <- BA9_SFTLD_weighted_df_summary %>% dplyr::select(celltype, optimal_genes)
        result <- result %>% distinct()
        
        write.csv(result, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_lambda_0.00025_gene_set_BA9_SFTLD_narval_2.csv'))

        ## Plot
        ggplot(BA9_SFTLD_weighted_df_summary, aes(x = n_genes, y = mean_accuracy_by_group, colour = method, group = method)) + 
        theme_bw() + 
        geom_vline(aes(xintercept = optimal_genes), linetype = "dashed", color = "black") +
        geom_point(size = 1) +
        geom_line() +
        geom_hline(yintercept = 0.8, colour = "red") +
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
        ylab("Balanced accuracy") +
        xlab("N genes") +
        scale_y_continuous(limits = c(0.5, 1)) +
        scale_colour_manual(values = c("#F18F01", "#2E4057"))
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA9_SFTLD_permutation_dotplot.pdf', height = 2, width = 15)
    ##

    ###########################
    ## C9FTLD BA9 -- Lambda approach and plot -- DONE
    ###########################
    ## code
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
    
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
        
        BA9_C9FTLD_weighted_df_summary <- BA9_C9FTLD_weighted_df %>%
            dplyr::group_by(group, celltype, n_genes, method) %>%
            dplyr::summarise(mean_accuracy_by_group = mean(mean_accuracy, na.rm = TRUE))

        BA9_C9FTLD_weighted_df_summary <- subset(BA9_C9FTLD_weighted_df_summary, n_genes >= 25)

        
        result <- subset(result, lambda == 0.00025)
        result <- result %>% dplyr::select(celltype, n_genes)
        
        colnames(result) <- c("celltype", "optimal_genes")
        
        BA9_C9FTLD_weighted_df_summary <- merge(BA9_C9FTLD_weighted_df_summary, result, by = "celltype")

        ## modify so that if the threshold value is below we take the number of genes corresponding to 80. 
        for(celltype2 in celltype_list){
        if (BA9_C9FTLD_weighted_df_summary$mean_accuracy_by_group[BA9_C9FTLD_weighted_df_summary$celltype == celltype2 &  BA9_C9FTLD_weighted_df_summary$method == "LIME" & BA9_C9FTLD_weighted_df_summary$n_genes == BA9_C9FTLD_weighted_df_summary$optimal_genes] < 0.80 & max(BA9_C9FTLD_weighted_df_summary$mean_accuracy_by_group[BA9_C9FTLD_weighted_df_summary$celltype == celltype2 &  BA9_C9FTLD_weighted_df_summary$method == "LIME"]) >= 0.8 ) {
    
            # Find the smallest n_genes where mean_accuracy_by_group >= 0.80
            new_optimal_genes <- min(BA9_C9FTLD_weighted_df_summary$n_genes[BA9_C9FTLD_weighted_df_summary$mean_accuracy_by_group >= 0.80 & BA9_C9FTLD_weighted_df_summary$celltype == celltype2 & BA9_C9FTLD_weighted_df_summary$method == "LIME"])
            
            # Update the optimal_genes column with the new value
            BA9_C9FTLD_weighted_df_summary$optimal_genes[BA9_C9FTLD_weighted_df_summary$celltype == celltype2] <- new_optimal_genes
            }
        }

        ## custom label
        BA9_C9FTLD_weighted_df_summary <- BA9_C9FTLD_weighted_df_summary %>%
            mutate(label = paste(celltype, "\n(", optimal_genes, " genes)", sep = ""))

        BA9_C9FTLD_weighted_df_summary$method <- factor(BA9_C9FTLD_weighted_df_summary$method, levels = c("LIME", "random_HVG", "random_gene"))
        BA9_C9FTLD_weighted_df_summary$celltype <- factor(BA9_C9FTLD_weighted_df_summary$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))
        
        BA9_C9FTLD_weighted_df_summary <- BA9_C9FTLD_weighted_df_summary %>%
            arrange(celltype)  # Arrange by the 'celltype' factor

        BA9_C9FTLD_weighted_df_summary$label <- factor(BA9_C9FTLD_weighted_df_summary$label, levels = unique(BA9_C9FTLD_weighted_df_summary$label))
        
    
        ## prin the number of genes corresponding to the optimal for each cell type
        result <- BA9_C9FTLD_weighted_df_summary %>% dplyr::select(celltype, optimal_genes)
        result <- result %>% distinct()
        
        write.csv(result, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_lambda_0.00025_gene_set_BA9_C9FTLD_narval_2.csv'))

        ## Plot
        ggplot(BA9_C9FTLD_weighted_df_summary, aes(x = n_genes, y = mean_accuracy_by_group, colour = method, group = method)) + 
        theme_bw() + 
        geom_vline(aes(xintercept = optimal_genes), linetype = "dashed", color = "black") +
        geom_point(size = 1) +
        geom_line() +
        geom_hline(yintercept = 0.8, colour = "red") +
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
        ylab("Balanced accuracy") +
        xlab("N genes") +
        scale_y_continuous(limits = c(0.5, 1)) +
        scale_colour_manual(values = c("#F18F01", "#2E4057"))
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA9_C9FTLD_permutation_dotplot.pdf', height = 2, width = 15)
    ##

##

############################################################################ Print file containing optimal gene set genes for each DNN group --> we will use this in downstream analyses -- First for 5CV
############################################################################
############################################################################
############################################################################ 
## code
    ###########################
    ## SALS BA4
    ###########################
    ## code
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        #celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",  "OPC",  "Micro", "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA4"
        par_status = "SALS"
        par_prep = "CombatSeq"

        for(celltype2 in unique(celltype_list)){
            
            ## read in optimal number of genes, obtained from the above code
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_lambda_0.00025_gene_set_',par_brain_region,'_',par_status,'_narval_2.csv'), sep = ",")
            df_lim <- subset(df, celltype == celltype2 )

            ## read in the LIME weighted percentile rank file which was created by and used in the permutation scripts
            df_rank = read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_percentile_rank_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval_2.csv'), sep = ",")

            # Extract the optimal_genes value from df_lim
            n <- df_lim$optimal_genes[1]

            # Subset the first n columns of df_rank
            df_rank_subset <- df_rank[ 1:n,]

            # View the result
            print(nrow(df_rank_subset) == n)

            write.csv(df_rank_subset, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_lambda_0.00025_optimal_genes_',par_brain_region,'_',par_status,'_',celltype2,'_narval_2.csv'))
            
        }
    ##

    ###########################
    ## C9ALS BA4
    ###########################
    ## code
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",   "Mural",   "Endo",    "Fibro", "L5")
        #celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",   "Astro",   "OPC",  "Micro", "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA4"
        par_status = "C9ALS"
        par_prep = "CombatSeq"

        for(celltype2 in unique(celltype_list)){
            
            ## read in optimal number of genes, obtained from the above code
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_lambda_0.00025_gene_set_',par_brain_region,'_',par_status,'_narval_2.csv'), sep = ",")
            df_lim <- subset(df, celltype == celltype2 )

            ## read in the LIME weighted percentile rank file which was created by and used in the permutation scripts
            df_rank = read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_percentile_rank_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval_2.csv'), sep = ",")

            # Extract the optimal_genes value from df_lim
            n <- df_lim$optimal_genes[1]

            # Subset the first n columns of df_rank
            df_rank_subset <- df_rank[ 1:n,]

            # View the result
            print(nrow(df_rank_subset) == n)

            write.csv(df_rank_subset, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_lambda_0.00025_optimal_genes_',par_brain_region,'_',par_status,'_',celltype2,'_narval_2.csv'))
            
        }
    ##

    ###########################
    ## SFTLD BA4
    ###########################
    ## code
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro", "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA4"
        par_status = "SFTLD"
        par_prep = "CombatSeq"

        for(celltype2 in unique(celltype_list)){
            
            ## read in optimal number of genes, obtained from the above code
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_lambda_0.00025_gene_set_',par_brain_region,'_',par_status,'_narval_2.csv'), sep = ",")
            df_lim <- subset(df, celltype == celltype2 )

            ## read in the LIME weighted percentile rank file which was created by and used in the permutation scripts
            df_rank = read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_percentile_rank_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval_2.csv'), sep = ",")

            # Extract the optimal_genes value from df_lim
            n <- df_lim$optimal_genes[1]

            # Subset the first n columns of df_rank
            df_rank_subset <- df_rank[ 1:n,]

            # View the result
            print(nrow(df_rank_subset) == n)

            write.csv(df_rank_subset, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_lambda_0.00025_optimal_genes_',par_brain_region,'_',par_status,'_',celltype2,'_narval_2.csv'))
            
        }
    ##

    ###########################
    ## C9FTLD BA4
    ###########################
    ## code
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA4"
        par_status = "C9FTLD"
        par_prep = "CombatSeq"

        for(celltype2 in unique(celltype_list)){
            
            ## read in optimal number of genes, obtained from the above code
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_lambda_0.00025_gene_set_',par_brain_region,'_',par_status,'_narval_2.csv'), sep = ",")
            df_lim <- subset(df, celltype == celltype2 )

            ## read in the LIME weighted percentile rank file which was created by and used in the permutation scripts
            df_rank = read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_percentile_rank_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval_2.csv'), sep = ",")

            # Extract the optimal_genes value from df_lim
            n <- df_lim$optimal_genes[1]

            # Subset the first n columns of df_rank
            df_rank_subset <- df_rank[ 1:n,]

            # View the result
            print(nrow(df_rank_subset) == n)

            write.csv(df_rank_subset, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_lambda_0.00025_optimal_genes_',par_brain_region,'_',par_status,'_',celltype2,'_narval_2.csv'))
            
        }
    ##

    ###########################
    ## SALS BA9
    ###########################
    ## code
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA9"
        par_status = "SALS"
        par_prep = "CombatSeq"

        for(celltype2 in unique(celltype_list)){
            
            ## read in optimal number of genes, obtained from the above code
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_lambda_0.00025_gene_set_',par_brain_region,'_',par_status,'_narval_2.csv'), sep = ",")
            df_lim <- subset(df, celltype == celltype2 )

            ## read in the LIME weighted percentile rank file which was created by and used in the permutation scripts
            df_rank = read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_percentile_rank_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval_2.csv'), sep = ",")

            # Extract the optimal_genes value from df_lim
            n <- df_lim$optimal_genes[1]

            # Subset the first n columns of df_rank
            df_rank_subset <- df_rank[ 1:n,]

            # View the result
            print(nrow(df_rank_subset) == n)

            write.csv(df_rank_subset, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_lambda_0.00025_optimal_genes_',par_brain_region,'_',par_status,'_',celltype2,'_narval_2.csv'))
            
        }
    ##

    ###########################
    ## C9ALS BA9
    ###########################
    ## code
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA9"
        par_status = "C9ALS"
        par_prep = "CombatSeq"

        for(celltype2 in unique(celltype_list)){
            
            ## read in optimal number of genes, obtained from the above code
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_lambda_0.00025_gene_set_',par_brain_region,'_',par_status,'_narval_2.csv'), sep = ",")
            df_lim <- subset(df, celltype == celltype2 )

            ## read in the LIME weighted percentile rank file which was created by and used in the permutation scripts
            df_rank = read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_percentile_rank_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval_2.csv'), sep = ",")

            # Extract the optimal_genes value from df_lim
            n <- df_lim$optimal_genes[1]

            # Subset the first n columns of df_rank
            df_rank_subset <- df_rank[ 1:n,]

            # View the result
            print(nrow(df_rank_subset) == n)

            write.csv(df_rank_subset, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_lambda_0.00025_optimal_genes_',par_brain_region,'_',par_status,'_',celltype2,'_narval_2.csv'))
            
        }
    ##

    ###########################
    ## SFTLD BA9
    ###########################
    ## code
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA9"
        par_status = "SFTLD"
        par_prep = "CombatSeq"

        for(celltype2 in unique(celltype_list)){
            
            ## read in optimal number of genes, obtained from the above code
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_lambda_0.00025_gene_set_',par_brain_region,'_',par_status,'_narval_2.csv'), sep = ",")
            df_lim <- subset(df, celltype == celltype2 )

            ## read in the LIME weighted percentile rank file which was created by and used in the permutation scripts
            df_rank = read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_percentile_rank_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval_2.csv'), sep = ",")

            # Extract the optimal_genes value from df_lim
            n <- df_lim$optimal_genes[1]

            # Subset the first n columns of df_rank
            df_rank_subset <- df_rank[ 1:n,]

            # View the result
            print(nrow(df_rank_subset) == n)

            write.csv(df_rank_subset, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_lambda_0.00025_optimal_genes_',par_brain_region,'_',par_status,'_',celltype2,'_narval_2.csv'))
            
        }
    ##

    ###########################
    ## C9FTLD BA9
    ###########################
    ## code
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro", "Mural",   "Endo",    "Fibro", "L5")

        par_brain_region = "BA9"
        par_status = "C9FTLD"
        par_prep = "CombatSeq"

        for(celltype2 in unique(celltype_list)){
            
            ## read in optimal number of genes, obtained from the above code
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_lambda_0.00025_gene_set_',par_brain_region,'_',par_status,'_narval_2.csv'), sep = ",")
            df_lim <- subset(df, celltype == celltype2 )

            ## read in the LIME weighted percentile rank file which was created by and used in the permutation scripts
            df_rank = read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_percentile_rank_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval_2.csv'), sep = ",")

            # Extract the optimal_genes value from df_lim
            n <- df_lim$optimal_genes[1]

            # Subset the first n columns of df_rank
            df_rank_subset <- df_rank[ 1:n,]

            # View the result
            print(nrow(df_rank_subset) == n)

            write.csv(df_rank_subset, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_lambda_0.00025_optimal_genes_',par_brain_region,'_',par_status,'_',celltype2,'_narval_2.csv'))
            
        }
    ##
##

############################################################################ Read in all of the files describing the optimal gene sets and create one file set
############################################################################
############################################################################
############################################################################ 
## code
    # Set variables
    celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",  "PV",      "5HT3aR",  "Rosehip", "SOM",  "OPC",  "Micro", "Mural",   "Endo",    "Fibro", "L5")
    region_list = c("BA4", "BA9")
    status_list = c("SALS", "C9ALS", "SFTLD", "C9FTLD")
    par_prep = "CombatSeq"

    # Create an empty list to store data frames
    all_data <- list()

    # Loop through each celltype in the list
    for(par_brain_region in region_list){
        for(par_status in status_list){
            for (celltype in celltype_list) {
            # Define the file path using the current celltype
            file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_lambda_0.00025_optimal_genes_', 
                                par_brain_region, '_', par_status, '_', celltype, '_narval_2.csv')
            
            # Read the CSV file into a data frame
            df <- read.csv(file_path)
            
            # Append the data frame to the list
            all_data[[paste0(par_brain_region, '_', par_status, '_',celltype)]] <- df
            }
        }
    }

    # Combine all the data frames by row
    combined_data <- do.call(rbind, all_data)

    # Add the category
    combined_data$category <- sub("\\..*", "", rownames(combined_data))
    length(unique(combined_data$feature))

    # View the combined data
    head(combined_data)

    # Create a dataframe and print it so that it can be used for downstream analyses.
    write.csv(combined_data, '/home/fiorini9/scratch/machine_learning_ALS/LIME_gene_set_outputs/temp_pineda_LIME_lambda_0.00025_optimal_gene_sets.csv')

    ## combined_data_lim <- subset(combined_data, feature == "DNAJC7")
    ## only retain protein coding
    ## remove duplicate AS
    ## remove LINC
    ## Remove MT genes
    ## remove ribosomal proteins
    
##



############################################################################ Plot 5fold optimal accuracy for 5CV workflow -- run this after permutation testing (e.g. Optimal_5CV_SALS_MCx_DNN_LIME_narval_weighted_2.py )
############################################################################
############################################################################
############################################################################ 
## code
    
    ## BA4
    
    ###########################
    ## SALS BA4
    ###########################
    ## code
        ## define variables
        condition = "SALS" #one ofSALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_5CV_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            test_column <- c(rep("LIME", 5), rep("Random", 25))
            result_df <- data.frame(accuracy = accuracies, test = test_column)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region
            result_df$n_genes <- unique(temp$n_genes)

            result_list[[i]] <- result_df
        }

        final_df_SALS_BA4 <- do.call(rbind, result_list)
    ##

    ###########################
    ## C9ALS BA4
    ###########################
    ## code
        ## define variables
        condition = "C9ALS" #one ofSALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_5CV_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            test_column <- c(rep("LIME", 5), rep("Random", 25))
            result_df <- data.frame(accuracy = accuracies, test = test_column)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region
            result_df$n_genes <- unique(temp$n_genes)

            result_list[[i]] <- result_df
        }

        final_df_C9ALS_BA4 <- do.call(rbind, result_list)
    ##

    ###########################
    ## SFTLD BA4
    ###########################
    ## code
        ## define variables
        condition = "SFTLD" #one ofSALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_5CV_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            test_column <- c(rep("LIME", 5), rep("Random", 25))
            result_df <- data.frame(accuracy = accuracies, test = test_column)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region
            result_df$n_genes <- unique(temp$n_genes)

            result_list[[i]] <- result_df
        }

        final_df_SFTLD_BA4 <- do.call(rbind, result_list)
    ##

    ###########################
    ## C9FTLD BA4
    ###########################
    ## code
        ## define variables
        condition = "C9FTLD" #one ofSALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_5CV_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            test_column <- c(rep("LIME", 5), rep("Random", 25))
            result_df <- data.frame(accuracy = accuracies, test = test_column)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region
            result_df$n_genes <- unique(temp$n_genes)

            result_list[[i]] <- result_df
        }

        final_df_C9FTLD_BA4 <- do.call(rbind, result_list)
    ##

    ###########################
    ## merge and plot BA4
    ###########################

    ## code for barplot
        BA4_bind_all <- rbind(final_df_SALS_BA4, final_df_C9ALS_BA4, final_df_SFTLD_BA4, final_df_C9FTLD_BA4)

        n_gene_df <- BA4_bind_all %>% dplyr::select(celltype, condition,  n_genes)
        
        ## Calculate the mean and standard deviation
        summary_stats <- BA4_bind_all %>%
            group_by(celltype, condition, test) %>%
            summarise(
                median_accuracy = median(accuracy), 
                mean_accuracy = mean(accuracy),           # Calculate the median of X0
                sd_accuracy = sd(accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(accuracy - median(accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )
        
        summary_stats <- merge(summary_stats, n_gene_df, by = c('celltype', 'condition') )

        ## set factor levels
        summary_stats$condition <- factor(summary_stats$condition, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD"))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))

        summary_stats$n_genes[summary_stats$test == "Random"] <- NA

        summary_stats <- summary_stats %>% distinct()
        
        ## plot
        ggplot(summary_stats, aes(x = mean_accuracy, y = celltype, fill = test, label = n_genes)) + 
        theme_bw() +
        geom_text(aes(x = 1.15), size = 3)+
        geom_vline(xintercept = 1.05) +
        geom_vline(xintercept = 0.8, colour = "red") +
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
        xlab("Balanced accuracy") +
        scale_fill_manual(values = c("#F18F01", "#2E4057")) +
        coord_cartesian(xlim=c(0.5,1.2)) +
        facet_grid(~condition) +
        scale_x_continuous(
            breaks = c(0.5, 0.75, 1.00),          # Set the tick positions
            labels = c("0.5", "0.75", "1.0")      # Optional: Customize the tick labels
        )
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/Optimal_5CV_for_5CV_worflow_barplot_BA4.pdf', height = 4, width = 7)
    ##

    ## BA9
    
    ###########################
    ## SALS BA9
    ###########################
    ## code
        ## define variables
        condition = "SALS" #one ofSALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_5CV_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            test_column <- c(rep("LIME", 5), rep("Random", 25))
            result_df <- data.frame(accuracy = accuracies, test = test_column)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region
            result_df$n_genes <- unique(temp$n_genes)

            result_list[[i]] <- result_df
        }

        final_df_SALS_BA9 <- do.call(rbind, result_list)
    ##

    ###########################
    ## C9ALS BA9
    ###########################
    ## code
        ## define variables
        condition = "C9ALS" #one ofSALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_5CV_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            test_column <- c(rep("LIME", 5), rep("Random", 25))
            result_df <- data.frame(accuracy = accuracies, test = test_column)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region
            result_df$n_genes <- unique(temp$n_genes)

            result_list[[i]] <- result_df
        }

        final_df_C9ALS_BA9 <- do.call(rbind, result_list)
    ##

    ###########################
    ## SFTLD BA9
    ###########################
    ## code
        ## define variables
        condition = "SFTLD" #one ofSALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_5CV_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            test_column <- c(rep("LIME", 5), rep("Random", 25))
            result_df <- data.frame(accuracy = accuracies, test = test_column)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region
            result_df$n_genes <- unique(temp$n_genes)

            result_list[[i]] <- result_df
        }

        final_df_SFTLD_BA9 <- do.call(rbind, result_list)
    ##

    ###########################
    ## C9FTLD BA9
    ###########################
    ## code
        ## define variables
        condition = "C9FTLD" #one ofSALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_5CV_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            test_column <- c(rep("LIME", 5), rep("Random", 25))
            result_df <- data.frame(accuracy = accuracies, test = test_column)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region
            result_df$n_genes <- unique(temp$n_genes)

            result_list[[i]] <- result_df
        }

        final_df_C9FTLD_BA9 <- do.call(rbind, result_list)
    ##

    ###########################
    ## merge and plot BA9
    ###########################

    ## code for barplot
        BA9_bind_all <- rbind(final_df_SALS_BA9, final_df_C9ALS_BA9, final_df_SFTLD_BA9, final_df_C9FTLD_BA9)

        n_gene_df <- BA9_bind_all %>% dplyr::select(celltype, condition,  n_genes)
        
        ## Calculate the mean and standard deviation
        summary_stats <- BA9_bind_all %>%
            group_by(celltype, condition, test) %>%
            summarise(
                median_accuracy = median(accuracy), 
                mean_accuracy = mean(accuracy),           # Calculate the median of X0
                sd_accuracy = sd(accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(accuracy - median(accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )
        
        summary_stats <- merge(summary_stats, n_gene_df, by = c('celltype', 'condition') )

        ## set factor levels
        summary_stats$condition <- factor(summary_stats$condition, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD"))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))

        summary_stats$n_genes[summary_stats$test == "Random"] <- NA

        summary_stats <- summary_stats %>% distinct()
        
        ## plot
        ggplot(summary_stats, aes(x = mean_accuracy, y = celltype, fill = test, label = n_genes)) + 
        theme_bw() +
        geom_text(aes(x = 1.15), size = 3)+
        geom_vline(xintercept = 1.05) +
        geom_vline(xintercept = 0.8, colour = "red") +
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
        xlab("Balanced accuracy") +
        scale_fill_manual(values = c("#F18F01", "#2E4057")) +
        coord_cartesian(xlim=c(0.5,1.2)) +
        facet_grid(~condition) +
        scale_x_continuous(
            breaks = c(0.5, 0.75, 1.00),          # Set the tick positions
            labels = c("0.5", "0.75", "1.0")      # Optional: Customize the tick labels
        )
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/Optimal_5CV_for_5CV_worflow_barplot_BA9.pdf', height = 4, width = 7)
    ##
        
        
        
        
    
##





############################################################################ Plot LOO optimal accuracy for 5CV workflow
############################################################################
############################################################################
############################################################################ 

## code
    
    ## BA4
    
    ###########################
    ## SALS BA4
    ###########################
    ## code
        ## define variables
        condition = "SALS" #one ofSALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_LOO_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            
            result_df <- temp %>%
            group_by(donor, method, n_genes) %>%
            summarise(
                mean_accuracy = mean(test_accuracy),           # Calculate the median of X0
                .groups = "drop"                 # Drop the grouping after summarising
                )

            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df
        }

        final_df_SALS_BA4 <- do.call(rbind, result_list)
    ##

    ###########################
    ## C9ALS BA4
    ###########################
    ## code
        ## define variables
        condition = "C9ALS" #one ofSALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_LOO_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            
            result_df <- temp %>%
            group_by(donor, method, n_genes) %>%
            summarise(
                mean_accuracy = mean(test_accuracy),           # Calculate the median of X0
                .groups = "drop"                 # Drop the grouping after summarising
                )

            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df
        }

        final_df_C9ALS_BA4 <- do.call(rbind, result_list)
    ##

    ###########################
    ## SFTLD BA4
    ###########################
    ## code
        ## define variables
        condition = "SFTLD" #one ofSALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_LOO_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            
            result_df <- temp %>%
            group_by(donor, method, n_genes) %>%
            summarise(
                mean_accuracy = mean(test_accuracy),           # Calculate the median of X0
                .groups = "drop"                 # Drop the grouping after summarising
                )

            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df
        }

        final_df_SFTLD_BA4 <- do.call(rbind, result_list)
    ##

    ###########################
    ## C9FTLD BA4
    ###########################
    ## code
        ## define variables
        condition = "C9FTLD" #one ofSALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_LOO_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            
            result_df <- temp %>%
            group_by(donor, method, n_genes) %>%
            summarise(
                mean_accuracy = mean(test_accuracy),           # Calculate the median of X0
                .groups = "drop"                 # Drop the grouping after summarising
                )

            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df
        }

        final_df_C9FTLD_BA4 <- do.call(rbind, result_list)
    ##

    ###########################
    ## merge and plot BA4
    ###########################

    ## code for barplot
        BA4_bind_all <- rbind(final_df_SALS_BA4, final_df_C9ALS_BA4, final_df_SFTLD_BA4, final_df_C9FTLD_BA4)

        n_gene_df <- BA4_bind_all %>% dplyr::select(celltype, condition,  n_genes)
        
        ## Calculate the mean and standard deviation
        summary_stats <- BA4_bind_all %>%
            group_by(celltype, condition, method) %>%
            summarise(
                median_accuracy = median(mean_accuracy), 
                mean2_accuracy = mean(mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(mean_accuracy - median(mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )
        
        summary_stats <- merge(summary_stats, n_gene_df, by = c('celltype', 'condition') )

        ## set factor levels
        summary_stats$condition <- factor(summary_stats$condition, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD"))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))

        summary_stats$n_genes[summary_stats$method == "random"] <- NA

        summary_stats <- summary_stats %>% distinct()
        
        ## plot
        ggplot(summary_stats, aes(x = median_accuracy, y = celltype, fill = method, label = n_genes)) + 
        theme_bw() +
        geom_text(aes(x = 1.15), size = 3)+
        geom_vline(xintercept = 1.05) +
        geom_vline(xintercept = 0.8, colour = "red") +
        geom_bar(stat = 'identity', position = "dodge") +
        geom_errorbar(aes(xmin=median_accuracy-0, xmax=median_accuracy+0), width=.2,
                position=position_dodge(.9)) +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(colour = "black"),
            axis.text.y = element_text(colour = "black"),
            axis.title.y = element_blank(),
            axis.title.x = element_text(face = "bold")
        ) +
        xlab("Balanced accuracy") +
        scale_fill_manual(values = c("#F18F01", "#2E4057")) +
        coord_cartesian(xlim=c(0.5,1.2)) +
        facet_grid(~condition) +
        scale_x_continuous(
            breaks = c(0.5, 0.75, 1.00),          # Set the tick positions
            labels = c("0.5", "0.75", "1.0")      # Optional: Customize the tick labels
        )
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/Optimal_LOO_for_5CV_worflow_barplot_BA4.pdf', height = 4, width = 7)
    ##

    ## BA9
    
    ###########################
    ## SALS BA9
    ###########################
    ## code
        ## define variables
        condition = "SALS" #one ofSALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_LOO_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            
            result_df <- temp %>%
            group_by(donor, method, n_genes) %>%
            summarise(
                mean_accuracy = mean(test_accuracy),           # Calculate the median of X0
                .groups = "drop"                 # Drop the grouping after summarising
                )

            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df
        }

        final_df_SALS_BA9 <- do.call(rbind, result_list)
    ##

    ###########################
    ## C9ALS BA9
    ###########################
    ## code
        ## define variables
        condition = "C9ALS" #one ofSALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_LOO_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            
            result_df <- temp %>%
            group_by(donor, method, n_genes) %>%
            summarise(
                mean_accuracy = mean(test_accuracy),           # Calculate the median of X0
                .groups = "drop"                 # Drop the grouping after summarising
                )

            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df
        }

        final_df_C9ALS_BA9 <- do.call(rbind, result_list)
    ##

    ###########################
    ## SFTLD BA9
    ###########################
    ## code
        ## define variables
        condition = "SFTLD" #one ofSALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_LOO_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            
            result_df <- temp %>%
            group_by(donor, method, n_genes) %>%
            summarise(
                mean_accuracy = mean(test_accuracy),           # Calculate the median of X0
                .groups = "drop"                 # Drop the grouping after summarising
                )

            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df
        }

        final_df_SFTLD_BA9 <- do.call(rbind, result_list)
    ##

    ###########################
    ## C9FTLD BA9
    ###########################
    ## code
        ## define variables
        condition = "C9FTLD" #one ofSALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_LOO_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            
            result_df <- temp %>%
            group_by(donor, method, n_genes) %>%
            summarise(
                mean_accuracy = mean(test_accuracy),           # Calculate the median of X0
                .groups = "drop"                 # Drop the grouping after summarising
                )

            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df
        }

        final_df_C9FTLD_BA9 <- do.call(rbind, result_list)
    ##

    ###########################
    ## merge and plot BA9
    ###########################

    ## code for barplot
        BA9_bind_all <- rbind(final_df_SALS_BA9, final_df_C9ALS_BA9, final_df_SFTLD_BA9, final_df_C9FTLD_BA9)

        n_gene_df <- BA9_bind_all %>% dplyr::select(celltype, condition,  n_genes)
        
        ## Calculate the mean and standard deviation
        summary_stats <- BA9_bind_all %>%
            group_by(celltype, condition, method) %>%
            summarise(
                median_accuracy = median(mean_accuracy), 
                mean2_accuracy = mean(mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(mean_accuracy - median(mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )
        
        summary_stats <- merge(summary_stats, n_gene_df, by = c('celltype', 'condition') )

        ## set factor levels
        summary_stats$condition <- factor(summary_stats$condition, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD"))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))

        summary_stats$n_genes[summary_stats$method == "random"] <- NA

        summary_stats <- summary_stats %>% distinct()
        
        ## plot
        ggplot(summary_stats, aes(x = median_accuracy, y = celltype, fill = method, label = n_genes)) + 
        theme_bw() +
        geom_text(aes(x = 1.15), size = 3)+
        geom_vline(xintercept = 1.05) +
        geom_vline(xintercept = 0.8, colour = "red") +
        geom_bar(stat = 'identity', position = "dodge") +
        geom_errorbar(aes(xmin=median_accuracy-0, xmax=median_accuracy+0), width=.2,
                position=position_dodge(.9)) +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(colour = "black"),
            axis.text.y = element_text(colour = "black"),
            axis.title.y = element_blank(),
            axis.title.x = element_text(face = "bold")
        ) +
        xlab("Balanced accuracy") +
        scale_fill_manual(values = c("#F18F01", "#2E4057")) +
        coord_cartesian(xlim=c(0.5,1.2)) +
        facet_grid(~condition) +
        scale_x_continuous(
            breaks = c(0.5, 0.75, 1.00),          # Set the tick positions
            labels = c("0.5", "0.75", "1.0")      # Optional: Customize the tick labels
        )
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/Optimal_LOO_for_5CV_worflow_barplot_BA9.pdf', height = 4, width = 7)
    ##
