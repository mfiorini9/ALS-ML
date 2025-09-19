
## We are running this in Narval
salloc -A def-tdurcan --time=0-4 -c 1 --mem=40g

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

#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
# The following initializes usage of Bioc devel
#BiocManager::install(version='devel')
#BiocManager::install("MAST")
#library(MAST, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
#install.packages("MAST", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")

################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ 5 CV workflow: Number of highly variable genes
## code
    
    #BA4
    
    ###########################
    ## SALS BA4
    ###########################
    ## code
        ## define variables
        condition = "SALS" #one ofSALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            HVGs <- as.numeric(temp$n_HVGs)
            result_df <- data.frame(HVGs = HVGs)
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
        condition = "C9ALS" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        #cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            HVGs <- as.numeric(temp$n_HVGs)
            result_df <- data.frame(HVGs = HVGs)
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
        condition = "SFTLD" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        #cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            HVGs <- as.numeric(temp$n_HVGs)
            result_df <- data.frame(HVGs = HVGs)
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
        condition = "C9FTLD" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        #cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            HVGs <- as.numeric(temp$n_HVGs)
            result_df <- data.frame(HVGs = HVGs)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df

        }

        final_df_C9FTLD_BA4 <- do.call(rbind, result_list)
    ##

    ###########################
    ## all_conditions BA4
    ###########################
    ## code
        ## define variables
        condition = "all_conditions" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        #cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')


        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            HVGs <- as.numeric(temp$n_HVGs)
            result_df <- data.frame(HVGs = HVGs)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df

        }

        final_df_all_conditions_BA4 <- do.call(rbind, result_list)
    ##

    ###########################
    ## merge and plot BA4
    ###########################
    ## code
        BA4_bind_all <- rbind(final_df_SALS_BA4, final_df_C9ALS_BA4, final_df_SFTLD_BA4, final_df_C9FTLD_BA4, final_df_all_conditions_BA4)
        BA4_bind_all$condition[BA4_bind_all$condition == "all_conditions"] <- "All conditions"

        ## set factor levels
        BA4_bind_all$condition <- factor(BA4_bind_all$condition, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All conditions"))
        #BA4_bind_all$celltype <- factor(BA4_bind_all$celltype, levels = c('L4_L6','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Astro', 'OPC','T_Cell'  ))
        BA4_bind_all$celltype <- factor(BA4_bind_all$celltype, levels = c('L4_L6','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Astro', 'OPC'))


        ggplot(BA4_bind_all, aes(x = condition, y = HVGs, fill = condition, colour = condition)) + 
        theme_bw() + 
        geom_bar(stat = 'identity') +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~celltype, scales = "free_x") +
        ylab("n HVGs") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "black"))
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA4_nHVGs.pdf', height = 2, width = 13)
    ##


    ## BA9
    
    
    ###########################
    ## SALS BA9
    ###########################
    ## code
        ## define variables
        condition = "SALS" #one ofSALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        #cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            HVGs <- as.numeric(temp$n_HVGs)
            result_df <- data.frame(HVGs = HVGs)
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
        condition = "C9ALS" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        #cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            HVGs <- as.numeric(temp$n_HVGs)
            result_df <- data.frame(HVGs = HVGs)
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
        condition = "SFTLD" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        #cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            HVGs <- as.numeric(temp$n_HVGs)
            result_df <- data.frame(HVGs = HVGs)
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
        condition = "C9FTLD" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        #cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            HVGs <- as.numeric(temp$n_HVGs)
            result_df <- data.frame(HVGs = HVGs)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df

        }

        final_df_C9FTLD_BA9 <- do.call(rbind, result_list)
    ##

    ###########################
    ## all_conditions BA9
    ###########################
    ## code
        ## define variables
        condition = "all_conditions" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        #cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            HVGs <- as.numeric(temp$n_HVGs)
            result_df <- data.frame(HVGs = HVGs)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df

        }

        final_df_all_conditions_BA9 <- do.call(rbind, result_list)
    ##

    ###########################
    ## merge and plot BA9
    ###########################
    ## code
        BA9_bind_all <- rbind(final_df_SALS_BA9, final_df_C9ALS_BA9, final_df_SFTLD_BA9, final_df_C9FTLD_BA9, final_df_all_conditions_BA9)
        BA9_bind_all$condition[BA9_bind_all$condition == "all_conditions"] <- "All conditions"

        ## set factor levels
        BA9_bind_all$condition <- factor(BA9_bind_all$condition, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All conditions"))
        #BA9_bind_all$celltype <- factor(BA9_bind_all$celltype, levels = c('L4_L6','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Astro', 'OPC','T_Cell'  ))
        BA9_bind_all$celltype <- factor(BA9_bind_all$celltype, levels = c('L4_L6','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Astro', 'OPC'))


        ggplot(BA9_bind_all, aes(x = condition, y = HVGs, fill = condition, colour = condition)) + 
        theme_bw() + 
        geom_bar(stat = 'identity') +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~celltype, scales = "free_x") +
        ylab("n HVGs") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "black"))
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA9_nHVGs.pdf', height = 2, width = 13)
    ##
##


################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ 5 CV workflow: Number of cells
## code
    ###########################
    ## SALS BA4
    ###########################
    ## code
        ## define variables
        condition = "SALS" #one ofSALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        cell_type_list <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        result_list <- list()

        for(i in unique(cell_type_list)){
            file_1 <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/test_meta_",condition,"_",region,"_",i,"_pineda_narval.csv")
            temp_1 <- read.delim(file_1, header = T, sep = ",")
            length_1 <- nrow(temp_1)

            file_2 <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/train_meta_",condition,"_",region,"_",i,"_pineda_narval.csv")
            temp_2 <- read.delim(file_2, header = T, sep = ",")
            length_2 <- nrow(temp_2)

            n_cells = length_1 + length_2
            result_df <- data.frame(n_cells = n_cells)
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
        condition = "C9ALS" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        cell_type_list <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        result_list <- list()

        for(i in unique(cell_type_list)){
            file_1 <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/test_meta_",condition,"_",region,"_",i,"_pineda_narval.csv")
            temp_1 <- read.delim(file_1, header = T, sep = ",")
            length_1 <- nrow(temp_1)

            file_2 <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/train_meta_",condition,"_",region,"_",i,"_pineda_narval.csv")
            temp_2 <- read.delim(file_2, header = T, sep = ",")
            length_2 <- nrow(temp_2)

            n_cells = length_1 + length_2
                        result_df <- data.frame(n_cells = n_cells)
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
        condition = "SFTLD" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        cell_type_list <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        result_list <- list()

        for(i in unique(cell_type_list)){
            file_1 <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/test_meta_",condition,"_",region,"_",i,"_pineda_narval.csv")
            temp_1 <- read.delim(file_1, header = T, sep = ",")
            length_1 <- nrow(temp_1)

            file_2 <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/train_meta_",condition,"_",region,"_",i,"_pineda_narval.csv")
            temp_2 <- read.delim(file_2, header = T, sep = ",")
            length_2 <- nrow(temp_2)

            n_cells = length_1 + length_2
                        result_df <- data.frame(n_cells = n_cells)
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
        condition = "C9FTLD" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        cell_type_list <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        result_list <- list()

        for(i in unique(cell_type_list)){
            file_1 <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/test_meta_",condition,"_",region,"_",i,"_pineda_narval.csv")
            temp_1 <- read.delim(file_1, header = T, sep = ",")
            length_1 <- nrow(temp_1)

            file_2 <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/train_meta_",condition,"_",region,"_",i,"_pineda_narval.csv")
            temp_2 <- read.delim(file_2, header = T, sep = ",")
            length_2 <- nrow(temp_2)

            n_cells = length_1 + length_2
                        result_df <- data.frame(n_cells = n_cells)
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
    ## code
        BA4_bind_all <- rbind(final_df_SALS_BA4, final_df_C9ALS_BA4, final_df_SFTLD_BA4, final_df_C9FTLD_BA4)

        ## set factor levels
        BA4_bind_all$condition <- factor(BA4_bind_all$condition, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD"))
        BA4_bind_all$celltype <- factor(BA4_bind_all$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))


        ggplot(BA4_bind_all, aes(x = condition, y = n_cells, fill = condition, colour = condition)) + 
        theme_bw() + 
        geom_bar(stat = 'identity') +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~celltype, scales = "free_x") +
        ylab("n cells") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple"))
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA4_nCells.pdf', height = 2, width = 13)
    ##

    ###########################
    ## SALS BA9
    ###########################
    ## code
        ## define variables
        condition = "SALS" #one ofSALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        cell_type_list <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        result_list <- list()

        for(i in unique(cell_type_list)){
            file_1 <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/test_meta_",condition,"_",region,"_",i,"_pineda_narval.csv")
            temp_1 <- read.delim(file_1, header = T, sep = ",")
            length_1 <- nrow(temp_1)

            file_2 <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/train_meta_",condition,"_",region,"_",i,"_pineda_narval.csv")
            temp_2 <- read.delim(file_2, header = T, sep = ",")
            length_2 <- nrow(temp_2)

            n_cells = length_1 + length_2
            result_df <- data.frame(n_cells = n_cells)
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
        condition = "C9ALS" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        cell_type_list <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        result_list <- list()

        for(i in unique(cell_type_list)){
            file_1 <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/test_meta_",condition,"_",region,"_",i,"_pineda_narval.csv")
            temp_1 <- read.delim(file_1, header = T, sep = ",")
            length_1 <- nrow(temp_1)

            file_2 <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/train_meta_",condition,"_",region,"_",i,"_pineda_narval.csv")
            temp_2 <- read.delim(file_2, header = T, sep = ",")
            length_2 <- nrow(temp_2)

            n_cells = length_1 + length_2
                        result_df <- data.frame(n_cells = n_cells)
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
        condition = "SFTLD" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        cell_type_list <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        result_list <- list()

        for(i in unique(cell_type_list)){
            file_1 <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/test_meta_",condition,"_",region,"_",i,"_pineda_narval.csv")
            temp_1 <- read.delim(file_1, header = T, sep = ",")
            length_1 <- nrow(temp_1)

            file_2 <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/train_meta_",condition,"_",region,"_",i,"_pineda_narval.csv")
            temp_2 <- read.delim(file_2, header = T, sep = ",")
            length_2 <- nrow(temp_2)

            n_cells = length_1 + length_2
                        result_df <- data.frame(n_cells = n_cells)
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
        condition = "C9FTLD" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        cell_type_list <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        result_list <- list()

        for(i in unique(cell_type_list)){
            file_1 <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/test_meta_",condition,"_",region,"_",i,"_pineda_narval.csv")
            temp_1 <- read.delim(file_1, header = T, sep = ",")
            length_1 <- nrow(temp_1)

            file_2 <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/train_meta_",condition,"_",region,"_",i,"_pineda_narval.csv")
            temp_2 <- read.delim(file_2, header = T, sep = ",")
            length_2 <- nrow(temp_2)

            n_cells = length_1 + length_2
                        result_df <- data.frame(n_cells = n_cells)
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
    ## code
        BA9_bind_all <- rbind(final_df_SALS_BA9, final_df_C9ALS_BA9, final_df_SFTLD_BA9, final_df_C9FTLD_BA9)

        ## set factor levels
        BA9_bind_all$condition <- factor(BA9_bind_all$condition, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD"))
        BA9_bind_all$celltype <- factor(BA9_bind_all$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))


        ggplot(BA9_bind_all, aes(x = condition, y = n_cells, fill = condition, colour = condition)) + 
        theme_bw() + 
        geom_bar(stat = 'identity') +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~celltype, scales = "free_x") +
        ylab("n cells") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple"))
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA9_nCells.pdf', height = 2, width = 13)
    ##

##


################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ 5 CV workflow: DNN 5CV classification balanced accuracies (Supplemental)
## code
    
    ## BA4
    
    ###########################
    ## SALS BA4
    ###########################
    ## code
        ## define variables
        condition = "SALS" #one ofSALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        #cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')


        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            test_column <- c(rep("validation", length(accuracies)), "test")
            result_df <- data.frame(accuracy = all_accuracies, test = test_column)
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
        condition = "C9ALS" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        #cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            test_column <- c(rep("validation", length(accuracies)), "test")
            result_df <- data.frame(accuracy = all_accuracies, test = test_column)
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
        condition = "SFTLD" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        #cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            test_column <- c(rep("validation", length(accuracies)), "test")
            result_df <- data.frame(accuracy = all_accuracies, test = test_column)
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
        condition = "C9FTLD" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        #cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            test_column <- c(rep("validation", length(accuracies)), "test")
            result_df <- data.frame(accuracy = all_accuracies, test = test_column)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df

        }

        final_df_C9FTLD_BA4 <- do.call(rbind, result_list)
    ##

    ###########################
    ## all_conditions BA4
    ###########################
    ## code
        ## define variables
        condition = "all_conditions" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        #cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            test_column <- c(rep("validation", length(accuracies)), "test")
            result_df <- data.frame(accuracy = all_accuracies, test = test_column)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df

        }

        final_df_all_conditions_BA4 <- do.call(rbind, result_list)
    ##

    ###########################
    ## merge and plot BA4
    ###########################
    ## code for boxplot
        BA4_bind_all <- rbind(final_df_SALS_BA4, final_df_C9ALS_BA4, final_df_SFTLD_BA4, final_df_C9FTLD_BA4, final_df_all_conditions_BA4)
        BA4_bind_all$condition[BA4_bind_all$condition == "all_conditions"] <- "All conditions"

        ## Only retain the 5CV values, remove test (we will plot test in the main figure)
        BA4_bind_all <- subset(BA4_bind_all, test == "validation" )
        
        ## set factor levels
        BA4_bind_all$condition <- factor(BA4_bind_all$condition, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All conditions"))
        BA4_bind_all$celltype <- factor(BA4_bind_all$celltype, levels = c('L4_L6','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Astro', 'OPC','T_Cell'  ))


        ggplot(BA4_bind_all, aes(x = condition, y = accuracy, fill = condition, colour = condition)) + 
        theme_bw() + 
        geom_boxplot() +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~celltype, scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "black"))
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA4_5CV_balanced_acc_boxplot.pdf', height = 2, width = 13)
    ##

    ## code for barplot
        BA4_bind_all <- rbind(final_df_SALS_BA4, final_df_C9ALS_BA4, final_df_SFTLD_BA4, final_df_C9FTLD_BA4, final_df_all_conditions_BA4)
        BA4_bind_all$condition[BA4_bind_all$condition == "all_conditions"] <- "All conditions"

        ## Only retain the 5CV values, remove test (we will plot test in the main figure)
        BA4_bind_all <- subset(BA4_bind_all, test == "validation" )
        
        ## Calculate the mean and standard deviation
        summary_stats <- BA4_bind_all %>%
            group_by(celltype, condition) %>%
            summarise(
                median_accuracy = median(accuracy), 
                mean_accuracy = mean(accuracy),           # Calculate the median of X0
                sd_accuracy = sd(accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(accuracy - median(accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## set factor levels
        summary_stats$condition <- factor(summary_stats$condition, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All conditions"))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L4_L6','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Astro', 'OPC','T_Cell'  ))


        ## Plot
        ggplot(summary_stats, aes(x = condition, y = mean_accuracy, fill = condition, colour = condition)) + 
        geom_errorbar(aes(ymin=mean_accuracy - sd_accuracy, ymax=mean_accuracy + sd_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~celltype, scales = "free_x") +
        ylab("Balanced\naccuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "black")) #+
        #coord_cartesian(ylim=c(0.5,1.0)) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA4_5CV_balanced_acc_barplot.pdf', height = 2, width = 13)
    ##

    ## BA9
    
    ###########################
    ## SALS BA9
    ###########################
    ## code
        ## define variables
        condition = "SALS" #one ofSALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        #cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            test_column <- c(rep("validation", length(accuracies)), "test")
            result_df <- data.frame(accuracy = all_accuracies, test = test_column)
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
        condition = "C9ALS" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        #cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            test_column <- c(rep("validation", length(accuracies)), "test")
            result_df <- data.frame(accuracy = all_accuracies, test = test_column)
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
        condition = "SFTLD" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        #cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            test_column <- c(rep("validation", length(accuracies)), "test")
            result_df <- data.frame(accuracy = all_accuracies, test = test_column)
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
        condition = "C9FTLD" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        #cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            test_column <- c(rep("validation", length(accuracies)), "test")
            result_df <- data.frame(accuracy = all_accuracies, test = test_column)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df

        }

        final_df_C9FTLD_BA9 <- do.call(rbind, result_list)
    ##

    ###########################
    ## all_conditions BA9
    ###########################
    ## code
        ## define variables
        condition = "all_conditions" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        #cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')


        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            test_column <- c(rep("validation", length(accuracies)), "test")
            result_df <- data.frame(accuracy = all_accuracies, test = test_column)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df

        }

        final_df_all_conditions_BA9 <- do.call(rbind, result_list)
    ##

    ###########################
    ## merge and plot BA9
    ###########################
    ## code for boxplot
        BA9_bind_all <- rbind(final_df_SALS_BA9, final_df_C9ALS_BA9, final_df_SFTLD_BA9, final_df_C9FTLD_BA9, final_df_all_conditions_BA9)
        BA9_bind_all$condition[BA9_bind_all$condition == "all_conditions"] <- "All conditions"

        ## Only retain the 5CV values, remove test (we will plot test in the main figure)
        BA9_bind_all <- subset(BA9_bind_all, test == "validation" )
        
        ## set factor levels
        BA9_bind_all$condition <- factor(BA9_bind_all$condition, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All conditions"))
        BA9_bind_all$celltype <- factor(BA9_bind_all$celltype, levels = c('L4_L6','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Astro', 'OPC','T_Cell'  ))


        ggplot(BA9_bind_all, aes(x = condition, y = accuracy, fill = condition, colour = condition)) + 
        theme_bw() + 
        geom_boxplot() +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~celltype, scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "black"))
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA9_5CV_balanced_acc_boxplot.pdf', height = 2, width = 13)
    ##

    ## code for barplot
        BA9_bind_all <- rbind(final_df_SALS_BA9, final_df_C9ALS_BA9, final_df_SFTLD_BA9, final_df_C9FTLD_BA9, final_df_all_conditions_BA9)
        BA9_bind_all$condition[BA9_bind_all$condition == "all_conditions"] <- "All conditions"

        ## Only retain the 5CV values, remove test (we will plot test in the main figure)
        BA9_bind_all <- subset(BA9_bind_all, test == "validation" )
        
        ## Calculate the mean and standard deviation
        summary_stats <- BA9_bind_all %>%
            group_by(celltype, condition) %>%
            summarise(
                median_accuracy = median(accuracy), 
                mean_accuracy = mean(accuracy),           # Calculate the median of X0
                sd_accuracy = sd(accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(accuracy - median(accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## set factor levels
        summary_stats$condition <- factor(summary_stats$condition, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All conditions"))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L4_L6','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Astro', 'OPC','T_Cell'  ))


        ## Plot
        ggplot(summary_stats, aes(x = condition, y = mean_accuracy, fill = condition, colour = condition)) + 
        geom_errorbar(aes(ymin=mean_accuracy - sd_accuracy, ymax=mean_accuracy + sd_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~celltype, scales = "free_x") +
        ylab("Balanced\naccuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "black")) #+
        #coord_cartesian(ylim=c(0.5,1.0)) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA9_5CV_balanced_acc_barplot.pdf', height = 2, width = 13)
    ##
##


################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ 5 CV workflow: DNN LOO classification balanced accuracies (Supplemental)
## code
    
    ## BA4
    
    ###########################
    ## SALS BA4
    ###########################
    ## code
        # List of cell types
        #cell_type <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')
        cell_type <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')

        status <- 'SALS'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2_LOO.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_SALS_BA4 <- do.call(rbind, df_list)
    ##

    ###########################
    ## C9ALS BA4
    ###########################
    ## code
        # List of cell types
        #cell_type <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')
        cell_type <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')
        status <- 'C9ALS'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2_LOO.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_C9ALS_BA4 <- do.call(rbind, df_list)
    ##

    ###########################
    ## SFTLD BA4
    ###########################
    ## code
        # List of cell types
        #cell_type <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')
        cell_type <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')
        
        status <- 'SFTLD'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2_LOO.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_SFTLD_BA4 <- do.call(rbind, df_list)
    ##

    ###########################
    ## C9FTLD BA4
    ###########################
    ## code
        # List of cell types
        #cell_type <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')
        cell_type <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')
        
        status <- 'C9FTLD'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2_LOO.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_C9FTLD_BA4 <- do.call(rbind, df_list)
    ##

    ###########################
    ## all_conditions BA4
    ###########################
    ## code
        # List of cell types
        #cell_type <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')
        cell_type <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')

        status <- 'all_conditions'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2_LOO.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_all_conditions_BA4 <- do.call(rbind, df_list)
    ##

    ###########################
    ## merge and plot BA4
    ###########################
    ## code for barplot
        BA4_bind_all <- rbind(final_SALS_BA4, final_C9ALS_BA4, final_SFTLD_BA4, final_C9FTLD_BA4, final_all_conditions_BA4)
        BA4_bind_all$group[BA4_bind_all$group == "all_conditions"] <- "All conditions"
        
        ## Calculate the mean and standard deviation
        summary_stats <- BA4_bind_all %>%
            group_by(celltype, group) %>%
            summarise(
                median_accuracy = median(test_accuracy), 
                mean_accuracy = mean(test_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(test_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(test_accuracy - median(test_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## set factor levels
        summary_stats$group <- factor(summary_stats$group, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All conditions"))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L4_L6','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Astro', 'OPC','T_Cell'  ))


        ## Plot
        ggplot(summary_stats, aes(x = group, y = median_accuracy, fill = group, colour = group)) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~celltype, scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "black")) #+
        #coord_cartesian(ylim=c(0.5,1.0)) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA4_LOO_balanced_acc_barplot.pdf', height = 2, width = 13)
    ##

    ## BA9
    
    ###########################
    ## SALS BA9
    ###########################
    ## code
        # List of cell types
        #cell_type <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')
        cell_type <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')

        status <- 'SALS'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2_LOO.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_SALS_BA9 <- do.call(rbind, df_list)
    ##

    ###########################
    ## C9ALS BA9
    ###########################
    ## code
        # List of cell types
        #cell_type <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')
        cell_type <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')
        status <- 'C9ALS'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2_LOO.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_C9ALS_BA9 <- do.call(rbind, df_list)
    ##

    ###########################
    ## SFTLD BA9
    ###########################
    ## code
        # List of cell types
        #cell_type <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')
        cell_type <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')
        
        status <- 'SFTLD'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2_LOO.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_SFTLD_BA9 <- do.call(rbind, df_list)
    ##

    ###########################
    ## C9FTLD BA9
    ###########################
    ## code
        # List of cell types
        #cell_type <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')
        cell_type <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')
        
        status <- 'C9FTLD'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2_LOO.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_C9FTLD_BA9 <- do.call(rbind, df_list)
    ##

    ###########################
    ## all_conditions BA9
    ###########################
    ## code
        # List of cell types
        #cell_type <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')
        cell_type <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')

        status <- 'all_conditions'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2_LOO.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_all_conditions_BA9 <- do.call(rbind, df_list)
    ##

    ###########################
    ## merge and plot BA9
    ###########################
    ## code for barplot
        BA9_bind_all <- rbind(final_SALS_BA9, final_C9ALS_BA9, final_SFTLD_BA9, final_C9FTLD_BA9, final_all_conditions_BA9)
        BA9_bind_all$group[BA9_bind_all$group == "all_conditions"] <- "All conditions"
        
        ## Calculate the mean and standard deviation
        summary_stats <- BA9_bind_all %>%
            group_by(celltype, group) %>%
            summarise(
                median_accuracy = median(test_accuracy), 
                mean_accuracy = mean(test_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(test_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(test_accuracy - median(test_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## set factor levels
        summary_stats$group <- factor(summary_stats$group, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All conditions"))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L4_L6','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Astro', 'OPC','T_Cell'  ))


        ## Plot
        ggplot(summary_stats, aes(x = group, y = median_accuracy, fill = group, colour = group)) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~celltype, scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "black")) #+
        #coord_cartesian(ylim=c(0.5,1.0)) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA9_LOO_balanced_acc_barplot.pdf', height = 2, width = 13)
    ##
##


################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ Optimal: LOO workflow 0.9995 delta threshold: Number of sig genes -- not using this
## code
    
    #BA4
    
    ###########################
    ## SALS BA4
    ###########################
    ## code
        ## define variables
        par_brain_region = "BA4"
        par_status = "SALS"
        par_prep = "CombatSeq"
            
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
            "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
            "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        df_list <- list()
        
        for (celltype2 in celltype_list) {
        file_path <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_",
                            par_prep, "_", par_status, "_", par_brain_region, "_", 
                            celltype2, "_.9995_abs_case_control_narval_2.csv")

        temp_df <- read_csv(file_path, col_types = cols(.default = col_character())) %>%
            mutate(celltype = celltype2)
        n_genes = length((unique(temp_df$gene)))
        HVGs <- as.numeric(n_genes)
        result_df <- data.frame(HVGs = HVGs)
        result_df$celltype <- celltype2
        result_df$condition <- par_status
        result_df$region <- par_brain_region
        result_list[[celltype2]] <- result_df
        }

        final_df_SALS_BA4 <- do.call(rbind, result_list) 
    ##

    ###########################
    ## C9ALS BA4
    ###########################
    ## code
        ## define variables
        par_brain_region = "BA4"
        par_status = "C9ALS"
        par_prep = "CombatSeq"
            
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
            "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
            "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        df_list <- list()
        
        for (celltype2 in celltype_list) {
        file_path <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_",
                            par_prep, "_", par_status, "_", par_brain_region, "_", 
                            celltype2, "_.9995_abs_case_control_narval_2.csv")

        temp_df <- read_csv(file_path, col_types = cols(.default = col_character())) %>%
            mutate(celltype = celltype2)
        n_genes = length((unique(temp_df$gene)))
        HVGs <- as.numeric(n_genes)
        result_df <- data.frame(HVGs = HVGs)
        result_df$celltype <- celltype2
        result_df$condition <- par_status
        result_df$region <- par_brain_region
        result_list[[celltype2]] <- result_df
        }

        final_df_C9ALS_BA4 <- do.call(rbind, result_list) 
    ##

    ###########################
    ## SFTLD BA4
    ###########################
    ## code
        ## define variables
        par_brain_region = "BA4"
        par_status = "SFTLD"
        par_prep = "CombatSeq"
            
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
            "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
            "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        df_list <- list()
        
        for (celltype2 in celltype_list) {
        file_path <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_",
                            par_prep, "_", par_status, "_", par_brain_region, "_", 
                            celltype2, "_.9995_abs_case_control_narval_2.csv")

        temp_df <- read_csv(file_path, col_types = cols(.default = col_character())) %>%
            mutate(celltype = celltype2)
        n_genes = length((unique(temp_df$gene)))
        HVGs <- as.numeric(n_genes)
        result_df <- data.frame(HVGs = HVGs)
        result_df$celltype <- celltype2
        result_df$condition <- par_status
        result_df$region <- par_brain_region
        result_list[[celltype2]] <- result_df
        }

        final_df_SFTLD_BA4 <- do.call(rbind, result_list) 
    ##

    ###########################
    ## C9FTLD BA4
    ###########################
    ## code
        ## define variables
        par_brain_region = "BA4"
        par_status = "C9FTLD"
        par_prep = "CombatSeq"
            
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
            "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
            "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        df_list <- list()
        
        for (celltype2 in celltype_list) {
        file_path <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_",
                            par_prep, "_", par_status, "_", par_brain_region, "_", 
                            celltype2, "_.9995_abs_case_control_narval_2.csv")

        temp_df <- read_csv(file_path, col_types = cols(.default = col_character())) %>%
            mutate(celltype = celltype2)
        n_genes = length((unique(temp_df$gene)))
        HVGs <- as.numeric(n_genes)
        result_df <- data.frame(HVGs = HVGs)
        result_df$celltype <- celltype2
        result_df$condition <- par_status
        result_df$region <- par_brain_region
        result_list[[celltype2]] <- result_df
        }

        final_df_C9FTLD_BA4 <- do.call(rbind, result_list) 
    ##

    
    ###########################
    ## merge and plot BA4
    ###########################
    ## code
        BA4_bind_all <- rbind(final_df_SALS_BA4, final_df_C9ALS_BA4, final_df_SFTLD_BA4, final_df_C9FTLD_BA4)

        ## set factor levels
        BA4_bind_all$condition <- factor(BA4_bind_all$condition, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD"))
        BA4_bind_all$celltype <- factor(BA4_bind_all$celltype, levels = c('L4_L6','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Astro', 'OPC'))


        ggplot(BA4_bind_all, aes(x = condition, y = HVGs, fill = condition, colour = condition)) + 
        theme_bw() + 
        geom_bar(stat = 'identity') +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~celltype, scales = "free_x") +
        ylab("n HVGs") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "black"))
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA4_n_sig_genes_0.9995.pdf', height = 2, width = 13)
    ##


    ## BA9

    ###########################
    ## SALS BA9
    ###########################
    ## code
        ## define variables
        par_brain_region = "BA9"
        par_status = "SALS"
        par_prep = "CombatSeq"
            
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
            "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
            "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        df_list <- list()
        
        for (celltype2 in celltype_list) {
        file_path <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_",
                            par_prep, "_", par_status, "_", par_brain_region, "_", 
                            celltype2, "_.9995_abs_case_control_narval_2.csv")

        temp_df <- read_csv(file_path, col_types = cols(.default = col_character())) %>%
            mutate(celltype = celltype2)
        n_genes = length((unique(temp_df$gene)))
        HVGs <- as.numeric(n_genes)
        result_df <- data.frame(HVGs = HVGs)
        result_df$celltype <- celltype2
        result_df$condition <- par_status
        result_df$region <- par_brain_region
        result_list[[celltype2]] <- result_df
        }

        final_df_SALS_BA9 <- do.call(rbind, result_list) 
    ##

    ###########################
    ## C9ALS BA9
    ###########################
    ## code
        ## define variables
        par_brain_region = "BA9"
        par_status = "C9ALS"
        par_prep = "CombatSeq"
            
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
            "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
            "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        df_list <- list()
        
        for (celltype2 in celltype_list) {
        file_path <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_",
                            par_prep, "_", par_status, "_", par_brain_region, "_", 
                            celltype2, "_.9995_abs_case_control_narval_2.csv")

        temp_df <- read_csv(file_path, col_types = cols(.default = col_character())) %>%
            mutate(celltype = celltype2)
        n_genes = length((unique(temp_df$gene)))
        HVGs <- as.numeric(n_genes)
        result_df <- data.frame(HVGs = HVGs)
        result_df$celltype <- celltype2
        result_df$condition <- par_status
        result_df$region <- par_brain_region
        result_list[[celltype2]] <- result_df
        }

        final_df_C9ALS_BA9 <- do.call(rbind, result_list) 
    ##

    ###########################
    ## SFTLD BA9
    ###########################
    ## code
        ## define variables
        par_brain_region = "BA9"
        par_status = "SFTLD"
        par_prep = "CombatSeq"
            
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
            "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
            "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        df_list <- list()
        
        for (celltype2 in celltype_list) {
        file_path <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_",
                            par_prep, "_", par_status, "_", par_brain_region, "_", 
                            celltype2, "_.9995_abs_case_control_narval_2.csv")

        temp_df <- read_csv(file_path, col_types = cols(.default = col_character())) %>%
            mutate(celltype = celltype2)
        n_genes = length((unique(temp_df$gene)))
        HVGs <- as.numeric(n_genes)
        result_df <- data.frame(HVGs = HVGs)
        result_df$celltype <- celltype2
        result_df$condition <- par_status
        result_df$region <- par_brain_region
        result_list[[celltype2]] <- result_df
        }

        final_df_SFTLD_BA9 <- do.call(rbind, result_list) 
    ##

    ###########################
    ## C9FTLD BA9
    ###########################
    ## code
        ## define variables
        par_brain_region = "BA9"
        par_status = "C9FTLD"
        par_prep = "CombatSeq"
            
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
            "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
            "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        df_list <- list()
        
        for (celltype2 in celltype_list) {
        file_path <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_",
                            par_prep, "_", par_status, "_", par_brain_region, "_", 
                            celltype2, "_.9995_abs_case_control_narval_2.csv")

        temp_df <- read_csv(file_path, col_types = cols(.default = col_character())) %>%
            mutate(celltype = celltype2)
        n_genes = length((unique(temp_df$gene)))
        HVGs <- as.numeric(n_genes)
        result_df <- data.frame(HVGs = HVGs)
        result_df$celltype <- celltype2
        result_df$condition <- par_status
        result_df$region <- par_brain_region
        result_list[[celltype2]] <- result_df
        }

        final_df_C9FTLD_BA9 <- do.call(rbind, result_list) 
    ##

    
    ###########################
    ## merge and plot BA9
    ###########################
    ## code
        BA9_bind_all <- rbind(final_df_SALS_BA9, final_df_C9ALS_BA9, final_df_SFTLD_BA9, final_df_C9FTLD_BA9)

        ## set factor levels
        BA9_bind_all$condition <- factor(BA9_bind_all$condition, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD"))
        BA9_bind_all$celltype <- factor(BA9_bind_all$celltype, levels = c('L4_L6','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Astro', 'OPC'))


        ggplot(BA9_bind_all, aes(x = condition, y = HVGs, fill = condition, colour = condition)) + 
        theme_bw() + 
        geom_bar(stat = 'identity') +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~celltype, scales = "free_x") +
        ylab("n HVGs") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "black"))
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA9_n_sig_genes_0.9995.pdf', height = 2, width = 13)
    ##





##

################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ Optimal: LOO workflow 0 delta threshold: Number of sig genes
## code
    
    #BA4
    
    ###########################
    ## SALS BA4
    ###########################
    ## code
        ## define variables
        par_brain_region = "BA4"
        par_status = "SALS"
        par_prep = "CombatSeq"
            
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
            "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
            "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        df_list <- list()
        
        for (celltype2 in celltype_list) {
        file_path <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_",
                            par_prep, "_", par_status, "_", par_brain_region, "_", 
                            celltype2, "_0_abs_case_control_narval_2.csv")

        temp_df <- read_csv(file_path, col_types = cols(.default = col_character())) %>%
            mutate(celltype = celltype2)
        n_genes = length((unique(temp_df$gene)))
        HVGs <- as.numeric(n_genes)
        result_df <- data.frame(HVGs = HVGs)
        result_df$celltype <- celltype2
        result_df$condition <- par_status
        result_df$region <- par_brain_region
        result_list[[celltype2]] <- result_df
        }

        final_df_SALS_BA4 <- do.call(rbind, result_list) 
    ##

    ###########################
    ## C9ALS BA4
    ###########################
    ## code
        ## define variables
        par_brain_region = "BA4"
        par_status = "C9ALS"
        par_prep = "CombatSeq"
            
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
            "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
            "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        df_list <- list()
        
        for (celltype2 in celltype_list) {
        file_path <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_",
                            par_prep, "_", par_status, "_", par_brain_region, "_", 
                            celltype2, "_0_abs_case_control_narval_2.csv")

        temp_df <- read_csv(file_path, col_types = cols(.default = col_character())) %>%
            mutate(celltype = celltype2)
        n_genes = length((unique(temp_df$gene)))
        HVGs <- as.numeric(n_genes)
        result_df <- data.frame(HVGs = HVGs)
        result_df$celltype <- celltype2
        result_df$condition <- par_status
        result_df$region <- par_brain_region
        result_list[[celltype2]] <- result_df
        }

        final_df_C9ALS_BA4 <- do.call(rbind, result_list) 
    ##

    ###########################
    ## SFTLD BA4
    ###########################
    ## code
        ## define variables
        par_brain_region = "BA4"
        par_status = "SFTLD"
        par_prep = "CombatSeq"
            
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
            "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
            "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        df_list <- list()
        
        for (celltype2 in celltype_list) {
        file_path <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_",
                            par_prep, "_", par_status, "_", par_brain_region, "_", 
                            celltype2, "_0_abs_case_control_narval_2.csv")

        temp_df <- read_csv(file_path, col_types = cols(.default = col_character())) %>%
            mutate(celltype = celltype2)
        n_genes = length((unique(temp_df$gene)))
        HVGs <- as.numeric(n_genes)
        result_df <- data.frame(HVGs = HVGs)
        result_df$celltype <- celltype2
        result_df$condition <- par_status
        result_df$region <- par_brain_region
        result_list[[celltype2]] <- result_df
        }

        final_df_SFTLD_BA4 <- do.call(rbind, result_list) 
    ##

    ###########################
    ## C9FTLD BA4
    ###########################
    ## code
        ## define variables
        par_brain_region = "BA4"
        par_status = "C9FTLD"
        par_prep = "CombatSeq"
            
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
            "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
            "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        df_list <- list()
        
        for (celltype2 in celltype_list) {
        file_path <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_",
                            par_prep, "_", par_status, "_", par_brain_region, "_", 
                            celltype2, "_0_abs_case_control_narval_2.csv")

        temp_df <- read_csv(file_path, col_types = cols(.default = col_character())) %>%
            mutate(celltype = celltype2)
        n_genes = length((unique(temp_df$gene)))
        HVGs <- as.numeric(n_genes)
        result_df <- data.frame(HVGs = HVGs)
        result_df$celltype <- celltype2
        result_df$condition <- par_status
        result_df$region <- par_brain_region
        result_list[[celltype2]] <- result_df
        }

        final_df_C9FTLD_BA4 <- do.call(rbind, result_list) 
    ##

    
    ###########################
    ## merge and plot BA4
    ###########################
    ## code
        BA4_bind_all <- rbind(final_df_SALS_BA4, final_df_C9ALS_BA4, final_df_SFTLD_BA4, final_df_C9FTLD_BA4)

        ## set factor levels
        BA4_bind_all$condition <- factor(BA4_bind_all$condition, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD"))
        BA4_bind_all$celltype <- factor(BA4_bind_all$celltype, levels = c('L4_L6','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Astro', 'OPC'))


        ggplot(BA4_bind_all, aes(x = condition, y = HVGs, fill = condition, colour = condition)) + 
        theme_bw() + 
        geom_bar(stat = 'identity') +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~celltype, scales = "free_x") +
        ylab("n HVGs") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "black"))
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA4_n_sig_genes_0.pdf', height = 2, width = 13)
    ##


    ## BA9

    ###########################
    ## SALS BA9
    ###########################
    ## code
        ## define variables
        par_brain_region = "BA9"
        par_status = "SALS"
        par_prep = "CombatSeq"
            
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
            "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
            "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        df_list <- list()
        
        for (celltype2 in celltype_list) {
        file_path <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_",
                            par_prep, "_", par_status, "_", par_brain_region, "_", 
                            celltype2, "_0_abs_case_control_narval_2.csv")

        temp_df <- read_csv(file_path, col_types = cols(.default = col_character())) %>%
            mutate(celltype = celltype2)
        n_genes = length((unique(temp_df$gene)))
        HVGs <- as.numeric(n_genes)
        result_df <- data.frame(HVGs = HVGs)
        result_df$celltype <- celltype2
        result_df$condition <- par_status
        result_df$region <- par_brain_region
        result_list[[celltype2]] <- result_df
        }

        final_df_SALS_BA9 <- do.call(rbind, result_list) 
    ##

    ###########################
    ## C9ALS BA9
    ###########################
    ## code
        ## define variables
        par_brain_region = "BA9"
        par_status = "C9ALS"
        par_prep = "CombatSeq"
            
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
            "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
            "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        df_list <- list()
        
        for (celltype2 in celltype_list) {
        file_path <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_",
                            par_prep, "_", par_status, "_", par_brain_region, "_", 
                            celltype2, "_0_abs_case_control_narval_2.csv")

        temp_df <- read_csv(file_path, col_types = cols(.default = col_character())) %>%
            mutate(celltype = celltype2)
        n_genes = length((unique(temp_df$gene)))
        HVGs <- as.numeric(n_genes)
        result_df <- data.frame(HVGs = HVGs)
        result_df$celltype <- celltype2
        result_df$condition <- par_status
        result_df$region <- par_brain_region
        result_list[[celltype2]] <- result_df
        }

        final_df_C9ALS_BA9 <- do.call(rbind, result_list) 
    ##

    ###########################
    ## SFTLD BA9
    ###########################
    ## code
        ## define variables
        par_brain_region = "BA9"
        par_status = "SFTLD"
        par_prep = "CombatSeq"
            
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
            "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
            "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        df_list <- list()
        
        for (celltype2 in celltype_list) {
        file_path <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_",
                            par_prep, "_", par_status, "_", par_brain_region, "_", 
                            celltype2, "_0_abs_case_control_narval_2.csv")

        temp_df <- read_csv(file_path, col_types = cols(.default = col_character())) %>%
            mutate(celltype = celltype2)
        n_genes = length((unique(temp_df$gene)))
        HVGs <- as.numeric(n_genes)
        result_df <- data.frame(HVGs = HVGs)
        result_df$celltype <- celltype2
        result_df$condition <- par_status
        result_df$region <- par_brain_region
        result_list[[celltype2]] <- result_df
        }

        final_df_SFTLD_BA9 <- do.call(rbind, result_list) 
    ##

    ###########################
    ## C9FTLD BA9
    ###########################
    ## code
        ## define variables
        par_brain_region = "BA9"
        par_status = "C9FTLD"
        par_prep = "CombatSeq"
            
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
            "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
            "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        df_list <- list()
        
        for (celltype2 in celltype_list) {
        file_path <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_",
                            par_prep, "_", par_status, "_", par_brain_region, "_", 
                            celltype2, "_0_abs_case_control_narval_2.csv")

        temp_df <- read_csv(file_path, col_types = cols(.default = col_character())) %>%
            mutate(celltype = celltype2)
        n_genes = length((unique(temp_df$gene)))
        HVGs <- as.numeric(n_genes)
        result_df <- data.frame(HVGs = HVGs)
        result_df$celltype <- celltype2
        result_df$condition <- par_status
        result_df$region <- par_brain_region
        result_list[[celltype2]] <- result_df
        }

        final_df_C9FTLD_BA9 <- do.call(rbind, result_list) 
    ##

    
    ###########################
    ## merge and plot BA9
    ###########################
    ## code
        BA9_bind_all <- rbind(final_df_SALS_BA9, final_df_C9ALS_BA9, final_df_SFTLD_BA9, final_df_C9FTLD_BA9)

        ## set factor levels
        BA9_bind_all$condition <- factor(BA9_bind_all$condition, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD"))
        BA9_bind_all$celltype <- factor(BA9_bind_all$celltype, levels = c('L4_L6','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Astro', 'OPC'))


        ggplot(BA9_bind_all, aes(x = condition, y = HVGs, fill = condition, colour = condition)) + 
        theme_bw() + 
        geom_bar(stat = 'identity') +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~celltype, scales = "free_x") +
        ylab("n HVGs") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "black"))
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA9_n_sig_genes_0.pdf', height = 2, width = 13)
    ##


##

################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ Optimal: LOO workflow 0.9995 delta threshold: Accuracy -- not using this
## code

    ## BA4
    
    ###########################
    ## SALS BA4
    ###########################
    ## code
        # List of cell types
        cell_type <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')
        status <- 'SALS'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_LOO_LIME_0.9995_abs_case_control_report_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_SALS_BA4 <- do.call(rbind, df_list)
    ##

    ###########################
    ## C9ALS BA4
    ###########################
    ## code
        # List of cell types
        cell_type <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')
        status <- 'C9ALS'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_LOO_LIME_0.9995_abs_case_control_report_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_C9ALS_BA4 <- do.call(rbind, df_list)
    ##

    ###########################
    ## SFTLD BA4
    ###########################
    ## code
        # List of cell types
        cell_type <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')
        status <- 'SFTLD'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_LOO_LIME_0.9995_abs_case_control_report_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_SFTLD_BA4 <- do.call(rbind, df_list)
    ##

    ###########################
    ## C9FTLD BA4
    ###########################
    ## code
        # List of cell types
        cell_type <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')
        status <- 'C9FTLD'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_LOO_LIME_0.9995_abs_case_control_report_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_C9FTLD_BA4 <- do.call(rbind, df_list)
    ##

    ###########################
    ## merge and plot BA4
    ###########################
    ## code for barplot
        BA4_bind_all <- rbind(final_SALS_BA4, final_C9ALS_BA4, final_SFTLD_BA4, final_C9FTLD_BA4)
        
        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, group) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy),        
                .groups = "drop"               
            )

        ## Calculate the mean and standard deviation
        summary_stats <- BA4_bind_all %>%
            group_by(celltype, group) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## set factor levels
        summary_stats$group <- factor(summary_stats$group, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD"))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L4_L6','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Astro', 'OPC','T_Cell'  ))


        ## Plot
        ggplot(summary_stats, aes(x = group, y = median_accuracy, fill = group, colour = group)) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~celltype, scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "black")) #+
        #coord_cartesian(ylim=c(0.5,1.0)) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA4_LOO_balanced_acc_barplot_0.9995.pdf', height = 2, width = 13)
    ##

    ## BA9
    
    ###########################
    ## SALS BA9
    ###########################
    ## code
        # List of cell types
        cell_type <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')
        status <- 'SALS'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_LOO_LIME_0.9995_abs_case_control_report_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_SALS_BA9 <- do.call(rbind, df_list)
    ##

    ###########################
    ## C9ALS BA9
    ###########################
    ## code
        # List of cell types
        cell_type <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')
        status <- 'C9ALS'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_LOO_LIME_0.9995_abs_case_control_report_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_C9ALS_BA9 <- do.call(rbind, df_list)
    ##

    ###########################
    ## SFTLD BA9
    ###########################
    ## code
        # List of cell types
        cell_type <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')
        status <- 'SFTLD'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_LOO_LIME_0.9995_abs_case_control_report_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_SFTLD_BA9 <- do.call(rbind, df_list)
    ##

    ###########################
    ## C9FTLD BA9
    ###########################
    ## code
        # List of cell types
        cell_type <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')
        status <- 'C9FTLD'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_LOO_LIME_0.9995_abs_case_control_report_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_C9FTLD_BA9 <- do.call(rbind, df_list)
    ##

    ###########################
    ## merge and plot BA9
    ###########################
    ## code for barplot
        BA9_bind_all <- rbind(final_SALS_BA9, final_C9ALS_BA9, final_SFTLD_BA9, final_C9FTLD_BA9)
        
        ## Calculate sample wise mean
        BA9_bind_all <- BA9_bind_all %>%
            group_by(donor, celltype, group) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy),        
                .groups = "drop"               
            )

        ## Calculate the mean and standard deviation
        summary_stats <- BA9_bind_all %>%
            group_by(celltype, group) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## set factor levels
        summary_stats$group <- factor(summary_stats$group, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD"))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L4_L6','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Astro', 'OPC','T_Cell'  ))


        ## Plot
        ggplot(summary_stats, aes(x = group, y = median_accuracy, fill = group, colour = group)) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~celltype, scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "black")) #+
        #coord_cartesian(ylim=c(0.5,1.0)) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA9_LOO_balanced_acc_barplot_0.9995.pdf', height = 2, width = 13)
    ##
##



################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ Optimal: LOO workflow 0 delta threshold: per sample accuracy visualization
## code
    ## BA4
    
    ###########################
    ## SALS BA4
    ###########################
    ## code
        # List of cell types
        cell_type <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')
        status <- 'SALS'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_SALS_BA4 <- do.call(rbind, df_list)

        ## Calculate sample wise mean
        final_SALS_BA4 <- final_SALS_BA4 %>%
            group_by(donor, celltype, group, n_cells) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy),        
                .groups = "drop"               
            )

        final_SALS_BA4$label <- paste0(round(final_SALS_BA4$sample_mean_accuracy, 2),'\n(',final_SALS_BA4$n_cells, ')' )

        final_SALS_BA4$celltype <- factor(final_SALS_BA4$celltype, levels = c('L4_L6','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Astro', 'OPC','T_Cell'  ))

        final_SALS_BA4 <- subset(final_SALS_BA4, n_cells > 10)
        
        ## Plot
        BA4_SALS_heat <- ggplot(final_SALS_BA4, aes(x = celltype, y = donor, fill = sample_mean_accuracy, label = label)) + 
        geom_tile() +
        geom_text(size = 1.75) +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title = element_blank()
        ) +
        scale_fill_gradientn(colors = brewer.pal(9, "RdYlGn"))   # Use the "RdYlBu" palette for continuous data
        #coord_cartesian(ylim=c(0.5,1.0)) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 8, width = 7)
        


    ##

    ###########################
    ## C9ALS BA4
    ###########################
    ## code
        # List of cell types
        cell_type <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')
        status <- 'C9ALS'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_C9ALS_BA4 <- do.call(rbind, df_list)

        ## Calculate sample wise mean
        final_C9ALS_BA4 <- final_C9ALS_BA4 %>%
            group_by(donor, celltype, group, n_cells) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy),        
                .groups = "drop"               
            )

        final_C9ALS_BA4$label <- paste0(round(final_C9ALS_BA4$sample_mean_accuracy, 2),'\n(',final_C9ALS_BA4$n_cells, ')' )

        final_C9ALS_BA4$celltype <- factor(final_C9ALS_BA4$celltype, levels = c('L4_L6','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Astro', 'OPC','T_Cell'  ))

        #final_C9ALS_BA4 <- subset(final_C9ALS_BA4, n_cells > 10)
        
        ## Plot
        BA4_C9ALS_heat <- ggplot(final_C9ALS_BA4, aes(x = celltype, y = donor, fill = sample_mean_accuracy, label = label)) + 
        geom_tile() +
        geom_text(size = 1.75) +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title = element_blank()
        ) +
        scale_fill_gradientn(colors = brewer.pal(9, "RdYlGn"))   # Use the "RdYlBu" palette for continuous data
        #coord_cartesian(ylim=c(0.5,1.0)) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 8, width = 7)
    ##

    ###########################
    ## SFTLD BA4
    ###########################
    ## code
        # List of cell types
        cell_type <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')
        status <- 'SFTLD'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_SFTLD_BA4 <- do.call(rbind, df_list)

        ## Calculate sample wise mean
        final_SFTLD_BA4 <- final_SFTLD_BA4 %>%
            group_by(donor, celltype, group, n_cells) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy),        
                .groups = "drop"               
            )

        final_SFTLD_BA4$label <- paste0(round(final_SFTLD_BA4$sample_mean_accuracy, 2),'\n(',final_SFTLD_BA4$n_cells, ')' )

        final_SFTLD_BA4$celltype <- factor(final_SFTLD_BA4$celltype, levels = c('L4_L6','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Astro', 'OPC','T_Cell'  ))

        #final_SFTLD_BA4 <- subset(final_SFTLD_BA4, n_cells > 10)
        
        ## Plot
        BA4_SFTLD_heat <- ggplot(final_SFTLD_BA4, aes(x = celltype, y = donor, fill = sample_mean_accuracy, label = label)) + 
        geom_tile() +
        geom_text(size = 1.75) +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title = element_blank()
        ) +
        scale_fill_gradientn(colors = brewer.pal(9, "RdYlGn"))   # Use the "RdYlBu" palette for continuous data
        #coord_cartesian(ylim=c(0.5,1.0)) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 8, width = 7)
        
    ##

    ###########################
    ## C9FTLD BA4
    ###########################
    ## code
        # List of cell types
        cell_type <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')
        status <- 'SFTLD'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_SFTLD_BA4 <- do.call(rbind, df_list)

        ## Calculate sample wise mean
        final_SFTLD_BA4 <- final_SFTLD_BA4 %>%
            group_by(donor, celltype, group, n_cells) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy),        
                .groups = "drop"               
            )

        final_SFTLD_BA4$label <- paste0(round(final_SFTLD_BA4$sample_mean_accuracy, 2),'\n(',final_SFTLD_BA4$n_cells, ')' )

        final_SFTLD_BA4$celltype <- factor(final_SFTLD_BA4$celltype, levels = c('L4_L6','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Astro', 'OPC','T_Cell'  ))

        #final_SFTLD_BA4 <- subset(final_SFTLD_BA4, n_cells > 10)
        
        ## Plot
        BA4_SFTLD_heat <- ggplot(final_SFTLD_BA4, aes(x = celltype, y = donor, fill = sample_mean_accuracy, label = label)) + 
        geom_tile() +
        geom_text(size = 1.75) +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title = element_blank()
        ) +
        scale_fill_gradientn(colors = brewer.pal(9, "RdYlGn"))   # Use the "RdYlBu" palette for continuous data
        #coord_cartesian(ylim=c(0.5,1.0)) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 8, width = 7)
        
    ##

    ## BA9
    
    ###########################
    ## SALS BA9
    ###########################
    ## code
        # List of cell types
        cell_type <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')
        status <- 'SALS'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_SALS_BA9 <- do.call(rbind, df_list)

        ## Calculate sample wise mean
        final_SALS_BA9 <- final_SALS_BA9 %>%
            group_by(donor, celltype, group, n_cells) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy),        
                .groups = "drop"               
            )

        final_SALS_BA9$label <- paste0(round(final_SALS_BA9$sample_mean_accuracy, 2),'\n(',final_SALS_BA9$n_cells, ')' )

        final_SALS_BA9$celltype <- factor(final_SALS_BA9$celltype, levels = c('L4_L6','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Astro', 'OPC','T_Cell'  ))

        final_SALS_BA9 <- subset(final_SALS_BA9, n_cells > 10)
        
        ## Plot
        BA9_SALS_heat <- ggplot(final_SALS_BA9, aes(x = celltype, y = donor, fill = sample_mean_accuracy, label = label)) + 
        geom_tile() +
        geom_text(size = 1.75) +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title = element_blank()
        ) +
        scale_fill_gradientn(colors = brewer.pal(9, "RdYlGn"))   # Use the "RdYlBu" palette for continuous data
        #coord_cartesian(ylim=c(0.5,1.0)) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 8, width = 7)
        


    ##

    ###########################
    ## C9ALS BA9
    ###########################
    ## code
        # List of cell types
        cell_type <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')
        status <- 'C9ALS'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_C9ALS_BA9 <- do.call(rbind, df_list)

        ## Calculate sample wise mean
        final_C9ALS_BA9 <- final_C9ALS_BA9 %>%
            group_by(donor, celltype, group, n_cells) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy),        
                .groups = "drop"               
            )

        final_C9ALS_BA9$label <- paste0(round(final_C9ALS_BA9$sample_mean_accuracy, 2),'\n(',final_C9ALS_BA9$n_cells, ')' )

        final_C9ALS_BA9$celltype <- factor(final_C9ALS_BA9$celltype, levels = c('L4_L6','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Astro', 'OPC','T_Cell'  ))

        #final_C9ALS_BA9 <- subset(final_C9ALS_BA9, n_cells > 10)
        
        ## Plot
        BA9_C9ALS_heat <- ggplot(final_C9ALS_BA9, aes(x = celltype, y = donor, fill = sample_mean_accuracy, label = label)) + 
        geom_tile() +
        geom_text(size = 1.75) +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title = element_blank()
        ) +
        scale_fill_gradientn(colors = brewer.pal(9, "RdYlGn"))   # Use the "RdYlBu" palette for continuous data
        #coord_cartesian(ylim=c(0.5,1.0)) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 8, width = 7)
    ##

    ###########################
    ## SFTLD BA9
    ###########################
    ## code
        # List of cell types
        cell_type <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')
        status <- 'SFTLD'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_SFTLD_BA9 <- do.call(rbind, df_list)

        ## Calculate sample wise mean
        final_SFTLD_BA9 <- final_SFTLD_BA9 %>%
            group_by(donor, celltype, group, n_cells) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy),        
                .groups = "drop"               
            )

        final_SFTLD_BA9$label <- paste0(round(final_SFTLD_BA9$sample_mean_accuracy, 2),'\n(',final_SFTLD_BA9$n_cells, ')' )

        final_SFTLD_BA9$celltype <- factor(final_SFTLD_BA9$celltype, levels = c('L4_L6','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Astro', 'OPC','T_Cell'  ))

        #final_SFTLD_BA9 <- subset(final_SFTLD_BA9, n_cells > 10)
        
        ## Plot
        BA9_SFTLD_heat <- ggplot(final_SFTLD_BA9, aes(x = celltype, y = donor, fill = sample_mean_accuracy, label = label)) + 
        geom_tile() +
        geom_text(size = 1.75) +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title = element_blank()
        ) +
        scale_fill_gradientn(colors = brewer.pal(9, "RdYlGn"))   # Use the "RdYlBu" palette for continuous data
        #coord_cartesian(ylim=c(0.5,1.0)) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 8, width = 7)
        
    ##

    ###########################
    ## C9FTLD BA9
    ###########################
    ## code
        # List of cell types
        cell_type <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')
        status <- 'C9FTLD'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_C9FTLD_BA9 <- do.call(rbind, df_list)

        ## Calculate sample wise mean
        final_C9FTLD_BA9 <- final_C9FTLD_BA9 %>%
            group_by(donor, celltype, group, n_cells) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy),        
                .groups = "drop"               
            )

        final_C9FTLD_BA9$label <- paste0(round(final_C9FTLD_BA9$sample_mean_accuracy, 2),'\n(',final_C9FTLD_BA9$n_cells, ')' )

        final_C9FTLD_BA9$celltype <- factor(final_C9FTLD_BA9$celltype, levels = c('L4_L6','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Astro', 'OPC','T_Cell'  ))

        #final_C9FTLD_BA9 <- subset(final_C9FTLD_BA9, n_cells > 10)
        
        ## Plot
        BA9_C9FTLD_heat <- ggplot(final_C9FTLD_BA9, aes(x = celltype, y = donor, fill = sample_mean_accuracy, label = label)) + 
        geom_tile() +
        geom_text(size = 1.75) +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title = element_blank()
        ) +
        scale_fill_gradientn(colors = brewer.pal(9, "RdYlGn"))   # Use the "RdYlBu" palette for continuous data
        #coord_cartesian(ylim=c(0.5,1.0)) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 8, width = 7)
        
    ##
##



################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ All HVGs: LOO workflow ABS scaling: Accuracy
## code

    ## BA4
    
    ###########################
    ## SALS BA4
    ###########################
    ## code
        # List of cell types
        cell_type <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')
        status <- 'SALS'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_LOO_report_CombatSeq_',status,'_',brain_region,'_',i,'_abs_scaling_narval_2.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_SALS_BA4 <- do.call(rbind, df_list)
    ##

    ###########################
    ## C9ALS BA4
    ###########################
    ## code
        # List of cell types
        cell_type <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')
        status <- 'C9ALS'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_LOO_report_CombatSeq_',status,'_',brain_region,'_',i,'_abs_scaling_narval_2.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_C9ALS_BA4 <- do.call(rbind, df_list)
    ##

    ###########################
    ## SFTLD BA4
    ###########################
    ## code
        # List of cell types
        cell_type <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')
        status <- 'SFTLD'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_LOO_report_CombatSeq_',status,'_',brain_region,'_',i,'_abs_scaling_narval_2.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_SFTLD_BA4 <- do.call(rbind, df_list)
    ##

    ###########################
    ## C9FTLD BA4
    ###########################
    ## code
        # List of cell types
        cell_type <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')
        status <- 'C9FTLD'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_LOO_report_CombatSeq_',status,'_',brain_region,'_',i,'_abs_scaling_narval_2.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_C9FTLD_BA4 <- do.call(rbind, df_list)
    ##

    ###########################
    ## All FTLD BA4
    ###########################
    ## code
        # List of cell types
        cell_type <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')
        status <- 'All FTLD'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_LOO_report_CombatSeq_',status,'_',brain_region,'_',i,'_abs_scaling_narval_2.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_All_FTLD_BA4 <- do.call(rbind, df_list)
    ##

    ###########################
    ## All ALS BA4
    ###########################
    ## code
        # List of cell types
        cell_type <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')
        status <- 'All ALS'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_LOO_report_CombatSeq_',status,'_',brain_region,'_',i,'_abs_scaling_narval_2.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_All_ALS_BA4 <- do.call(rbind, df_list)
    ##

    ###########################
    ## merge and plot BA4
    ###########################
    ## code for barplot
        BA4_bind_all <- rbind(final_SALS_BA4, final_C9ALS_BA4, final_SFTLD_BA4, final_C9FTLD_BA4, final_All_ALS_BA4, final_All_FTLD_BA4)
        
        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, group) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy),        
                .groups = "drop"               
            )

        ## Calculate the mean and standard deviation
        summary_stats <- BA4_bind_all %>%
            group_by(celltype, group) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## set factor levels
        summary_stats$group <- factor(summary_stats$group, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All ALS", "All FTLD"))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L4_L6','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Astro', 'OPC','T_Cell'  ))


        ## Plot
        ggplot(summary_stats, aes(x = group, y = median_accuracy, fill = group, colour = group)) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~celltype, scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) #+
        #coord_cartesian(ylim=c(0.5,1.0)) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA4_LOO_balanced_acc_barplot_HVGs_abs_scaling.pdf', height = 2, width = 13)
    ##

    ## BA9
    
    ###########################
    ## SALS BA9
    ###########################
    ## code
        # List of cell types
        cell_type <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')
        status <- 'SALS'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_LOO_report_CombatSeq_',status,'_',brain_region,'_',i,'_abs_scaling_narval_2.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_SALS_BA9 <- do.call(rbind, df_list)
    ##

    ###########################
    ## C9ALS BA9
    ###########################
    ## code
        # List of cell types
        cell_type <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')
        status <- 'C9ALS'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_LOO_report_CombatSeq_',status,'_',brain_region,'_',i,'_abs_scaling_narval_2.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_C9ALS_BA9 <- do.call(rbind, df_list)
    ##

    ###########################
    ## SFTLD BA9
    ###########################
    ## code
        # List of cell types
        cell_type <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')
        status <- 'SFTLD'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_LOO_report_CombatSeq_',status,'_',brain_region,'_',i,'_abs_scaling_narval_2.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_SFTLD_BA9 <- do.call(rbind, df_list)
    ##

    ###########################
    ## C9FTLD BA9
    ###########################
    ## code
        # List of cell types
        cell_type <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')
        status <- 'C9FTLD'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_LOO_report_CombatSeq_',status,'_',brain_region,'_',i,'_abs_scaling_narval_2.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_C9FTLD_BA9 <- do.call(rbind, df_list)
    ##

    ###########################
    ## All FTLD BA9
    ###########################
    ## code
        # List of cell types
        cell_type <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')
        status <- 'All FTLD'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_LOO_report_CombatSeq_',status,'_',brain_region,'_',i,'_abs_scaling_narval_2.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_All_FTLD_BA9 <- do.call(rbind, df_list)
    ##

    ###########################
    ## All ALS BA9
    ###########################
    ## code
        # List of cell types
        cell_type <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')
        status <- 'All ALS'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_LOO_report_CombatSeq_',status,'_',brain_region,'_',i,'_abs_scaling_narval_2.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_All_ALS_BA9 <- do.call(rbind, df_list)
    ##

    ###########################
    ## merge and plot BA9
    ###########################
    ## code for barplot
        BA9_bind_all <- rbind(final_SALS_BA9, final_C9ALS_BA9, final_SFTLD_BA9, final_C9FTLD_BA9, final_All_ALS_BA9, final_All_FTLD_BA9)
        
        ## Calculate sample wise mean
        BA9_bind_all <- BA9_bind_all %>%
            group_by(donor, celltype, group) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy),        
                .groups = "drop"               
            )

        ## Calculate the mean and standard deviation
        summary_stats <- BA9_bind_all %>%
            group_by(celltype, group) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## set factor levels
        summary_stats$group <- factor(summary_stats$group, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All ALS", "All FTLD"))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L4_L6','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Astro', 'OPC','T_Cell'  ))


        ## Plot
        ggplot(summary_stats, aes(x = group, y = median_accuracy, fill = group, colour = group)) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~celltype, scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) #+
        #coord_cartesian(ylim=c(0.5,1.0)) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA9_LOO_balanced_acc_barplot_HVGs_abs_scaling.pdf', height = 2, width = 13)
    ##

##

############################################################################################################################################################################### HVG methods

################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ HVGs: Number of highly variable genes
## code
    
    #BA4
    
    ###########################
    ## SALS BA4
    ###########################
    ## code
        ## define variables
        condition = "SALS" #one ofSALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            HVGs <- as.numeric(temp$n_HVGs)
            result_df <- data.frame(HVGs = HVGs)
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
        condition = "C9ALS" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        #cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            HVGs <- as.numeric(temp$n_HVGs)
            result_df <- data.frame(HVGs = HVGs)
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
        condition = "SFTLD" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        #cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            HVGs <- as.numeric(temp$n_HVGs)
            result_df <- data.frame(HVGs = HVGs)
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
        condition = "C9FTLD" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        #cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            HVGs <- as.numeric(temp$n_HVGs)
            result_df <- data.frame(HVGs = HVGs)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df

        }

        final_df_C9FTLD_BA4 <- do.call(rbind, result_list)
    ##

    ###########################
    ## all_conditions BA4
    ###########################
    ## code
        ## define variables
        condition = "all_conditions" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        #cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')


        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            HVGs <- as.numeric(temp$n_HVGs)
            result_df <- data.frame(HVGs = HVGs)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df

        }

        final_df_all_conditions_BA4 <- do.call(rbind, result_list)
    ##

    ###########################
    ## All ALS BA4
    ###########################
    ## code
        ## define variables
        condition = "All ALS" #one ofSALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_LOO_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            HVGs <- as.numeric(temp$n_genes)
            result_df <- data.frame(HVGs = HVGs)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df
        }

        final_df_All_ALS_BA4 <- do.call(rbind, result_list)
    ##

    ###########################
    ## All FTLD BA4
    ###########################
    ## code
        ## define variables
        condition = "All FTLD" #one ofSALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')

        result_list <- list()
        
        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_LOO_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            HVGs <- as.numeric(temp$n_genes)
            result_df <- data.frame(HVGs = HVGs)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df
        }

        final_df_All_FTLD_BA4 <- do.call(rbind, result_list)
    ##

    ###########################
    ## merge and plot BA4
    ###########################
    ## code
        BA4_bind_all <- rbind(final_df_SALS_BA4, final_df_C9ALS_BA4, final_df_SFTLD_BA4, final_df_C9FTLD_BA4, final_df_All_ALS_BA4, final_df_All_FTLD_BA4)
        BA4_bind_all <- unique(BA4_bind_all)
        
        ## set factor levels
        BA4_bind_all$condition <- factor(BA4_bind_all$condition, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All ALS", "All FTLD"))
        BA4_bind_all$celltype <- factor(BA4_bind_all$celltype, levels = c('L4_L6','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Astro', 'OPC'))


        ggplot(BA4_bind_all, aes(x = condition, y = HVGs, fill = condition, colour = condition)) + 
        theme_bw() + 
        geom_bar(stat = 'identity') +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~celltype, scales = "free_x") +
        ylab("n HVGs") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue"))
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA4_nHVGs.pdf', height = 2, width = 13)
    ##

    #BA9
    
    ###########################
    ## SALS BA9
    ###########################
    ## code
        ## define variables
        condition = "SALS" #one ofSALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            HVGs <- as.numeric(temp$n_HVGs)
            result_df <- data.frame(HVGs = HVGs)
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
        condition = "C9ALS" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        #cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            HVGs <- as.numeric(temp$n_HVGs)
            result_df <- data.frame(HVGs = HVGs)
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
        condition = "SFTLD" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        #cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            HVGs <- as.numeric(temp$n_HVGs)
            result_df <- data.frame(HVGs = HVGs)
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
        condition = "C9FTLD" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        #cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            HVGs <- as.numeric(temp$n_HVGs)
            result_df <- data.frame(HVGs = HVGs)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df

        }

        final_df_C9FTLD_BA9 <- do.call(rbind, result_list)
    ##

    ###########################
    ## all_conditions BA9
    ###########################
    ## code
        ## define variables
        condition = "all_conditions" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        #cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')


        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            HVGs <- as.numeric(temp$n_HVGs)
            result_df <- data.frame(HVGs = HVGs)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df

        }

        final_df_all_conditions_BA9 <- do.call(rbind, result_list)
    ##

    ###########################
    ## All ALS BA9
    ###########################
    ## code
        ## define variables
        condition = "All ALS" #one ofSALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_LOO_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            HVGs <- as.numeric(temp$n_genes)
            result_df <- data.frame(HVGs = HVGs)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df
        }

        final_df_All_ALS_BA9 <- do.call(rbind, result_list)
    ##

    ###########################
    ## All FTLD BA9
    ###########################
    ## code
        ## define variables
        condition = "All FTLD" #one ofSALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM')

        result_list <- list()
        
        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_LOO_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            HVGs <- as.numeric(temp$n_genes)
            result_df <- data.frame(HVGs = HVGs)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df
        }

        final_df_All_FTLD_BA9 <- do.call(rbind, result_list)
    ##

    ###########################
    ## merge and plot BA9
    ###########################
    ## code
        BA9_bind_all <- rbind(final_df_SALS_BA9, final_df_C9ALS_BA9, final_df_SFTLD_BA9, final_df_C9FTLD_BA9, final_df_All_ALS_BA9, final_df_All_FTLD_BA9)
        BA9_bind_all <- unique(BA9_bind_all)
        
        ## set factor levels
        BA9_bind_all$condition <- factor(BA9_bind_all$condition, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All ALS", "All FTLD"))
        BA9_bind_all$celltype <- factor(BA9_bind_all$celltype, levels = c('L4_L6','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Astro', 'OPC'))


        ggplot(BA9_bind_all, aes(x = condition, y = HVGs, fill = condition, colour = condition)) + 
        theme_bw() + 
        geom_bar(stat = 'identity') +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~celltype, scales = "free_x") +
        ylab("n HVGs") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue"))
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA9_nHVGs.pdf', height = 2, width = 13)
    ##

##


################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ HVGs: LOO workflow: Accuracy
## code
    
    ## BA4
    
    ###########################
    ## SALS BA4
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        status <- 'SALS'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2_LOO.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_SALS_BA4 <- do.call(rbind, df_list)
    ##

    ###########################
    ## C9ALS BA4
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'C9ALS'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2_LOO.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_C9ALS_BA4 <- do.call(rbind, df_list)
    ##

    ###########################
    ## SFTLD BA4
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        
        status <- 'SFTLD'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2_LOO.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_SFTLD_BA4 <- do.call(rbind, df_list)
    ##

    ###########################
    ## C9FTLD BA4
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        
        status <- 'C9FTLD'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2_LOO.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_C9FTLD_BA4 <- do.call(rbind, df_list)
    ##

    ###########################
    ## All ALS BA4
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        
        status <- 'All ALS'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_LOO_report_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_All_ALS_BA4 <- do.call(rbind, df_list)
    ##

    ###########################
    ## All FTLD BA4
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        
        status <- 'All FTLD'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_LOO_report_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_All_FTLD_BA4 <- do.call(rbind, df_list)
    ##

    ###########################
    ## all_conditions BA4
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        status <- 'all_conditions'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2_LOO.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_all_conditions_BA4 <- do.call(rbind, df_list)
    ##

    ###########################
    ## merge and plot BA4
    ###########################
    ## code for barplot
        final_All_ALS_BA4 <- final_All_ALS_BA4 %>% dplyr::select(prep, donor, region, group, celltype, n_genes, n_cells, learning_rate, batch_size, test_accuracy, dataset)
        colnames(final_All_ALS_BA4) <- c('prep', 'donor', 'region', 'group', 'celltype', 'n_HVGs', 'n_cells', 'learning_rate', 'batch_size', 'test_accuracy', 'dataset')
        final_All_FTLD_BA4 <- final_All_FTLD_BA4 %>% dplyr::select(prep, donor, region, group, celltype, n_genes, n_cells, learning_rate, batch_size, test_accuracy, dataset)
        colnames(final_All_FTLD_BA4) <- c('prep', 'donor', 'region', 'group', 'celltype', 'n_HVGs', 'n_cells', 'learning_rate', 'batch_size', 'test_accuracy', 'dataset')
        
        BA4_bind_all <- rbind(final_SALS_BA4, final_C9ALS_BA4, final_SFTLD_BA4, final_C9FTLD_BA4, final_All_ALS_BA4, final_All_FTLD_BA4)
        write.csv(BA4_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_LOO_HVGs.csv')

        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, group) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy),        
                .groups = "drop"               
            )
        
        ## Calculate the mean and standard deviation
        summary_stats <- BA4_bind_all %>%
            group_by(celltype, group) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## set factor levels
        summary_stats$group <- factor(summary_stats$group, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All ALS", "All FTLD"))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC' ))


        ## Plot
        ggplot(summary_stats, aes(x = group, y = median_accuracy, fill = group, colour = group)) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~celltype, scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) #+
        #coord_cartesian(ylim=c(0.5,1.0)) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA4_LOO_balanced_acc_barplot.pdf', height = 2, width = 13)
    ##

    ## BA9
    
    ###########################
    ## SALS BA9
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")


        status <- 'SALS'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2_LOO.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_SALS_BA9 <- do.call(rbind, df_list)
    ##

    ###########################
    ## C9ALS BA9
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        status <- 'C9ALS'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2_LOO.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_C9ALS_BA9 <- do.call(rbind, df_list)
    ##

    ###########################
    ## SFTLD BA9
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        
        status <- 'SFTLD'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2_LOO.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_SFTLD_BA9 <- do.call(rbind, df_list)
    ##

    ###########################
    ## C9FTLD BA9
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        
        status <- 'C9FTLD'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2_LOO.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_C9FTLD_BA9 <- do.call(rbind, df_list)
    ##

    ###########################
    ## All ALS BA9
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        
        status <- 'All ALS'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_LOO_report_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_All_ALS_BA9 <- do.call(rbind, df_list)
    ##

    ###########################
    ## All FTLD BA9
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        
        status <- 'All FTLD'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_LOO_report_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_All_FTLD_BA9 <- do.call(rbind, df_list)
    ##

    ###########################
    ## all_conditions BA9
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")


        status <- 'all_conditions'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2_LOO.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_all_conditions_BA9 <- do.call(rbind, df_list)
    ##

    ###########################
    ## merge and plot BA9
    ###########################
    ## code for barplot
        final_All_ALS_BA9 <- final_All_ALS_BA9 %>% dplyr::select(prep, donor, region, group, celltype, n_genes, n_cells, learning_rate, batch_size, test_accuracy, dataset)
        colnames(final_All_ALS_BA9) <- c('prep', 'donor', 'region', 'group', 'celltype', 'n_HVGs', 'n_cells', 'learning_rate', 'batch_size', 'test_accuracy', 'dataset')
        final_All_FTLD_BA9 <- final_All_FTLD_BA9 %>% dplyr::select(prep, donor, region, group, celltype, n_genes, n_cells, learning_rate, batch_size, test_accuracy, dataset)
        colnames(final_All_FTLD_BA9) <- c('prep', 'donor', 'region', 'group', 'celltype', 'n_HVGs', 'n_cells', 'learning_rate', 'batch_size', 'test_accuracy', 'dataset')
        
        BA9_bind_all <- rbind(final_SALS_BA9, final_C9ALS_BA9, final_SFTLD_BA9, final_C9FTLD_BA9, final_All_ALS_BA9, final_All_FTLD_BA9)
        write.csv(BA9_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_LOO_HVGs.csv')

        BA9_bind_all <- BA9_bind_all %>%
            group_by(donor, celltype, group) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy),        
                .groups = "drop"               
            )
        
        ## Calculate the mean and standard deviation
        summary_stats <- BA9_bind_all %>%
            group_by(celltype, group) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## set factor levels
        summary_stats$group <- factor(summary_stats$group, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All ALS", "All FTLD"))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC' ))


        ## Plot
        ggplot(summary_stats, aes(x = group, y = median_accuracy, fill = group, colour = group)) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~celltype, scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) #+
        #coord_cartesian(ylim=c(0.5,1.0)) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA9_LOO_balanced_acc_barplot.pdf', height = 2, width = 13)
    ##

 
##

################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ HVGs: LOO workflow ABS scaling: Accuracy
## code

    ## BA4
    
    ###########################
    ## SALS BA4
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'SALS'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_LOO_report_CombatSeq_',status,'_',brain_region,'_',i,'_abs_scaling_narval_2.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_SALS_BA4 <- do.call(rbind, df_list)
    ##

    ###########################
    ## C9ALS BA4
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'C9ALS'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_LOO_report_CombatSeq_',status,'_',brain_region,'_',i,'_abs_scaling_narval_2.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_C9ALS_BA4 <- do.call(rbind, df_list)
    ##

    ###########################
    ## SFTLD BA4
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'SFTLD'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_LOO_report_CombatSeq_',status,'_',brain_region,'_',i,'_abs_scaling_narval_2.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_SFTLD_BA4 <- do.call(rbind, df_list)
    ##

    ###########################
    ## C9FTLD BA4
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'C9FTLD'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_LOO_report_CombatSeq_',status,'_',brain_region,'_',i,'_abs_scaling_narval_2.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_C9FTLD_BA4 <- do.call(rbind, df_list)
    ##

    ###########################
    ## All FTLD BA4
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'All FTLD'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_LOO_report_CombatSeq_',status,'_',brain_region,'_',i,'_abs_scaling_narval_2.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_All_FTLD_BA4 <- do.call(rbind, df_list)
    ##

    ###########################
    ## All ALS BA4
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'All ALS'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_LOO_report_CombatSeq_',status,'_',brain_region,'_',i,'_abs_scaling_narval_2.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_All_ALS_BA4 <- do.call(rbind, df_list)
    ##

    ###########################
    ## merge and plot BA4
    ###########################
    ## code for barplot
        BA4_bind_all <- rbind(final_SALS_BA4, final_C9ALS_BA4, final_SFTLD_BA4, final_C9FTLD_BA4, final_All_ALS_BA4, final_All_FTLD_BA4)
        write.csv(BA4_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_LOO_HVGs_ABS_scaling.csv')

        
        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, group) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy),        
                .groups = "drop"               
            )

        ## Calculate the mean and standard deviation
        summary_stats <- BA4_bind_all %>%
            group_by(celltype, group) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## set factor levels
        summary_stats$group <- factor(summary_stats$group, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All ALS", "All FTLD"))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC' ))


        ## Plot
        ggplot(summary_stats, aes(x = group, y = median_accuracy, fill = group, colour = group)) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~celltype, scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) #+
        #coord_cartesian(ylim=c(0.5,1.0)) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA4_LOO_balanced_acc_barplot_HVGs_abs_scaling.pdf', height = 2, width = 13)
    ##

    ## BA9
    
    ###########################
    ## SALS BA9
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'SALS'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_LOO_report_CombatSeq_',status,'_',brain_region,'_',i,'_abs_scaling_narval_2.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_SALS_BA9 <- do.call(rbind, df_list)
    ##

    ###########################
    ## C9ALS BA9
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'C9ALS'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_LOO_report_CombatSeq_',status,'_',brain_region,'_',i,'_abs_scaling_narval_2.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_C9ALS_BA9 <- do.call(rbind, df_list)
    ##

    ###########################
    ## SFTLD BA9
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'SFTLD'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_LOO_report_CombatSeq_',status,'_',brain_region,'_',i,'_abs_scaling_narval_2.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_SFTLD_BA9 <- do.call(rbind, df_list)
    ##

    ###########################
    ## C9FTLD BA9
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'C9FTLD'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_LOO_report_CombatSeq_',status,'_',brain_region,'_',i,'_abs_scaling_narval_2.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_C9FTLD_BA9 <- do.call(rbind, df_list)
    ##

    ###########################
    ## All FTLD BA9
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'All FTLD'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_LOO_report_CombatSeq_',status,'_',brain_region,'_',i,'_abs_scaling_narval_2.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_All_FTLD_BA9 <- do.call(rbind, df_list)
    ##

    ###########################
    ## All ALS BA9
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'All ALS'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_LOO_report_CombatSeq_',status,'_',brain_region,'_',i,'_abs_scaling_narval_2.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_All_ALS_BA9 <- do.call(rbind, df_list)
    ##

    ###########################
    ## merge and plot BA9
    ###########################
    ## code for barplot
        BA9_bind_all <- rbind(final_SALS_BA9, final_C9ALS_BA9, final_SFTLD_BA9, final_C9FTLD_BA9, final_All_ALS_BA9, final_All_FTLD_BA9)
        write.csv(BA9_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_LOO_HVGs_ABS_scaling.csv')
        
        ## Calculate sample wise mean
        BA9_bind_all <- BA9_bind_all %>%
            group_by(donor, celltype, group) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy),        
                .groups = "drop"               
            )

        ## Calculate the mean and standard deviation
        summary_stats <- BA9_bind_all %>%
            group_by(celltype, group) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## set factor levels
        summary_stats$group <- factor(summary_stats$group, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All ALS", "All FTLD"))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC' ))


        ## Plot
        ggplot(summary_stats, aes(x = group, y = median_accuracy, fill = group, colour = group)) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~celltype, scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) #+
        #coord_cartesian(ylim=c(0.5,1.0)) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA9_LOO_balanced_acc_barplot_HVGs_abs_scaling.pdf', height = 2, width = 13)
    ##

##

################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ HVGs: LOO workflow ABS scaling KNN Bayesian: Accuracy
## code
    ## BA4
    
    ###########################
    ## SALS BA4
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'SALS'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVGs_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_',status,'_',brain_region,'_',i,'_abs_scaling_sweep_confidence_threshold_narval_2.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_SALS_BA4 <- do.call(rbind, df_list)
    ##

    ###########################
    ## C9ALS BA4
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'C9ALS'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVGs_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_',status,'_',brain_region,'_',i,'_abs_scaling_sweep_confidence_threshold_narval_2.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_C9ALS_BA4 <- do.call(rbind, df_list)
    ##

    ###########################
    ## SFTLD BA4
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'SFTLD'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVGs_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_',status,'_',brain_region,'_',i,'_abs_scaling_sweep_confidence_threshold_narval_2.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_SFTLD_BA4 <- do.call(rbind, df_list)
    ##

    ###########################
    ## C9FTLD BA4
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'C9FTLD'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVGs_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_',status,'_',brain_region,'_',i,'_abs_scaling_sweep_confidence_threshold_narval_2.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_C9FTLD_BA4 <- do.call(rbind, df_list)
    ##

    ###########################
    ## All ALS BA4
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'All ALS'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVGs_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_',status,'_',brain_region,'_',i,'_abs_scaling_sweep_confidence_threshold_narval_2.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_All_ALS_BA4 <- do.call(rbind, df_list)
    ##

    ###########################
    ## All FTLD BA4
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'All FTLD'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVGs_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_',status,'_',brain_region,'_',i,'_abs_scaling_sweep_confidence_threshold_narval_2.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_All_FTLD_BA4 <- do.call(rbind, df_list)
    ##

    ###########################
    ## merge and plot BA4
    ###########################
    ## code for barplot
        BA4_bind_all <- rbind(final_SALS_BA4, final_C9ALS_BA4, final_SFTLD_BA4, final_C9FTLD_BA4, final_All_ALS_BA4, final_All_FTLD_BA4)
        write.csv(BA4_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_LOO_HVGs_ABS_scaling_KNN_Bayesian.csv')

        
        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, group, kNN_thresh) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )

        ## Calculate the mean and standard deviation
        summary_stats <- BA4_bind_all %>%
            group_by(celltype, group, kNN_thresh) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## set factor levels
        summary_stats$group <- factor(summary_stats$group, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All ALS", "All FTLD"))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC' ))

        ## Plot
        ggplot(summary_stats, aes(x = group, y = median_accuracy, fill = group, colour = group)) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(kNN_thresh~celltype, scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) #+
        #coord_cartesian(ylim=c(0.5,1.0)) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA4_LOO_balanced_acc_barplot_HVGs_abs_scaling_KNN_bayesian.pdf', height = 8, width = 13)

        ## Threshold wise median
        summary_stats_thresh <- BA4_bind_all %>%
            group_by(kNN_thresh) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )
        ggplot(summary_stats_thresh, aes(x = kNN_thresh, y = median_accuracy, label = round(median_accuracy, 4))) + 
        geom_text(vjust = -0.5, size = 2) +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        scale_fill_manual(values = "black") +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        ylab("Accuracy") 
        #coord_cartesian(ylim=c(0.5,1.0)) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA4_LOO_balanced_acc_barplot_HVGs_threshold_comparison.pdf', height = 4, width = 4)

    ##

    ## BA9
    
    ###########################
    ## SALS BA9
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'SALS'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVGs_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_',status,'_',brain_region,'_',i,'_abs_scaling_sweep_confidence_threshold_narval_2.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_SALS_BA9 <- do.call(rbind, df_list)
    ##

    ###########################
    ## C9ALS BA9
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'C9ALS'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVGs_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_',status,'_',brain_region,'_',i,'_abs_scaling_sweep_confidence_threshold_narval_2.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_C9ALS_BA9 <- do.call(rbind, df_list)
    ##

    ###########################
    ## SFTLD BA9
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'SFTLD'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVGs_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_',status,'_',brain_region,'_',i,'_abs_scaling_sweep_confidence_threshold_narval_2.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_SFTLD_BA9 <- do.call(rbind, df_list)
    ##

    ###########################
    ## C9FTLD BA9
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'C9FTLD'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVGs_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_',status,'_',brain_region,'_',i,'_abs_scaling_sweep_confidence_threshold_narval_2.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_C9FTLD_BA9 <- do.call(rbind, df_list)
    ##

    ###########################
    ## All ALS BA9
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'All ALS'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVGs_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_',status,'_',brain_region,'_',i,'_abs_scaling_sweep_confidence_threshold_narval_2.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_All_ALS_BA9 <- do.call(rbind, df_list)
    ##

    ###########################
    ## All FTLD BA9
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'All FTLD'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVGs_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_',status,'_',brain_region,'_',i,'_abs_scaling_sweep_confidence_threshold_narval_2.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_All_FTLD_BA9 <- do.call(rbind, df_list)
    ##

    ###########################
    ## merge and plot BA9
    ###########################
    ## code for barplot
        BA9_bind_all <- rbind(final_SALS_BA9, final_C9ALS_BA9, final_SFTLD_BA9, final_C9FTLD_BA9, final_All_ALS_BA9, final_All_FTLD_BA9)
        write.csv(BA9_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_LOO_HVGs_ABS_scaling_KNN_Bayesian.csv')
        
        ## Calculate sample wise mean
        BA9_bind_all <- BA9_bind_all %>%
            group_by(donor, celltype, group, kNN_thresh) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )

        ## Calculate the mean and standard deviation
        summary_stats <- BA9_bind_all %>%
            group_by(celltype, group, kNN_thresh) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## set factor levels
        summary_stats$group <- factor(summary_stats$group, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All ALS", "All FTLD"))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC' ))

        ## Plot
        ggplot(summary_stats, aes(x = group, y = median_accuracy, fill = group, colour = group)) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(kNN_thresh~celltype, scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) #+
        #coord_cartesian(ylim=c(0.5,1.0)) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA9_LOO_balanced_acc_barplot_HVGs_abs_scaling_KNN_bayesian.pdf', height = 8, width = 13)

        ## Threshold wise median
        summary_stats_thresh <- BA9_bind_all %>%
            group_by(kNN_thresh) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )
        ggplot(summary_stats_thresh, aes(x = kNN_thresh, y = median_accuracy, label = round(median_accuracy, 4))) + 
        geom_text(vjust = -0.5, size = 2) +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        scale_fill_manual(values = "black") +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        ylab("Accuracy") 
        #coord_cartesian(ylim=c(0.5,1.0)) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA9_LOO_balanced_acc_barplot_HVGs_threshold_comparison.pdf', height = 4, width = 4)

    ##

##


############################################################################################################################################################################### Optimal methods


################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ Optimal: LOO workflow 0 delta threshold: Number of sig genes
## code
    
    #BA4
    
    ###########################
    ## SALS BA4
    ###########################
    ## code
        ## define variables
        par_brain_region = "BA4"
        par_status = "SALS"
        par_prep = "CombatSeq"
            
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
            "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
            "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        df_list <- list()
        
        for (celltype2 in celltype_list) {
        file_path <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_",
                            par_prep, "_", par_status, "_", par_brain_region, "_", 
                            celltype2, "_0_abs_case_control_narval_2.csv")

        temp_df <- read_csv(file_path, col_types = cols(.default = col_character())) %>%
            mutate(celltype = celltype2)
        n_genes = length((unique(temp_df$gene)))
        HVGs <- as.numeric(n_genes)
        result_df <- data.frame(HVGs = HVGs)
        result_df$celltype <- celltype2
        result_df$condition <- par_status
        result_df$region <- par_brain_region
        df_list[[celltype2]] <- result_df
        }

        final_df_SALS_BA4 <- do.call(rbind, df_list) 
    ##

    ###########################
    ## C9ALS BA4
    ###########################
    ## code
        ## define variables
        par_brain_region = "BA4"
        par_status = "C9ALS"
        par_prep = "CombatSeq"
            
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
            "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
            "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        df_list <- list()
        
        for (celltype2 in celltype_list) {
        file_path <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_",
                            par_prep, "_", par_status, "_", par_brain_region, "_", 
                            celltype2, "_0_abs_case_control_narval_2.csv")

        temp_df <- read_csv(file_path, col_types = cols(.default = col_character())) %>%
            mutate(celltype = celltype2)
        n_genes = length((unique(temp_df$gene)))
        HVGs <- as.numeric(n_genes)
        result_df <- data.frame(HVGs = HVGs)
        result_df$celltype <- celltype2
        result_df$condition <- par_status
        result_df$region <- par_brain_region
        df_list[[celltype2]] <- result_df
        }

        final_df_C9ALS_BA4 <- do.call(rbind, df_list) 
    ##

    ###########################
    ## SFTLD BA4
    ###########################
    ## code
        ## define variables
        par_brain_region = "BA4"
        par_status = "SFTLD"
        par_prep = "CombatSeq"
            
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
            "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
            "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        df_list <- list()
        
        for (celltype2 in celltype_list) {
        file_path <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_",
                            par_prep, "_", par_status, "_", par_brain_region, "_", 
                            celltype2, "_0_abs_case_control_narval_2.csv")

        temp_df <- read_csv(file_path, col_types = cols(.default = col_character())) %>%
            mutate(celltype = celltype2)
        n_genes = length((unique(temp_df$gene)))
        HVGs <- as.numeric(n_genes)
        result_df <- data.frame(HVGs = HVGs)
        result_df$celltype <- celltype2
        result_df$condition <- par_status
        result_df$region <- par_brain_region
        df_list[[celltype2]] <- result_df
        }

        final_df_SFTLD_BA4 <- do.call(rbind, df_list) 
    ##

    ###########################
    ## C9FTLD BA4
    ###########################
    ## code
        ## define variables
        par_brain_region = "BA4"
        par_status = "C9FTLD"
        par_prep = "CombatSeq"
            
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
            "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
            "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        df_list <- list()
        
        for (celltype2 in celltype_list) {
        file_path <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_",
                            par_prep, "_", par_status, "_", par_brain_region, "_", 
                            celltype2, "_0_abs_case_control_narval_2.csv")

        temp_df <- read_csv(file_path, col_types = cols(.default = col_character())) %>%
            mutate(celltype = celltype2)
        n_genes = length((unique(temp_df$gene)))
        HVGs <- as.numeric(n_genes)
        result_df <- data.frame(HVGs = HVGs)
        result_df$celltype <- celltype2
        result_df$condition <- par_status
        result_df$region <- par_brain_region
        df_list[[celltype2]] <- result_df
        }

        final_df_C9FTLD_BA4 <- do.call(rbind, df_list) 
    ##

    ###########################
    ## All ALS BA4
    ###########################
    ## code
        ## define variables
        par_brain_region = "BA4"
        par_status = "All ALS"
        par_prep = "CombatSeq"
            
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
            "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
            "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        df_list <- list()
        
        for (celltype2 in celltype_list) {
        file_path <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_",
                            par_prep, "_", par_status, "_", par_brain_region, "_", 
                            celltype2, "_0_abs_case_control_narval_2.csv")

        temp_df <- read_csv(file_path, col_types = cols(.default = col_character())) %>%
            mutate(celltype = celltype2)
        n_genes = length((unique(temp_df$gene)))
        HVGs <- as.numeric(n_genes)
        result_df <- data.frame(HVGs = HVGs)
        result_df$celltype <- celltype2
        result_df$condition <- par_status
        result_df$region <- par_brain_region
        df_list[[celltype2]] <- result_df
        }

        final_df_All_ALS_BA4 <- do.call(rbind, df_list) 
    ##

    ###########################
    ## All FTLD BA4
    ###########################
    ## code
        ## define variables
        par_brain_region = "BA4"
        par_status = "All FTLD"
        par_prep = "CombatSeq"
            
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
            "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
            "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        df_list <- list()
        
        for (celltype2 in celltype_list) {
        file_path <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_",
                            par_prep, "_", par_status, "_", par_brain_region, "_", 
                            celltype2, "_0_abs_case_control_narval_2.csv")

        temp_df <- read_csv(file_path, col_types = cols(.default = col_character())) %>%
            mutate(celltype = celltype2)
        n_genes = length((unique(temp_df$gene)))
        HVGs <- as.numeric(n_genes)
        result_df <- data.frame(HVGs = HVGs)
        result_df$celltype <- celltype2
        result_df$condition <- par_status
        result_df$region <- par_brain_region
        df_list[[celltype2]] <- result_df
        }

        final_df_All_FTLD_BA4 <- do.call(rbind, df_list) 
    ##

    
    ###########################
    ## merge and plot BA4
    ###########################
    ## code
        BA4_bind_all <- rbind(final_df_SALS_BA4, final_df_C9ALS_BA4, final_df_SFTLD_BA4, final_df_C9FTLD_BA4, final_df_All_ALS_BA4, final_df_All_FTLD_BA4)

        ## set factor levels
        BA4_bind_all$condition <- factor(BA4_bind_all$condition, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All ALS", "All FTLD"))
        BA4_bind_all$celltype <- factor(BA4_bind_all$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))


        ggplot(BA4_bind_all, aes(x = condition, y = HVGs, fill = condition, colour = condition)) + 
        theme_bw() + 
        geom_bar(stat = 'identity') +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~celltype, scales = "free_x") +
        ylab("n HVGs") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue"))
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA4_n_sig_genes_0.pdf', height = 2, width = 13)
    ##

    #BA9
    
    ###########################
    ## SALS BA9
    ###########################
    ## code
        ## define variables
        par_brain_region = "BA9"
        par_status = "SALS"
        par_prep = "CombatSeq"
            
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
            "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
            "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        df_list <- list()
        
        for (celltype2 in celltype_list) {
        file_path <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_",
                            par_prep, "_", par_status, "_", par_brain_region, "_", 
                            celltype2, "_0_abs_case_control_narval_2.csv")

        temp_df <- read_csv(file_path, col_types = cols(.default = col_character())) %>%
            mutate(celltype = celltype2)
        n_genes = length((unique(temp_df$gene)))
        HVGs <- as.numeric(n_genes)
        result_df <- data.frame(HVGs = HVGs)
        result_df$celltype <- celltype2
        result_df$condition <- par_status
        result_df$region <- par_brain_region
        df_list[[celltype2]] <- result_df
        }

        final_df_SALS_BA9 <- do.call(rbind, df_list) 
    ##

    ###########################
    ## C9ALS BA9
    ###########################
    ## code
        ## define variables
        par_brain_region = "BA9"
        par_status = "C9ALS"
        par_prep = "CombatSeq"
            
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
            "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
            "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        df_list <- list()
        
        for (celltype2 in celltype_list) {
        file_path <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_",
                            par_prep, "_", par_status, "_", par_brain_region, "_", 
                            celltype2, "_0_abs_case_control_narval_2.csv")

        temp_df <- read_csv(file_path, col_types = cols(.default = col_character())) %>%
            mutate(celltype = celltype2)
        n_genes = length((unique(temp_df$gene)))
        HVGs <- as.numeric(n_genes)
        result_df <- data.frame(HVGs = HVGs)
        result_df$celltype <- celltype2
        result_df$condition <- par_status
        result_df$region <- par_brain_region
        df_list[[celltype2]] <- result_df
        }

        final_df_C9ALS_BA9 <- do.call(rbind, df_list) 
    ##

    ###########################
    ## SFTLD BA9
    ###########################
    ## code
        ## define variables
        par_brain_region = "BA9"
        par_status = "SFTLD"
        par_prep = "CombatSeq"
            
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
            "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
            "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        df_list <- list()
        
        for (celltype2 in celltype_list) {
        file_path <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_",
                            par_prep, "_", par_status, "_", par_brain_region, "_", 
                            celltype2, "_0_abs_case_control_narval_2.csv")

        temp_df <- read_csv(file_path, col_types = cols(.default = col_character())) %>%
            mutate(celltype = celltype2)
        n_genes = length((unique(temp_df$gene)))
        HVGs <- as.numeric(n_genes)
        result_df <- data.frame(HVGs = HVGs)
        result_df$celltype <- celltype2
        result_df$condition <- par_status
        result_df$region <- par_brain_region
        df_list[[celltype2]] <- result_df
        }

        final_df_SFTLD_BA9 <- do.call(rbind, df_list) 
    ##

    ###########################
    ## C9FTLD BA9
    ###########################
    ## code
        ## define variables
        par_brain_region = "BA9"
        par_status = "C9FTLD"
        par_prep = "CombatSeq"
            
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
            "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
            "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        df_list <- list()
        
        for (celltype2 in celltype_list) {
        file_path <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_",
                            par_prep, "_", par_status, "_", par_brain_region, "_", 
                            celltype2, "_0_abs_case_control_narval_2.csv")

        temp_df <- read_csv(file_path, col_types = cols(.default = col_character())) %>%
            mutate(celltype = celltype2)
        n_genes = length((unique(temp_df$gene)))
        HVGs <- as.numeric(n_genes)
        result_df <- data.frame(HVGs = HVGs)
        result_df$celltype <- celltype2
        result_df$condition <- par_status
        result_df$region <- par_brain_region
        df_list[[celltype2]] <- result_df
        }

        final_df_C9FTLD_BA9 <- do.call(rbind, df_list) 
    ##

    ###########################
    ## All ALS BA9
    ###########################
    ## code
        ## define variables
        par_brain_region = "BA9"
        par_status = "All ALS"
        par_prep = "CombatSeq"
            
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
            "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
            "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        df_list <- list()
        
        for (celltype2 in celltype_list) {
        file_path <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_",
                            par_prep, "_", par_status, "_", par_brain_region, "_", 
                            celltype2, "_0_abs_case_control_narval_2.csv")

        temp_df <- read_csv(file_path, col_types = cols(.default = col_character())) %>%
            mutate(celltype = celltype2)
        n_genes = length((unique(temp_df$gene)))
        HVGs <- as.numeric(n_genes)
        result_df <- data.frame(HVGs = HVGs)
        result_df$celltype <- celltype2
        result_df$condition <- par_status
        result_df$region <- par_brain_region
        df_list[[celltype2]] <- result_df
        }

        final_df_All_ALS_BA9 <- do.call(rbind, df_list) 
    ##

    ###########################
    ## All FTLD BA9
    ###########################
    ## code
        ## define variables
        par_brain_region = "BA9"
        par_status = "All FTLD"
        par_prep = "CombatSeq"
            
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
            "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
            "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        df_list <- list()
        
        for (celltype2 in celltype_list) {
        file_path <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_",
                            par_prep, "_", par_status, "_", par_brain_region, "_", 
                            celltype2, "_0_abs_case_control_narval_2.csv")

        temp_df <- read_csv(file_path, col_types = cols(.default = col_character())) %>%
            mutate(celltype = celltype2)
        n_genes = length((unique(temp_df$gene)))
        HVGs <- as.numeric(n_genes)
        result_df <- data.frame(HVGs = HVGs)
        result_df$celltype <- celltype2
        result_df$condition <- par_status
        result_df$region <- par_brain_region
        df_list[[celltype2]] <- result_df
        }

        final_df_All_FTLD_BA9 <- do.call(rbind, df_list) 
    ##

    
    ###########################
    ## merge and plot BA9
    ###########################
    ## code
        BA9_bind_all <- rbind(final_df_SALS_BA9, final_df_C9ALS_BA9, final_df_SFTLD_BA9, final_df_C9FTLD_BA9, final_df_All_ALS_BA9, final_df_All_FTLD_BA9)

        ## set factor levels
        BA9_bind_all$condition <- factor(BA9_bind_all$condition, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All ALS", "All FTLD"))
        BA9_bind_all$celltype <- factor(BA9_bind_all$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))


        ggplot(BA9_bind_all, aes(x = condition, y = HVGs, fill = condition, colour = condition)) + 
        theme_bw() + 
        geom_bar(stat = 'identity') +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~celltype, scales = "free_x") +
        ylab("n HVGs") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue"))
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA9_n_sig_genes_0.pdf', height = 2, width = 13)
    ##

##


################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ Optimal: LOO workflow 0 delta threshold: Accuracy
## code

    ## BA4
    
    ###########################
    ## SALS BA4
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'SALS'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_SALS_BA4 <- do.call(rbind, df_list)
    ##

    ###########################
    ## C9ALS BA4
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'C9ALS'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_C9ALS_BA4 <- do.call(rbind, df_list)
    ##

    ###########################
    ## SFTLD BA4
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'SFTLD'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_SFTLD_BA4 <- do.call(rbind, df_list)
    ##

    ###########################
    ## C9FTLD BA4
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'C9FTLD'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_C9FTLD_BA4 <- do.call(rbind, df_list)
    ##

    ###########################
    ## All ALS BA4
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'All ALS'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/balanced_optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_All_ALS_BA4 <- do.call(rbind, df_list)
    ##

    ###########################
    ## All FTLD BA4
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'All FTLD'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/balanced_optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_All_FTLD_BA4 <- do.call(rbind, df_list)
    ##

    ###########################
    ## merge and plot BA4
    ###########################
    ## code for barplot
        BA4_bind_all <- rbind(final_SALS_BA4, final_C9ALS_BA4, final_SFTLD_BA4, final_C9FTLD_BA4, final_All_ALS_BA4, final_All_FTLD_BA4)
        write.csv(BA4_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_LOO_optimal.csv')
        
        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, group) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy),        
                .groups = "drop"               
            )

        ## Calculate the mean and standard deviation
        summary_stats <- BA4_bind_all %>%
            group_by(celltype, group) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## set factor levels
        summary_stats$group <- factor(summary_stats$group, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All ALS", "All FTLD"))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC' ))


        ## Plot
        ggplot(summary_stats, aes(x = group, y = median_accuracy, fill = group, colour = group)) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~celltype, scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple",  "darkred", "darkblue")) #+
        #coord_cartesian(ylim=c(0.5,1.0)) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA4_LOO_balanced_acc_barplot_0.pdf', height = 2, width = 13)
    ##

    ## BA9
    
    ###########################
    ## SALS BA9
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'SALS'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_SALS_BA9 <- do.call(rbind, df_list)
    ##

    ###########################
    ## C9ALS BA9
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'C9ALS'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_C9ALS_BA9 <- do.call(rbind, df_list)
    ##

    ###########################
    ## SFTLD BA9
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'SFTLD'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_SFTLD_BA9 <- do.call(rbind, df_list)
    ##

    ###########################
    ## C9FTLD BA9
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'C9FTLD'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_C9FTLD_BA9 <- do.call(rbind, df_list)
    ##

    ###########################
    ## All ALS BA9
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'All ALS'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/balanced_optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_All_ALS_BA9 <- do.call(rbind, df_list)
    ##

    ###########################
    ## All FTLD BA9
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'All FTLD'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/balanced_optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_All_FTLD_BA9 <- do.call(rbind, df_list)
    ##

    ###########################
    ## merge and plot BA9
    ###########################
    ## code for barplot
        BA9_bind_all <- rbind(final_SALS_BA9, final_C9ALS_BA9, final_SFTLD_BA9, final_C9FTLD_BA9, final_All_ALS_BA9, final_All_FTLD_BA9)
        write.csv(BA9_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_LOO_optimal.csv')
        
        ## Calculate sample wise mean
        BA9_bind_all <- BA9_bind_all %>%
            group_by(donor, celltype, group) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy),        
                .groups = "drop"               
            )

        ## Calculate the mean and standard deviation
        summary_stats <- BA9_bind_all %>%
            group_by(celltype, group) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## set factor levels
        summary_stats$group <- factor(summary_stats$group, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All ALS", "All FTLD"))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC' ))


        ## Plot
        ggplot(summary_stats, aes(x = group, y = median_accuracy, fill = group, colour = group)) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~celltype, scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple",  "darkred", "darkblue")) #+
        #coord_cartesian(ylim=c(0.5,1.0)) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA9_LOO_balanced_acc_barplot_0.pdf', height = 2, width = 13)
    ##

    
##

################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ Optimal: LOO workflow 0 delta threshold ABS scaling: Accuracy 
## code

    ## BA4
    
    ###########################
    ## SALS BA4
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'SALS'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_',status,'_',brain_region,'_',i,'_abs_scaling_narval_2.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_SALS_BA4 <- do.call(rbind, df_list)
    ##

    ###########################
    ## C9ALS BA4
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'C9ALS'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_',status,'_',brain_region,'_',i,'_abs_scaling_narval_2.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_C9ALS_BA4 <- do.call(rbind, df_list)
    ##

    ###########################
    ## SFTLD BA4
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'SFTLD'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_',status,'_',brain_region,'_',i,'_abs_scaling_narval_2.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_SFTLD_BA4 <- do.call(rbind, df_list)
    ##

    ###########################
    ## C9FTLD BA4
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'C9FTLD'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_',status,'_',brain_region,'_',i,'_abs_scaling_narval_2.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_C9FTLD_BA4 <- do.call(rbind, df_list)
    ##

    ###########################
    ## All ALS BA4
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'All ALS'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/balanced_optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_',status,'_',brain_region,'_',i,'_abs_scaling_narval_2.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_All_ALS_BA4 <- do.call(rbind, df_list)
    ##

    ###########################
    ## All FTLD BA4
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'All FTLD'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/balanced_optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_',status,'_',brain_region,'_',i,'_abs_scaling_narval_2.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_All_FTLD_BA4 <- do.call(rbind, df_list)
    ##

    ###########################
    ## merge and plot BA4
    ###########################
    ## code for barplot
        BA4_bind_all <- rbind(final_SALS_BA4, final_C9ALS_BA4, final_SFTLD_BA4, final_C9FTLD_BA4, final_All_ALS_BA4, final_All_FTLD_BA4)
        write.csv(BA4_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_LOO_optimal_ABS_scaling.csv')
        
        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, group) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy),        
                .groups = "drop"               
            )

        ## Calculate the mean and standard deviation
        summary_stats <- BA4_bind_all %>%
            group_by(celltype, group) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## set factor levels
        summary_stats$group <- factor(summary_stats$group, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All ALS", "All FTLD"))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC' ))


        ## Plot
        ggplot(summary_stats, aes(x = group, y = median_accuracy, fill = group, colour = group)) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~celltype, scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) #+
        #coord_cartesian(ylim=c(0.5,1.0)) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA4_LOO_balanced_acc_barplot_0_abs_scaling.pdf', height = 2, width = 13)
    ##

    ## BA9
    
    ###########################
    ## SALS BA9
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'SALS'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_',status,'_',brain_region,'_',i,'_abs_scaling_narval_2.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_SALS_BA9 <- do.call(rbind, df_list)
    ##

    ###########################
    ## C9ALS BA9
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'C9ALS'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_',status,'_',brain_region,'_',i,'_abs_scaling_narval_2.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_C9ALS_BA9 <- do.call(rbind, df_list)
    ##

    ###########################
    ## SFTLD BA9
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'SFTLD'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_',status,'_',brain_region,'_',i,'_abs_scaling_narval_2.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_SFTLD_BA9 <- do.call(rbind, df_list)
    ##

    ###########################
    ## C9FTLD BA9
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'C9FTLD'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_',status,'_',brain_region,'_',i,'_abs_scaling_narval_2.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_C9FTLD_BA9 <- do.call(rbind, df_list)
    ##

    ###########################
    ## All ALS BA9
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'All ALS'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/balanced_optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_',status,'_',brain_region,'_',i,'_abs_scaling_narval_2.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_All_ALS_BA9 <- do.call(rbind, df_list)
    ##

    ###########################
    ## All FTLD BA9
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'All FTLD'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/balanced_optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_',status,'_',brain_region,'_',i,'_abs_scaling_narval_2.csv')

        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_All_FTLD_BA9 <- do.call(rbind, df_list)
    ##

    ###########################
    ## merge and plot BA9
    ###########################
    ## code for barplot
        BA9_bind_all <- rbind(final_SALS_BA9, final_C9ALS_BA9, final_SFTLD_BA9, final_C9FTLD_BA9, final_All_ALS_BA9, final_All_FTLD_BA9)
        write.csv(BA9_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_LOO_optimal_ABS_scaling.csv')
        
        ## Calculate sample wise mean
        BA9_bind_all <- BA9_bind_all %>%
            group_by(donor, celltype, group) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy),        
                .groups = "drop"               
            )

        ## Calculate the mean and standard deviation
        summary_stats <- BA9_bind_all %>%
            group_by(celltype, group) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## set factor levels
        summary_stats$group <- factor(summary_stats$group, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All ALS", "All FTLD"))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC' ))


        ## Plot
        ggplot(summary_stats, aes(x = group, y = median_accuracy, fill = group, colour = group)) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~celltype, scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) #+
        #coord_cartesian(ylim=c(0.5,1.0)) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA9_LOO_balanced_acc_barplot_0_abs_scaling.pdf', height = 2, width = 13)
    ##

    

##

################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ Optimal: LOO workflow 0 delta threshold ABS scaling KNN Bayesian: Accuracy
## code
    ## BA4
    
    ###########################
    ## SALS BA4
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'SALS'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_',status,'_',brain_region,'_',i,'_abs_scaling_sweep_confidence_threshold_narval_2.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_SALS_BA4 <- do.call(rbind, df_list)
    ##

    ###########################
    ## C9ALS BA4
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'C9ALS'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_',status,'_',brain_region,'_',i,'_abs_scaling_sweep_confidence_threshold_narval_2.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_C9ALS_BA4 <- do.call(rbind, df_list)
    ##

    ###########################
    ## SFTLD BA4
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'SFTLD'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_',status,'_',brain_region,'_',i,'_abs_scaling_sweep_confidence_threshold_narval_2.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_SFTLD_BA4 <- do.call(rbind, df_list)
    ##

    ###########################
    ## C9FTLD BA4
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'C9FTLD'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_',status,'_',brain_region,'_',i,'_abs_scaling_sweep_confidence_threshold_narval_2.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_C9FTLD_BA4 <- do.call(rbind, df_list)
    ##

    ###########################
    ## All ALS BA4
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'All ALS'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/balanced_optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_',status,'_',brain_region,'_',i,'_abs_scaling_sweep_confidence_threshold_narval_2.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_All_ALS_BA4 <- do.call(rbind, df_list)
    ##

    ###########################
    ## All FTLD BA4
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'All FTLD'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/balanced_optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_',status,'_',brain_region,'_',i,'_abs_scaling_sweep_confidence_threshold_narval_2.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_All_FTLD_BA4 <- do.call(rbind, df_list)
    ##

    ###########################
    ## merge and plot BA4
    ###########################
    ## code for barplot
        BA4_bind_all <- rbind(final_SALS_BA4, final_C9ALS_BA4, final_SFTLD_BA4, final_C9FTLD_BA4, final_All_ALS_BA4, final_All_FTLD_BA4)
        write.csv(BA4_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_LOO_optimal_ABS_scaling_KNN_Bayesian.csv')
        
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_LOO_optimal_ABS_scaling_KNN_Bayesian.csv')
        
        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, group, kNN_thresh) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )

        ## Calculate the mean and standard deviation
        summary_stats <- BA4_bind_all %>%
            group_by(celltype, group, kNN_thresh) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## set factor levels
        summary_stats$group <- factor(summary_stats$group, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All ALS", "All FTLD"))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))


        ## Plot
        ggplot(summary_stats, aes(x = group, y = median_accuracy, fill = group, colour = group)) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(kNN_thresh~celltype, scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA4_LOO_balanced_acc_barplot_optimal_abs_scaling_KNN_bayesian.pdf', height = 8, width = 13)

        ## Plot only Bayes 0.9 thresh
        summary_stats_lim <- subset(summary_stats, kNN_thresh == 0.9 )
        ggplot(summary_stats_lim, aes(x = group, y = median_accuracy, fill = group, colour = group)) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~celltype, scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA4_LOO_balanced_acc_barplot_optimal_abs_scaling_KNN_bayesian_09.pdf', height = 2, width = 13)

    ##

    ###########################
    ## merge and plot BA4 at subject level
    ###########################
    ## code for barplot        
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_LOO_optimal_ABS_scaling_KNN_Bayesian.csv')
        BA4_bind_all_lim <- subset(BA4_bind_all, kNN_thresh == 0.9)
        
        ## Calculate sample wise mean
        BA4_bind_all_lim <- BA4_bind_all_lim %>%
            group_by(donor, celltype, group, kNN_thresh) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )
        BA4_bind_all_lim$sample_level <- 0
        BA4_bind_all_lim$sample_level[BA4_bind_all_lim$sample_mean_accuracy >= 0.5] <- 1

        summary_stats <-  BA4_bind_all_lim %>%
            group_by(celltype, group) %>%
            summarise(
                prop_sample_level_1 = mean(sample_level == 1),
                .groups = 'drop'
            )



        ## set factor levels
        summary_stats$group <- factor(summary_stats$group, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All ALS", "All FTLD"))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))


        ## Plot
        ggplot(summary_stats, aes(x = group, y = prop_sample_level_1, fill = group, colour = group)) + 
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~celltype, scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 2, width = 13)




        

    ##

    ## BA9
    
    ###########################
    ## SALS BA9
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'SALS'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_',status,'_',brain_region,'_',i,'_abs_scaling_sweep_confidence_threshold_narval_2.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_SALS_BA9 <- do.call(rbind, df_list)
    ##

    ###########################
    ## C9ALS BA9
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'C9ALS'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_',status,'_',brain_region,'_',i,'_abs_scaling_sweep_confidence_threshold_narval_2.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_C9ALS_BA9 <- do.call(rbind, df_list)
    ##

    ###########################
    ## SFTLD BA9
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'SFTLD'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_',status,'_',brain_region,'_',i,'_abs_scaling_sweep_confidence_threshold_narval_2.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_SFTLD_BA9 <- do.call(rbind, df_list)
    ##

    ###########################
    ## C9FTLD BA9
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'C9FTLD'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_',status,'_',brain_region,'_',i,'_abs_scaling_sweep_confidence_threshold_narval_2.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_C9FTLD_BA9 <- do.call(rbind, df_list)
    ##

    ###########################
    ## All ALS BA9
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'All ALS'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/balanced_optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_',status,'_',brain_region,'_',i,'_abs_scaling_sweep_confidence_threshold_narval_2.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_All_ALS_BA9 <- do.call(rbind, df_list)
    ##

    ###########################
    ## All FTLD BA9
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'All FTLD'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/balanced_optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_',status,'_',brain_region,'_',i,'_abs_scaling_sweep_confidence_threshold_narval_2.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_All_FTLD_BA9 <- do.call(rbind, df_list)
    ##

    ###########################
    ## merge and plot BA9
    ###########################
    ## code for barplot
        BA9_bind_all <- rbind(final_SALS_BA9, final_C9ALS_BA9, final_SFTLD_BA9, final_C9FTLD_BA9, final_All_ALS_BA9, final_All_FTLD_BA9)
        write.csv(BA9_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_LOO_optimal_ABS_scaling_KNN_Bayesian.csv')
        
        BA9_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_LOO_optimal_ABS_scaling_KNN_Bayesian.csv')
        
        ## Calculate sample wise mean
        BA9_bind_all <- BA9_bind_all %>%
            group_by(donor, celltype, group, kNN_thresh) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )

        ## Calculate the mean and standard deviation
        summary_stats <- BA9_bind_all %>%
            group_by(celltype, group, kNN_thresh) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## set factor levels
        summary_stats$group <- factor(summary_stats$group, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All ALS", "All FTLD"))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))


        ## Plot
        ggplot(summary_stats, aes(x = group, y = median_accuracy, fill = group, colour = group)) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(kNN_thresh~celltype, scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA9_LOO_balanced_acc_barplot_optimal_abs_scaling_KNN_bayesian.pdf', height = 8, width = 13)

        ## Plot only Bayes 0.9 thresh
        summary_stats_lim <- subset(summary_stats, kNN_thresh == 0.9 )
        ggplot(summary_stats_lim, aes(x = group, y = median_accuracy, fill = group, colour = group)) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~celltype, scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA9_LOO_balanced_acc_barplot_optimal_abs_scaling_KNN_bayesian_09.pdf', height = 2, width = 13)


    ##

    ###########################
    ## merge and plot BA9 at subject level
    ###########################
    ## code for barplot
        
        BA9_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_LOO_optimal_ABS_scaling_KNN_Bayesian.csv')
        BA9_bind_all_lim <- subset(BA9_bind_all, kNN_thresh == 0.9)
        
        ## Calculate sample wise mean
        BA9_bind_all_lim <- BA9_bind_all_lim %>%
            group_by(donor, celltype, group, kNN_thresh) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )
        BA9_bind_all_lim$sample_level <- 0
        BA9_bind_all_lim$sample_level[BA9_bind_all_lim$sample_mean_accuracy >= 0.5] <- 1

        BA9_bind_all_lim <- data.frame(BA9_bind_all_lim)
        
        summary_stats <-  BA9_bind_all_lim %>%
            group_by(celltype, group) %>%
            summarise(
                prop_sample_level_1 = mean(sample_level == 1),
                .groups = 'drop'
            )
        summary_stats <- write.csv(BA9_bind_all_lim, '/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.csv')

        ## set factor levels
        summary_stats$group <- factor(summary_stats$group, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All ALS", "All FTLD"))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))


        ## Plot
        ggplot(summary_stats, aes(x = group, y = prop_sample_level_1, fill = group, colour = group)) + 
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~celltype, scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 2, width = 13)

        mean(summary_stats$prop_sample_level_1)




        

    ##
    
    
    ###########################
    ## BA4 and BA9 0.9 Bayes plots for Sali pres. 
    ###########################
    ## Code
        ###################
        ## BA4
        ###################
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_LOO_optimal_ABS_scaling_KNN_Bayesian.csv')
        
        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, group, kNN_thresh) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )

        ## Calculate the mean and standard deviation
        summary_stats_BA4 <- BA4_bind_all %>%
            group_by(celltype, group, kNN_thresh) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## set factor levels
        summary_stats_BA4$group <- factor(summary_stats_BA4$group, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All ALS", "All FTLD"))
        summary_stats_BA4$celltype <- factor(summary_stats_BA4$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))

        ## Plot only Bayes 0.9 thresh
        summary_stats_BA4 <- subset(summary_stats_BA4, kNN_thresh == 0.9 )
        summary_stats_BA4$region <- "BA4"

        ###################
        ## BA9
        ###################
        BA9_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_LOO_optimal_ABS_scaling_KNN_Bayesian.csv')
        
        ## Calculate sample wise mean
        BA9_bind_all <- BA9_bind_all %>%
            group_by(donor, celltype, group, kNN_thresh) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )

        ## Calculate the mean and standard deviation
        summary_stats_BA9 <- BA9_bind_all %>%
            group_by(celltype, group, kNN_thresh) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## set factor levels
        summary_stats_BA9$group <- factor(summary_stats_BA9$group, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All ALS", "All FTLD"))
        summary_stats_BA9$celltype <- factor(summary_stats_BA9$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))

        ## Plot only Bayes 0.9 thresh
        summary_stats_BA9 <- subset(summary_stats_BA9, kNN_thresh == 0.9 )
        summary_stats_BA9$region <- "BA9"

        ###################
        ## Merge and plot
        ###################
        bind_all_test <- rbind(summary_stats_BA4, summary_stats_BA9)

        ## Calculate mean of groups across brain regions
        summary_stats_all <- bind_all_test %>%
            group_by(group, region) %>%
            summarise(
                mean_accuracy = mean(median_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(median_accuracy),                  # Calculate the standard deviation of X0
                .groups = "drop"                 # Drop the grouping after summarising
            )
        
        ## Plot
        ggplot(summary_stats_all, aes(x = group, y = mean_accuracy, fill = group, colour = group, label = round(mean_accuracy, 2) )) + 
        geom_errorbar(aes(ymin=mean_accuracy - sd_accuracy, ymax=mean_accuracy + sd_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_text(colour = "white", size = 3, vjust = 1) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(region ~ ., scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 4, width = 3)


        ## By major celltype group
        bind_all_test <- data.frame(bind_all_test)
        bind_all_test$major_group[bind_all_test$celltype == "L2_L3"] <- "Excitatory"
        bind_all_test$major_group[bind_all_test$celltype == "L3_L5"] <- "Excitatory"
        bind_all_test$major_group[bind_all_test$celltype == "L4_L5"] <- "Excitatory"
        bind_all_test$major_group[bind_all_test$celltype == "L4_L6"] <- "Excitatory"
        bind_all_test$major_group[bind_all_test$celltype == "L5"] <- "Excitatory"
        bind_all_test$major_group[bind_all_test$celltype == "L5_L6"] <- "Excitatory"
        bind_all_test$major_group[bind_all_test$celltype == "L6"] <- "Excitatory"

        bind_all_test$major_group[bind_all_test$celltype == "5HT3aR"] <- "Inhibitory"
        bind_all_test$major_group[bind_all_test$celltype == "PV"] <- "Inhibitory"
        bind_all_test$major_group[bind_all_test$celltype == "Rosehip"] <- "Inhibitory"
        bind_all_test$major_group[bind_all_test$celltype == "SOM"] <- "Inhibitory"

        bind_all_test$major_group[bind_all_test$celltype == "Astro"] <- "non-neuronal"
        bind_all_test$major_group[bind_all_test$celltype == "Endo"] <- "non-neuronal"
        bind_all_test$major_group[bind_all_test$celltype == "Fibro"] <- "non-neuronal"
        bind_all_test$major_group[bind_all_test$celltype == "Micro"] <- "non-neuronal"
        bind_all_test$major_group[bind_all_test$celltype == "Mural"] <- "non-neuronal"
        bind_all_test$major_group[bind_all_test$celltype == "Oligo"] <- "non-neuronal"
        bind_all_test$major_group[bind_all_test$celltype == "OPC"] <- "non-neuronal"

        ## Calculate mean of groups across brain regions
        summary_stats_all <- bind_all_test %>%
            group_by(group, region, major_group) %>%
            summarise(
                mean_accuracy = mean(median_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(median_accuracy),                  # Calculate the standard deviation of X0
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## Plot
        ggplot(summary_stats_all, aes(x = group, y = mean_accuracy, fill = group, colour = group, label = round(mean_accuracy, 2) )) + 
        geom_errorbar(aes(ymin=mean_accuracy - sd_accuracy, ymax=mean_accuracy + sd_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_text(colour = "white", size = 3, vjust = 1) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(region ~ major_group, scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 4, width = 6)

        ## Calculate mean of groups across brain regions
        summary_stats_all <- bind_all_test %>%
            group_by(major_group) %>%
            summarise(
                mean_accuracy = mean(median_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(median_accuracy),                  # Calculate the standard deviation of X0
                .groups = "drop"                 # Drop the grouping after summarising
            )
        
    ##


##



############################################################################################################################################################################### Compare methods



################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ HVGs: compare all methods
## code

    ## BA4

    ##########################
    ## BA4 cell type specific
    ##########################
    ## code
        BA4_LOO <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_LOO_HVGs.csv')
        BA4_LOO <- BA4_LOO %>%
            group_by(donor, celltype, group) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy),        
                .groups = "drop"               
            )
        BA4_LOO <- BA4_LOO %>%
            group_by(celltype) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )
        BA4_LOO$method <- "Base"

        ## AbsScaling
        BA4_LOO_abs <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_LOO_HVGs_ABS_scaling.csv')
        BA4_LOO_abs <- BA4_LOO_abs %>%
            group_by(donor, celltype, group) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy),        
                .groups = "drop"               
            )
        BA4_LOO_abs <- BA4_LOO_abs %>%
            group_by(celltype) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )
        BA4_LOO_abs$method <- "AbsScaling"

        ## AbsScaling + KNN Bayesian
        BA4_LOO_abs_KNN <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_LOO_HVGs_ABS_scaling_KNN_Bayesian.csv')
        BA4_LOO_abs_KNN <- BA4_LOO_abs_KNN %>%
            group_by(donor, celltype, group, kNN_thresh) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )
        BA4_LOO_abs_KNN <- BA4_LOO_abs_KNN %>%
            group_by(celltype, kNN_thresh) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        BA4_LOO_abs_KNN$method <- "AbsScaling + KNN Bayesian"
        BA4_LOO_abs_KNN$method <- paste0(BA4_LOO_abs_KNN$method, " (", BA4_LOO_abs_KNN$kNN_thresh, ")" )
        BA4_LOO_abs_KNN <- BA4_LOO_abs_KNN %>% dplyr::select(celltype, median_accuracy, mean_accuracy, sd_accuracy, mad_accuracy, method)

        ## Bind all
        bind_total <- rbind(BA4_LOO, BA4_LOO_abs, BA4_LOO_abs_KNN)

        ## set factor levels
        bind_total$method <- factor(bind_total$method, levels = c("Base", "AbsScaling", "AbsScaling + KNN Bayesian (0.8)",  "AbsScaling + KNN Bayesian (0.85)", "AbsScaling + KNN Bayesian (0.9)", "AbsScaling + KNN Bayesian (0.95)", "AbsScaling + KNN Bayesian (0.99)"   ))
        bind_total$celltype <- factor(bind_total$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'  ))

        ## Plot
        ggplot(bind_total, aes(x = method, y = median_accuracy, fill = method, colour = method)) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_wrap(~celltype, scales = "free_x", ncol = 18) +
        ylab("Accuracy") 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA4_LOO_cell_type_specific_all_methods_bar_HVGs.pdf', height = 4, width = 17)
    ##

    ##########################
    ## BA4 across all cell types
    ##########################
    ## code
        BA4_LOO <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_LOO_HVGs.csv')
        BA4_LOO <- BA4_LOO %>%
            group_by(donor, celltype, group) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy),        
                .groups = "drop"               
            )
        BA4_LOO <- BA4_LOO %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )
        BA4_LOO$method <- "Base"

        ## AbsScaling
        BA4_LOO_abs <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_LOO_HVGs_ABS_scaling.csv')
        BA4_LOO_abs <- BA4_LOO_abs %>%
            group_by(donor, celltype, group) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy),        
                .groups = "drop"               
            )
        BA4_LOO_abs <- BA4_LOO_abs %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )
        BA4_LOO_abs$method <- "AbsScaling"

        ## AbsScaling + KNN Bayesian
        BA4_LOO_abs_KNN <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_LOO_HVGs_ABS_scaling_KNN_Bayesian.csv')
        BA4_LOO_abs_KNN <- BA4_LOO_abs_KNN %>%
            group_by(donor, celltype, group, kNN_thresh) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )
        BA4_LOO_abs_KNN <- BA4_LOO_abs_KNN %>%
            group_by(kNN_thresh) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        BA4_LOO_abs_KNN$method <- "AbsScaling + KNN Bayesian"
        BA4_LOO_abs_KNN$method <- paste0(BA4_LOO_abs_KNN$method, " (", BA4_LOO_abs_KNN$kNN_thresh, ")" )
        BA4_LOO_abs_KNN <- BA4_LOO_abs_KNN %>% dplyr::select(median_accuracy, mean_accuracy, sd_accuracy, mad_accuracy, method)

        ## Bind all
        bind_total <- rbind(BA4_LOO, BA4_LOO_abs, BA4_LOO_abs_KNN)

        ## set factor levels
        bind_total$method <- factor(bind_total$method, levels = c("Base", "AbsScaling", "AbsScaling + KNN Bayesian (0.8)",  "AbsScaling + KNN Bayesian (0.85)", "AbsScaling + KNN Bayesian (0.9)", "AbsScaling + KNN Bayesian (0.95)", "AbsScaling + KNN Bayesian (0.99)"   ))

        ## Plot
        ggplot(bind_total, aes(x = method, y = median_accuracy, fill = method, colour = method, label = round(median_accuracy, 3))) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_text(colour = "black", size = 1.5, vjust = 1) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        ylab("Accuracy") 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA4_LOO_all_methods_bar_HVGs.pdf', height = 3.5, width = 2)
    ##

    ## BA9

    ##########################
    ## BA9 cell type specific
    ##########################
    ## code
        BA9_LOO <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_LOO_HVGs.csv')
        BA9_LOO <- BA9_LOO %>%
            group_by(donor, celltype, group) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy),        
                .groups = "drop"               
            )
        BA9_LOO <- BA9_LOO %>%
            group_by(celltype) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )
        BA9_LOO$method <- "Base"

        ## AbsScaling
        BA9_LOO_abs <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_LOO_HVGs_ABS_scaling.csv')
        BA9_LOO_abs <- BA9_LOO_abs %>%
            group_by(donor, celltype, group) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy),        
                .groups = "drop"               
            )
        BA9_LOO_abs <- BA9_LOO_abs %>%
            group_by(celltype) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )
        BA9_LOO_abs$method <- "AbsScaling"

        ## AbsScaling + KNN Bayesian
        BA9_LOO_abs_KNN <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_LOO_HVGs_ABS_scaling_KNN_Bayesian.csv')
        BA9_LOO_abs_KNN <- BA9_LOO_abs_KNN %>%
            group_by(donor, celltype, group, kNN_thresh) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )
        BA9_LOO_abs_KNN <- BA9_LOO_abs_KNN %>%
            group_by(celltype, kNN_thresh) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        BA9_LOO_abs_KNN$method <- "AbsScaling + KNN Bayesian"
        BA9_LOO_abs_KNN$method <- paste0(BA9_LOO_abs_KNN$method, " (", BA9_LOO_abs_KNN$kNN_thresh, ")" )
        BA9_LOO_abs_KNN <- BA9_LOO_abs_KNN %>% dplyr::select(celltype, median_accuracy, mean_accuracy, sd_accuracy, mad_accuracy, method)

        ## Bind all
        bind_total <- rbind(BA9_LOO, BA9_LOO_abs, BA9_LOO_abs_KNN)

        ## set factor levels
        bind_total$method <- factor(bind_total$method, levels = c("Base", "AbsScaling", "AbsScaling + KNN Bayesian (0.8)",  "AbsScaling + KNN Bayesian (0.85)", "AbsScaling + KNN Bayesian (0.9)", "AbsScaling + KNN Bayesian (0.95)", "AbsScaling + KNN Bayesian (0.99)"   ))
        bind_total$celltype <- factor(bind_total$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'  ))

        ## Plot
        ggplot(bind_total, aes(x = method, y = median_accuracy, fill = method, colour = method)) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_wrap(~celltype, scales = "free_x", ncol = 18) +
        ylab("Accuracy") 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA9_LOO_cell_type_specific_all_methods_bar_HVGs.pdf', height = 4, width = 17)
    ##

    ##########################
    ## BA9 across all cell types
    ##########################
    ## code
        BA9_LOO <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_LOO_HVGs.csv')
        BA9_LOO <- BA9_LOO %>%
            group_by(donor, celltype, group) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy),        
                .groups = "drop"               
            )
        BA9_LOO <- BA9_LOO %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )
        BA9_LOO$method <- "Base"

        ## AbsScaling
        BA9_LOO_abs <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_LOO_HVGs_ABS_scaling.csv')
        BA9_LOO_abs <- BA9_LOO_abs %>%
            group_by(donor, celltype, group) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy),        
                .groups = "drop"               
            )
        BA9_LOO_abs <- BA9_LOO_abs %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )
        BA9_LOO_abs$method <- "AbsScaling"

        ## AbsScaling + KNN Bayesian
        BA9_LOO_abs_KNN <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_LOO_HVGs_ABS_scaling_KNN_Bayesian.csv')
        BA9_LOO_abs_KNN <- BA9_LOO_abs_KNN %>%
            group_by(donor, celltype, group, kNN_thresh) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )
        BA9_LOO_abs_KNN <- BA9_LOO_abs_KNN %>%
            group_by(kNN_thresh) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        BA9_LOO_abs_KNN$method <- "AbsScaling + KNN Bayesian"
        BA9_LOO_abs_KNN$method <- paste0(BA9_LOO_abs_KNN$method, " (", BA9_LOO_abs_KNN$kNN_thresh, ")" )
        BA9_LOO_abs_KNN <- BA9_LOO_abs_KNN %>% dplyr::select(median_accuracy, mean_accuracy, sd_accuracy, mad_accuracy, method)

        ## Bind all
        bind_total <- rbind(BA9_LOO, BA9_LOO_abs, BA9_LOO_abs_KNN)

        ## set factor levels
        bind_total$method <- factor(bind_total$method, levels = c("Base", "AbsScaling", "AbsScaling + KNN Bayesian (0.8)",  "AbsScaling + KNN Bayesian (0.85)", "AbsScaling + KNN Bayesian (0.9)", "AbsScaling + KNN Bayesian (0.95)", "AbsScaling + KNN Bayesian (0.99)"   ))

        ## Plot
        ggplot(bind_total, aes(x = method, y = median_accuracy, fill = method, colour = method, label = round(median_accuracy, 3))) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_text(colour = "black", size = 1.5, vjust = 1) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        scale_y_continuous(limits = c(0,1)) +
        ylab("Accuracy") 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA9_LOO_all_methods_bar_HVGs.pdf', height = 3.5, width = 2)
    ##


##


################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ Optimal: compare all methods
## code

    ## BA4

    ##########################
    ## BA4 cell type specific
    ##########################
    ## code
        BA4_LOO <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_LOO_optimal.csv')
        BA4_LOO <- BA4_LOO %>%
            group_by(donor, celltype, group) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy),        
                .groups = "drop"               
            )
        BA4_LOO <- BA4_LOO %>%
            group_by(celltype) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )
        BA4_LOO$method <- "Base"

        ## AbsScaling
        BA4_LOO_abs <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_LOO_optimal_ABS_scaling.csv')
        BA4_LOO_abs <- BA4_LOO_abs %>%
            group_by(donor, celltype, group) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy),        
                .groups = "drop"               
            )
        BA4_LOO_abs <- BA4_LOO_abs %>%
            group_by(celltype) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )
        BA4_LOO_abs$method <- "AbsScaling"

        ## AbsScaling + KNN Bayesian
        BA4_LOO_abs_KNN <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_LOO_optimal_ABS_scaling_KNN_Bayesian.csv')
        BA4_LOO_abs_KNN <- BA4_LOO_abs_KNN %>%
            group_by(donor, celltype, group, kNN_thresh) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )
        BA4_LOO_abs_KNN <- BA4_LOO_abs_KNN %>%
            group_by(celltype, kNN_thresh) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        BA4_LOO_abs_KNN$method <- "AbsScaling + KNN Bayesian"
        BA4_LOO_abs_KNN$method <- paste0(BA4_LOO_abs_KNN$method, " (", BA4_LOO_abs_KNN$kNN_thresh, ")" )
        BA4_LOO_abs_KNN <- BA4_LOO_abs_KNN %>% dplyr::select(celltype, median_accuracy, mean_accuracy, sd_accuracy, mad_accuracy, method)

        ## Bind all
        bind_total <- rbind(BA4_LOO, BA4_LOO_abs, BA4_LOO_abs_KNN)

        ## set factor levels
        bind_total$method <- factor(bind_total$method, levels = c("Base", "AbsScaling", "AbsScaling + KNN Bayesian (0.8)",  "AbsScaling + KNN Bayesian (0.85)", "AbsScaling + KNN Bayesian (0.9)", "AbsScaling + KNN Bayesian (0.95)", "AbsScaling + KNN Bayesian (0.99)"   ))
        bind_total$celltype <- factor(bind_total$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'  ))

        ## Plot
        ggplot(bind_total, aes(x = method, y = median_accuracy, fill = method, colour = method)) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_wrap(~celltype, scales = "free_x", ncol = 18) +
        ylab("Accuracy") 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA4_LOO_cell_type_specific_all_methods_bar.pdf', height = 4, width = 17)
    ##

    ##########################
    ## BA4 across all cell types
    ##########################
    ## code
        BA4_LOO <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_LOO_optimal.csv')
        BA4_LOO <- BA4_LOO %>%
            group_by(donor, celltype, group) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy),        
                .groups = "drop"               
            )
        BA4_LOO <- BA4_LOO %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )
        BA4_LOO$method <- "Base"

        ## AbsScaling
        BA4_LOO_abs <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_LOO_optimal_ABS_scaling.csv')
        BA4_LOO_abs <- BA4_LOO_abs %>%
            group_by(donor, celltype, group) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy),        
                .groups = "drop"               
            )
        BA4_LOO_abs <- BA4_LOO_abs %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )
        BA4_LOO_abs$method <- "AbsScaling"

        ## AbsScaling + KNN Bayesian
        BA4_LOO_abs_KNN <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_LOO_optimal_ABS_scaling_KNN_Bayesian.csv')
        BA4_LOO_abs_KNN <- BA4_LOO_abs_KNN %>%
            group_by(donor, celltype, group, kNN_thresh) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )
        BA4_LOO_abs_KNN <- BA4_LOO_abs_KNN %>%
            group_by(kNN_thresh) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        BA4_LOO_abs_KNN$method <- "AbsScaling + KNN Bayesian"
        BA4_LOO_abs_KNN$method <- paste0(BA4_LOO_abs_KNN$method, " (", BA4_LOO_abs_KNN$kNN_thresh, ")" )
        BA4_LOO_abs_KNN <- BA4_LOO_abs_KNN %>% dplyr::select(median_accuracy, mean_accuracy, sd_accuracy, mad_accuracy, method)

        ## Bind all
        bind_total <- rbind(BA4_LOO, BA4_LOO_abs, BA4_LOO_abs_KNN)

        ## set factor levels
        bind_total$method <- factor(bind_total$method, levels = c("Base", "AbsScaling", "AbsScaling + KNN Bayesian (0.8)",  "AbsScaling + KNN Bayesian (0.85)", "AbsScaling + KNN Bayesian (0.9)", "AbsScaling + KNN Bayesian (0.95)", "AbsScaling + KNN Bayesian (0.99)"   ))

        ## Plot
        ggplot(bind_total, aes(x = method, y = median_accuracy, fill = method, colour = method, label = round(median_accuracy, 3))) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_text(colour = "black", size = 1.5, vjust = 1) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        ylab("Accuracy") 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA4_LOO_all_methods_bar.pdf', height = 3.5, width = 2)
    ##

    ## BA9

    ##########################
    ## BA9 cell type specific
    ##########################
    ## code
        BA9_LOO <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_LOO_optimal.csv')
        BA9_LOO <- BA9_LOO %>%
            group_by(donor, celltype, group) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy),        
                .groups = "drop"               
            )
        BA9_LOO <- BA9_LOO %>%
            group_by(celltype) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )
        BA9_LOO$method <- "Base"

        ## AbsScaling
        BA9_LOO_abs <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_LOO_optimal_ABS_scaling.csv')
        BA9_LOO_abs <- BA9_LOO_abs %>%
            group_by(donor, celltype, group) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy),        
                .groups = "drop"               
            )
        BA9_LOO_abs <- BA9_LOO_abs %>%
            group_by(celltype) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )
        BA9_LOO_abs$method <- "AbsScaling"

        ## AbsScaling + KNN Bayesian
        BA9_LOO_abs_KNN <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_LOO_optimal_ABS_scaling_KNN_Bayesian.csv')
        BA9_LOO_abs_KNN <- BA9_LOO_abs_KNN %>%
            group_by(donor, celltype, group, kNN_thresh) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )
        BA9_LOO_abs_KNN <- BA9_LOO_abs_KNN %>%
            group_by(celltype, kNN_thresh) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        BA9_LOO_abs_KNN$method <- "AbsScaling + KNN Bayesian"
        BA9_LOO_abs_KNN$method <- paste0(BA9_LOO_abs_KNN$method, " (", BA9_LOO_abs_KNN$kNN_thresh, ")" )
        BA9_LOO_abs_KNN <- BA9_LOO_abs_KNN %>% dplyr::select(celltype, median_accuracy, mean_accuracy, sd_accuracy, mad_accuracy, method)

        ## Bind all
        bind_total <- rbind(BA9_LOO, BA9_LOO_abs, BA9_LOO_abs_KNN)

        ## set factor levels
        bind_total$method <- factor(bind_total$method, levels = c("Base", "AbsScaling", "AbsScaling + KNN Bayesian (0.8)",  "AbsScaling + KNN Bayesian (0.85)", "AbsScaling + KNN Bayesian (0.9)", "AbsScaling + KNN Bayesian (0.95)", "AbsScaling + KNN Bayesian (0.99)"   ))
        bind_total$celltype <- factor(bind_total$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'  ))

        ## Plot
        ggplot(bind_total, aes(x = method, y = median_accuracy, fill = method, colour = method)) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_wrap(~celltype, scales = "free_x", ncol = 18) +
        ylab("Accuracy") 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA9_LOO_cell_type_specific_all_methods_bar.pdf', height = 4, width = 17)
    ##

    ##########################
    ## BA9 across all cell types
    ##########################
    ## code
        BA9_LOO <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_LOO_optimal.csv')
        BA9_LOO <- BA9_LOO %>%
            group_by(donor, celltype, group) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy),        
                .groups = "drop"               
            )
        BA9_LOO <- BA9_LOO %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )
        BA9_LOO$method <- "Base"

        ## AbsScaling
        BA9_LOO_abs <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_LOO_optimal_ABS_scaling.csv')
        BA9_LOO_abs <- BA9_LOO_abs %>%
            group_by(donor, celltype, group) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy),        
                .groups = "drop"               
            )
        BA9_LOO_abs <- BA9_LOO_abs %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )
        BA9_LOO_abs$method <- "AbsScaling"

        ## AbsScaling + KNN Bayesian
        BA9_LOO_abs_KNN <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_LOO_optimal_ABS_scaling_KNN_Bayesian.csv')
        BA9_LOO_abs_KNN <- BA9_LOO_abs_KNN %>%
            group_by(donor, celltype, group, kNN_thresh) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )
        BA9_LOO_abs_KNN <- BA9_LOO_abs_KNN %>%
            group_by(kNN_thresh) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        BA9_LOO_abs_KNN$method <- "AbsScaling + KNN Bayesian"
        BA9_LOO_abs_KNN$method <- paste0(BA9_LOO_abs_KNN$method, " (", BA9_LOO_abs_KNN$kNN_thresh, ")" )
        BA9_LOO_abs_KNN <- BA9_LOO_abs_KNN %>% dplyr::select(median_accuracy, mean_accuracy, sd_accuracy, mad_accuracy, method)

        ## Bind all
        bind_total <- rbind(BA9_LOO, BA9_LOO_abs, BA9_LOO_abs_KNN)

        ## set factor levels
        bind_total$method <- factor(bind_total$method, levels = c("Base", "AbsScaling", "AbsScaling + KNN Bayesian (0.8)",  "AbsScaling + KNN Bayesian (0.85)", "AbsScaling + KNN Bayesian (0.9)", "AbsScaling + KNN Bayesian (0.95)", "AbsScaling + KNN Bayesian (0.99)"   ))

        ## Plot
        ggplot(bind_total, aes(x = method, y = median_accuracy, fill = method, colour = method, label = round(median_accuracy, 3))) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_text(colour = "black", size = 1.5, vjust = 1) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        ylab("Accuracy") 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA9_LOO_all_methods_bar.pdf', height = 3.5, width = 2)
    ##
##

write.csv(BA4_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_LOO_optimal.csv')
write.csv(BA9_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_LOO_optimal.csv')
write.csv(BA4_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_LOO_optimal_ABS_scaling.csv')
write.csv(BA9_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_LOO_optimal_ABS_scaling.csv')
write.csv(BA4_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_LOO_optimal_ABS_scaling_KNN_Bayesian.csv')
write.csv(BA4_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_LOO_optimal_ABS_scaling_KNN_Bayesian.csv')

############################################################################################################################################################################### ETM explore

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################ Topic-epoch sweep explore: reg - no bayesian -- THIS IS WHERE WE PRINT THE FILES WITH THE OPTIMAL TOPICS AND EPOCHS.
## code 
    ###########################
    ## SALS BA4
    ###########################
    ## code
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'SALS'
        brain_region <- 'BA4'
        topic_list <- c(100, 200, 300)
        epoch_list <- c(35, 50, 100, 200, 300, 400, 500)


        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
            print(i)
            for (topic in topic_list) {
                for (epoch in epoch_list) {
                    # Construct the file path for the current cell type
                    file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_',topic,'_',epoch,'_narval_2.csv')

                    # Check if the file exists
                    if (file.exists(file_path)) {
                        # Read the CSV file
                        current_df <- read.csv(file_path)
                        df_list[[paste0(i, '_',topic,'_',epoch)]] <- current_df
                    } else {
                        warning(paste("File does not exist:", file_path))
                    }
                }
            }
        }
        

        # Merge all dataframes into one
        final_SALS_BA4 <- do.call(rbind, df_list)

        ## Calculate sample wise mean
        final_SALS_BA4 <- final_SALS_BA4 %>%
            group_by(donor, celltype, group, n_genes, n_epochs) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy),        
                .groups = "drop"               
            )

        ## Calculate the mean and standard deviation
        summary_stats <- final_SALS_BA4 %>%
            group_by(celltype, group, n_genes, n_epochs) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## set factor levels
        #summary_stats$group <- factor(summary_stats$group, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD"))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC' ))


        ## Plot
        ggplot(summary_stats, aes(x = as.factor(n_epochs), y = median_accuracy, fill = as.factor(n_epochs), colour = as.factor(n_epochs))) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(n_genes~celltype, scales = "free_x") +
        ylab("Accuracy") 
        #coord_cartesian(ylim=c(0.5,1.0)) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/bar_epoch_sweep_SALS_BA4.pdf', height = 4, width = 13)

        ## Save the complete file
        write.csv(summary_stats, paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_summary_CombatSeq_',status,'_',brain_region,'_narval_2.csv'))

        ## Define the optimal topic-epoch combination by considering both mean and median accuracy
        summary_stats_optimal <- summary_stats %>%
            mutate(composite_score = (median_accuracy + mean_accuracy) / 2) %>%
            group_by(celltype) %>%
            arrange(desc(composite_score)) %>%
            slice(1) %>%
            ungroup()

        table(summary_stats_optimal$celltype)

        summary_stats_optimal <- summary_stats_optimal %>%
            group_by(celltype) %>%
            slice_min(sd_accuracy, n = 1, with_ties = FALSE) %>%
            ungroup()
        
        table(summary_stats_optimal$celltype)

        write.csv(summary_stats_optimal, paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_optimal_CombatSeq_',status,'_',brain_region,'_narval_2.csv'))
    ##

    ###########################
    ## SALS BA9
    ###########################
    ## code
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'SALS'
        brain_region <- 'BA9'
        topic_list <- c(100, 200, 300)
        epoch_list <- c(35, 50, 100, 200, 300, 400, 500)


        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
            print(i)
            for (topic in topic_list) {
                for (epoch in epoch_list) {
                    # Construct the file path for the current cell type
                    file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_',topic,'_',epoch,'_narval_2.csv')

                    # Check if the file exists
                    if (file.exists(file_path)) {
                        # Read the CSV file
                        current_df <- read.csv(file_path)
                        df_list[[paste0(i, '_',topic,'_',epoch)]] <- current_df
                    } else {
                        warning(paste("File does not exist:", file_path))
                    }
                }
            }
        }
        
        for (name in names(df_list)) {
        cat("----", name, "----\n")
        print(dim(df_list[[name]]))
        #print(colnames(df_list[[name]]))
        }

        for (name in names(df_list)) {
        df <- df_list[[name]]
        if (ncol(df) == 11) {
            cat("----", name, "----\n")
            print(dim(df))
            # Optionally print column names
            # print(colnames(df))
        }
        }
        
        df_list <- df_list[sapply(df_list, function(df) ncol(df) != 11)]

        for (name in names(df_list)) {
        df <- df_list[[name]]
        if (ncol(df) == 11) {
            cat("----", name, "----\n")
            print(dim(df))
            # Optionally print column names
            # print(colnames(df))
        }
        }

        # Merge all dataframes into one
        final_SALS_BA4 <- do.call(rbind, df_list)

        ## Calculate sample wise mean
        final_SALS_BA4 <- final_SALS_BA4 %>%
            group_by(donor, celltype, group, n_genes, n_epochs) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy),        
                .groups = "drop"               
            )

        ## Calculate the mean and standard deviation
        summary_stats <- final_SALS_BA4 %>%
            group_by(celltype, group, n_genes, n_epochs) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## set factor levels
        #summary_stats$group <- factor(summary_stats$group, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD"))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC' ))


        ## Plot
        ggplot(summary_stats, aes(x = as.factor(n_epochs), y = median_accuracy, fill = as.factor(n_epochs), colour = as.factor(n_epochs))) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(n_genes~celltype, scales = "free_x") +
        ylab("Accuracy") 
        #coord_cartesian(ylim=c(0.5,1.0)) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/bar_epoch_sweep_SALS_BA9.pdf', height = 4, width = 13)

        ## Save the complete file
        write.csv(summary_stats, paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_summary_CombatSeq_',status,'_',brain_region,'_narval_2.csv'))

        ## Define the optimal topic-epoch combination by considering both mean and median accuracy
        summary_stats_optimal <- summary_stats %>%
            mutate(composite_score = (median_accuracy + mean_accuracy) / 2) %>%
            group_by(celltype) %>%
            arrange(desc(composite_score)) %>%
            slice(1) %>%
            ungroup()

        table(summary_stats_optimal$celltype)

        summary_stats_optimal <- summary_stats_optimal %>%
            group_by(celltype) %>%
            slice_min(sd_accuracy, n = 1, with_ties = FALSE) %>%
            ungroup()
        
        table(summary_stats_optimal$celltype)

        write.csv(summary_stats_optimal, paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_optimal_CombatSeq_',status,'_',brain_region,'_narval_2.csv'))
    ##

    ###########################
    ## C9ALS BA4
    ###########################
    ## code
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'C9ALS'
        brain_region <- 'BA4'
        topic_list <- c(100, 200, 300)
        epoch_list <- c(35, 50, 100, 200, 300, 400, 500)


        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
            print(i)
            for (topic in topic_list) {
                for (epoch in epoch_list) {
                    # Construct the file path for the current cell type
                    file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs_cedar/model_outs/regETM_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_',topic,'_',epoch,'_narval_2.csv')

                    # Check if the file exists
                    if (file.exists(file_path)) {
                        # Read the CSV file
                        current_df <- read.csv(file_path)
                        df_list[[paste0(i, '_',topic,'_',epoch)]] <- current_df
                    } else {
                        warning(paste("File does not exist:", file_path))
                    }
                }
            }
        }
        

        # Merge all dataframes into one
        final_SALS_BA4 <- do.call(rbind, df_list)

        ## Calculate sample wise mean
        final_SALS_BA4 <- final_SALS_BA4 %>%
            group_by(donor, celltype, group, n_genes, n_epochs) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy),        
                .groups = "drop"               
            )

        ## Calculate the mean and standard deviation
        summary_stats <- final_SALS_BA4 %>%
            group_by(celltype, group, n_genes, n_epochs) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## set factor levels
        #summary_stats$group <- factor(summary_stats$group, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD"))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC' ))


        ## Plot
        ggplot(summary_stats, aes(x = as.factor(n_epochs), y = median_accuracy, fill = as.factor(n_epochs), colour = as.factor(n_epochs))) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(n_genes~celltype, scales = "free_x") +
        ylab("Accuracy") 
        #coord_cartesian(ylim=c(0.5,1.0)) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/bar_epoch_sweep_C9ALS_BA4.pdf', height = 4, width = 13)

        ## Save the complete file
        write.csv(summary_stats, paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_summary_CombatSeq_',status,'_',brain_region,'_narval_2.csv'))

        ## Define the optimal topic-epoch combination by considering both mean and median accuracy
        summary_stats_optimal <- summary_stats %>%
            mutate(composite_score = (median_accuracy + mean_accuracy) / 2) %>%
            group_by(celltype) %>%
            arrange(desc(composite_score)) %>%
            slice(1) %>%
            ungroup()

        table(summary_stats_optimal$celltype)

        summary_stats_optimal <- summary_stats_optimal %>%
            group_by(celltype) %>%
            slice_min(sd_accuracy, n = 1, with_ties = FALSE) %>%
            ungroup()
        
        table(summary_stats_optimal$celltype)

        write.csv(summary_stats_optimal, paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_optimal_CombatSeq_',status,'_',brain_region,'_narval_2.csv'))
    ##

    ###########################
    ## C9ALS BA9
    ###########################
    ## code
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'C9ALS'
        brain_region <- 'BA9'
        topic_list <- c(100, 200, 300)
        epoch_list <- c(35, 50, 100, 200, 300, 400, 500)


        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
            print(i)
            for (topic in topic_list) {
                for (epoch in epoch_list) {
                    # Construct the file path for the current cell type
                    file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs_graham/model_outs/regETM_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_',topic,'_',epoch,'_narval_2.csv')

                    # Check if the file exists
                    if (file.exists(file_path)) {
                        # Read the CSV file
                        current_df <- read.csv(file_path)
                        df_list[[paste0(i, '_',topic,'_',epoch)]] <- current_df
                    } else {
                        warning(paste("File does not exist:", file_path))
                    }
                }
            }
        }
        

        # Merge all dataframes into one
        final_SALS_BA4 <- do.call(rbind, df_list)

        ## Calculate sample wise mean
        final_SALS_BA4 <- final_SALS_BA4 %>%
            group_by(donor, celltype, group, n_genes, n_epochs) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy),        
                .groups = "drop"               
            )

        ## Calculate the mean and standard deviation
        summary_stats <- final_SALS_BA4 %>%
            group_by(celltype, group, n_genes, n_epochs) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## set factor levels
        #summary_stats$group <- factor(summary_stats$group, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD"))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC' ))


        ## Plot
        ggplot(summary_stats, aes(x = as.factor(n_epochs), y = median_accuracy, fill = as.factor(n_epochs), colour = as.factor(n_epochs))) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(n_genes~celltype, scales = "free_x") +
        ylab("Accuracy") 
        #coord_cartesian(ylim=c(0.5,1.0)) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/bar_epoch_sweep_C9ALS_BA9.pdf', height = 4, width = 13)

        ## Save the complete file
        write.csv(summary_stats, paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_summary_CombatSeq_',status,'_',brain_region,'_narval_2.csv'))

        ## Define the optimal topic-epoch combination by considering both mean and median accuracy
        summary_stats_optimal <- summary_stats %>%
            mutate(composite_score = (median_accuracy + mean_accuracy) / 2) %>%
            group_by(celltype) %>%
            arrange(desc(composite_score)) %>%
            slice(1) %>%
            ungroup()

        table(summary_stats_optimal$celltype)

        summary_stats_optimal <- summary_stats_optimal %>%
            group_by(celltype) %>%
            slice_min(sd_accuracy, n = 1, with_ties = FALSE) %>%
            ungroup()
        
        table(summary_stats_optimal$celltype)

        write.csv(summary_stats_optimal, paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_optimal_CombatSeq_',status,'_',brain_region,'_narval_2.csv'))
    ##

    ###########################
    ## All ALS BA4
    ###########################
    ## code
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'All ALS'
        brain_region <- 'BA4'
        topic_list <- c(100, 200, 300)
        epoch_list <- c(35, 50, 100, 200, 300, 400, 500)


        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
            print(i)
            for (topic in topic_list) {
                for (epoch in epoch_list) {
                    # Construct the file path for the current cell type
                    file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_',topic,'_',epoch,'_narval_2.csv')

                    # Check if the file exists
                    if (file.exists(file_path)) {
                        # Read the CSV file
                        current_df <- read.csv(file_path)
                        df_list[[paste0(i, '_',topic,'_',epoch)]] <- current_df
                    } else {
                        warning(paste("File does not exist:", file_path))
                    }
                }
            }
        }
        

        # Merge all dataframes into one
        final_SALS_BA4 <- do.call(rbind, df_list)

        ## Calculate sample wise mean
        final_SALS_BA4 <- final_SALS_BA4 %>%
            group_by(donor, celltype, group, n_genes, n_epochs) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy),        
                .groups = "drop"               
            )

        ## Calculate the mean and standard deviation
        summary_stats <- final_SALS_BA4 %>%
            group_by(celltype, group, n_genes, n_epochs) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## set factor levels
        #summary_stats$group <- factor(summary_stats$group, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD"))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC' ))


        ## Plot
        ggplot(summary_stats, aes(x = as.factor(n_epochs), y = median_accuracy, fill = as.factor(n_epochs), colour = as.factor(n_epochs))) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(n_genes~celltype, scales = "free_x") +
        ylab("Accuracy") 
        #coord_cartesian(ylim=c(0.5,1.0)) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/bar_epoch_sweep_All_ALS_BA4.pdf', height = 4, width = 13)

        ## Save the complete file
        write.csv(summary_stats, paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_summary_CombatSeq_',status,'_',brain_region,'_narval_2.csv'))

        ## Define the optimal topic-epoch combination by considering both mean and median accuracy
        summary_stats_optimal <- summary_stats %>%
            mutate(composite_score = (median_accuracy + mean_accuracy) / 2) %>%
            group_by(celltype) %>%
            arrange(desc(composite_score)) %>%
            slice(1) %>%
            ungroup()

        table(summary_stats_optimal$celltype)

        summary_stats_optimal <- summary_stats_optimal %>%
            group_by(celltype) %>%
            slice_min(sd_accuracy, n = 1, with_ties = FALSE) %>%
            ungroup()
        
        table(summary_stats_optimal$celltype)

        write.csv(summary_stats_optimal, paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_optimal_CombatSeq_',status,'_',brain_region,'_narval_2.csv'))
    ##

    ###########################
    ## All ALS BA9
    ###########################
    ## code
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'All ALS'
        brain_region <- 'BA9'
        topic_list <- c(100, 200, 300)
        epoch_list <- c(35, 50, 100, 200, 300, 400, 500)


        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
            print(i)
            for (topic in topic_list) {
                for (epoch in epoch_list) {
                    # Construct the file path for the current cell type
                    file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs_beluga/model_outs/regETM_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_',topic,'_',epoch,'_narval_2.csv')

                    # Check if the file exists
                    if (file.exists(file_path)) {
                        # Read the CSV file
                        current_df <- read.csv(file_path)
                        df_list[[paste0(i, '_',topic,'_',epoch)]] <- current_df
                    } else {
                        warning(paste("File does not exist:", file_path))
                    }
                }
            }
        }
        

        # Merge all dataframes into one
        final_SALS_BA4 <- do.call(rbind, df_list)

        ## Calculate sample wise mean
        final_SALS_BA4 <- final_SALS_BA4 %>%
            group_by(donor, celltype, group, n_genes, n_epochs) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy),        
                .groups = "drop"               
            )

        ## Calculate the mean and standard deviation
        summary_stats <- final_SALS_BA4 %>%
            group_by(celltype, group, n_genes, n_epochs) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## set factor levels
        #summary_stats$group <- factor(summary_stats$group, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD"))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC' ))


        ## Plot
        ggplot(summary_stats, aes(x = as.factor(n_epochs), y = median_accuracy, fill = as.factor(n_epochs), colour = as.factor(n_epochs))) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(n_genes~celltype, scales = "free_x") +
        ylab("Accuracy") 
        #coord_cartesian(ylim=c(0.5,1.0)) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/bar_epoch_sweep_All_ALS_BA9.pdf', height = 4, width = 13)

        ## Save the complete file
        write.csv(summary_stats, paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_summary_CombatSeq_',status,'_',brain_region,'_narval_2.csv'))

        ## Define the optimal topic-epoch combination by considering both mean and median accuracy
        summary_stats_optimal <- summary_stats %>%
            mutate(composite_score = (median_accuracy + mean_accuracy) / 2) %>%
            group_by(celltype) %>%
            arrange(desc(composite_score)) %>%
            slice(1) %>%
            ungroup()

        table(summary_stats_optimal$celltype)

        summary_stats_optimal <- summary_stats_optimal %>%
            group_by(celltype) %>%
            slice_min(sd_accuracy, n = 1, with_ties = FALSE) %>%
            ungroup()
        
        table(summary_stats_optimal$celltype)

        write.csv(summary_stats_optimal, paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_optimal_CombatSeq_',status,'_',brain_region,'_narval_2.csv'))
    ##

    ###########################
    ## SFTLD BA4
    ###########################
    ## code
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'SFTLD'
        brain_region <- 'BA4'
        topic_list <- c(100, 200, 300)
        epoch_list <- c(35, 50, 100, 200, 300, 400, 500)


        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
            print(i)
            for (topic in topic_list) {
                for (epoch in epoch_list) {
                    # Construct the file path for the current cell type
                    file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_',topic,'_',epoch,'_narval_2.csv')

                    # Check if the file exists
                    if (file.exists(file_path)) {
                        # Read the CSV file
                        current_df <- read.csv(file_path)
                        df_list[[paste0(i, '_',topic,'_',epoch)]] <- current_df
                    } else {
                        warning(paste("File does not exist:", file_path))
                    }
                }
            }
        }
        

        # Merge all dataframes into one
        final_SALS_BA4 <- do.call(rbind, df_list)

        ## Calculate sample wise mean
        final_SALS_BA4 <- final_SALS_BA4 %>%
            group_by(donor, celltype, group, n_genes, n_epochs) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy),        
                .groups = "drop"               
            )

        ## Calculate the mean and standard deviation
        summary_stats <- final_SALS_BA4 %>%
            group_by(celltype, group, n_genes, n_epochs) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## set factor levels
        #summary_stats$group <- factor(summary_stats$group, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD"))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC' ))


        ## Plot
        ggplot(summary_stats, aes(x = as.factor(n_epochs), y = median_accuracy, fill = as.factor(n_epochs), colour = as.factor(n_epochs))) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(n_genes~celltype, scales = "free_x") +
        ylab("Accuracy") 
        #coord_cartesian(ylim=c(0.5,1.0)) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/bar_epoch_sweep_SFTLD_BA4.pdf', height = 4, width = 13)

        ## Save the complete file
        write.csv(summary_stats, paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_summary_CombatSeq_',status,'_',brain_region,'_narval_2.csv'))

        ## Define the optimal topic-epoch combination by considering both mean and median accuracy
        summary_stats_optimal <- summary_stats %>%
            mutate(composite_score = (median_accuracy + mean_accuracy) / 2) %>%
            group_by(celltype) %>%
            arrange(desc(composite_score)) %>%
            slice(1) %>%
            ungroup()

        table(summary_stats_optimal$celltype)

        summary_stats_optimal <- summary_stats_optimal %>%
            group_by(celltype) %>%
            slice_min(sd_accuracy, n = 1, with_ties = FALSE) %>%
            ungroup()
        
        table(summary_stats_optimal$celltype)

        write.csv(summary_stats_optimal, paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_optimal_CombatSeq_',status,'_',brain_region,'_narval_2.csv'))
    ##

    ###########################
    ## SFTLD BA9
    ###########################
    ## code
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'SFTLD'
        brain_region <- 'BA9'
        topic_list <- c(100, 200, 300)
        epoch_list <- c(35, 50, 100, 200, 300, 400, 500)


        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
            print(i)
            for (topic in topic_list) {
                for (epoch in epoch_list) {
                    # Construct the file path for the current cell type
                    file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_',topic,'_',epoch,'_narval_2.csv')

                    # Check if the file exists
                    if (file.exists(file_path)) {
                        # Read the CSV file
                        current_df <- read.csv(file_path)
                        df_list[[paste0(i, '_',topic,'_',epoch)]] <- current_df
                    } else {
                        warning(paste("File does not exist:", file_path))
                    }
                }
            }
        }
        

        # Merge all dataframes into one
        final_SALS_BA4 <- do.call(rbind, df_list)

        ## Calculate sample wise mean
        final_SALS_BA4 <- final_SALS_BA4 %>%
            group_by(donor, celltype, group, n_genes, n_epochs) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy),        
                .groups = "drop"               
            )

        ## Calculate the mean and standard deviation
        summary_stats <- final_SALS_BA4 %>%
            group_by(celltype, group, n_genes, n_epochs) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## set factor levels
        #summary_stats$group <- factor(summary_stats$group, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD"))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC' ))


        ## Plot
        ggplot(summary_stats, aes(x = as.factor(n_epochs), y = median_accuracy, fill = as.factor(n_epochs), colour = as.factor(n_epochs))) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(n_genes~celltype, scales = "free_x") +
        ylab("Accuracy") 
        #coord_cartesian(ylim=c(0.5,1.0)) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/bar_epoch_sweep_SFTLD_BA9.pdf', height = 4, width = 13)

        ## Save the complete file
        write.csv(summary_stats, paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_summary_CombatSeq_',status,'_',brain_region,'_narval_2.csv'))

        ## Define the optimal topic-epoch combination by considering both mean and median accuracy
        summary_stats_optimal <- summary_stats %>%
            mutate(composite_score = (median_accuracy + mean_accuracy) / 2) %>%
            group_by(celltype) %>%
            arrange(desc(composite_score)) %>%
            slice(1) %>%
            ungroup()

        table(summary_stats_optimal$celltype)

        summary_stats_optimal <- summary_stats_optimal %>%
            group_by(celltype) %>%
            slice_min(sd_accuracy, n = 1, with_ties = FALSE) %>%
            ungroup()
        
        table(summary_stats_optimal$celltype)

        write.csv(summary_stats_optimal, paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_optimal_CombatSeq_',status,'_',brain_region,'_narval_2.csv'))
    ##

    ###########################
    ## C9FTLD BA4
    ###########################
    ## code
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'C9FTLD'
        brain_region <- 'BA4'
        topic_list <- c(100, 200, 300)
        epoch_list <- c(35, 50, 100, 200, 300, 400, 500)


        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
            print(i)
            for (topic in topic_list) {
                for (epoch in epoch_list) {
                    # Construct the file path for the current cell type
                    file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_',topic,'_',epoch,'_narval_2.csv')

                    # Check if the file exists
                    if (file.exists(file_path)) {
                        # Read the CSV file
                        current_df <- read.csv(file_path)
                        df_list[[paste0(i, '_',topic,'_',epoch)]] <- current_df
                    } else {
                        warning(paste("File does not exist:", file_path))
                    }
                }
            }
        }
        
        # Merge all dataframes into one
        final_SALS_BA4 <- do.call(rbind, df_list)

        ## Calculate sample wise mean
        final_SALS_BA4 <- final_SALS_BA4 %>%
            group_by(donor, celltype, group, n_genes, n_epochs) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy),        
                .groups = "drop"               
            )

        ## Calculate the mean and standard deviation
        summary_stats <- final_SALS_BA4 %>%
            group_by(celltype, group, n_genes, n_epochs) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## set factor levels
        #summary_stats$group <- factor(summary_stats$group, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD"))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC' ))


        ## Plot
        ggplot(summary_stats, aes(x = as.factor(n_epochs), y = median_accuracy, fill = as.factor(n_epochs), colour = as.factor(n_epochs))) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(n_genes~celltype, scales = "free_x") +
        ylab("Accuracy") 
        #coord_cartesian(ylim=c(0.5,1.0)) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/bar_epoch_sweep_C9FTLD_BA4.pdf', height = 4, width = 13)

        ## Save the complete file
        write.csv(summary_stats, paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_summary_CombatSeq_',status,'_',brain_region,'_narval_2.csv'))

        ## Define the optimal topic-epoch combination by considering both mean and median accuracy
        summary_stats_optimal <- summary_stats %>%
            mutate(composite_score = (median_accuracy + mean_accuracy) / 2) %>%
            group_by(celltype) %>%
            arrange(desc(composite_score)) %>%
            slice(1) %>%
            ungroup()

        table(summary_stats_optimal$celltype)

        summary_stats_optimal <- summary_stats_optimal %>%
            group_by(celltype) %>%
            slice_min(sd_accuracy, n = 1, with_ties = FALSE) %>%
            ungroup()
        
        table(summary_stats_optimal$celltype)

        write.csv(summary_stats_optimal, paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_optimal_CombatSeq_',status,'_',brain_region,'_narval_2.csv'))
    ##

    ###########################
    ## C9FTLD BA9
    ###########################
    ## code
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'C9FTLD'
        brain_region <- 'BA9'
        topic_list <- c(100, 200, 300)
        epoch_list <- c(35, 50, 100, 200, 300, 400, 500)


        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
            print(i)
            for (topic in topic_list) {
                for (epoch in epoch_list) {
                    # Construct the file path for the current cell type
                    file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_',topic,'_',epoch,'_narval_2.csv')

                    # Check if the file exists
                    if (file.exists(file_path)) {
                        # Read the CSV file
                        current_df <- read.csv(file_path)
                        df_list[[paste0(i, '_',topic,'_',epoch)]] <- current_df
                    } else {
                        warning(paste("File does not exist:", file_path))
                    }
                }
            }
        }
        

        for (name in names(df_list)) {
        df <- df_list[[name]]
        if (ncol(df) == 11) {
            cat("----", name, "----\n")
            print(dim(df))
            # Optionally print column names
            # print(colnames(df))
        }
        }
        
        df_list <- df_list[sapply(df_list, function(df) ncol(df) != 11)]

        for (name in names(df_list)) {
        df <- df_list[[name]]
        if (ncol(df) == 11) {
            cat("----", name, "----\n")
            print(dim(df))
            # Optionally print column names
            # print(colnames(df))
        }
        }

        # Merge all dataframes into one
        final_SALS_BA4 <- do.call(rbind, df_list)

        ## Calculate sample wise mean
        final_SALS_BA4 <- final_SALS_BA4 %>%
            group_by(donor, celltype, group, n_genes, n_epochs) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy),        
                .groups = "drop"               
            )

        ## Calculate the mean and standard deviation
        summary_stats <- final_SALS_BA4 %>%
            group_by(celltype, group, n_genes, n_epochs) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## set factor levels
        #summary_stats$group <- factor(summary_stats$group, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD"))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC' ))


        ## Plot
        ggplot(summary_stats, aes(x = as.factor(n_epochs), y = median_accuracy, fill = as.factor(n_epochs), colour = as.factor(n_epochs))) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(n_genes~celltype, scales = "free_x") +
        ylab("Accuracy") 
        #coord_cartesian(ylim=c(0.5,1.0)) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/bar_epoch_sweep_C9FTLD_BA9.pdf', height = 4, width = 13)

        ## Save the complete file
        write.csv(summary_stats, paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_summary_CombatSeq_',status,'_',brain_region,'_narval_2.csv'))

        ## Define the optimal topic-epoch combination by considering both mean and median accuracy
        summary_stats_optimal <- summary_stats %>%
            mutate(composite_score = (median_accuracy + mean_accuracy) / 2) %>%
            group_by(celltype) %>%
            arrange(desc(composite_score)) %>%
            slice(1) %>%
            ungroup()

        table(summary_stats_optimal$celltype)

        summary_stats_optimal <- summary_stats_optimal %>%
            group_by(celltype) %>%
            slice_min(sd_accuracy, n = 1, with_ties = FALSE) %>%
            ungroup()
        
        table(summary_stats_optimal$celltype)

        write.csv(summary_stats_optimal, paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_optimal_CombatSeq_',status,'_',brain_region,'_narval_2.csv'))
    ##


    ###########################
    ## All FTLD BA4
    ###########################
    ## code
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'All FTLD'
        brain_region <- 'BA4'
        topic_list <- c(100, 200, 300)
        epoch_list <- c(35, 50, 100, 200, 300, 400, 500)


        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
            print(i)
            for (topic in topic_list) {
                for (epoch in epoch_list) {
                    # Construct the file path for the current cell type
                    file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs_beluga/model_outs/regETM_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_',topic,'_',epoch,'_narval_2.csv')

                    # Check if the file exists
                    if (file.exists(file_path)) {
                        # Read the CSV file
                        current_df <- read.csv(file_path)
                        df_list[[paste0(i, '_',topic,'_',epoch)]] <- current_df
                    } else {
                        warning(paste("File does not exist:", file_path))
                    }
                }
            }
        }
        

        # Merge all dataframes into one
        final_SALS_BA4 <- do.call(rbind, df_list)

        ## Calculate sample wise mean
        final_SALS_BA4 <- final_SALS_BA4 %>%
            group_by(donor, celltype, group, n_genes, n_epochs) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy),        
                .groups = "drop"               
            )

        ## Calculate the mean and standard deviation
        summary_stats <- final_SALS_BA4 %>%
            group_by(celltype, group, n_genes, n_epochs) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## set factor levels
        #summary_stats$group <- factor(summary_stats$group, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD"))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC' ))


        ## Plot
        ggplot(summary_stats, aes(x = as.factor(n_epochs), y = median_accuracy, fill = as.factor(n_epochs), colour = as.factor(n_epochs))) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(n_genes~celltype, scales = "free_x") +
        ylab("Accuracy") 
        #coord_cartesian(ylim=c(0.5,1.0)) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/bar_epoch_sweep_All_FTLD_BA4.pdf', height = 4, width = 13)

        ## Save the complete file
        write.csv(summary_stats, paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_summary_CombatSeq_',status,'_',brain_region,'_narval_2.csv'))

        ## Define the optimal topic-epoch combination by considering both mean and median accuracy
        summary_stats_optimal <- summary_stats %>%
            mutate(composite_score = (median_accuracy + mean_accuracy) / 2) %>%
            group_by(celltype) %>%
            arrange(desc(composite_score)) %>%
            slice(1) %>%
            ungroup()

        table(summary_stats_optimal$celltype)

        summary_stats_optimal <- summary_stats_optimal %>%
            group_by(celltype) %>%
            slice_min(sd_accuracy, n = 1, with_ties = FALSE) %>%
            ungroup()
        
        table(summary_stats_optimal$celltype)

        write.csv(summary_stats_optimal, paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_optimal_CombatSeq_',status,'_',brain_region,'_narval_2.csv'))
    ##


    ###########################
    ## All FTLD BA9
    ###########################
    ## code
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'All FTLD'
        brain_region <- 'BA9'
        topic_list <- c(100, 200, 300)
        epoch_list <- c(35, 50, 100, 200, 300, 400, 500)


        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
            print(i)
            for (topic in topic_list) {
                for (epoch in epoch_list) {
                    # Construct the file path for the current cell type
                    file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs_beluga/model_outs/regETM_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_',topic,'_',epoch,'_narval_2.csv')

                    # Check if the file exists
                    if (file.exists(file_path)) {
                        # Read the CSV file
                        current_df <- read.csv(file_path)
                        df_list[[paste0(i, '_',topic,'_',epoch)]] <- current_df
                    } else {
                        warning(paste("File does not exist:", file_path))
                    }
                }
            }
        }
        

        # Merge all dataframes into one
        final_SALS_BA4 <- do.call(rbind, df_list)

        ## Calculate sample wise mean
        final_SALS_BA4 <- final_SALS_BA4 %>%
            group_by(donor, celltype, group, n_genes, n_epochs) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy),        
                .groups = "drop"               
            )

        ## Calculate the mean and standard deviation
        summary_stats <- final_SALS_BA4 %>%
            group_by(celltype, group, n_genes, n_epochs) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## set factor levels
        #summary_stats$group <- factor(summary_stats$group, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD"))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC' ))


        ## Plot
        ggplot(summary_stats, aes(x = as.factor(n_epochs), y = median_accuracy, fill = as.factor(n_epochs), colour = as.factor(n_epochs))) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(n_genes~celltype, scales = "free_x") +
        ylab("Accuracy") 
        #coord_cartesian(ylim=c(0.5,1.0)) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/bar_epoch_sweep_All_FTLD_BA9.pdf', height = 4, width = 13)

        ## Save the complete file
        write.csv(summary_stats, paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_summary_CombatSeq_',status,'_',brain_region,'_narval_2.csv'))

        ## Define the optimal topic-epoch combination by considering both mean and median accuracy
        summary_stats_optimal <- summary_stats %>%
            mutate(composite_score = (median_accuracy + mean_accuracy) / 2) %>%
            group_by(celltype) %>%
            arrange(desc(composite_score)) %>%
            slice(1) %>%
            ungroup()

        table(summary_stats_optimal$celltype)

        summary_stats_optimal <- summary_stats_optimal %>%
            group_by(celltype) %>%
            slice_min(sd_accuracy, n = 1, with_ties = FALSE) %>%
            ungroup()
        
        table(summary_stats_optimal$celltype)

        write.csv(summary_stats_optimal, paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_optimal_CombatSeq_',status,'_',brain_region,'_narval_2.csv'))
    ##

##

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################ Optimal Topic-epoch: reg - no bayesian
## code
    ###########################
    ## BA4 and BA9
    ###########################
    ## code
        ## load in files
        SALS_BA4 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_optimal_CombatSeq_SALS_BA4_narval_2.csv')
        C9ALS_BA4 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_optimal_CombatSeq_C9ALS_BA4_narval_2.csv')
        SFTLD_BA4 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_optimal_CombatSeq_SFTLD_BA4_narval_2.csv')
        C9FTLD_BA4 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_optimal_CombatSeq_C9FTLD_BA4_narval_2.csv')
        All_ALS_BA4 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_optimal_CombatSeq_All ALS_BA4_narval_2.csv')
        All_FTLD_BA4 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_optimal_CombatSeq_All FTLD_BA4_narval_2.csv')

        ## bind all
        BA4_bind_all <- rbind(SALS_BA4, C9ALS_BA4, SFTLD_BA4, C9FTLD_BA4, All_ALS_BA4, All_FTLD_BA4 )

        ## add brain region
        BA4_bind_all$brain_region <- "BA4"

        ## load in files
        SALS_BA9 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_optimal_CombatSeq_SALS_BA9_narval_2.csv')
        C9ALS_BA9 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_optimal_CombatSeq_C9ALS_BA9_narval_2.csv')
        SFTLD_BA9 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_optimal_CombatSeq_SFTLD_BA9_narval_2.csv')
        C9FTLD_BA9 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_optimal_CombatSeq_C9FTLD_BA9_narval_2.csv')
        All_ALS_BA9 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_optimal_CombatSeq_All ALS_BA9_narval_2.csv')
        All_FTLD_BA9 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_optimal_CombatSeq_All FTLD_BA9_narval_2.csv')

        ## bind all
        BA9_bind_all <- rbind(SALS_BA9, C9ALS_BA9, SFTLD_BA9, C9FTLD_BA9, All_ALS_BA9, All_FTLD_BA9)

        ## add brain region
        BA9_bind_all$brain_region <- "BA9"

        ## Bind all
        bind_all <- rbind(BA4_bind_all, BA9_bind_all)

        ## set factor levels
        bind_all$group <- factor(bind_all$group, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All ALS", "All FTLD"))
        bind_all$celltype <- factor(bind_all$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))
        
        ## Plot
        ggplot(bind_all, aes(x = group, y = median_accuracy, fill = group, colour = group)) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(brain_region~celltype, scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) #+
        #coord_cartesian(ylim=c(0.5,1.0)) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/Bar_optimal_ETM_no_bayes.pdf', height = 3, width = 13)
    
        mean(bind_all$median_accuracy)

    ##

    ###########################
    ## Compute mean across all models
    ###########################
    ## code 
        ## load in files
        SALS_BA4 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_optimal_CombatSeq_SALS_BA4_narval_2.csv')
        C9ALS_BA4 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_optimal_CombatSeq_C9ALS_BA4_narval_2.csv')
        SFTLD_BA4 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_optimal_CombatSeq_SFTLD_BA4_narval_2.csv')
        C9FTLD_BA4 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_optimal_CombatSeq_C9FTLD_BA4_narval_2.csv')
        All_ALS_BA4 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_optimal_CombatSeq_All ALS_BA4_narval_2.csv')
        All_FTLD_BA4 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_optimal_CombatSeq_All FTLD_BA4_narval_2.csv')

        ## bind all
        BA4_bind_all <- rbind(SALS_BA4, C9ALS_BA4, SFTLD_BA4, C9FTLD_BA4, All_ALS_BA4, All_FTLD_BA4 )

        ## add brain region
        BA4_bind_all$brain_region <- "BA4"
        
        ## load in files
        SALS_BA9 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_optimal_CombatSeq_SALS_BA9_narval_2.csv')
        C9ALS_BA9 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_optimal_CombatSeq_C9ALS_BA9_narval_2.csv')
        SFTLD_BA9 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_optimal_CombatSeq_SFTLD_BA9_narval_2.csv')
        C9FTLD_BA9 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_optimal_CombatSeq_C9FTLD_BA9_narval_2.csv')
        All_ALS_BA9 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_optimal_CombatSeq_All ALS_BA9_narval_2.csv')
        All_FTLD_BA9 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_optimal_CombatSeq_All FTLD_BA9_narval_2.csv')

        ## bind all
        BA9_bind_all <- rbind(SALS_BA9, C9ALS_BA9, SFTLD_BA9, C9FTLD_BA9, All_ALS_BA9, All_FTLD_BA9 )

        ## add brain region
        BA9_bind_all$brain_region <- "BA9"

        ## bind all and compute mean
        bind_all <- rbind(BA4_bind_all, BA9_bind_all)

        ## compute overall mean
        mean(bind_all$median_accuracy)

        bind_all$group <- factor(bind_all$group, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All ALS", "All FTLD"))
        bind_all$celltype <- factor(bind_all$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))

        ## Calculate mean of groups across brain regions
        summary_stats_all <- bind_all %>%
            group_by(group, brain_region) %>%
            summarise(
                mean_accuracy = mean(median_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(median_accuracy),                  # Calculate the standard deviation of X0
                .groups = "drop"                 # Drop the grouping after summarising
            )
        
        ## Plot
        ggplot(summary_stats_all, aes(x = group, y = mean_accuracy, fill = group, colour = group, label = round(mean_accuracy, 2) )) + 
        geom_errorbar(aes(ymin=mean_accuracy - sd_accuracy, ymax=mean_accuracy + sd_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_text(colour = "white", size = 3, vjust = 1) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(brain_region ~ ., scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 4, width = 3)

        ## By major celltype group
        bind_all <- data.frame(bind_all)
        bind_all$major_group[bind_all$celltype == "L2_L3"] <- "Excitatory"
        bind_all$major_group[bind_all$celltype == "L3_L5"] <- "Excitatory"
        bind_all$major_group[bind_all$celltype == "L4_L5"] <- "Excitatory"
        bind_all$major_group[bind_all$celltype == "L4_L6"] <- "Excitatory"
        bind_all$major_group[bind_all$celltype == "L5"] <- "Excitatory"
        bind_all$major_group[bind_all$celltype == "L5_L6"] <- "Excitatory"
        bind_all$major_group[bind_all$celltype == "L6"] <- "Excitatory"

        bind_all$major_group[bind_all$celltype == "5HT3aR"] <- "Inhibitory"
        bind_all$major_group[bind_all$celltype == "PV"] <- "Inhibitory"
        bind_all$major_group[bind_all$celltype == "Rosehip"] <- "Inhibitory"
        bind_all$major_group[bind_all$celltype == "SOM"] <- "Inhibitory"

        bind_all$major_group[bind_all$celltype == "Astro"] <- "non-neuronal"
        bind_all$major_group[bind_all$celltype == "Endo"] <- "non-neuronal"
        bind_all$major_group[bind_all$celltype == "Fibro"] <- "non-neuronal"
        bind_all$major_group[bind_all$celltype == "Micro"] <- "non-neuronal"
        bind_all$major_group[bind_all$celltype == "Mural"] <- "non-neuronal"
        bind_all$major_group[bind_all$celltype == "Oligo"] <- "non-neuronal"
        bind_all$major_group[bind_all$celltype == "OPC"] <- "non-neuronal"

        ## Calculate mean of groups across brain regions
        summary_stats_all <- bind_all %>%
            group_by(group, brain_region, major_group) %>%
            summarise(
                mean_accuracy = mean(median_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(median_accuracy),                  # Calculate the standard deviation of X0
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## Plot
        ggplot(summary_stats_all, aes(x = group, y = mean_accuracy, fill = group, colour = group, label = round(mean_accuracy, 2) )) + 
        geom_errorbar(aes(ymin=mean_accuracy - sd_accuracy, ymax=mean_accuracy + sd_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_text(colour = "white", size = 3, vjust = 1) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(brain_region ~ major_group, scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 4, width = 6)

        ## Calculate mean of major cell groups
        summary_stats_all <- bind_all %>%
            group_by(major_group) %>%
            summarise(
                mean_accuracy = mean(median_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(median_accuracy),                  # Calculate the standard deviation of X0
                .groups = "drop"                 # Drop the grouping after summarising
            )
##


################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################ Optimal Topic-epoch: with bayesian 
## code 

    ## BA4
    ###########################
    ## SALS BA4 -- done
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'SALS'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
            # Construct the file path for the current cell type
            info_df <- read.csv(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_optimal_CombatSeq_',status,'_',brain_region,'_narval_2.csv'))
            info_df <- subset(info_df, celltype == i)
            par_n_topics <- info_df$n_genes
            par_n_epochs <- info_df$n_epochs
            file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_',par_n_topics,'_',par_n_epochs,'_narval_2_KNN_Bayesian.csv')
            # Check if the file exists
            if (file.exists(file_path)) {
                # Read the CSV file
                current_df <- read.csv(file_path)
                
                # Optionally, add the cell type as a new column to identify the source of each row
                current_df$dataset <- paste0(status, "_", brain_region)
                
                # Add the dataframe to the list
                df_list[[i]] <- current_df
            } else {
                warning(paste("File does not exist:", file_path))
            }
        }

        # Merge all dataframes into one
        final_SALS_BA4 <- do.call(rbind, df_list)
    ##    

    ###########################
    ## C9ALS BA4 -- not done
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'C9ALS'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
            # Construct the file path for the current cell type
            info_df <- read.csv(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_optimal_CombatSeq_',status,'_',brain_region,'_narval_2.csv'))
            info_df <- subset(info_df, celltype == i)
            par_n_topics <- info_df$n_genes
            par_n_epochs <- info_df$n_epochs
            file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_',par_n_topics,'_',par_n_epochs,'_narval_2_KNN_Bayesian.csv')
            # Check if the file exists
            if (file.exists(file_path)) {
                # Read the CSV file
                current_df <- read.csv(file_path)
                
                # Optionally, add the cell type as a new column to identify the source of each row
                current_df$dataset <- paste0(status, "_", brain_region)
                
                # Add the dataframe to the list
                df_list[[i]] <- current_df
            } else {
                warning(paste("File does not exist:", file_path))
            }
        }

        # Merge all dataframes into one
        final_C9ALS_BA4 <- do.call(rbind, df_list)
    ##  

    ###########################
    ## SFTLD BA4 -- done
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'SFTLD'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
            # Construct the file path for the current cell type
            info_df <- read.csv(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_optimal_CombatSeq_',status,'_',brain_region,'_narval_2.csv'))
            info_df <- subset(info_df, celltype == i)
            par_n_topics <- info_df$n_genes
            par_n_epochs <- info_df$n_epochs
            file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_',par_n_topics,'_',par_n_epochs,'_narval_2_KNN_Bayesian.csv')
            # Check if the file exists
            if (file.exists(file_path)) {
                # Read the CSV file
                current_df <- read.csv(file_path)
                
                # Optionally, add the cell type as a new column to identify the source of each row
                current_df$dataset <- paste0(status, "_", brain_region)
                
                # Add the dataframe to the list
                df_list[[i]] <- current_df
            } else {
                warning(paste("File does not exist:", file_path))
            }
        }

        # Merge all dataframes into one
        final_SFTLD_BA4 <- do.call(rbind, df_list)
    ##  

    ###########################
    ## C9FTLD BA4
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'C9FTLD'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
            # Construct the file path for the current cell type
            info_df <- read.csv(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_optimal_CombatSeq_',status,'_',brain_region,'_narval_2.csv'))
            info_df <- subset(info_df, celltype == i)
            par_n_topics <- info_df$n_genes
            par_n_epochs <- info_df$n_epochs
            file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_',par_n_topics,'_',par_n_epochs,'_narval_2_KNN_Bayesian.csv')
            # Check if the file exists
            if (file.exists(file_path)) {
                # Read the CSV file
                current_df <- read.csv(file_path)
                
                # Optionally, add the cell type as a new column to identify the source of each row
                current_df$dataset <- paste0(status, "_", brain_region)
                
                # Add the dataframe to the list
                df_list[[i]] <- current_df
            } else {
                warning(paste("File does not exist:", file_path))
            }
        }

        # Merge all dataframes into one
        final_C9FTLD_BA4 <- do.call(rbind, df_list)
    ##  

    ###########################
    ## All ALS BA4 
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'All ALS'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
            # Construct the file path for the current cell type
            info_df <- read.csv(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_optimal_CombatSeq_',status,'_',brain_region,'_narval_2.csv'))
            info_df <- subset(info_df, celltype == i)
            par_n_topics <- info_df$n_genes
            par_n_epochs <- info_df$n_epochs
            file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_',par_n_topics,'_',par_n_epochs,'_narval_2_KNN_Bayesian.csv')
            # Check if the file exists
            if (file.exists(file_path)) {
                # Read the CSV file
                current_df <- read.csv(file_path)
                
                # Optionally, add the cell type as a new column to identify the source of each row
                current_df$dataset <- paste0(status, "_", brain_region)
                
                # Add the dataframe to the list
                df_list[[i]] <- current_df
            } else {
                warning(paste("File does not exist:", file_path))
            }
        }

        # Merge all dataframes into one
        final_All_ALS_BA4 <- do.call(rbind, df_list)
    ##  

    ###########################
    ## All FTLD BA4 
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'All FTLD'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
            # Construct the file path for the current cell type
            info_df <- read.csv(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_optimal_CombatSeq_',status,'_',brain_region,'_narval_2.csv'))
            info_df <- subset(info_df, celltype == i)
            par_n_topics <- info_df$n_genes
            par_n_epochs <- info_df$n_epochs
            file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_',par_n_topics,'_',par_n_epochs,'_narval_2_KNN_Bayesian.csv')
            # Check if the file exists
            if (file.exists(file_path)) {
                # Read the CSV file
                current_df <- read.csv(file_path)
                
                # Optionally, add the cell type as a new column to identify the source of each row
                current_df$dataset <- paste0(status, "_", brain_region)
                
                # Add the dataframe to the list
                df_list[[i]] <- current_df
            } else {
                warning(paste("File does not exist:", file_path))
            }
        }

        # Merge all dataframes into one
        final_All_FTLD_BA4 <- do.call(rbind, df_list)
    ##  

    ###########################
    ## merge and plot BA4
    ###########################
    ## code for barplot
        BA4_bind_all <- rbind(final_SALS_BA4, final_C9ALS_BA4, final_SFTLD_BA4, final_C9FTLD_BA4, final_All_ALS_BA4, final_All_FTLD_BA4)
        write.csv(BA4_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_LOO_optimal_ABS_scaling_KNN_Bayesian_ETM.csv')
        
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_LOO_optimal_ABS_scaling_KNN_Bayesian_ETM.csv')
        
        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, group, kNN_thresh) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )

        ## Calculate the mean and standard deviation
        summary_stats <- BA4_bind_all %>%
            group_by(celltype, group, kNN_thresh) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## set factor levels
        summary_stats$group <- factor(summary_stats$group, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All ALS", "All FTLD"))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))

        ## Plot
        ggplot(summary_stats, aes(x = group, y = median_accuracy, fill = group, colour = group)) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(kNN_thresh~celltype, scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/ETM_BA4_LOO_balanced_acc_barplot_optimal_abs_scaling_KNN_bayesian.pdf', height = 8, width = 13)

        ## Plot only Bayes 0.9 thresh
        summary_stats_lim <- subset(summary_stats, kNN_thresh == 0.9 )
        ggplot(summary_stats_lim, aes(x = group, y = median_accuracy, fill = group, colour = group)) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~celltype, scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA4_LOO_balanced_acc_barplot_optimal_abs_scaling_KNN_bayesian_09.pdf', height = 2, width = 13)

        summary_stats_lim_BA4 <- summary_stats_lim

    ##

    ## BA9
    ###########################
    ## SALS BA9
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'SALS'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
            # Construct the file path for the current cell type
            info_df <- read.csv(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_optimal_CombatSeq_',status,'_',brain_region,'_narval_2.csv'))
            info_df <- subset(info_df, celltype == i)
            par_n_topics <- info_df$n_genes
            par_n_epochs <- info_df$n_epochs
            file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_',par_n_topics,'_',par_n_epochs,'_narval_2_KNN_Bayesian.csv')
            # Check if the file exists
            if (file.exists(file_path)) {
                # Read the CSV file
                current_df <- read.csv(file_path)
                
                # Optionally, add the cell type as a new column to identify the source of each row
                current_df$dataset <- paste0(status, "_", brain_region)
                
                # Add the dataframe to the list
                df_list[[i]] <- current_df
            } else {
                warning(paste("File does not exist:", file_path))
            }
        }

        # Merge all dataframes into one
        final_SALS_BA9 <- do.call(rbind, df_list)
    ##    

    ###########################
    ## C9ALS BA9
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'C9ALS'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
            # Construct the file path for the current cell type
            info_df <- read.csv(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_optimal_CombatSeq_',status,'_',brain_region,'_narval_2.csv'))
            info_df <- subset(info_df, celltype == i)
            par_n_topics <- info_df$n_genes
            par_n_epochs <- info_df$n_epochs
            file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_',par_n_topics,'_',par_n_epochs,'_narval_2_KNN_Bayesian.csv')
            # Check if the file exists
            if (file.exists(file_path)) {
                # Read the CSV file
                current_df <- read.csv(file_path)
                
                # Optionally, add the cell type as a new column to identify the source of each row
                current_df$dataset <- paste0(status, "_", brain_region)
                
                # Add the dataframe to the list
                df_list[[i]] <- current_df
            } else {
                warning(paste("File does not exist:", file_path))
            }
        }

        # Merge all dataframes into one
        final_C9ALS_BA9 <- do.call(rbind, df_list)
    ##  

    ###########################
    ## SFTLD BA9
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'SFTLD'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
            # Construct the file path for the current cell type
            info_df <- read.csv(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_optimal_CombatSeq_',status,'_',brain_region,'_narval_2.csv'))
            info_df <- subset(info_df, celltype == i)
            par_n_topics <- info_df$n_genes
            par_n_epochs <- info_df$n_epochs
            file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_',par_n_topics,'_',par_n_epochs,'_narval_2_KNN_Bayesian.csv')
            # Check if the file exists
            if (file.exists(file_path)) {
                # Read the CSV file
                current_df <- read.csv(file_path)
                
                # Optionally, add the cell type as a new column to identify the source of each row
                current_df$dataset <- paste0(status, "_", brain_region)
                
                # Add the dataframe to the list
                df_list[[i]] <- current_df
            } else {
                warning(paste("File does not exist:", file_path))
            }
        }

        # Merge all dataframes into one
        final_SFTLD_BA9 <- do.call(rbind, df_list)
    ##  

    ###########################
    ## C9FTLD BA9
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'C9FTLD'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
            # Construct the file path for the current cell type
            info_df <- read.csv(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_optimal_CombatSeq_',status,'_',brain_region,'_narval_2.csv'))
            info_df <- subset(info_df, celltype == i)
            par_n_topics <- info_df$n_genes
            par_n_epochs <- info_df$n_epochs
            file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_',par_n_topics,'_',par_n_epochs,'_narval_2_KNN_Bayesian.csv')
            # Check if the file exists
            if (file.exists(file_path)) {
                # Read the CSV file
                current_df <- read.csv(file_path)
                
                # Optionally, add the cell type as a new column to identify the source of each row
                current_df$dataset <- paste0(status, "_", brain_region)
                
                # Add the dataframe to the list
                df_list[[i]] <- current_df
            } else {
                warning(paste("File does not exist:", file_path))
            }
        }

        # Merge all dataframes into one
        final_C9FTLD_BA9 <- do.call(rbind, df_list)
    ##  

    ###########################
    ## All ALS BA9 -- 
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'All ALS'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
            # Construct the file path for the current cell type
            info_df <- read.csv(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_optimal_CombatSeq_',status,'_',brain_region,'_narval_2.csv'))
            info_df <- subset(info_df, celltype == i)
            par_n_topics <- info_df$n_genes
            par_n_epochs <- info_df$n_epochs
            file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_',par_n_topics,'_',par_n_epochs,'_narval_2_KNN_Bayesian.csv')
            # Check if the file exists
            if (file.exists(file_path)) {
                # Read the CSV file
                current_df <- read.csv(file_path)
                
                # Optionally, add the cell type as a new column to identify the source of each row
                current_df$dataset <- paste0(status, "_", brain_region)
                
                # Add the dataframe to the list
                df_list[[i]] <- current_df
            } else {
                warning(paste("File does not exist:", file_path))
            }
        }

        # Merge all dataframes into one
        final_All_ALS_BA9 <- do.call(rbind, df_list)
    ##  

    ###########################
    ## All FTLD BA9 -- 
    ###########################
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'All FTLD'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
            # Construct the file path for the current cell type
            info_df <- read.csv(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_optimal_CombatSeq_',status,'_',brain_region,'_narval_2.csv'))
            info_df <- subset(info_df, celltype == i)
            par_n_topics <- info_df$n_genes
            par_n_epochs <- info_df$n_epochs
            file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_',par_n_topics,'_',par_n_epochs,'_narval_2_KNN_Bayesian.csv')
            # Check if the file exists
            if (file.exists(file_path)) {
                # Read the CSV file
                current_df <- read.csv(file_path)
                
                # Optionally, add the cell type as a new column to identify the source of each row
                current_df$dataset <- paste0(status, "_", brain_region)
                
                # Add the dataframe to the list
                df_list[[i]] <- current_df
            } else {
                warning(paste("File does not exist:", file_path))
            }
        }

        # Merge all dataframes into one
        final_All_FTLD_BA9 <- do.call(rbind, df_list)
    ##  

    ###########################
    ## merge and plot BA9
    ###########################
    ## code for barplot
        BA9_bind_all <- rbind(final_SALS_BA9, final_C9ALS_BA9, final_SFTLD_BA9, final_C9FTLD_BA9, final_All_ALS_BA9, final_All_FTLD_BA9)
        write.csv(BA9_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_LOO_optimal_ABS_scaling_KNN_Bayesian_ETM.csv')
        
        BA9_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_LOO_optimal_ABS_scaling_KNN_Bayesian_ETM.csv')
        
        ## Calculate sample wise mean
        BA9_bind_all <- BA9_bind_all %>%
            group_by(donor, celltype, group, kNN_thresh) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )

        ## Calculate the mean and standard deviation
        summary_stats <- BA9_bind_all %>%
            group_by(celltype, group, kNN_thresh) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## set factor levels
        summary_stats$group <- factor(summary_stats$group, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All ALS", "All FTLD"))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))

        ## Plot
        ggplot(summary_stats, aes(x = group, y = median_accuracy, fill = group, colour = group)) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(kNN_thresh~celltype, scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/ETM_BA9_LOO_balanced_acc_barplot_optimal_abs_scaling_KNN_bayesian.pdf', height = 8, width = 13)

        
    
        ## Plot only Bayes 0.9 thresh
        summary_stats_lim <- subset(summary_stats, kNN_thresh == 0.9 )
        ggplot(summary_stats_lim, aes(x = group, y = median_accuracy, fill = group, colour = group)) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~celltype, scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA9_LOO_balanced_acc_barplot_optimal_abs_scaling_KNN_bayesian_09.pdf', height = 2, width = 13)

        summary_stats_lim_BA9 <- summary_stats_lim

    ##

    ###########################
    ## BA4 and BA9 0.9 Bayes plots 
    ## By cell type, group, and region
    ## By region
    ## Br major group
    ###########################
    ## Code
        ###################
        ## BA4
        ###################
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_LOO_optimal_ABS_scaling_KNN_Bayesian_ETM.csv')
        
        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, group, kNN_thresh) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )

        ## Calculate the mean and standard deviation
        summary_stats_BA4 <- BA4_bind_all %>%
            group_by(celltype, group, kNN_thresh) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## set factor levels
        summary_stats_BA4$group <- factor(summary_stats_BA4$group, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All ALS", "All FTLD"))
        summary_stats_BA4$celltype <- factor(summary_stats_BA4$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))

        ## Plot only Bayes 0.9 thresh
        summary_stats_BA4 <- subset(summary_stats_BA4, kNN_thresh == 0.9 )
        summary_stats_BA4$region <- "BA4"

        ###################
        ## BA9
        ###################
        BA9_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_LOO_optimal_ABS_scaling_KNN_Bayesian_ETM.csv')
        
        ## Calculate sample wise mean
        BA9_bind_all <- BA9_bind_all %>%
            group_by(donor, celltype, group, kNN_thresh) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )

        ## Calculate the mean and standard deviation
        summary_stats_BA9 <- BA9_bind_all %>%
            group_by(celltype, group, kNN_thresh) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## set factor levels
        summary_stats_BA9$group <- factor(summary_stats_BA9$group, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All ALS", "All FTLD"))
        summary_stats_BA9$celltype <- factor(summary_stats_BA9$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))

        ## Plot only Bayes 0.9 thresh
        summary_stats_BA9 <- subset(summary_stats_BA9, kNN_thresh == 0.9 )
        summary_stats_BA9$region <- "BA9"

        ###################
        ## Merge and plot
        ###################
        bind_all_test <- rbind(summary_stats_BA4, summary_stats_BA9)

        ## Plot by cell type
        ggplot(bind_all_test, aes(x = group, y = median_accuracy, fill = group, colour = group)) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(region ~ celltype, scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA4_BA9_ETM_LOO_balanced_acc_barplot_optimal_abs_scaling_KNN_bayesian_09.pdf', height = 3, width = 13)

        mean(bind_all_test$median_accuracy)
        
        
        ## Calculate mean of groups across brain regions
        summary_stats_all <- bind_all_test %>%
            group_by(group, region) %>%
            summarise(
                mean_accuracy = mean(median_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(median_accuracy),                  # Calculate the standard deviation of X0
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## Plot by cell type
        ggplot(summary_stats_all, aes(x = group, y = mean_accuracy, fill = group, colour = group, label = round(mean_accuracy, 2) )) + 
        geom_errorbar(aes(ymin=mean_accuracy - sd_accuracy, ymax=mean_accuracy + sd_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_text(colour = "white", size = 3, vjust = 1) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(region ~ ., scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 4, width = 3)

        
        ## By major celltype group
        bind_all_test <- data.frame(bind_all_test)
        bind_all_test$major_group[bind_all_test$celltype == "L2_L3"] <- "Excitatory"
        bind_all_test$major_group[bind_all_test$celltype == "L3_L5"] <- "Excitatory"
        bind_all_test$major_group[bind_all_test$celltype == "L4_L5"] <- "Excitatory"
        bind_all_test$major_group[bind_all_test$celltype == "L4_L6"] <- "Excitatory"
        bind_all_test$major_group[bind_all_test$celltype == "L5"] <- "Excitatory"
        bind_all_test$major_group[bind_all_test$celltype == "L5_L6"] <- "Excitatory"
        bind_all_test$major_group[bind_all_test$celltype == "L6"] <- "Excitatory"

        bind_all_test$major_group[bind_all_test$celltype == "5HT3aR"] <- "Inhibitory"
        bind_all_test$major_group[bind_all_test$celltype == "PV"] <- "Inhibitory"
        bind_all_test$major_group[bind_all_test$celltype == "Rosehip"] <- "Inhibitory"
        bind_all_test$major_group[bind_all_test$celltype == "SOM"] <- "Inhibitory"

        bind_all_test$major_group[bind_all_test$celltype == "Astro"] <- "non-neuronal"
        bind_all_test$major_group[bind_all_test$celltype == "Endo"] <- "non-neuronal"
        bind_all_test$major_group[bind_all_test$celltype == "Fibro"] <- "non-neuronal"
        bind_all_test$major_group[bind_all_test$celltype == "Micro"] <- "non-neuronal"
        bind_all_test$major_group[bind_all_test$celltype == "Mural"] <- "non-neuronal"
        bind_all_test$major_group[bind_all_test$celltype == "Oligo"] <- "non-neuronal"
        bind_all_test$major_group[bind_all_test$celltype == "OPC"] <- "non-neuronal"

        ## Calculate mean of groups across brain regions
        summary_stats_all <- bind_all_test %>%
            group_by(group, region, major_group) %>%
            summarise(
                mean_accuracy = mean(median_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(median_accuracy),                  # Calculate the standard deviation of X0
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## Plot
        ggplot(summary_stats_all, aes(x = group, y = mean_accuracy, fill = group, colour = group, label = round(mean_accuracy, 2) )) + 
        geom_errorbar(aes(ymin=mean_accuracy - sd_accuracy, ymax=mean_accuracy + sd_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_text(colour = "white", size = 3, vjust = 1) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(region ~ major_group, scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 4, width = 6)

        ## Calculate mean of groups across brain regions
        summary_stats_all <- bind_all_test %>%
            group_by(major_group) %>%
            summarise(
                mean_accuracy = mean(median_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(median_accuracy),                  # Calculate the standard deviation of X0
                .groups = "drop"                 # Drop the grouping after summarising
            )
        
    ##


##  

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################ Optimal Topic-epoch and LIME optimal topic: with bayesian 

## code
    ###########
    ## SALS BA4
    ###########
    ## Code
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'SALS'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
            # Construct the file path for the current cell type
            file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2_KNN_Bayesian_LIME_optimal.csv')
            # Check if the file exists
            if (file.exists(file_path)) {
                # Read the CSV file
                current_df <- read.csv(file_path)
                
                # Optionally, add the cell type as a new column to identify the source of each row
                current_df$dataset <- paste0(status, "_", brain_region)
                
                # Add the dataframe to the list
                df_list[[i]] <- current_df
            } else {
                warning(paste("File does not exist:", file_path))
            }
        }

        # Merge all dataframes into one
        final_SALS_BA4 <- do.call(rbind, df_list)
    ##

    ###########
    ## C9ALS BA4
    ###########
    ## Code
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'C9ALS'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
            # Construct the file path for the current cell type
            file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2_KNN_Bayesian_LIME_optimal.csv')
            # Check if the file exists
            if (file.exists(file_path)) {
                # Read the CSV file
                current_df <- read.csv(file_path)
                
                # Optionally, add the cell type as a new column to identify the source of each row
                current_df$dataset <- paste0(status, "_", brain_region)
                
                # Add the dataframe to the list
                df_list[[i]] <- current_df
            } else {
                warning(paste("File does not exist:", file_path))
            }
        }

        # Merge all dataframes into one
        final_C9ALS_BA4 <- do.call(rbind, df_list)
    ##

    ###########
    ## SFTLD BA4
    ###########
    ## Code
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'SFTLD'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
            # Construct the file path for the current cell type
            file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2_KNN_Bayesian_LIME_optimal.csv')
            # Check if the file exists
            if (file.exists(file_path)) {
                # Read the CSV file
                current_df <- read.csv(file_path)
                
                # Optionally, add the cell type as a new column to identify the source of each row
                current_df$dataset <- paste0(status, "_", brain_region)
                
                # Add the dataframe to the list
                df_list[[i]] <- current_df
            } else {
                warning(paste("File does not exist:", file_path))
            }
        }

        # Merge all dataframes into one
        final_SFTLD_BA4 <- do.call(rbind, df_list)
    ##

    ###########
    ## C9FTLD BA4
    ###########
    ## Code
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'C9FTLD'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
            # Construct the file path for the current cell type
            file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2_KNN_Bayesian_LIME_optimal.csv')
            # Check if the file exists
            if (file.exists(file_path)) {
                # Read the CSV file
                current_df <- read.csv(file_path)
                
                # Optionally, add the cell type as a new column to identify the source of each row
                current_df$dataset <- paste0(status, "_", brain_region)
                
                # Add the dataframe to the list
                df_list[[i]] <- current_df
            } else {
                warning(paste("File does not exist:", file_path))
            }
        }

        # Merge all dataframes into one
        final_C9FTLD_BA4 <- do.call(rbind, df_list)
    ##

    ###########
    ## All ALS BA4
    ###########
    ## Code
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'All ALS'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
            # Construct the file path for the current cell type
            file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2_KNN_Bayesian_LIME_optimal.csv')
            # Check if the file exists
            if (file.exists(file_path)) {
                # Read the CSV file
                current_df <- read.csv(file_path)
                
                # Optionally, add the cell type as a new column to identify the source of each row
                current_df$dataset <- paste0(status, "_", brain_region)
                
                # Add the dataframe to the list
                df_list[[i]] <- current_df
            } else {
                warning(paste("File does not exist:", file_path))
            }
        }

        # Merge all dataframes into one
        final_All_ALS_BA4 <- do.call(rbind, df_list)
    ##

    ###########
    ## All FTLD BA4
    ###########
    ## Code
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'All FTLD'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
            # Construct the file path for the current cell type
            file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2_KNN_Bayesian_LIME_optimal.csv')
            # Check if the file exists
            if (file.exists(file_path)) {
                # Read the CSV file
                current_df <- read.csv(file_path)
                
                # Optionally, add the cell type as a new column to identify the source of each row
                current_df$dataset <- paste0(status, "_", brain_region)
                
                # Add the dataframe to the list
                df_list[[i]] <- current_df
            } else {
                warning(paste("File does not exist:", file_path))
            }
        }

        # Merge all dataframes into one
        final_All_FTLD_BA4 <- do.call(rbind, df_list)
    ##

    ###########################
    ## merge and plot BA4
    ###########################
    ## code for barplot
        BA4_bind_all <- rbind(final_SALS_BA4, final_C9ALS_BA4, final_SFTLD_BA4, final_C9FTLD_BA4, final_All_ALS_BA4, final_All_FTLD_BA4)
        write.csv(BA4_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_LOO_optimal_ABS_scaling_KNN_Bayesian_ETM_LIME_optimal.csv')
        
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_LOO_optimal_ABS_scaling_KNN_Bayesian_ETM_LIME_optimal.csv')
        
        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, group, kNN_thresh) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )

        ## Calculate the mean and standard deviation
        summary_stats <- BA4_bind_all %>%
            group_by(celltype, group, kNN_thresh) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## set factor levels
        summary_stats$group <- factor(summary_stats$group, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All ALS", "All FTLD"))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))

        ## Plot
        ggplot(summary_stats, aes(x = group, y = median_accuracy, fill = group, colour = group)) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(kNN_thresh~celltype, scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/ETM_BA4_LOO_balanced_acc_barplot_optimal_abs_scaling_KNN_bayesian_LIME_optimal.pdf', height = 8, width = 13)

        ## Plot only Bayes 0.9 thresh
        summary_stats_lim <- subset(summary_stats, kNN_thresh == 0.9 )
        ggplot(summary_stats_lim, aes(x = group, y = median_accuracy, fill = group, colour = group)) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~celltype, scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA4_LOO_balanced_acc_barplot_optimal_abs_scaling_KNN_bayesian_09_LIME_optimal.pdf', height = 2, width = 13)

        summary_stats_lim_BA4 <- summary_stats_lim

    ##

    ###########
    ## SALS BA9
    ###########
    ## Code
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'SALS'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
            # Construct the file path for the current cell type
            file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2_KNN_Bayesian_LIME_optimal.csv')
            # Check if the file exists
            if (file.exists(file_path)) {
                # Read the CSV file
                current_df <- read.csv(file_path)
                
                # Optionally, add the cell type as a new column to identify the source of each row
                current_df$dataset <- paste0(status, "_", brain_region)
                
                # Add the dataframe to the list
                df_list[[i]] <- current_df
            } else {
                warning(paste("File does not exist:", file_path))
            }
        }

        # Merge all dataframes into one
        final_SALS_BA9 <- do.call(rbind, df_list)
    ##

    ###########
    ## C9ALS BA9
    ###########
    ## Code
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'C9ALS'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
            # Construct the file path for the current cell type
            file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2_KNN_Bayesian_LIME_optimal.csv')
            # Check if the file exists
            if (file.exists(file_path)) {
                # Read the CSV file
                current_df <- read.csv(file_path)
                
                # Optionally, add the cell type as a new column to identify the source of each row
                current_df$dataset <- paste0(status, "_", brain_region)
                
                # Add the dataframe to the list
                df_list[[i]] <- current_df
            } else {
                warning(paste("File does not exist:", file_path))
            }
        }

        # Merge all dataframes into one
        final_C9ALS_BA9 <- do.call(rbind, df_list)
    ##

    ###########
    ## SFTLD BA9
    ###########
    ## Code
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'SFTLD'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
            # Construct the file path for the current cell type
            file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2_KNN_Bayesian_LIME_optimal.csv')
            # Check if the file exists
            if (file.exists(file_path)) {
                # Read the CSV file
                current_df <- read.csv(file_path)
                
                # Optionally, add the cell type as a new column to identify the source of each row
                current_df$dataset <- paste0(status, "_", brain_region)
                
                # Add the dataframe to the list
                df_list[[i]] <- current_df
            } else {
                warning(paste("File does not exist:", file_path))
            }
        }

        # Merge all dataframes into one
        final_SFTLD_BA9 <- do.call(rbind, df_list)
    ##

    ###########
    ## C9FTLD BA9
    ###########
    ## Code
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'C9FTLD'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
            # Construct the file path for the current cell type
            file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2_KNN_Bayesian_LIME_optimal.csv')
            # Check if the file exists
            if (file.exists(file_path)) {
                # Read the CSV file
                current_df <- read.csv(file_path)
                
                # Optionally, add the cell type as a new column to identify the source of each row
                current_df$dataset <- paste0(status, "_", brain_region)
                
                # Add the dataframe to the list
                df_list[[i]] <- current_df
            } else {
                warning(paste("File does not exist:", file_path))
            }
        }

        # Merge all dataframes into one
        final_C9FTLD_BA9 <- do.call(rbind, df_list)
    ##

    ###########
    ## All ALS BA9
    ###########
    ## Code
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'All ALS'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
            # Construct the file path for the current cell type
            file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2_KNN_Bayesian_LIME_optimal.csv')
            # Check if the file exists
            if (file.exists(file_path)) {
                # Read the CSV file
                current_df <- read.csv(file_path)
                
                # Optionally, add the cell type as a new column to identify the source of each row
                current_df$dataset <- paste0(status, "_", brain_region)
                
                # Add the dataframe to the list
                df_list[[i]] <- current_df
            } else {
                warning(paste("File does not exist:", file_path))
            }
        }

        # Merge all dataframes into one
        final_All_ALS_BA9 <- do.call(rbind, df_list)
    ##

    ###########
    ## All FTLD BA9
    ###########
    ## Code
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'All FTLD'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()

        # Loop through each cell type
        for (i in cell_type) {
            # Construct the file path for the current cell type
            file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2_KNN_Bayesian_LIME_optimal.csv')
            # Check if the file exists
            if (file.exists(file_path)) {
                # Read the CSV file
                current_df <- read.csv(file_path)
                
                # Optionally, add the cell type as a new column to identify the source of each row
                current_df$dataset <- paste0(status, "_", brain_region)
                
                # Add the dataframe to the list
                df_list[[i]] <- current_df
            } else {
                warning(paste("File does not exist:", file_path))
            }
        }

        # Merge all dataframes into one
        final_All_FTLD_BA9 <- do.call(rbind, df_list)
    ##

    ###########################
    ## merge and plot BA9
    ###########################
    ## code for barplot
        #BA9_bind_all <- rbind(final_SALS_BA9, final_C9ALS_BA9, final_SFTLD_BA9, final_C9FTLD_BA9, final_All_ALS_BA9, final_All_FTLD_BA9)
        BA9_bind_all <- rbind(final_SALS_BA9, final_C9ALS_BA9, final_SFTLD_BA9, final_C9FTLD_BA9)
        write.csv(BA9_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_LOO_optimal_ABS_scaling_KNN_Bayesian_ETM_LIME_optimal.csv')
        
        BA9_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_LOO_optimal_ABS_scaling_KNN_Bayesian_ETM_LIME_optimal.csv')
        
        ## Calculate sample wise mean
        BA9_bind_all <- BA9_bind_all %>%
            group_by(donor, celltype, group, kNN_thresh) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )

        ## Calculate the mean and standard deviation
        summary_stats <- BA9_bind_all %>%
            group_by(celltype, group, kNN_thresh) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## set factor levels
        summary_stats$group <- factor(summary_stats$group, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All ALS", "All FTLD"))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))

        ## Plot
        ggplot(summary_stats, aes(x = group, y = median_accuracy, fill = group, colour = group)) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(kNN_thresh~celltype, scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/ETM_BA9_LOO_balanced_acc_barplot_optimal_abs_scaling_KNN_bayesian_LIME_optimal.pdf', height = 8, width = 13)

        ## Plot only Bayes 0.9 thresh
        summary_stats_lim <- subset(summary_stats, kNN_thresh == 0.9 )
        ggplot(summary_stats_lim, aes(x = group, y = median_accuracy, fill = group, colour = group)) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~celltype, scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA9_LOO_balanced_acc_barplot_optimal_abs_scaling_KNN_bayesian_09_LIME_optimal.pdf', height = 2, width = 13)

        summary_stats_lim_BA9 <- summary_stats_lim
    ##

    ###########################
    ## BA4 and BA9 0.9 Bayes plots 
    ## By cell type, group, and region
    ## By region
    ## Br major group
    ###########################
    ## Code
        ###################
        ## BA4
        ###################
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_LOO_optimal_ABS_scaling_KNN_Bayesian_ETM_LIME_optimal.csv')
        
        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, group, kNN_thresh) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )

        ## Calculate the mean and standard deviation
        summary_stats_BA4 <- BA4_bind_all %>%
            group_by(celltype, group, kNN_thresh) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## set factor levels
        summary_stats_BA4$group <- factor(summary_stats_BA4$group, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All ALS", "All FTLD"))
        summary_stats_BA4$celltype <- factor(summary_stats_BA4$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))

        ## Plot only Bayes 0.9 thresh
        summary_stats_BA4 <- subset(summary_stats_BA4, kNN_thresh == 0.9 )
        summary_stats_BA4$region <- "BA4"

        ###################
        ## BA9
        ###################
        BA9_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_LOO_optimal_ABS_scaling_KNN_Bayesian_ETM_LIME_optimal.csv')
        
        ## Calculate sample wise mean
        BA9_bind_all <- BA9_bind_all %>%
            group_by(donor, celltype, group, kNN_thresh) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )

        ## Calculate the mean and standard deviation
        summary_stats_BA9 <- BA9_bind_all %>%
            group_by(celltype, group, kNN_thresh) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## set factor levels
        summary_stats_BA9$group <- factor(summary_stats_BA9$group, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All ALS", "All FTLD"))
        summary_stats_BA9$celltype <- factor(summary_stats_BA9$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))

        ## Plot only Bayes 0.9 thresh
        summary_stats_BA9 <- subset(summary_stats_BA9, kNN_thresh == 0.9 )
        summary_stats_BA9$region <- "BA9"

        ###################
        ## Merge and plot
        ###################
        bind_all_test <- rbind(summary_stats_BA4, summary_stats_BA9)

        ## Plot by cell type
        ggplot(bind_all_test, aes(x = group, y = median_accuracy, fill = group, colour = group)) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(region ~ celltype, scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA4_BA9_ETM_LOO_balanced_acc_barplot_optimal_abs_scaling_KNN_bayesian_09_LIME_optimal.pdf', height = 3, width = 13)

        mean(bind_all_test$median_accuracy)
        
        
        ## Calculate mean of groups across brain regions
        summary_stats_all <- bind_all_test %>%
            group_by(group, region) %>%
            summarise(
                mean_accuracy = mean(median_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(median_accuracy),                  # Calculate the standard deviation of X0
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## Plot by cell type
        ggplot(summary_stats_all, aes(x = group, y = mean_accuracy, fill = group, colour = group, label = round(mean_accuracy, 2) )) + 
        geom_errorbar(aes(ymin=mean_accuracy - sd_accuracy, ymax=mean_accuracy + sd_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_text(colour = "white", size = 3, vjust = 1) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(region ~ ., scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 4, width = 3)

        
        ## By major celltype group
        bind_all_test <- data.frame(bind_all_test)
        bind_all_test$major_group[bind_all_test$celltype == "L2_L3"] <- "Excitatory"
        bind_all_test$major_group[bind_all_test$celltype == "L3_L5"] <- "Excitatory"
        bind_all_test$major_group[bind_all_test$celltype == "L4_L5"] <- "Excitatory"
        bind_all_test$major_group[bind_all_test$celltype == "L4_L6"] <- "Excitatory"
        bind_all_test$major_group[bind_all_test$celltype == "L5"] <- "Excitatory"
        bind_all_test$major_group[bind_all_test$celltype == "L5_L6"] <- "Excitatory"
        bind_all_test$major_group[bind_all_test$celltype == "L6"] <- "Excitatory"

        bind_all_test$major_group[bind_all_test$celltype == "5HT3aR"] <- "Inhibitory"
        bind_all_test$major_group[bind_all_test$celltype == "PV"] <- "Inhibitory"
        bind_all_test$major_group[bind_all_test$celltype == "Rosehip"] <- "Inhibitory"
        bind_all_test$major_group[bind_all_test$celltype == "SOM"] <- "Inhibitory"

        bind_all_test$major_group[bind_all_test$celltype == "Astro"] <- "non-neuronal"
        bind_all_test$major_group[bind_all_test$celltype == "Endo"] <- "non-neuronal"
        bind_all_test$major_group[bind_all_test$celltype == "Fibro"] <- "non-neuronal"
        bind_all_test$major_group[bind_all_test$celltype == "Micro"] <- "non-neuronal"
        bind_all_test$major_group[bind_all_test$celltype == "Mural"] <- "non-neuronal"
        bind_all_test$major_group[bind_all_test$celltype == "Oligo"] <- "non-neuronal"
        bind_all_test$major_group[bind_all_test$celltype == "OPC"] <- "non-neuronal"

        ## Calculate mean of groups across brain regions
        summary_stats_all <- bind_all_test %>%
            group_by(group, region, major_group) %>%
            summarise(
                mean_accuracy = mean(median_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(median_accuracy),                  # Calculate the standard deviation of X0
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## Plot
        ggplot(summary_stats_all, aes(x = group, y = mean_accuracy, fill = group, colour = group, label = round(mean_accuracy, 2) )) + 
        geom_errorbar(aes(ymin=mean_accuracy - sd_accuracy, ymax=mean_accuracy + sd_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_text(colour = "white", size = 3, vjust = 1) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(region ~ major_group, scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 4, width = 6)

        ## Calculate mean of groups across brain regions
        summary_stats_all <- bind_all_test %>%
            group_by(major_group) %>%
            summarise(
                mean_accuracy = mean(median_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(median_accuracy),                  # Calculate the standard deviation of X0
                .groups = "drop"                 # Drop the grouping after summarising
            )
        
    ##

    ###########################
    ## Not using this
    ###########################
    ## code
        ## Calculate sample wise mean
        BA4_bind_all <- final_SALS_BA4 %>%
            group_by(donor, celltype, group, kNN_thresh) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )

        ## Calculate the mean and standard deviation
        summary_stats <- BA4_bind_all %>%
            group_by(celltype, group, kNN_thresh) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## set factor levels
        summary_stats$group <- factor(summary_stats$group, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All ALS", "All FTLD"))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))
        res_LIME <- summary_stats
        res_LIME$region = "BA4"

        res_LIME <- subset(res_LIME,kNN_thresh == 0.9 )
        
        ####################
        ## Merge and plot
        ####################
        mean(res_LIME$median_accuracy)
        mean(res$median_accuracy)
        
        res$method <- "no_LIME"
        res_LIME$method <- "LIME"

        bind <- rbind(res, res_LIME)

        bind <- data.frame(bind)

        ## Plot by cell type
        ggplot(bind, aes(x = method, y = median_accuracy, fill = method, colour = method)) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "bottom",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(region ~ celltype, scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red")) +
        scale_colour_manual(values = c("orange", "red")) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp10.pdf', height = 3, width = 13)

        mean(bind_all_test$median_accuracy)


        ####################
        ## Explore LIME sets
        ####################

        i = "L4_L6"
        temp <- read.csv(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_topic_gene_info_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2_KNN_Bayesian_LIME_optimal.csv'))
        temp$topic <- paste0(temp$variable, "_", temp$donor)

        ## Compute jaccard simmilarity
        #gene_sets <- temp %>%
        #    group_by(topic) %>%
        #    summarise(genes = list(unique(value))) %>%
        #    deframe()

        gene_sets <- temp %>%
            group_by(topic) %>%
            slice_head(n = 100) %>%             # keep the first 100 rows per group
            summarise(genes = list(unique(value))) %>%
            deframe()

        length(gene_sets)


        #jaccard <- function(a, b) {
        #ci <- length(intersect(a, b))
        #cu <- length(union(a, b))
        #ci / cu
        #}

        #vars <- names(gene_sets)
        #mat <- outer(vars, vars, Vectorize(function(x, y) {
        #jaccard(gene_sets[[x]], gene_sets[[y]])
        #}))
        #rownames(mat) <- colnames(mat) <- vars

        overlap_coef <- function(a, b) {
        ci <- length(intersect(a, b))
        mi <- min(length(a), length(b))
        ci / mi
        }

        vars <- names(gene_sets)
        mat_overlap <- outer(vars, vars, Vectorize(function(x, y) {
        overlap_coef(gene_sets[[x]], gene_sets[[y]])
        }))
        rownames(mat_overlap) <- colnames(mat_overlap) <- vars

        mat_long <- melt(mat_overlap, varnames = c("Var1", "Var2"), value.name = "value")


        ggplot(mat_long, aes(x = Var2, y = Var1, fill = value)) +
        geom_tile(color = NA) +  
        theme_minimal(base_size = 12) +  
        theme(axis.text = element_blank(),
        axis.ticks = element_blank())     +   
        coord_fixed() +                          # keep tiles square
        scale_fill_gradient2(
            low = "blue", mid = "white", high = "red", midpoint = median(mat_long$value),
            name = "Similarity"
        )                                     # custom diverging palette
        
        
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp10.pdf', height = 13, width = 13)
    ##

##


############################################################################################################################################################################### HVGs and ETM 

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################ Gene and ETM with Bayesian

## code
    #########
    ## SALS BA4
    #########
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'SALS'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/genes_ETM_KNN_bayesian_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_SALS_BA4 <- do.call(rbind, df_list)
    ##

    #########
    ## C9ALS BA4
    #########
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'C9ALS'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/genes_ETM_KNN_bayesian_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_C9ALS_BA4 <- do.call(rbind, df_list)
    ##

    #########
    ## SFTLD BA4
    #########
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'SFTLD'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/genes_ETM_KNN_bayesian_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_SFTLD_BA4 <- do.call(rbind, df_list)
    ##

    #########
    ## C9FTLD BA4
    #########
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'C9FTLD'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/genes_ETM_KNN_bayesian_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_C9FTLD_BA4 <- do.call(rbind, df_list)
    ##

    #########
    ## All ALS BA4
    #########
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'All ALS'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/genes_ETM_KNN_bayesian_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_All_ALS_BA4 <- do.call(rbind, df_list)
    ##

    #########
    ## All FTLD BA4
    #########
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'All FTLD'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/genes_ETM_KNN_bayesian_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_All_FTLD_BA4 <- do.call(rbind, df_list)
    ##

    ###########################
    ## merge and plot BA4
    ###########################
    ## code for barplot
        BA4_bind_all <- rbind(final_SALS_BA4, final_C9ALS_BA4, final_SFTLD_BA4, final_C9FTLD_BA4, final_All_ALS_BA4, final_All_FTLD_BA4)
        write.csv(BA4_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_LOO_optimal_ABS_scaling_KNN_Bayesian_gene_ETM.csv')
        
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_LOO_optimal_ABS_scaling_KNN_Bayesian_gene_ETM.csv')
        
        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, group, kNN_thresh) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )

        ## Calculate the mean and standard deviation
        summary_stats <- BA4_bind_all %>%
            group_by(celltype, group, kNN_thresh) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## set factor levels
        summary_stats$group <- factor(summary_stats$group, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All ALS", "All FTLD"))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))

        ## Plot
        ggplot(summary_stats, aes(x = group, y = median_accuracy, fill = group, colour = group)) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(kNN_thresh~celltype, scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/ETM_BA4_LOO_balanced_acc_barplot_optimal_abs_scaling_KNN_bayesian_LIME_optimal.pdf', height = 8, width = 13)

        ## Plot only Bayes 0.9 thresh
        summary_stats_lim <- subset(summary_stats, kNN_thresh == 0.9 )
        ggplot(summary_stats_lim, aes(x = group, y = median_accuracy, fill = group, colour = group)) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~celltype, scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA4_LOO_balanced_acc_barplot_optimal_abs_scaling_KNN_bayesian_09_LIME_optimal.pdf', height = 2, width = 13)

        summary_stats_lim_BA4 <- summary_stats_lim
    ##

    #########
    ## SALS BA9
    #########
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'SALS'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/genes_ETM_KNN_bayesian_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_SALS_BA9 <- do.call(rbind, df_list)
    ##

    #########
    ## C9ALS BA9
    #########
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'C9ALS'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/genes_ETM_KNN_bayesian_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_C9ALS_BA9 <- do.call(rbind, df_list)
    ##

    #########
    ## SFTLD BA9
    #########
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'SFTLD'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/genes_ETM_KNN_bayesian_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_SFTLD_BA9 <- do.call(rbind, df_list)
    ##

    #########
    ## C9FTLD BA9
    #########
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'C9FTLD'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/genes_ETM_KNN_bayesian_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_C9FTLD_BA9 <- do.call(rbind, df_list)
    ##

    #########
    ## All ALS BA9
    #########
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'All ALS'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/genes_ETM_KNN_bayesian_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_All_ALS_BA9 <- do.call(rbind, df_list)
    ##

    #########
    ## All FTLD BA9
    #########
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'All FTLD'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/genes_ETM_KNN_bayesian_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_All_FTLD_BA9 <- do.call(rbind, df_list)
    ##

    ###########################
    ## merge and plot BA9
    ###########################
    ## code for barplot
        #BA9_bind_all <- rbind(final_SALS_BA9, final_C9ALS_BA9, final_SFTLD_BA9, final_C9FTLD_BA9, final_All_ALS_BA9, final_All_FTLD_BA9)
        BA9_bind_all <- rbind(final_SALS_BA9, final_C9ALS_BA9, final_SFTLD_BA9, final_C9FTLD_BA9)
        write.csv(BA9_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_LOO_optimal_ABS_scaling_KNN_Bayesian_gene_ETM.csv')
        
        BA9_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_LOO_optimal_ABS_scaling_KNN_Bayesian_gene_ETM.csv')
        
        ## Calculate sample wise mean
        BA9_bind_all <- BA9_bind_all %>%
            group_by(donor, celltype, group, kNN_thresh) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )

        ## Calculate the mean and standard deviation
        summary_stats <- BA9_bind_all %>%
            group_by(celltype, group, kNN_thresh) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## set factor levels
        summary_stats$group <- factor(summary_stats$group, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All ALS", "All FTLD"))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))

        ## Plot
        ggplot(summary_stats, aes(x = group, y = median_accuracy, fill = group, colour = group)) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(kNN_thresh~celltype, scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/ETM_BA9_LOO_balanced_acc_barplot_optimal_abs_scaling_KNN_bayesian_LIME_optimal.pdf', height = 8, width = 13)

        ## Plot only Bayes 0.9 thresh
        summary_stats_lim <- subset(summary_stats, kNN_thresh == 0.9 )
        ggplot(summary_stats_lim, aes(x = group, y = median_accuracy, fill = group, colour = group)) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~celltype, scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA9_LOO_balanced_acc_barplot_optimal_abs_scaling_KNN_bayesian_09_LIME_optimal.pdf', height = 2, width = 13)

        summary_stats_lim_BA9 <- summary_stats_lim
    ##

    ###########################
    ## BA4 and BA9 0.9 Bayes plots 
    ## By cell type, group, and region
    ## By region
    ## Br major group
    ###########################
    ## Code
        ###################
        ## BA4
        ###################
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_LOO_optimal_ABS_scaling_KNN_Bayesian_gene_ETM.csv')
        
        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, group, kNN_thresh) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )

        ## Calculate the mean and standard deviation
        summary_stats_BA4 <- BA4_bind_all %>%
            group_by(celltype, group, kNN_thresh) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## set factor levels
        summary_stats_BA4$group <- factor(summary_stats_BA4$group, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All ALS", "All FTLD"))
        summary_stats_BA4$celltype <- factor(summary_stats_BA4$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))

        ## Plot only Bayes 0.9 thresh
        summary_stats_BA4 <- subset(summary_stats_BA4, kNN_thresh == 0.9 )
        summary_stats_BA4$region <- "BA4"

        ###################
        ## BA9
        ###################
        BA9_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_LOO_optimal_ABS_scaling_KNN_Bayesian_gene_ETM.csv')
        
        ## Calculate sample wise mean
        BA9_bind_all <- BA9_bind_all %>%
            group_by(donor, celltype, group, kNN_thresh) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )

        ## Calculate the mean and standard deviation
        summary_stats_BA9 <- BA9_bind_all %>%
            group_by(celltype, group, kNN_thresh) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## set factor levels
        summary_stats_BA9$group <- factor(summary_stats_BA9$group, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All ALS", "All FTLD"))
        summary_stats_BA9$celltype <- factor(summary_stats_BA9$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))

        ## Plot only Bayes 0.9 thresh
        summary_stats_BA9 <- subset(summary_stats_BA9, kNN_thresh == 0.9 )
        summary_stats_BA9$region <- "BA9"

        ###################
        ## Merge and plot
        ###################
        bind_all_test <- rbind(summary_stats_BA4, summary_stats_BA9)

        ## Plot by cell type
        ggplot(bind_all_test, aes(x = group, y = median_accuracy, fill = group, colour = group)) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(region ~ celltype, scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA4_BA9_ETM_gene_LOO_balanced_acc_barplot_optimal_abs_scaling_KNN_bayesian_09.pdf', height = 3, width = 13)

        mean(bind_all_test$median_accuracy)
        
        
        ## Calculate mean of groups across brain regions
        summary_stats_all <- bind_all_test %>%
            group_by(group, region) %>%
            summarise(
                mean_accuracy = mean(median_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(median_accuracy),                  # Calculate the standard deviation of X0
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## Plot by cell type
        ggplot(summary_stats_all, aes(x = group, y = mean_accuracy, fill = group, colour = group, label = round(mean_accuracy, 2) )) + 
        geom_errorbar(aes(ymin=mean_accuracy - sd_accuracy, ymax=mean_accuracy + sd_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_text(colour = "white", size = 3, vjust = 1) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(region ~ ., scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 4, width = 3)

        
        ## By major celltype group
        bind_all_test <- data.frame(bind_all_test)
        bind_all_test$major_group[bind_all_test$celltype == "L2_L3"] <- "Excitatory"
        bind_all_test$major_group[bind_all_test$celltype == "L3_L5"] <- "Excitatory"
        bind_all_test$major_group[bind_all_test$celltype == "L4_L5"] <- "Excitatory"
        bind_all_test$major_group[bind_all_test$celltype == "L4_L6"] <- "Excitatory"
        bind_all_test$major_group[bind_all_test$celltype == "L5"] <- "Excitatory"
        bind_all_test$major_group[bind_all_test$celltype == "L5_L6"] <- "Excitatory"
        bind_all_test$major_group[bind_all_test$celltype == "L6"] <- "Excitatory"

        bind_all_test$major_group[bind_all_test$celltype == "5HT3aR"] <- "Inhibitory"
        bind_all_test$major_group[bind_all_test$celltype == "PV"] <- "Inhibitory"
        bind_all_test$major_group[bind_all_test$celltype == "Rosehip"] <- "Inhibitory"
        bind_all_test$major_group[bind_all_test$celltype == "SOM"] <- "Inhibitory"

        bind_all_test$major_group[bind_all_test$celltype == "Astro"] <- "non-neuronal"
        bind_all_test$major_group[bind_all_test$celltype == "Endo"] <- "non-neuronal"
        bind_all_test$major_group[bind_all_test$celltype == "Fibro"] <- "non-neuronal"
        bind_all_test$major_group[bind_all_test$celltype == "Micro"] <- "non-neuronal"
        bind_all_test$major_group[bind_all_test$celltype == "Mural"] <- "non-neuronal"
        bind_all_test$major_group[bind_all_test$celltype == "Oligo"] <- "non-neuronal"
        bind_all_test$major_group[bind_all_test$celltype == "OPC"] <- "non-neuronal"

        ## Calculate mean of groups across brain regions
        summary_stats_all <- bind_all_test %>%
            group_by(group, region, major_group) %>%
            summarise(
                mean_accuracy = mean(median_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(median_accuracy),                  # Calculate the standard deviation of X0
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## Plot
        ggplot(summary_stats_all, aes(x = group, y = mean_accuracy, fill = group, colour = group, label = round(mean_accuracy, 2) )) + 
        geom_errorbar(aes(ymin=mean_accuracy - sd_accuracy, ymax=mean_accuracy + sd_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_text(colour = "white", size = 3, vjust = 1) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(region ~ major_group, scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 4, width = 6)

        ## Calculate mean of groups across brain regions
        summary_stats_all <- bind_all_test %>%
            group_by(major_group) %>%
            summarise(
                mean_accuracy = mean(median_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(median_accuracy),                  # Calculate the standard deviation of X0
                .groups = "drop"                 # Drop the grouping after summarising
            )
        
    ##
##

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################ Gene and LIME ETM with Bayesian

## code

    #########
    ## SALS BA4
    #########
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'SALS'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_optimal_genes_ETM_KNN_bayesian_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_SALS_BA4 <- do.call(rbind, df_list)
    ##

    #########
    ## C9ALS BA4
    #########
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'C9ALS'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_optimal_genes_ETM_KNN_bayesian_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_C9ALS_BA4 <- do.call(rbind, df_list)
    ##

    #########
    ## SFTLD BA4
    #########
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'SFTLD'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_optimal_genes_ETM_KNN_bayesian_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_SFTLD_BA4 <- do.call(rbind, df_list)

        ## read in fibro
        i = "Fibro"
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/genes_ETM_KNN_bayesian_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')
        current_df <- read.csv(file_path)
        current_df$dataset <- paste0(status, "_", brain_region)

        ##Bind all
        final_SFTLD_BA4 <- rbind(final_SFTLD_BA4, current_df)

        print(length(unique(final_SFTLD_BA4$celltype)))
        


    ##
    ## missing fibro

    #########
    ## C9FTLD BA4
    #########
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'C9FTLD'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_optimal_genes_ETM_KNN_bayesian_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_C9FTLD_BA4 <- do.call(rbind, df_list)

        ## read in Oligo
        i = "Oligo"
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/genes_ETM_KNN_bayesian_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')
        current_df <- read.csv(file_path)
        current_df$dataset <- paste0(status, "_", brain_region)

        ##Bind all
        final_C9FTLD_BA4 <- rbind(final_C9FTLD_BA4, current_df)

        print(length(unique(final_C9FTLD_BA4$celltype)))
    ##
    ## missing oligo

    #########
    ## All ALS BA4
    #########
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'All ALS'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_optimal_genes_ETM_KNN_bayesian_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_All_ALS_BA4 <- do.call(rbind, df_list)
    ##

    #########
    ## All FTLD BA4
    #########
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'All FTLD'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_optimal_genes_ETM_KNN_bayesian_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_All_FTLD_BA4 <- do.call(rbind, df_list)

        ## read in Mural
        i = "Mural"
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/genes_ETM_KNN_bayesian_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')
        current_df <- read.csv(file_path)
        current_df$dataset <- paste0(status, "_", brain_region)

        ##Bind all
        final_All_FTLD_BA4 <- rbind(final_All_FTLD_BA4, current_df)

        print(length(unique(final_All_FTLD_BA4$celltype)))
    ##
    ## missing mural

    ###########################
    ## merge and plot BA4
    ###########################
    ## code for barplot
        BA4_bind_all <- rbind(final_SALS_BA4, final_C9ALS_BA4, final_SFTLD_BA4, final_C9FTLD_BA4, final_All_ALS_BA4, final_All_FTLD_BA4)
        write.csv(BA4_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_LOO_optimal_ABS_scaling_KNN_Bayesian_gene_ETM_LIME_optimal.csv')
        
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_LOO_optimal_ABS_scaling_KNN_Bayesian_gene_ETM_LIME_optimal.csv')
        
        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, group, kNN_thresh) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )

        ## Calculate the mean and standard deviation
        summary_stats <- BA4_bind_all %>%
            group_by(celltype, group, kNN_thresh) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## set factor levels
        summary_stats$group <- factor(summary_stats$group, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All ALS", "All FTLD"))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))

        ## Plot
        ggplot(summary_stats, aes(x = group, y = median_accuracy, fill = group, colour = group)) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(kNN_thresh~celltype, scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/ETM_BA4_LOO_balanced_acc_barplot_optimal_abs_scaling_KNN_bayesian_LIME_optimal.pdf', height = 8, width = 13)

        ## Plot only Bayes 0.9 thresh
        summary_stats_lim <- subset(summary_stats, kNN_thresh == 0.9 )
        ggplot(summary_stats_lim, aes(x = group, y = median_accuracy, fill = group, colour = group)) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~celltype, scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA4_LOO_balanced_acc_barplot_optimal_abs_scaling_KNN_bayesian_09_LIME_optimal.pdf', height = 2, width = 13)

        summary_stats_lim_BA4 <- summary_stats_lim
    ##

    #########
    ## SALS BA9
    #########
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'SALS'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_optimal_genes_ETM_KNN_bayesian_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_SALS_BA9 <- do.call(rbind, df_list)

        ## read in L3_L5
        i = "L3_L5"
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/genes_ETM_KNN_bayesian_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')
        current_df <- read.csv(file_path)
        current_df$dataset <- paste0(status, "_", brain_region)

        ##Bind all
        final_SALS_BA9 <- rbind(final_SALS_BA9, current_df)

        ## read in L2_L3
        i = "L2_L3"
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/genes_ETM_KNN_bayesian_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')
        current_df <- read.csv(file_path)
        current_df$dataset <- paste0(status, "_", brain_region)

        ##Bind all
        final_SALS_BA9 <- rbind(final_SALS_BA9, current_df)

        ## read in L5
        i = "L5"
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/genes_ETM_KNN_bayesian_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')
        current_df <- read.csv(file_path)
        current_df$dataset <- paste0(status, "_", brain_region)

        ##Bind all
        final_SALS_BA9 <- rbind(final_SALS_BA9, current_df)

        print(length(unique(final_SALS_BA9$celltype)))
    ##
    ## misisng L3_L5
    ## misisng L2_L3
    ## misisng L5

    #########
    ## C9ALS BA9
    #########
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'C9ALS'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_optimal_genes_ETM_KNN_bayesian_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_C9ALS_BA9 <- do.call(rbind, df_list)

        ## read in L3_L5
        i = "L3_L5"
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/genes_ETM_KNN_bayesian_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')
        current_df <- read.csv(file_path)
        current_df$dataset <- paste0(status, "_", brain_region)

        ##Bind all
        final_C9ALS_BA9 <- rbind(final_C9ALS_BA9, current_df)

        ## read in L2_L3
        i = "L2_L3"
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/genes_ETM_KNN_bayesian_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')
        current_df <- read.csv(file_path)
        current_df$dataset <- paste0(status, "_", brain_region)

        ##Bind all
        final_C9ALS_BA9 <- rbind(final_C9ALS_BA9, current_df)


        print(length(unique(final_C9ALS_BA9$celltype)))
    ##
    ## misisng L3_L5
    ## misisng L2_L3


    #########
    ## SFTLD BA9
    #########
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'SFTLD'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_optimal_genes_ETM_KNN_bayesian_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_SFTLD_BA9 <- do.call(rbind, df_list)

        ## read in L5
        i = "L5"
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/genes_ETM_KNN_bayesian_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')
        current_df <- read.csv(file_path)
        current_df$dataset <- paste0(status, "_", brain_region)

        ##Bind all
        final_SFTLD_BA9 <- rbind(final_SFTLD_BA9, current_df)

        print(length(unique(final_SFTLD_BA9$celltype)))
    ##
    ## missing L5

    #########
    ## C9FTLD BA9
    #########
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'C9FTLD'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_optimal_genes_ETM_KNN_bayesian_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_C9FTLD_BA9 <- do.call(rbind, df_list)
    ##

    #########
    ## All ALS BA9
    #########
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'All ALS'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_optimal_genes_ETM_KNN_bayesian_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_All_ALS_BA9 <- do.call(rbind, df_list)

        ## read in L3_L5
        i = "L3_L5"
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/genes_ETM_KNN_bayesian_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')
        current_df <- read.csv(file_path)
        current_df$dataset <- paste0(status, "_", brain_region)

        ##Bind all
        final_All_ALS_BA9 <- rbind(final_All_ALS_BA9, current_df)

        ## read in L2_L3
        i = "L2_L3"
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/genes_ETM_KNN_bayesian_model_report2_CombatSeq_',status,'_BA4_',i,'_narval_2.csv')
        current_df <- read.csv(file_path)
        current_df$dataset <- paste0(status, "_", brain_region)

        ##Bind all
        final_All_ALS_BA9 <- rbind(final_All_ALS_BA9, current_df)

        print(length(unique(final_All_ALS_BA9$celltype)))
    ##
    ## missing L3_L5
    ## missing L2_L3 -- NEED TO FIX.

    #########
    ## All FTLD BA9
    #########
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'All FTLD'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_optimal_genes_ETM_KNN_bayesian_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_All_FTLD_BA9 <- do.call(rbind, df_list)
    ##

    ###########################
    ## merge and plot BA9
    ###########################
    ## code for barplot
        BA9_bind_all <- rbind(final_SALS_BA9, final_C9ALS_BA9, final_SFTLD_BA9, final_C9FTLD_BA9, final_All_ALS_BA9, final_All_FTLD_BA9)
        write.csv(BA9_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_LOO_optimal_ABS_scaling_KNN_Bayesian_gene_ETM_LIME_optimal.csv')
        
        BA9_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_LOO_optimal_ABS_scaling_KNN_Bayesian_gene_ETM_LIME_optimal.csv')
        
        ## Calculate sample wise mean
        BA9_bind_all <- BA9_bind_all %>%
            group_by(donor, celltype, group, kNN_thresh) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )

        ## Calculate the mean and standard deviation
        summary_stats <- BA9_bind_all %>%
            group_by(celltype, group, kNN_thresh) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## set factor levels
        summary_stats$group <- factor(summary_stats$group, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All ALS", "All FTLD"))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))

        ## Plot
        ggplot(summary_stats, aes(x = group, y = median_accuracy, fill = group, colour = group)) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(kNN_thresh~celltype, scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/ETM_BA9_LOO_balanced_acc_barplot_optimal_abs_scaling_KNN_bayesian_LIME_optimal.pdf', height = 8, width = 13)

        ## Plot only Bayes 0.9 thresh
        summary_stats_lim <- subset(summary_stats, kNN_thresh == 0.9 )
        ggplot(summary_stats_lim, aes(x = group, y = median_accuracy, fill = group, colour = group)) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~celltype, scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA9_LOO_balanced_acc_barplot_optimal_abs_scaling_KNN_bayesian_09_LIME_optimal.pdf', height = 2, width = 13)

        summary_stats_lim_BA9 <- summary_stats_lim
    ##

    

    ###########################
    ## BA4 and BA9 0.9 Bayes plots 
    ## By cell type, group, and region
    ## By region
    ## Br major group
    ###########################
    ## Code
        ###################
        ## BA4
        ###################
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_LOO_optimal_ABS_scaling_KNN_Bayesian_gene_ETM_LIME_optimal.csv')
        unique(BA4_bind_all$kNN_thresh)

        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, group, kNN_thresh) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )

        ## Calculate the mean and standard deviation
        summary_stats_BA4 <- BA4_bind_all %>%
            group_by(celltype, group, kNN_thresh) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## set factor levels
        summary_stats_BA4$group <- factor(summary_stats_BA4$group, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All ALS", "All FTLD"))
        summary_stats_BA4$celltype <- factor(summary_stats_BA4$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))

        ## Plot only Bayes 0.9 thresh
        summary_stats_BA4 <- subset(summary_stats_BA4, kNN_thresh == 0.90 )
        summary_stats_BA4$region <- "BA4"

        ###################
        ## BA9
        ###################
        BA9_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_LOO_optimal_ABS_scaling_KNN_Bayesian_gene_ETM_LIME_optimal.csv')
        
        ## Calculate sample wise mean
        BA9_bind_all <- BA9_bind_all %>%
            group_by(donor, celltype, group, kNN_thresh) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )

        ## Calculate the mean and standard deviation
        summary_stats_BA9 <- BA9_bind_all %>%
            group_by(celltype, group, kNN_thresh) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## set factor levels
        summary_stats_BA9$group <- factor(summary_stats_BA9$group, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All ALS", "All FTLD"))
        summary_stats_BA9$celltype <- factor(summary_stats_BA9$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))

        ## Plot only Bayes 0.9 thresh
        summary_stats_BA9 <- subset(summary_stats_BA9, kNN_thresh == 0.90 )
        summary_stats_BA9$region <- "BA9"

        ###################
        ## Merge and plot
        ###################
        bind_all_test <- rbind(summary_stats_BA4, summary_stats_BA9)

        ## Plot by cell type
        ggplot(bind_all_test, aes(x = group, y = median_accuracy, fill = group, colour = group)) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(region ~ celltype, scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA4_BA9_ETM_gene_LOO_balanced_acc_barplot_optimal_abs_scaling_KNN_bayesian_09_LIME_optimal.pdf', height = 3, width = 13)

        mean(bind_all_test$median_accuracy)
        
        ## Calculate mean of groups across brain regions
        summary_stats_all <- bind_all_test %>%
            group_by(group, region) %>%
            summarise(
                mean_accuracy = mean(median_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(median_accuracy),                  # Calculate the standard deviation of X0
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## Plot by cell type
        ggplot(summary_stats_all, aes(x = group, y = mean_accuracy, fill = group, colour = group, label = round(mean_accuracy, 2) )) + 
        geom_errorbar(aes(ymin=mean_accuracy - sd_accuracy, ymax=mean_accuracy + sd_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_text(colour = "white", size = 3, vjust = 1) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(region ~ ., scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 4, width = 3)

        
        ## By major celltype group
        bind_all_test <- data.frame(bind_all_test)
        bind_all_test$major_group[bind_all_test$celltype == "L2_L3"] <- "Excitatory"
        bind_all_test$major_group[bind_all_test$celltype == "L3_L5"] <- "Excitatory"
        bind_all_test$major_group[bind_all_test$celltype == "L4_L5"] <- "Excitatory"
        bind_all_test$major_group[bind_all_test$celltype == "L4_L6"] <- "Excitatory"
        bind_all_test$major_group[bind_all_test$celltype == "L5"] <- "Excitatory"
        bind_all_test$major_group[bind_all_test$celltype == "L5_L6"] <- "Excitatory"
        bind_all_test$major_group[bind_all_test$celltype == "L6"] <- "Excitatory"

        bind_all_test$major_group[bind_all_test$celltype == "5HT3aR"] <- "Inhibitory"
        bind_all_test$major_group[bind_all_test$celltype == "PV"] <- "Inhibitory"
        bind_all_test$major_group[bind_all_test$celltype == "Rosehip"] <- "Inhibitory"
        bind_all_test$major_group[bind_all_test$celltype == "SOM"] <- "Inhibitory"

        bind_all_test$major_group[bind_all_test$celltype == "Astro"] <- "non-neuronal"
        bind_all_test$major_group[bind_all_test$celltype == "Endo"] <- "non-neuronal"
        bind_all_test$major_group[bind_all_test$celltype == "Fibro"] <- "non-neuronal"
        bind_all_test$major_group[bind_all_test$celltype == "Micro"] <- "non-neuronal"
        bind_all_test$major_group[bind_all_test$celltype == "Mural"] <- "non-neuronal"
        bind_all_test$major_group[bind_all_test$celltype == "Oligo"] <- "non-neuronal"
        bind_all_test$major_group[bind_all_test$celltype == "OPC"] <- "non-neuronal"

        ## Calculate mean of groups across brain regions
        summary_stats_all <- bind_all_test %>%
            group_by(group, region, major_group) %>%
            summarise(
                mean_accuracy = mean(median_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(median_accuracy),                  # Calculate the standard deviation of X0
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## Plot
        ggplot(summary_stats_all, aes(x = group, y = mean_accuracy, fill = group, colour = group, label = round(mean_accuracy, 2) )) + 
        geom_errorbar(aes(ymin=mean_accuracy - sd_accuracy, ymax=mean_accuracy + sd_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_text(colour = "white", size = 3, vjust = 1) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(region ~ major_group, scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 4, width = 6)

        ## Calculate mean of groups across brain regions
        summary_stats_all <- bind_all_test %>%
            group_by(major_group) %>%
            summarise(
                mean_accuracy = mean(median_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(median_accuracy),                  # Calculate the standard deviation of X0
                .groups = "drop"                 # Drop the grouping after summarising
            )
        
    ##
##


############################################################################################################################################################################### HVGs and NMF

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################ Gene and ETM with Bayesian

## code
    #########
    ## SALS BA4
    #########
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'SALS'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/genes_NMF_KNN_bayesian_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_SALS_BA4 <- do.call(rbind, df_list)
    ##

    #########
    ## C9ALS BA4
    #########
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'C9ALS'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/genes_NMF_KNN_bayesian_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_C9ALS_BA4 <- do.call(rbind, df_list)
    ##

    #########
    ## SFTLD BA4
    #########
    ## Missing fibro
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'SFTLD'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/genes_NMF_KNN_bayesian_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_SFTLD_BA4 <- do.call(rbind, df_list)
    ##

    #########
    ## C9FTLD BA4
    #########
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'C9FTLD'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/genes_NMF_KNN_bayesian_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_C9FTLD_BA4 <- do.call(rbind, df_list)
    ##

    #########
    ## All ALS BA4
    #########
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'All ALS'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/genes_NMF_KNN_bayesian_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_All_ALS_BA4 <- do.call(rbind, df_list)
    ##

    #########
    ## All FTLD BA4
    #########
    ## missing fibro
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'All FTLD'
        brain_region <- 'BA4'

        # Create an empty list to store dataframes
        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/genes_NMF_KNN_bayesian_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_All_FTLD_BA4 <- do.call(rbind, df_list)
    ##

    ###########################
    ## merge and plot BA4
    ###########################
    ## code for barplot
        BA4_bind_all <- rbind(final_SALS_BA4, final_C9ALS_BA4, final_SFTLD_BA4, final_C9FTLD_BA4, final_All_ALS_BA4, final_All_FTLD_BA4)
        write.csv(BA4_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_LOO_optimal_ABS_scaling_KNN_Bayesian_gene_NMF.csv')
        
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_LOO_optimal_ABS_scaling_KNN_Bayesian_gene_NMF.csv')
        
        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, group, kNN_thresh) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )

        ## Calculate the mean and standard deviation
        summary_stats <- BA4_bind_all %>%
            group_by(celltype, group, kNN_thresh) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## set factor levels
        summary_stats$group <- factor(summary_stats$group, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All ALS", "All FTLD"))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))

        ## Plot
        ggplot(summary_stats, aes(x = group, y = median_accuracy, fill = group, colour = group)) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(kNN_thresh~celltype, scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/NMF_BA4_LOO_balanced_acc_barplot_optimal_abs_scaling_KNN_bayesian_LIME_optimal.pdf', height = 8, width = 13)

        ## Plot only Bayes 0.9 thresh
        summary_stats_lim <- subset(summary_stats, kNN_thresh == 0.9 )
        ggplot(summary_stats_lim, aes(x = group, y = median_accuracy, fill = group, colour = group)) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~celltype, scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA4_LOO_balanced_acc_barplot_optimal_abs_scaling_KNN_bayesian_09_LIME_optimal.pdf', height = 2, width = 13)

        summary_stats_lim_BA4 <- summary_stats_lim
    ##

    #########
    ## SALS BA9
    #########
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'SALS'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/genes_NMF_KNN_bayesian_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_SALS_BA9 <- do.call(rbind, df_list)
    ##

    #########
    ## C9ALS BA9
    #########
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'C9ALS'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/genes_NMF_KNN_bayesian_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_C9ALS_BA9 <- do.call(rbind, df_list)
    ##

    #########
    ## SFTLD BA9
    #########
    ## missing L5
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'SFTLD'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/genes_NMF_KNN_bayesian_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_SFTLD_BA9 <- do.call(rbind, df_list)
    ##

    #########
    ## C9FTLD BA9
    #########
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'C9FTLD'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/genes_NMF_KNN_bayesian_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_C9FTLD_BA9 <- do.call(rbind, df_list)
    ##

    #########
    ## All ALS BA9
    #########
    ## Missing L2_L3; L3_L5
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'All ALS'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/genes_NMF_KNN_bayesian_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_All_ALS_BA9 <- do.call(rbind, df_list)
    ##

    #########
    ## All FTLD BA9
    #########
    ## Missing L5
    ## code
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        status <- 'All FTLD'
        brain_region <- 'BA9'

        # Create an empty list to store dataframes
        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/genes_NMF_KNN_bayesian_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            
            # Optionally, add the cell type as a new column to identify the source of each row
            current_df$dataset <- paste0(status, "_", brain_region)
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        final_All_FTLD_BA9 <- do.call(rbind, df_list)
    ##

    ###########################
    ## merge and plot BA9
    ###########################
    ## code for barplot
        BA9_bind_all <- rbind(final_SALS_BA9, final_C9ALS_BA9, final_SFTLD_BA9, final_C9FTLD_BA9, final_All_ALS_BA9, final_All_FTLD_BA9)
        #BA9_bind_all <- rbind(final_SALS_BA9, final_C9ALS_BA9, final_SFTLD_BA9, final_C9FTLD_BA9, final_All_ALS_BA9)
        write.csv(BA9_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_LOO_optimal_ABS_scaling_KNN_Bayesian_gene_NMF.csv')
        
        BA9_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_LOO_optimal_ABS_scaling_KNN_Bayesian_gene_NMF.csv')
        
        ## Calculate sample wise mean
        BA9_bind_all <- BA9_bind_all %>%
            group_by(donor, celltype, group, kNN_thresh) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )

        ## Calculate the mean and standard deviation
        summary_stats <- BA9_bind_all %>%
            group_by(celltype, group, kNN_thresh) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## set factor levels
        summary_stats$group <- factor(summary_stats$group, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All ALS", "All FTLD"))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))

        ## Plot
        ggplot(summary_stats, aes(x = group, y = median_accuracy, fill = group, colour = group)) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(kNN_thresh~celltype, scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/NMF_BA9_LOO_balanced_acc_barplot_optimal_abs_scaling_KNN_bayesian_LIME_optimal.pdf', height = 8, width = 13)

        ## Plot only Bayes 0.9 thresh
        summary_stats_lim <- subset(summary_stats, kNN_thresh == 0.9 )
        ggplot(summary_stats_lim, aes(x = group, y = median_accuracy, fill = group, colour = group)) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~celltype, scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA9_LOO_balanced_acc_barplot_optimal_abs_scaling_KNN_bayesian_09_LIME_optimal.pdf', height = 2, width = 13)

        summary_stats_lim_BA9 <- summary_stats_lim
    ##

    ###########################
    ## BA4 and BA9 0.9 Bayes plots 
    ## By cell type, group, and region
    ## By region
    ## Br major group
    ###########################
    ## Code
        ###################
        ## BA4
        ###################
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_LOO_optimal_ABS_scaling_KNN_Bayesian_gene_NMF.csv')
        
        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, group, kNN_thresh) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )

        ## Calculate the mean and standard deviation
        summary_stats_BA4 <- BA4_bind_all %>%
            group_by(celltype, group, kNN_thresh) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## set factor levels
        summary_stats_BA4$group <- factor(summary_stats_BA4$group, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All ALS", "All FTLD"))
        summary_stats_BA4$celltype <- factor(summary_stats_BA4$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))

        ## Plot only Bayes 0.9 thresh
        summary_stats_BA4 <- subset(summary_stats_BA4, kNN_thresh == 0.9 )
        summary_stats_BA4$region <- "BA4"

        ###################
        ## BA9
        ###################
        BA9_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_LOO_optimal_ABS_scaling_KNN_Bayesian_gene_NMF.csv')
        
        ## Calculate sample wise mean
        BA9_bind_all <- BA9_bind_all %>%
            group_by(donor, celltype, group, kNN_thresh) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )

        ## Calculate the mean and standard deviation
        summary_stats_BA9 <- BA9_bind_all %>%
            group_by(celltype, group, kNN_thresh) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## set factor levels
        summary_stats_BA9$group <- factor(summary_stats_BA9$group, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All ALS", "All FTLD"))
        summary_stats_BA9$celltype <- factor(summary_stats_BA9$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))

        ## Plot only Bayes 0.9 thresh
        summary_stats_BA9 <- subset(summary_stats_BA9, kNN_thresh == 0.9 )
        summary_stats_BA9$region <- "BA9"

        ###################
        ## Merge and plot
        ###################
        bind_all_test <- rbind(summary_stats_BA4, summary_stats_BA9)

        ## Plot by cell type
        ggplot(bind_all_test, aes(x = group, y = median_accuracy, fill = group, colour = group)) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(region ~ celltype, scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA4_BA9_NMF_gene_LOO_balanced_acc_barplot_optimal_abs_scaling_KNN_bayesian_09.pdf', height = 3, width = 13)

        mean(bind_all_test$median_accuracy)
        
        
        ## Calculate mean of groups across brain regions
        summary_stats_all <- bind_all_test %>%
            group_by(group, region) %>%
            summarise(
                mean_accuracy = mean(median_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(median_accuracy),                  # Calculate the standard deviation of X0
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## Plot by cell type
        ggplot(summary_stats_all, aes(x = group, y = mean_accuracy, fill = group, colour = group, label = round(mean_accuracy, 2) )) + 
        geom_errorbar(aes(ymin=mean_accuracy - sd_accuracy, ymax=mean_accuracy + sd_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_text(colour = "white", size = 3, vjust = 1) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(region ~ ., scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 4, width = 3)

        
        ## By major celltype group
        bind_all_test <- data.frame(bind_all_test)
        bind_all_test$major_group[bind_all_test$celltype == "L2_L3"] <- "Excitatory"
        bind_all_test$major_group[bind_all_test$celltype == "L3_L5"] <- "Excitatory"
        bind_all_test$major_group[bind_all_test$celltype == "L4_L5"] <- "Excitatory"
        bind_all_test$major_group[bind_all_test$celltype == "L4_L6"] <- "Excitatory"
        bind_all_test$major_group[bind_all_test$celltype == "L5"] <- "Excitatory"
        bind_all_test$major_group[bind_all_test$celltype == "L5_L6"] <- "Excitatory"
        bind_all_test$major_group[bind_all_test$celltype == "L6"] <- "Excitatory"

        bind_all_test$major_group[bind_all_test$celltype == "5HT3aR"] <- "Inhibitory"
        bind_all_test$major_group[bind_all_test$celltype == "PV"] <- "Inhibitory"
        bind_all_test$major_group[bind_all_test$celltype == "Rosehip"] <- "Inhibitory"
        bind_all_test$major_group[bind_all_test$celltype == "SOM"] <- "Inhibitory"

        bind_all_test$major_group[bind_all_test$celltype == "Astro"] <- "non-neuronal"
        bind_all_test$major_group[bind_all_test$celltype == "Endo"] <- "non-neuronal"
        bind_all_test$major_group[bind_all_test$celltype == "Fibro"] <- "non-neuronal"
        bind_all_test$major_group[bind_all_test$celltype == "Micro"] <- "non-neuronal"
        bind_all_test$major_group[bind_all_test$celltype == "Mural"] <- "non-neuronal"
        bind_all_test$major_group[bind_all_test$celltype == "Oligo"] <- "non-neuronal"
        bind_all_test$major_group[bind_all_test$celltype == "OPC"] <- "non-neuronal"

        ## Calculate mean of groups across brain regions
        summary_stats_all <- bind_all_test %>%
            group_by(group, region, major_group) %>%
            summarise(
                mean_accuracy = mean(median_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(median_accuracy),                  # Calculate the standard deviation of X0
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## Plot
        ggplot(summary_stats_all, aes(x = group, y = mean_accuracy, fill = group, colour = group, label = round(mean_accuracy, 2) )) + 
        geom_errorbar(aes(ymin=mean_accuracy - sd_accuracy, ymax=mean_accuracy + sd_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_text(colour = "white", size = 3, vjust = 1) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(region ~ major_group, scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "darkred", "darkblue")) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 4, width = 6)

        ## Calculate mean of groups across brain regions
        summary_stats_all <- bind_all_test %>%
            group_by(major_group) %>%
            summarise(
                mean_accuracy = mean(median_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(median_accuracy),                  # Calculate the standard deviation of X0
                .groups = "drop"                 # Drop the grouping after summarising
            )
        
    ##
##





########################################################################################################################################################## Not using beyond this point, for now. 





################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################ 35 epoch 100 topics 10000 genes
## code 
    
    ## BA4
    
    ###########################
    ## SALS BA4
    ###########################
    ## code
        cell_type <- c("L3_L5", "L4_L6", "L4_L6", "L4_L5", "L5_L6", "L6", "PV", "5HT3aR", 
        "Rosehip", "SOM", "Astro", "Astro", "OPC", "Micro",
        "Mural", "Endo", "Fibro", "L5")

        status <- 'SALS'
        brain_region <- 'BA4'

        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
                # Construct the file path for the current cell type
                file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_100_35_narval_2.csv')

                # Check if the file exists
                if (file.exists(file_path)) {
                    # Read the CSV file
                    current_df <- read.csv(file_path)
                    df_list[[i]] <- current_df
                } else {
                    warning(paste("File does not exist:", file_path))
                }
            }
        

        # Merge all dataframes into one
        final_SALS_BA4 <- do.call(rbind, df_list)
    ##

    ###########################
    ## C9ALS BA4
    ###########################
    ## code
        cell_type <- c("L3_L5", "L4_L6", "L4_L6", "L4_L5", "L5_L6", "L6", "PV", "5HT3aR", 
        "Rosehip", "SOM", "Astro", "Astro", "OPC", "Micro",
        "Mural", "Endo", "Fibro", "L5")

        status <- 'C9ALS'
        brain_region <- 'BA4'

        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
                # Construct the file path for the current cell type
                file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_100_35_narval_2.csv')

                # Check if the file exists
                if (file.exists(file_path)) {
                    # Read the CSV file
                    current_df <- read.csv(file_path)
                    df_list[[i]] <- current_df
                } else {
                    warning(paste("File does not exist:", file_path))
                }
            }

        # Merge all dataframes into one
        final_C9ALS_BA4 <- do.call(rbind, df_list)
    ##

    ###########################
    ## SFTLD BA4
    ###########################
    ## code
        cell_type <- c("L3_L5", "L4_L6", "L4_L6", "L4_L5", "L5_L6", "L6", "PV", "5HT3aR", 
        "Rosehip", "SOM", "Astro", "Astro", "OPC", "Micro",
        "Mural", "Endo", "Fibro", "L5")

        status <- 'SFTLD'
        brain_region <- 'BA4'

        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
                # Construct the file path for the current cell type
                file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_100_35_narval_2.csv')

                # Check if the file exists
                if (file.exists(file_path)) {
                    # Read the CSV file
                    current_df <- read.csv(file_path)
                    df_list[[i]] <- current_df
                } else {
                    warning(paste("File does not exist:", file_path))
                }
            }

        # Merge all dataframes into one
        final_SFTLD_BA4 <- do.call(rbind, df_list)
    ##

    ###########################
    ## C9FTLD BA4
    ###########################
    ## code
        cell_type <- c("L3_L5", "L4_L6", "L4_L6", "L4_L5", "L5_L6", "L6", "PV", "5HT3aR", 
        "Rosehip", "SOM", "Astro", "Astro", "OPC", "Micro",
        "Mural", "Endo", "Fibro", "L5")

        status <- 'C9FTLD'
        brain_region <- 'BA4'

        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
                # Construct the file path for the current cell type
                file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_100_35_narval_2.csv')

                # Check if the file exists
                if (file.exists(file_path)) {
                    # Read the CSV file
                    current_df <- read.csv(file_path)
                    df_list[[i]] <- current_df
                } else {
                    warning(paste("File does not exist:", file_path))
                }
            }

        # Merge all dataframes into one
        final_C9FTLD_BA4 <- do.call(rbind, df_list)
    ##

    ###########################
    ## merge and plot BA4
    ###########################
    ## code for barplot
        
        BA4_bind_all <- rbind(final_SALS_BA4, final_C9ALS_BA4, final_SFTLD_BA4, final_C9FTLD_BA4)
        
        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, group) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy),        
                .groups = "drop"               
            )

        ## Calculate the mean and standard deviation
        summary_stats <- BA4_bind_all %>%
            group_by(celltype, group) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## set factor levels
        summary_stats$group <- factor(summary_stats$group, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD"))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L4_L6','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Astro', 'OPC','T_Cell'  ))


        ## Plot
        ggplot(summary_stats, aes(x = group, y = median_accuracy, fill = group, colour = group)) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~celltype, scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "black")) #+
        #coord_cartesian(ylim=c(0.5,1.0)) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 2, width = 13)
    ##
    
    ## BA9
    
    ###########################
    ## SALS BA9
    ###########################
    ## code
        cell_type <- c("L3_L5", "L4_L6", "L4_L6", "L4_L5", "L5_L6", "L6", "PV", "5HT3aR", 
        "Rosehip", "SOM", "Astro", "Astro", "OPC", "Micro",
        "Mural", "Endo", "Fibro", "L5")

        status <- 'SALS'
        brain_region <- 'BA9'

        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
                # Construct the file path for the current cell type
                file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_100_35_narval_2.csv')

                # Check if the file exists
                if (file.exists(file_path)) {
                    # Read the CSV file
                    current_df <- read.csv(file_path)
                    df_list[[i]] <- current_df
                } else {
                    warning(paste("File does not exist:", file_path))
                }
            }
        

        # Merge all dataframes into one
        final_SALS_BA9 <- do.call(rbind, df_list)
    ##

    ###########################
    ## C9ALS BA9
    ###########################
    ## code
        cell_type <- c("L3_L5", "L4_L6", "L4_L6", "L4_L5", "L5_L6", "L6", "PV", "5HT3aR", 
        "Rosehip", "SOM", "Astro", "Astro", "OPC", "Micro",
        "Mural", "Endo", "Fibro", "L5")

        status <- 'C9ALS'
        brain_region <- 'BA9'

        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
                # Construct the file path for the current cell type
                file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_100_35_narval_2.csv')

                # Check if the file exists
                if (file.exists(file_path)) {
                    # Read the CSV file
                    current_df <- read.csv(file_path)
                    df_list[[i]] <- current_df
                } else {
                    warning(paste("File does not exist:", file_path))
                }
            }

        # Merge all dataframes into one
        final_C9ALS_BA9 <- do.call(rbind, df_list)
    ##

    ###########################
    ## SFTLD BA9
    ###########################
    ## code
        cell_type <- c("L3_L5", "L4_L6", "L4_L6", "L4_L5", "L5_L6", "L6", "PV", "5HT3aR", 
        "Rosehip", "SOM", "Astro", "Astro", "OPC", "Micro",
        "Mural", "Endo", "Fibro", "L5")

        status <- 'SFTLD'
        brain_region <- 'BA9'

        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
                # Construct the file path for the current cell type
                file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_100_35_narval_2.csv')

                # Check if the file exists
                if (file.exists(file_path)) {
                    # Read the CSV file
                    current_df <- read.csv(file_path)
                    df_list[[i]] <- current_df
                } else {
                    warning(paste("File does not exist:", file_path))
                }
            }

        # Merge all dataframes into one
        final_SFTLD_BA9 <- do.call(rbind, df_list)
    ##

    ###########################
    ## C9FTLD BA9
    ###########################
    ## code
        cell_type <- c("L3_L5", "L4_L6", "L4_L6", "L4_L5", "L5_L6", "L6", "PV", "5HT3aR", 
        "Rosehip", "SOM", "Astro", "Astro", "OPC", "Micro",
        "Mural", "Endo", "Fibro", "L5")

        status <- 'C9FTLD'
        brain_region <- 'BA9'

        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
                # Construct the file path for the current cell type
                file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_100_35_narval_2.csv')

                # Check if the file exists
                if (file.exists(file_path)) {
                    # Read the CSV file
                    current_df <- read.csv(file_path)
                    df_list[[i]] <- current_df
                } else {
                    warning(paste("File does not exist:", file_path))
                }
            }

        # Merge all dataframes into one
        final_C9FTLD_BA9 <- do.call(rbind, df_list)
    ##

    ###########################
    ## merge and plot BA9
    ###########################
    ## code for barplot
        BA9_bind_all <- rbind(final_SALS_BA9, final_C9ALS_BA9, final_SFTLD_BA9, final_C9FTLD_BA9)
        
        ## Calculate sample wise mean
        BA9_bind_all <- BA9_bind_all %>%
            group_by(donor, celltype, group) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy),        
                .groups = "drop"               
            )

        ## Calculate the mean and standard deviation
        summary_stats <- BA9_bind_all %>%
            group_by(celltype, group) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## set factor levels
        summary_stats$group <- factor(summary_stats$group, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD"))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L4_L6','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Astro', 'OPC','T_Cell'  ))


        ## Plot
        ggplot(summary_stats, aes(x = group, y = median_accuracy, fill = group, colour = group)) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~celltype, scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "black")) #+
        #coord_cartesian(ylim=c(0.5,1.0)) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 2, width = 13)
    ##
##

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################ 35 epoch 100 topics 5000 genes
## code 
    
    ## BA4
    
    ###########################
    ## SALS BA4
    ###########################
    ## code
        cell_type <- c("L3_L5", "L4_L6", "L4_L6", "L4_L5", "L5_L6", "L6", "PV", "5HT3aR", 
        "Rosehip", "SOM", "Astro", "Astro", "OPC", "Micro",
        "Mural", "Endo", "Fibro", "L5")

        status <- 'SALS'
        brain_region <- 'BA4'

        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
                # Construct the file path for the current cell type
                file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_100_35_narval_2_5000.csv')

                # Check if the file exists
                if (file.exists(file_path)) {
                    # Read the CSV file
                    current_df <- read.csv(file_path)
                    df_list[[i]] <- current_df
                } else {
                    warning(paste("File does not exist:", file_path))
                }
            }
        

        # Merge all dataframes into one
        final_SALS_BA4 <- do.call(rbind, df_list)
    ##

    ###########################
    ## C9ALS BA4
    ###########################
    ## code
        cell_type <- c("L3_L5", "L4_L6", "L4_L6", "L4_L5", "L5_L6", "L6", "PV", "5HT3aR", 
        "Rosehip", "SOM", "Astro", "Astro", "OPC", "Micro",
        "Mural", "Endo", "Fibro", "L5")

        status <- 'C9ALS'
        brain_region <- 'BA4'

        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
                # Construct the file path for the current cell type
                file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_100_35_narval_2_5000.csv')

                # Check if the file exists
                if (file.exists(file_path)) {
                    # Read the CSV file
                    current_df <- read.csv(file_path)
                    df_list[[i]] <- current_df
                } else {
                    warning(paste("File does not exist:", file_path))
                }
            }

        # Merge all dataframes into one
        final_C9ALS_BA4 <- do.call(rbind, df_list)
    ##

    ###########################
    ## SFTLD BA4
    ###########################
    ## code
        cell_type <- c("L3_L5", "L4_L6", "L4_L6", "L4_L5", "L5_L6", "L6", "PV", "5HT3aR", 
        "Rosehip", "SOM", "Astro", "Astro", "OPC", "Micro",
        "Mural", "Endo", "Fibro", "L5")

        status <- 'SFTLD'
        brain_region <- 'BA4'

        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
                # Construct the file path for the current cell type
                file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_100_35_narval_2_5000.csv')

                # Check if the file exists
                if (file.exists(file_path)) {
                    # Read the CSV file
                    current_df <- read.csv(file_path)
                    df_list[[i]] <- current_df
                } else {
                    warning(paste("File does not exist:", file_path))
                }
            }

        # Merge all dataframes into one
        final_SFTLD_BA4 <- do.call(rbind, df_list)
    ##

    ###########################
    ## C9FTLD BA4
    ###########################
    ## code
        cell_type <- c("L3_L5", "L4_L6", "L4_L6", "L4_L5", "L5_L6", "L6", "PV", "5HT3aR", 
        "Rosehip", "SOM", "Astro", "Astro", "OPC", "Micro",
        "Mural", "Endo", "Fibro", "L5")

        status <- 'C9FTLD'
        brain_region <- 'BA4'

        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
                # Construct the file path for the current cell type
                file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_100_35_narval_2_5000.csv')

                # Check if the file exists
                if (file.exists(file_path)) {
                    # Read the CSV file
                    current_df <- read.csv(file_path)
                    df_list[[i]] <- current_df
                } else {
                    warning(paste("File does not exist:", file_path))
                }
            }

        # Merge all dataframes into one
        final_C9FTLD_BA4 <- do.call(rbind, df_list)
    ##

    ###########################
    ## merge and plot BA4
    ###########################
    ## code for barplot
        
        BA4_bind_all <- rbind(final_SALS_BA4, final_C9ALS_BA4, final_SFTLD_BA4, final_C9FTLD_BA4)
        
        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, group) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy),        
                .groups = "drop"               
            )

        ## Calculate the mean and standard deviation
        summary_stats <- BA4_bind_all %>%
            group_by(celltype, group) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## set factor levels
        summary_stats$group <- factor(summary_stats$group, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD"))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L4_L6','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Astro', 'OPC','T_Cell'  ))


        ## Plot
        ggplot(summary_stats, aes(x = group, y = median_accuracy, fill = group, colour = group)) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~celltype, scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "black")) #+
        #coord_cartesian(ylim=c(0.5,1.0)) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 2, width = 13)
    ##
    
    ## BA9
    
    ###########################
    ## SALS BA9
    ###########################
    ## code
        cell_type <- c("L3_L5", "L4_L6", "L4_L6", "L4_L5", "L5_L6", "L6", "PV", "5HT3aR", 
        "Rosehip", "SOM", "Astro", "Astro", "OPC", "Micro",
        "Mural", "Endo", "Fibro", "L5")

        status <- 'SALS'
        brain_region <- 'BA9'

        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
                # Construct the file path for the current cell type
                file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_100_35_narval_2.csv')

                # Check if the file exists
                if (file.exists(file_path)) {
                    # Read the CSV file
                    current_df <- read.csv(file_path)
                    df_list[[i]] <- current_df
                } else {
                    warning(paste("File does not exist:", file_path))
                }
            }
        

        # Merge all dataframes into one
        final_SALS_BA9 <- do.call(rbind, df_list)
    ##

    ###########################
    ## C9ALS BA9
    ###########################
    ## code
        cell_type <- c("L3_L5", "L4_L6", "L4_L6", "L4_L5", "L5_L6", "L6", "PV", "5HT3aR", 
        "Rosehip", "SOM", "Astro", "Astro", "OPC", "Micro",
        "Mural", "Endo", "Fibro", "L5")

        status <- 'C9ALS'
        brain_region <- 'BA9'

        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
                # Construct the file path for the current cell type
                file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_100_35_narval_2.csv')

                # Check if the file exists
                if (file.exists(file_path)) {
                    # Read the CSV file
                    current_df <- read.csv(file_path)
                    df_list[[i]] <- current_df
                } else {
                    warning(paste("File does not exist:", file_path))
                }
            }

        # Merge all dataframes into one
        final_C9ALS_BA9 <- do.call(rbind, df_list)
    ##

    ###########################
    ## SFTLD BA9
    ###########################
    ## code
        cell_type <- c("L3_L5", "L4_L6", "L4_L6", "L4_L5", "L5_L6", "L6", "PV", "5HT3aR", 
        "Rosehip", "SOM", "Astro", "Astro", "OPC", "Micro",
        "Mural", "Endo", "Fibro", "L5")

        status <- 'SFTLD'
        brain_region <- 'BA9'

        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
                # Construct the file path for the current cell type
                file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_100_35_narval_2.csv')

                # Check if the file exists
                if (file.exists(file_path)) {
                    # Read the CSV file
                    current_df <- read.csv(file_path)
                    df_list[[i]] <- current_df
                } else {
                    warning(paste("File does not exist:", file_path))
                }
            }

        # Merge all dataframes into one
        final_SFTLD_BA9 <- do.call(rbind, df_list)
    ##

    ###########################
    ## C9FTLD BA9
    ###########################
    ## code
        cell_type <- c("L3_L5", "L4_L6", "L4_L6", "L4_L5", "L5_L6", "L6", "PV", "5HT3aR", 
        "Rosehip", "SOM", "Astro", "Astro", "OPC", "Micro",
        "Mural", "Endo", "Fibro", "L5")

        status <- 'C9FTLD'
        brain_region <- 'BA9'

        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
                # Construct the file path for the current cell type
                file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_100_35_narval_2.csv')

                # Check if the file exists
                if (file.exists(file_path)) {
                    # Read the CSV file
                    current_df <- read.csv(file_path)
                    df_list[[i]] <- current_df
                } else {
                    warning(paste("File does not exist:", file_path))
                }
            }

        # Merge all dataframes into one
        final_C9FTLD_BA9 <- do.call(rbind, df_list)
    ##

    ###########################
    ## merge and plot BA9
    ###########################
    ## code for barplot
        BA9_bind_all <- rbind(final_SALS_BA9, final_C9ALS_BA9, final_SFTLD_BA9, final_C9FTLD_BA9)
        
        ## Calculate sample wise mean
        BA9_bind_all <- BA9_bind_all %>%
            group_by(donor, celltype, group) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy),        
                .groups = "drop"               
            )

        ## Calculate the mean and standard deviation
        summary_stats <- BA9_bind_all %>%
            group_by(celltype, group) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## set factor levels
        summary_stats$group <- factor(summary_stats$group, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD"))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L4_L6','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Astro', 'OPC','T_Cell'  ))


        ## Plot
        ggplot(summary_stats, aes(x = group, y = median_accuracy, fill = group, colour = group)) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~celltype, scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "black")) #+
        #coord_cartesian(ylim=c(0.5,1.0)) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 2, width = 13)
    ##
##





############################################################################################################################################################################### ETM process documentation


################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ ETM SFTLD L4_L6 BA9 parameter sweep-- Need to re do this. 
## code 
    # Create an empty list to store dataframes
    #topic_list = c(100, 150, 200, 250, 300)
    #epoch_list = c(35)
    #out_path = f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_report2_{par_prep}_{par_status}_{par_brain_region}_{par_keep_cell_type}_{par_n_topics}_{par_n_epochs}_narval_2.csv"
    #i = "L5_L6"
    #file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_report2_CombatSeq_SALS_BA4_',i,'_100_35_narval_2.csv')
    #current_df <- read.csv(file_path)

    ###########################
    ## 35-100-10000
    ###########################
    
    cell_type <- c("L3_L5", "L4_L6", "L4_L6", "L4_L5", "L5_L6", "L6", "PV", "5HT3aR", 
    "Rosehip", "SOM", "Astro", "Astro", "OPC", "Micro",
    "Mural", "Endo", "Fibro", "L5")

    status <- 'C9FTLD'
    brain_region <- 'BA9'

    df_list <- list()
    
    # Loop through each cell type
    for (i in cell_type) {
            # Construct the file path for the current cell type
            file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_100_35_narval_2.csv')

            # Check if the file exists
            if (file.exists(file_path)) {
                # Read the CSV file
                current_df <- read.csv(file_path)
                df_list[[i]] <- current_df
            } else {
                warning(paste("File does not exist:", file_path))
            }
        }
    

    # Merge all dataframes into one
    final_SALS_BA4 <- do.call(rbind, df_list)

    temp <- final_SALS_BA4 %>% dplyr::select(donor, celltype, test_accuracy)


    ###########################
    ## 35-300-10000
    ###########################
    cell_type <- c("L3_L5", "L4_L6", "L4_L6", "L4_L5", "L5_L6", "L6", "PV", "5HT3aR", 
    "Rosehip", "SOM", "Astro", "Astro", "OPC", "Micro",
    "Mural", "Endo", "Fibro", "L5")

    status <- 'C9ALS'
    brain_region <- 'BA9'

    df_list <- list()
    
    # Loop through each cell type
    for (i in cell_type) {
            # Construct the file path for the current cell type
            file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_300_35_narval_2.csv')

            # Check if the file exists
            if (file.exists(file_path)) {
                # Read the CSV file
                current_df <- read.csv(file_path)
                df_list[[i]] <- current_df
            } else {
                warning(paste("File does not exist:", file_path))
            }
        }


    ###########################
    ## 35-100-5000
    ###########################
    cell_type <- c("L3_L5", "L4_L6", "L4_L6", "L4_L5", "L5_L6", "L6", "PV", "5HT3aR", 
    "Rosehip", "SOM", "Astro", "Astro", "OPC", "Micro",
    "Mural", "Endo", "Fibro", "L5")

    status <- 'C9ALS'
    brain_region <- 'BA9'

    df_list <- list()
    
    # Loop through each cell type
    for (i in cell_type) {
            # Construct the file path for the current cell type
            file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_100_35_narval_2_5000.csv')

            # Check if the file exists
            if (file.exists(file_path)) {
                # Read the CSV file
                current_df <- read.csv(file_path)
                df_list[[i]] <- current_df
            } else {
                warning(paste("File does not exist:", file_path))
            }
        }

    ###########################
    ## 35-300-5000
    ###########################
    cell_type <- c("L3_L5", "L4_L6", "L4_L6", "L4_L5", "L5_L6", "L6", "PV", "5HT3aR", 
    "Rosehip", "SOM", "Astro", "Astro", "OPC", "Micro",
    "Mural", "Endo", "Fibro", "L5")

    status <- 'SFTLD'
    brain_region <- 'BA9'

    df_list <- list()
    
    # Loop through each cell type
    for (i in cell_type) {
            # Construct the file path for the current cell type
            file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_report2_CombatSeq_',status,'_',brain_region,'_',i,'_300_35_narval_2_5000.csv')

            # Check if the file exists
            if (file.exists(file_path)) {
                # Read the CSV file
                current_df <- read.csv(file_path)
                df_list[[i]] <- current_df
            } else {
                warning(paste("File does not exist:", file_path))
            }
        }
    

    # Merge all dataframes into one
    final_SALS_BA4 <- do.call(rbind, df_list)

    temp <- final_SALS_BA4 %>% dplyr::select(donor, celltype, test_accuracy)


    ###########################
    ## LOO bayesian
    ###########################
    cell_type <- c("L3_L5", "L4_L6", "L4_L6", "L4_L5", "L5_L6", "L6", "PV", "5HT3aR", 
    "Rosehip", "SOM", "Astro", "Astro", "OPC", "Micro",
    "Mural", "Endo", "Fibro", "L5")

    status <- 'C9FTLD'
    brain_region <- 'BA9'

    df_list <- list()
    
    # Loop through each cell type
    for (i in cell_type) {
            # Construct the file path for the current cell type
            file_path <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_CombatSeq_All ALS_BA9_",i,"_0_abs_case_control_narval_2.csv")

            # Check if the file exists
            if (file.exists(file_path)) {
                # Read the CSV file
                current_df <- read.csv(file_path)
                df_list[[i]] <- current_df
            } else {
                warning(paste("File does not exist:", file_path))
            }
        }
    

    # Merge all dataframes into one
    final_SALS_BA4 <- do.call(rbind, df_list)

    temp <- final_SALS_BA4 %>% dplyr::select(donor, celltype, test_accuracy)




    write.csv(mean_deltas_sorted, file = paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_CombatSeq_All ALS_BA9_",celltype2,"_0_abs_case_control_narval_2.csv"), sep = ",")





    ## code for barplot
        
        ## Calculate sample wise mean
        final_SALS_BA4 <- final_SALS_BA4 %>%
            group_by(donor, celltype) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy),        
                .groups = "drop"               
            )

        ## Calculate the mean and standard deviation
        summary_stats <- final_SALS_BA4 %>%
            group_by(celltype) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L4_L6','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Astro', 'OPC','T_Cell'  ))

        test <- subset(final_SALS_BA4, celltype == "L4_L5")
        test <- data.frame(test)
        
        ## Plot
        ggplot(summary_stats, aes(x = celltype, y = median_accuracy, fill = "black", colour = "black", label = round(median_accuracy, 4))) + 
        geom_errorbar(aes(ymin=median_accuracy - mad_accuracy, ymax=median_accuracy + mad_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_text(colour = "black")+
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.y = element_text(face = "bold")
        ) +
        #facet_grid(~n_topic, scales = "free_x") +
        ylab("Accuracy") + xlab("n_epoch")
        #scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) +
        #scale_colour_manual(values = c("orange", "red", "blue", "purple", "black")) #+
        #coord_cartesian(ylim=c(0.5,1.0)) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 4, width = 4)
##


############################################################################################### END supplemental figures. 



################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ DNN test classification accuracies (Main figure)
## code
    ###########################
    ## USE THIS ONE: Main figure accuracy heatmap: TEST only
    ###########################
    ## code
        BA4_bind_all_test <- subset(BA4_bind_all, test == "test")
        BA9_bind_all_test <- subset(BA9_bind_all, test == "test")
        bind_all_test <- rbind(BA4_bind_all_test, BA9_bind_all_test)

        ggplot(bind_all_test, aes(x = celltype, y = condition, fill = accuracy)) +
        theme_bw() + 
        geom_tile() +
        geom_text(aes(label = round(accuracy, 2)), size = 2.5, colour = "black") +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
            axis.text.y = element_text(colour = c("orange", "red", "blue", "purple", "black")),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            strip.background =element_rect(fill="white", colour = "white")
        ) +
        facet_grid(region ~ .) +
        ylab("Accuracy") +
        scale_x_discrete(expand = c(0,0)) +
        scale_y_discrete(expand = c(0,0)) +
        scale_fill_gradientn(colors = brewer.pal(9, "RdYlGn"))   # Use the "RdYlBu" palette for continuous data
        #scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) +
        #scale_colour_manual(values = c("orange", "red", "blue", "purple", "black"))
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 2.5, width = 9)
    ##

    ###########################
    ## DO NOT USE THIS ONE: Main figure accuracy heatmap: MEAN across test and validation only
    ###########################
    ## code
        ## Bind and process BA4
        summary_df_BA4 <- BA4_bind_all %>%
            group_by(celltype, condition, region) %>%
            summarize(mean_accuracy = mean(accuracy, na.rm = TRUE)) %>% as.data.frame()

        summary_df_BA9 <- BA9_bind_all %>%
            group_by(celltype, condition, region) %>%
            summarize(mean_accuracy = mean(accuracy, na.rm = TRUE)) %>% as.data.frame()

        bind_all_test <- rbind(summary_df_BA4, summary_df_BA9)

        ggplot(bind_all_test, aes(x = celltype, y = condition, fill = mean_accuracy)) +
        theme_bw() + 
        geom_tile() +
        geom_text(aes(label = round(mean_accuracy, 2)), size = 2.5, colour = "black") +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
            axis.text.y = element_text(colour = c("orange", "red", "blue", "purple", "black")),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            strip.background =element_rect(fill="white", colour = "white")
        ) +
        facet_grid(region ~ .) +
        scale_x_discrete(expand = c(0,0)) +
        scale_y_discrete(expand = c(0,0)) +
        scale_fill_gradientn(colors = brewer.pal(9, "RdYlGn"))   # Use the "RdYlBu" palette for continuous data
        #scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) +
        #scale_colour_manual(values = c("orange", "red", "blue", "purple", "black"))
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 2.5, width = 9)
    ##

    ###########################
    ## DO NOT USE THIS ONE: Main figure accuracy heatmap: MEDIAN across test and validation only
    ###########################
    ## code
        ## Bind and process BA4
        summary_df_BA4 <- BA4_bind_all %>%
            group_by(celltype, condition, region) %>%
            summarize(median_accuracy = median(accuracy, na.rm = TRUE)) %>% as.data.frame()

        summary_df_BA9 <- BA9_bind_all %>%
            group_by(celltype, condition, region) %>%
            summarize(median_accuracy = median(accuracy, na.rm = TRUE)) %>% as.data.frame()

        bind_all_test <- rbind(summary_df_BA4, summary_df_BA9)

        ggplot(bind_all_test, aes(x = celltype, y = condition, fill = median_accuracy)) +
        theme_bw() + 
        geom_tile() +
        geom_text(aes(label = round(median_accuracy, 2)), size = 2.5, colour = "black") +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
            axis.text.y = element_text(colour = c("orange", "red", "blue", "purple", "black")),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            strip.background =element_rect(fill="white", colour = "white")
        ) +
        facet_grid(region ~ .) +
        scale_x_discrete(expand = c(0,0)) +
        scale_y_discrete(expand = c(0,0)) +
        scale_fill_gradientn(colors = brewer.pal(9, "RdYlGn"))   # Use the "RdYlBu" palette for continuous data
        #scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) +
        #scale_colour_manual(values = c("orange", "red", "blue", "purple", "black"))
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 2.5, width = 9)
    ##
##



################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ DNN classification accuracy correlations

## code
    ###########################
    ## SALS BA4
    ###########################
    ## code
        ## define variables
        condition = "SALS" #one ofSALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')

        result_list <- list()

        i = "5HT3aR"

        for(i in unique(cell_type_list)){
            
            file_1 <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/test_meta_",condition,"_",region,"_",i,"_pineda_narval.csv")
            temp_1 <- read.delim(file_1, header = T, sep = ",")
            length_1 <- nrow(temp_1)
            file_2 <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/train_meta_",condition,"_",region,"_",i,"_pineda_narval.csv")
            temp_2 <- read.delim(file_2, header = T, sep = ",")
            length_2 <- nrow(temp_2)
            n_cells = length_1 + length_2
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval.csv')
            temp <- read.delim(file, header = T, sep = ",")
            HVGs <- as.numeric(temp$n_HVGs)
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            all_accuracies <- mean(all_accuracies)
            result_df <- data.frame(accuracy = all_accuracies, n_cells = n_cells, HVGs = HVGs)
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
        condition = "C9ALS" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file_1 <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/test_meta_",condition,"_",region,"_",i,"_pineda_narval.csv")
            temp_1 <- read.delim(file_1, header = T, sep = ",")
            length_1 <- nrow(temp_1)
            file_2 <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/train_meta_",condition,"_",region,"_",i,"_pineda_narval.csv")
            temp_2 <- read.delim(file_2, header = T, sep = ",")
            length_2 <- nrow(temp_2)
            n_cells = length_1 + length_2
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval.csv')
            temp <- read.delim(file, header = T, sep = ",")
            HVGs <- as.numeric(temp$n_HVGs)
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            all_accuracies <- mean(all_accuracies)
            result_df <- data.frame(accuracy = all_accuracies, n_cells = n_cells, HVGs = HVGs)
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
        condition = "SFTLD" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file_1 <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/test_meta_",condition,"_",region,"_",i,"_pineda_narval.csv")
            temp_1 <- read.delim(file_1, header = T, sep = ",")
            length_1 <- nrow(temp_1)
            file_2 <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/train_meta_",condition,"_",region,"_",i,"_pineda_narval.csv")
            temp_2 <- read.delim(file_2, header = T, sep = ",")
            length_2 <- nrow(temp_2)
            n_cells = length_1 + length_2
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval.csv')
            temp <- read.delim(file, header = T, sep = ",")
            HVGs <- as.numeric(temp$n_HVGs)
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            all_accuracies <- mean(all_accuracies)
            result_df <- data.frame(accuracy = all_accuracies, n_cells = n_cells, HVGs = HVGs)
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
        condition = "C9FTLD" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file_1 <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/test_meta_",condition,"_",region,"_",i,"_pineda_narval.csv")
            temp_1 <- read.delim(file_1, header = T, sep = ",")
            length_1 <- nrow(temp_1)
            file_2 <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/train_meta_",condition,"_",region,"_",i,"_pineda_narval.csv")
            temp_2 <- read.delim(file_2, header = T, sep = ",")
            length_2 <- nrow(temp_2)
            n_cells = length_1 + length_2
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval.csv')
            temp <- read.delim(file, header = T, sep = ",")
            HVGs <- as.numeric(temp$n_HVGs)
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            all_accuracies <- mean(all_accuracies)
            result_df <- data.frame(accuracy = all_accuracies, n_cells = n_cells, HVGs = HVGs)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region
            result_list[[i]] <- result_df

        }

        final_df_C9FTLD_BA4 <- do.call(rbind, result_list)
    ##

    ###########################
    ## all_conditions BA4
    ###########################
    ## code
        ## define variables
        condition = "all_conditions" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')
        #cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file_1 <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/test_meta_",condition,"_",region,"_",i,"_pineda_narval.csv")
            temp_1 <- read.delim(file_1, header = T, sep = ",")
            length_1 <- nrow(temp_1)
            file_2 <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/train_meta_",condition,"_",region,"_",i,"_pineda_narval.csv")
            temp_2 <- read.delim(file_2, header = T, sep = ",")
            length_2 <- nrow(temp_2)
            n_cells = length_1 + length_2
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval.csv')
            temp <- read.delim(file, header = T, sep = ",")
            HVGs <- as.numeric(temp$n_HVGs)
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            all_accuracies <- mean(all_accuracies)
            result_df <- data.frame(accuracy = all_accuracies, n_cells = n_cells, HVGs = HVGs)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region
            result_list[[i]] <- result_df

        }

        final_df_all_conditions_BA4 <- do.call(rbind, result_list)
    ##

    ###########################
    ## merge and plot BA4
    ###########################
    ## code
        BA4_bind_all <- rbind(final_df_SALS_BA4, final_df_C9ALS_BA4, final_df_SFTLD_BA4, final_df_C9FTLD_BA4, final_df_all_conditions_BA4)
        BA4_bind_all$condition[BA4_bind_all$condition == "all_conditions"] <- "All conditions"

        ## set factor levels
        BA4_bind_all$condition <- factor(BA4_bind_all$condition, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All conditions"))
        BA4_bind_all$celltype <- factor(BA4_bind_all$celltype, levels = c('L4_L6','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Astro', 'OPC','T_Cell'  ))


        ## correlation with N HVGs
        ggplot(BA4_bind_all, aes(x = HVGs, y = accuracy)) + 
        theme_bw() + 
        geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "solid") + # Add line of best fit
        geom_point(colour = "darkgrey") +
        stat_cor(method = "pearson", colour = "black", label.y = 1,) +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(colour = "black"),
            axis.text.y = element_text(colour = "black"),
            axis.title.x = element_text(face = "bold"),
            axis.title.y = element_text(face = "bold")
        ) +
        ylab("Accuracy") + 
        xlab("n HVGs") 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 3, width = 3)

        ## correlation with N Cells
        ggplot(BA4_bind_all, aes(x = n_cells, y = accuracy)) + 
        theme_bw() + 
        geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "solid") + # Add line of best fit
        geom_point(colour = "darkgrey") +
        stat_cor(method = "pearson", colour = "black", label.y = 1,) +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(colour = "black"),
            axis.text.y = element_text(colour = "black"),
            axis.title.x = element_text(face = "bold"),
            axis.title.y = element_text(face = "bold")
        ) +
        ylab("Accuracy") + 
        xlab("n cells") 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 3, width = 3)
    ##


    ###########################
    ## SALS BA9
    ###########################
    ## code
        ## define variables
        condition = "SALS" #one ofSALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file_1 <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/test_meta_",condition,"_",region,"_",i,"_pineda_narval.csv")
            temp_1 <- read.delim(file_1, header = T, sep = ",")
            length_1 <- nrow(temp_1)
            file_2 <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/train_meta_",condition,"_",region,"_",i,"_pineda_narval.csv")
            temp_2 <- read.delim(file_2, header = T, sep = ",")
            length_2 <- nrow(temp_2)
            n_cells = length_1 + length_2
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval.csv')
            temp <- read.delim(file, header = T, sep = ",")
            HVGs <- as.numeric(temp$n_HVGs)
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            all_accuracies <- mean(all_accuracies)
            result_df <- data.frame(accuracy = all_accuracies, n_cells = n_cells, HVGs = HVGs)
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
        condition = "C9ALS" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file_1 <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/test_meta_",condition,"_",region,"_",i,"_pineda_narval.csv")
            temp_1 <- read.delim(file_1, header = T, sep = ",")
            length_1 <- nrow(temp_1)
            file_2 <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/train_meta_",condition,"_",region,"_",i,"_pineda_narval.csv")
            temp_2 <- read.delim(file_2, header = T, sep = ",")
            length_2 <- nrow(temp_2)
            n_cells = length_1 + length_2
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval.csv')
            temp <- read.delim(file, header = T, sep = ",")
            HVGs <- as.numeric(temp$n_HVGs)
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            all_accuracies <- mean(all_accuracies)
            result_df <- data.frame(accuracy = all_accuracies, n_cells = n_cells, HVGs = HVGs)
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
        condition = "SFTLD" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file_1 <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/test_meta_",condition,"_",region,"_",i,"_pineda_narval.csv")
            temp_1 <- read.delim(file_1, header = T, sep = ",")
            length_1 <- nrow(temp_1)
            file_2 <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/train_meta_",condition,"_",region,"_",i,"_pineda_narval.csv")
            temp_2 <- read.delim(file_2, header = T, sep = ",")
            length_2 <- nrow(temp_2)
            n_cells = length_1 + length_2
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval.csv')
            temp <- read.delim(file, header = T, sep = ",")
            HVGs <- as.numeric(temp$n_HVGs)
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            all_accuracies <- mean(all_accuracies)
            result_df <- data.frame(accuracy = all_accuracies, n_cells = n_cells, HVGs = HVGs)
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
        condition = "C9FTLD" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file_1 <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/test_meta_",condition,"_",region,"_",i,"_pineda_narval.csv")
            temp_1 <- read.delim(file_1, header = T, sep = ",")
            length_1 <- nrow(temp_1)
            file_2 <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/train_meta_",condition,"_",region,"_",i,"_pineda_narval.csv")
            temp_2 <- read.delim(file_2, header = T, sep = ",")
            length_2 <- nrow(temp_2)
            n_cells = length_1 + length_2
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval.csv')
            temp <- read.delim(file, header = T, sep = ",")
            HVGs <- as.numeric(temp$n_HVGs)
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            all_accuracies <- mean(all_accuracies)
            result_df <- data.frame(accuracy = all_accuracies, n_cells = n_cells, HVGs = HVGs)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region
            result_list[[i]] <- result_df

        }

        final_df_C9FTLD_BA9 <- do.call(rbind, result_list)
    ##

    ###########################
    ## all_conditions BA9
    ###########################
    ## code
        ## define variables
        condition = "all_conditions" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L4_L6', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Astro', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file_1 <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/test_meta_",condition,"_",region,"_",i,"_pineda_narval.csv")
            temp_1 <- read.delim(file_1, header = T, sep = ",")
            length_1 <- nrow(temp_1)
            file_2 <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/train_meta_",condition,"_",region,"_",i,"_pineda_narval.csv")
            temp_2 <- read.delim(file_2, header = T, sep = ",")
            length_2 <- nrow(temp_2)
            n_cells = length_1 + length_2
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval.csv')
            temp <- read.delim(file, header = T, sep = ",")
            HVGs <- as.numeric(temp$n_HVGs)
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            all_accuracies <- mean(all_accuracies)
            result_df <- data.frame(accuracy = all_accuracies, n_cells = n_cells, HVGs = HVGs)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region
            result_list[[i]] <- result_df

        }

        final_df_all_conditions_BA9 <- do.call(rbind, result_list)
    ##

    ###########################
    ## merge and plot BA9
    ###########################
    ## code
        BA9_bind_all <- rbind(final_df_SALS_BA9, final_df_C9ALS_BA9, final_df_SFTLD_BA9, final_df_C9FTLD_BA9, final_df_all_conditions_BA9)
        BA9_bind_all$condition[BA9_bind_all$condition == "all_conditions"] <- "All conditions"

        ## set factor levels
        BA9_bind_all$condition <- factor(BA9_bind_all$condition, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All conditions"))
        BA9_bind_all$celltype <- factor(BA9_bind_all$celltype, levels = c('L4_L6','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Astro', 'OPC','T_Cell'  ))


        ## correlation with N HVGs
        ggplot(BA9_bind_all, aes(x = HVGs, y = accuracy)) + 
        theme_bw() + 
        geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "solid") + # Add line of best fit
        geom_point(colour = "darkgrey") +
        stat_cor(method = "pearson", colour = "black", label.y = 1,) +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(colour = "black"),
            axis.text.y = element_text(colour = "black"),
            axis.title.x = element_text(face = "bold"),
            axis.title.y = element_text(face = "bold")
        ) +
        ylab("Accuracy") + 
        xlab("n HVGs") 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 3, width = 3)

        ## correlation with N Cells
        ggplot(BA9_bind_all, aes(x = n_cells, y = accuracy)) + 
        theme_bw() + 
        geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "solid") + # Add line of best fit
        geom_point(colour = "darkgrey") +
        stat_cor(method = "pearson", colour = "black", label.y = 1,) +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(colour = "black"),
            axis.text.y = element_text(colour = "black"),
            axis.title.x = element_text(face = "bold"),
            axis.title.y = element_text(face = "bold")
        ) +
        ylab("Accuracy") + 
        xlab("n cells") 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 3, width = 3)
    ##

##


################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ DNN classification accuracies for major groups (Supplemental)

## code
    ###########################
    ## SALS BA4
    ###########################
    ## code
        ## define variables
        condition = "SALS" #one ofSALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        cell_type_list <- c("Ex", 'In', 'Vasc', 'Glia')

        result_list <- list()
        
        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            test_column <- c(rep("validation", length(accuracies)), "test")
            result_df <- data.frame(accuracy = all_accuracies, test = test_column)
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
        condition = "C9ALS" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        cell_type_list <- c("Ex", 'In', 'Vasc', 'Glia')
        
        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            test_column <- c(rep("validation", length(accuracies)), "test")
            result_df <- data.frame(accuracy = all_accuracies, test = test_column)
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
        condition = "SFTLD" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        cell_type_list <- c("Ex", 'In', 'Vasc', 'Glia')


        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            test_column <- c(rep("validation", length(accuracies)), "test")
            result_df <- data.frame(accuracy = all_accuracies, test = test_column)
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
        condition = "C9FTLD" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        cell_type_list <- c("Ex", 'In', 'Vasc', 'Glia')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            test_column <- c(rep("validation", length(accuracies)), "test")
            result_df <- data.frame(accuracy = all_accuracies, test = test_column)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df

        }

        final_df_C9FTLD_BA4 <- do.call(rbind, result_list)
    ##

    ###########################
    ## all_conditions BA4 -- skip for now. 
    ###########################
    ## code
        ## define variables
        condition = "all_conditions" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        cell_type_list <- c("Ex", 'In', 'Vasc', 'Glia')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            test_column <- c(rep("validation", length(accuracies)), "test")
            result_df <- data.frame(accuracy = all_accuracies, test = test_column)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df

        }

        final_df_all_conditions_BA4 <- do.call(rbind, result_list)
    ##

    ###########################
    ## merge and plot BA4
    ###########################
    ## code
        #BA4_bind_all <- rbind(final_df_SALS_BA4, final_df_C9ALS_BA4, final_df_SFTLD_BA4, final_df_C9FTLD_BA4, final_df_all_conditions_BA4)
        BA4_bind_all <- rbind(final_df_SALS_BA4, final_df_C9ALS_BA4, final_df_SFTLD_BA4, final_df_C9FTLD_BA4)
        
        
        BA4_bind_all$condition[BA4_bind_all$condition == "all_conditions"] <- "All conditions"

        ## set factor levels
        BA4_bind_all$condition <- factor(BA4_bind_all$condition, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All conditions"))
        BA4_bind_all$celltype <- factor(BA4_bind_all$celltype, levels = c('Ex', 'In', 'Glia', 'Vasc'  ))

        ggplot(BA4_bind_all, aes(x = condition, y = accuracy, fill = condition, colour = condition)) + 
        theme_bw() + 
        geom_boxplot() +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~celltype, scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "black"))
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 3, width = 6)
    ##

    ###########################
    ## SALS BA9
    ###########################
    ## code
        ## define variables
        condition = "SALS" #one ofSALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        cell_type_list <- c("Ex", 'In', 'Vasc', 'Glia')

        result_list <- list()
        
        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            test_column <- c(rep("validation", length(accuracies)), "test")
            result_df <- data.frame(accuracy = all_accuracies, test = test_column)
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
        condition = "C9ALS" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        cell_type_list <- c("Ex", 'In', 'Vasc', 'Glia')
        
        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            test_column <- c(rep("validation", length(accuracies)), "test")
            result_df <- data.frame(accuracy = all_accuracies, test = test_column)
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
        condition = "SFTLD" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        cell_type_list <- c("Ex", 'In', 'Vasc', 'Glia')


        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            test_column <- c(rep("validation", length(accuracies)), "test")
            result_df <- data.frame(accuracy = all_accuracies, test = test_column)
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
        condition = "C9FTLD" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        cell_type_list <- c("Ex", 'In', 'Vasc', 'Glia')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            test_column <- c(rep("validation", length(accuracies)), "test")
            result_df <- data.frame(accuracy = all_accuracies, test = test_column)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df

        }

        final_df_C9FTLD_BA9 <- do.call(rbind, result_list)
    ##

    ###########################
    ## all_conditions BA9 -- skip for now. 
    ###########################
    ## code
        ## define variables
        condition = "all_conditions" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        cell_type_list <- c("Ex", 'In', 'Vasc', 'Glia')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            test_column <- c(rep("validation", length(accuracies)), "test")
            result_df <- data.frame(accuracy = all_accuracies, test = test_column)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df

        }

        final_df_all_conditions_BA9 <- do.call(rbind, result_list)
    ##

    ###########################
    ## merge and plot BA9
    ###########################
    ## code
        #BA9_bind_all <- rbind(final_df_SALS_BA9, final_df_C9ALS_BA9, final_df_SFTLD_BA9, final_df_C9FTLD_BA9, final_df_all_conditions_BA9)
        BA9_bind_all <- rbind(final_df_SALS_BA9, final_df_C9ALS_BA9, final_df_SFTLD_BA9, final_df_C9FTLD_BA9)
        
        
        BA9_bind_all$condition[BA9_bind_all$condition == "all_conditions"] <- "All conditions"

        ## set factor levels
        BA9_bind_all$condition <- factor(BA9_bind_all$condition, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All conditions"))
        BA9_bind_all$celltype <- factor(BA9_bind_all$celltype, levels = c('Ex', 'In', 'Glia', 'Vasc'  ))

        ggplot(BA9_bind_all, aes(x = condition, y = accuracy, fill = condition, colour = condition)) + 
        theme_bw() + 
        geom_boxplot() +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~celltype, scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "black"))
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 3, width = 6)
    ##

##


################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ DNN classification accuracies for major groups with equal number of cells (Supplemental)

## code
    ###########################
    ## SALS BA4
    ###########################
    ## code
        ## define variables
        condition = "SALS" #one ofSALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        cell_type_list <- c("Ex", 'In', 'Vasc', 'Glia')

        result_list <- list()
        
        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_standard_n_CombatSeq_',condition,'_',region,'_',i,'_narval.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            test_column <- c(rep("validation", length(accuracies)), "test")
            result_df <- data.frame(accuracy = all_accuracies, test = test_column)
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
        condition = "C9ALS" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        cell_type_list <- c("Ex", 'In', 'Vasc', 'Glia')
        
        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_standard_n_CombatSeq_',condition,'_',region,'_',i,'_narval.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            test_column <- c(rep("validation", length(accuracies)), "test")
            result_df <- data.frame(accuracy = all_accuracies, test = test_column)
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
        condition = "SFTLD" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        cell_type_list <- c("Ex", 'In', 'Vasc', 'Glia')


        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_standard_n_CombatSeq_',condition,'_',region,'_',i,'_narval.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            test_column <- c(rep("validation", length(accuracies)), "test")
            result_df <- data.frame(accuracy = all_accuracies, test = test_column)
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
        condition = "C9FTLD" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        cell_type_list <- c("Ex", 'In', 'Vasc', 'Glia')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_standard_n_CombatSeq_',condition,'_',region,'_',i,'_narval.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            test_column <- c(rep("validation", length(accuracies)), "test")
            result_df <- data.frame(accuracy = all_accuracies, test = test_column)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df

        }

        final_df_C9FTLD_BA4 <- do.call(rbind, result_list)
    ##

    ###########################
    ## all_conditions BA4 -- skip for now. 
    ###########################
    ## code
        ## define variables
        condition = "all_conditions" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        cell_type_list <- c("Ex", 'In', 'Vasc', 'Glia')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_standard_n_CombatSeq_',condition,'_',region,'_',i,'_narval.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            test_column <- c(rep("validation", length(accuracies)), "test")
            result_df <- data.frame(accuracy = all_accuracies, test = test_column)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df

        }

        final_df_all_conditions_BA4 <- do.call(rbind, result_list)
    ##

    ###########################
    ## merge and plot BA4
    ###########################
    ## code
        #BA4_bind_all <- rbind(final_df_SALS_BA4, final_df_C9ALS_BA4, final_df_SFTLD_BA4, final_df_C9FTLD_BA4, final_df_all_conditions_BA4)
        BA4_bind_all <- rbind(final_df_SALS_BA4, final_df_C9ALS_BA4, final_df_SFTLD_BA4, final_df_C9FTLD_BA4)
        
        
        BA4_bind_all$condition[BA4_bind_all$condition == "all_conditions"] <- "All conditions"

        ## set factor levels
        BA4_bind_all$condition <- factor(BA4_bind_all$condition, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All conditions"))
        BA4_bind_all$celltype <- factor(BA4_bind_all$celltype, levels = c('Ex', 'In', 'Glia', 'Vasc'  ))

        ggplot(BA4_bind_all, aes(x = condition, y = accuracy, fill = condition, colour = condition)) + 
        theme_bw() + 
        geom_boxplot() +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~celltype, scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "black"))
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 3, width = 6)
    ##

    ###########################
    ## SALS BA9
    ###########################
    ## code
        ## define variables
        condition = "SALS" #one ofSALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        cell_type_list <- c("Ex", 'In', 'Vasc', 'Glia')

        result_list <- list()
        
        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_standard_n_CombatSeq_',condition,'_',region,'_',i,'_narval.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            test_column <- c(rep("validation", length(accuracies)), "test")
            result_df <- data.frame(accuracy = all_accuracies, test = test_column)
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
        condition = "C9ALS" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        cell_type_list <- c("Ex", 'In', 'Vasc', 'Glia')
        
        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_standard_n_CombatSeq_',condition,'_',region,'_',i,'_narval.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            test_column <- c(rep("validation", length(accuracies)), "test")
            result_df <- data.frame(accuracy = all_accuracies, test = test_column)
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
        condition = "SFTLD" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        cell_type_list <- c("Ex", 'In', 'Vasc', 'Glia')


        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_standard_n_CombatSeq_',condition,'_',region,'_',i,'_narval.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            test_column <- c(rep("validation", length(accuracies)), "test")
            result_df <- data.frame(accuracy = all_accuracies, test = test_column)
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
        condition = "C9FTLD" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        cell_type_list <- c("Ex", 'In', 'Vasc', 'Glia')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_standard_n_CombatSeq_',condition,'_',region,'_',i,'_narval.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            test_column <- c(rep("validation", length(accuracies)), "test")
            result_df <- data.frame(accuracy = all_accuracies, test = test_column)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df

        }

        final_df_C9FTLD_BA9 <- do.call(rbind, result_list)
    ##

    ###########################
    ## all_conditions BA9 -- skip for now. 
    ###########################
    ## code
        ## define variables
        condition = "all_conditions" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        cell_type_list <- c("Ex", 'In', 'Vasc', 'Glia')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_standard_n_CombatSeq_',condition,'_',region,'_',i,'_narval.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            test_column <- c(rep("validation", length(accuracies)), "test")
            result_df <- data.frame(accuracy = all_accuracies, test = test_column)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df

        }

        final_df_all_conditions_BA9 <- do.call(rbind, result_list)
    ##

    ###########################
    ## merge and plot BA9
    ###########################
    ## code
        #BA9_bind_all <- rbind(final_df_SALS_BA9, final_df_C9ALS_BA9, final_df_SFTLD_BA9, final_df_C9FTLD_BA9, final_df_all_conditions_BA9)
        BA9_bind_all <- rbind(final_df_SALS_BA9, final_df_C9ALS_BA9, final_df_SFTLD_BA9, final_df_C9FTLD_BA9)
        
        
        BA9_bind_all$condition[BA9_bind_all$condition == "all_conditions"] <- "All conditions"

        ## set factor levels
        BA9_bind_all$condition <- factor(BA9_bind_all$condition, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All conditions"))
        BA9_bind_all$celltype <- factor(BA9_bind_all$celltype, levels = c('Ex', 'In', 'Glia', 'Vasc'  ))

        ggplot(BA9_bind_all, aes(x = condition, y = accuracy, fill = condition, colour = condition)) + 
        theme_bw() + 
        geom_boxplot() +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~celltype, scales = "free_x") +
        ylab("Accuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "black"))
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 3, width = 6)
    ##

##



