## We are running this in Narval
salloc -A def-grouleau --time=0-8 -c 1 --mem=50g

module load StdEnv/2023
module load r/4.4.0
R


## Load libraries
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
#install.packages("RColorBrewer", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
library(RColorBrewer, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")

################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ Number of highly variable genes
## code
    ###########################
    ## SALS BA4
    ###########################
    ## code
        ## define variables
        condition = "SALS" #one ofSALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval.csv')
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
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval.csv')
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
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval.csv')
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
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval.csv')
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
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval.csv')
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
        BA4_bind_all$celltype <- factor(BA4_bind_all$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))


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
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 2, width = 13)
    ##


    ###########################
    ## SALS BA9
    ###########################
    ## code
        ## define variables
        condition = "SALS" #one ofSALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval.csv')
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
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval.csv')
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
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval.csv')
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
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval.csv')
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
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval.csv')
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
        BA9_bind_all$celltype <- factor(BA9_bind_all$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))


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
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 2, width = 13)
    ##
##


################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ Number of cells
## code
    ###########################
    ## SALS BA4
    ###########################
    ## code
        ## define variables
        condition = "SALS" #one ofSALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')

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
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')

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
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')

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
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')

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
    ## all_conditions BA4
    ###########################
    ## code
        ## define variables
        condition = "all_conditions" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')

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
        BA4_bind_all$celltype <- factor(BA4_bind_all$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))


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
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "black"))
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 2, width = 13)
    ##


    ###########################
    ## SALS BA9
    ###########################
    ## code
        ## define variables
        condition = "SALS" #one ofSALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')

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
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')

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
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')

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
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')

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
    ## all_conditions BA9
    ###########################
    ## code
        ## define variables
        condition = "all_conditions" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')

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
        BA9_bind_all$celltype <- factor(BA9_bind_all$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))


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
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "black"))
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 2, width = 13)
    ##
##



################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ DNN classification accuracies (Main Figure and Supplemental)

## code
    ###########################
    ## SALS BA4
    ###########################
    ## code
        ## define variables
        condition = "SALS" #one ofSALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')

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
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')

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
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')

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
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')

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
    ## all_conditions BA4
    ###########################
    ## code
        ## define variables
        condition = "all_conditions" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')
        #cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')

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
        BA4_bind_all <- rbind(final_df_SALS_BA4, final_df_C9ALS_BA4, final_df_SFTLD_BA4, final_df_C9FTLD_BA4, final_df_all_conditions_BA4)
        BA4_bind_all$condition[BA4_bind_all$condition == "all_conditions"] <- "All conditions"

        ## set factor levels
        BA4_bind_all$condition <- factor(BA4_bind_all$condition, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All conditions"))
        BA4_bind_all$celltype <- factor(BA4_bind_all$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))


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
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 2, width = 13)
    ##

    ###########################
    ## SALS BA9
    ###########################
    ## code
        ## define variables
        condition = "SALS" #one ofSALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')

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
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')

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
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')

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
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')

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
    ## all_conditions BA9
    ###########################
    ## code
        ## define variables
        condition = "all_conditions" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')

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
        BA9_bind_all <- rbind(final_df_SALS_BA9, final_df_C9ALS_BA9, final_df_SFTLD_BA9, final_df_C9FTLD_BA9, final_df_all_conditions_BA9)
        BA9_bind_all$condition[BA9_bind_all$condition == "all_conditions"] <- "All conditions"

        ## set factor levels
        BA9_bind_all$condition <- factor(BA9_bind_all$condition, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All conditions"))
        BA9_bind_all$celltype <- factor(BA9_bind_all$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))


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
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 2, width = 13)
    ##

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
    ## DONT USE THIS ONE: Main figure accuracy heatmap: MEAN across test and validation only
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
    ## DONT USE THIS ONE: Main figure accuracy heatmap: MEDIAN across test and validation only
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
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')

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
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')

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
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')

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
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')

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
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')
        #cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')

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
        BA4_bind_all$celltype <- factor(BA4_bind_all$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))


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
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')

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
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')

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
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')

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
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')

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
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')

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
        BA9_bind_all$celltype <- factor(BA9_bind_all$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  ))


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



